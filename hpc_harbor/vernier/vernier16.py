#!/usr/bin/env python3
"""
VERNIER-16 -- does deconfusion feed the wheel? EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier8.

A REAL SEARCH, not a Fisher rank. Inject -> search -> recover -> subtract.
F-statistic over a sky grid (amplitude/phase/psi/inclination maximised
analytically as two complex coefficients), greedy CLEAN subtraction, then explicit
classification of every recovery as GENUINE / SPLIT / SPURIOUS.
"""
import os
import numpy as np
import vernier_toy as V

REPORTS = V.REPORTS


def fib_sphere(n):
    i = np.arange(n) + 0.5
    ph = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.sin(ph) * np.cos(th), np.sin(ph) * np.sin(th),
                     np.cos(ph)], axis=1)


def antenna_grid(P, G):
    """F+, Fx for every (pulsar, grid direction). Shapes (n_psr, n_grid)."""
    z = np.array([0.0, 0.0, 1.0])
    u = np.cross(z[None, :], G)
    nn = np.linalg.norm(u, axis=1, keepdims=True)
    u = np.where(nn > 1e-12, u / np.maximum(nn, 1e-30), np.array([1.0, 0, 0]))
    v = np.cross(G, u)
    pu, pv = P @ u.T, P @ v.T
    d = 1.0 + P @ G.T
    d = np.where(np.abs(d) < 1e-6, 1e-6, d)
    return (pu ** 2 - pv ** 2) / (2 * d), (pu * pv) / d


def fstat(d, Fp, Fc):
    """F-statistic over the whole grid at once."""
    Muu = np.sum(Fp * Fp, axis=0)
    Muv = np.sum(Fp * Fc, axis=0)
    Mvv = np.sum(Fc * Fc, axis=0)
    bu = Fp.T @ d
    bv = Fc.T @ d
    det = Muu * Mvv - Muv ** 2
    det = np.where(np.abs(det) < 1e-12, 1e-12, det)
    return (Mvv * np.abs(bu) ** 2 - 2 * Muv * np.real(np.conj(bu) * bv)
            + Muu * np.abs(bv) ** 2) / det


def fit_at(d, up, uc):
    """Least-squares complex coefficients of one source at one sky position."""
    A = np.stack([up, uc], axis=1)
    M = A.T @ A
    b = A.T @ d
    try:
        a = np.linalg.solve(M + 1e-12 * np.eye(2) * np.trace(M), b)
    except np.linalg.LinAlgError:
        a = np.zeros(2, complex)
    return a, A @ a


def inject(rng, P, G_true, sig):
    """Complex data: sum of M sources plus complex Gaussian noise."""
    n_psr = P.shape[0]
    Fp, Fc = antenna_grid(P, G_true)
    d = np.zeros(n_psr, complex)
    for i in range(G_true.shape[0]):
        a1 = (rng.normal() + 1j * rng.normal())
        a2 = (rng.normal() + 1j * rng.normal())
        d = d + a1 * Fp[:, i] + a2 * Fc[:, i]
    d = d + sig * (rng.normal(size=n_psr) + 1j * rng.normal(size=n_psr)) / np.sqrt(2)
    return d


def clean(d, G, Fp, Fc, thresh, max_iter=40):
    """Greedy CLEAN: find the loudest sky position, fit, subtract, repeat."""
    r = d.copy()
    found = []
    for _ in range(max_iter):
        F = fstat(r, Fp, Fc)
        j = int(np.argmax(F))
        if F[j] < thresh:
            break
        a, model = fit_at(r, Fp[:, j], Fc[:, j])
        if not np.all(np.isfinite(model)):
            break
        r = r - model
        found.append((j, float(F[j]), float(np.sum(np.abs(model) ** 2))))
    return found, r


def classify(found, G, G_true, tol=0.35):
    """GENUINE (first match to a true source) / SPLIT (repeat) / SPURIOUS (no match)."""
    claimed = set()
    gen = spl = spu = 0
    for j, _, _ in found:
        cosang = G_true @ G[j]
        i = int(np.argmax(cosang))
        if np.arccos(np.clip(cosang[i], -1, 1)) <= tol:
            if i in claimed:
                spl += 1
            else:
                claimed.add(i); gen += 1
        else:
            spu += 1
    return gen, spl, spu
