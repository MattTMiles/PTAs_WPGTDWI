#!/usr/bin/env python3
"""
VERNIER-18 -- deconfusion on the MEASURED brightness ladder.
EXPLORATORY, NOT pre-registered. Clean room: imports only vernier16 / vernier_toy.

VERNIER-16 injected every source at equal SNR -- the worst possible case for a
greedy method, since nothing stands out and the first template lands between
sources. Here the ladder is built from PHOENIX's measured lambda_fore and bin
occupancy: h_k^2 ~ k^-alpha with alpha solved so h_1^2/sum h^2 = lambda_fore.

Two algorithms:
  GREEDY  -- fit one source, subtract, repeat (VERNIER-16's)
  JOINT   -- after each detection, RE-FIT ALL found amplitudes simultaneously at
             the current positions and subtract the full model. Linear in the
             amplitudes, so cheap, and it removes the bias from fitting each
             source in the presence of the others.
"""
import os
import numpy as np
from scipy.optimize import brentq

import vernier_toy as V
import vernier16 as S

REPORTS = V.REPORTS

BINS = [("2.00-2.85", 0.263, 156), ("2.85-3.71", 0.579, 46),
        ("3.71-4.56", 0.722, 22), ("4.56-5.41", 0.636, 11),
        ("5.41-6.27", 0.811, 5)]


def ladder(lam, N):
    """h_k / h_1 reproducing the measured lambda_fore at the measured occupancy."""
    f = lambda a: np.sum(np.arange(1, N + 1) ** -a) - 1.0 / lam
    a = brentq(f, 0.01, 25.0)
    return np.sqrt(np.arange(1, N + 1) ** -a), a


def inject_ladder(rng, P, Gt, amps, snr1):
    """Inject sources with the given relative amplitudes, top one at snr1."""
    n = P.shape[0]
    fp, fc = S.antenna_grid(P, Gt)
    d = np.zeros(n, complex)
    for i in range(Gt.shape[0]):
        a1 = rng.normal() + 1j * rng.normal()
        a2 = rng.normal() + 1j * rng.normal()
        m = a1 * fp[:, i] + a2 * fc[:, i]
        nn = np.sqrt(np.sum(np.abs(m) ** 2))
        d = d + m / max(nn, 1e-30) * snr1 * amps[i]
    return d + (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2)


def clean_joint(d, G, Fp, Fc, thresh, max_iter=40):
    """CLEAN with a simultaneous amplitude refit of all found sources each cycle."""
    found = []
    for _ in range(max_iter):
        if found:
            A = np.concatenate([np.stack([Fp[:, j], Fc[:, j]], axis=1)
                                for j in found], axis=1)
            M = A.T @ A
            try:
                coef = np.linalg.solve(M + 1e-10 * np.eye(M.shape[0]) * np.trace(M),
                                       A.T @ d)
            except np.linalg.LinAlgError:
                break
            r = d - A @ coef
        else:
            r = d.copy()
        F = S.fstat(r, Fp, Fc)
        j = int(np.argmax(F))
        if F[j] < thresh:
            break
        found.append(j)
    return [(j, 0.0, 0.0) for j in found]


def depth_of_clean_recovery(found, G, G_true, tol=0.15):
    """
    How far DOWN the ladder does clean recovery reach before blending sets in?
    Returns the largest k such that true sources 1..k are ALL recovered genuinely.
    """
    hit = set()
    for j, _, _ in found:
        ca = G_true @ G[j]
        i = int(np.argmax(ca))
        if np.arccos(np.clip(ca[i], -1, 1)) <= tol:
            hit.add(i)
    k = 0
    while k < G_true.shape[0] and k in hit:
        k += 1
    return k
