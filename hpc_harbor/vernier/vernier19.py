#!/usr/bin/env python3
"""
VERNIER-19 -- does the deconfusion wheel turn? EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier16 / vernier18.

The loop: find the top source, subtract THE FITTED MODEL (not the truth -- that is
the whole point), recompute, repeat. Run on a HETEROGENEOUS array (VERNIER-13's
log-normal weight spread) since subtraction and heterogeneity have never been
tested together.

Safety: every turn is scored genuine/blend against the injected truth, and the
whole loop is run on a pure-noise null.
"""
import os
import numpy as np

import vernier_toy as V
import vernier16 as S
import vernier18 as W

REPORTS = V.REPORTS
SIG_LOG_W = 0.40           # VERNIER-13 heterogeneous weight spread


def het_array(rng, n_psr, sig_log=SIG_LOG_W):
    P = V.iso_dirs(n_psr, rng)
    w = 10.0 ** rng.normal(0.0, sig_log, n_psr)
    return P, w / np.sqrt(np.mean(w ** 2))       # unit mean power


def whiten(Fp, Fc, w):
    return Fp * w[:, None], Fc * w[:, None]


def build(rng, P, w, G_all, lam, N, snr1):
    """Inject a full bin population on the measured ladder. Returns data + truth."""
    amps, alpha = W.ladder(lam, N)
    Gt = G_all[rng.choice(G_all.shape[0], N, replace=False)]
    fp, fc = S.antenna_grid(P, Gt)
    sig = np.zeros(P.shape[0], complex)
    parts = []
    for i in range(N):
        a1 = rng.normal() + 1j * rng.normal()
        a2 = rng.normal() + 1j * rng.normal()
        m = a1 * fp[:, i] + a2 * fc[:, i]
        # Normalise in the WHITENED metric: the F-statistic operates on d*w, so
        # the SNR that matters is ||m*w||, not ||m||. Normalising the unwhitened
        # norm makes the effective SNR geometry-dependent and systematically low.
        m = m / max(np.sqrt(np.sum(np.abs(m * w) ** 2)), 1e-30) * snr1 * amps[i]
        parts.append(m); sig = sig + m
    noise = (rng.normal(size=P.shape[0]) + 1j * rng.normal(size=P.shape[0])) / np.sqrt(2)
    return sig + noise / np.maximum(w, 1e-9), Gt, np.array(parts), amps


def loop(rng, P, w, G, Fpw, Fcw, thresh, lam, N, snr1, max_turn=8, tol=0.15,
         null=False):
    """The deconfusion loop. Subtracts the FITTED model each turn."""
    if null:
        d = (rng.normal(size=P.shape[0]) + 1j * rng.normal(size=P.shape[0])) / np.sqrt(2)
        d = d / np.maximum(w, 1e-9)
        Gt = np.zeros((0, 3)); parts = np.zeros((0, P.shape[0]), complex)
    else:
        d, Gt, parts, amps = build(rng, P, w, S.fib_sphere(6000), lam, N, snr1)
    P0 = float(np.sum(np.abs(parts.sum(axis=0)) ** 2)) if parts.size else 1.0
    r = d.copy()
    hist = []
    removed = set()
    for t in range(max_turn):
        F = S.fstat(r * w, Fpw, Fcw)
        j = int(np.argmax(F))
        if F[j] < thresh:
            break
        a, model = S.fit_at(r * w, Fpw[:, j], Fcw[:, j])
        model = model / np.maximum(w, 1e-9)          # back to unwhitened units
        gen = False; which = -1
        if Gt.shape[0]:
            ca = Gt @ G[j]; i = int(np.argmax(ca))
            if np.arccos(np.clip(ca[i], -1, 1)) <= tol and i not in removed:
                gen = True; which = i; removed.add(i)
        r = r - model
        # residual SIGNAL power (exclude noise): recompute from truth minus removals
        if parts.size:
            keep = [k for k in range(parts.shape[0]) if k not in removed]
            Pres = float(np.sum(np.abs(parts[keep].sum(axis=0)) ** 2)) if keep else 1e-12
            lam_now = (float(np.sum(np.abs(parts[keep[0]]) ** 2)) / Pres
                       if keep else 1.0)
        else:
            Pres, lam_now = 1.0, 0.0
        hist.append(dict(turn=t + 1, genuine=gen, which=which,
                         gain=0.5 * np.log10(P0 / max(Pres, 1e-12)),
                         lam=lam_now, F=float(F[j])))
    return hist
