#!/usr/bin/env python3
"""
VERNIER-20 -- run them together, honestly. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier16 / vernier18 / vernier19.

VERNIER-19 ADDED 0.244 (heterogeneity) + 0.308 (subtraction) = 0.552 dex and
refused to stand behind it. Here the three mechanisms are simulated TOGETHER in
one self-consistent noise budget, so the interactions are measured:

  target k's SNR:   rho_k^2 = sum_a R_ak^2 (1 + r^2 c_a) / (1 + sum_{i!=k,unsub} R_ai^2)
  certification:    pulsar a certifies iff its own pulsar-term SNR >= 2, which
                    itself depends on the confusion -- so subtraction RAISES
                    f_cert, and heterogeneity changes WHICH pulsars clear it.

Nothing is added by hand; every configuration is evaluated from the same budget.
"""
import os
import numpy as np

import vernier_toy as V
import vernier16 as S
import vernier18 as W
import vernier19 as L

REPORTS = V.REPORTS
R_AMP = 1.196              # pulsar-term / Earth-term residual ratio (VERNIER-8)
CERT_BAR = 2.0             # per-pulsar rho_P for certification (VERNIER-1)


def scenario(rng, n_psr, lam, N, snr1, hetero, sig_log=0.40):
    """One bin: array weights, per-(pulsar,source) response amplitudes."""
    P = V.iso_dirs(n_psr, rng)
    w = (10.0 ** rng.normal(0, sig_log, n_psr)) if hetero else np.ones(n_psr)
    w = w / np.sqrt(np.mean(w ** 2))
    amps, _ = W.ladder(lam, N)
    Gt = S.fib_sphere(6000)[rng.choice(6000, N, replace=False)]
    fp, fc = S.antenna_grid(P, Gt)
    ant = np.sqrt(fp ** 2 + fc ** 2)                     # (n_psr, N)
    R = ant * w[:, None] * amps[None, :]
    R = R / np.sqrt(np.sum(R[:, 0] ** 2)) * snr1          # top source at snr1
    return R


def snr_of(R, k, subtracted, use_pterm, joint=True, n_psr_eff=None):
    """
    Self-consistent SNR of source k given which sources are subtracted and
    whether certified pulsar terms are used.
    """
    n_psr, N = R.shape
    others = [i for i in range(N) if i != k and i not in subtracted]
    conf = R[:, others] ** 2 if others else np.zeros((n_psr, 1))
    noise = 1.0 + conf.sum(axis=1)                        # white(=1) + confusion
    if not use_pterm:
        return float(np.sqrt(np.sum(R[:, k] ** 2 / noise))), 0.0
    # per-pulsar certification, evaluated in the SAME budget
    rho_p = R_AMP * R[:, k] / np.sqrt(noise)
    cert = rho_p >= CERT_BAR
    if not joint:
        cert = cert & False                               # per-pulsar treatment fails
    fc = (np.sum(R[cert, k] ** 2) / max(np.sum(R[:, k] ** 2), 1e-30))
    s2 = np.sum(R[:, k] ** 2 / noise) + R_AMP ** 2 * np.sum(R[cert, k] ** 2 / noise[cert])
    return float(np.sqrt(s2)), float(fc)
