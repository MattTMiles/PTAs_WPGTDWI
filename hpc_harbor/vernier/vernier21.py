#!/usr/bin/env python3
"""
VERNIER-21 -- is lambda a resolution parameter? EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier16 / vernier18.

HYPOTHESIS. Certified distances let the pulsar term resolve frequency WITHIN a
Fourier bin (the Earth term has no frequency dependence inside a bin at all), so
the effective resolution element is narrower, fewer sources fall in it, and the
effective lambda rises.

TESTED DIRECTLY, not argued: lambda_eff is computed from the actual TEMPLATE
OVERLAPS. An interloper only confuses the target to the extent its template
projects onto the target's:

    lambda_eff = h_k^2 / ( h_k^2 + sum_{i!=k} h_i^2 * |<t_k|t_i>|^2 )

with the overlap normalised. Earth-only vs Earth+pulsar-term templates.
"""
import os
import numpy as np

import vernier_toy as V
import vernier16 as S
import vernier18 as W

REPORTS = V.REPORTS
T_OBS = 20 * V.YR
DF_BIN = 1.0 / T_OBS


def templates(rng, P, L, Gt, freqs, use_pterm):
    """Complex template per source, stacked over pulsars (+ pulsar-term block)."""
    fp, fc = S.antenna_grid(P, Gt)
    psi = rng.uniform(0, np.pi, Gt.shape[0])
    ci = rng.uniform(-1, 1, Gt.shape[0])
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    Fp = fp * c2[None, :] + fc * s2[None, :]
    Fc = -fp * s2[None, :] + fc * c2[None, :]
    E = Fp * ((1 + ci ** 2) / 2)[None, :] + 1j * Fc * ci[None, :]
    if not use_pterm:
        return E
    cmu = P @ Gt.T
    omc = np.clip(1.0 - cmu, 5e-2, None)
    tau = V.tau_of(L[:, None], omc)                     # (n_psr, n_src)
    Pt = E * np.exp(-2j * np.pi * freqs[None, :] * tau)
    return np.concatenate([E, Pt], axis=0)


def lambda_eff(rng, n_psr, N, lam, use_pterm, f0=5.8e-9, n_real=40):
    """Effective lambda from measured template overlaps."""
    amps, _ = W.ladder(lam, N)
    out = []
    for _ in range(n_real):
        P = V.iso_dirs(n_psr, rng)
        L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), n_psr)
        Gt = S.fib_sphere(6000)[rng.choice(6000, N, replace=False)]
        freqs = f0 + rng.uniform(-0.5, 0.5, N) * DF_BIN   # spread across ONE bin
        Tm = templates(rng, P, L, Gt, freqs, use_pterm)
        nrm = np.sqrt(np.sum(np.abs(Tm) ** 2, axis=0))
        Tn = Tm / np.maximum(nrm, 1e-30)
        ov = np.abs(Tn.conj().T @ Tn) ** 2               # (N, N) overlap^2
        h2 = amps ** 2
        num = h2[0]
        den = h2[0] + np.sum(h2[1:] * ov[0, 1:])
        out.append(num / den)
    return float(np.median(out))
