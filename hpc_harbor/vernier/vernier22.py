#!/usr/bin/env python3
"""
VERNIER-22 -- the small array with a soft bar. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier6b / vernier16 / vernier18 / vernier21.
"""
import os
import numpy as np

import vernier_toy as V
import vernier6b as B
import vernier16 as S
import vernier18 as W
import vernier21 as R21

REPORTS = V.REPORTS
R_AMP = 1.196


def lambda_eff_searched(rng, n_psr, N, lam, n_real=60, f0=5.8e-9):
    """
    lambda_eff when the target's own position is SEARCHED rather than known.

    V21 evaluated overlaps at the true positions, which is an upper bound. Here the
    target template is placed at the best-fit sky position found by maximising the
    F-statistic on the full (target + interlopers + noise) data. The purity of that
    recovery -- what fraction of the fitted template's power is genuinely the target
    -- is the operational lambda_eff a real search sees.
    """
    amps, _ = W.ladder(lam, N)
    G = S.fib_sphere(2000)
    out = []
    for _ in range(n_real):
        P = V.iso_dirs(n_psr, rng)
        Fp, Fc = S.antenna_grid(P, G)
        Gt = S.fib_sphere(6000)[rng.choice(6000, N, replace=False)]
        fpt, fct = S.antenna_grid(P, Gt)
        parts = []
        d = np.zeros(n_psr, complex)
        for i in range(N):
            a1 = rng.normal() + 1j * rng.normal()
            a2 = rng.normal() + 1j * rng.normal()
            m = a1 * fpt[:, i] + a2 * fct[:, i]
            m = m / max(np.sqrt(np.sum(np.abs(m) ** 2)), 1e-30) * amps[i]
            parts.append(m); d = d + m
        d = d + (rng.normal(size=n_psr) + 1j * rng.normal(size=n_psr)) / np.sqrt(2) * 0.05
        j = int(np.argmax(S.fstat(d, Fp, Fc)))
        a, model = S.fit_at(d, Fp[:, j], Fc[:, j])
        pm = float(np.sum(np.abs(model) ** 2))
        if pm < 1e-30:
            continue
        # how much of the FITTED template is genuinely the top source?
        proj = abs(np.vdot(model, parts[0])) ** 2 / (pm * max(float(np.sum(np.abs(parts[0]) ** 2)), 1e-30))
        out.append(proj)
    return float(np.median(out))


def soft_wheel(seed, n_psr, lam, N, snr_per_psr_std=1.0, prior_pc=30.0,
               npts=6000, max_turn=6, scramble=False, subtract=True):
    """
    The loop with a SOFT bar: nothing is thresholded. Pulsar state is the full
    fringe weight function; coherence is the resultant length; progress is K_eff.
    Per-pulsar SNR comes from a self-consistent confusion budget with the top
    source subtracted.
    """
    rng = np.random.default_rng(seed)
    amps, _ = W.ladder(lam, N)
    P = V.iso_dirs(n_psr, rng)
    L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), n_psr)
    Gt = S.fib_sphere(6000)[rng.choice(6000, N, replace=False)]
    fp, fc = S.antenna_grid(P, Gt)
    ant = np.sqrt(fp ** 2 + fc ** 2)
    Rmat = ant * amps[None, :]
    Rmat = Rmat / np.sqrt(np.mean(Rmat[:, 0] ** 2)) * snr_per_psr_std
    f0 = 5.8e-9
    tc = V.t_merge(f0, 9.5)
    omc = np.clip(1.0 - P @ Gt.T, 5e-2, None)
    det = [0]
    H = {"K_eff": [], "w_true": [], "R": [], "f_cert": [], "rho_p": []}
    for t in range(max_turn):
        sub = set(det[:1]) if subtract else set()
        # per-pulsar pulsar-term SNR of the next target, in the confusion budget
        k = min(len(det), N - 1)
        others = [i for i in range(N) if i != k and i not in sub]
        noise = 1.0 + (Rmat[:, others] ** 2).sum(axis=1) if others else np.ones(n_psr)
        rho_p = R_AMP * Rmat[:, k] / np.sqrt(noise)
        # SOFT weights: full fringe posterior per pulsar, nothing thresholded
        Ke, wt, Rc = [], [], []
        for a in range(n_psr):
            Lstar = L[a]
            Lg = Lstar + np.linspace(-3 * prior_pc, 3 * prior_pc, npts)
            off = rng.uniform(-500, 500) if scramble else 0.0
            comps = [B.X.Comp(f0 / 2, tc, omc[a, j], 2, rho_p[a], 20 * V.YR, j)
                     for j in det]
            lnl = B.X.lnL_curve(Lg, max(Lstar + off, 50.0), comps, rng)
            w = B.soft_state(lnl, Lg, prior_pc, Lstar)
            dL = V.dL_of(f0, omc[a, k])
            ke, wtr, _ = B.keff_wtrue(w, Lg, Lstar, dL)
            phi = B.X.phase_env(Lg, Lstar, f0 / 2, tc, omc[a, k], 20 * V.YR, 2)[0]
            Ke.append(ke); wt.append(wtr); Rc.append(B.coherence(w, phi, Lg))
        Rc = np.array(Rc)
        fcert = float(np.sum(Rc ** 2 * rho_p ** 2) / max(np.sum(rho_p ** 2), 1e-30))
        H["K_eff"].append(float(np.median(Ke)))
        H["w_true"].append(float(np.median(wt)))
        H["R"].append(float(np.median(Rc)))
        H["f_cert"].append(fcert)
        H["rho_p"].append(float(np.median(rho_p)))
        if k >= N - 1:
            break
        det.append(k)
    return {kk: np.array(vv) for kk, vv in H.items()}
