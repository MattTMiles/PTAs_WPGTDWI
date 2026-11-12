#!/usr/bin/env python3
"""
VERNIER-9 -- raise the ceiling. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier6 / vernier6b / vernier7.

VERNIER-8 established the budget is exactly log10 sqrt(1 + r^2) with
    r = (f_P/f_E)^(-1/3)                (residual amplitude ratio)
and with the exact Peters chirp f_P/f_E = (1 + tau/t_c)^(-3/8), so

    *** r = (1 + tau/t_c)^(1/8) ***

tau = L (1 - cos mu)/c scales with PULSAR DISTANCE, so r is a per-pulsar quantity,
not an array constant. Distant pulsars sit further up the chirp and carry more
pulsar-term residual. The array ceiling is the rho_E^2-weighted quantity

    ceiling = log10 sqrt( 1 + <r^2>_w ),   <r^2>_w = sum_a r_a^2 rho_E,a^2 / sum_a rho_E,a^2
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V
import vernier6b as B
import vernier7 as S
import vernier8 as E

YR, C_PC_S, REPORTS = V.YR, V.C_PC_S, V.REPORTS
V.PTERM_STATIC = False

L_LO, L_HI = 300.0, 5000.0          # pc -- IPTA-like spread, 0.3-5 kpc
MC_FID = 9.5
F_FID = 25e-9


def r_of(L_pc, omc, f_E, log10mc):
    """Per-pulsar residual amplitude ratio r = (1 + tau/t_c)^(1/8)."""
    tc = V.t_merge(f_E, log10mc)
    tau = V.tau_of(L_pc, omc)
    return (1.0 + tau / tc) ** 0.125, tau, tc


def draw_array(rng, n_psr, log_uniform=True):
    L = (10 ** rng.uniform(np.log10(L_LO), np.log10(L_HI), n_psr) if log_uniform
         else rng.uniform(L_LO, L_HI, n_psr))
    p = V.iso_dirs(n_psr, rng)
    return p, L


def array_ceiling(rng, f_E, log10mc=MC_FID, n_psr=200, n_real=60):
    """rho_E^2-weighted array ceiling, and the single-distance comparison."""
    out = []
    for _ in range(n_real):
        p, L = draw_array(rng, n_psr)
        s = V.iso_dirs(1, rng)[0]
        cmu = p @ s
        omc = np.clip(1.0 - cmu, 5e-2, None)
        ant2 = np.clip(1.0 - cmu ** 2, 1e-6, None)      # rho_E^2 weight
        r, tau, tc = r_of(L, omc, f_E, log10mc)
        w = ant2 / ant2.sum()
        r2w = float(np.sum(w * r ** 2))
        # contribution of the most distant quartile
        q = L >= np.percentile(L, 75)
        frac_q = float(np.sum(w[q] * r[q] ** 2) / max(r2w, 1e-12))
        out.append([r2w, np.log10(np.sqrt(1 + r2w)), float(np.median(r)),
                    float(r.min()), float(r.max()), frac_q, float(np.median(tau / tc))])
    return np.array(out)


def run_wheel_dist(seed, cfg, f_E=F_FID, rho_P=6.0, npts=7000, max_turn=30,
                   scramble=False):
    """
    VERNIER-7's wheel with a PER-PULSAR distance (and hence per-pulsar tau, r).
    Each pulsar gets its own true distance and its own likelihood grid.
    """
    rng = np.random.default_rng(seed)
    n_psr, n_src = cfg["n_psr"], cfg["n_src"]
    p, Ltrue = draw_array(rng, n_psr)
    s = V.iso_dirs(n_src, rng)
    f_gw = f_E * 10 ** rng.uniform(-0.15, 0.0, n_src)
    tc = V.t_merge(f_gw, MC_FID)
    cmu = p @ s.T
    omc = np.clip(1.0 - cmu, 5e-2, None)
    ant = np.sqrt(np.clip(1.0 - cmu ** 2, 1e-6, None))
    amp = 10.0 ** S.LAD[:n_src]
    r_a = (1.0 + V.tau_of(Ltrue[:, None], omc) / tc[None, :]) ** 0.125   # per (psr, src)

    grids = [Ltrue[a] + np.linspace(-3 * cfg["prior_pc"], 3 * cfg["prior_pc"], npts)
             for a in range(n_psr)]
    dL_ref = float(np.median([V.dL_of(2 * f_gw[0] / 2, omc[a, 0]) for a in range(n_psr)]))
    off = {k: (rng.uniform(-500.0, 500.0) if scramble else 0.0) for k in range(n_src)}
    noise = {(a, k): (rng.normal(), rng.normal())
             for a in range(n_psr) for k in range(n_src)}
    Lc = {}

    def lnl_one(a, k):
        g = grids[a]
        tot = np.zeros_like(g)
        ns = np.arange(1, cfg["n_harm"] + 1)
        gw = np.array([V.pm_g(int(n), cfg["ecc"]) for n in ns]); gw /= gw.sum()
        rho = rho_P * ant[a, k] / np.mean(ant) * amp[k]
        for n, wq in zip(ns, gw):
            if wq < 1e-3:
                continue
            # clip: a scramble offset must not push the trial distance negative
            Ls = max(Ltrue[a] + off[k], 50.0)
            phi, env = B.X.phase_env(g, Ls, f_gw[k] / 2.0, tc[k],
                                     omc[a, k], cfg["T"], int(n))
            rr = rho * np.sqrt(wq)
            nc, nsn = noise[(a, k)]
            tot += rr ** 2 * (env * np.cos(phi) - 1.0) \
                + rr * env * (nc * np.cos(phi) + nsn * np.sin(phi))
        return tot

    def acc(a, det):
        g = grids[a]
        tot = np.zeros_like(g)
        for k in det:
            if (a, k) not in Lc:
                Lc[(a, k)] = lnl_one(a, k)
            tot = tot + Lc[(a, k)]
        return tot

    rhoE = ant * amp[None, :]
    nrm = cfg["snr0"] / np.sqrt(np.sum(rhoE[:, 0] ** 2))
    rhoE = rhoE * nrm
    rhoP = rhoE * r_a                              # PER-PULSAR r

    det = [0]
    H = {"rank": [], "gain": [], "R": [], "w_true": [], "K_eff": [], "step": []}
    for t in range(max_turn):
        W = [B.soft_state(acc(a, det), grids[a], cfg["prior_pc"], Ltrue[a])
             for a in range(n_psr)]
        kw = [B.keff_wtrue(W[a], grids[a], Ltrue[a], dL_ref) for a in range(n_psr)]
        H["K_eff"].append(float(np.median([x[0] for x in kw])))
        H["w_true"].append(float(np.median([x[1] for x in kw])))
        avail = [k for k in range(n_src) if k not in det]
        if not avail:
            break
        k = avail[0]
        Rk = np.array([B.coherence(W[a],
                                   B.X.phase_env(grids[a], Ltrue[a], f_gw[k] / 2.0,
                                                 tc[k], omc[a, k], cfg["T"], 2)[0],
                                   grids[a]) for a in range(n_psr)])
        s0 = np.sqrt(np.sum(rhoE[:, k] ** 2))
        sn = np.sqrt(s0 ** 2 + np.sum(Rk ** 2 * rhoP[:, k] ** 2))
        H["R"].append(float(np.median(Rk)))
        H["gain"].append(float(np.log10(sn / s0)))
        H["rank"].append(k + 1)
        H["step"].append(E.step_needed(k + 1))
        if sn < cfg["det_thresh"]:
            break
        det.append(k)
    H = {kk: np.array(vv) for kk, vv in H.items()}
    H["final_rank"] = int(max(det) + 1)
    return H
