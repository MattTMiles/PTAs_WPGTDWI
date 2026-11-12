#!/usr/bin/env python3
"""
VERNIER-12 -- pairwise, properly. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier6 / vernier6b / vernier7 / vernier8.

MECHANISM. Every pulsar's distance inference depends on the SAME shared source
parameters (f, Mc, sky). An error df in the shared frequency shifts pulsar a's
pulsar-term phase by
        d(phi_a) = 2 pi df tau_a
-- proportional to that pulsar's tau, but driven by ONE common df. Handled
per-pulsar the error is marginalised independently N times and the smear is paid
in full; handled jointly across a group of m pulsars the common df is constrained
m-fold and the residual smear is sigma_phi / sqrt(m).

SHARP PREDICTION: the pairwise benefit must scale with sigma_f and VANISH when the
source parameters are known exactly. If it does not vanish, the 105 nats came from
somewhere else.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier_toy as V
import vernier6 as X
import vernier6b as B
import vernier7 as S
import vernier8 as E

YR, REPORTS = V.YR, V.REPORTS
V.PTERM_STATIC = False


def sigma_phi_of(cfg, sigma_f_rel, group_m):
    """Residual shared-parameter phase smear when m pulsars are fitted jointly."""
    tau = V.tau_of(cfg["L_kpc"] * 1000.0, 1.0)
    s = 2 * np.pi * sigma_f_rel * cfg["f_gw_top"] * tau
    return s / np.sqrt(max(group_m, 1))


def contrast_cell(rng, cfg, sigma_f_rel, group_m, n_trial=60, seeds=None):
    """Contrast and its scrambled null at a given (sigma_f, group size)."""
    sp = sigma_phi_of(cfg, sigma_f_rel, group_m)
    if seeds is None:
        seeds = list(rng.integers(0, 2 ** 31, n_trial))
    cs, cn, qs, qn = [], [], [], []
    for sd in seeds:
        r2 = np.random.default_rng(sd)
        comps, _, dL = X.build_comps(r2, cfg)
        Lstar = cfg["L_kpc"] * 1000.0
        Lg = Lstar + np.linspace(-3 * cfg["prior_pc"], 3 * cfg["prior_pc"], 16000)
        for scr, C, Q in ((False, cs, qs), (True, cn, qn)):
            r3 = np.random.default_rng(sd + (10 ** 7 if scr else 0))
            lnl = X.lnL_curve(Lg, Lstar, comps, r3, scramble=scr, sigma_phi=sp)
            m = X.metrics(lnl, Lg, Lstar, dL, cfg["prior_pc"])
            C.append(m["contrast"]); Q.append(m["q_true"])
    f = lambda v: float(np.median(np.array(v)[np.isfinite(v)]))
    return dict(sig=f(cs), null=f(cn), q_sig=f(qs), q_null=f(qn), sigma_phi=sp)


def run_wheel_joint(seed, cfg, f_E=20e-9, rho_P=6.0, sigma_f_rel=3e-4, group_m=None,
                    npts=7000, max_turn=30, scramble=False):
    """
    VERNIER-9's per-pulsar-distance wheel WITH the shared-source-parameter smear
    that 6b/7/8/9/10 silently dropped. group_m = pulsars fitted jointly
    (1 = per-pulsar, n_psr = full joint).
    """
    import vernier9 as N
    rng = np.random.default_rng(seed)
    n_psr, n_src = cfg["n_psr"], cfg["n_src"]
    if group_m is None:
        group_m = n_psr
    p, Ltrue = N.draw_array(rng, n_psr)
    s_ = V.iso_dirs(n_src, rng)
    f_gw = f_E * 10 ** rng.uniform(-0.15, 0.0, n_src)
    tc = V.t_merge(f_gw, N.MC_FID)
    cmu = p @ s_.T
    omc = np.clip(1.0 - cmu, 5e-2, None)
    ant = np.sqrt(np.clip(1.0 - cmu ** 2, 1e-6, None))
    amp = 10.0 ** S.LAD[:n_src]
    r_a = (1.0 + V.tau_of(Ltrue[:, None], omc) / tc[None, :]) ** 0.125
    # per-pulsar shared-parameter smear, reduced by the joint group size
    sig_phi = 2 * np.pi * sigma_f_rel * f_E * V.tau_of(Ltrue[:, None], omc) \
        / np.sqrt(max(group_m, 1))
    dec = np.exp(-0.5 * sig_phi ** 2)
    grids = [Ltrue[a] + np.linspace(-3 * cfg["prior_pc"], 3 * cfg["prior_pc"], npts)
             for a in range(n_psr)]
    dL_ref = float(np.median([V.dL_of(f_gw[0], omc[a, 0]) for a in range(n_psr)]))
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
            Ls = max(Ltrue[a] + off[k], 50.0)
            phi, env = X.phase_env(g, Ls, f_gw[k] / 2.0, tc[k], omc[a, k],
                                   cfg["T"], int(n))
            env = env * dec[a, k]                    # <-- the smear, now present
            rr = rho * np.sqrt(wq)
            nc, nsn = noise[(a, k)]
            tot += rr ** 2 * (env * np.cos(phi) - 1.0) \
                + rr * env * (nc * np.cos(phi) + nsn * np.sin(phi))
        return tot

    def acc(a, det):
        tot = np.zeros_like(grids[a])
        for k in det:
            if (a, k) not in Lc:
                Lc[(a, k)] = lnl_one(a, k)
            tot = tot + Lc[(a, k)]
        return tot

    rhoE = ant * amp[None, :]
    rhoE = rhoE * cfg["snr0"] / np.sqrt(np.sum(rhoE[:, 0] ** 2))
    rhoP = rhoE * r_a
    det = [0]
    H = {"w_true": [], "K_eff": [], "gain": []}
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
        Rk = np.array([B.coherence(W[a], X.phase_env(grids[a], Ltrue[a], f_gw[k] / 2.0,
                                                     tc[k], omc[a, k], cfg["T"], 2)[0],
                                   grids[a]) for a in range(n_psr)])
        s0 = np.sqrt(np.sum(rhoE[:, k] ** 2))
        sn = np.sqrt(s0 ** 2 + np.sum(Rk ** 2 * rhoP[:, k] ** 2))
        H["gain"].append(float(np.log10(sn / s0)))
        if sn < cfg["det_thresh"]:
            break
        det.append(k)
    H = {kk: np.array(vv) for kk, vv in H.items()}
    H["final_rank"] = int(max(det) + 1)
    return H
