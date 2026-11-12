#!/usr/bin/env python3
"""
VERNIER-14 -- reproduce them on their own terms. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier6.

Xiao et al. configuration, used exactly as stated:
  85 pulsars, 20 yr, 2-week cadence, sigma_n = 50 / 100 ns
  Gaussian distance priors, sigma = 10 pc at 1 kpc
  Mc = 5e9 Msun, d_L = 1 Gpc, f_GW in 3-15 nHz, sky positions FIXED
  N = 5 sources combined; Mc, f, d_L marginalised
"""
import os
import numpy as np
import vernier_toy as V
import vernier6 as X

YR, C_PC_S, REPORTS = V.YR, V.C_PC_S, V.REPORTS
V.PTERM_STATIC = False
T_SUN = 4.925490947e-6
PC = 3.0856775814913673e16
C = 299792458.0

XIAO = dict(n_psr=85, T=20 * YR, cad=14.0, sig_n=50e-9, prior_pc=10.0,
            mc=5e9, dL_Gpc=1.0, n_src=5, f_lo=3e-9, f_hi=15e-9, L_kpc=1.0)
YUPAN = dict(n_psr=20, T=30 * YR, cad=14.0, sig_n=20e-9, prior_pc=10.0,
             mc=5e9, dL_Gpc=1.0, n_src=8, f_lo=3e-9, f_hi=15e-9, L_kpc=1.0)


def h0_of(mc, dL_Gpc, f):
    return 2.0 * (T_SUN * mc) ** (5 / 3.) * (np.pi * f) ** (2 / 3.) \
        / (dL_Gpc * 1e9 * PC / C)


def rho_of(cfg, f):
    """Array and per-pulsar matched-filter SNR for one source."""
    R = h0_of(cfg["mc"], cfg["dL_Gpc"], f) / (2 * np.pi * f)
    Sw = 2 * cfg["sig_n"] ** 2 * cfg["cad"] * 86400.0
    ra = np.sqrt((1 / 3.) * cfg["n_psr"] * R ** 2 * cfg["T"] / (2 * Sw))
    return ra, ra / np.sqrt(cfg["n_psr"])


def sigma_phi(cfg, f, group_m=1):
    """Shared-frequency-error phase smear; sky is FIXED so only f (and Mc) smear."""
    ra, _ = rho_of(cfg, f)
    sf = 1.0 / (2 * np.pi * cfg["T"] * ra)
    tau = V.tau_of(cfg["L_kpc"] * 1000.0, 1.0)
    return 2 * np.pi * sf * tau / np.sqrt(max(group_m, 1))


def pair_posterior(rng, cfg, f, sig_phi_pair, rho_p, D_true, npts=4000, n_src=5):
    """
    Posterior over ONE pulsar's distance from a PAIR analysis, with fringe structure,
    noise, and the residual shared-parameter smear. Returns (D_est, half_width_68).
    """
    half = 3.0 * cfg["prior_pc"]
    D = D_true + np.linspace(-half, half, npts)
    tc = V.t_merge(f, 9.7)
    lnl = np.zeros_like(D)
    for _ in range(n_src):
        fk = f * 10 ** rng.uniform(-0.15, 0.15)
        phi, env = X.phase_env(D, D_true, fk / 2.0, tc, 1.0, cfg["T"], 2)
        env = env * np.exp(-0.5 * sig_phi_pair ** 2)
        nc, ns = rng.normal(), rng.normal()
        lnl += rho_p ** 2 * (env * np.cos(phi) - 1.0) \
            + rho_p * env * (nc * np.cos(phi) + ns * np.sin(phi))
    pr = np.exp(-0.5 * ((D - D_true) / cfg["prior_pc"]) ** 2)
    post = np.exp(lnl - lnl.max()) * pr
    post /= np.trapezoid(post, D)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (post[1:] + post[:-1]) * np.diff(D))])
    cdf /= cdf[-1]
    lo, med, hi = np.interp([0.16, 0.5, 0.84], cdf, D)
    return float(med), float(0.5 * (hi - lo))


def selection_trial(rng, cfg, f, D_true, n_partner=84, scramble=False):
    """Xiao's procedure: build all pair posteriors, ADOPT THE NARROWEST."""
    ra, rp = rho_of(cfg, f)
    ests, hws = [], []
    for _ in range(n_partner):
        # partner quality varies with its antenna coupling -> pair constrains f
        # better or worse -> different residual smear and different rho
        q = np.clip(rng.uniform(0.3, 1.0), 0.05, None)
        sp = sigma_phi(cfg, f, 2) / q
        Dt = D_true + (rng.uniform(-30, 30) if scramble else 0.0)
        e, h = pair_posterior(rng, cfg, f, sp, rp * q, Dt)
        ests.append(e); hws.append(h)
    ests, hws = np.array(ests), np.array(hws)
    j = int(np.argmin(hws))                      # <-- the selection
    return dict(sel_est=ests[j], sel_hw=hws[j],
                med_hw=float(np.median(hws)),
                sel_err=abs(ests[j] - D_true),
                med_err=float(np.median(np.abs(ests - D_true))),
                cov_sel=bool(abs(ests[j] - D_true) <= hws[j]),
                cov_typ=bool(abs(ests[0] - D_true) <= hws[0]))
