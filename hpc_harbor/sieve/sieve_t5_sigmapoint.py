"""SIEVE-C T5 -- SIGMA-POINT E-STEP, GATE HALF.  Zero-noise / Asimov only.

Question: the E-step currently evaluates the per-pulsar FRINGE POSTERIOR at the
source POINT ESTIMATE.  Should it instead be AVERAGED over the source belief?

  gate    : as the sigma-point width -> 0 the sigma-point E-step must reproduce
            the point E-step bit-comparably.
  readout : the fringe-posterior shift at REALISTIC-class source uncertainty --
            the STORY.md S4.2.3/S4.2.4 widths (sigma(log10 Mc) = 0.501 / 0.864 /
            0.866 dex; log10 f pinned by the Earth term to 162-389x the prior).

WHY THE SOURCE MUST CHIRP HERE.  S4.2.4: *95-100 % of the pulsar-term phase
wander is log10_mc*.  A non-evolving toy has no mc at all, so this test is run
on a GR-driven circular binary with the EXACT closed-form evolution

    f(t) = [ f0^(-8/3) - k t ]^(-3/8),  k = (256/5) pi^(8/3) (G Mc / c^3)^(5/3)

evaluated at the Earth time t and at the RETARDED time t - dt_p.  dt_p ~ 1e11 s
(~3 kyr), so the pulsar term is a genuinely earlier, lower-frequency epoch of the
same binary and the pulsar-term phase carries mc.

SUB-FRINGE ARTEFACT GUARD (the banked hazard: a fringe-CENTRE grid pinned to
u = 0 fabricates fringe evidence).  Every fringe cell here is scored at its
PROFILED sub-fringe offset -- a 1-D maximisation inside the cell -- and the
fringe ladder is anchored on the PRIOR MEAN, never on the injected truth.

DEFERRED TO SIEVE-A (stated, not attempted here): calibration-against-truth over
NOISY realisations.  cronus is RECOMPUTED-UNSAFE for noise draws; this file
contains no `rng` call that touches a datum.
"""

import os
import sys
import json
import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sieve_toy as S

OUT = os.path.join(S.ROOT, "reports", "sieve_t5_sigmapoint.npz")

G = 6.67430e-11
MSUN = 1.98892e30

# the injected source (one source, Asimov, loud)
SRC = dict(costheta=0.35, phi=1.90, log10_h=-13.0, log10_f=-8.00,
           cosi=0.40, psi=0.70, phi0=1.10, log10_mc=9.0)

# STORY.md S4.2.3 / S4.2.4 realistic-class belief widths
SIG_MC_DEX = [0.501, 0.864, 0.866]
SIG_F_DEX = 0.144 / 389.0           # Earth-term-pinned (gain 389x on the U(-8,-7.5) prior)
NFRINGE = 8                          # +-8 fringes about the prior mean


def kchirp(log10_mc):
    mc = (10.0 ** log10_mc) * MSUN
    return (256.0 / 5.0) * np.pi ** (8.0 / 3.0) * (G * mc / S.C_LIGHT ** 3) ** (5.0 / 3.0)


def chirp_phase(t, f0, k, phi0):
    """Phi(t) = 2 pi int_0^t f dt' for f = [f0^{-8/3} - k t]^{-3/8}."""
    A = f0 ** (-8.0 / 3.0)
    arg = np.clip(A - k * t, 1e-30, None)
    return phi0 + (16.0 * np.pi / (5.0 * k)) * (A ** (5.0 / 8.0) - arg ** (5.0 / 8.0))


def chirp_freq(t, f0, k):
    A = f0 ** (-8.0 / 3.0)
    return np.clip(A - k * t, 1e-30, None) ** (-3.0 / 8.0)


def residuals_chirp(spine, src, L_m=None):
    """Zero-noise CW residuals for a GR-driven circular binary, Earth + pulsar
    term, with the pulsar term at the EXACT retarded epoch."""
    t = spine["toas"]
    phat = spine["phat"]
    L = spine["L_m"] if L_m is None else np.asarray(L_m, float)
    shat = S.sky_unit(src["costheta"], src["phi"])
    Fp, Fx = S.antenna(shat, phat)
    f0 = 10.0 ** src["log10_f"]
    k = kchirp(src["log10_mc"])
    h0 = 10.0 ** src["log10_h"]
    cosi, psi, phi0 = src["cosi"], src["psi"], src["phi0"]
    Ap = h0 * (1.0 + cosi ** 2) / 2.0
    Ax = h0 * cosi
    c2, s2 = np.cos(2.0 * psi), np.sin(2.0 * psi)
    Fp_r = Fp * c2 + Fx * s2
    Fx_r = -Fp * s2 + Fx * c2
    dt = (L / S.C_LIGHT) * (1.0 - phat @ shat)

    def wave(tt):
        ph = chirp_phase(tt, f0, k, phi0)
        fe = chirp_freq(tt, f0, k)
        return (Ap / (2 * np.pi * fe)) * np.sin(ph), -(Ax / (2 * np.pi * fe)) * np.cos(ph)

    spE, sxE = wave(t[None, :])
    spP, sxP = wave(t[None, :] - dt[:, None])
    return Fp_r[:, None] * (spE - spP) + Fx_r[:, None] * (sxE - sxP)


def psr_lnL(spine, data, src, L_m):
    r = residuals_chirp(spine, src, L_m=L_m)
    d = (data - r) / spine["sigma"][:, None]
    return -0.5 * np.sum(d * d, axis=1)          # (npsr,)  -- exactly separable


# ------------------------------------------------------------------
def fringe_posterior(spine, data, src, L_ref, dLf, guard=True):
    """Per-pulsar posterior over the NFRINGE ladder, anchored on L_ref (the PRIOR
    MEAN).  Each cell is scored at its PROFILED sub-fringe offset (the artefact
    guard); with guard=False the cell is scored at its centre -- the banked
    fabrication mode, kept as a switch so the difference is measurable."""
    npsr = spine["npsr"]
    js = np.arange(-NFRINGE, NFRINGE + 1)
    lp = np.zeros((npsr, len(js)))
    for i, j in enumerate(js):
        Lc = L_ref + j * dLf
        if not guard:
            lp[:, i] = psr_lnL(spine, data, src, Lc)
            continue
        for p in range(npsr):                    # profile the sub-fringe offset
            def f1(u):
                LL = Lc.copy(); LL[p] = Lc[p] + u * dLf[p]
                return -psr_lnL(spine, data, src, LL)[p]
            r = minimize_scalar(f1, bounds=(-0.5, 0.5), method="bounded",
                                options=dict(xatol=1e-10))
            lp[p, i] = -float(r.fun)
    # Gaussian distance prior on the ladder
    Lg = L_ref[:, None] + js[None, :] * dLf[:, None]
    lp = lp - 0.5 * ((Lg - spine["L_m"][:, None]) / spine["Lerr_m"][:, None]) ** 2
    lp -= lp.max(axis=1, keepdims=True)
    w = np.exp(lp)
    return w / w.sum(axis=1, keepdims=True), js


def sigma_points(mu, sig, scale, kappa=0.0):
    """Unscented 2n+1 = 5 points for the 2-D belief (log10_mc, log10_f)."""
    n = len(mu)
    lam = 3.0 - n + kappa
    P = np.diag((np.asarray(sig) * scale) ** 2)
    Sq = np.linalg.cholesky((n + lam) * P) if scale > 0 else np.zeros((n, n))
    pts = [mu] + [mu + Sq[:, i] for i in range(n)] + [mu - Sq[:, i] for i in range(n)]
    wts = [lam / (n + lam)] + [1.0 / (2 * (n + lam))] * (2 * n)
    return np.array(pts), np.array(wts)


def main():
    spine = S.build_spine(npsr=16, T_yr=18.0, cadence_per_yr=20.0)
    data = residuals_chirp(spine, SRC)
    snr = np.sqrt(np.sum((data / spine["sigma"][:, None]) ** 2))
    shat = S.sky_unit(SRC["costheta"], SRC["phi"])
    f0 = 10.0 ** SRC["log10_f"]
    dLf = S.C_LIGHT / (f0 * (1.0 - spine["phat"] @ shat))
    L_ref = spine["L_m"].copy()                  # the PRIOR MEAN -- never the truth
    dt = (spine["L_m"] / S.C_LIGHT) * (1.0 - spine["phat"] @ shat)
    print(f"[T5] SNR={snr:.1f}; retardation dt = {dt.min()/S.YR:.0f}-{dt.max()/S.YR:.0f} yr; "
          f"fringe {dLf.min()/S.PC:.3f}-{dLf.max()/S.PC:.3f} pc")

    # ---- how much of the pulsar-term wander is mc?  (the S4.2.4 cross-check) --
    def wander(axis, d):
        s2 = dict(SRC); s2[axis] = SRC[axis] + d
        k1, k2 = kchirp(SRC["log10_mc"]), kchirp(s2["log10_mc"])
        f1_, f2_ = 10.0 ** SRC["log10_f"], 10.0 ** s2["log10_f"]
        p1 = chirp_phase(-dt, f1_, k1, 0.0)
        p2 = chirp_phase(-dt, f2_, k2, 0.0)
        return np.abs(p2 - p1)
    w_mc = wander("log10_mc", 1e-3) / 1e-3
    w_f = wander("log10_f", 1e-3) / 1e-3
    # S4.2.4's "95-100 % of the wander is mc" is a FISHER-WEIGHTED statement --
    # the gradient times the axis's POSTERIOR sigma, not the bare per-dex
    # gradient.  Quoted both ways so the comparison is unambiguous.
    share_dex = w_mc ** 2 / (w_mc ** 2 + w_f ** 2)
    gm, gf = w_mc * SIG_MC_DEX[0], w_f * SIG_F_DEX
    share = gm ** 2 / (gm ** 2 + gf ** 2)
    print(f"[T5] pulsar-term wander per dex: |dphi/dlog10mc| = "
          f"{w_mc.min():.3e}-{w_mc.max():.3e} rad, |dphi/dlog10f| = "
          f"{w_f.min():.3e}-{w_f.max():.3e} rad")
    print(f"[T5] mc share of the wander VARIANCE: per-dex "
          f"{100*share_dex.min():.2f}-{100*share_dex.max():.2f} %  ->  "
          f"FISHER-WEIGHTED (x posterior sigma) "
          f"{100*share.min():.2f}-{100*share.max():.2f} %   "
          f"(S4.2.4 banked: 95-100 %)")

    # ---- the point E-step -------------------------------------------------
    w_pt, js = fringe_posterior(spine, data, SRC, L_ref, dLf)
    w_pt_nog, _ = fringe_posterior(spine, data, SRC, L_ref, dLf, guard=False)
    tv_guard = 0.5 * np.abs(w_pt - w_pt_nog).sum(axis=1)
    print(f"[T5] sub-fringe artefact guard: TV(profiled vs centre-pinned) = "
          f"{tv_guard.min():.3e}-{tv_guard.max():.3e} (median {np.median(tv_guard):.3e})")

    mu = np.array([SRC["log10_mc"], SRC["log10_f"]])

    def sp_estep(sig, scale):
        pts, wts = sigma_points(mu, sig, scale)
        acc = 0.0
        for (mc, lf), om in zip(pts, wts):
            s2 = dict(SRC); s2["log10_mc"] = mc; s2["log10_f"] = lf
            w, _ = fringe_posterior(spine, data, s2, L_ref, dLf)
            acc = acc + om * w
        return acc / acc.sum(axis=1, keepdims=True)

    # ---- GATE: width -> 0 reproduces the point E-step ----------------------
    print("\n[T5-GATE] sigma-point width -> 0 vs the point E-step")
    sig0 = np.array([SIG_MC_DEX[0], SIG_F_DEX])
    gate = []
    for sc in [1e-0, 1e-2, 1e-4, 1e-6, 1e-8, 0.0]:
        w = sp_estep(sig0, sc)
        tv = 0.5 * np.abs(w - w_pt).sum(axis=1)
        mx = np.abs(w - w_pt).max()
        gate.append(dict(scale=sc, tv_max=float(tv.max()), absmax=float(mx)))
        print(f"   scale={sc:8.1e}  max TV = {tv.max():.6e}   max|dw| = {mx:.6e}")

    # ---- READOUT: the shift at realistic-class widths ----------------------
    print("\n[T5-READOUT] fringe-posterior shift at S4.2.3/S4.2.4 widths")
    read = []
    jt = np.argmin(np.abs(js))          # the truth fringe (L_ref == injected L)
    for smc in SIG_MC_DEX:
        w = sp_estep(np.array([smc, SIG_F_DEX]), 1.0)
        tv = 0.5 * np.abs(w - w_pt).sum(axis=1)
        map_pt = js[np.argmax(w_pt, axis=1)]
        map_sp = js[np.argmax(w, axis=1)]
        ptrue_pt, ptrue_sp = w_pt[:, jt], w[:, jt]
        nflip = int(np.sum(map_pt != map_sp))
        read.append(dict(sig_mc=smc, tv=tv.tolist(),
                         ptrue_point=ptrue_pt.tolist(),
                         ptrue_sigma=ptrue_sp.tolist(),
                         map_point=map_pt.tolist(), map_sigma=map_sp.tolist()))
        print(f"   sigma(log10 Mc) = {smc:.3f} dex:  TV median {np.median(tv):.4f} "
              f"(min {tv.min():.4f}, max {tv.max():.4f});  "
              f"MAP-fringe flips {nflip}/{spine['npsr']};  "
              f"p(true fringe) {ptrue_pt.mean():.4f} -> {ptrue_sp.mean():.4f} "
              f"(mean over pulsars)")
        # per-pulsar wander in rad at this width, for the report's mechanism line
        print(f"      pulsar-term wander at this width: "
              f"{(w_mc*smc).min():.2e}-{(w_mc*smc).max():.2e} rad "
              f"= {(w_mc*smc/(2*np.pi)).min():.2e}-{(w_mc*smc/(2*np.pi)).max():.2e} fringes")

    np.savez(OUT, gate=json.dumps(gate), read=json.dumps(read),
             names=np.array(spine["names"]), js=js, w_point=w_pt,
             w_point_nogurad=w_pt_nog, tv_guard=tv_guard,
             wander_mc=w_mc, wander_f=w_f, mc_share=share, mc_share_perdex=share_dex,
             snr=snr, fringe_pc=dLf / S.PC, dt_yr=dt / S.YR,
             src=json.dumps(SRC), zero_noise=True)
    print(f"\n[T5] banked -> {OUT}")
    print("[T5] HANDOFF: calibration-against-truth over NOISY realisations is "
          "DEFERRED to SIEVE-A on banked data (cronus is RECOMPUTED-UNSAFE).")


if __name__ == "__main__":
    main()
