"""Track B — B1 STEP 1: the (f, mc) registration ladder + wide-lane cascade feasibility.

F2 asked, for the SKY axis: does the per-(pulsar, source) registration-tolerance ladder span
continuously from the blind float down to the certification tolerance? Answer was NO -- zero
rungs above the float, so the GPS-RTK wide-lane cascade could never fix its first rung.

B1's pilot showed the TARGETED scenario (sky supplied by an EM counterpart) leaves TWO
registration axes, not one: log10_fgw and log10_mc. F2 was structurally blind to the second
(`trackB_F2_ladder.py:47` uses the chirp-free phase 2 pi f L (1-cos mu)/c, whose gradient in
log10_mc is identically zero). This module asks F2's question again on the joint (f, mc) axes,
where the pilot found the loosest mc rung at 60 SCALED -- far ABOVE the Earth-term float,
unlike sky. So a first rung may exist here.

The two halves
--------------
A. LADDER (pure geometry, exact-phase autodiff; no likelihood, no GPU).
   Pulsar a's pulsar-term phase wanders by
       dPhi_a = sum_k [ g_f,ak * d(log10_fgw)_k + g_mc,ak * d(log10_mc)_k ]     (scaled units)
   over the current (f, mc) box. Its fringe is registrable when the rms wander < 1 rad. So
       R_a = 1 / sqrt( sum_k [ (sig_f,k * g_f,ak)^2 + (sig_mc,k * g_mc,ak)^2 ] )
   is the factor by which the box must SHRINK before pulsar a registers. R_a >= 1 means a
   registers at the current box: the cascade has a first rung. The ladder is the sorted R_a
   spectrum; it must reach the census-3 (the pulsars that actually certify).

B. LOOP GAIN (likelihood; the half F2 never got to run).
   Fixing a subset S of fringes and re-solving (f, mc) conditionally is EXACTLY a masked
   likelihood: pulsar terms ON for S, OFF elsewhere. `trackB_b1_core.MaskedDelay` makes the
   mask a runtime argument, so the conditional Fisher

       F_S = -H[ lnL(theta ; data_S, mask_S) ]   over the 6 params (log10_fgw, log10_mc) x 3

   is one call per subset with NO optimiser and NO pull-in problem. Asimov data is injected
   with the SAME mask, so residual(truth) = 0 and H = -F exactly (the D2 identity).
   Cascade converges iff sigma_S (from F_S + prior) is smaller than the registration tolerance
   of the NEXT rung down. That ratio IS the loop gain, per rung.

Everything here is TRUTH-PLACED Asimov geometry (labelled). Env: jug-gpu, autotune off, x64.
Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_b1_core as C
import trackB_b1_pilot as PL
import trackB_estimator as TE

CWT = "/home/mattm/projects/HSYMT/CW_transition"
I_MC, I_FGW = PL.I_MC, PL.I_FGW
# certification tolerances, SCALED, measured by the pilot's M2 (last delta keeping >=2 census
# pulsars at P(true) > 0.9). Quoted, not re-derived.
CERT_TOL = {"log10_fgw": 3e-5, "log10_mc": 1e-3}


# ============================================================
# A. the geometry ladder
# ============================================================
def build_gradients(P, t_rel):
    """g_f[a,k], g_mc[a,k] = d(pulsar-term phase)/d(scaled param), radians. Exact phase."""
    src = P["theta_truth"][P["n_dist"]:]
    scale = TE.phi_scale(P)
    npsr = P["npsr"]
    g_f = np.zeros((npsr, C.N_LOUD)); g_mc = np.zeros((npsr, C.N_LOUD))
    for k in range(C.N_LOUD):
        p4 = jnp.array([src[8*k+PL.I_COSGT], src[8*k+PL.I_GWPHI], src[8*k+I_MC], src[8*k+I_FGW]])
        for a in range(npsr):
            g = np.asarray(PL._grad_dphi(p4, jnp.asarray(P["disco_psrs"][a].pos),
                                         float(P["L0"][a]), t_rel))
            g_mc[a, k] = g[2] * scale[I_MC]
            g_f[a, k] = g[3] * scale[I_FGW]
    C._finite("phase gradients", g_f, g_mc)
    return g_f, g_mc


def registration_radius(g_f, g_mc, sig_f, sig_mc):
    """R_a = box-shrink factor at which pulsar a's pulsar-term phase wander falls to 1 rad."""
    var = np.sum((sig_f[None, :] * g_f) ** 2 + (sig_mc[None, :] * g_mc) ** 2, axis=1)
    return 1.0 / np.sqrt(np.clip(var, 1e-300, None))


def mode_ladder(args):
    print("=== B1 STEP 1A: the joint (f, mc) registration ladder ===", flush=True)
    P = C.build_b1_problem()
    tref = 86400.0 * 51544.5
    t_rel = float(np.mean(np.concatenate([np.asarray(p.toas) for p in P["disco_psrs"]])) - tref)
    names = P["names"]

    m3 = np.load(f"{CWT}/b1_pilot_m3.npz")
    sig = m3["sigma_scaled"]                       # Earth-term posterior, SCALED
    sig_f = np.array([sig[8 * k + I_FGW] for k in range(C.N_LOUD)])
    sig_mc = np.array([sig[8 * k + I_MC] for k in range(C.N_LOUD)])
    print(f"  FLOAT (Earth-term posterior, scaled): sigma_f = {np.round(sig_f,5)}", flush=True)
    print(f"                                        sigma_mc = {np.round(sig_mc,4)} "
          f"(= the prior; the Earth term carries no mc information)", flush=True)

    g_f, g_mc = build_gradients(P, t_rel)
    R = registration_radius(g_f, g_mc, sig_f, sig_mc)
    order = np.argsort(-R)                         # loosest (largest R) first

    # what fraction of the phase-wander variance each axis contributes, per pulsar
    var_f = np.sum((sig_f[None, :] * g_f) ** 2, axis=1)
    var_mc = np.sum((sig_mc[None, :] * g_mc) ** 2, axis=1)
    frac_mc = var_mc / np.clip(var_f + var_mc, 1e-300, None)

    cen = np.array([names.index(n) for n in C.CENSUS])
    anc = P["anchor_idx"]
    print(f"\n  registration radius R_a (box-shrink factor to reach 1 rad; R>=1 = registers NOW):")
    print(f"    max {R.max():.3e} ({names[order[0]]})   median {np.median(R):.3e}   "
          f"min {R.min():.3e} ({names[order[-1]]})", flush=True)
    n_ign = int(np.sum(R >= 1.0))
    print(f"    pulsars with R >= 1 (IGNITION set)      : {n_ign}", flush=True)
    for thr in (0.3, 0.1, 0.03, 0.01):
        print(f"    pulsars with R >= {thr:<5}                 : {int(np.sum(R>=thr))}", flush=True)
    print(f"\n  census-3 / anchor:", flush=True)
    for a in list(cen) + [anc]:
        print(f"    {names[a]:14s} R = {R[a]:.3e}  (needs box shrink {1/R[a]:.3e}x)   "
              f"mc share of wander {frac_mc[a]*100:5.1f}%", flush=True)
    print(f"\n  top-10 loosest rungs:", flush=True)
    print(f"    {'pulsar':14s} {'R_a':>10s} {'L (kpc)':>9s} {'mc share':>9s}", flush=True)
    for a in order[:10]:
        print(f"    {names[a]:14s} {R[a]:10.3e} {P['L0'][a]:9.3f} {frac_mc[a]*100:8.1f}%", flush=True)

    # ---- ladder-span test, F2's criterion, on the joint axes ----
    S_target = 1.0 / R[cen].min()                  # shrink needed for the LAST census pulsar
    print(f"\n  TARGET: census-3 all register at box shrink {S_target:.3e}x "
          f"(binding: {names[cen[np.argmin(R[cen])]]})", flush=True)
    band = np.sort(R[(R <= 1.0) & (R >= 1.0 / S_target)])[::-1]
    ignite = R.max() >= 1.0
    print(f"  IGNITION (a rung at the float, R >= 1): {'YES' if ignite else 'NO'} "
          f"-- {n_ign} pulsars", flush=True)
    if len(band) > 1:
        lr = np.log10(band)
        gaps = -np.diff(lr)
        spans = gaps.max() < 0.7
        print(f"  in-band rungs between the float and the census target: {len(band)}; "
              f"max log10 gap {gaps.max():.3f} dex at R ~ {10**lr[np.argmax(gaps)]:.2e}", flush=True)
        print(f"  LADDER SPANS (no gap > 0.7 dex, F2's criterion): {spans}", flush=True)
    else:
        spans = False
        print(f"  in-band rungs < 2 -> ladder does NOT span", flush=True)

    np.savez(f"{CWT}/b1_ladder.npz", R=R, g_f=g_f, g_mc=g_mc, sig_f=sig_f, sig_mc=sig_mc,
             frac_mc=frac_mc, names=np.array(names), order=order, S_target=S_target,
             ignite=ignite, spans=spans, L0=P["L0"])
    print(f"\n  saved b1_ladder.npz", flush=True)
    return P, R, order, sig_f, sig_mc, g_f, g_mc


# ============================================================
# B. loop gain: conditional (f, mc) precision given a FIXED fringe subset
# ============================================================
def make_cond_fisher(P, sel):
    """x -> Hessian of lnL(theta[sel]=x ; data, mask) over the len(sel) selected theta slots.
    data and mask are RUNTIME args, so ONE compile serves every subset. HVP assembly (D1)."""
    tt = jnp.asarray(P["theta_truth"])
    sel_j = jnp.asarray(sel)
    logL = P["amo"]["logL"]

    def f_of_x(x, data, mask):
        return logL(tt.at[sel_j].set(x), data, mask)

    g = jax.grad(f_of_x, argnums=0)

    def hvp(x, v, data, mask):
        return jax.jvp(lambda z: g(z, data, mask), (x,), (v,))[1]

    return jax.jit(jax.vmap(hvp, in_axes=(None, 0, None, None)))


def mode_loopgain(args):
    P, R, order, sig_f, sig_mc, g_f, g_mc = mode_ladder(args)
    print("\n\n=== B1 STEP 1B: cascade loop gain (conditional (f,mc) Fisher at fixed fringes) ===",
          flush=True)
    nd = P["n_dist"]; names = P["names"]
    sel = np.array([nd + 8 * k + j for k in range(C.N_LOUD) for j in (I_FGW, I_MC)])
    lab = [f"loud{k}_{'f' if j==I_FGW else 'mc'}" for k in range(C.N_LOUD) for j in (I_FGW, I_MC)]
    scale = TE.phi_scale(P)
    sc = np.array([scale[j] for k in range(C.N_LOUD) for j in (I_FGW, I_MC)])

    lo, hi = TE.phi_bounds(P)
    var_prior = np.array([((hi[8*k+j] - lo[8*k+j]) / scale[j]) ** 2 / 12.0
                          for k in range(C.N_LOUD) for j in (I_FGW, I_MC)])
    Pi = np.diag(1.0 / var_prior)

    hvp = make_cond_fisher(P, sel)
    x0 = jnp.asarray(P["theta_truth"][sel])
    V = jnp.asarray(np.eye(len(sel)))

    # nested subsets, loosest rung first (the cascade order), plus reference sets
    subsets = [("earth only (mask=0)", np.array([], int))]
    for n in (1, 3, 6, 12, 24, 48, 116):
        subsets.append((f"top-{n} rungs", order[:n]))
    subsets.append(("census-3 only", np.array([names.index(n) for n in C.CENSUS])))

    print(f"\n  conditional 1-sigma on the 6 registration params, given the SUBSET's fringes "
          f"fixed at truth", flush=True)
    print(f"  (Asimov data injected with the SAME mask -> residual(truth)=0 -> H = -Fisher exactly)",
          flush=True)
    print(f"\n  {'subset':22s} {'sig_f (med)':>12s} {'sig_mc (med)':>13s} {'sig_f/tol_f':>12s} "
          f"{'sig_mc/tol_mc':>14s} {'next rung R':>12s} {'GAIN':>8s}", flush=True)
    rows = []
    t0 = time.time()
    for lab_s, idx in subsets:
        mask = np.zeros(P["npsr"]); mask[idx] = 1.0
        mj = jnp.asarray(mask)
        data = P["amo"]["inject_delay"](jnp.asarray(P["theta_truth"]), mj)
        H = np.asarray(hvp(x0, V, data, mj))
        C._finite(f"Hessian {lab_s}", H)
        F = -0.5 * (H + H.T)
        Fs = F * sc[:, None] * sc[None, :]                 # scaled units
        Cov = np.linalg.inv(Fs + Pi)
        s = np.sqrt(np.clip(np.diag(Cov), 0, None))
        s_f = np.median(s[0::2]); s_mc = np.median(s[1::2])
        # after this subset is fixed, what box do we have -> which pulsar registers next?
        Rn = registration_radius(g_f, g_mc, s[0::2], s[1::2])
        free = np.setdiff1d(np.arange(P["npsr"]), idx)
        nxt = Rn[free].max() if len(free) else np.nan
        gain = nxt                                          # >=1 -> another rung ignites
        rows.append((lab_s, s_f, s_mc, s_f / CERT_TOL["log10_fgw"], s_mc / CERT_TOL["log10_mc"],
                     nxt, int(np.sum(Rn[free] >= 1.0))))
        print(f"  {lab_s:22s} {s_f:12.3e} {s_mc:13.3e} {s_f/CERT_TOL['log10_fgw']:12.1f} "
              f"{s_mc/CERT_TOL['log10_mc']:14.1f} {nxt:12.3e} {gain:8.2f}", flush=True)
    print(f"\n  ({time.time()-t0:.0f}s)", flush=True)

    # ---- cascade simulation, seeded from EACH conditioning tier's float ----
    allhist = {}
    for tier in ("A", "B", "C"):
        sf, smc, prov = tier_sigmas(P, tier)
        print(f"\n  CASCADE from tier {tier} float ({prov})", flush=True)
        print(f"    sigma_f {np.round(sf,5)}  sigma_mc {np.round(smc,4)} (scaled)", flush=True)
        allhist[tier] = _cascade(P, hvp, x0, V, sc, Pi, g_f, g_mc, sf, smc, names)
    np.savez(f"{CWT}/b1_cascade_tiers.npz",
             **{f"hist_{t}": np.array(v[0]) for t, v in allhist.items()},
             **{f"fixed_{t}": v[1] for t, v in allhist.items()},
             **{f"conv_{t}": v[2] for t, v in allhist.items()})
    print(f"\n  saved b1_cascade_tiers.npz", flush=True)
    return 0


def tier_sigmas(P, tier, sig_p_over_p=1e-3):
    """1-sigma (f, mc) float widths, SCALED, for the three pre-registered conditioning tiers.
    Tier C's mc comes from log10 h = (5/3)log10 Mc + (2/3)log10 f - log10 D_L + const."""
    scale = TE.phi_scale(P)
    m3 = np.load(f"{CWT}/b1_pilot_m3.npz")
    sig = m3["sigma_scaled"]
    sig_f_e = np.array([sig[8 * k + I_FGW] for k in range(C.N_LOUD)])
    sig_h_e = np.array([sig[8 * k + 5] for k in range(C.N_LOUD)]) * scale[5]      # dex
    sig_mc_prior = np.array([sig[8 * k + I_MC] for k in range(C.N_LOUD)])
    if tier == "A":
        return sig_f_e, sig_mc_prior, "sky only: Earth-term f, prior mc"
    sig_f_dex = sig_p_over_p / np.log(10.0)
    sig_f = np.minimum(np.full(C.N_LOUD, sig_f_dex / scale[I_FGW]), sig_f_e)
    if tier == "B":
        return sig_f, sig_mc_prior, f"sky + EM period (sigma_P/P={sig_p_over_p:.0e})"
    sig_mc_dex = 0.6 * np.sqrt(sig_h_e ** 2 + ((2.0 / 3.0) * sig_f_dex) ** 2 + 0.005 ** 2)
    sig_mc = np.minimum(sig_mc_dex / scale[I_MC], sig_mc_prior)
    return sig_f, sig_mc, f"sky + period + D_L: sigma_mc {np.round(sig_mc_dex,4)} dex"


def _cascade(P, hvp, x0, V, sc, Pi, g_f, g_mc, s_f_cur, s_mc_cur, names, nmax=14):
    fixed = np.zeros(P["npsr"], bool)
    hist = []
    for it in range(nmax):
        Rn = registration_radius(g_f, g_mc, s_f_cur, s_mc_cur)
        newly = (Rn >= 1.0) & (~fixed)
        if not newly.any():
            print(f"    round {it}: no new pulsar reaches R>=1 (best free R = "
                  f"{Rn[~fixed].max() if (~fixed).any() else np.nan:.3e}) -> CASCADE STALLS", flush=True)
            break
        fixed |= newly
        mask = fixed.astype(float)
        mj = jnp.asarray(mask)
        data = P["amo"]["inject_delay"](jnp.asarray(P["theta_truth"]), mj)
        H = np.asarray(hvp(x0, V, data, mj))
        F = -0.5 * (H + H.T)
        Fs = F * sc[:, None] * sc[None, :]
        s = np.sqrt(np.clip(np.diag(np.linalg.inv(Fs + Pi)), 0, None))
        s_f_cur, s_mc_cur = s[0::2], s[1::2]
        ncen = int(np.sum(fixed[[names.index(n) for n in C.CENSUS]]))
        hist.append((it, int(fixed.sum()), float(np.median(s_f_cur)), float(np.median(s_mc_cur)), ncen))
        print(f"    round {it}: +{int(newly.sum())} registered (total {int(fixed.sum())}/116); "
              f"sigma_f {np.median(s_f_cur):.3e} sigma_mc {np.median(s_mc_cur):.3e}; "
              f"census-3 registered {ncen}/3", flush=True)
        if fixed.sum() == P["npsr"]:
            break
    conv = bool(np.all(fixed[[names.index(n) for n in C.CENSUS]]))
    print(f"    >>> reaches the census-3: {conv}; registered {int(fixed.sum())}/116; "
          f"final sigma_f {np.median(s_f_cur):.3e} (cert tol {CERT_TOL['log10_fgw']:.0e}), "
          f"sigma_mc {np.median(s_mc_cur):.3e} (cert tol {CERT_TOL['log10_mc']:.0e})", flush=True)
    return hist, fixed, conv


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="loopgain", choices=["ladder", "loopgain"])
    a = ap.parse_args()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    if a.mode == "ladder":
        mode_ladder(a); sys.exit(0)
    sys.exit(mode_loopgain(a))
