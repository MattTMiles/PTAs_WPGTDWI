"""Track B — B1 STEP 2C: the BREAK-EVEN RESPONSE CURVE Z_box(lambda_mc).

WHY A CURVE AND NOT A POINT. The Tier-C break-even was computed as
`lambda_mc = exp((d - logit 0.95)/3)`, i.e. assuming `Z_box ∝ lambda_mc^3` (the uniform-plateau-
density approximation R used for theta* = 0.188 deg). **Tier A's first SMC seed falsified that
assumption.** Tier A's 6-D box is e^{10.1} ~ 24000x larger in volume than Tier C's, yet
    ln Z_box(A) = 405632.432   vs   ln Z_box(C) = 405633.035
agree inside the 0.86-nat sampler scatter. Z_box is an INTEGRAL, so Z_box(A) >= Z_box(C)
necessarily; measuring them equal means the plateau's evidence has SATURATED -- it is confined well
inside Tier C's box, and enlarging the box adds volume carrying negligible likelihood. So Z_box is
insensitive to the box over precisely the range the old break-even extrapolated across.

Measuring a single shrunken box would repeat the error class (one point + an assumed scaling).
Instead: measure the RESPONSE CURVE and read the break-even off it.

  lambda_mc in {1.0 (= Tier C, BANKED), 0.3, 0.12, 0.05}, 2 seeds each.
  The mc half-widths scale by lambda_mc; the f half-widths stay at Tier C's (period-supplied).

DELIVERABLES
  (a) SATURATION SCALE: the lambda at which Z_box starts responding = the plateau's own mc extent.
      A NEW measured quantity -- it is the width of the wrong-fringe plateau in chirp mass.
  (b) TRUE BREAK-EVEN lambda*: where f = Z_needle/(Z_needle + Z_box(lambda)) = 0.95, read off the
      curve. If lambda* lies below the smallest measured lambda, report a BOUND, not a value.
  (c) CORRECTED DEFICIT = 1 / lambda*, replacing the suspended 8.29x. No proportionality anywhere.

CAVEAT ENFORCED IN CODE: as the box shrinks toward the needle, the needle sub-box occupies a growing
fraction of it. Particles within 5 needle-sigma are excised (as R2 did) and the excised count is
reported per lambda; once excision is non-negligible, Z_box is no longer a clean "plateau" evidence
and the point is FLAGGED rather than silently used.

Env: jug-gpu, autotune off, x64. Matt commits; never commit from here.
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
import trackB_b1_referendum as REF
import trackB_estimator as TE

CWT = "/home/mattm/projects/HSYMT/CW_transition"
I_MC, I_FGW = 3, 4
LOGIT95 = np.log(0.95 / 0.05)

# Tier C, banked (4 seeds, frozen protocol)
LNZN = 405629.6337                 # ln Z_needle, tier-independent to 0.003 nat (A vs C)
LNZBOX_C = 405633.035              # lambda_mc = 1.0
LNZBOX_C_SEM = 0.429


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lams", default="0.3,0.12,0.05")
    ap.add_argument("--N", type=int, default=192)
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--mcmc", type=int, default=3)
    ap.add_argument("--max_mcmc", type=int, default=20)
    ap.add_argument("--target_acc", type=float, default=0.25)
    ap.add_argument("--adapt_acc", type=float, default=0.35)
    a = ap.parse_args()
    if a.adapt_acc <= a.target_acc:
        sys.exit("CONFIG ERROR: --adapt_acc must exceed --target_acc")
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    print("=== B1 STEP 2C: break-even RESPONSE CURVE Z_box(lambda_mc) ===", flush=True)

    P = C.build_b1_problem()
    tt = P["theta_truth"]; mask1 = P["mask_one"]
    data = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))
    G = REF.TargetedMarg(P, data, mask1)
    G.E.gate_fast_full(G.src_from_x(np.zeros((1, G.D))))

    half_f, half_mc, prov = REF.tier_boxes(P, "C")
    print(f"  Tier-C reference box: {prov}", flush=True)
    scale_mc = TE.phi_scale(P)[I_MC]

    # needle sigma for excision, from the banked Tier-C needle diagnostics if present
    try:
        z = np.load(f"{CWT}/b1_referendum_tierC.npz", allow_pickle=True)
        needle_sig = np.asarray(z["sig_needle"])
    except Exception:
        needle_sig = np.full(G.D, 1e-5)

    lams = [float(x) for x in a.lams.split(",")]
    rows = [dict(lam=1.0, lnZ=LNZBOX_C, sem=LNZBOX_C_SEM, seeds=4, excised=0, banked=True,
                 mc_dex=float(np.median(half_mc)) * scale_mc)]
    for lam in lams:
        bh = np.empty(G.D)
        for k in range(C.N_LOUD):
            bh[REF.active_index(k, I_FGW)] = half_f[k]
            bh[REF.active_index(k, I_MC)] = half_mc[k] * lam
        lnV = float(np.sum(np.log(2 * bh)))
        print(f"\n  --- lambda_mc = {lam:g}: mc half-width {np.round(half_mc*lam,5)} scaled "
              f"= {np.round(half_mc*lam*scale_mc,6)} dex; ln V_box = {lnV:.4f} ---", flush=True)
        lz, exc, accs = [], [], []
        for sd in range(a.seeds):
            t0 = time.time()
            r = REF.z_box(G, bh, a.N, sd, needle_sig, n_mcmc=a.mcmc, tag=f"lam{lam:g}",
                          target_acc=a.target_acc, max_mcmc=a.max_mcmc, adapt_acc=a.adapt_acc)
            lz.append(r["logZ"]); exc.append(r["excised"]); accs.append(r["acc_hi"])
            print(f"     [seed {sd}] ln Z_box = {r['logZ']:.3f} (density {r['logZ_density']:.3f} "
                  f"+ lnV {r['lnVbox']:.3f}); excised {r['excised']}/{a.N} acc_hi {r['acc_hi']:.3f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        lz = np.array(lz)
        sem = float(lz.std(ddof=1) / np.sqrt(len(lz))) if len(lz) > 1 else np.inf
        exc_frac = float(np.mean(exc)) / a.N
        rows.append(dict(lam=lam, lnZ=float(lz.mean()), sem=sem, seeds=len(lz),
                         excised=exc_frac, banked=False,
                         mc_dex=float(np.median(half_mc)) * lam * scale_mc,
                         acc_hi=float(min(accs))))
        print(f"     ln Z_box(lam={lam:g}) = {lz.mean():.3f} +- {sem:.3f}; "
              f"needle excision {exc_frac*100:.1f}% of particles"
              f"{'  *** FLAG: excision non-negligible, not a clean plateau evidence ***' if exc_frac > 0.02 else ''}",
              flush=True)

    # ---- read the curve ----
    rows.sort(key=lambda r: -r["lam"])
    print(f"\n=== RESPONSE CURVE ===", flush=True)
    print(f"  {'lambda':>8s} {'mc box (dex)':>13s} {'ln Z_box':>12s} {'+-sem':>7s} "
          f"{'d=lnZn-lnZbox':>14s} {'f':>10s} {'excised':>8s}", flush=True)
    for r in rows:
        d = LNZN - r["lnZ"]; f = 1.0 / (1.0 + np.exp(-d))
        r["d"] = d; r["f"] = f
        print(f"  {r['lam']:8.3g} {r['mc_dex']:13.5f} {r['lnZ']:12.3f} {r['sem']:7.3f} "
              f"{d:+14.3f} {f:10.5g} {r['excised']*100:7.1f}%", flush=True)

    # (a) saturation scale: the largest lambda at which ln Z_box has dropped by > 2 sem from lam=1
    ref = rows[0]["lnZ"]
    sat = None
    for r in rows[1:]:
        thresh = 2.0 * max(r["sem"], rows[0]["sem"])
        if (ref - r["lnZ"]) > thresh:
            sat = r["lam"]
            break
    if sat is None:
        print(f"\n  (a) SATURATION SCALE: Z_box has NOT started responding down to lambda = "
              f"{rows[-1]['lam']:g} (mc box {rows[-1]['mc_dex']:.5f} dex). The plateau's mc extent "
              f"is SMALLER than the smallest box measured -> report as an UPPER BOUND.", flush=True)
    else:
        print(f"\n  (a) SATURATION SCALE: Z_box starts responding at lambda ~ {sat:g} "
              f"(mc box {[r['mc_dex'] for r in rows if r['lam']==sat][0]:.5f} dex). "
              f"That is the PLATEAU'S OWN mc EXTENT -- a newly measured quantity.", flush=True)

    # (b) true break-even: f = 0.95  <=>  ln Z_box = LNZN - logit(0.95)
    target = LNZN - LOGIT95
    lams_arr = np.array([r["lam"] for r in rows]); lnz = np.array([r["lnZ"] for r in rows])
    lam_star, bound = None, False
    for i in range(len(rows) - 1):
        y0, y1 = lnz[i], lnz[i + 1]
        if (y0 - target) * (y1 - target) <= 0 and y0 != y1:
            x0, x1 = np.log(lams_arr[i]), np.log(lams_arr[i + 1])
            lam_star = float(np.exp(x0 + (target - y0) * (x1 - x0) / (y1 - y0)))
            break
    if lam_star is None:
        lam_star = float(lams_arr.min()); bound = True
    if bound:
        print(f"  (b) TRUE BREAK-EVEN: ln Z_box never falls to {target:.3f} within the measured "
              f"range -> lambda* < {lam_star:g}. Report as a BOUND, not a value.", flush=True)
        print(f"  (c) CORRECTED DEFICIT: > {1.0/lam_star:.3g}x  (BOUND)", flush=True)
    else:
        print(f"  (b) TRUE BREAK-EVEN: lambda* = {lam_star:.4g} "
              f"(mc box {np.median(half_mc)*lam_star*scale_mc:.6f} dex), by interpolation in log "
              f"lambda on the MEASURED curve -- no proportionality assumed.", flush=True)
        print(f"  (c) CORRECTED DEFICIT = 1/lambda* = {1.0/lam_star:.3g}x "
              f"(replaces the suspended 8.29x)", flush=True)

    np.savez(f"{CWT}/b1_breakeven_curve.npz",
             lam=lams_arr, lnZbox=lnz, sem=np.array([r["sem"] for r in rows]),
             f=np.array([r["f"] for r in rows]), d=np.array([r["d"] for r in rows]),
             excised=np.array([r["excised"] for r in rows]),
             mc_dex=np.array([r["mc_dex"] for r in rows]),
             lnZn=LNZN, lam_star=lam_star, is_bound=bound, deficit=1.0 / lam_star,
             saturation_lam=(sat if sat is not None else np.nan))
    print(f"\n  saved b1_breakeven_curve.npz", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
