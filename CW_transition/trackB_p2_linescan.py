"""Track B — P2 STEP 1a: static soft-objective line scan through truth.

The "why continuous fails" proof for the paper. Along loud0's sky params (gwphi,
cos_gwtheta) through truth, with the OTHER distances at their SELF-CONSISTENT MAP (the
real algorithm's conditioning), evaluate:
  HARD_marg = sum_p max_n   [lnL_p(n) - offs^2/2sig^2]   (true-fringe cusp; peaks at truth)
  SOFT_marg = sum_p logsumexp_n[same]                     (the smoothed objective the soft
                                                           M-step climbs)
and locate each argmax. Claim: the SOFT argmax is STRUCTURALLY DISPLACED ~0.01-0.05 scaled
from the cusp -> a continuous optimiser of the smoothed objective lands off the needle.

Run (jug-gpu, background): python trackB_p2_linescan.py
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
from trackB_p1_map import batched_scan

CWT = "/home/mattm/projects/HSYMT/CW_transition"
PASSES = 3
SCALED = np.array([-0.08, -0.06, -0.045, -0.03, -0.02, -0.012, -0.006, -0.003,
                   0.0, 0.003, 0.006, 0.012, 0.02, 0.03, 0.045, 0.06, 0.08])


def _lse(x):
    m = np.max(x); return m + np.log(np.sum(np.exp(x - m)))


def surfaces_selfconsistent(P, pidx, val, prior):
    """Set loud0 source param pidx=val, self-consistent E-step (others->MAP), return
    HARD/SOFT marginal + sum-QTRUE (weight on true fringe k=0, diagnostic)."""
    ndist = P["n_dist"]; tt = P["theta_truth"]
    src = tt[ndist:].copy(); src[pidx] = val
    dist = tt[:ndist].copy()
    post = None
    for _ in range(PASSES):
        th = np.concatenate([dist, src])
        LNL, _ = batched_scan(P, [th])
        post = TE.fringe_posterior(P, LNL[0], None, prior, 1.0)
        newd = post["map_evalL"].copy()
        if np.max(np.abs(newd - dist)) < 1e-9:
            dist = newd; break
        dist = newd
    hard = soft = qtrue = 0.0
    for p in range(P["npsr"]):
        uk, offs_u, lnL_u, w = post["qlist"][p]
        sig = prior[p]
        lp = lnL_u - offs_u ** 2 / (2.0 * sig ** 2)
        hard += float(np.max(lp)); soft += float(_lse(lp))
        qtrue += float(w[uk == 0].sum()) if (uk == 0).any() else 0.0
    return hard, soft, qtrue


def run_axis(P, axis_name, pidx, prior):
    scale = TE.phi_scale(P)[pidx]
    v0 = P["theta_truth"][P["n_dist"] + pidx]
    H = np.zeros(len(SCALED)); S = np.zeros(len(SCALED)); Q = np.zeros(len(SCALED))
    t0 = time.time()
    for i, sc in enumerate(SCALED):
        H[i], S[i], Q[i] = surfaces_selfconsistent(P, pidx, v0 + sc * scale, prior)
    H -= H.max(); S -= S.max()
    iH = int(np.argmax(H)); iS = int(np.argmax(S))
    print(f"\n=== {axis_name} (pidx={pidx}, {time.time()-t0:.0f}s) ===", flush=True)
    print(f"  {'scaled':>8s} {'HARD-max':>10s} {'SOFT-max':>10s} {'QTRUE':>8s}", flush=True)
    for i, sc in enumerate(SCALED):
        mk = "  <-HARD*" if i == iH else ("  <-SOFT*" if i == iS else "")
        print(f"  {sc:8.4f} {H[i]:10.2f} {S[i]:10.2f} {Q[i]:8.2f}{mk}", flush=True)
    print(f"  HARD argmax at scaled {SCALED[iH]:+.4f}; SOFT argmax at scaled {SCALED[iS]:+.4f} "
          f"-> soft displaced {abs(SCALED[iS]-SCALED[iH]):.4f} scaled from the cusp", flush=True)
    return dict(scaled=SCALED, H=H, S=S, Q=Q, hard_argmax=SCALED[iH], soft_argmax=SCALED[iS])


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop"); prior = P["priors"]["lit"]
    print(f"  build_problem {time.time()-t0:.0f}s", flush=True)
    res = {}
    res["gwphi"] = run_axis(P, "gwphi", 1, prior)
    res["cos_gwtheta"] = run_axis(P, "cos_gwtheta", 0, prior)
    np.savez(f"{CWT}/trackB_p2_linescan.npz",
             scaled=SCALED,
             **{f"{a}_{k}": res[a][k] for a in res for k in ("H", "S", "Q")},
             gwphi_hardmax=res["gwphi"]["hard_argmax"], gwphi_softmax=res["gwphi"]["soft_argmax"],
             cos_hardmax=res["cos_gwtheta"]["hard_argmax"], cos_softmax=res["cos_gwtheta"]["soft_argmax"])
    print(f"\n  saved trackB_p2_linescan.npz ; total {time.time()-t0:.0f}s", flush=True)
