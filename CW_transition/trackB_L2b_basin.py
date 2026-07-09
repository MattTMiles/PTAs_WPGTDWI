"""Track B — LAMBDA probe, L2 supplement: CONDITIONAL RE-SOLVE PULL-IN RADIUS.

The RTK claim is: with the fringe integers FIXED at truth, the source likelihood becomes a
smooth quadratic whose max is truth -> Newton locks to mas. This tests that DIRECTLY,
decoupled from the (diverged) float: fix ALL distances at truth (n_true=0), perturb the loud
sources by delta (scaled), LBFGS-ascend the FULL lnL over source, report how close it locks
and whether lnL(result) exceeds lnL(truth). Two variants:
  (i)  optimise LOUD source params only (faint frozen at truth) -- clean basin test
  (ii) optimise ALL source params (faint free) -- as the pipeline re-solve does
Establishes the pull-in radius = the largest delta from which the re-solve returns to truth.
Run: python trackB_L2b_basin.py
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
import jax.numpy as jnp
import scipy.optimize as sopt
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
from trackB_p1_map import ang_deg
CWT = "/home/mattm/projects/HSYMT/CW_transition"


def main():
    print("=== L2b CONDITIONAL RE-SOLVE PULL-IN RADIUS (integers fixed at truth) ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False)
    fl = P["amo_full"]["fisher_logL"]; data = P["data_obs"]
    ndist = P["n_dist"]; theta_t = P["theta_truth"].copy()
    base_L = np.asarray(theta_t[:ndist]).copy()               # distances FIXED at truth (n=0)
    baseLj = jnp.asarray(base_L)
    src_t = np.asarray(theta_t[ndist:]).copy()
    scale = TE.phi_scale(P); plo, phi = TE.phi_bounds(P)
    n_loud = TE.CONFIGS["pop"]["population"][0]; loud_sl = slice(0, 8 * n_loud)
    lnL_truth = float(fl(jnp.asarray(theta_t), data))
    print(f"  lnL(truth) = {lnL_truth:.3f}; built {time.time()-t0:.0f}s", flush=True)

    # value_and_grad wrt full source; for LOUD-only we zero the faint grads + keep faint fixed
    vg_full = jax.jit(lambda s: jax.value_and_grad(
        lambda ss: -fl(jnp.concatenate([baseLj, ss]), data))(s))

    def resolve(src0, loud_only):
        lo = np.asarray(plo, float).copy(); hi = np.asarray(phi, float).copy()
        if loud_only:                                          # pin faint at truth via tight bounds
            lo[8 * n_loud:] = src_t[8 * n_loud:] - 1e-12; hi[8 * n_loud:] = src_t[8 * n_loud:] + 1e-12
        def fun(s):
            v, g = vg_full(jnp.asarray(s)); return float(v), np.asarray(g, float)
        r = sopt.minimize(fun, np.asarray(src0, float), jac=True, method="L-BFGS-B",
                          bounds=list(zip(lo, hi)), options=dict(maxiter=400, ftol=1e-13, gtol=1e-10))
        s = r.x
        skyd = max(ang_deg(s[8*k], s[8*k+1], src_t[8*k], src_t[8*k+1]) for k in range(n_loud))
        off = np.abs(s[loud_sl] - src_t[loud_sl]) / scale[loud_sl]
        lnl = -float(r.fun)
        return skyd, float(off.max()), lnl

    rng = np.random.default_rng(3)
    deltas = [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
    print(f"\n  {'delta':>7s} | {'LOUD-only sky/off/dlnL':>34s} | {'ALL-src sky/off/dlnL':>34s}", flush=True)
    rows = []
    for d in deltas:
        pert = rng.normal(size=8 * n_loud) * scale[loud_sl] * d
        src0 = src_t.copy(); src0[loud_sl] = src_t[loud_sl] + pert
        sk_l, of_l, ll_l = resolve(src0, True)
        sk_a, of_a, ll_a = resolve(src0, False)
        rows.append((d, sk_l, of_l, ll_l, sk_a, of_a, ll_a))
        print(f"  {d:7.0e} | sky={sk_l:8.4f} off={of_l:7.4f} dlnL={ll_l-lnL_truth:+9.2f} | "
              f"sky={sk_a:8.4f} off={of_a:7.4f} dlnL={ll_a-lnL_truth:+9.2f}  ({time.time()-t0:.0f}s)", flush=True)

    locks_l = [d for d, sk_l, *_ in rows if sk_l < 0.006]
    locks_a = [d for d, *_ , sk_a, of_a, ll_a in [(r[0], r[1], r[2], r[3], r[4], r[5], r[6]) for r in rows] if sk_a < 0.006]
    pull_l = max(locks_l) if locks_l else 0.0
    pull_a = max([r[0] for r in rows if r[4] < 0.006], default=0.0)
    print(f"\n  PULL-IN RADIUS (sky<0.006deg): loud-only <= {pull_l:.0e} ; all-src <= {pull_a:.0e}", flush=True)
    print(f"  cert tolerance ~1e-4 (scaled); achievable float precision ~O(1) scaled", flush=True)
    np.savez(f"{CWT}/trackB_L2b_basin.npz", deltas=np.array(deltas),
             rows=np.array([[r[0], r[1], r[2], r[3], r[4], r[5], r[6]] for r in rows]),
             lnL_truth=lnL_truth, pull_loud=pull_l, pull_all=pull_a)
    print(f"  saved trackB_L2b_basin.npz ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}", flush=True)
    sys.exit(main())
