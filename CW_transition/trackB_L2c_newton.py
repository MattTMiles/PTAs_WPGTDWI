"""Track B — LAMBDA probe, L2 supplement 2: NEWTON conditional re-solve.

The RTK claim is Newton-specific: integers fixed -> source objective is a smooth quadratic ->
ONE (or a few damped) Newton steps lock to truth regardless of conditioning. LBFGS (L2b)
stalled on the cusp shoulder, but that can be pure ill-conditioning (~1e11: sky sharp, faint
flat), NOT a refutation of the quadratic claim. This tests Newton DIRECTLY on the 24 loud
params (distances + 13 faint sources FIXED at truth), with backtracking + Levenberg damping.
  - delta=0: gradient norm + Hessian eigenvalue sign at EXACT truth (is truth a stationary max?)
  - delta>0: damped-Newton from truth+delta -> does it lock to sky<0.006deg / recover lnL(truth)?
If Newton locks from >~1e-3 but LBFGS did not -> the objective IS quadratic, failure was the
optimizer; the LAMBDA verdict then rests solely on the float being too far. If Newton ALSO
fails from ~1e-4 -> the cusp is genuinely non-quadratic and the smooth-quadratic premise is false.
Run: python trackB_L2c_newton.py
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
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_estimator as TE
from trackB_p1_map import ang_deg
CWT = "/home/mattm/projects/HSYMT/CW_transition"


def main():
    print("=== L2c NEWTON conditional re-solve (loud params, integers+faint fixed at truth) ===", flush=True)
    t0 = time.time()
    P = TE.build_problem("pop", build_earth=False)
    fl = P["amo_full"]["fisher_logL"]; data = P["data_obs"]
    ndist = P["n_dist"]; theta_t = np.asarray(P["theta_truth"]).copy()
    src_t = theta_t[ndist:].copy(); scale = TE.phi_scale(P)
    n_loud = TE.CONFIGS["pop"]["population"][0]; nlp = 8 * n_loud
    lnL_truth = float(fl(jnp.asarray(theta_t), data))
    fixed_tail = jnp.asarray(theta_t[ndist + nlp:])       # faint source params, fixed at truth
    base_dist = jnp.asarray(theta_t[:ndist])              # distances fixed at truth (n=0)

    # lnL as a function of the 24 loud params only
    def lnl_loud(x):
        theta = jnp.concatenate([base_dist, x, fixed_tail])
        return fl(theta, data)
    g_fn = jax.jit(jax.grad(lnl_loud))
    H_fn = jax.jit(jax.hessian(lnl_loud))

    x_t = jnp.asarray(src_t[:nlp])
    # delta=0 diagnostics: gradient + Hessian spectrum at exact truth
    g0 = np.asarray(g_fn(x_t)); H0 = np.asarray(H_fn(x_t))
    ev = np.linalg.eigvalsh(H0)
    print(f"  lnL(truth)={lnL_truth:.3f}; built {time.time()-t0:.0f}s", flush=True)
    print(f"  AT TRUTH: |grad|_inf={np.abs(g0).max():.3e}; Hessian eig min/max="
          f"{ev.min():.3e}/{ev.max():.3e} (neg-def => local max); cond="
          f"{abs(ev.max()/ev.min()) if ev.min()!=0 else np.inf:.2e}", flush=True)

    def newton(x0, nsteps=25):
        x = np.asarray(x0, float).copy()
        f = float(lnl_loud(jnp.asarray(x)))
        for _ in range(nsteps):
            g = np.asarray(g_fn(jnp.asarray(x))); H = np.asarray(H_fn(jnp.asarray(x)))
            # ascend: step = -H^-1 g on the CONCAVE objective; damp if H not neg-def
            lam = 0.0
            for _try in range(40):
                Hd = H - lam * np.eye(len(x))            # push toward neg-def
                try:
                    step = -np.linalg.solve(Hd, g)
                except np.linalg.LinAlgError:
                    lam = max(1e-6, lam * 10); continue
                xn = x + step
                fn = float(lnl_loud(jnp.asarray(xn)))
                if np.isfinite(fn) and fn > f:           # accept ascent
                    x, f = xn, fn; break
                lam = max(1e-6, lam * 10)
            else:
                break                                     # no ascent -> stop
        return x, f

    print(f"\n  {'delta':>7s} {'start_sky':>10s} {'newton_sky':>11s} {'off':>8s} {'dlnL_vs_truth':>14s} {'lock':>5s}", flush=True)
    rng = np.random.default_rng(3)
    rows = []
    for d in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
        pert = rng.normal(size=nlp) * scale[:nlp] * d
        x0 = src_t[:nlp] + pert
        sky0 = max(ang_deg(x0[8*k], x0[8*k+1], src_t[8*k], src_t[8*k+1]) for k in range(n_loud))
        xN, fN = newton(x0)
        skyN = max(ang_deg(xN[8*k], xN[8*k+1], src_t[8*k], src_t[8*k+1]) for k in range(n_loud))
        off = float(np.abs(xN - src_t[:nlp]).max() / scale[:nlp].max())
        lock = skyN < 0.006
        rows.append((d, sky0, skyN, off, fN - lnL_truth, lock))
        print(f"  {d:7.0e} {sky0:10.4f} {skyN:11.5f} {off:8.4f} {fN-lnL_truth:+14.3f} {str(lock):>5s} "
              f"({time.time()-t0:.0f}s)", flush=True)

    locks = [d for d, *_ , lock in rows if lock]
    print(f"\n  NEWTON PULL-IN RADIUS (sky<0.006deg): <= {max(locks) if locks else 0.0:.0e}", flush=True)
    np.savez(f"{CWT}/trackB_L2c_newton.npz",
             grad_truth=g0, hess_eig=ev, lnL_truth=lnL_truth,
             rows=np.array([[r[0], r[1], r[2], r[3], r[4], float(r[5])] for r in rows]),
             pull_newton=(max(locks) if locks else 0.0))
    print(f"  saved trackB_L2c_newton.npz ({time.time()-t0:.0f}s)", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}", flush=True)
    sys.exit(main())
