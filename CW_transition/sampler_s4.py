"""SAMPLER g4 — S-4 REVISITED: does an honest marginal posterior restore coverage with a GWB?

THE STANDING STOP-POINT (RING S-4, reports/RING_q1_modernized.md §5.2/§7)
------------------------------------------------------------------------
With a GWB in the data (scenario C), the sky credible region UNDERCOVERS even at EXACT
distances: `inside90` = 0.40-0.50 against a nominal 0.90, while scenarios A (white) and B
(white+RN) hold 0.90-1.00 at the same tier. RING could not separate the two candidate causes
and forbade quoting the number until a sampler had been run:

  (i)  PHYSICS   -- a correlated background shifts which fringe/peak the MAP selects;
  (ii) ESTIMATOR -- RING's "posterior" is a PROFILE SLICE on a sky grid (every other parameter
                    held at truth), not a marginal posterior, and its credible regions are only
                    approximately calibrated exactly where the region is smallest (~1e-4 deg2).

WHAT THIS MODULE DOES, AND WHAT IT DOES NOT
-------------------------------------------
It runs the S-4 QUESTION in the B1 harness, which differs from RING's and is stated as such:
116 real pulsars (not the 30-pulsar ring), the pop-3000 loud sources, and -- the point --
a likelihood in which the GWB is an HD-correlated GP MARGINALISED ANALYTICALLY AND EXACTLY
inside the covariance (frozen hyperparameters), plus a genuine MARGINAL posterior from NUTS
over every free source parameter. So the estimator leg (ii) is removed by construction:
  * no profile slice   -- all 8 params/loud source are sampled, the sky among them;
  * no Laplace         -- the credible region is read off the posterior draws;
  * no grid            -- RING S-1's wrong-peak refinement cannot occur.
Distances are FIXED AT TRUTH in injection and recovery (RING's tier3 / "exact" tier), so there
are no fringes to marginalise and lnL is the plain array likelihood.

If coverage returns to nominal -> RING's undercoverage was (ii), the estimator.
If it does not                 -> it is (i), the physics, and the number stands as a property
                                  of a GWB-contaminated sky posterior.

THE STATISTIC (chosen before looking, and it is not `inside90`)
--------------------------------------------------------------
At N = 5 realisations `inside90` can only take the values 0/5..5/5 and is nearly uninformative.
The per-realisation CREDIBLE LEVEL AT TRUTH is recorded instead:
    c = fraction of posterior draws whose 2-D sky density is BELOW the density at truth
      = the smallest credible level whose region contains truth.
Calibrated inference => c ~ U(0,1) across realisations; systematic undercoverage => c piles up
near 1. `inside90 = (c < 0.90)` is reported beside it, so the RING comparison survives.
Density estimated by a Gaussian KDE on the 2-D sky marginal of each loud source (RING's object).

Env: jug-gpu, x64. Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, argparse
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp
from scipy.stats import gaussian_kde

import numpyro
from numpyro.infer import MCMC, NUTS

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import sampler_core as S
import sampler_run as R
import trackB_b1_core as C
import trackB_estimator as TE

RES = S.RES
ALL_AXES = list(range(8))          # sky + everything: the full loud-source vector


class ExactDistTarget:
    """Plain array lnL at FIXED (true) distances -- RING's exact tier. Differentiable, cheap."""

    def __init__(self, P, data, mask, theta_truth, active_idx):
        self.P = P
        self.active_idx = np.asarray(active_idx, int)
        self.ndim = len(self.active_idx)
        self.theta0 = np.asarray(theta_truth, float).copy()      # distances = TRUE distances
        self.ndist = P["n_dist"]
        logL = P["amo"]["logL"]
        th0 = jnp.asarray(self.theta0)
        idx = jnp.asarray(self.active_idx + self.ndist)          # into the full theta vector
        d, m = data, jnp.asarray(mask)

        def f(x):
            theta = th0.at[idx].set(x)
            return logL(theta, d, m)
        self.lnL = jax.jit(f)
        self.lnL_vg = jax.jit(jax.value_and_grad(f))

    def value(self, x):
        return float(self.lnL(jnp.asarray(x, float)))

    def x_truth(self):
        return self.theta0[self.ndist:][self.active_idx]


def sky_credible_level(X, truth_x, loud, axes):
    """c = P(density < density(truth)) on each loud source's 2-D sky marginal (KDE)."""
    out = []
    for k, kk in enumerate(loud):
        j0 = k * len(axes) + axes.index(S.I_COSGT)
        j1 = k * len(axes) + axes.index(S.I_GWPHI)
        pts = X[:, [j0, j1]].T                                   # (2, N)
        kde = gaussian_kde(pts)
        d_s = kde(pts)
        d_t = float(kde(np.array([[truth_x[j0]], [truth_x[j1]]]))[0])
        c = float(np.mean(d_s < d_t))                            # credible level at truth
        out.append(dict(src=kk, c=c, inside90=bool(c < 0.90),
                        sky_sd=(float(np.std(X[:, j0])), float(np.std(X[:, j1]))),
                        truth=(truth_x[j0], truth_x[j1]),
                        med=(float(np.median(X[:, j0])), float(np.median(X[:, j1])))))
    return out


def run_one(P, seed, with_gwb, n_warm, n_samp, loud, max_tree_depth, verbose=True):
    tt = P["theta_truth"].copy()                    # arm A: true distances == prior means
    mask1 = P["mask_one"]
    comps = ("white", "rn", "gwb") if with_gwb else ("white", "rn")
    data_clean = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))
    nz = P["nd"].draw(seed=seed, components=comps)
    data = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(data_clean, nz))
    C._finite("s4 data", *[np.asarray(d) for d in data])

    idx = S.active_index(loud, ALL_AXES, P["ncw"])
    T = ExactDistTarget(P, data, mask1, tt, idx)
    lo, hi = S.box_bounds(idx, P["ncw"])
    prior = S.BoxPrior(lo, hi)
    x_tr = T.x_truth()

    pot = jax.jit(lambda z: -(T.lnL(prior.to_x(z)[0]) + prior.to_x(z)[1]))
    logp = jax.jit(lambda z: T.lnL(prior.to_x(z)[0]) + prior.to_x(z)[1])
    z0 = np.asarray(prior.to_z(jnp.asarray(x_tr)))
    # Hessian by FD ON THE GRADIENT, not jax.hessian. jax.hessian here is forward-over-reverse
    # through the WHOLE 116-pulsar array likelihood and its host-side compile footprint
    # OOM-KILLED the box (kernel: "Out of memory: Killed process ... anon-rss:34 GB"). The FD
    # form costs 2*ndim gradient evals and a few hundred MB.
    inv_mass, eig, cond = R.hessian_precond(logp, z0, len(z0), verbose=verbose, fd=True, h=1e-3)

    kern = NUTS(potential_fn=pot, dense_mass=True, max_tree_depth=max_tree_depth,
                target_accept_prob=0.8, inverse_mass_matrix=inv_mass)
    mcmc = MCMC(kern, num_warmup=n_warm, num_samples=n_samp, num_chains=1, progress_bar=False)
    t0 = time.time()
    mcmc.run(jax.random.PRNGKey(seed), init_params=jnp.asarray(z0),
             extra_fields=("diverging", "num_steps"))
    Z = np.asarray(mcmc.get_samples())
    ex = mcmc.get_extra_fields()
    X = np.array([np.asarray(prior.to_x(jnp.asarray(z))[0]) for z in Z])
    C._finite("s4 samples", X)
    div = int(np.sum(np.asarray(ex["diverging"])))
    lev = sky_credible_level(X, x_tr, loud, ALL_AXES)
    if verbose:
        for L in lev:
            print(f"     src{L['src']}: c(truth) = {L['c']:.3f}  inside90 = {L['inside90']}  "
                  f"sky sd = ({L['sky_sd'][0]:.2e}, {L['sky_sd'][1]:.2e})", flush=True)
        print(f"     {n_samp} draws, {time.time()-t0:.0f}s, divergences {div}, "
              f"cond(-H) {cond:.1e}", flush=True)
    return dict(X=X, levels=lev, div=div, cond=cond, wall=time.time() - t0,
                x_truth=x_tr, seed=seed, with_gwb=with_gwb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nreal", type=int, default=5)
    ap.add_argument("--warm", type=int, default=400)
    ap.add_argument("--samp", type=int, default=1500)
    ap.add_argument("--depth", type=int, default=8)
    ap.add_argument("--nloud", type=int, default=1, help="loud sources sampled (8 dims each)")
    ap.add_argument("--seed0", type=int, default=770001)
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    print("=== g4 / S-4 REVISITED: sky coverage with a properly marginalised GWB ===", flush=True)

    P = C.build_b1_problem()
    loud = list(range(a.nloud))
    tt = P["theta_truth"]
    nd = P["n_dist"]
    for k in loud:
        print(f"  loud{k}: log10_fgw = {tt[nd+8*k+S.I_FGW]:.3f}, log10_h = "
              f"{tt[nd+8*k+S.I_H]:.3f}, log10_mc = {tt[nd+8*k+S.I_MC]:.3f}", flush=True)

    out = {}
    for arm, gwb in [("C", True), ("B", False)]:
        print(f"\n-- scenario {arm} ({'white+RN+GWB' if gwb else 'white+RN'}), "
              f"{a.nreal} realisations, EXACT distances --", flush=True)
        recs = []
        for i in range(a.nreal):
            seed = a.seed0 + 1000 * i + (0 if gwb else 500)
            print(f"   [{arm}{i}] seed {seed}", flush=True)
            r = run_one(P, seed, gwb, a.warm, a.samp, loud, a.depth)
            recs.append(r)
            cs = np.array([[L["c"] for L in rr["levels"]] for rr in recs])
            np.savez(f"{RES}/s4_scen{arm}.npz",
                     c=cs, inside90=(cs < 0.90),
                     div=np.array([rr["div"] for rr in recs]),
                     cond=np.array([rr["cond"] for rr in recs]),
                     wall=np.array([rr["wall"] for rr in recs]),
                     seeds=np.array([rr["seed"] for rr in recs]),
                     X_last=recs[-1]["X"], x_truth=recs[-1]["x_truth"])
        out[arm] = recs

    print("\n=== S-4 SCORECARD (credible level at truth; calibrated => U(0,1)) ===", flush=True)
    for arm in ("C", "B"):
        cs = np.array([[L["c"] for L in r["levels"]] for r in out[arm]]).ravel()
        ins = float(np.mean(cs < 0.90))
        print(f"  scenario {arm}: c = [{', '.join(f'{c:.3f}' for c in cs)}]  "
              f"median {np.median(cs):.3f}  inside90 = {ins:.2f} (nominal 0.90)", flush=True)
    print("\n  RING S-4 for reference: scenario C tier3 inside90 = 0.40-0.50; "
          "scenarios A/B tier3 = 0.90-1.00", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
