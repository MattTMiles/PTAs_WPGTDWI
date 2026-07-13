"""SAMPLER — NUTS on the fringe-marginalised conditional posterior + the gate ladder.

Modes
-----
  g12       zero-noise Asimov: does the posterior concentrate at truth (g1), and do two
            independent chains/seeds agree (g2, R-hat < 1.01 on every dim)?
  g3        ONE noisy realisation at the VLBI cell: coverage of truth + the per-pulsar
            fringe posteriors from POSTERIOR-AVERAGED sources vs the E-step's at fixed source.
  sbc       simulation-based calibration: prior draw -> simulate -> posterior -> rank stats.

Conventions
-----------
  * The target is `sampler_core.MargJax` (gated value-identical to `trackB_b1.B1Marg`).
  * Sampling happens in a logit-transformed unit box (the spec's physical ranges), so no
    proposal ever leaves the box; the log-Jacobian is carried exactly.
  * THE MASS MATRIX IS SEEDED FROM THE TARGET'S OWN HESSIAN at the init point (dense), because
    the registration dims (log10_fgw, log10_mc) carry the pulsar-term lever arm and the rest do
    not: the condition number is ~1e10-class and window-adapted diagonal mass will not find it
    inside any affordable warmup. What it takes is REPORTED, not hidden.
  * Background runs: nohup + per-chunk npz checkpoints + finiteness-on-values.
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

import numpyro
from numpyro.infer import MCMC, NUTS
from numpyro.diagnostics import split_gelman_rubin, effective_sample_size

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import sampler_core as S
import trackB_b1 as B1
import trackB_b1_core as C
import trackB_estimator as TE

RES = S.RES


# ============================================================
# preconditioner
# ============================================================
def hessian_precond(logp_z, z0, ndim, floor_rel=1e-8, verbose=True, fd=True, h=1e-4):
    """Dense inverse mass matrix from the target's own Hessian at z0.

    inv_mass = (-H)^-1 (the Laplace covariance), eigen-floored to stay positive definite.
    Returns (inv_mass, eigenvalues of -H, condition number).

    `fd=True` builds H by central differences OF THE GRADIENT (2*ndim gradient evals). This
    is required for the fringe-marginalised target, whose lnL_marg is an XLA-opaque primitive
    with a custom first-order VJP (MargJax.as_primitive) -- jax.hessian cannot traverse it.
    The gradient it differences is the one gated against FD, so the Hessian inherits that gate.
    """
    t0 = time.time()
    if fd:
        g = jax.jit(jax.grad(logp_z))
        H = np.zeros((ndim, ndim))
        for j in range(ndim):
            e = np.zeros(ndim); e[j] = h
            gp = np.asarray(g(jnp.asarray(z0 + e)))
            gm = np.asarray(g(jnp.asarray(z0 - e)))
            H[:, j] = (gp - gm) / (2 * h)
    else:
        H = np.asarray(jax.hessian(logp_z)(jnp.asarray(z0)))
    H = 0.5 * (H + H.T)
    A = -H                                            # curvature (pos-def at a max)
    w, V = np.linalg.eigh(A)
    cond_raw = float(np.max(np.abs(w)) / max(np.min(np.abs(w)), 1e-300))
    wpos = np.clip(w, floor_rel * np.max(np.abs(w)), None)
    inv_mass = (V * (1.0 / wpos)) @ V.T
    inv_mass = 0.5 * (inv_mass + inv_mass.T)
    if verbose:
        print(f"   Hessian at init: {time.time()-t0:.0f}s; eig(-H) in "
              f"[{w.min():.3e}, {w.max():.3e}]; cond = {cond_raw:.2e}; "
              f"neg eigs {int(np.sum(w <= 0))}/{ndim}", flush=True)
    return inv_mass, w, cond_raw


# ============================================================
# NUTS driver (chunked, checkpointed)
# ============================================================
def run_chain(M, prior, x_init, seed, n_warm, n_samp, chunk, tag, max_tree_depth=8,
              dense_mass=True, adapt_mass=True, precond=True, target_accept=0.8):
    pot, logp = S.make_potential(M, prior)
    z0 = np.asarray(prior.to_z(jnp.asarray(x_init)))
    ndim = len(z0)

    inv_mass = None
    eig = None
    cond = np.nan
    if precond:
        inv_mass, eig, cond = hessian_precond(logp, z0, ndim)

    kern = NUTS(potential_fn=pot, dense_mass=dense_mass, max_tree_depth=max_tree_depth,
                target_accept_prob=target_accept, adapt_mass_matrix=adapt_mass,
                inverse_mass_matrix=inv_mass, adapt_step_size=True)
    mcmc = MCMC(kern, num_warmup=n_warm, num_samples=chunk, num_chains=1,
                progress_bar=False, jit_model_args=True)

    key = jax.random.PRNGKey(seed)
    Z = []
    extras = []
    t0 = time.time()
    got = 0
    first = True
    while got < n_samp:
        if first:
            mcmc.run(key, init_params=jnp.asarray(z0), extra_fields=("diverging",
                                                                     "num_steps",
                                                                     "adapt_state.step_size"))
            first = False
            print(f"   [{tag}] warmup ({n_warm}) + chunk done in {time.time()-t0:.0f}s",
                  flush=True)
        else:
            mcmc.post_warmup_state = mcmc.last_state
            mcmc.run(mcmc.post_warmup_state.rng_key, extra_fields=("diverging", "num_steps",
                                                                   "adapt_state.step_size"))
        z = np.asarray(mcmc.get_samples())
        ex = mcmc.get_extra_fields()
        Z.append(z)
        extras.append(ex)
        got += z.shape[0]
        X = np.array([np.asarray(prior.to_x(jnp.asarray(zz))[0]) for zz in np.concatenate(Z)])
        C._finite(f"{tag} samples", X)
        div = int(np.sum(np.concatenate([np.asarray(e["diverging"]) for e in extras])))
        nst = np.concatenate([np.asarray(e["num_steps"]) for e in extras])
        np.savez(f"{RES}/chain_{tag}.npz", X=X, Z=np.concatenate(Z), seed=seed,
                 n_warm=n_warm, div=div, num_steps=nst, cond=cond,
                 eig=eig if eig is not None else np.zeros(1),
                 step_size=float(np.asarray(ex["adapt_state.step_size"])[-1]),
                 inv_mass=inv_mass if inv_mass is not None else np.zeros((1, 1)),
                 wall=time.time() - t0)
        print(f"   [{tag}] {got}/{n_samp} samples; {time.time()-t0:.0f}s; "
              f"divergences {div}; median leapfrog {np.median(nst):.0f}; "
              f"step {float(np.asarray(ex['adapt_state.step_size'])[-1]):.3e}", flush=True)
    Zc = np.concatenate(Z)
    X = np.array([np.asarray(prior.to_x(jnp.asarray(zz))[0]) for zz in Zc])
    lnl = np.array([M.value(x) for x in X[::max(1, len(X) // 200)]])
    return dict(X=X, Z=Zc, div=div, num_steps=np.concatenate(
        [np.asarray(e["num_steps"]) for e in extras]), cond=cond, eig=eig,
        wall=time.time() - t0, lnl_sub=lnl, inv_mass=inv_mass)


def diagnose(chains, labels, x_truth):
    """R-hat / ESS per dim over >=2 chains; concentration vs truth."""
    Xs = np.stack([c["X"] for c in chains])                # (nchain, nsamp, ndim)
    rh = np.array([float(split_gelman_rubin(Xs[:, :, j])) for j in range(Xs.shape[2])])
    ess = np.array([float(effective_sample_size(Xs[:, :, j])) for j in range(Xs.shape[2])])
    med = np.median(Xs.reshape(-1, Xs.shape[2]), axis=0)
    sd = np.std(Xs.reshape(-1, Xs.shape[2]), axis=0)
    z = (med - x_truth) / np.where(sd > 0, sd, np.inf)
    print(f"\n   {'param':>22s} {'truth':>12s} {'post med':>12s} {'post sd':>10s} "
          f"{'z':>7s} {'R-hat':>7s} {'ESS':>7s}", flush=True)
    for j, lab in enumerate(labels):
        print(f"   {lab:>22s} {x_truth[j]:12.6f} {med[j]:12.6f} {sd[j]:10.3e} "
              f"{z[j]:7.2f} {rh[j]:7.4f} {ess[j]:7.0f}", flush=True)
    return dict(rhat=rh, ess=ess, med=med, sd=sd, z=z)


# ============================================================
# post-hoc per-pulsar fringe / distance posterior at posterior draws
# ============================================================
def fringe_posterior_at_draws(P, B1M, M, X, n_draw=40, seed=0):
    """q_bar_p(n) = mean over posterior draws of the per-pulsar fringe posterior.

    This is how distances are recovered: never sampled, always marginalised, then read off
    the fringe posteriors AT posterior draws (mission convention).
    """
    rng = np.random.default_rng(seed)
    pick = rng.choice(len(X), size=min(n_draw, len(X)), replace=False)
    qs, qmax, mapk = [], [], []
    for i, s in enumerate(pick):
        src = M.src_fixed.copy()
        src[M.active_idx] = X[s]
        post, LNL, lnref = B1M.posterior(src)
        qs.append([np.asarray(q) for q in B1M.FT.q])
        qmax.append(post["q_max"]); mapk.append(post["map_k"])
    # average the per-pulsar fringe posteriors over draws (same fringe lattice for all draws)
    npsr = P["npsr"]
    qbar = []
    for a in range(npsr):
        Q = np.stack([q[a] for q in qs])
        qbar.append(Q.mean(axis=0))
    return dict(qbar=qbar, qmax=np.array(qmax), mapk=np.array(mapk), pick=pick)


def cert_from_qbar(P, B1M, qbar, n_true):
    """q_max / MAP fringe / P(true fringe) from the posterior-averaged fringe posteriors."""
    npsr = P["npsr"]
    qm = np.zeros(npsr); mk = np.zeros(npsr, int); pt = np.full(npsr, np.nan)
    for a in range(npsr):
        j = int(np.argmax(qbar[a]))
        qm[a] = qbar[a][j]
        mk[a] = B1M.FT.uk[a][j]
        hit = np.where(B1M.FT.uk[a] == n_true[a])[0]
        if len(hit):
            pt[a] = qbar[a][hit[0]]
    return qm, mk, pt


# ============================================================
# g1 / g2 : zero-noise concentration + two-chain agreement
# ============================================================
def mode_g12(a):
    print("=== g1/g2: ZERO-NOISE Asimov posterior on the conditional target ===", flush=True)
    cell = S.build_cell(seed=None, tier=a.tier, noise=False, arm="A")
    P = cell["P"]
    src_truth = cell["theta_truth"][P["n_dist"]:].copy()
    loud = list(range(C.N_LOUD))
    idx = S.active_index(loud, S.ACTIVE_AXES, P["ncw"])
    labels = S.active_labels(loud, S.ACTIVE_AXES)
    M = S.MargJax(P, cell["data"], cell["mask"], idx, src_truth, tier=a.tier)
    lo, hi = S.box_bounds(idx, P["ncw"])
    prior = S.BoxPrior(lo, hi)
    x_tr = src_truth[idx]
    scale = TE.phi_scale(P)[idx]

    t0 = time.time(); v, g = M.value_and_grad(x_tr); t_g = time.time() - t0
    print(f"  lnL_marg(truth) = {v:.4f}; warm grad {t_g:.2f}s; ndim {len(idx)}", flush=True)

    chains = []
    for c in range(a.chains):
        # chain 0 starts AT truth; chain 1+ start displaced (an independent seed AND an
        # independent start -- g2 must not be a two-chain rerun of the same trajectory)
        rng = np.random.default_rng(1000 + c)
        x0 = x_tr if c == 0 else np.clip(
            x_tr + a.disp * scale * rng.normal(size=len(idx)), lo + 1e-6, hi - 1e-6)
        print(f"\n-- chain {c} (seed {a.seed0+c}, start offset "
              f"{0.0 if c == 0 else a.disp:.0e} scaled) --", flush=True)
        ch = run_chain(M, prior, x0, seed=a.seed0 + c, n_warm=a.warm, n_samp=a.samp,
                       chunk=a.chunk, tag=f"g12_c{c}", max_tree_depth=a.depth)
        chains.append(ch)

    d = diagnose(chains, labels, x_tr)
    # g1: sigma(log10 mc) vs IGNITE-2's banked profile widths (1e-5..1e-4 dex class)
    mci = [j for j, l in enumerate(labels) if l.endswith("log10_mc")]
    mcs = ", ".join("%.2e" % d["sd"][j] for j in mci)
    print(f"\n  [g1] sigma(log10_mc) per loud source: {mcs} dex", flush=True)
    print(f"       IGNITE-2 soft-loop PROFILE widths (different object): 1e-5..1e-4 dex",
          flush=True)
    print(f"  [g1] |median - truth| / sigma per dim: max {np.max(np.abs(d['z'])):.2f}",
          flush=True)
    rmax = float(np.max(d["rhat"]))
    print(f"  [g2] max R-hat = {rmax:.4f} -> {'PASS' if rmax < 1.01 else 'FAIL'} (<1.01)",
          flush=True)
    np.savez(f"{RES}/g12_summary.npz", rhat=d["rhat"], ess=d["ess"], med=d["med"],
             sd=d["sd"], z=d["z"], labels=np.array(labels), x_truth=x_tr,
             cond=np.array([c["cond"] for c in chains]),
             div=np.array([c["div"] for c in chains]),
             wall=np.array([c["wall"] for c in chains]),
             num_steps=np.concatenate([c["num_steps"] for c in chains]))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["g12"])
    ap.add_argument("--tier", default="lit")
    ap.add_argument("--warm", type=int, default=150)
    ap.add_argument("--samp", type=int, default=250)
    ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--chains", type=int, default=2)
    ap.add_argument("--depth", type=int, default=6)
    ap.add_argument("--disp", type=float, default=1e-3,
                    help="chain-1 start displacement, scaled units")
    ap.add_argument("--seed0", type=int, default=990001)
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    if a.mode == "g12":
        return mode_g12(a)
    return 2


if __name__ == "__main__":
    sys.exit(main())
