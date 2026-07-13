"""SAMPLER — the SBC harness (simulation-based calibration) on the fringe-marginalised posterior.

THE TEST THE PROGRAMME HAS NEVER RUN (STORY S10.1, deliverable (ii))
--------------------------------------------------------------------
SBC is the only test that can validate `q_max` as a PROBABILITY rather than as a score:

    for i in 1..N:
        theta_i   ~ prior                        (the spec box, ACTIVE dims; rest at truth)
        n_true_i  ~ prior over fringe centres    (the arm-B distance draw)
        d_i       ~ p(d | theta_i, n_true_i)     (CW + white + RN + HD-GWB, all drawn)
        {theta_s} ~ p(theta | d_i)               (NUTS on lnL_marg -- fringes marginalised)
        r_i[j]    = #{s : theta_s[j] < theta_i[j]}          (rank of truth, per dim)

If the posterior is calibrated, each r_i is Uniform{0..S}. Deviations diagnose:
    - a U-shaped rank histogram   -> posterior too NARROW (overconfident)
    - a central-peaked histogram  -> posterior too WIDE
    - a sloped histogram          -> biased.

AND THE FRINGE-INTEGER RANK, which is the one that matters here: for each pulsar, the
posterior-averaged fringe posterior q_bar_p(n) is compared with the TRUE fringe n_true_p via
the probability-integral-transform of a DISCRETE distribution (randomised PIT, Czado et al.):
    u_p = P(n < n_true) + U * P(n = n_true),   U ~ U(0,1)
which is Uniform(0,1) under calibration. THIS is the direct test of whether `q_max` means what
the certification criterion has always assumed it means.

DEV PHASE (this module): N = 10 smoke test on the 4090 -- shapes, cost, plumbing, and the
rank machinery. PRODUCTION (N >= 200): ACCRE, per the handoff block in SAMPLER_dev.md.

Env: jug-gpu, x64. nohup + per-realisation npz checkpoints + finiteness-on-values.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, argparse, glob
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import sampler_core as S
import sampler_run as R
import trackB_b1 as B1
import trackB_b1_core as C
import trackB_estimator as TE

RES = S.RES


def draw_prior_source(P, src_truth, active_idx, seed):
    """theta_act ~ uniform on the spec box; everything else held at truth."""
    rng = np.random.default_rng(seed)
    lo, hi = S.box_bounds(active_idx, P["ncw"])
    src = np.asarray(src_truth, float).copy()
    src[active_idx] = lo + (hi - lo) * rng.uniform(size=len(active_idx))
    return src


def rank_of_truth(X, x_true):
    """r[j] = #{s: X[s,j] < x_true[j]}  (0..S), the SBC rank statistic."""
    return (X < x_true[None, :]).sum(axis=0)


def fringe_pit(qbar, uk, n_true, rng):
    """Randomised PIT of the TRUE fringe under the posterior-averaged fringe posterior."""
    out = np.full(len(qbar), np.nan)
    for a in range(len(qbar)):
        hit = np.where(uk[a] == n_true[a])[0]
        if not len(hit):
            continue                                  # true fringe outside the sampled window
        j = int(hit[0])
        below = float(np.sum(qbar[a][:j])) if j > 0 else 0.0
        at = float(qbar[a][j])
        out[a] = below + rng.uniform() * at
    return out


def one_realisation(P, B1M_cache, i, seed, active_idx, src_truth, args):
    t0 = time.time()
    # ---- prior draw (source) + arm-B fringe truth (distances) ----
    src_i = draw_prior_source(P, src_truth, active_idx, seed=seed)
    L_true, n_true, u_true = C.draw_true_distances(P, seed=seed + 7)
    theta_i = np.concatenate([L_true, src_i])

    # ---- simulate ----
    mask1 = P["mask_one"]
    clean = P["amo"]["inject_delay"](jnp.asarray(theta_i), jnp.asarray(mask1))
    nz = P["nd"].draw(seed=seed + 11)
    data = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(clean, nz))
    C._finite(f"sbc{i} data", *[np.asarray(d) for d in data])

    # ---- posterior ----
    M = S.MargJax(P, data, mask1, active_idx, src_i, tier=args.tier)
    lo, hi = S.box_bounds(active_idx, P["ncw"])
    prior = S.BoxPrior(lo, hi)
    x_true = src_i[active_idx]
    # init AT the truth of this realisation (dev-phase convention: SBC tests the posterior's
    # calibration, not the sampler's global search; a cold init would confound the two and the
    # referendum (f = 6.9e-7) already says global search fails on the data alone)
    ch = R.run_chain(M, prior, x_true, seed=seed, n_warm=args.warm, n_samp=args.samp,
                     chunk=args.samp, tag=f"sbc{i}", max_tree_depth=args.depth)
    X = ch["X"]

    # ---- rank statistics ----
    r = rank_of_truth(X, x_true)
    B1M = B1M_cache[0]
    fp = R.fringe_posterior_at_draws(P, B1M, M, X, n_draw=args.ndraw, seed=seed)
    rng = np.random.default_rng(seed + 13)
    pit = fringe_pit(fp["qbar"], B1M.FT.uk, n_true, rng)
    qm, mk, pt = R.cert_from_qbar(P, B1M, fp["qbar"], n_true)

    rec = dict(i=i, seed=seed, x_true=x_true, rank=r, S=len(X),
               div=ch["div"], cond=ch["cond"], wall=time.time() - t0,
               pit=pit, qmax=qm, mapk=mk, ptrue=pt, n_true=n_true,
               post_med=np.median(X, axis=0), post_sd=np.std(X, axis=0),
               num_steps=ch["num_steps"])
    np.savez(f"{RES}/sbc_{i:03d}.npz", **{k: np.asarray(v) for k, v in rec.items()})
    print(f"  [sbc {i}] {time.time()-t0:.0f}s; ranks {r}; div {ch['div']}; "
          f"n_cert(qbar>0.9) {int(np.sum(qm > 0.9))}", flush=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--warm", type=int, default=200)
    ap.add_argument("--samp", type=int, default=300)
    ap.add_argument("--depth", type=int, default=7)
    ap.add_argument("--ndraw", type=int, default=24)
    ap.add_argument("--tier", default="lit")
    ap.add_argument("--seed0", type=int, default=880001)
    args = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    print(f"=== SBC SMOKE TEST, N = {args.n} (production N>=200 -> ACCRE) ===", flush=True)

    P = C.build_b1_problem()
    src_truth = P["theta_truth"][P["n_dist"]:].copy()
    idx = S.active_index(list(range(C.N_LOUD)), S.ACTIVE_AXES, P["ncw"])
    labels = S.active_labels(list(range(C.N_LOUD)), S.ACTIVE_AXES)

    # one B1Marg for the fringe-posterior read-out (data swapped per realisation is NOT
    # possible -- B1Marg bakes data; rebuild per realisation is the cost we pay, cached below)
    recs = []
    for i in range(args.n):
        seed = args.seed0 + 100 * i
        done = f"{RES}/sbc_{i:03d}.npz"
        if os.path.exists(done):
            print(f"  [sbc {i}] already banked -- skip", flush=True)
            recs.append(dict(np.load(done, allow_pickle=True)))
            continue
        # B1Marg must see THIS realisation's data
        src_i = draw_prior_source(P, src_truth, idx, seed=seed)
        L_true, n_true, _ = C.draw_true_distances(P, seed=seed + 7)
        theta_i = np.concatenate([L_true, src_i])
        clean = P["amo"]["inject_delay"](jnp.asarray(theta_i), jnp.asarray(P["mask_one"]))
        nz = P["nd"].draw(seed=seed + 11)
        data = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(clean, nz))
        B1M = B1.B1Marg(P, data, P["mask_one"])
        recs.append(one_realisation(P, [B1M], i, seed, idx, src_truth, args))

    # ---- SBC summary ----
    Rk = np.array([np.asarray(r["rank"]) for r in recs])
    Sm = int(np.asarray(recs[0]["S"]))
    print(f"\n=== SBC RANKS (N={len(recs)}, S={Sm}); uniform under calibration ===", flush=True)
    print(f"   {'param':>22s} {'ranks':>40s}", flush=True)
    for j, lab in enumerate(labels):
        print(f"   {lab:>22s} {str(list(Rk[:, j])):>40s}", flush=True)
    pit = np.concatenate([np.asarray(r["pit"])[np.isfinite(np.asarray(r["pit"]))]
                          for r in recs])
    print(f"\n   fringe PIT (all pulsars x realisations, N = {len(pit)}): "
          f"mean {pit.mean():.3f} (0.5), sd {pit.std():.3f} (0.289); "
          f"KS-vs-U(0,1) D = {ks_uniform(pit):.3f}", flush=True)
    np.savez(f"{RES}/sbc_summary.npz", rank=Rk, S=Sm, pit=pit, labels=np.array(labels),
             wall=np.array([float(np.asarray(r["wall"])) for r in recs]))
    return 0


def ks_uniform(u):
    u = np.sort(np.asarray(u))
    n = len(u)
    if n == 0:
        return np.nan
    i = np.arange(1, n + 1)
    return float(max(np.max(i / n - u), np.max(u - (i - 1) / n)))


if __name__ == "__main__":
    sys.exit(main())
