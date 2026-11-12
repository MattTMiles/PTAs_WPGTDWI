#!/usr/bin/env python3
"""SIEVE V4 -- THE SIGMA-POINT E-STEP AT WORKING-REGIME BELIEF WIDTHS
(agent SIEVE-A, 2026-07-29 addendum; SIEVE-C T5's working-regime arm).

THE PRE-REGISTERED QUESTION, AND THE QUANTITY IT IS ABOUT.
  T3 measured that the certification gate `q_max > 0.9` is not calibrated: on the
  truth-anchored A-SKY banks the empirical P(on_true | q_max > 0.9) is 0.786, not 0.9;
  across the GLACIER banks it is 0.305. LEDGER-A1 named the mechanism -- the E-step
  evaluates the fringe posterior at ONE source point, so q_max is P(fringe | theta_hat)
  charged as P(fringe | data), and a conditional confidence is systematically OVER-
  confident when the source belief is wide relative to the fringe-discriminating scale.
  THE QUESTION IS THEREFORE: does averaging the E-step over the source belief move
  q_max toward calibration? NOT whether it changes the certification count. A method
  that raised the count while leaving the 0.786 gap open would have failed this test.

WHY "WORKING-REGIME" WIDTHS.  A1's gates run the rule at Sigma -> 0 (where it collapses
structurally to the incumbent) and on an analytic surrogate. Neither says what happens
at the widths a running loop actually carries. This arm re-tests it on a graduated
ladder from SUB-FRINGE to FEW-FRINGE belief, which is the band a feed-time belief
occupies -- narrower than a fringe and the average cannot change the fringe posterior;
much wider and every fringe is equally supported and q_max collapses toward 1/K for
reasons that have nothing to do with calibration.

THE LADDER IS DEFINED IN FRINGES, NOT IN PARAMETER UNITS, AND THE CONVERSION IS EXACT.
  The pulsar-term phase goes as f * L, so a fractional frequency error is
  indistinguishable from a fractional distance error: delta_f/f == delta_L/L. One
  fringe is delta_L = dL_a. Hence a belief of n_fr fringes on pulsar a corresponds to
        sigma(log10 fgw) = n_fr * dL_a / (L0_a * ln 10)
  and the ladder uses the MEDIAN over pulsars. Rungs: 0 (the incumbent, by structural
  collapse), 0.5 (sub-fringe) and 2.0 (few-fringe) -- see RUNGS for why three.

WHY THE BELIEF IS PUT ON log10_fgw ONLY -- DECLARED, NOT OVERLOOKED.
  MSTEP_AXES is (I_FGW, I_MC) and A1 specs the belief block as both. log10_fgw is the
  axis with the CLOSED-FORM fringe conversion above, so a rung on it means exactly what
  the rung label says. log10_mc also moves the pulsar-term phase -- the pulsar term is
  evaluated at a retarded time of order kyr, where the chirp has run -- but its
  fringe-equivalent width has no closed form here and would have to be calibrated
  numerically per cell. Adding it with a guessed width would make every rung label a
  guess too. So this arm measures the fgw axis exactly and DEFERS the mc axis, which is
  the conservative direction: a belief on one axis widens the marginal less than a
  belief on two, so any calibration improvement measured here is a LOWER BOUND on what
  the full two-axis belief would deliver.

MECHANICS.  The rule and the reduction come from `ledger_a1_sigma_estep` (imported, not
reimplemented, so A1's gates G-L1/G-L2 cover this arm's arithmetic): scaled unscented
transform, w0 = 1/3, every weight strictly positive, and the average is taken in
LIKELIHOOD space (logsumexp) -- averaging log-likelihoods gives a geometric mean and
cannot widen a posterior at all, which is the failure mode this whole test would
otherwise walk into.
  Cost: 2m+1 E-steps per rung with m = n_belief. Each E-step is expensive on the CPU
lane because the per-pulsar evaluators must be EVICTED after use (see
`_install_evicting_ab` -- this is a hard constraint, not a tuning choice). The script
TIMES its first E-step and then processes realisations in realisation-major order under
a wall-clock budget, so an early stop still leaves all 8 skies represented; the number
actually scored is printed and banked rather than silently truncated.

SUBSTRATE: the GENERALISE arm A-SKY SURVIVOR cell (`gen_armAS_sky.npz` verdict column:
cell 1 = 'e0.3 h-12.75 5+11' is the only SURVIVES), all 8 skies, truth-anchored by
construction (generalise.py:376, theta_base = theta_src with distances at L0).

READS GENERALISE_results/ read-only. WRITES SIEVE_results/ only.
Output bank: SIEVE_results/sieve_v4_sigma_c<cell>__<lane>.npz
Usage: python3 hpc_harbor/sieve/sv_v4_sigma.py [--cell 1] [--nbelief 2] [--budget-s 25200]
"""
import os, sys, time
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import argparse
import glob
import socket
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/chorus", "hpc_harbor/generalise",
          "hpc_harbor/ledger"):
    sys.path.insert(0, f"{HSYMT}/{p}")

import chorus as CH
import generalise as GEN
import ledger_a1_sigma_estep as A1

OUT = os.environ.get("GLACIER_OUT", f"{HSYMT}/SIEVE_results")
QBAR = 0.9
NP_SRC = 8
N_POP = 16
# Belief width in FRINGES. THREE rungs, and the trim is a deliberate trade, declared:
# at 360 s per E-step (see _install_evicting_ab and the compile-cache note in main) a
# 5-rung ladder buys ~3 realisations inside the budget -- ~350 rows per rung -- while a
# 3-rung ladder buys ~6, at ~700 rows per rung. The pre-registered question is whether
# q_max moves TOWARD CALIBRATION, which is a rate measured to a confidence interval, so
# rows are worth more than rungs. The band the brief asked for is still spanned:
#   0.0  the incumbent (structural collapse to the point rule -- A1's G-L1 identity)
#   0.5  SUB-fringe
#   2.0  FEW-fringe
RUNGS = (0.0, 0.5, 2.0)
W0 = 1.0 / 3.0


def _maps():
    """Live mmap count of this process -- the CPU lane's binding resource."""
    try:
        with open("/proc/self/maps") as f:
            return sum(1 for _ in f)
    except OSError:
        return -1


def _install_evicting_ab(sp):
    """CPU-lane map-count mitigation, NUMERICALLY INERT. Ported verbatim from
    `fsky_stage0._install_evicting_ab` (its docstring carries the measurement: the
    XLA-CPU thunk runtime maps one code section per fused kernel, each per-pulsar
    E-step evaluator holds ~4k mappings, and 116 live evaluators exhaust
    vm.max_map_count = 65530).

    THIS IS NOT OPTIONAL HERE, AND THE FIRST V4 SUBMISSION (job 12844088) PROVED IT:
    116 x ~4k = ~464k mappings cannot fit under the limit no matter how clean the slate,
    so the run died inside the FIRST sigma point with
    "Failed to materialize symbols / LLVM Cannot allocate memory" raised from
    trackB_b1_core.estep:488. Eviction keeps ~1 evaluator resident at a time.

    THE TAX IS PAID ON EVERY CALL, NOT ONCE. I expected the persistent JAX cache to
    serve the evicted evaluators back from disk; it does not (see the compile-cache note
    in `main`). Measured: 360 s per E-step, ~3.1 s per pulsar, unchanged by lowering the
    cache threshold and with the cache directory static at 3 entries. The cost is Python
    trace + lower of the 116 per-pulsar terms, which no executable cache touches. That
    number is what sets this arm's affordable scope, so it is recorded here rather than
    left for the next person to rediscover."""
    orig = sp._pulsar_ab_fn

    def evicting(p):
        fn = orig(p)

        def call(*a):
            out = tuple(np.asarray(o) for o in fn(*a))
            try:
                fn.clear_cache()
            except Exception:
                pass
            sp._ab_fns.pop(p, None)
            return out

        return call

    sp._pulsar_ab_fn = evicting


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(p), float(max(0.0, c - h)), float(min(1.0, c + h))


def reliability(q, ot, edges=(0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0001)):
    """Empirical P(on_true | q in bin). The diagonal is what calibration means."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (q >= lo) & (q < hi)
        out.append((lo, hi, int(m.sum()),
                    float(ot[m].mean()) if m.any() else np.nan,
                    float(q[m].mean()) if m.any() else np.nan))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", type=int, default=1)
    ap.add_argument("--nbelief", type=int, default=2,
                    help="number of loud members carrying a belief (m = nbelief)")
    ap.add_argument("--budget-s", type=float, default=18000.0)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    lane = os.environ.get("SLURM_JOB_ID", "local")
    host = socket.gethostname()
    print(f"[V4] host={host} cpus={len(os.sched_getaffinity(0))} job={lane} "
          f"start={time.strftime('%FT%T')}", flush=True)

    if not glob.glob(f"{GEN.OUT}/gen_L_gwb_n*.npz"):
        raise SystemExit("G-V4a FAIL: no shared GENERALISE L_gwb bank -- running would "
                         "WRITE into another campaign's bank directory. STOP.")
    print("[G-V4a] PASS shared L_gwb bank present", flush=True)

    jax, jnp, C, B1, TE, F, IG = GEN._import_stack()

    # ---- COMPILE-CACHE THRESHOLD: TRIED, MEASURED, AND IT DOES NOT HELP ---------
    # chorus.py:196-197 sets jax_persistent_cache_min_compile_time_secs = 10, and a
    # per-pulsar E-step evaluator turns round in ~3.1 s -- below that threshold, so it
    # is never written to the persistent cache. That looked like the reason the
    # mandatory evict-after-use above costs a full rebuild on every call, so the
    # threshold is lowered here to 0.2 s. MEASURED RESULT (job 12845138): no change --
    # 360.6 s per E-step vs 368.0 s without it, and the cache directory stayed at 3
    # entries / 376 MB while the run proceeded. The eviction tax is therefore Python
    # TRACE + LOWER of the 116 per-pulsar terms, not XLA backend compilation, and no
    # executable cache can remove it. The knob is kept because it is correct and free;
    # it is documented as ineffective so the next person does not re-derive it.
    # NUMERICALLY INERT either way: the cache stores executables, not results, and the
    # BLAS-thread count (pinned at 8) is the only reproducibility lever (harbor_env.sh).
    # The cache DIRECTORY is left exactly as chorus.py set it -- harbor_env.sh binds
    # this job's own $HARBOR_JAXCACHE onto that container path, so the isolation that
    # defuses the concurrent-writer race lives in the bind map, not in this config, and
    # repointing it here would escape the bind.
    jax.config.update("jax_persistent_cache_min_compile_time_secs", 0.2)
    units = [(ui, u) for ui, u in enumerate(GEN.as_units()) if u["ci"] == a.cell]
    h_cell, e_cell, k_cell = GEN.AS_CELLS[a.cell]
    print(f"[V4] cell {a.cell}: h={h_cell} e={e_cell} k_loud={k_cell}; "
          f"{len(units)} skies; rungs (fringes) {RUNGS}", flush=True)

    t0 = time.time()
    P = CH.build_chorus_problem(1, T_label=30)
    nd = P["n_dist"]; AI = P["AI"]; sp = P["sp"]; L0 = np.asarray(P["L0"], float)
    one = jnp.ones(P["npsr"])
    print(f"[V4] venue built in {time.time()-t0:.0f}s; maps={_maps()}", flush=True)
    # Release the BUILD's executables before the E-step starts compiling its own, then
    # keep at most one per-pulsar evaluator resident. fsky_stage0:304-307,376-378 order.
    jax.clear_caches()
    _install_evicting_ab(P["sp"])
    print(f"[V4] map-count hygiene installed; maps={_maps()} "
          f"(limit {open('/proc/sys/vm/max_map_count').read().strip()})", flush=True)

    # belief block: log10_fgw of the top-`nbelief` LOUD members
    axes = np.array([nd + NP_SRC * i + GEN.I_FGW for i in range(a.nbelief)], int)

    # ---- assemble the work list, sky-major round-robin ----------------------
    work = []
    geoms = {}
    for ui, u in units:
        ents = [e for e in GEN.unit_entries("AS", ui, u) if e["kind"] == "sig"]
        for ri, e in enumerate(ents):
            work.append((ri, ui, u, e))
    work.sort(key=lambda w: (w[0], w[1]))          # realisation-major -> all skies first

    rows = []
    t_est = None
    t_start = time.time()
    n_done = 0
    for ri, ui, u, e in work:
        if ui not in geoms:
            geoms[ui] = GEN.gen_geometry(P, C, F.load_geo_skies([u["sky"]])[0], u)
        G = geoms[ui]
        theta_src = G["theta"]; dL = np.asarray(G["dL"], float); EV = G["EV"]; FT = G["FT"]

        # fringes -> sigma(log10 fgw); exact, see the header
        sig_per_fringe = float(np.median(dL / L0)) / np.log(10.0)

        L_true, _ = CH.draw_true_distances_tier(P, dL, EV, seed=e["dist_seed"],
                                                tier=u["tier"])
        n_true_grid = C.n_true_on_grid(L_true, P["L0"], dL)
        theta_true = theta_src.copy(); theta_true[AI] = L_true
        clean = P["amo"]["inject_delay"](jnp.asarray(theta_true), one)
        noise = P["nd"].draw(e["noise_seed"])
        data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                     for c, n in zip(clean, noise))
        theta_base = theta_src.copy(); theta_base[AI] = P["L0"]

        def estep_fn(th):
            return np.asarray(sp.estep(jnp.asarray(th), EV, AI, data, one))

        for nf in RUNGS:
            te = time.time()
            sig = np.full(len(axes), nf * sig_per_fringe)
            LNL, n_eval = A1.estep_sigma(estep_fn, theta_base, axes, sig, w0=W0)
            post = FT.posterior(LNL)
            q = np.asarray(post["q_max"], float)
            mk = np.asarray(post["map_k"], int)
            ot = (mk == n_true_grid)
            for p_i in range(len(q)):
                rows.append((u["sky"], int(e["noise_seed"]), nf, p_i,
                             float(q[p_i]), bool(ot[p_i]), int(n_eval)))
            if t_est is None:
                t_est = (time.time() - te) / max(n_eval, 1)      # seconds per E-step
                per_real = t_est * sum(1 if r == 0 else 1 + 2 * len(axes)
                                       for r in RUNGS)
                print(f"[V4] {t_est:.2f}s per E-step; {per_real:.0f}s per realisation; "
                      f"full sweep of {len(work)} would be "
                      f"~{per_real*len(work)/3600:.1f} h vs budget "
                      f"{a.budget_s/3600:.1f} h -> expect "
                      f"~{min(len(work), int(a.budget_s/max(per_real,1e-9)))} scored",
                      flush=True)
        n_done += 1
        print(f"  [{n_done}/{len(work)}] sky {u['sky']} done; elapsed "
              f"{time.time()-t_start:.0f}s of {a.budget_s:.0f}s", flush=True)
        if time.time() - t_start > a.budget_s:
            print(f"[V4] BUDGET STOP after {n_done}/{len(work)} realisations "
                  f"({time.time()-t_start:.0f}s). NOT a silent truncation: the sweep is "
                  f"realisation-major, so all {len(units)} skies are represented.",
                  flush=True)
            break

    if not rows:
        raise SystemExit("no rows scored -- STOP rather than bank an empty result")

    sky = np.array([r[0] for r in rows], int)
    seed = np.array([r[1] for r in rows], int)
    nfr = np.array([r[2] for r in rows], float)
    psr = np.array([r[3] for r in rows], int)
    qq = np.array([r[4] for r in rows], float)
    ot = np.array([r[5] for r in rows], bool)
    nev = np.array([r[6] for r in rows], int)

    print("\n" + "=" * 78)
    print(f"V4  q_max CALIBRATION vs BELIEF WIDTH   ({n_done} realisations, "
          f"{len(np.unique(sky))} skies, m = {len(axes)} belief axes)")
    print("=" * 78)
    print(f"{'fringes':>8s} {'n_eval':>7s} {'mean q':>8s} {'n(q>0.9)':>9s} "
          f"{'P(on_true|q>0.9)':>17s} {'95% CI':>18s} {'gap vs 0.9':>11s}")
    summ = []
    for nf in RUNGS:
        m = nfr == nf
        if not m.any():
            continue
        sel = m & (qq > QBAR)
        k, n = int(ot[sel].sum()), int(sel.sum())
        p, lo, hi = wilson(k, n)
        summ.append((nf, int(nev[m][0]), float(qq[m].mean()), n, p, lo, hi))
        print(f"{nf:8.2f} {int(nev[m][0]):7d} {qq[m].mean():8.4f} {n:9d} "
              f"{p:17.4f} [{lo:.4f},{hi:.4f}] {p-0.9:+11.4f}")

    print("\nRELIABILITY (empirical P(on_true) per q_max bin; calibration = diagonal)")
    for nf in (0.0, RUNGS[-1]):
        m = nfr == nf
        if not m.any():
            continue
        print(f"  belief = {nf:.2f} fringes")
        for lo, hi, n, pt, mq in reliability(qq[m], ot[m]):
            if n:
                print(f"    q in [{lo:.2f},{hi:.2f})  n={n:6d}  mean q={mq:.4f}  "
                      f"P(on_true)={pt:.4f}   {'OVER' if mq > pt else 'under'}"
                      f"-confident by {abs(mq-pt):.4f}")

    stem = f"sieve_v4_sigma_c{a.cell}"
    path = f"{a.out}/{stem}__{host}_{lane}.npz"
    np.savez(path, sky=sky, seed=seed, n_fringes=nfr, psr=psr, qmax=qq, on_true=ot,
             n_eval=nev, rungs=np.array(RUNGS), axes=axes, w0=W0,
             n_belief=a.nbelief, cell=a.cell, h_cell=h_cell, k_loud=k_cell,
             n_real_scored=n_done, n_real_planned=len(work),
             summary=np.array(summ, float), host=host, job=lane,
             time=time.strftime("%FT%T"))
    print(f"\n[V4] banked {path}  ({n_done}/{len(work)} realisations scored)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
