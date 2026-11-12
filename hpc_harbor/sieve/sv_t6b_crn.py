#!/usr/bin/env python3
"""SIEVE T6b -- COMMON RANDOM NUMBERS: PAIRED vs INDEPENDENT frozen-vs-live.

REPORT-ONLY. Measures ONE number: how much of the variance of an ARM DIFFERENCE is
removed by sharing the noise seed across the two arms. That number is the sample-size
savings factor for every future two-arm comparison in this project.

    Var_indep(D) / Var_paired(D)  =  1 / (1 - rho)   when the arms have equal variance,
    rho = corr(Y_frozen(s), Y_live(s)) across seeds.

WHAT IS AND IS NOT BORROWED FROM PHOENIX. The freeze DIAL is reimplemented here,
verbatim in behaviour (per-slot freeze-after-first-fit, hpc_harbor/frozen/frozen_mstep.py
`mstep_frozen`), but this driver NEVER imports that module, never writes to
FROZEN_results/, and never uses PHOENIX's stems, rung, seeds or lane. It is a toy at
n_src = 16, T = 30, 3 iterations, on its own seed base, banking to SIEVE_results/.
The CRN variance ratio is a property of the ESTIMATOR, not of the rung; measuring it on
a toy is the whole point of costing it at 0.5 GPU-hr rather than at campaign scale.

WHAT IS DELIBERATELY OMITTED FROM THE TOY (declared, so the scope is auditable):
  * null-offender FLOORS. They are the expensive part of run_cell step (f) and they do
    not enter the arm difference under measurement. floor = 0 throughout; the scoreboard
    columns are still cut, so n_cert here is a floorless count and is NOT comparable to
    any banked campaign count. It is used only as one of several arm metrics.
  * the STOP machinery, banking per iteration, and checkpoint/resume. This is not a
    campaign cell and produces no campaign bank.

THE DESIGN (12 runs from ONE venue build):
    6 noise seeds s_1..s_6, each run through BOTH arms.
    PAIRED    D_s   = Y_frozen(s)   - Y_live(s)          -> Var over the 6 diffs
    INDEP (a) analytic: Var(Y_frozen) + Var(Y_live)      -> what unpaired sampling costs
    INDEP (b) crossed:  D_ij = Y_f(s_i) - Y_l(s_j), i!=j -> all 30 mismatched pairs
    Both independent estimates are reported; (a) is the estimator-theory answer, (b) is
    the empirical check on the same 12 runs (its 30 values are not independent, so it is
    quoted as a cross-check, not as the headline).

METRICS Y (one ratio reported per metric):
    a_bg      the drain amplitude at the last iteration -- what PHOENIX's frozen-vs-live
              comparison actually reads (S4.15.1 item 2: "0.03-0.30 dex below baseline")
    a_bg_sig  the drain width at that fit
    logL_end  the joint log-likelihood at the final template
    fgw_hat   the M-step's log10 f_gw for the first-fed slot (the dial's own axis)
    mc_hat    the M-step's log10 M_c for the first-fed slot (the dial's other axis)
    n_res     resolved (fed) member count

  THE SCOREBOARD IS NOT CALLED. `CertScoreboard.columns` routes to
  `spark3.estep_per_target`, which is the documented XLA-CPU `vm.max_map_count`
  exhaustion hazard (fsky_stage0._install_evicting_ab) -- and the certification columns
  do not enter the arm difference this test measures. The drain IS the frozen-vs-live
  readout. Dropping the E-step is what keeps T6b on the general CPU lane and therefore
  off every claimed GPU lane; the GPU-hr allocated to it is not spent.

Output bank: SIEVE_results/sieve_t6b_crn__<lane>.npz
Usage: python3 hpc_harbor/sieve/sv_t6b_crn.py [--nsrc 16] [--iters 3] [--seeds 6]
"""
import os, sys, time
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import argparse
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/spark", "hpc_harbor/glacier"):
    sys.path.insert(0, f"{HSYMT}/{p}")

import glacier_loop as GL
import glacier_pop as GP

# SIEVE's OWN seed base -- disjoint from GL.SEED_NOISE_BASE (9_500_000) and from every
# campaign stem, so a toy run can never be mistaken for or collide with a banked cell.
SEED_BASE = 6_100_000
T = 30
ARM = "none"
SKY = 0
METRICS = ("a_bg", "a_bg_sig", "logL_end", "fgw_hat", "mc_hat", "n_res")


def run_arm(IG, TE, C, G, ndraw, jnp, noise_seed, frozen, n_iter):
    """One toy cell. `frozen=True` applies per-slot freeze-after-first-fit to the
    M-step and nothing else; every other stage is stock glacier_loop code."""
    amo = G["amo"]; nd = amo["n_dist"]; slots = G["slots"]
    ones = jnp.ones(amo["npsr"])
    # the arm-B distance truth, drawn exactly as CertScoreboard.draw_truth does but
    # without constructing the scoreboard (whose E-step is the CPU-lane hazard)
    from cw_helpers import compute_mode_spacing
    import stagec_anchor_a2 as A2
    ent = G["ent_psrs"]
    names = [pe.name for pe in ent]
    L0 = np.array([pe.pdist[0] for pe in ent])
    dL = np.array([min(compute_mode_spacing(s[GL.I_COSGT], s[GL.I_GWPHI], s[GL.I_FGW],
                                            G["disco_psrs"][ai].pos)
                       for s in slots) for ai in range(len(ent))])
    Pd = dict(npsr=len(ent), L0=L0, priors=A2.load_real_priors(names), ent_psrs=ent)
    EV = C.build_EV(Pd, dL)
    L_true, _ = IG.draw_true_distances_tier(Pd, dL, EV,
                                            seed=noise_seed + 10_000_000, tier=GL.TIER)
    theta_true = np.asarray(amo["theta_truth"], float).copy()
    theta_true[np.arange(nd)] = L_true
    clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
    noise = ndraw.draw(noise_seed, components=("white", "rn"))
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                 for c, n in zip(clean, noise))

    theta_rec = theta_true.copy()
    led = GL.PromoteLedger(slots)
    bf = GL.BackgroundFit(amo, band_log10f=GL.BAND_CAMPAIGN)
    box_sigma = np.array([1.0, np.pi, 1.0, 0.5, 0.25, 0.25, np.pi,
                          0.5 * np.pi]) / np.sqrt(3.0)
    a_grid = np.linspace(GP.A_TARGET_LOG10 - 1.0, GP.A_TARGET_LOG10 + 1.0, 41)
    scale = TE.phi_scale({"ncw": 1})
    prev_cert_idx, prev_q = np.zeros(0, int), np.zeros(amo["npsr"])
    held = {}                       # slot -> frozen 8-vector (the dial's state)
    out = {}

    for it in range(n_iter):
        carried = np.where(~led.fed)[0]
        sig_opt, sig_pes, F_ii = GL.fisher_conditional(amo, jnp, theta_rec, carried, nd,
                                                       box_sigma, data, ones)
        ratio = np.max(sig_opt / sig_pes, axis=1)
        new = [int(k) for k in np.where((ratio < GL.GATE_RATIO) & ~led.fed)[0]]
        # frontier-v2 data-support term, stock form
        if new:
            th_off = GL.theta_with_absent(theta_rec, nd, carried)
            for k in list(new):
                th_on = GL.theta_with_absent(theta_rec, nd, np.setdiff1d(carried, [k]))
                ll = np.asarray(amo["logL_batch_theta"](
                    jnp.asarray(np.stack([th_on, th_off])), data, ones))
                if not float(ll[0] - ll[1]) > 0.0:
                    new.remove(k)
        for k in new:
            led.promote(k, slots[k], iteration=it)
        fed_idx = np.where(led.fed)[0]

        prof = bf.profile(theta_rec, data, ones, led.fed.astype(float), a_grid)

        if len(fed_idx):
            carried_now = np.where(~led.fed)[0]
            th_eval = GL.theta_with_absent(theta_rec, nd, carried_now)

            def marg_fn(srcs):
                ths = np.tile(th_eval, (len(srcs), 1))
                ths[:, nd:] = srcs
                return np.asarray(amo["logL_batch_theta"](jnp.asarray(ths), data, ones))

            if frozen:
                fresh = [k for k in fed_idx if int(k) not in held]
                if fresh:
                    sh, _, _ = GL.mstep_quadratic(marg_fn, th_eval[nd:].copy(), fresh,
                                                  scale)
                    for k in fresh:
                        held[int(k)] = sh[GL.NP_SRC * k: GL.NP_SRC * (k + 1)].copy()
                src_hat = th_eval[nd:].copy()
                for k in fed_idx:
                    src_hat[GL.NP_SRC * k: GL.NP_SRC * (k + 1)] = held[int(k)]
            else:
                src_hat, _, _ = GL.mstep_quadratic(marg_fn, th_eval[nd:].copy(),
                                                   fed_idx, scale)
            theta_rec = theta_rec.copy()
            for k in fed_idx:
                theta_rec[nd + GL.NP_SRC * k: nd + GL.NP_SRC * (k + 1)] = \
                    src_hat[GL.NP_SRC * k: GL.NP_SRC * (k + 1)]

        th_end = GL.theta_with_absent(theta_rec, nd, np.where(~led.fed)[0])
        ll_end = float(np.asarray(amo["logL_batch_theta"](
            jnp.asarray(th_end[None, :]), data, ones))[0])
        k0 = int(fed_idx[0]) if len(fed_idx) else -1
        out = dict(a_bg=float(prof["ahat"]), a_bg_sig=float(prof["sig"]),
                   logL_end=ll_end,
                   fgw_hat=(float(theta_rec[nd + GL.NP_SRC * k0 + GL.I_FGW])
                            if k0 >= 0 else np.nan),
                   mc_hat=(float(theta_rec[nd + GL.NP_SRC * k0 + GL.I_MC])
                           if k0 >= 0 else np.nan),
                   n_res=float(led.n_resolved()))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nsrc", type=int, default=16)
    ap.add_argument("--iters", type=int, default=3)
    ap.add_argument("--seeds", type=int, default=6)
    a = ap.parse_args()

    GL.check_affinity()
    t0 = time.time()
    jax, jnp, C, B1, TE, IG, F, FL = GL._stack()
    print(f"[T6b] stack {time.time()-t0:.1f}s", flush=True)

    pop = GL.draw_population(GL.SEED_POP_BASE + SKY, n_src=a.nsrc,
                             band_log10f=GL.BAND_CAMPAIGN)
    slots, n_harm, active, chan, n_clip = GL.embed_igniter(
        pop, GL.E_IGNITER[ARM], GL.VENUE_SPAN_S[T])
    pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
    t = time.time()
    G = GP.build_glacier_problem(T, pop_s, verbose=True)
    G["slots"] = slots
    print(f"[T6b] build {time.time()-t:.1f}s  n_slot={len(slots)}", flush=True)
    amo = G["amo"]
    # B1Split is built directly (the NoiseDrawer needs it); CertScoreboard is NOT
    # constructed -- see the module docstring on the CPU-lane E-step hazard.
    t = time.time()
    theta_truth = np.asarray(amo["theta_truth"], float)
    data_zero = amo["inject_delay"](jnp.asarray(theta_truth), jnp.ones(amo["npsr"]))
    sp = C.B1Split(amo, theta_truth, data_zero, np.ones(amo["npsr"]))
    ndraw = C.NoiseDrawer(sp, amo, lgwb_path=GL.LGWB_BANKS[T], strict=True)
    print(f"[T6b] split+drawer {time.time()-t:.1f}s  L_gwb {ndraw.lgwb_prov}",
          flush=True)

    seeds = [SEED_BASE + 1000 * i for i in range(a.seeds)]
    Yf, Yl = [], []
    for i, s in enumerate(seeds):
        for frozen, store in ((True, Yf), (False, Yl)):
            t = time.time()
            r = run_arm(IG, TE, C, G, ndraw, jnp, s, frozen, a.iters)
            store.append(r)
            print(f"[T6b] seed {i} ({s}) {'FROZEN' if frozen else 'LIVE  '}: "
                  + " ".join(f"{k}={r[k]:.5g}" for k in METRICS)
                  + f"  [{time.time()-t:.0f}s]", flush=True)

    print("\n[T6b] ===== CRN VARIANCE RATIO =====")
    print(f"   {a.seeds} seed pairs, n_src={a.nsrc}(+{n_harm} harm), "
          f"T={T}, {a.iters} iterations, arm={ARM}")
    print("   metric       n  Var(Yf)     Var(Yl)   Var_paired  Var_indep(a)  "
          "Var_indep(b)     rho   ratio(a)  ratio(b)")
    res = []
    for m in METRICS:
        yf = np.array([r[m] for r in Yf], float)
        yl = np.array([r[m] for r in Yl], float)
        # A seed that feeds NOTHING has no fitted member, so fgw_hat/mc_hat are NaN for
        # that seed in BOTH arms. Drop those pairs rather than let one NaN erase the
        # metric; the surviving pair count is printed so the reader sees the n.
        ok = np.isfinite(yf) & np.isfinite(yl)
        n_pair = int(ok.sum())
        if n_pair < 3:
            print(f"   {m:10s} only {n_pair} finite seed pairs -- not reported")
            res.append([m, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                        np.nan, np.nan])
            continue
        yf, yl = yf[ok], yl[ok]
        d = yf - yl
        vf, vl = float(yf.var(ddof=1)), float(yl.var(ddof=1))
        vp = float(d.var(ddof=1))
        via = vf + vl
        cross = np.array([yf[i] - yl[j] for i in range(len(yf))
                          for j in range(len(yl)) if i != j], float)
        vib = float(cross.var(ddof=1))
        if vf > 0 and vl > 0:
            rho = float(np.corrcoef(yf, yl)[0, 1])
        else:
            rho = float("nan")
        ra = via / vp if vp > 0 else float("inf")
        rb = vib / vp if vp > 0 else float("inf")
        print(f"   {m:10s} {n_pair:3d} {vf:10.3e} {vl:10.3e} {vp:11.3e} {via:12.3e} "
              f"{vib:12.3e} {rho:7.3f} {ra:9.2f} {rb:9.2f}")
        res.append([m, vf, vl, vp, via, vib, rho, ra, rb, n_pair])

    fin = [r for r in res if np.isfinite(r[7]) and np.isfinite(r[3]) and r[3] > 0]
    if fin:
        best = max(r[7] for r in fin)
        med = float(np.median([r[7] for r in fin]))
        print(f"\n   median ratio(a) over metrics with non-degenerate paired variance: "
              f"{med:.2f}x   (max {best:.2f}x)")
        print(f"   FLAG TO SUMMIT SS2 CONVENTIONS: "
              f"{'YES -- ratio > 1.5x' if med > 1.5 else 'no -- ratio <= 1.5x'}")
    else:
        med = float("nan")
        print("\n   every metric has zero paired variance: the arms are IDENTICAL on "
              "every seed at this toy scale -- the ratio is unmeasurable here and the "
              "toy must be made larger (more iterations, or a rung that late-feeds) "
              "before a savings factor can be quoted.")

    path = GP.bank_npz(
        "sieve_t6b_crn",
        note=("SIEVE T6b: common-random-numbers demo. PAIRED (shared noise seed across "
              "the frozen and live arms) vs INDEPENDENT sampling; the variance ratio of "
              "the arm difference is the sample-size savings factor for future two-arm "
              "comparisons. Toy scale, floorless, SIEVE's own seed base; the freeze dial "
              "is reimplemented (behaviour-verbatim) and PHOENIX's module, stems, lane "
              "and banks are untouched. REPORT-ONLY."),
        n_seeds=a.seeds, n_src=a.nsrc, n_harm=n_harm, n_slot=len(slots),
        t_label=T, n_iter=a.iters, arm=ARM, sky=SKY, seed_base=SEED_BASE,
        seeds=np.array(seeds), metrics=np.array(METRICS),
        yf=np.array([[r[m] for m in METRICS] for r in Yf]),
        yl=np.array([[r[m] for m in METRICS] for r in Yl]),
        result=np.array(res, dtype=object),
        result_keys=np.array(["metric", "var_frozen", "var_live", "var_paired",
                              "var_indep_analytic", "var_indep_crossed", "rho",
                              "ratio_analytic", "ratio_crossed", "n_pairs"]),
        median_ratio=med)
    print(f"\nbanked -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
