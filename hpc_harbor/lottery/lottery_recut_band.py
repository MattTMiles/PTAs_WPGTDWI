"""LOTTERY re-cut: the CHORUS m1 e=0.3 (lit) grade as a BAND over a multi-basis ensemble.

WHY. LOTTERY tier-2's basis control showed that m1e03_lit -- RECUT's headline "a single e=0.3
member does NOT switch the count on", the cell called "the single most expensive consequence of
the floor defect in the whole repo" -- FLIPS below -> CONFIRMED when its floor is re-estimated
from a fresh 100-null draw on a different host (floor 11.30 -> 9.91, Delta -1.40 nat ~ 1 sigma of
its own bootstrap error). A grade that reverses under a nuisance knob is not a measurement. This
driver replaces the point grade with a BAND over an ensemble of GWB-basis states, and applies a
PRE-REGISTERED two-sided reading.

THE BASIS KNOB. NoiseDrawer's L_gwb is an eigendecomposition of a near-degenerate Phi_gwb. Its
eigenvectors are only defined up to rotation within degenerate blocks, and that rotation depends
on the BLAS thread count and on the host's BLAS kernel dispatch (see the standing
noise-draw-thread-count hazard: cpus-per-task=8 is NECESSARY but NOT SUFFICIENT for bit-repro --
the host matters too). Rotating the basis leaves the physics invariant but redraws the noise
realisations, so each basis state is an independent-but-legitimate estimate of the same floor.
The spread ACROSS basis states is the systematic the point grade omitted.

ENSEMBLE (what is actually reachable on the general queue):
  B0  dgx      /  8 thr  -- CHORUS's banked cut (recut_chorus.npz). READ ONLY; the reserved
                            dgx_iacc share is never touched by this campaign.
  B1  hgx03    /  8 thr  -- LOTTERY tier-2's own cut (already banked, prefix `lot`).
  B2  hgx03    /  4 thr  -- thread-count basis variation.
  B3  hgx03    / 16 thr  -- thread-count basis variation.
  B4  hgx03    /  2 thr  -- thread-count basis variation (strongest perturbation).
>= 3 basis states as briefed. THREE distinct hosts were NOT reachable: the QOS
`taylor_group_account_batch_gpu` grants only gres h200=8 and gh200=2; the only h200 node is
hgx03, and the gh200 nodes are aarch64 while the harbor container + jug-gpu env are x86-64
(jaxlib/nvidia-*-cu12 wheels). So the brief's stated fallback -- "host + thread-count
variations" -- is the branch taken, with the dgx bank supplying the one genuine second host.

CELLS. The target m1e03_lit, and m2e03_lit as a FAR-FROM-BAR CONTRAST: if the band rule is sane
it must leave m2e03 (count 2.8, ~2 nat clear of the bar) CONFIRMED in every basis state while the
near-bar cell straddles. The contrast is what separates "this cell is fragile" from "the
estimator is noisy".

SEEDS ARE IDENTICAL ACROSS BASIS STATES -- the work plan is imported verbatim from
lottery_tier2._plan() and filtered to the control cells, so the ONLY thing that differs between
states is the basis. Files are keyed by prefix (`lot`, `lotb04`, ...) so states never collide and
skip-on-exist resume works per state.

MODES:
  run      : GPU. env LOT_BASIS=<tag> selects the prefix; cpus-per-task sets the thread count.
             Writes a provenance sidecar LOTTERY_results/basis_<tag>.json (host, cpus, BLAS
             threads, lgwb provenance string) -- the convention this campaign proposes adopting
             repo-wide.
  analyze  : CPU. Refits the floor per (cell, basis state) under the RECUT adopt convention,
             applies the pre-registered band rule, emits reports/lottery_recut_band.npz.

Usage (inside harbor container):
  LOT_BASIS=b04 harbor_py hpc_harbor/lottery/lottery_recut_band.py run
  harbor_py hpc_harbor/lottery/lottery_recut_band.py analyze
"""
import os, sys, glob, json, time, socket
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/hpc_harbor/chorus")
sys.path.insert(0, f"{HSYMT}/hpc_harbor/atlas")
sys.path.insert(0, f"{HSYMT}/hpc_harbor/lottery")
sys.path.insert(0, f"{HSYMT}/CW_transition")

import chorus as CH
import lottery_tier2 as T2

LOT_OUT = T2.LOT_OUT
REPORTS = T2.REPORTS
TIER = T2.TIER

# ---- the ensemble ----------------------------------------------------------
# tag -> (prefix, host_expected, cpus, kind). B0 is read from the CHORUS bank, never re-run.
BASIS = [
    dict(tag="B0", prefix=None,     host="dgx (banked)", cpus=8,  src="recut_chorus.npz"),
    dict(tag="B1", prefix="lot",    host="hgx03",        cpus=8,  src="LOTTERY tier-2"),
    dict(tag="B2", prefix="lotb04", host="hgx03",        cpus=4,  src="band ensemble"),
    dict(tag="B3", prefix="lotb16", host="hgx03",        cpus=16, src="band ensemble"),
    dict(tag="B4", prefix="lotb02", host="hgx03",        cpus=2,  src="band ensemble"),
]
BY_PREFIX = {b["prefix"]: b for b in BASIS if b["prefix"]}
CELLS = [(1, 0.30, "e03"), (2, 0.30, "e03")]      # target + far-from-bar contrast
TARGET_TAG = CH.cell_tag(1, "e03", TIER)

# ---- RECUT floor convention (reused verbatim from tier 2) ------------------
adopt = T2.adopt
_offender = T2._offender
QBAR = T2.QBAR


# ============================================================================
# GPU run mode
# ============================================================================
def _plan_for(prefix):
    """lottery_tier2's plan, filtered to the control cells. Seeds are IDENTICAL to tier 2 --
    the basis is the only free variable."""
    keep = {(k, lab) for (k, _e, lab) in CELLS}
    return [e for e in T2._plan() if (e["n_ecc"], e["edist"]) in keep]

def _path_of(prefix, e):
    return (f"{LOT_OUT}/{prefix}_{e['kind']}_"
            f"{CH.cell_tag(e['n_ecc'], e['edist'], e['tier'])}_g{e['geo_id']}_n{e['noise_seed']}.npz")

def mode_run(btag=None):
    # argv wins over env: harbor_py runs apptainer with --cleanenv, which strips LOT_BASIS from
    # the container, so the sbatch passes the basis tag positionally.
    btag = (btag or os.environ.get("LOT_BASIS", "")).strip()
    prefix = f"lot{btag}" if btag else None
    if prefix not in BY_PREFIX:
        print(f"FATAL: LOT_BASIS={btag!r} -> prefix {prefix!r} not in ensemble "
              f"{sorted(BY_PREFIX)}", flush=True)
        return 2
    spec = BY_PREFIX[prefix]

    CH.E_FIXED.update(T2.E_LABELS)
    CH.OUT = LOT_OUT
    os.makedirs(LOT_OUT, exist_ok=True)
    host = socket.gethostname()
    # apptainer runs with --cleanenv, so SLURM_CPUS_PER_TASK is NOT visible in here. The quantity
    # that actually sets the basis is the CPU count OpenBLAS sees = this process's affinity mask
    # (cgroup-limited by cpus-per-task). Read that directly and require it to match the ensemble
    # spec -- a mislabelled basis must never enter the bank.
    ncpu = len(os.sched_getaffinity(0))
    print(f"[BAND {spec['tag']}] prefix={prefix} host={host} affinity_cpus={ncpu} "
          f"(expected {spec['cpus']}) start={time.strftime('%FT%T')}", flush=True)
    if ncpu != spec["cpus"]:
        print(f"  !! REFUSING: affinity CPUs {ncpu} != ensemble spec {spec['cpus']}; the basis "
              f"label would be wrong.", flush=True)
        return 3

    jax, jnp, C, B1, TE, F, IG = CH._import_stack()
    print(f"jax {jax.__version__} {jax.devices()}", flush=True)

    # Same forced-recompute convention as tier 2: the frozen b1_L_gwb bank is the T=15 baseline
    # (3248^2) and shape-mismatches the T=30 extended GWB (5336^2). Forcing recompute is what
    # makes the BLAS thread count an actual basis knob -- the eigh runs here, in this job.
    _orig = C.load_or_build_L_gwb
    def _force_recompute(Phi_gwb, path=None, strict=False):
        return _orig(Phi_gwb, path="/nonexistent/lottery_forces_recompute.npz", strict=False)
    C.load_or_build_L_gwb = _force_recompute

    plan = _plan_for(prefix)
    todo = [e for e in plan if not os.path.exists(_path_of(prefix, e))]
    print(f"plan: {len(plan)}; already banked: {len(plan)-len(todo)}; to run: {len(todo)}", flush=True)
    if not todo:
        print("NOTHING TO DO (all banked)", flush=True)
        return 0

    skies = {g: (None if g < 0 else F.load_geo_skies([g])[0]) for g in set(CH.SKIES)}
    by_shape = {}
    for e in todo:
        by_shape.setdefault(e["n_ecc"], []).append(e)
    prov = {}
    for n_ecc in sorted(by_shape):
        P = CH.build_chorus_problem(n_ecc)
        prov[f"n_ecc={n_ecc}"] = str(P["nd"].lgwb_prov)
        print(f"  [NoiseDrawer] {P['nd'].lgwb_prov}", flush=True)
        for e in by_shape[n_ecc]:
            t0 = time.time()
            tag, msg, _ = CH.run_realisation(
                P, C, F, jnp, e["kind"], e["n_ecc"], e["edist"], e["tier"],
                e["geo_id"], skies[e["geo_id"]], e["noise_seed"], e["dist_seed"],
                no_cw=e.get("no_cw", False), prefix=prefix)
            print(f"  {tag:40s} {msg}  ({time.time()-t0:.1f}s)", flush=True)
        del P
        import gc; gc.collect()

    # ---- provenance sidecar: the convention this campaign proposes -----------
    side = f"{LOT_OUT}/basis_{spec['tag']}.json"
    with open(side, "w") as fh:
        json.dump(dict(basis_tag=spec["tag"], prefix=prefix, host=host,
                       affinity_cpus=ncpu, slurm_job=os.environ.get("SLURM_JOB_ID"),
                       lgwb_prov=prov, jax=str(jax.__version__),
                       devices=[str(d) for d in jax.devices()],
                       written=time.strftime("%FT%T")), fh, indent=2)
    print(f"provenance -> {side}", flush=True)
    print(f"RUN COMPLETE {time.strftime('%FT%T')}", flush=True)
    return 0


# ============================================================================
# CPU analyze mode
# ============================================================================
def _load_nulls(prefix, tag):
    fs = sorted(glob.glob(f"{LOT_OUT}/{prefix}_fnullN_{tag}_g*_n*.npz"))
    out = []
    for f in fs:
        d = np.load(f)
        out.append(_offender(d["dlnL_det"], d["lnK"], d["qmax"]))
    return np.array(out)

def _load_sig(prefix, tag):
    fs = sorted(glob.glob(f"{LOT_OUT}/{prefix}_sig_{tag}_g*_n*.npz"))
    D, K, Q, M, N = [], [], [], [], []
    for f in fs:
        d = np.load(f, allow_pickle=True)
        D.append(np.nan_to_num(np.asarray(d["dlnL_det"]), posinf=1e30))
        K.append(np.asarray(d["lnK"])); Q.append(np.asarray(d["qmax"]))
        M.append(np.asarray(d["mapk"])); N.append(np.asarray(d["n_true_grid"]))
    return dict(dlnL=np.array(D), lnK=np.array(K), qmax=np.array(Q),
                mapk=np.array(M), ntg=np.array(N), n=len(fs))

def _count_at(S, floor):
    cert = (S["dlnL"] > np.maximum(S["lnK"], floor)) & (S["qmax"] > QBAR)
    corr = cert & (S["mapk"] == S["ntg"])
    return float(corr.sum()) / S["n"]

# ---- THE PRE-REGISTERED BAND RULE ------------------------------------------
# The standing RECUT convention is ONE-SIDED: CONFIRMED iff count(floor+err) > 1, MARGINAL iff
# count(floor) > 1, else below. It has no way to say "the refutation is inside the error", so a
# cell whose floor moves by ~1 sigma flips its headline. The band rule closes that by reading the
# floor error in BOTH directions:
#     CONFIRMED_b  iff count(floor_b + err_b) >  1     (survives the pessimistic floor)
#     below_b      iff count(floor_b - err_b) <= 1     (fails even the optimistic floor)
#     MARGINAL_b   otherwise                            (the bar is inside 1 sigma -> undecided)
# Ensemble grade = unanimous label if all basis states agree, else MARGINAL (band), quoting
# [min,max] of floor and count across states. A cell within 1 sigma of its own bootstrap error
# NEITHER CONFIRMS NOR REFUTES.
def grade_band(count_lo, count_hi):
    """count_lo = count at (floor+err) [pessimistic], count_hi = count at (floor-err)."""
    if count_lo > 1.0:
        return "CONFIRMED"
    if count_hi <= 1.0:
        return "below"
    return "MARGINAL"

def mode_analyze():
    print("=" * 104)
    print("LOTTERY re-cut -- m1 e=0.3 (lit) as a BAND over a multi-basis ensemble")
    print("=" * 104)

    # B0: CHORUS's banked dgx cut
    z = np.load(f"{REPORTS}/recut_chorus.npz", allow_pickle=True)
    cols = [str(c) for c in z["cols"]]; tab = z["table"]
    tags = [str(t) for t in z["tags"]]; grades = [str(g) for g in z["grade"]]

    rows = []
    for (k, e, lab) in CELLS:
        ctag = CH.cell_tag(k, lab, TIER)
        for spec in BASIS:
            if spec["prefix"] is None:            # banked dgx
                if ctag not in tags:
                    continue
                i = tags.index(ctag)
                r = {c: float(tab[i][j]) for j, c in enumerate(cols)}
                floor, err = r["floor"], r["err"]
                # recut_chorus banks corr and corr_lo (=count at floor+err) but not the
                # optimistic side; recover it from the banked slope if unavailable.
                c_mid, c_lo = r["corr"], r["corr_lo"]
                c_hi = c_mid + (c_mid - c_lo)      # symmetric linearisation of count(floor)
                approx_hi = True
                nn = int(r.get("n", 0))
            else:
                oN = _load_nulls(spec["prefix"], ctag)
                S = _load_sig(spec["prefix"], ctag)
                if len(oN) == 0 or S["n"] == 0:
                    print(f"  {ctag} [{spec['tag']}]: INCOMPLETE (nulls={len(oN)}, sig={S['n']})",
                          flush=True)
                    continue
                zf = float(np.mean(oN == 0.0))
                floor, err, est = adopt(oN, zf)
                c_mid = _count_at(S, floor)
                c_lo = _count_at(S, floor + err)
                c_hi = _count_at(S, max(floor - err, 0.0))
                approx_hi = False
                nn = S["n"]
            g_band = grade_band(c_lo, c_hi)
            g_old = ("CONFIRMED" if c_lo > 1.0 else ("MARGINAL" if c_mid > 1.0 else "below"))
            rows.append(dict(cell=ctag, n_ecc=k, e=e, basis=spec["tag"], host=spec["host"],
                             cpus=spec["cpus"], src=spec["src"], floor=floor, err=err,
                             count=c_mid, count_lo=c_lo, count_hi=c_hi, n=nn,
                             grade_band=g_band, grade_pointwise=g_old, approx_hi=approx_hi))

    # ---- per-cell ensemble summary -----------------------------------------
    summary = []
    for (k, e, lab) in CELLS:
        ctag = CH.cell_tag(k, lab, TIER)
        rs = [r for r in rows if r["cell"] == ctag]
        if not rs:
            continue
        fl = np.array([r["floor"] for r in rs]); ct = np.array([r["count"] for r in rs])
        er = np.array([r["err"] for r in rs])
        labels = sorted({r["grade_band"] for r in rs})
        ens = labels[0] if len(labels) == 1 else "MARGINAL (band)"
        pw = sorted({r["grade_pointwise"] for r in rs})
        summary.append(dict(cell=ctag, n_states=len(rs),
                            floor_min=float(fl.min()), floor_max=float(fl.max()),
                            floor_spread=float(fl.max() - fl.min()),
                            err_mean=float(er.mean()),
                            count_min=float(ct.min()), count_max=float(ct.max()),
                            straddles_bar=bool(ct.min() <= 1.0 < ct.max()),
                            band_labels=labels, ensemble_grade=ens,
                            pointwise_labels=pw,
                            pointwise_unstable=bool(len(pw) > 1),
                            spread_over_err=float((fl.max() - fl.min()) / er.mean())))

    # ---- VARIANCE DECOMPOSITION: thread-count rotation vs host change -------
    # B1..B4 are the SAME host (hgx03) at 2/4/8/16 BLAS threads -> four genuinely different
    # eigenbases (distinct printed L_gwb fingerprints), same LAPACK build. B0 is a different
    # host. If the floor were sensitive only to the BASIS ROTATION, the within-host scatter
    # would be comparable to the host-to-host step. It is not -- which localises the systematic.
    decomp = []
    for (k, e, lab) in CELLS:
        ctag = CH.cell_tag(k, lab, TIER)
        rs = [r for r in rows if r["cell"] == ctag]
        wi = [r for r in rs if r["basis"] != "B0"]          # hgx03, 4 thread states
        b0 = next((r for r in rs if r["basis"] == "B0"), None)
        if len(wi) < 2 or b0 is None:
            continue
        f = np.array([r["floor"] for r in wi]); c = np.array([r["count"] for r in wi])
        errm = float(np.mean([r["err"] for r in wi]))
        d = dict(cell=ctag,
                 within_host_floor_mean=float(f.mean()), within_host_floor_sd=float(f.std(ddof=1)),
                 within_host_floor_range=float(f.max() - f.min()),
                 within_host_count_sd=float(c.std(ddof=1)),
                 within_host_count_vals=[float(x) for x in c],
                 cross_host_dfloor=float(f.mean() - b0["floor"]),
                 cross_host_dfloor_pct=float(100.0 * (f.mean() - b0["floor"]) / b0["floor"]),
                 cross_host_dcount=float(c.mean() - b0["count"]),
                 mean_bootstrap_err=errm,
                 within_over_err=float(f.std(ddof=1) / errm),
                 cross_over_err=float(abs(f.mean() - b0["floor"]) / errm))
        decomp.append(d)
    print("\n" + "=" * 104)
    print("VARIANCE DECOMPOSITION -- is the systematic the BASIS ROTATION (threads) or the HOST?")
    for d in decomp:
        print(f"  {d['cell']:12s} within-host (2/4/8/16 thr, 4 distinct bases): "
              f"floor {d['within_host_floor_mean']:.2f} sd={d['within_host_floor_sd']:.3f} "
              f"range={d['within_host_floor_range']:.2f} | counts={d['within_host_count_vals']}")
        print(f"  {'':12s} cross-host (hgx03 - dgx): Dfloor={d['cross_host_dfloor']:+.2f} nat "
              f"({d['cross_host_dfloor_pct']:+.1f}%)  Dcount={d['cross_host_dcount']:+.2f}")
        print(f"  {'':12s} vs mean bootstrap err {d['mean_bootstrap_err']:.2f}: "
              f"within/err={d['within_over_err']:.2f}  cross/err={d['cross_over_err']:.2f}  "
              f"-> {'HOST dominates (rotation is benign)' if d['cross_over_err'] > 2*max(d['within_over_err'],1e-9) else 'comparable'}")

    # ---- print --------------------------------------------------------------
    print(f"\n{'cell':12s} {'basis':5s} {'host':14s} {'thr':>3} | {'floor':>7} {'±':>5} "
          f"{'cnt(-e)':>8} {'count':>7} {'cnt(+e)':>8} | {'BAND':>10} {'pointwise':>10}")
    print("-" * 104)
    for ctag in [CH.cell_tag(k, lab, TIER) for (k, _e, lab) in CELLS]:
        for r in [x for x in rows if x["cell"] == ctag]:
            star = "~" if r["approx_hi"] else " "
            print(f"{r['cell']:12s} {r['basis']:5s} {r['host']:14s} {r['cpus']:3d} | "
                  f"{r['floor']:7.2f} {r['err']:5.2f} {r['count_hi']:7.2f}{star} "
                  f"{r['count']:7.2f} {r['count_lo']:8.2f} | "
                  f"{r['grade_band']:>10} {r['grade_pointwise']:>10}")
        print("-" * 104)
    print("  (~ = optimistic-side count linearised from the banked corr/corr_lo pair)")

    print(f"\n{'cell':12s} {'states':>6} {'floor band':>16} {'spread/err':>11} "
          f"{'count band':>14} {'ENSEMBLE GRADE':>18} {'pointwise':>22}")
    for s in summary:
        print(f"{s['cell']:12s} {s['n_states']:6d} "
              f"{s['floor_min']:7.2f}-{s['floor_max']:<8.2f} {s['spread_over_err']:11.2f} "
              f"{s['count_min']:6.2f}-{s['count_max']:<7.2f} {s['ensemble_grade']:>18} "
              f"{'/'.join(s['pointwise_labels']):>22}"
              + ("  <-- UNSTABLE" if s["pointwise_unstable"] else ""))

    os.makedirs(REPORTS, exist_ok=True)
    np.savez(f"{REPORTS}/lottery_recut_band.npz",
             rows_json=json.dumps(rows), summary_json=json.dumps(summary),
             decomp_json=json.dumps(decomp), ensemble_json=json.dumps(BASIS),
             rule=json.dumps(dict(
                 name="two-sided band rule (pre-registered)",
                 CONFIRMED="count(floor+err) > 1",
                 below="count(floor-err) <= 1",
                 MARGINAL="otherwise -- the >1 bar lies inside 1 sigma of the floor error",
                 ensemble="unanimous label, else MARGINAL (band) with [min,max] floor+count",
                 external_quote="unchanged: e=0.5 single / e=0.3 pair (conservative, band-safe)")))
    print(f"\nbanked -> {REPORTS}/lottery_recut_band.npz")
    return 0


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    rest = sys.argv[2:]
    sys.exit(mode_run(*rest) if mode == "run" else mode_analyze())
