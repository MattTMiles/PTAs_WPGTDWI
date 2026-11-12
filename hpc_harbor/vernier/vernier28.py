#!/usr/bin/env python3
"""
VERNIER-28 -- BLEND RECOVERY: the one lever left.
EXPLORATORY, NOT pre-registered. Toy, cronus, CPU/numpy, multiprocess.

Clean room: imports vernier28_lib (new, this instalment) and through it
vernier27_lib / vernier25 / vernier24 / vernier16 / vernier6 / vernier_toy.
None of those files is edited on disk. One imported CONSTANT is overridden at
run time and it is announced in the header: vernier27_lib.SUB_SPAN -> 0.0
(see vernier28_lib's "departure 1" -- handing R3 a sub-bin frequency the toy's
search cannot measure would flatter it).

WHY
  VERNIER-25 measured the conversion from AVAILABLE detections to RECOVERED
  combs at ~0.29-0.5: at gamma = 0.87 ten sources clear the floor, 2.9 are
  genuinely recovered, and 9.2 of twelve turns land on blends. The array
  forecasts (~4 resolvable MeerKAT-era, 30-40 SKA-era) count AVAILABLE sources;
  the certification bar is denominated in RECOVERED combs. Conversion is the
  bridge between the two, and it is the only methods lever left that can move
  the threshold rather than just measure it.

THE THREE LADDERS, on RNG-matched universes
  R1 BASELINE            sequential/greedy per-source recovery, exactly as
                         V24/V25/V27 run it. Control. Not improved.
  R2 JOINT-BIN           all sources in a bin fitted simultaneously, positions
                         re-solved against each other, M scanned 1..12 and
                         chosen against a null-calibrated floor.
  R3 JOINT-BIN + PULSAR  as R2 with the pulsar-term block as a second
                         observable and the pulsar distances marginalised over
                         their 10 pc priors. VERNIER-23 Part 1's degree-of-
                         freedom surplus, spent on recovery instead of on
                         certification.

THE CELLS  (K_available MEASURED, not asserted; A_equivalent labelled)
  C5    gamma 1.56 +  2 seeds     K_av ~  5
  C10   gamma 0.87 +  0 seeds     K_av ~ 11   -- V25's flat-slope cell
  C20   gamma 0.87 + 10 seeds     K_av ~ 19
  C10b  gamma 1.56 +  8 seeds     K_av ~ 11   -- construction ablation: same
                                    availability reached by seeding instead of
                                    by flattening, so conversion-vs-K_available
                                    is not confounded with gamma.

SAFETY KIT
  V18 genuine/split/blend accounting against the AVAILABLE pool with its chance
  floor printed; a blend counts as recovered only if it resolves into its true
  members. Pure-noise null universes through all three ladders (the phantom
  ledger). Source-level scrambled null on EVERY certification of every ladder.
  Wrong-fringe certifications reported per ladder. Null-calibrated model
  selection, floors recut by V24's own estimator.

Code hpc_harbor/vernier/vernier28.py + vernier28_lib.py . data reports/vernier28.npz
Staged, not committed.
"""
import os
import sys
import time
import numpy as np
import multiprocessing as mp

import vernier27_lib as L27

L27.SUB_SPAN = 0.0                     # announced override, see header

import vernier28_lib as Z              # noqa: E402
import vernier25 as V25                # noqa: E402
import vernier24 as W                  # noqa: E402

REPORTS = W.REPORTS

CELLS = [
    ("C5",   1.56,  2, "gamma 1.56 +  2 seeds"),
    ("C10",  0.87,  0, "gamma 0.87, unseeded"),
    ("C20",  0.87, 10, "gamma 0.87 + 10 seeds"),
    ("C10b", 1.56,  8, "gamma 1.56 +  8 seeds (construction ablation)"),
]
N_REAL = int(os.environ.get("V28_NREAL", 12))
N_PROC = int(os.environ.get("V28_NPROC", 6))
LADDER_I = {"R1": 1, "R2": 2, "R3": 3}


# ==================================================================== one run
def run_universe(task):
    ci, r = task
    name, gamma, K_seed, _ = CELLS[ci]
    t0 = time.time()
    arr = W.make_array(np.random.default_rng(2500 + r))
    seed = 28000 + 977 * r + 13 * ci
    ns = 880000 + 7919 * r + 101 * ci
    U = Z.universe(seed, arr, gamma, K_seed, ns)
    Un = Z.null_twin(U, arr)

    # --- null-calibrated model-selection floors, one per bin, shared by the
    #     signal arm and the null arm (they are properties of the bin, not of
    #     the data) and by R2 and R3 as far as the Earth block goes.
    rb = np.random.default_rng(ns + 310000)
    ctx = Z.BandCtx(arr, U)
    bars = [Z.null_curve_earth(rb, B) for B in U["bins"]]
    thrT = Z.null_dT(rb, U, arr, ctx)

    out = dict(cell=ci, real=r, K_av=U["K_av"],
               K_av_bin=[len(a) for a in U["avail"]])
    # A_equivalent label of the cell
    B0 = U["bins"][0]
    out["a_eq"] = (float(np.sqrt(B0["R2"][0] / B0["R2"][B0["n_seed"]]))
                   if B0["n_seed"] else 1.0)
    out["lam_fore"] = float(B0["diag"]["lam_fore"])
    out["A_eq"] = float(V25.a_equivalent(B0["diag"]["lam_fore"]))

    ladders = {}
    R1 = Z.ladder_R1(U, arr)
    R2 = Z.ladder_R2(U, arr, bars)
    R3 = Z.ladder_R3(U, arr, ctx, R2, thrT)
    ladders["R1"] = (R1["pos"], R1)
    ladders["R2"] = (R2["pos"], R2)
    ladders["R3"] = (R3["pos"], R3)
    out["R3_drop"], out["R3_add"] = R3["n_drop"], R3["n_add"]

    # the null arm reuses ctx: its comb kernels are data-independent, and the
    # null twin shares the same sky, so nothing about them changes.
    R1n = Z.ladder_R1(Un, arr)
    R2n = Z.ladder_R2(Un, arr, bars)
    R3n = Z.ladder_R3(Un, arr, ctx, R2n, thrT)
    nulls = {"R1": R1n["pos"], "R2": R2n["pos"], "R3": R3n["pos"]}

    for lab, (pos, raw) in ladders.items():
        n, per_bin, kinds = Z.score(U, pos)
        crng = np.random.default_rng(ns + 610000 + 1009 * LADDER_I[lab])
        rec = Z.rec_order(U, per_bin)
        hist = Z.certify_pass(crng, arr, U, rec)
        es = Z.end_state(hist)
        nn, _, _ = Z.score(Un, nulls[lab])
        out[lab] = dict(
            claims=n["claims"], genuine=n["genuine"], split=n["split"],
            faint=n["faint"], phantom=n["phantom"], blend=n["blend"],
            rec=len(rec), **{f"es_{k}": v for k, v in es.items()},
            null_claims=nn["claims"], null_phantom=nn["phantom"],
            null_genuine=nn["genuine"],
            Mhat=[int(raw["Mhat"].get(b, 0)) for b in range(len(U["bins"]))]
            if "Mhat" in raw else
            [sum(1 for h in raw["hist"] if h["bin"] == b)
             for b in range(len(U["bins"]))],
            gen_bin=[len(per_bin.get(b, [])) for b in range(len(U["bins"]))],
        )
    out["R1"]["native_genuine"] = R1["native_genuine"]
    out["R1"]["native_blend"] = R1["native_blend"]
    out["R1"]["ignite"] = float(R1["ignite"])
    out["R1"]["turns"] = R1["turns"]
    out["R1"]["int_cert"] = R1["n_cert_tot"]
    out["R1"]["int_wrong"] = R1["n_wrong_tot"]
    out["R1"]["int_null"] = R1["n_null_tot"]
    out["secs"] = time.time() - t0
    return out


# ============================================================ R1 == the control
def validate_R1(n=3):
    """
    The control has to BE the control: assert that vernier28_lib's instrumented
    copy reproduces vernier27_lib.chain turn for turn on the same universe.
    """
    print("=" * 110)
    print("CONTROL VALIDATION -- instrumented R1 vs the imported vernier27_lib.chain")
    print("=" * 110)
    ok = True
    for r in range(n):
        arr = W.make_array(np.random.default_rng(2500 + r))
        seed, ns = 28000 + 977 * r, 880000 + 7919 * r
        U = Z.universe(seed, arr, 1.56, 8, ns, mc_fix=None)
        mine = Z.ladder_R1(U, arr)
        theirs = L27.chain(seed, arr, 1.56, n_seed=8, max_turn=Z.MAX_TURN,
                           comb_cap=Z.COMB_CAP, noise_seed=ns)
        same = (mine["turns"] == theirs["turns"]
                and mine["native_genuine"] == theirs["genuine"]
                and mine["native_blend"] == theirs["blend"]
                and mine["n_cert_tot"] == theirs["n_cert_tot"])
        ok &= same
        print(f"  real {r}: turns {mine['turns']}/{theirs['turns']}  "
              f"genuine {mine['native_genuine']}/{theirs['genuine']}  "
              f"blend {mine['native_blend']}/{theirs['blend']}  "
              f"cert {mine['n_cert_tot']}/{theirs['n_cert_tot']}   "
              f"{'MATCH' if same else 'DIFFERS'}")
    print(f"  -> control {'reproduced exactly' if ok else 'DOES NOT MATCH -- stop'}")
    print()
    return ok


# ================================================================ deliverables
def agg(res, ci, lab, key):
    return np.array([q[lab][key] for q in res if q["cell"] == ci], float)


def d1_conversion(out, res):
    print()
    print("=" * 118)
    print("DELIVERABLE 1 -- CONVERSION EFFICIENCY  K_recovered / K_available")
    print("=" * 118)
    print("  GENUINE = a fitted component inside 0.15 rad of an AVAILABLE true source")
    print("  that no other component of the same fit has claimed. SPLIT = a second")
    print("  component on an already-claimed source. FAINT = a real source that was")
    print("  not available. PHANTOM = nothing within tolerance. No credit is given for")
    print("  fitting one component to two truths.")
    print()
    hdr = (f"  {'cell':>5} {'construction':>34} {'K_av':>6} {'A_eq/A15':>9} "
           f"{'ladder':>7} {'claims':>7} {'GEN':>6} {'split':>6} {'faint':>6} "
           f"{'phantom':>8} {'K_rec/K_av':>11} {'chance':>7}")
    print(hdr)
    rows = []
    for ci, (name, gamma, K_seed, desc) in enumerate(CELLS):
        sub = [q for q in res if q["cell"] == ci]
        if not sub:
            continue
        kav = np.array([q["K_av"] for q in sub], float)
        aeq = np.median([q["A_eq"] for q in sub]) / W.A_GWB
        cf = float(np.mean([Z.chance_floor(q["K_av"] / len(q["K_av_bin"]))
                            for q in sub]))
        for lab in ("R1", "R2", "R3"):
            g = agg(res, ci, lab, "genuine")
            c = agg(res, ci, lab, "claims")
            eff = g / np.maximum(kav, 1e-9)
            rows.append([ci, K_seed, float(np.median(kav)), aeq,
                         {"R1": 1, "R2": 2, "R3": 3}[lab],
                         float(np.mean(c)), float(np.mean(g)),
                         float(np.mean(agg(res, ci, lab, "split"))),
                         float(np.mean(agg(res, ci, lab, "faint"))),
                         float(np.mean(agg(res, ci, lab, "phantom"))),
                         float(np.mean(eff)), float(np.std(eff, ddof=1)
                                                    / max(len(eff) ** 0.5, 1)),
                         cf])
            print(f"  {name if lab == 'R1' else '':>5} "
                  f"{desc if lab == 'R1' else '':>34} "
                  f"{np.median(kav) if lab == 'R1' else 0:6.1f} "
                  f"{aeq if lab == 'R1' else 0:9.2f} "
                  f"{lab:>7} {np.mean(c):7.2f} {np.mean(g):6.2f} "
                  f"{np.mean(agg(res, ci, lab, 'split')):6.2f} "
                  f"{np.mean(agg(res, ci, lab, 'faint')):6.2f} "
                  f"{np.mean(agg(res, ci, lab, 'phantom')):8.2f} "
                  f"{np.mean(eff):8.3f}+-{np.std(eff, ddof=1)/max(len(eff)**0.5,1):.3f} "
                  f"{cf:7.3f}")
        print()
    out["conv"] = np.array(rows)
    return np.array(rows)


def d2_separability(out, res):
    print()
    print("=" * 118)
    print("DELIVERABLE 2 -- REALIZED PER-BIN SEPARABILITY vs VERNIER-24's 5-10")
    print("=" * 118)
    print("  M_hat  = components the data supports in a bin (null-calibrated).")
    print("  G_bin  = of those, how many are GENUINE. The counting bound for 116")
    print(f"  pulsars is 3N/7 = {3*116/7:.0f}; V24's REALIZED frontier was 5-10.")
    print()
    print(f"  {'cell':>5} {'ladder':>7} {'M_hat med':>10} {'M_hat max':>10} "
          f"{'G_bin med':>10} {'G_bin max':>10} {'M_hat|G=M':>10} "
          f"{'bins with G>=2':>15}")
    rows = []
    for ci, (name, *_rest) in enumerate(CELLS):
        sub = [q for q in res if q["cell"] == ci]
        if not sub:
            continue
        for lab in ("R1", "R2", "R3"):
            Mh = np.array([m for q in sub for m in q[lab]["Mhat"]], float)
            Gb = np.array([g for q in sub for g in q[lab]["gen_bin"]], float)
            clean = Mh[(Mh > 0) & (Gb == Mh)]
            rows.append([ci, {"R1": 1, "R2": 2, "R3": 3}[lab],
                         float(np.median(Mh[Mh > 0]) if (Mh > 0).any() else 0),
                         float(Mh.max()), float(np.median(Gb[Gb > 0])
                                                if (Gb > 0).any() else 0),
                         float(Gb.max()),
                         float(clean.max() if clean.size else 0),
                         float(np.mean(Gb >= 2))])
            print(f"  {name if lab == 'R1' else '':>5} {lab:>7} "
                  f"{rows[-1][2]:10.1f} {rows[-1][3]:10.0f} "
                  f"{rows[-1][4]:10.1f} {rows[-1][5]:10.0f} "
                  f"{rows[-1][6]:10.0f} {rows[-1][7]:15.2f}")
        print()
    out["sep"] = np.array(rows)
    return np.array(rows)


def d3_kstar(out, res, conv):
    print()
    print("=" * 118)
    print("DELIVERABLE 3 -- THE CALENDAR:  K* RE-EXPRESSED IN AVAILABLE-SOURCE UNITS")
    print("=" * 118)
    print("  The invariant is the number of ACCUMULATED COMBS the fringe contrast")
    print("  needs to clear its trials bar. K* in available units is that number")
    print("  divided by the ladder's conversion efficiency.")
    print()
    print(f"  {'cell':>5} {'ladder':>7} {'combs':>6} {'S_end':>8} {'bar':>6} "
          f"{'S/bar':>7} {'cross':>6} {'CERT':>6} {'WRONG':>6} {'null':>6}")
    pts = []
    for ci, (name, *_rest) in enumerate(CELLS):
        sub = [q for q in res if q["cell"] == ci]
        if not sub:
            continue
        for lab in ("R1", "R2", "R3"):
            nc = agg(res, ci, lab, "es_combs")
            S = agg(res, ci, lab, "es_S_end")
            bar = agg(res, ci, lab, "es_bar")
            cr = agg(res, ci, lab, "es_cross")
            ce = agg(res, ci, lab, "es_cert")
            wr = agg(res, ci, lab, "es_wrong")
            nu = agg(res, ci, lab, "es_null")
            bm = float(np.nanmedian(bar))
            print(f"  {name if lab == 'R1' else '':>5} {lab:>7} "
                  f"{np.median(nc):6.1f} {np.median(S):8.3f} {bm:6.2f} "
                  f"{np.median(S)/max(bm,1e-9):7.3f} {np.mean(cr):6.2f} "
                  f"{np.mean(ce):6.2f} {np.mean(wr):6.2f} {np.mean(nu):6.2f}")
            pts.append([np.median(nc), np.median(S), bm])
        print()
    pts = np.array(pts)
    out["contrast"] = pts

    # --- combs* : where does the contrast curve meet the bar
    m = pts[:, 0] > 0
    combs_star = np.nan
    if m.sum() >= 3:
        sl, ic = np.polyfit(pts[m, 0], pts[m, 1], 1)
        bar = float(np.nanmedian(pts[m, 2]))
        if sl > 0:
            combs_star = (bar - ic) / sl
        print(f"  contrast vs combs (pooled, linear): S = {sl:+.3f} x combs "
              f"{ic:+.3f} nat;  bar {bar:.2f} nat")
        print(f"  -> combs* (contrast reaches the bar) = {combs_star:.1f} "
              f"accumulated combs")
        out["fit"] = np.array([sl, ic, bar, combs_star])
    ref = 10.5           # V25's extrapolated combs*, used when this run's own
    use = combs_star if np.isfinite(combs_star) else ref
    print()
    print(f"  combs* used below: {use:.1f} "
          f"({'measured here' if np.isfinite(combs_star) else 'V25 reference, this run did not resolve one'})")
    print(f"  {'ladder':>7} {'conversion (pooled)':>21} {'K*_available':>14} "
          f"{'K* at V25 combs*=10.5':>23} {'era':>30}")
    rows = []
    for lab, li in (("R1", 1), ("R2", 2), ("R3", 3)):
        eff = float(np.mean(conv[conv[:, 4] == li, 10]))
        ka = use / max(eff, 1e-9)
        era = ("MeerKAT-era (<= 4 available)" if ka <= 4 else
               "SKA 30 yr (30-40 available)" if ka <= 40 else
               "beyond the SKA 30 yr forecast")
        rows.append([li, eff, ka, ref / max(eff, 1e-9)])
        print(f"  {lab:>7} {eff:21.3f} {ka:14.1f} "
              f"{ref/max(eff,1e-9):23.1f} {era:>30}")
    out["kstar"] = np.array(rows)
    print()
    print("  NOTE on V25's number. V25 quoted K* ~ 20-25 as SIMULTANEOUS DETECTIONS,")
    print("  reached by dividing ~10-11 needed combs by a conversion of ~0.5 -- so")
    print("  20-25 was already in AVAILABLE/SEEDED units, not recovered ones.")
    print("  Re-expressing it a second time through the conversion would double-count.")
    print("  The invariant quoted above is combs*, and K*_available = combs*/conversion.")
    return rows


def d4_phantoms(out, res):
    print()
    print("=" * 118)
    print("DELIVERABLE 4 -- THE WRONG-RECOVERY LEDGER")
    print("=" * 118)
    print("  PHANTOM = fitted component matching no true source at all.")
    print("  FAINT   = matches a real source that was not AVAILABLE. At 300 sources a")
    print(f"            bin, a random sky point lands within 0.15 rad of one of them")
    print(f"            {Z.chance_floor(300):.3f} of the time, so the faint column is")
    print("            mostly chance coincidence and is NOT a recovery.")
    print("  WRONG   = the V18 blend rate, (faint + phantom)/claims: the fraction of")
    print("            everything a ladder banks that is not a clean recovery.")
    print("  Null arm: the SAME sky, same noise draws, signal removed from BOTH blocks")
    print("  -- every component it claims is false. Certification columns carry the")
    print("  source-level scrambled null, run at every step of every ladder.")
    print()
    print(f"  {'cell':>5} {'ladder':>7} {'phantom/run':>12} {'phantom/claim':>14} "
          f"{'WRONG/claim':>12} {'NULL claims':>12} {'NULL phantom':>13} "
          f"{'wrong cert':>11} {'scrambled cert':>15}")
    rows = []
    for ci, (name, *_rest) in enumerate(CELLS):
        sub = [q for q in res if q["cell"] == ci]
        if not sub:
            continue
        for lab in ("R1", "R2", "R3"):
            ph = agg(res, ci, lab, "phantom")
            fa = agg(res, ci, lab, "faint")
            cl = agg(res, ci, lab, "claims")
            nc = agg(res, ci, lab, "null_claims")
            np_ = agg(res, ci, lab, "null_phantom")
            wr = agg(res, ci, lab, "es_wrong")
            nu = agg(res, ci, lab, "es_null")
            bl = float(np.sum(ph + fa) / max(np.sum(cl), 1e-9))
            rows.append([ci, {"R1": 1, "R2": 2, "R3": 3}[lab],
                         float(np.mean(ph)),
                         float(np.sum(ph) / max(np.sum(cl), 1e-9)), bl,
                         float(np.mean(nc)), float(np.mean(np_)),
                         float(np.mean(wr)), float(np.mean(nu))])
            print(f"  {name if lab == 'R1' else '':>5} {lab:>7} "
                  f"{np.mean(ph):12.2f} "
                  f"{np.sum(ph)/max(np.sum(cl),1e-9):14.3f} {bl:12.3f} "
                  f"{np.mean(nc):12.2f} {np.mean(np_):13.2f} "
                  f"{np.mean(wr):11.2f} {np.mean(nu):15.2f}")
        print()
    out["ledger"] = np.array(rows)
    return rows


def d5_fisher(out, n_real=3):
    """
    The Fisher companion: V24's separability(), which is where the dof surplus
    is COUNTED, run beside the ladders that try to SPEND it. Unmodified call.
    """
    print()
    print("=" * 118)
    print("DELIVERABLE 2b -- THE FISHER COMPANION (V24.separability, unmodified)")
    print("=" * 118)
    print("  resolved = sigma_sky inside half the nearest-neighbour separation AND")
    print("  sigma(A)/A < 0.5, from the INVERSE marginalised Fisher.")
    print()
    print("  Called at V24's identifiability-limited working point, rho_1 = 30 fixed,")
    print("  which is the column V24's 'realized frontier 5-10' was read from.")
    print()
    print(f"  {'M':>4} {'Earth only':>22} {'+ pulsar, marginalised':>26}")
    print(f"  {'':>4} {'resolved  top_ok':>22} {'resolved  top_ok':>26}")
    rows = []
    rng = np.random.default_rng(28061)
    for M in (2, 5, 10, 17, 25, 40):
        a = [W.separability(rng, W.N_PSR, M, 15e-9, use_pterm=False,
                            rho_fixed=30.0) for _ in range(n_real)]
        b = [W.separability(rng, W.N_PSR, M, 15e-9, use_pterm=True,
                            marginalise=True, rho_fixed=30.0)
             for _ in range(n_real)]
        ra, rb = np.median([x["resolved"] for x in a]), \
            np.median([x["resolved"] for x in b])
        ta, tb = np.mean([x["top_ok"] for x in a]), \
            np.mean([x["top_ok"] for x in b])
        rows.append([M, ra, ta, rb, tb])
        print(f"  {M:4d} {ra:12.1f} {ta:9.2f} {rb:16.1f} {tb:9.2f}")
    out["fisher"] = np.array(rows)


# ======================================================================= main
def main():
    print()
    print("#" * 118)
    print("# VERNIER-28 -- BLEND RECOVERY: THE ONE LEVER LEFT")
    print("# EXPLORATORY, NOT PRE-REGISTERED")
    print("#" * 118)
    print(f"  cells {[c[0] for c in CELLS]}, {N_REAL} realisations each, "
          f"{N_PROC} processes")
    print(f"  M_MAX {Z.M_MAX}  SAGE sweeps {Z.N_SWEEP}  candidates {Z.N_CAND}  "
          f"L-grid {Z.NPTS_L} (even)  model-selection nulls {Z.NJ_NULL}")
    print(f"  COMB_CAP {Z.COMB_CAP}  max_turn {Z.MAX_TURN}  "
          f"vernier27_lib.SUB_SPAN overridden to {L27.SUB_SPAN}")
    print()
    if not validate_R1():
        sys.exit(1)

    tasks = [(ci, r) for ci in range(len(CELLS)) for r in range(N_REAL)]
    t0 = time.time()
    if N_PROC > 1:
        with mp.Pool(N_PROC) as pool:
            res = []
            for q in pool.imap_unordered(run_universe, tasks):
                res.append(q)
                print(f"    [{CELLS[q['cell']][0]:>5} r{q['real']:02d}] "
                      f"K_av {q['K_av']:2d}  gen R1/R2/R3 "
                      f"{q['R1']['genuine']:2d}/{q['R2']['genuine']:2d}/"
                      f"{q['R3']['genuine']:2d}  "
                      f"phantom {q['R1']['phantom']:2d}/{q['R2']['phantom']:2d}/"
                      f"{q['R3']['phantom']:2d}  "
                      f"S_end {q['R1']['es_S_end']:.2f}/"
                      f"{q['R2']['es_S_end']:.2f}/{q['R3']['es_S_end']:.2f}  "
                      f"({q['secs']:.0f}s, {len(res)}/{len(tasks)}, "
                      f"{time.time()-t0:.0f}s elapsed)", flush=True)
    else:
        res = [run_universe(t) for t in tasks]
    print(f"\n  all universes done in {time.time()-t0:.0f} s")

    out = {}
    conv = d1_conversion(out, res)
    d2_separability(out, res)
    d5_fisher(out)
    d3_kstar(out, res, conv)
    d4_phantoms(out, res)
    out["n_real"] = np.array([N_REAL])
    out["cells"] = np.array([c[0] for c in CELLS])
    np.savez(os.path.join(REPORTS, "vernier28.npz"), **out)
    print(f"\n  banked -> {os.path.join(REPORTS, 'vernier28.npz')}")


if __name__ == "__main__":
    main()
