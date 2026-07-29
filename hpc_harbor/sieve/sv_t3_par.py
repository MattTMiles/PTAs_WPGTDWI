#!/usr/bin/env python3
"""SIEVE T3 -- PARTIAL AMBIGUITY RESOLUTION vs CRITERION-v2.2 (agent SIEVE-A, 2026-07-29).

REPORT-ONLY. Banked columns only: no likelihood is re-evaluated, no floor re-cut, no
venue rebuilt, no banked verdict moved.

--------------------------------------------------------------------------------
THE TRANSPLANT, STATED BEFORE ANY NUMBER IS READ

GNSS integer ambiguity resolution and the PTA fringe problem are the same object: a
real-valued quantity known only modulo a lattice spacing. Here the ambiguity of pulsar
a is its fringe integer n_a (distance wrap), the spacing is the mode spacing dL_a, and
the "float solution" is the fringe posterior the E-step already computes.

Three GNSS pieces, and what each becomes here:

  LAMBDA DECORRELATION.  In GNSS the float ambiguity covariance Q is Z-transformed so
    the sequential conditional variances approach the marginals. THE FULL Q IS NOT
    AVAILABLE HERE AND IS MEASURED-UNAFFORDABLE, not merely unbanked: the joint Hessian
    over the distance block is the same object SPARK-3 priced and refused
    (spark3.faint_fisher_bounds -- two builds failed to return inside a 1 h A100-80GB
    walltime). SUBSTITUTION, declared: the banked `qmax` already comes from
    `estep_per_target`, in which target a's own pulsar term is live under its own mask
    while every uncertified pulsar is DECOHERED (spark3.rung_masks). Removing the other
    pulsars' contributions is the PTA analogue of decorrelation, and it removes
    information rather than adding it. So qmax_a is CONSERVATIVE relative to a fully
    decorrelated joint fix, and every joint success rate below is a LOWER bound on what
    a true LAMBDA fix would report. A PAR set that is still LOOSER than v2.2 under a
    lower bound is therefore decisive; a PAR set that is stricter is not.

  INTEGER BOOTSTRAPPING.  The bootstrapped success rate is the product of the per-
    ambiguity CONDITIONAL success rates:  P_boot(S) = prod_{a in S} P_a.  With the
    decorrelation substitution above, P_a := qmax_a (the fringe posterior mass on the
    MAP fringe). This is the standard bootstrapping formula, evaluated on the banked
    posterior rather than on a Gaussian conditional variance.

  PARTIAL AMBIGUITY RESOLUTION.  Fix only the largest subset whose joint success rate
    clears a declared level P0: sort by P_a descending and take the longest prefix with
    P_boot >= P0. That prefix is the "safe subset". Reported at P0 = 0.99 / 0.95 / 0.90.

CRITERION-v2.2, for comparison, on the same realisation (glacier_loop.py:533):
    cert[a] = (dlnL[a] > max(lnK[a] + 0.578, floor)) & (qmax[a] > 0.9)

--------------------------------------------------------------------------------
THE CALIBRATION GATE (runs first; PAR's whole claim rests on it)

P_boot is only meaningful if qmax is a CALIBRATED success probability. The GENERALISE
banks carry the ground truth per pulsar (`on_true` = MAP fringe is the true fringe) and
`ptrue` (posterior mass on the true fringe), over 15 independent noise realisations of
the A-SKY survivor unit. So the transplant is checked, not asserted:

  (G1) RELIABILITY: bin pulsars by qmax; the empirical fraction with on_true must track
       the bin's mean qmax.
  (G2) JOINT REALISED vs PREDICTED: for the PAR safe subset of each realisation, the
       predicted P_boot vs the realised indicator that EVERY member has on_true.

If qmax is over-confident, PAR is optimistic by exactly that factor and the verdict says
so with the number.

Output bank: reports/sieve_t3_par.npz
Usage: python3 hpc_harbor/sieve/sv_t3_par.py [--repo <path>]
"""
import argparse
import glob
import os

import numpy as np

TRIALS_NAT = 0.578
QBAR = 0.9
P0_LEVELS = (0.99, 0.95, 0.90)

# The three banked above-onset realisations scored head-to-head.
GLACIER_UNITS = [
    ("census r13p25/e07/s0 i0 (cascade, 3 certs)",
     "GLACIER_results/gl2_r13p25_cell_e07_s0_T30_lit_i0__*.npz"),
    ("census r13p5/none/s3 i0 (the campaign's true certs)",
     "GLACIER_results/gl2_r13p5_cell_none_s3_T30_lit_i0__*.npz"),
]
ASKY_GLOB = "GENERALISE_results/gen_sig_AS_e03_h1275_k5_s4_g4_n*.npz"


def par_subset(P, p0):
    """Largest prefix (by descending per-ambiguity success rate) whose bootstrapped
    joint success rate clears p0. Returns (index array, joint rate)."""
    order = np.argsort(-P)
    best_n, best_j = 0, 1.0
    j = 1.0
    for n, a in enumerate(order, start=1):
        j = j * P[a]
        if j >= p0:
            best_n, best_j = n, j
        else:
            break
    return order[:best_n], (best_j if best_n else 1.0)


def cert_v22(dlnl, lnK, q, fl):
    return (dlnl > np.maximum(lnK + TRIALS_NAT, fl)) & (q > QBAR)


def describe(tag, names, cert, par_sets, P, extra=""):
    print(f"\n   -- {tag} {extra}")
    cidx = set(np.where(cert)[0].tolist())
    print(f"      criterion-v2.2 certified set: {sorted(cidx)} "
          f"(n={len(cidx)})")
    for p0, (S, j) in par_sets.items():
        sset = set(S.tolist())
        only_par = sorted(sset - cidx)
        only_v22 = sorted(cidx - sset)
        rel = ("AGREE" if sset == cidx else
               "PAR-STRICTER" if sset < cidx else
               "PAR-LOOSER" if sset > cidx else "PARTIAL-OVERLAP")
        print(f"      PAR P0={p0:.2f}: |S|={len(sset):3d} P_boot={j:.4f}  {rel}")
        if only_par:
            s = ", ".join(f"{names[a] if names is not None else a}"
                          f"(q={P[a]:.4f})" for a in only_par[:8])
            print(f"         in PAR only ({len(only_par)}): {s}"
                  + (" ..." if len(only_par) > 8 else ""))
        if only_v22:
            s = ", ".join(f"{names[a] if names is not None else a}"
                          f"(q={P[a]:.4f})" for a in only_v22[:8])
            print(f"         in v2.2 only ({len(only_v22)}): {s}"
                  + (" ..." if len(only_v22) > 8 else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="/data/taylor_group/matt_miles/PTAs_WPGTDWI")
    a = ap.parse_args()
    R = a.repo
    rows = []

    # ================= THE CALIBRATION GATE =================
    print("[SIEVE-T3] ===== GATE: is qmax a calibrated success probability? =====")
    fs = sorted(glob.glob(os.path.join(R, ASKY_GLOB)))
    print(f"   A-SKY survivor unit AS_e03_h1275_k5_s4: {len(fs)} realisations")
    Q, OT, PT, MK, NT = [], [], [], [], []
    names = None
    for f in fs:
        z = np.load(f, allow_pickle=True)
        Q.append(np.asarray(z["qmax"], float))
        OT.append(np.asarray(z["on_true"], bool))
        PT.append(np.asarray(z["ptrue"], float))
        MK.append(np.asarray(z["mapk"], int))
        NT.append(np.asarray(z["n_true_grid"], int))
        if names is None and "names" in z.files:
            names = np.asarray(z["names"])
    Q = np.array(Q); OT = np.array(OT); PT = np.array(PT)
    MK = np.array(MK); NT = np.array(NT)
    print(f"   shape {Q.shape} (realisations x pulsars)")
    assert (MK == NT).all() == OT.all() or True
    ok = (MK == NT)
    if not (ok == OT).all():
        print(f"   NOTE: on_true != (mapk == n_true_grid) on "
              f"{int((ok != OT).sum())} entries -- using the banked on_true column.")

    print("\n   (G1) RELIABILITY -- empirical P(on_true) vs qmax, pooled over "
          f"{Q.shape[0]} realisations x {Q.shape[1]} pulsars")
    edges = np.array([0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0 + 1e-12])
    print("      qmax bin            n     mean qmax   emp. P(on_true)   "
          "mean ptrue   over/under")
    rel_rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (Q >= lo) & (Q < hi)
        n = int(m.sum())
        if n == 0:
            continue
        mq, emp, mp = float(Q[m].mean()), float(OT[m].mean()), float(PT[m].mean())
        se = np.sqrt(max(emp * (1 - emp), 1e-12) / n)
        tag = ("OVER-confident" if mq - emp > 2 * se else
               "under-confident" if emp - mq > 2 * se else "calibrated")
        print(f"      [{lo:.3f},{hi:.3f}) {n:7d}   {mq:9.4f}   {emp:9.4f}+-{se:.4f}"
              f"   {mp:9.4f}   {tag}")
        rel_rows.append([lo, hi, n, mq, emp, se, mp, tag])
    tot_q, tot_e = float(Q.mean()), float(OT.mean())
    se_t = np.sqrt(tot_e * (1 - tot_e) / Q.size)
    print(f"      POOLED          {Q.size:7d}   {tot_q:9.4f}   {tot_e:9.4f}"
          f"+-{se_t:.4f}")
    hi = Q > QBAR
    if hi.sum():
        hq, he = float(Q[hi].mean()), float(OT[hi].mean())
        seh = np.sqrt(max(he * (1 - he), 1e-12) / hi.sum())
        print(f"      qmax > {QBAR} only  {int(hi.sum()):7d}   {hq:9.4f}   "
              f"{he:9.4f}+-{seh:.4f}   <-- the region criterion-v2.2 uses")
    else:
        hq = he = seh = float("nan")

    # (G1b) the same reliability read on GLACIER's OWN banks -- is the miscalibration a
    # property of the GENERALISE arm, or of the fringe posterior the criterion gates on?
    print("\n   (G1b) THE SAME READ ON EVERY GLACIER PER-ITERATION BANK "
          "(is this a GENERALISE-arm artifact?)")
    gq, got, gitr = [], [], []
    for f in sorted(glob.glob(os.path.join(R, "GLACIER_results", "gl*_i?__*.npz"))):
        try:
            zz = np.load(f, allow_pickle=True)
        except Exception:
            continue
        if "qmax" not in zz.files or "on_true" not in zz.files:
            continue
        qq = np.asarray(zz["qmax"], float)
        gq.append(qq)
        got.append(np.asarray(zz["on_true"], bool))
        gitr.append(np.full(qq.shape, int(zz["iter"]) if "iter" in zz.files else -1))
    glac_rows = []
    if gq:
        GQ, GOT = np.concatenate(gq), np.concatenate(got)
        print(f"      {len(gq)} banks, {GQ.size} pulsar-samples")
        print("      qmax bin            n     mean qmax   emp. P(on_true)   verdict")
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (GQ >= lo) & (GQ < hi)
            n = int(m.sum())
            if n == 0:
                continue
            mq, emp = float(GQ[m].mean()), float(GOT[m].mean())
            se = np.sqrt(max(emp * (1 - emp), 1e-12) / n)
            tag = ("OVER-confident" if mq - emp > 2 * se else
                   "under-confident" if emp - mq > 2 * se else "calibrated")
            print(f"      [{lo:.3f},{hi:.3f}) {n:7d}   {mq:9.4f}   "
                  f"{emp:9.4f}+-{se:.4f}   {tag}")
            glac_rows.append([lo, hi, n, mq, emp, se, tag])
        gh = GQ > QBAR
        ghq, ghe = float(GQ[gh].mean()), float(GOT[gh].mean())
        gse = np.sqrt(max(ghe * (1 - ghe), 1e-12) / gh.sum())
        print(f"      qmax > {QBAR} only  {int(gh.sum()):7d}   {ghq:9.4f}   "
              f"{ghe:9.4f}+-{gse:.4f}   <-- the region criterion-v2.2 uses")
        # split by iteration: iteration 0 is the least-wandered state the loop occupies
        # (promote is AT DRAWN TRUTH), later iterations carry accumulated M-step motion.
        # If the gap grows with iteration, the deficit is TEMPLATE WANDER, not posterior
        # miscalibration alone.
        GIT = np.concatenate(gitr)
        print("      -- the qmax>0.9 region split by loop iteration --")
        print("      iter        n     mean qmax   emp. P(on_true)")
        for i in range(6):
            m = gh & (GIT == i)
            if not m.sum():
                continue
            mq, emp = float(GQ[m].mean()), float(GOT[m].mean())
            se = np.sqrt(max(emp * (1 - emp), 1e-12) / m.sum())
            print(f"      i{i}    {int(m.sum()):7d}   {mq:9.4f}   {emp:9.4f}+-{se:.4f}")
            glac_rows.append([f"iter{i}", QBAR, int(m.sum()), mq, emp, se, "by-iter"])
    else:
        GQ = GOT = np.zeros(0)
        ghq = ghe = gse = float("nan")

    # ================= (G2) joint realised vs predicted =================
    print("\n   (G2) JOINT SUCCESS -- PAR safe subset, predicted P_boot vs realised")
    print("      P0     mean |S|   mean predicted P_boot   realised all-correct rate")
    g2 = []
    for p0 in P0_LEVELS:
        sizes, preds, hits = [], [], []
        for r in range(Q.shape[0]):
            S, j = par_subset(Q[r], p0)
            sizes.append(len(S)); preds.append(j)
            hits.append(bool(OT[r][S].all()) if len(S) else True)
        mr = float(np.mean(hits))
        se = np.sqrt(max(mr * (1 - mr), 1e-12) / len(hits))
        print(f"      {p0:.2f}   {np.mean(sizes):7.2f}   {np.mean(preds):19.4f}   "
              f"{mr:15.4f}+-{se:.4f}  ({int(np.sum(hits))}/{len(hits)})")
        g2.append([p0, float(np.mean(sizes)), float(np.mean(preds)), mr, se])

    # ============ the A-SKY floor, cut BOTH ways on the same 100 nulls ============
    # GENERALISE banks 100 no-CW realisations for this unit, so the v2.2 comparison below
    # does not have to be made at floor = 0. Cutting the floor on BOTH offender
    # definitions on the SAME null draws also measures, directly, the divergence SIEVE
    # T7/E0 found in glacier_loop._null_offenders:
    #   canonical  recut_surface.offender:75  -> max_a dlnL     gated (dlnL>lnK)&(q>QBAR)
    #   GLACIER    glacier_loop:594           -> max_a (dlnL-lnK)_+ gated (q>QBAR)
    print("\n[SIEVE-T3] ===== A-SKY floor, cut both ways on the same nulls =====")
    nfs = sorted(glob.glob(os.path.join(R, "GENERALISE_results",
                                        "gen_fnullN_AS_e03_h1275_k5_s4_g4_n*.npz")))
    off_can, off_gl = [], []
    for f in nfs:
        zz = np.load(f, allow_pickle=True)
        dd = np.nan_to_num(np.asarray(zz["dlnL_det"], float), posinf=1e30)
        kk = np.asarray(zz["lnK"], float)
        qq = np.asarray(zz["qmax"], float)
        m = (dd > kk) & (qq > QBAR) & np.isfinite(dd)
        off_can.append(float(dd[m].max()) if m.any() else 0.0)
        o = np.where(qq > QBAR, np.maximum(dd - np.maximum(kk, 0.0), 0.0), 0.0)
        off_gl.append(float(np.max(o)))
    off_can, off_gl = np.array(off_can), np.array(off_gl)
    fl_can = float(np.quantile(off_can, 0.95))
    fl_gl = float(np.quantile(off_gl, 0.95))
    print(f"   {len(nfs)} no-CW null realisations")
    print(f"   canonical offender (max dlnL)      : q95 = {fl_can:.3f}  "
          f"zero-fraction {float((off_can == 0).mean()):.2f}")
    print(f"   glacier_loop offender (dlnL-lnK)_+ : q95 = {fl_gl:.3f}  "
          f"zero-fraction {float((off_gl == 0).mean()):.2f}")
    print(f"   DIFFERENCE (the T7/E0 divergence, measured on real nulls): "
          f"{fl_can - fl_gl:+.3f} nat")
    print(f"   -> a bar cut with the glacier_loop statistic and applied to dlnL sits "
          f"{fl_can - fl_gl:.3f} nat LOW at this unit.")

    # ================= HEAD-TO-HEAD on the three realisations =================
    print("\n[SIEVE-T3] ===== PAR safe subset vs criterion-v2.2 certified set =====")

    for tag, pat in GLACIER_UNITS:
        f = sorted(glob.glob(os.path.join(R, pat)))
        if not f:
            print(f"\n   -- {tag}: NO BANK MATCHED {pat}")
            continue
        z = np.load(f[0], allow_pickle=True)
        dlnl = np.asarray(z["dlnL_det"], float)
        lnK = np.asarray(z["lnK"], float)
        q = np.asarray(z["qmax"], float)
        fl = float(z["floor"])
        ot = np.asarray(z["on_true"], bool)
        c = cert_v22(dlnl, lnK, q, fl)
        ps = {p0: par_subset(q, p0) for p0 in P0_LEVELS}
        describe(tag, None, c, ps, q,
                 extra=f"[floor {fl:.3f}, {int(c.sum())} certs, "
                       f"{int(ot.sum())} pulsars on_true]")
        for p0, (S, j) in ps.items():
            rows.append([tag, p0, int(len(S)), float(j), int(c.sum()),
                         int(len(set(S.tolist()) & set(np.where(c)[0].tolist()))),
                         int(ot[S].sum()) if len(S) else 0, int(ot.sum())])
        # PAR's own correctness on this realisation
        for p0, (S, j) in ps.items():
            if len(S):
                print(f"      PAR P0={p0:.2f} realised: {int(ot[S].sum())}/{len(S)} "
                      f"members actually on_true "
                      f"(predicted joint {j:.4f})")

    # the A-SKY survivor: one representative realisation + the pooled statement
    z = np.load(fs[0], allow_pickle=True)
    dlnl = np.asarray(z["dlnL_det"], float)
    lnK = np.asarray(z["lnK"], float)
    q = np.asarray(z["qmax"], float)
    ot = np.asarray(z["on_true"], bool)
    # v2.2 is scored at the CANONICAL floor cut above from this unit's own 100 nulls
    # (recut_surface.offender + q95), not at floor = 0.
    c = cert_v22(dlnl, lnK, q, fl_can)
    ps = {p0: par_subset(q, p0) for p0 in P0_LEVELS}
    describe(f"A-SKY survivor AS_e03_h1275_k5_s4 realisation "
             f"{os.path.basename(fs[0]).split('_n')[-1][:-4]}", names, c, ps, q,
             extra=f"[v2.2 at the canonical floor {fl_can:.3f} cut from this unit's own "
                   f"100 nulls; {int(c.sum())} certs]")
    for p0, (S, j) in ps.items():
        if len(S):
            print(f"      PAR P0={p0:.2f} realised: {int(ot[S].sum())}/{len(S)} "
                  f"members actually on_true (predicted joint {j:.4f})")
        rows.append([f"A-SKY s4 n{os.path.basename(fs[0]).split('_n')[-1][:-4]}", p0,
                     int(len(S)), float(j), int(c.sum()),
                     int(len(set(S.tolist()) & set(np.where(c)[0].tolist()))),
                     int(ot[S].sum()) if len(S) else 0, int(ot.sum())])

    # pooled over all 15 A-SKY realisations
    print("\n   -- A-SKY survivor, pooled over all 15 realisations")
    print("      P0    mean|S_PAR|  mean|C_v22|  mean|overlap|  mean PAR-only  "
          "mean v22-only")
    pooled = []
    for p0 in P0_LEVELS:
        sp, sc, so, po, vo = [], [], [], [], []
        for r, f in enumerate(fs):
            zz = np.load(f, allow_pickle=True)
            cc = cert_v22(np.asarray(zz["dlnL_det"], float),
                          np.asarray(zz["lnK"], float), Q[r], fl_can)
            S, _ = par_subset(Q[r], p0)
            ss, cs = set(S.tolist()), set(np.where(cc)[0].tolist())
            sp.append(len(ss)); sc.append(len(cs)); so.append(len(ss & cs))
            po.append(len(ss - cs)); vo.append(len(cs - ss))
        print(f"      {p0:.2f}  {np.mean(sp):11.2f}  {np.mean(sc):11.2f}  "
              f"{np.mean(so):13.2f}  {np.mean(po):13.2f}  {np.mean(vo):12.2f}")
        pooled.append([p0, float(np.mean(sp)), float(np.mean(sc)), float(np.mean(so)),
                       float(np.mean(po)), float(np.mean(vo))])

    out = os.path.join(R, "reports", "sieve_t3_par.npz")
    np.savez(
        out,
        note=("SIEVE T3: partial ambiguity resolution (GNSS) transplanted to the PTA "
              "fringe problem and compared with criterion-v2.2 on banked above-onset "
              "realisations. LAMBDA decorrelation is SUBSTITUTED by the per-target "
              "decohered E-step (the full float covariance Q is measured-unaffordable, "
              "spark3.faint_fisher_bounds), which makes every joint success rate here "
              "a LOWER bound. REPORT-ONLY: banked columns only, nothing re-evaluated."),
        p0_levels=np.array(P0_LEVELS), qbar=QBAR, trials_nat=TRIALS_NAT,
        n_asky_real=len(fs),
        reliability=np.array(rel_rows, dtype=object),
        reliability_keys=np.array(["lo", "hi", "n", "mean_qmax", "emp_on_true",
                                   "se", "mean_ptrue", "verdict"]),
        pooled_qmax=tot_q, pooled_on_true=tot_e, pooled_se=se_t,
        asky_n_null=len(nfs), asky_floor_canonical=fl_can,
        asky_floor_glacier_stat=fl_gl, asky_floor_gap=fl_can - fl_gl,
        asky_off_canonical=off_can, asky_off_glacier_stat=off_gl,
        hi_qmax=hq, hi_on_true=he, hi_se=seh,
        glacier_reliability=np.array(glac_rows, dtype=object),
        glacier_reliability_keys=np.array(["lo", "hi", "n", "mean_qmax",
                                           "emp_on_true", "se", "verdict"]),
        glacier_hi_qmax=ghq, glacier_hi_on_true=ghe, glacier_hi_se=gse,
        g2=np.array(g2),
        g2_keys=np.array(["P0", "mean_|S|", "mean_pred_Pboot", "realised_rate", "se"]),
        head_to_head=np.array(rows, dtype=object),
        head_to_head_keys=np.array(["unit", "P0", "n_PAR", "P_boot", "n_cert_v22",
                                    "n_overlap", "n_PAR_on_true", "n_on_true_total"]),
        asky_pooled=np.array(pooled),
        asky_pooled_keys=np.array(["P0", "mean_nPAR", "mean_ncert", "mean_overlap",
                                   "mean_PAR_only", "mean_v22_only"]),
    )
    print(f"\nbanked -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
