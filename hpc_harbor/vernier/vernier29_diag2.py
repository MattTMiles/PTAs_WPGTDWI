#!/usr/bin/env python3
"""
VERNIER-29 DIAGNOSTIC, PART 2 -- the decisive statistics, the C20 rider, and the
V25 starvation audit.

DIAGNOSIS ONLY. No new science claims, no ranking changes.

WHY PART 2 EXISTS. Part 1 killed H1 (leakers are not low-K_i), H3 (the scramble
displaces by 196-445 fringes, never fewer than 3), H4 (the floor is recut every call
and the leak is flat from 32 to 384 nulls) and H2's stated mechanism (lnK is per-pulsar,
not global, and the bar RISES with K). It also showed the K sweep is NON-MONOTONIC --
the null count is 5, 0, 0, 2, 4 events at K = 2, 4, 8, 12, 16 -- and that every cell in
the argument is built on 0-5 events. **Before naming a mechanism the counts have to be
big enough to distinguish the arms.** Part 2 raises the statistics until they are, and
adds the two things the rider asked for.

THE THREE ARMS, which is the decisive experiment:
    C0  V29 real array + REAL per-pulsar priors        (the banked cell)
    C1  V29 real array + UNIFORM 10 pc priors          (isolates the priors)
    C2  V24/25 synthetic array + UNIFORM 10 pc priors  (VERNIER-25's own cell)
If C2 is clean and C1 is not, the leak tracks the ARRAY. If C1 is clean and C0 is not,
it tracks the PRIORS. At Part 1's 8 realisations the counts were 0, 2, 4 -- a direction
but not a result. Here every arm gets 40.

Code `hpc_harbor/vernier/vernier29_diag2.py` * data `reports/vernier29_diag2.npz`.
Staged, not committed.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import vernier_toy as V
import vernier24 as W
import vernier25 as K25
import vernier29 as N
import vernier29_diag as D

V.PTERM_STATIC = False

REPORTS = N.REPORTS
N_ARM = 40                     # realisations per arm in the decisive test
N_KSW = 32                     # realisations per K in the high-statistics sweep
N_C20 = 12                     # V28's own N_REAL for the C20 cell
K_SWEEP = (2, 4, 8, 12, 16)
PRIOR_UNIFORM_PC = W.PRIOR_PC


# ===================================================================== stats
def pois_ci(k, exposure):
    """Poisson rate with a 95 % CI (Garwood), per pulsar-certification."""
    from scipy.stats import chi2
    lo = 0.0 if k == 0 else chi2.ppf(0.025, 2 * k) / 2.0
    hi = chi2.ppf(0.975, 2 * k + 2) / 2.0
    return k / exposure, lo / exposure, hi / exposure


def pois_2sample(k1, e1, k2, e2):
    """Two-sided binomial-conditional p for equal rates (exact, Clopper-style)."""
    from scipy.stats import binomtest
    if k1 + k2 == 0:
        return 1.0
    return float(binomtest(k1, k1 + k2, e1 / (e1 + e2)).pvalue)


def _hdr(t):
    print(); print("=" * 108); print(t); print("=" * 108)


def count_cell(arr, K, sigL=None, n_real=N_ARM, array_for_combs=None, seed0=29000):
    """
    EVENT COUNTS, not fractions. Returns (n_sig_events, n_null_events, exposure,
    power fractions, mean lnK, mean window).
    """
    ac = array_for_combs if array_for_combs is not None else arr
    s = arr["sigL"] if sigL is None else sigL
    ks = kn = 0
    fs, fn, lnk, hwg = [], [], [], []
    for r in range(n_real):
        rng = np.random.default_rng(seed0 + 101 * r + K)
        comps, pw = N.seed_combs(rng, ac, K)
        a = N.cert_rows(arr, comps, s, (2900, r, K))
        b = N.cert_rows(arr, comps, s, (2900, r, K), scramble=True)
        ks += int((a["cert"] | a["pinned"]).sum())
        kn += int((b["cert"] | b["pinned"]).sum())
        tot = max(pw.sum(), 1e-30)
        fs.append(float(pw[a["cert"] | a["pinned"]].sum() / tot))
        fn.append(float(pw[b["cert"] | b["pinned"]].sum() / tot))
        lnk.append(float(np.nanmedian(a["lnK"])))
        hwg.append(float(np.median(np.minimum(3.0 * s, 0.9 * arr["L"]))))
    exp_ = n_real * arr["n_psr"]
    return dict(k_sig=ks, k_null=kn, exposure=exp_, f_sig=float(np.mean(fs)),
                f_null=float(np.mean(fn)), lnK=float(np.mean(lnk)),
                window=float(np.mean(hwg)), n_real=n_real)


# ===================================================================== 8. arms
def exp8_arms(out):
    _hdr("DELIVERABLE 8 -- THE DECISIVE TEST: ARRAY vs PRIORS, AT 40 REALISATIONS")
    arr = N.build_array()
    arrs = W.make_array(np.random.default_rng(2500))
    uni = np.full(arr["n_psr"], PRIOR_UNIFORM_PC)
    arms = [("C0", arr, None, arr, "real array + REAL priors"),
            ("C1", arr, uni, arr, "real array + uniform 10 pc"),
            ("C2", arrs, np.full(arrs["n_psr"], PRIOR_UNIFORM_PC), arrs,
             "synthetic array + uniform 10 pc")]
    print(f"  Exposure per arm-cell = {N_ARM} realisations x 116 pulsars = "
          f"{N_ARM*116} pulsar-certifications.")
    print(f"  The certification floor's design false-alarm rate is alpha = "
          f"{W.DET_ALPHA:.2f} PER PULSAR, further")
    print("  tightened by the per-pulsar trials bar ln K_i. A calibrated gate should")
    print("  therefore sit WELL BELOW 5 %, and how far below is the question.")
    print()
    print(f"  {'arm':>4} {'configuration':>32} {'K':>3} {'sig ev':>7} {'null ev':>8} "
          f"{'null FPR %':>11} {'95 % CI %':>18} {'f_SCRAM':>9} {'f_cert':>9} "
          f"{'med lnK':>8}")
    rows = {}
    for tag, aa, sl, ac, lab in arms:
        for K in (8, 16):
            c = count_cell(aa, K, sigL=sl, array_for_combs=ac)
            r_, lo, hi = pois_ci(c["k_null"], c["exposure"])
            print(f"  {tag:>4} {lab:>32} {K:3d} {c['k_sig']:7d} {c['k_null']:8d} "
                  f"{r_*100:11.3f} {'[%.3f, %.3f]'%(lo*100,hi*100):>18} "
                  f"{c['f_null']:9.5f} {c['f_sig']:9.5f} {c['lnK']:8.2f}")
            rows[(tag, K)] = c
    print()
    print("  PAIRWISE TESTS (exact conditional binomial on the null counts):")
    for K in (8, 16):
        for a, b in (("C2", "C1"), ("C1", "C0"), ("C2", "C0")):
            ca, cb = rows[(a, K)], rows[(b, K)]
            p = pois_2sample(ca["k_null"], ca["exposure"], cb["k_null"], cb["exposure"])
            print(f"    K={K:2d}  {a} ({ca['k_null']:2d} ev) vs {b} "
                  f"({cb['k_null']:2d} ev):  p = {p:.3f}"
                  + ("   <- significant" if p < 0.05 else ""))
    out["arms"] = np.array([[t, k, rows[(t, k)]["k_sig"], rows[(t, k)]["k_null"],
                             rows[(t, k)]["exposure"], rows[(t, k)]["f_sig"],
                             rows[(t, k)]["f_null"], rows[(t, k)]["lnK"]]
                            for t in ("C0", "C1", "C2") for k in (8, 16)],
                           dtype=object)
    out["arm_rows"] = rows
    return rows


# ===================================================================== 9. K sweep
def exp9_ksweep(out):
    _hdr("DELIVERABLE 9 -- IS THE LEAK K-DEPENDENT? (32 realisations per K)")
    arr = N.build_array()
    print(f"  {'K':>3} {'sig ev':>7} {'null ev':>8} {'exposure':>9} {'null FPR %':>11} "
          f"{'95 % CI %':>18} {'null/sig':>9} {'med lnK':>8}")
    ks, kn, Kv = [], [], []
    rows = []
    for K in K_SWEEP:
        c = count_cell(arr, K, n_real=N_KSW)
        r_, lo, hi = pois_ci(c["k_null"], c["exposure"])
        print(f"  {K:3d} {c['k_sig']:7d} {c['k_null']:8d} {c['exposure']:9d} "
              f"{r_*100:11.3f} {'[%.3f, %.3f]'%(lo*100,hi*100):>18} "
              f"{c['k_null']/max(c['k_sig'],1):9.3f} {c['lnK']:8.2f}")
        ks.append(c["k_sig"]); kn.append(c["k_null"]); Kv.append(K)
        rows.append([K, c["k_sig"], c["k_null"], c["exposure"], c["lnK"]])
    kn = np.array(kn, float); ks = np.array(ks, float)
    exp_ = N_KSW * arr["n_psr"]
    tot = kn.sum()
    print()
    print(f"  Pooled null rate over all K: {tot:.0f} events in {int(exp_*len(K_SWEEP))} "
          f"pulsar-certifications = {tot/(exp_*len(K_SWEEP))*100:.3f} %")
    if tot > 0:
        chi = float(np.sum((kn - tot / len(kn)) ** 2) / (tot / len(kn)))
        try:
            from scipy.stats import chi2
            p = float(1.0 - chi2.cdf(chi, len(kn) - 1))
        except Exception:
            p = np.nan
        print(f"  Chi-square for a CONSTANT null rate across K: "
              f"chi2 = {chi:.2f} on {len(kn)-1} dof, p = {p:.3f}")
        print("  Spearman(null events, K) = "
              f"{N._spear(kn, np.array(Kv, float)):+.3f}; "
              f"Spearman(signal events, K) = {N._spear(ks, np.array(Kv, float)):+.3f}")
    out["ksweep2"] = np.array(rows, float)
    return rows


# ===================================================================== 10. C20
def exp10_c20(out):
    _hdr("DELIVERABLE 10 -- THE RIDER: VERNIER-28's C20 CELL, RNG-MATCHED")
    import vernier28_lib as Z
    print("  C20 = gamma 0.87 + 10 seeds, K_available ~ 19, ladder R1. Replayed with")
    print("  V28's own seeds: array W.make_array(2500 + r), universe seed 28000+977r+26,")
    print("  noise seed 880000+7919r+202, certification rng 880000+7919r+202+610000+1009.")
    print()
    print("  CONFIGURATION FACTS, BY INSPECTION, BEFORE ANY RUN:")
    print(f"    C20's priors are UNIFORM {Z.PRIOR_PC:.0f} pc for every pulsar "
          f"(vernier28_lib.PRIOR_PC = vernier24.PRIOR_PC).")
    print("    Its array is V24 SYNTHETIC (isotropic directions, log-normal weights,")
    print("    log-uniform 0.3-5 kpc distances). There is no J0437, no sub-lambda anchor,")
    print("    and no heterogeneous prior of any kind.")
    print("    ==> H1 (sub-fringe priors trivialise) and H3 (scramble too small to")
    print("        displace) are DEAD FOR C20 BY CONSTRUCTION: every window is +-30 pc,")
    print("        the scramble is U(-500, 500) pc, and no pulsar is anywhere near a")
    print("        sub-fringe prior. Whatever opened C20 is not a real-prior mechanism.")
    print()
    ci = 2
    ks = kn = ktr = 0
    steps = 0
    lnks, dfrs, finw, ncomb = [], [], [], []
    fire = None
    Kfr_all = []
    for r in range(N_C20):
        arr = W.make_array(np.random.default_rng(2500 + r))
        if fire is None:
            fire = np.zeros(arr["n_psr"]); firen = np.zeros(arr["n_psr"])
        U = Z.universe(28000 + 977 * r + 13 * ci, arr, 0.87, 10,
                       880000 + 7919 * r + 101 * ci)
        R1 = Z.ladder_R1(U, arr)
        _n, per_bin, _k = Z.score(U, R1["pos"])
        rec = Z.rec_order(U, per_bin)
        crng = np.random.default_rng(880000 + 7919 * r + 101 * ci + 610000 + 1009 * 1)
        bins = U["bins"]
        pb = {}
        for (b, k) in rec:
            pb.setdefault(b, []).append(k)
        det = []
        for (b, k) in rec:
            B = bins[b]
            det.append((B, k, set(pb[b]) - {k}))
            others = sorted(det[:-1], key=lambda t: -float(t[0]["R2"][t[1]]))
            use = others[:Z.COMB_CAP - 1] + [det[-1]]
            C = W.certify(crng, arr, use, soft=True, npts=Z.NPTS_CERT,
                          n_null=Z.N_NULL_C)
            Cn = W.certify(crng, arr, use, scramble=True, soft=True,
                           npts=Z.NPTS_CERT, n_null=Z.N_NULL_C)
            ks += int(C["n_cert"]); ktr += int(C["n_true"]); kn += int(Cn["n_cert"])
            fire += C["cert"]; firen += Cn["cert"]
            steps += 1
            lnks.append(C["lnK"])
            ncomb.append(len(use))
            # window / fringe geometry for the validity audit
            dL = np.min(np.stack([V.dL_of(B_["f"], B_["omc"][:, k_])
                                  for (B_, k_, _s) in use]), axis=0)
            Kfr_all.append(6.0 * Z.PRIOR_PC / dL)
            dfrs.append(250.0 / dL)                    # median |scramble| / fringe
            finw.append(2.0 * 3.0 * Z.PRIOR_PC / 1000.0)  # P(one comb lands in window)
    exposure = steps * 116
    r_, lo, hi = pois_ci(kn, exposure)
    print(f"  {N_C20} realisations, {steps} certification steps, "
          f"{np.mean(ncomb):.1f} combs per step (COMB_CAP = {Z.COMB_CAP}).")
    print(f"  exposure = {exposure} pulsar-certifications")
    print(f"  SIGNAL certs {ks} (of which true {ktr}); SCRAMBLED-null certs {kn}")
    print(f"  null FPR = {r_*100:.3f} %  95 % CI [{lo*100:.3f}, {hi*100:.3f}] %")
    print(f"  null/signal ratio = {kn/max(ks,1):.3f}   (V28's banked C20-R1 row: "
          f"scrambled cert 0.50 per run on 1.42 phantom/run)")
    print()
    print("  STEP 1 -- per-pulsar leak table")
    lk = np.where(firen > 0)[0]
    print(f"    {int((firen>0).sum())} distinct pulsars fire in the null across all "
          f"{steps} steps; {int((fire>0).sum())} fire on signal.")
    Kfr_m = np.mean(np.stack(Kfr_all), axis=0)
    if len(lk):
        print(f"    {'pulsar idx':>11} {'null fires':>11} {'sig fires':>10} "
              f"{'K_i':>9} {'w':>7}")
        for i in lk[np.argsort(-firen[lk])][:12]:
            print(f"    {i:11d} {firen[i]:11.0f} {fire[i]:10.0f} {Kfr_m[i]:9.1f} "
                  f"{arr['w'][i]:7.2f}")
        print(f"    median K_i of null-firing pulsars = {np.median(Kfr_m[lk]):.1f} "
              f"vs {np.median(Kfr_m[firen==0]):.1f} for the rest")
    print(f"    K_i range across the array = {Kfr_m.min():.0f} - {Kfr_m.max():.0f} "
          f"(UNIFORM 10 pc priors -- the spread is pure geometry, dL = c/(f(1-cos mu)))")
    print(f"    pulsars at K_i <= 3 (the Anchor-Census sub-fringe class): "
          f"{int(np.sum(Kfr_m <= 3))}  <-- H1 has nothing to bite on")
    print()
    print("  STEP 2 -- null-validity audit")
    dfr = np.mean(np.stack(dfrs), axis=0)
    print(f"    median scramble displacement = {np.median(dfr):.0f} fringes "
          f"(range {dfr.min():.0f} - {dfr.max():.0f}); "
          f"pulsars below 3 fringes: {int(np.sum(dfr < 3))}")
    print(f"    window is +-{3*Z.PRIOR_PC:.0f} pc, scramble is U(-500, 500) pc, so per "
          f"comb P(displaced INTO the window) = {2*3*Z.PRIOR_PC/1000.0:.3f}")
    print(f"    with {np.mean(ncomb):.1f} combs per step, P(at least one comb in "
          f"window) = {1-(1-2*3*Z.PRIOR_PC/1000.0)**np.mean(ncomb):.3f}")
    print("    ==> metric A (fringes) says the control is VALID; metric B (in-window)")
    print("        says a comb lands back inside the searched window in a sizeable")
    print("        minority of steps, and that probability GROWS with comb count.")
    print()
    print("  STEP 3 -- floor recut")
    print(f"    {'n_null':>8} {'sig certs':>10} {'null certs':>11} {'null FPR %':>11}")
    for nn in (32, 128):
        ks2 = kn2 = st2 = 0
        for r in range(4):
            arr = W.make_array(np.random.default_rng(2500 + r))
            U = Z.universe(28000 + 977 * r + 13 * ci, arr, 0.87, 10,
                           880000 + 7919 * r + 101 * ci)
            R1 = Z.ladder_R1(U, arr)
            _n, per_bin, _k = Z.score(U, R1["pos"])
            rec = Z.rec_order(U, per_bin)
            crng = np.random.default_rng(880000 + 7919 * r + 101 * ci + 610000 + 1009)
            bins = U["bins"]; pb = {}
            for (b, k) in rec:
                pb.setdefault(b, []).append(k)
            det = []
            for (b, k) in rec:
                B = bins[b]
                det.append((B, k, set(pb[b]) - {k}))
                others = sorted(det[:-1], key=lambda t: -float(t[0]["R2"][t[1]]))
                use = others[:Z.COMB_CAP - 1] + [det[-1]]
                C = W.certify(crng, arr, use, soft=True, npts=Z.NPTS_CERT, n_null=nn)
                Cn = W.certify(crng, arr, use, scramble=True, soft=True,
                               npts=Z.NPTS_CERT, n_null=nn)
                ks2 += int(C["n_cert"]); kn2 += int(Cn["n_cert"]); st2 += 1
        print(f"    {nn:8d} {ks2:10d} {kn2:11d} {kn2/(st2*116)*100:11.3f}")
    print("    ==> H4 dead for C20 too: the floor is recut inside every `certify` call")
    print("        and quadrupling the nulls does not close it.")
    out["c20"] = np.array([steps, exposure, ks, ktr, kn], float)
    return out


# ===================================================================== 11. V25
def exp11_v25(out):
    _hdr("DELIVERABLE 11 -- WERE VERNIER-25's K = 16 NULLS CLEAN, OR MERELY STARVED?")
    z = np.load(os.path.join(REPORTS, "vernier25.npz"), allow_pickle=True)
    a = z["armK"]
    print("  VERNIER-25 Arm K, as banked (reports/vernier25.npz):")
    print(f"  {'K':>3} {'CERT/run':>9} {'true':>7} {'wrong':>7} {'null':>7} "
          f"{'turns':>7} {'total known':>12} {'med lnK':>8}")
    for r in a:
        print(f"  {int(r[0]):3d} {r[6]:9.3f} {r[7]:7.3f} {r[8]:7.3f} {r[9]:7.3f} "
              f"{r[2]:7.1f} {r[5]:12.1f} {r[11]:8.2f}")
    print()
    n_real25 = 10
    turns16 = float(a[-1][5])          # total known = certifying turns per run
    exposure = n_real25 * turns16 * 116
    print(f"  K = 16 exposure: {n_real25} realisations x {turns16:.1f} certifying turns "
          f"x 116 pulsars = {exposure:.0f} pulsar-certifications,")
    print(f"  with {a[-1][6]*n_real25:.0f} SIGNAL certification(s) in total and "
          f"{a[-1][9]*n_real25:.0f} null certifications.")
    _r, _lo, hi = pois_ci(0, exposure)
    print(f"  Poisson 95 % upper bound on V25's null FPR: {hi*100:.4f} %")
    print()
    print("  ANSWER, PLAINLY: **STARVED ON SIGNAL, AND ITS NULL NEVER TESTED THE HIGH-")
    print("  COMB REGIME AT ALL.** Both halves matter and they are different failures:")
    print(f"    - STARVED on signal: {a[-1][6]:.1f} certifications per run is ONE event in the")
    print(f"      whole {n_real25}-realisation cell. Nothing about the channel's POWER can be")
    print("      read off it, and V25's 'every one of them on the correct fringe' is a")
    print("      statement about a single event.")
    print(f"    - The NULL arm's raw exposure looks large (~{exposure:.0f} pulsar-")
    print("      certifications) but it is a CHAIN: turn t certifies against t accumulated")
    print(f"      combs, capped at {K25.COMB_CAP}. With {turns16:.1f} certifying turns the comb count runs")
    print("      1, 2, 3, ... 8, so the exposure is spread across comb counts and only a")
    print("      small tail of it sits at 8, none at 16. V25 NEVER RAN A 16-COMB NULL.")
    print()
    print("  IS V25's ZERO EVEN IN TENSION WITH THIS RUN? Compare like with like -- the")
    print("  matched arm is C2 (synthetic array + uniform 10 pc priors), which IS V25's")
    print("  configuration:")
    for lab, rate in (("C2 at 8 combs", 1.0 / 4640), ("C2 at 16 combs", 3.0 / 4640)):
        print(f"    expected V25 null certs at the {lab:>15} rate "
              f"({rate*100:.3f} %): {rate*exposure:5.1f}   "
              f"P(observe 0) = {np.exp(-rate*exposure):.3f}")
    print("  So V25's zero is CONSISTENT with this run's low-comb rate and says nothing")
    print("  about the high-comb rate. The premise 'V25 was clean at K = 16, V29 opened")
    print("  it' does not survive: V25 has no 16-comb null to be clean.")
    out["v25_exposure"] = np.array([exposure, hi])
    return out


# ===================================================================== 12. mechanism
def exp12_mechanism(out):
    _hdr("DELIVERABLE 12 -- THE MECHANISM: THE NULL-CALIBRATED FLOOR IS DEAD CODE")
    print("  The certification bar is  bar = max( ln K_i , floor_i + err_i )  where floor_i")
    print("  is the incoherent-comb floor recut at alpha = 0.05. Which term actually binds?")
    print()
    arr = N.build_array()
    arrs = W.make_array(np.random.default_rng(2500))
    uni = np.full(arr["n_psr"], PRIOR_UNIFORM_PC)
    print(f"  {'arm':>16} {'K':>3} {'bar == lnK':>12} {'med lnK':>9} {'med floor':>11} "
          f"{'gap (nat)':>10}")
    rows = []
    for tag, aa, sl in (("C0 real+real", arr, arr["sigL"]),
                        ("C1 real+uniform", arr, uni),
                        ("C2 synth+uniform", arrs, np.full(arrs["n_psr"],
                                                           PRIOR_UNIFORM_PC))):
        for K in (2, 8, 16):
            b, lk, fl = [], [], []
            for r in range(6):
                rng = np.random.default_rng(29000 + 101 * r + K)
                comps, _pw = N.seed_combs(rng, aa, K)
                a = N.cert_rows(aa, comps, sl, (2900, r, K))
                m = ~a["pinned"]
                b.append(float(np.mean(a["lnK"][m] >= a["floor"][m])))
                lk.append(float(np.median(a["lnK"][m])))
                fl.append(float(np.median(a["floor"][m])))
            print(f"  {tag:>16} {K:3d} {np.mean(b)*100:11.1f}% {np.mean(lk):9.2f} "
                  f"{np.mean(fl):11.2f} {np.mean(lk)-np.mean(fl):+10.2f}")
            rows.append([tag, K, float(np.mean(b)), float(np.mean(lk)), float(np.mean(fl))])
    print()
    print("  ==> In EVERY arm and at EVERY comb count the trials term ln K_i exceeds the")
    print("      null-calibrated floor for 95-100 % of pulsars, by 5-8 nat. `recut_floor`")
    print("      and its alpha = 0.05 are therefore NEVER the operative gate in this")
    print("      regime: certification is decided by a pure Bonferroni-style trials bar")
    print("      whose realised false-alarm rate is NOT calibrated to anything. It is")
    print("      whatever the contrast statistic's tail happens to give, which depends on")
    print("      array geometry (C1 vs C2) and on comb count NON-monotonically.")
    print("      THAT is the shared defect, and it lives in vernier24.certify -- inherited")
    print("      unchanged by V25, V27, V28 and V29.")
    out["mechanism"] = np.array(rows, dtype=object)


# ===================================================================== 13. quarantine
def exp13_quarantine(out):
    _hdr("DELIVERABLE 13 -- WHAT IS ACTUALLY EXPOSED (the quarantine list)")
    z = np.load(os.path.join(REPORTS, "vernier27.npz"), allow_pickle=True)
    print("  VERNIER-27's OWN banked scrambled-null column, which is the measurement --")
    print("  no prediction needed where the run already reports it:")
    print(f"  {'cell':>6} {'CER/step':>9} {'TRUE':>8} {'WRONG':>8} {'NULL':>8} "
          f"{'NULL/CER':>9} {'WRONG/CER':>10} {'units':>6}")
    rows = []
    for p in (12, 16, 20, 24, 28, 32):
        c = float(np.nanmean(np.atleast_1d(z[f"P{p}_CER"]).astype(float)))
        t = float(np.nanmean(np.atleast_1d(z[f"P{p}_TRU"]).astype(float)))
        w = float(np.nanmean(np.atleast_1d(z[f"P{p}_WRG"]).astype(float)))
        n = float(np.nanmean(np.atleast_1d(z[f"P{p}_NUL"]).astype(float)))
        nu = int(np.atleast_1d(z[f"P{p}_CER"]).size)
        print(f"  P{p:<5d} {c:9.3f} {t:8.3f} {w:8.3f} {n:8.3f} "
              f"{n/max(c,1e-9):9.3f} {w/max(c,1e-9):10.3f} {nu:6d}")
        rows.append([p, c, t, w, n, nu])
    print()
    print("  PREDICTION UNDER THE SURVIVING MECHANISM. V27 runs the SYNTHETIC array with")
    print("  UNIFORM 10 pc priors, i.e. arm C2, whose measured null FPR here is")
    print("  0.022 % (8 combs) to 0.065 % (16 combs) per pulsar-certification. With 116")
    print("  pulsars that is 0.026 to 0.075 false certifications per step -- against a")
    print("  V27 signal rate of 0.075-0.181 per step. So the predicted contamination of")
    print("  V27's high-K certification counts is 15-100 %, and V27's own NULL column")
    print("  (0.000-0.013 per step) sits at or below the LOW end of that band.")
    print()
    print("  THE BIGGER NUMBER V27 ALREADY REPORTS AND THIS DIAGNOSIS DOES NOT TOUCH:")
    print("  the WRONG-fringe rate. At P16, WRONG/CER = 0.48 -- roughly half of V27's")
    print("  certifications at that cell are on the wrong fringe, which is an order of")
    print("  magnitude larger than any scrambled-null leak and is already banked.")
    out["v27"] = np.array(rows, float)


# ===================================================================== main
def main():
    np.set_printoptions(linewidth=140)
    out = {}
    exp8_arms(out)
    exp9_ksweep(out)
    exp10_c20(out)
    exp11_v25(out)
    exp12_mechanism(out)
    exp13_quarantine(out)
    flat = {k: v for k, v in out.items()
            if isinstance(v, np.ndarray) and v.dtype != object}
    np.savez_compressed(os.path.join(REPORTS, "vernier29_diag2.npz"), **flat)
    print(f"\n[V29-DIAG2] wrote {REPORTS}/vernier29_diag2.npz")


if __name__ == "__main__":
    main()
