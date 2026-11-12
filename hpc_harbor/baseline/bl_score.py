#!/usr/bin/env python
"""BASELINE scorer -- CPU only, banked npz only, no new realisations.

Produces B3: one table (detection at matched FAP; certification purity vs the calibrated
criterion) and one panel. Every floor carries its estimator and its ZERO-FRACTION column
(criterion-v2.2). Every fraction carries a Wilson interval (the size-glance).

THE TWO THRESHOLDS ARE MATCHED BY CONSTRUCTION. Both are the (1 - alpha) = 0.95 point of a
statistic's distribution over THE SAME 100 pure-noise realisations:
  * FIELD      : the max-over-grid F_e statistic 2F  -> t95
  * CRITERION  : the D2 'offender' (max over pulsars of dlnL among pulsars passing
                 dlnL > lnK and q > 0.9; 0 if none) -> floor
so "detected at FAP = 0.05" means the same thing on both sides. The criterion's per-pulsar
gate then also carries its own trials factor lnK, which the field's has no analogue of --
that asymmetry is the thing being measured, not a defect of the comparison.

Run: bl_score.py
"""
import os, sys, glob, json
import numpy as np
from scipy.stats import gumbel_r

HSYMT = os.environ.get("BL_HSYMT", "/home/mattm/projects/HSYMT")
OUT = os.environ.get("BL_OUT", f"{HSYMT}/BASELINE_results")
REPORTS = os.environ.get("BL_REPORTS", f"{HSYMT}/reports")

ALPHA = 0.05
Z_ALPHA = 2.9701952521018403
C_SD = 2.80
ZF_GATE = 0.20
QBAR = 0.9
BOOT = 4000
BOOT_SEED = 20260729

SF_OUT = os.environ.get("BL_SF", f"{HSYMT}/SURFACE_results")

CELLS = [("C1_census_h1275_T30_lit_k3", "census 3+13, h=-12.75, T=30, lit"),
         ("C2_askysurv_e03_h1275_k5", "A-SKY survivor 5+11, e=0.3, h=-12.75, T=30, lit")]

# validated categorical slots 1/2/3 (dataviz reference palette, light mode)
CFIELD, CFIELD_D1, CCRIT = "#2a78d6", "#eb6834", "#1baf7a"


# ============================================================
# the criterion-v2.2 floor machinery (recut_surface.py, verbatim)
# ============================================================
def offender(dlnL, lnK, qmax, elig=None):
    """Largest dlnL among pulsars passing layers 1+3. MAX OVER PULSARS -> Gumbel domain.

    `elig` (the ELIGIBILITY mask, BASELINE-local and applied to BOTH arms identically) drops
    pulsars whose window holds a single fringe. On the criterion's comb that mask is all-True
    and this is the incumbent `offender` verbatim; on the FIELD's comb a mis-estimated sky can
    make dL_a diverge for a pulsar near the estimated source direction, leaving one fringe --
    a pulsar with no fringe ambiguity, whose `dlnL` is +inf by construction and whose q_max is
    trivially 1. Scoring that as an infinitely significant certification would be an artefact,
    and it would be an artefact that FLATTERS the field, so it is removed."""
    m = (dlnL > lnK) & (qmax > QBAR)
    if elig is not None:
        m = m & elig
    return float(dlnL[m].max()) if m.any() else 0.0


def gumbel_floor(x):
    x = np.asarray(x, float)
    mu, beta = gumbel_r.fit(x)
    return (float(mu + beta * Z_ALPHA), float(mu), float(beta),
            float(C_SD * beta / np.sqrt(len(x))), len(x))


def emp_quantile(off, alpha=ALPHA):
    return float(np.quantile(np.asarray(off, float), 1.0 - alpha))


def boot_sd(off, alpha=ALPHA, B=BOOT, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    off = np.asarray(off, float); n = len(off)
    return float(np.std([np.quantile(rng.choice(off, n, replace=True), 1.0 - alpha)
                         for _ in range(B)]))


def adopt(off):
    """(floor, err, estimator, zero_fraction) under ANCHOR S4 / criterion-v2.2."""
    off = np.asarray(off, float)
    zf = float(np.mean(off == 0.0))
    gu, _, _, gu_sd, _ = gumbel_floor(off)
    if zf <= ZF_GATE:
        return float(gu), float(gu_sd), "gumbel", zf
    return emp_quantile(off), boot_sd(off), "emp_q95", zf


def wilson(k, n, z=1.96):
    """Wilson score interval -- the size-glance on every fraction quoted here."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


# ============================================================
# load
# ============================================================
def load(cell, kind):
    fs = sorted(glob.glob(f"{OUT}/bl_{kind}_{cell}_g*_n*.npz"))
    if not fs:
        return None
    R = []
    for f in fs:
        d = np.load(f, allow_pickle=True)
        R.append({k: d[k] for k in d.files})
    return R


def col(R, k):
    return np.array([r[k] for r in R])


# ============================================================
# scoring
# ============================================================
def score_cell(cell, NUL, SIG):
    """All of B1 and B2 for one cell. Raw columns in, fractions + intervals out."""
    n_s, n_n = len(SIG), len(NUL)

    # ---------- B1: the FIELD's F_e threshold at matched FAP ----------
    twoF_n, twoF_s = col(NUL, "twoF"), col(SIG, "twoF")
    t95 = emp_quantile(twoF_n)
    t95_sd = boot_sd(twoF_n)
    gu_t, gu_mu, gu_beta, gu_sd, _ = gumbel_floor(twoF_n)
    zf_t = float(np.mean(twoF_n == 0.0))
    fe_det = twoF_s > t95
    fe_det_gu = twoF_s > gu_t

    # ---------- eligibility: >= 2 fringes in the window (both arms, same rule) ----------
    eN_c = col(NUL, "c_nfringe") >= 2
    eN_f = col(NUL, "f_nfringe") >= 2
    eS_c = col(SIG, "c_nfringe") >= 2
    eS_f = col(SIG, "f_nfringe") >= 2
    excl = dict(crit_null=int((~eN_c).sum()), field_null=int((~eN_f).sum()),
                crit_sig=int((~eS_c).sum()), field_sig=int((~eS_f).sum()),
                n_psr_null=int(eN_c.size), n_psr_sig=int(eS_c.size))

    # ---------- B1: the CRITERION's floor on the SAME nulls ----------
    off_c = np.array([offender(np.nan_to_num(r["c_dlnL_det"], posinf=1e30),
                               r["c_lnK"], r["c_qmax"], r["c_nfringe"] >= 2)
                      for r in NUL])
    fl_c, fl_c_sd, est_c, zf_c = adopt(off_c)
    cdet = np.array([(np.nan_to_num(r["c_dlnL_det"], posinf=1e30)
                      > np.maximum(r["c_lnK"], fl_c)) for r in SIG]) & eS_c
    ccert = cdet & (col(SIG, "c_qmax") > QBAR)
    ccorr = ccert & col(SIG, "c_on_true")
    cwrong = ccert & ~col(SIG, "c_on_true")
    cri_det_real = cdet.any(axis=1)          # realisation-level detection (matched to F_e)

    # ---------- the FIELD's own per-pulsar floor, refit on the SAME nulls ----------
    off_f = np.array([offender(np.nan_to_num(r["f_dlnL_det"], posinf=1e30),
                               r["f_lnK"], r["f_qmax"], r["f_nfringe"] >= 2)
                      for r in NUL])
    fl_f, fl_f_sd, est_f, zf_f = adopt(off_f)

    # ---------- B2: the field's per-pulsar picks, on the F_e-DETECTED realisations ----
    D = fe_det
    fq = col(SIG, "f_qmax"); fot = col(SIG, "f_on_true")
    fd = np.nan_to_num(col(SIG, "f_dlnL_det"), posinf=1e30); fk = col(SIG, "f_lnK")
    conf = (fq > QBAR) & eS_f                               # the FIELD as practised
    conf_d1 = conf & (fd > np.maximum(fk, fl_f))            # the FIELD + the criterion's D1
    res = {}

    def purity(mask, on_true, sub):
        m = mask[sub]; t = on_true[sub]
        k = int((m & t).sum()); n = int(m.sum())
        lo, hi = wilson(k, n)
        return dict(n_conf=n, n_corr=k, n_wrong=n - k,
                    purity=(k / n if n else np.nan), lo=lo, hi=hi,
                    per_real=(n / max(1, int(sub.sum()))),
                    corr_per_real=(k / max(1, int(sub.sum()))),
                    wrong_per_real=((n - k) / max(1, int(sub.sum()))))

    # `ptrue == -1` marks a pulsar whose TRUE fringe was never sampled by the 512-point
    # window. On the criterion's comb that cannot happen (`draw_true_distances_tier` draws
    # n_true only over sampled fringes). On the FIELD's comb it is not guaranteed, so it is
    # MEASURED. Such a pick is genuinely wrong -- the field states a distance that is not the
    # true one -- but the reason is a window budget, not the comb mechanism, so purity is
    # reported BOTH ways: headline counts them wrong, `field_inwin` restricts to pulsars
    # whose truth was reachable and isolates the mechanism.
    in_win = col(SIG, "f_ptrue") >= 0.0
    excl["field_truth_unsampled"] = int((~in_win).sum())

    res["field"] = purity(conf, fot, D)
    res["field_inwin"] = purity(conf & in_win, fot, D)
    res["field_d1"] = purity(conf_d1, fot, D)
    res["crit"] = purity(ccert, col(SIG, "c_on_true"), D)
    res["crit_all"] = purity(ccert, col(SIG, "c_on_true"), np.ones(n_s, bool))
    if "fmc_qmax" in SIG[0]:
        mq = col(SIG, "fmc_qmax"); mot = col(SIG, "fmc_on_true")
        res["field_mc"] = purity((mq > QBAR) & (col(SIG, "fmc_nfringe") >= 2), mot, D)
    if "for_qmax" in SIG[0]:
        # the ORACLE-SOURCE bound: the field with its search error set to ZERO.
        oq = col(SIG, "for_qmax"); oot = col(SIG, "for_on_true")
        res["field_oracle"] = purity((oq > QBAR) & (col(SIG, "for_nfringe") >= 2), oot, D)

    # ---------- detection fractions ----------
    def frac(m):
        k, n = int(m.sum()), len(m)
        lo, hi = wilson(k, n)
        return dict(k=k, n=n, f=k / n, lo=lo, hi=hi)

    return dict(
        cell=cell, n_sig=n_s, n_null=n_n,
        fe=dict(t95=t95, t95_sd=t95_sd, gumbel=gu_t, gumbel_sd=gu_sd, mu=gu_mu,
                beta=gu_beta, zero_frac=zf_t,
                null_med=float(np.median(twoF_n)), null_max=float(twoF_n.max()),
                sig_med=float(np.median(twoF_s)), sig_max=float(twoF_s.max()),
                det=frac(fe_det), det_gumbel=frac(fe_det_gu)),
        crit=dict(floor=fl_c, err=fl_c_sd, est=est_c, zero_frac=zf_c,
                  det=frac(cri_det_real),
                  cert_per_real=float(ccert.sum()) / n_s,
                  corr_per_real=float(ccorr.sum()) / n_s,
                  wrong_per_real=float(cwrong.sum()) / n_s),
        fieldfloor=dict(floor=fl_f, err=fl_f_sd, est=est_f, zero_frac=zf_f),
        pur=res, excl=excl,
        # diagnostics
        sep_deg=col(SIG, "sep_deg") if "sep_deg" in SIG[0] else None,
        f_dL_med=col(SIG, "f_dL_med"), c_dL_med=col(SIG, "c_dL_med"),
        polish=dict(
            # relative bar: the grid value comes from the VMAPPED profile and the polish
            # seed from a non-vmapped call, so they agree only to ~1e-6 relative.
            sig_rate=float(np.mean(twoF_s - col(SIG, "twoF_grid")
                                   > 1e-6 * np.maximum(1.0, np.abs(col(SIG, "twoF_grid"))))),
            null_rate=float(np.mean(twoF_n - col(NUL, "twoF_grid")
                                    > 1e-6 * np.maximum(1.0, np.abs(col(NUL, "twoF_grid"))))),
            sig_med=float(np.median(twoF_s - col(SIG, "twoF_grid"))),
            null_med=float(np.median(twoF_n - col(NUL, "twoF_grid"))),
            sig_max=float(np.max(twoF_s - col(SIG, "twoF_grid"))),
            worst=float(np.min(np.concatenate([twoF_s - col(SIG, "twoF_grid"),
                                               twoF_n - col(NUL, "twoF_grid")])))),
        twoF_n=twoF_n, twoF_s=twoF_s, off_c=off_c, off_f=off_f,
        fe_det=fe_det, skies=col(SIG, "sky"))


# ============================================================
# report
# ============================================================
def hr(t):
    print("\n" + "=" * 100 + f"\n{t}\n" + "=" * 100, flush=True)


def surface_census():
    """The census cell scored straight off SURFACE's OWN bank, with this file's estimator.

    An EXTERNAL consistency check, not an input to any number above. G-B2 measured that the
    C1 realisations re-drawn here carry a different GWB rotation from SURFACE's (max|dq| =
    0.51, max|d dlnL| = 4.1 nat at the same seed), so these two columns are NOT expected to
    agree realisation-by-realisation -- only distributionally. If they disagree by much more
    than the 30-realisation sampling scatter, something other than the rotation is moving."""
    sig = sorted(glob.glob(f"{SF_OUT}/sf_sig_h1275_T30_lit_k3_g*_n*.npz"))
    nul = sorted(glob.glob(f"{SF_OUT}/sf_nullN_h1275_T30_lit_k3_g*_n*.npz"))
    if not sig or len(nul) < 100:
        return None
    off = []
    for f in nul:
        d = np.load(f, allow_pickle=True)
        off.append(offender(np.nan_to_num(d["dlnL_det"], posinf=1e30), d["lnK"], d["qmax"]))
    fl, err, est, zf = adopt(np.array(off))
    cert = corr = wrong = det_real = 0
    for f in sig:
        d = np.load(f, allow_pickle=True)
        dd = np.nan_to_num(d["dlnL_det"], posinf=1e30)
        dt = dd > np.maximum(d["lnK"], fl)
        ct = dt & (d["qmax"] > QBAR)
        cert += int(ct.sum()); corr += int((ct & d["on_true"]).sum())
        wrong += int((ct & ~d["on_true"]).sum()); det_real += int(dt.any())
    n = len(sig)
    return dict(n=n, n_null=len(nul), floor=fl, err=err, est=est, zero_frac=zf,
                det_real=det_real / n, cert=cert / n, corr=corr / n, wrong=wrong / n,
                purity=(corr / cert if cert else np.nan))


def main():
    RES = {}
    for cell, label in CELLS:
        NUL, nsrc = load(cell, "null"), cell
        if NUL is None:                      # the POOLED null set -- licensed by G-B6 below
            nsrc = CELLS[0][0]
            NUL = load(nsrc, "null")
        SIG = load(cell, "sig")
        print(f"  [{cell}] nulls from {nsrc}"
              f"{'  (POOLED -- see G-B6)' if nsrc != cell else ''}")
        if not SIG or not NUL:
            print(f"  {cell}: incomplete (sig={len(SIG or [])}, null={len(NUL or [])}) "
                  "-- skipped", flush=True)
            continue
        RES[cell] = score_cell(cell, NUL, SIG)
        RES[cell]["label"] = label

    if not RES:
        print("no complete cell on disk"); return 1

    # ---- the cross-build pooling gate: same seeds, ncw=16 vs ncw=47 ----
    hr("GATE G-B6 -- the pure-noise null set is VENUE-ONLY (ncw=16 build == ncw=47 build)")
    XN = load("C2_askysurv_e03_h1275_k5", "xnull")
    N16 = load("C1_census_h1275_T30_lit_k3", "null")
    if XN and N16:
        m16 = {int(r["noise_seed"]): r for r in N16}
        d2f = d_dl = 0.0
        for r in XN:
            s = int(r["noise_seed"])
            if s in m16:
                d2f = max(d2f, abs(float(r["twoF"]) - float(m16[s]["twoF"])))
                a = np.nan_to_num(r["f_dlnL_det"], posinf=1e30)
                b = np.nan_to_num(m16[s]["f_dlnL_det"], posinf=1e30)
                d_dl = max(d_dl, float(np.max(np.abs(a - b))))
        ok = (d2f < 1e-6 and d_dl < 1e-6)
        print(f"  {len(XN)} seeds re-run in the ncw=47 build: max|d 2F| = {d2f:.3e}, "
              f"max|d dlnL_field| = {d_dl:.3e}  ->  {'PASS' if ok else '*** FAIL'}")
        print("  (a pure-noise draw contains no source, so the two builds must agree; "
              "this licenses ONE null set for both cells)")
    else:
        print("  xnull set absent -- pooling NOT licensed; each cell uses its own nulls")

    # ================= TABLE =================
    hr("BASELINE TABLE -- detection at matched FAP = 0.05, and certification purity")
    print("\n### B1 -- DETECTION (realisation level, both thresholds at FAP = 0.05 on the "
          "SAME nulls)\n")
    print(f"{'cell':46s} {'N_sig':>5s} {'N_null':>6s} | "
          f"{'Fe t95':>9s} {'est':>8s} {'zf':>5s} | {'Fe det':>16s} | "
          f"{'crit floor':>11s} {'est':>8s} {'zf':>5s} | {'crit det':>16s}")
    print("-" * 160)
    for c, R in RES.items():
        fe, cr = R["fe"], R["crit"]
        print(f"{R['label']:46s} {R['n_sig']:5d} {R['n_null']:6d} | "
              f"{fe['t95']:9.2f} {'emp_q95':>8s} {fe['zero_frac']:5.2f} | "
              f"{fe['det']['k']:3d}/{fe['det']['n']:<3d} = {fe['det']['f']:5.2f} "
              f"[{fe['det']['lo']:.2f},{fe['det']['hi']:.2f}] | "
              f"{cr['floor']:11.2f} {cr['est']:>8s} {cr['zero_frac']:5.2f} | "
              f"{cr['det']['k']:3d}/{cr['det']['n']:<3d} = {cr['det']['f']:5.2f} "
              f"[{cr['det']['lo']:.2f},{cr['det']['hi']:.2f}]")

    print("\n### B2 -- CERTIFICATION PURITY of confident picks (q_max > 0.9), on the "
          "Fe-DETECTED realisations\n")
    print(f"{'cell':46s} {'method':>22s} | {'conf/real':>9s} {'corr/real':>9s} "
          f"{'wrong/real':>10s} | {'purity':>7s} {'Wilson 95%':>16s} {'n_conf':>7s}")
    print("-" * 160)
    order = [("field", "FIELD (as practised)"), ("field_inwin", "FIELD, truth in window"),
             ("field_d1", "FIELD + D1 gate"),
             ("field_mc", "FIELD, oracle mc"),
             ("field_oracle", "FIELD, oracle source"),
             ("crit", "CRITERION v2.2"),
             ("crit_all", "CRITERION, ALL reals")]
    for c, R in RES.items():
        first = True
        for kk, nm in order:
            if kk not in R["pur"]:
                continue
            p = R["pur"][kk]
            pu = " n/a" if not np.isfinite(p["purity"]) else f"{p['purity']:4.2f}"
            ci = f"[{p['lo']:.2f},{p['hi']:.2f}]"
            print(f"{R['label'] if first else '':46s} {nm:>22s} | "
                  f"{p['per_real']:9.2f} {p['corr_per_real']:9.2f} "
                  f"{p['wrong_per_real']:10.2f} | {pu:>7s} {ci:>16s} {p['n_conf']:7d}")
            first = False
        print("-" * 160)
    print("  NOTE: 'CRITERION v2.2' is restricted to the Fe-DETECTED realisations, so it is "
          "read beside the field's on the SAME subset.\n"
          "        'CRITERION, ALL reals' is the criterion's own number over every "
          "realisation -- the one the programme quotes.")

    print("\n### floors, with the criterion-v2.2 validity column\n")
    for c, R in RES.items():
        ff = R["fieldfloor"]
        print(f"  {R['label']}")
        print(f"     FIELD per-pulsar floor  {ff['floor']:8.2f} +- {ff['err']:5.2f} nat  "
              f"[{ff['est']}, zero-fraction {ff['zero_frac']:.2f}]")
        print(f"     CRITERION floor         {R['crit']['floor']:8.2f} +- "
              f"{R['crit']['err']:5.2f} nat  [{R['crit']['est']}, zero-fraction "
              f"{R['crit']['zero_frac']:.2f}]")
        print(f"     FIELD Fe threshold      {R['fe']['t95']:8.2f} +- "
              f"{R['fe']['t95_sd']:5.2f} (2F)  [emp_q95, zero-fraction "
              f"{R['fe']['zero_frac']:.2f}; Gumbel {R['fe']['gumbel']:.2f} +- "
              f"{R['fe']['gumbel_sd']:.2f}]")
        print(f"     null 2F  median {R['fe']['null_med']:7.2f}  max {R['fe']['null_max']:7.2f}"
              f"   |  signal 2F  median {R['fe']['sig_med']:7.2f}  max "
              f"{R['fe']['sig_max']:7.2f}")
        if R["sep_deg"] is not None:
            d = R["sep_deg"][R["fe_det"]]
            if d.size:
                print(f"     Fe-detected sky error vs the nearest injected member: median "
                      f"{np.median(d):.1f} deg (min {d.min():.1f}, max {d.max():.1f})")
        print(f"     comb spacing (median over pulsars): FIELD "
              f"{np.median(R['f_dL_med']):.4f} kpc  vs  CRITERION "
              f"{np.median(R['c_dL_med']):.4f} kpc")
        po = R["polish"]
        print(f"     polish vs grid argmax: improved on {po['sig_rate']:.0%} of signal / "
              f"{po['null_rate']:.0%} of null reals; median gain "
              f"{po['sig_med']:+.2f} / {po['null_med']:+.2f} 2F, max {po['sig_max']:+.1f}; "
              f"WORST change {po['worst']:+.3e} (>= -1e-6 x 2F, the fusion tolerance)")
        ex = R["excl"]
        print(f"     single-fringe pulsars EXCLUDED (both arms, same rule): "
              f"FIELD {ex['field_sig']}/{ex['n_psr_sig']} signal, "
              f"{ex['field_null']}/{ex['n_psr_null']} null  |  "
              f"CRITERION {ex['crit_sig']}/{ex['n_psr_sig']} signal, "
              f"{ex['crit_null']}/{ex['n_psr_null']} null")
        print(f"     pulsar-realisations whose TRUE fringe was never sampled by the field's "
              f"512-pt window: {ex['field_truth_unsampled']}/{ex['n_psr_sig']}")
    # ---- external check: the census cell off SURFACE's own bank ----
    hr("EXTERNAL CHECK -- the census cell scored off SURFACE's OWN bank (different GWB "
       "rotation; distributional agreement only, see G-B2)")
    sc = surface_census()
    c1 = RES.get("C1_census_h1275_T30_lit_k3")
    if sc and c1:
        print(f"{'quantity':34s} {'BASELINE re-draw':>18s} {'SURFACE bank':>18s}")
        print("-" * 74)
        print(f"{'n signal / n null':34s} {c1['n_sig']:>8d} /{c1['n_null']:>7d} "
              f"{sc['n']:>10d} /{sc['n_null']:>6d}")
        print(f"{'criterion floor (nat)':34s} {c1['crit']['floor']:18.2f} "
              f"{sc['floor']:18.2f}")
        print(f"{'  its error / estimator / zero-f':34s} "
              f"{c1['crit']['err']:6.2f} {c1['crit']['est']:>7s} "
              f"{c1['crit']['zero_frac']:4.2f} "
              f"{sc['err']:8.2f} {sc['est']:>7s} {sc['zero_frac']:4.2f}")
        print(f"{'realisations detected':34s} {c1['crit']['det']['f']:18.2f} "
              f"{sc['det_real']:18.2f}")
        print(f"{'certifications / realisation':34s} "
              f"{c1['crit']['cert_per_real']:18.2f} {sc['cert']:18.2f}")
        print(f"{'CORRECT certs / realisation':34s} "
              f"{c1['crit']['corr_per_real']:18.2f} {sc['corr']:18.2f}")
        print(f"{'WRONG certs / realisation':34s} "
              f"{c1['crit']['wrong_per_real']:18.2f} {sc['wrong']:18.2f}")
        print("\n  These two columns are DIFFERENT REALISATIONS of the same cell (G-B2), so "
              "they are read for\n  distributional agreement only. Nothing above depends on "
              "this table.")
    else:
        print("  SURFACE bank not reachable, or C1 incomplete -- external check skipped")

    np.savez(f"{REPORTS}/baseline_score.npz",
             **{f"{c}__{k}": np.asarray(v, dtype=object) if isinstance(v, dict) else v
                for c, R in RES.items() for k, v in R.items() if v is not None})
    with open(f"{REPORTS}/BASELINE_tables.txt", "w") as fh:
        fh.write(json.dumps({c: {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                                 for k, v in R.items()
                                 if k in ("label", "n_sig", "n_null", "fe", "crit",
                                          "fieldfloor", "pur")}
                             for c, R in RES.items()}, indent=1, default=float))
    print(f"\nsaved {REPORTS}/baseline_score.npz and {REPORTS}/BASELINE_tables.txt")
    panel(RES)
    return 0


def panel(RES):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))
    cells = list(RES)
    x = np.arange(len(cells))

    # ---- A: detection at matched FAP ----
    ax = axes[0]
    series = [("Fe-statistic (field)", CFIELD, lambda R: R["fe"]["det"]),
              ("criterion-v2.2", CCRIT, lambda R: R["crit"]["det"])]
    for i, (nm, cl, get) in enumerate(series):
        dx = (i - (len(series) - 1) / 2) * 0.18
        v = [get(RES[c]) for c in cells]
        ax.errorbar(x + dx, [d["f"] for d in v],
                    yerr=[[d["f"] - d["lo"] for d in v], [d["hi"] - d["f"] for d in v]],
                    fmt="o", ms=8, lw=2, color=cl, capsize=4, label=nm, zorder=3)
        for xi, d in zip(x + dx, v):
            up = d["f"] < 0.85
            ax.annotate(f"{d['f']:.2f}", (xi, d["hi"] if up else d["lo"]),
                        textcoords="offset points", xytext=(0, 7 if up else -15),
                        ha="center", fontsize=8.5, color="0.25")
    ax.set_ylim(-0.06, 1.10)
    ax.set_ylabel("fraction of realisations detected")
    ax.set_title("A — detection at matched FAP = 0.05\n(same 100 pure-noise nulls set both "
                 "thresholds)", fontsize=10)

    # ---- B: purity of confident picks ----
    ax = axes[1]
    pser = [("FIELD (as practised)", CFIELD, "field"),
            ("FIELD + D1 gate", CFIELD_D1, "field_d1"),
            ("CRITERION v2.2", CCRIT, "crit")]
    for i, (nm, cl, kk) in enumerate(pser):
        dx = (i - (len(pser) - 1) / 2) * 0.16
        v = [RES[c]["pur"].get(kk) for c in cells]
        yy = [(p["purity"] if p and np.isfinite(p["purity"]) else np.nan) for p in v]
        lo = [(p["purity"] - p["lo"] if p and np.isfinite(p["purity"]) else 0) for p in v]
        hi = [(p["hi"] - p["purity"] if p and np.isfinite(p["purity"]) else 0) for p in v]
        ax.errorbar(x + dx, yy, yerr=[lo, hi], fmt="o", ms=8, lw=2, color=cl,
                    capsize=4, label=nm, zorder=3)
        for xi, p, y in zip(x + dx, v, yy):
            if p and np.isfinite(y):
                up = y < 0.80
                ax.annotate(f"{y:.2f}\nn={p['n_conf']}",
                            (xi, (p["hi"] if up else p["lo"])),
                            textcoords="offset points", xytext=(0, 7 if up else -22),
                            ha="center", fontsize=7.5, color="0.25")
    ax.axhline(1.0, color="0.35", lw=1.1, ls="--", zorder=1)
    ax.set_ylim(-0.06, 1.10)
    ax.set_ylabel("purity  =  correct / confident  (q$_{max}$ > 0.9)")
    ax.set_title("B — certification purity on Fe-detected realisations\n"
                 "(the same bar, applied to each method's own comb)", fontsize=10)

    # legends go BELOW the axes: every in-axes corner collides with a marker or a whisker
    # at some data value, and this figure has to survive whatever the fan returns.
    for ax in axes:
        ax.set_xticks(x, [RES[c]["label"].replace(", ", "\n", 1) for c in cells],
                      fontsize=8.5)
        ax.set_xlim(-0.55, len(cells) - 0.45)
        ax.grid(axis="y", color="0.9", lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8.5, loc="upper center",
                  bbox_to_anchor=(0.5, -0.17), ncol=3, columnspacing=1.2,
                  handletextpad=0.4)
    fig.text(0.5, -0.06, "dashed line in B: purity = 1 (every confident pick on the true "
             "fringe).  Error bars are Wilson 95 % intervals; n is the number of confident "
             "picks behind each point.", ha="center", fontsize=7.8, color="0.35")
    fig.tight_layout()
    p = f"{REPORTS}/BASELINE_panel.png"
    fig.savefig(p, dpi=180, bbox_inches="tight")
    print(f"  panel -> {p}")


if __name__ == "__main__":
    sys.exit(main())
