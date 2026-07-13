"""RECUT — CHORUS re-scored under the ADOPTED floor convention (ANCHOR §4).

All 26 CHORUS cells have a nullN zero-fraction above the 20% validity gate, so EVERY CHORUS
floor is invalid as a Gumbel and must be restated as the empirical (1-alpha) quantile with a
bootstrap error. FLOOR_FIX_provisional corrected the FLOORS (ch_floors banked its raw
offenders) but could not re-cut the COUNTS. The per-realisation signal banks
(CHORUS_results/ch_sig_*) are on this machine, so this script re-cuts them.

THE HEADLINE THIS SETTLES: does the e = 0.3 switch-on survive its floor rise?
  m1 e=0.3 lit : floor 7.39 -> 11.30 nat (+53%), banked count 1.57, margin +0.57 over the >1 bar
  m1 e=0.3 vlbi: floor 10.58 -> 10.78 nat (+2%), banked count 1.13, margin +0.13
The floor rises, so the count can only fall. Whether it falls THROUGH the bar is what the raw
banks decide, and nothing else can.

READBACK GATES (no corrected number is emitted until an uncorrected one is reproduced):
  GATE A  floors : banked raw offenders (offN_i) -> banked fN, fN_sd, zero_frac   [exact]
  GATE B  counts : this scorer at the OLD floors -> banked surface_corr/_wrong    [exact]

TWO BARS ARE REPORTED, and they are different questions:
  count > 1                      the campaign's own switch-on bar (CHORUS §5 "the >1 bar")
  count at (floor + err) > 1     the criterion-v2.1 ONSET grade: survives the floor's own
                                 calibration error. Under the convention `err` is the
                                 BOOTSTRAP sd of the empirical quantile.
A cell is CONFIRMED only if it clears the second. Clearing only the first is MARGINAL.

CPU only. No GPU, no jax, no new realisations. Run:
    python CW_transition/recut_chorus.py
"""
import os, sys, glob
import numpy as np
from scipy.stats import gumbel_r

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "CHORUS_results")
REPORTS = os.path.join(ROOT, "reports")

ALPHA = 0.05
Z_ALPHA = -np.log(-np.log(1.0 - ALPHA))
C_SD = 2.80
QBAR = 0.9
ZF_GATE = 0.20
BOOT = 4000
SEED = 20260713

E_ORDER = {"e00": 0.0, "e03": 0.3, "e05": 0.5, "e07": 0.7, "eU": -1.0}   # eU = uniform draw


def _cells():
    mixes = [(0, "e00")] + [(k, ed) for k in (1, 2, 3) for ed in ("e03", "e05", "e07", "eU")]
    return [(k, ed, t) for (k, ed) in mixes for t in ("lit", "vlbi")]


def cell_tag(k, ed, t):
    return f"m{k}{ed}_{t}"


def gumbel_floor(x):
    x = np.asarray(x, float)
    mu, beta = gumbel_r.fit(x)
    return (float(mu + beta * Z_ALPHA), float(mu), float(beta),
            float(C_SD * beta / np.sqrt(len(x))), len(x))


def emp_quantile(off, alpha=ALPHA):
    return float(np.quantile(np.asarray(off, float), 1.0 - alpha))


def boot_sd(off, alpha=ALPHA, B=BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    off = np.asarray(off, float)
    n = len(off)
    return float(np.std([np.quantile(rng.choice(off, n, replace=True), 1.0 - alpha)
                         for _ in range(B)]))


def adopt(gu, gu_sd, off, zf):
    if zf <= ZF_GATE:
        return float(gu), float(gu_sd), "gumbel"
    return emp_quantile(off), boot_sd(off), "emp_q95"


def load_sig(tag):
    """The per-realisation signal bank. CHORUS's correctness column is (mapk == n_true_grid)."""
    fs = sorted(glob.glob(f"{OUT}/ch_sig_{tag}_g*_n*.npz"))
    D, K, Q, M, N = [], [], [], [], []
    for f in fs:
        d = np.load(f, allow_pickle=True)
        D.append(np.nan_to_num(np.asarray(d["dlnL_det"]), posinf=1e30))
        K.append(np.asarray(d["lnK"])); Q.append(np.asarray(d["qmax"]))
        M.append(np.asarray(d["mapk"])); N.append(np.asarray(d["n_true_grid"]))
    return dict(dlnL=np.array(D), lnK=np.array(K), qmax=np.array(Q),
                mapk=np.array(M), ntg=np.array(N), n=len(fs))


def score(S, floor):
    """chorus_analysis.stage1_surface, verbatim: cert & (mapk == ntg), averaged per realisation."""
    cert = (S["dlnL"] > np.maximum(S["lnK"], floor)) & (S["qmax"] > QBAR)
    corr = cert & (S["mapk"] == S["ntg"])
    wrong = cert & (S["mapk"] != S["ntg"])
    return float(corr.sum()) / S["n"], float(wrong.sum()) / S["n"]


def hr(t):
    print("\n" + "=" * 108 + f"\n{t}\n" + "=" * 108, flush=True)


def main():
    F = np.load(os.path.join(REPORTS, "ch_floors.npz"), allow_pickle=True)
    A = np.load(os.path.join(REPORTS, "ch_analysis.npz"), allow_pickle=True)
    cells = [str(c) for c in F["cells"]]
    tags = [str(t) for t in A["surface_tags"]]
    b_corr, b_wrong = A["surface_corr"], A["surface_wrong"]

    as_tag = lambda c: (lambda p: f"m{p[0]}{p[1]}_{p[2]}")(c.split(","))
    assert [as_tag(c) for c in cells] == tags, "ch_floors / ch_analysis cell order disagrees"

    # ================= GATE A =================
    hr("GATE A — the banked raw offenders must reproduce the banked Gumbel floors")
    dfn = dsd = dzf = 0.0
    for i in range(len(cells)):
        o = F[f"offN_{i}"]
        fN, mu, beta, sd, n = gumbel_floor(o)
        dfn = max(dfn, abs(fN - float(F[f"fN_{i}"])))
        dsd = max(dsd, abs(sd - float(F[f"fN_sd_{i}"])))
        dzf = max(dzf, abs(float((o == 0.0).mean()) - float(F[f"zero_frac_{i}"])))
    print(f"  max |recomputed - banked|  fN        = {dfn:.3e}")
    print(f"  max |recomputed - banked|  fN_sd     = {dsd:.3e}")
    print(f"  max |recomputed - banked|  zero_frac = {dzf:.3e}")
    if max(dfn, dsd, dzf) > 1e-9:
        print("  *** STOP: banked offenders do not reproduce the banked floors.")
        return 1
    print("  GATE A PASSED.")

    # ================= GATE B =================
    hr("GATE B — this scorer at the OLD (Gumbel) floors must reproduce the banked counts")
    SIG = {}
    dc = dw = 0.0
    for i, c in enumerate(cells):
        k, ed, t = c.split(",")
        tag = cell_tag(int(k), ed, t)
        SIG[i] = load_sig(tag)
        cc, ww = score(SIG[i], float(F[f"fN_{i}"]))
        dc = max(dc, abs(cc - float(b_corr[i])))
        dw = max(dw, abs(ww - float(b_wrong[i])))
    print(f"  max |recomputed - banked|  corr  = {dc:.3e}")
    print(f"  max |recomputed - banked|  wrong = {dw:.3e}")
    if max(dc, dw) > 1e-9:
        print("  *** STOP: the scorer here is not the banked scorer.")
        return 1
    print("  GATE B PASSED — the scorer here IS the one that produced the published counts.")

    # ================= THE RE-CUT =================
    hr("THE RE-CUT — all 26 cells: Gumbel | empirical q95 + bootstrap | ADOPTED, counts re-scored")
    rows = []
    for i, c in enumerate(cells):
        k, ed, t = c.split(",")
        o = F[f"offN_{i}"]
        zf = float(F[f"zero_frac_{i}"])
        gu, gsd = float(F[f"fN_{i}"]), float(F[f"fN_sd_{i}"])
        emp, esd = emp_quantile(o), boot_sd(o)
        fl, err, est = adopt(gu, gsd, o, zf)
        cc, ww = score(SIG[i], fl)
        clo, _ = score(SIG[i], fl + err)          # the ONSET grade: survives the floor's error
        chi, _ = score(SIG[i], max(fl - err, 0.0))
        grade = "CONFIRMED" if clo > 1.0 else ("MARGINAL" if cc > 1.0 else "below")
        rows.append(dict(cell=c, tag=tags[i], k=int(k), ed=ed, tier=t, n=SIG[i]["n"], zf=zf,
                         gumbel=gu, gumbel_sd=gsd, emp=emp, emp_sd=esd, floor=fl, err=err,
                         est=est, corr=cc, corr_lo=clo, corr_hi=chi, wrong=ww,
                         old_corr=float(b_corr[i]), old_wrong=float(b_wrong[i]),
                         grade=grade))

    n_inv = sum(1 for r in rows if r["est"] != "gumbel")
    n_up = sum(1 for r in rows if r["floor"] > r["gumbel"])
    print(f"  zero-fraction > {ZF_GATE:.0%} -> Gumbel INVALID : {n_inv} / {len(rows)}")
    print(f"     of which the floor RISES (count can only fall) : {n_up}")
    print()
    print(f"  {'cell':>12} | {'zero-f':>6} | {'Gumbel':>7} {'±fit':>5} | {'ADOPTED':>8} {'±bs':>5} "
          f"{'est':>7} | {'ratio':>6} | {'old':>5} {'NEW':>5} {'@fl+bs':>6} {'wrong':>5} | grade")
    for r in rows:
        print(f"  {r['cell']:>12} | {r['zf']:6.2f} | {r['gumbel']:7.2f} {r['gumbel_sd']:5.2f} | "
              f"{r['floor']:8.2f} {r['err']:5.2f} {r['est']:>7} | "
              f"{r['floor'] / r['gumbel']:5.2f}x | {r['old_corr']:5.2f} {r['corr']:5.2f} "
              f"{r['corr_lo']:6.2f} {r['wrong']:5.2f} | {r['grade']}")

    # ================= THE HEADLINE =================
    hr("THE HEADLINE — does the e = 0.3 switch-on survive its floor rise?")
    for r in rows:
        if r["k"] == 1 and r["ed"] in ("e03", "e05"):
            surv = "SURVIVES" if r["corr"] > 1.0 else "DOES NOT SURVIVE"
            print(f"  {r['cell']:>12}: floor {r['gumbel']:6.2f} -> {r['floor']:6.2f} ± {r['err']:.2f} "
                  f"({r['floor'] / r['gumbel']:.2f}x, "
                  f"{abs(r['floor'] - r['gumbel']) / r['gumbel_sd']:.1f} sigma of its Gumbel fit err)")
            print(f"  {'':>12}  count {r['old_corr']:.2f} -> {r['corr']:.2f}  "
                  f"[at floor+bs: {r['corr_lo']:.2f}]   margin over the >1 bar "
                  f"{r['corr'] - 1.0:+.2f}   ==> {surv}  ({r['grade']})")

    # ---- the corrected switch-on threshold in e ----
    hr("THE CORRECTED SWITCH-ON THRESHOLD IN e (the externally quotable number)")
    print(f"  Rule: the switch is ON at eccentricity e for a given (n_ecc, tier) when the "
          f"re-cut count\n  exceeds the >1 bar. CONFIRMED additionally requires it at floor + "
          f"bootstrap error.\n")
    print(f"  {'n_ecc':>5} {'tier':>5} | " + " ".join(f"{e:>22}" for e in ("e=0.3", "e=0.5", "e=0.7")))
    thr = {}
    for k in (1, 2, 3):
        for t in ("lit", "vlbi"):
            cs = []
            on_e = None
            for ed in ("e03", "e05", "e07"):
                r = next(x for x in rows if (x["k"], x["ed"], x["tier"]) == (k, ed, t))
                mark = "ON " if r["corr"] > 1.0 else "off"
                cs.append(f"{mark} {r['corr']:5.2f}[{r['corr_lo']:5.2f}] {r['grade'][:4]:>4}")
                if r["corr"] > 1.0 and on_e is None:
                    on_e = E_ORDER[ed]
            thr[(k, t)] = on_e
            print(f"  {k:5d} {t:>5} | " + " ".join(f"{c:>22}" for c in cs))
    print()
    for k in (1, 2, 3):
        for t in ("lit", "vlbi"):
            e = thr[(k, t)]
            print(f"    n_ecc={k} {t:>4}: switch-on at e = "
                  f"{e if e is not None else 'NOT REACHED in {0.3,0.5,0.7}'}")

    e03 = {r["tier"]: r for r in rows if r["k"] == 1 and r["ed"] == "e03"}
    both = all(e03[t]["corr"] > 1.0 for t in ("lit", "vlbi"))
    conf = all(e03[t]["grade"] == "CONFIRMED" for t in ("lit", "vlbi"))
    print()
    if both and conf:
        print("  ==> THE SINGLE e = 0.3 MEMBER STILL SWITCHES THE COUNT ON, in BOTH tiers, and it")
        print("      survives the corrected floor AND that floor's bootstrap error.")
        print("      The corrected switch-on threshold is e = 0.3. The docs' interim binding to")
        print("      e = 0.5 is LIFTED.")
    elif both:
        print("  ==> The e = 0.3 switch-on SURVIVES the >1 bar in both tiers, but does not clear")
        print("      it at floor + bootstrap error in every tier (see grades). Quote e = 0.3 with")
        print("      the grade attached.")
    else:
        off = [t for t in ("lit", "vlbi") if e03[t]["corr"] <= 1.0]
        print(f"  ==> The e = 0.3 switch-on DOES NOT SURVIVE in tier(s): {off}.")
        print("      The corrected switch-on threshold is e = 0.5, and the docs' interim binding")
        print("      to e = 0.5 STANDS.")

    keys = ["k", "n", "zf", "gumbel", "gumbel_sd", "emp", "emp_sd", "floor", "err",
            "corr", "corr_lo", "corr_hi", "wrong", "old_corr", "old_wrong"]
    np.savez(
        os.path.join(REPORTS, "recut_chorus.npz"),
        alpha=ALPHA, zf_gate=ZF_GATE, boot=BOOT, seed=SEED, qbar=QBAR,
        cols=np.array(keys),
        table=np.array([[r[c] for c in keys] for r in rows], float),
        cells=np.array([r["cell"] for r in rows]),
        tags=np.array([r["tag"] for r in rows]),
        ecc=np.array([r["ed"] for r in rows]),
        tiers=np.array([r["tier"] for r in rows]),
        estimator=np.array([r["est"] for r in rows]),
        grade=np.array([r["grade"] for r in rows]),
        n_invalid=n_inv, n_floor_up=n_up,
        gateA_maxdev=float(max(dfn, dsd, dzf)), gateB_maxdev=float(max(dc, dw)),
        note=("counts RE-CUT from per-realisation ch_sig banks against the ADOPTED floor "
              "(empirical q95 + bootstrap sd wherever zero-fraction > 0.20 -- which is all 26 "
              "cells); grade CONFIRMED iff count at (floor + bootstrap err) > 1"))
    print(f"\nbanked -> reports/recut_chorus.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
