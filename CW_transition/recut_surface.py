"""RECUT — SURFACE re-scored under the ADOPTED floor convention (ANCHOR §4).

THE CONVENTION (adopted verbatim, ANCHOR 2026-07-13 §4):

    The D2 Gumbel-MLE floor is VALID ONLY where the nullN zero-fraction is <= 20%.
    Above that, the floor is the empirical (1 - alpha) quantile with a BOOTSTRAP error.
    The zero-fraction is a REQUIRED column beside every floor, not a caveat.

What FLOOR_FIX_provisional could NOT do, this does: it RE-CUTS THE COUNT. The per-realisation
signal banks (SURFACE_results/sf_sig_*) are on this machine, so every count below is scored
from the raw statistics against the ADOPTED floor -- not bounded, not interpolated.

THE READBACK GATES (this script refuses to emit a corrected number until it has reproduced
an uncorrected one from the same raw columns):

  GATE A  floors : recomputed nullN offenders  -> banked fN, fN_sd, fN_emp, fN_zerofrac
                   in surface_floors.npz                                    [must be exact]
  GATE B  counts : this scorer at the OLD floors -> banked corr, corr_lo, corr_hi, verdicts
                   in surface_analysis.npz                                  [must be exact]

Only if A and B pass is the adopted-floor re-cut trustworthy: they prove the offender
statistic, the floor estimator, and the scorer here are the SAME code paths that produced
the published surface.

Verdict rule is UNCHANGED from surface_analysis.py -- only the floor and its error change:
    ONSET    corr at (floor + err) > 1     (survives the floor's own 1-sigma calibration error)
    MARGINAL corr at  floor       > 1
    below    otherwise
Under the convention, for a touched cell (zf > 20%) `floor` is the empirical q95 and `err` is
its BOOTSTRAP sd -- so the onset test is made against a bootstrap error, exactly as required.

CPU only. No GPU, no jax, no new realisations. Run:
    python CW_transition/recut_surface.py
"""
import os, sys, glob
import numpy as np
from scipy.stats import gumbel_r

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "SURFACE_results")
REPORTS = os.path.join(ROOT, "reports")

# ---- criterion-v2.1 constants, copied verbatim from hpc_harbor/surface/surface.py ----
ALPHA = 0.05
Z_ALPHA = -np.log(-np.log(1.0 - ALPHA))      # 2.9702
C_SD = 2.80                                   # sd(floor_hat) = C_SD * beta / sqrt(N)   (G7)
QBAR = 0.9
QBAR_STRICT = 0.99
H_GRID = [-13.25, -13.0, -12.75, -12.5, -12.25, -12.0]
T_GRID = [30, 40, 50]
TIERS = ["lit", "vlbi"]
STRUCTS = [3, 5, 2]
STRUCT_LABEL = {3: "3+13", 5: "5+11", 2: "2+14"}

# ---- the convention ----
ZF_GATE = 0.20
BOOT = 4000
SEED = 20260713                               # fixed: the bootstrap is reproducible


def hkey(h):
    return f"{abs(h):.2f}".replace(".", "")


def cell_tag(h, T, tier, k):
    return f"h{hkey(h)}_T{T}_{tier}_k{k}"


def cells():
    return [(h, T, tier, k) for T in T_GRID for k in STRUCTS
            for tier in TIERS for h in H_GRID]


# ---- the offender statistic (surface.py, verbatim) ----
def offender(dlnL, lnK, qmax):
    """Largest dlnL among pulsars passing layers 1+3. MAX OVER PULSARS -> Gumbel domain."""
    m = (dlnL > lnK) & (qmax > QBAR)
    return float(dlnL[m].max()) if m.any() else 0.0


def gumbel_floor(x):
    x = np.asarray(x, float)
    mu, beta = gumbel_r.fit(x)
    return (float(mu + beta * Z_ALPHA), float(mu), float(beta),
            float(C_SD * beta / np.sqrt(len(x))), len(x))


def emp_quantile(off, alpha=ALPHA):
    return float(np.quantile(np.asarray(off, float), 1.0 - alpha))


def boot_sd(off, alpha=ALPHA, B=BOOT, seed=SEED):
    """Bootstrap standard error of the empirical quantile. Valid under a zero point mass."""
    rng = np.random.default_rng(seed)
    off = np.asarray(off, float)
    n = len(off)
    q = [np.quantile(rng.choice(off, n, replace=True), 1.0 - alpha) for _ in range(B)]
    return float(np.std(q))


def adopt(gu, gu_sd, off, zf):
    """(floor, err, estimator) under ANCHOR §4."""
    if zf <= ZF_GATE:
        return float(gu), float(gu_sd), "gumbel"
    return emp_quantile(off), boot_sd(off), "emp_q95"


# ---- the scorer (surface_analysis.py, verbatim) ----
def load_signal(h, T, tier, k):
    fs = sorted(glob.glob(f"{OUT}/sf_sig_{cell_tag(h, T, tier, k)}_g*_n*.npz"))
    if not fs:
        return None
    D, K, Q, OT = [], [], [], []
    names = None
    for f in fs:
        d = np.load(f, allow_pickle=True)
        D.append(np.nan_to_num(d["dlnL_det"], posinf=1e30))
        K.append(d["lnK"]); Q.append(d["qmax"]); OT.append(d["on_true"])
        names = [str(x) for x in d["names"]]
    return dict(dlnL=np.array(D), lnK=np.array(K), qmax=np.array(Q),
                on_true=np.array(OT), n=len(fs), names=names)


def score(S, floor, bar=QBAR):
    det = S["dlnL"] > np.maximum(S["lnK"], floor)
    cert = det & (S["qmax"] > bar)
    return det, cert, cert & S["on_true"], cert & ~S["on_true"]


def per_real(mask, n):
    return float(mask.sum()) / n


def verdict_of(corr, corr_lo):
    if corr_lo > 1.0:
        return "ONSET"
    return "MARGINAL" if corr > 1.0 else "below"


def hr(t):
    print("\n" + "=" * 104 + f"\n{t}\n" + "=" * 104, flush=True)


def main():
    # ================= recompute the nullN offenders from the raw banks =================
    hr("STEP 1 — recompute the nullN offender vectors from the raw SURFACE_results banks")
    OFF, meta = {}, []
    for (h, T, tier, k) in cells():
        tag = cell_tag(h, T, tier, k)
        fs = sorted(glob.glob(f"{OUT}/sf_nullN_{tag}_g*_n*.npz"))
        o = []
        for f in fs:
            d = np.load(f, allow_pickle=True)
            o.append(offender(np.nan_to_num(d["dlnL_det"], posinf=1e30), d["lnK"], d["qmax"]))
        OFF[(h, T, tier, k)] = np.array(o, float)
        meta.append((h, T, tier, k, len(o)))
    ns = np.array([m[4] for m in meta])
    print(f"  {len(OFF)} cells; nullN per cell: min {ns.min()}  max {ns.max()}  "
          f"(D2 rule requires N >= 100)")
    if ns.min() < 100:
        print("  *** STOP: a cell has fewer than 100 nullN realisations.")
        return 1

    # ================= GATE A: floors reproduce =================
    hr("GATE A — recomputed offenders must reproduce the BANKED floors (surface_floors.npz)")
    FB = np.load(os.path.join(REPORTS, "surface_floors.npz"), allow_pickle=True)
    fcols = [str(c) for c in FB["tab_cols"]]
    ftab, ftiers = FB["table"], [str(t) for t in FB["tiers"]]
    bank = {}
    for i in range(len(ftab)):
        r = {c: ftab[i, j] for j, c in enumerate(fcols)}
        bank[(float(r["h"]), int(r["T"]), ftiers[i], int(r["k"]))] = r

    dev = {"fN": 0.0, "fN_sd": 0.0, "fN_emp": 0.0, "fN_zerofrac": 0.0}
    for c in OFF:
        o = OFF[c]
        fN, mu, beta, sd, n = gumbel_floor(o)
        got = dict(fN=fN, fN_sd=sd, fN_emp=emp_quantile(o), fN_zerofrac=float((o == 0.0).mean()))
        for kk in dev:
            dev[kk] = max(dev[kk], abs(got[kk] - float(bank[c][kk])))
    for kk, v in dev.items():
        print(f"  max |recomputed - banked|  {kk:12s} = {v:.3e}")
    if max(dev.values()) > 1e-9:
        print("  *** STOP: recomputed floors do not reproduce the bank. The offender statistic "
              "or the estimator here is NOT the one that produced the published surface.")
        return 1
    print("  GATE A PASSED — the offender statistic and the D2 estimator here are the banked ones.")

    # ================= GATE B: counts reproduce under the OLD floors =================
    hr("GATE B — this scorer at the OLD (Gumbel) floors must reproduce the BANKED counts")
    AB = np.load(os.path.join(REPORTS, "surface_analysis.npz"), allow_pickle=True)
    acols = [str(c) for c in AB["cols"]]
    atab = AB["table"]
    atiers = np.array([str(t) for t in AB["tiers"]])
    averd = np.array([str(v) for v in AB["verdicts"]])
    g = lambda c: atab[:, acols.index(c)]
    key2row = {(float(atab[i, acols.index("h")]), int(atab[i, acols.index("T")]),
                atiers[i], int(atab[i, acols.index("k")])): i for i in range(len(atab))}

    SIG = {}
    dcorr = dlo = dhi = 0.0
    nv = 0
    for c in cells():
        S = load_signal(*c)
        if S is None:
            print(f"  *** STOP: no signal bank for {cell_tag(*c)}")
            return 1
        SIG[c] = S
        i = key2row[c]
        b = bank[c]
        fN, sd = float(b["fN"]), float(b["fN_sd"])
        got = {}
        for lab, fl in (("", fN), ("_lo", fN + sd), ("_hi", fN - sd)):
            _, _, corr, _ = score(S, fl)
            got[lab] = per_real(corr, S["n"])
        dcorr = max(dcorr, abs(got[""] - g("corr")[i]))
        dlo = max(dlo, abs(got["_lo"] - g("corr_lo")[i]))
        dhi = max(dhi, abs(got["_hi"] - g("corr_hi")[i]))
        nv += (verdict_of(got[""], got["_lo"]) == averd[i])
    print(f"  max |recomputed - banked|  corr    = {dcorr:.3e}")
    print(f"  max |recomputed - banked|  corr_lo = {dlo:.3e}")
    print(f"  max |recomputed - banked|  corr_hi = {dhi:.3e}")
    print(f"  verdicts reproduced: {nv} / {len(cells())}")
    if max(dcorr, dlo, dhi) > 1e-9 or nv != len(cells()):
        print("  *** STOP: the scorer here is not the banked scorer. No corrected count may be "
              "emitted.")
        return 1
    print("  GATE B PASSED — the scorer here IS the one that produced the published counts.")
    print(f"  (banked: ONSET {int((averd == 'ONSET').sum())}  "
          f"MARGINAL {int((averd == 'MARGINAL').sum())}  below {int((averd == 'below').sum())})")

    # ================= THE RE-CUT =================
    hr("STEP 2 — THE RE-CUT: adopted floor (Gumbel where zf <= 20%, empirical q95 + bootstrap "
       "above), counts re-scored from the raw signal banks")
    rows = []
    for c in cells():
        h, T, tier, k = c
        o, S, b, i = OFF[c], SIG[c], bank[c], key2row[c]
        zf = float(b["fN_zerofrac"])
        gu, gu_sd = float(b["fN"]), float(b["fN_sd"])
        fl, err, est = adopt(gu, gu_sd, o, zf)

        _, cert, corr, wrong = score(S, fl)
        _, _, corr_lo, _ = score(S, fl + err)
        _, _, corr_hi, _ = score(S, max(fl - err, 0.0))
        n = S["n"]
        cc, clo, chi = per_real(corr, n), per_real(corr_lo, n), per_real(corr_hi, n)
        vnew = verdict_of(cc, clo)
        rows.append(dict(
            h=h, T=T, tier=tier, k=k, n=n, zf=zf,
            gumbel=gu, gumbel_sd=gu_sd, emp=emp_quantile(o), emp_sd=boot_sd(o),
            floor=fl, err=err, est=est, touched=zf > ZF_GATE,
            corr=cc, corr_lo=clo, corr_hi=chi, wrong=per_real(wrong, n),
            corr_se=float(corr.sum(axis=1).std(ddof=1) / np.sqrt(n)),
            verdict=vnew,
            old_corr=float(g("corr")[i]), old_corr_lo=float(g("corr_lo")[i]),
            old_verdict=str(averd[i])))

    touched = [r for r in rows if r["touched"]]
    print(f"  cells touched (zero-fraction > {ZF_GATE:.0%}, Gumbel INVALID) : "
          f"{len(touched)} / {len(rows)}")
    print(f"  cells untouched (Gumbel valid, floor and count stand)         : "
          f"{len(rows) - len(touched)} / {len(rows)}")

    print("\n  THE 15 TOUCHED CELLS — floors and RE-CUT counts (the entire blast radius):")
    print(f"  {'h':>7} {'T':>3} {'tier':>5} {'struct':>6} | {'zero-f':>6} | {'Gumbel':>7} "
          f"{'EMPq95':>7} {'±bs':>5} {'ratio':>6} | {'old':>5} {'NEW':>5} {'@fl+bs':>6} | "
          f"{'was':>8} -> {'NOW':<8}")
    for r in sorted(touched, key=lambda x: -x["zf"]):
        flip = "" if r["verdict"] == r["old_verdict"] else "   <<< FLIP"
        print(f"  {r['h']:7.2f} {r['T']:3d} {r['tier']:>5} {STRUCT_LABEL[r['k']]:>6} | "
              f"{r['zf']:6.2f} | {r['gumbel']:7.2f} {r['emp']:7.2f} {r['emp_sd']:5.2f} "
              f"{r['emp'] / r['gumbel']:5.2f}x | {r['old_corr']:5.2f} {r['corr']:5.2f} "
              f"{r['corr_lo']:6.2f} | {r['old_verdict']:>8} -> {r['verdict']:<8}{flip}")

    # ---- the flips ----
    flips = [r for r in rows if r["verdict"] != r["old_verdict"]]
    on_new = [r for r in rows if r["verdict"] == "ONSET"]
    on_old = [r for r in rows if r["old_verdict"] == "ONSET"]
    lost = [r for r in rows if r["old_verdict"] == "ONSET" and r["verdict"] != "ONSET"]
    gained = [r for r in rows if r["old_verdict"] != "ONSET" and r["verdict"] == "ONSET"]

    hr("STEP 3 — THE VERDICT: the corrected onset count")
    print(f"  banked (pre-fix)   N_onset = {len(on_old)}")
    print(f"  provisional bound            57 <= N_onset <= 67   (FLOOR_FIX_provisional §2)")
    print(f"  RE-CUT             N_onset = {len(on_new)}")
    print(f"\n  onsets LOST  (were ONSET, now not) : {len(lost)}")
    for r in lost:
        print(f"     h={r['h']:6.2f} T={r['T']:2d} {r['tier']:4s} {STRUCT_LABEL[r['k']]:5s}: "
              f"floor {r['gumbel']:.2f} -> {r['floor']:.2f} ({r['floor'] / r['gumbel']:.2f}x), "
              f"count {r['old_corr']:.2f} -> {r['corr']:.2f} [lo {r['corr_lo']:.2f}] "
              f"=> {r['verdict']}")
    print(f"\n  onsets GAINED (were not, now ONSET): {len(gained)}")
    for r in gained:
        print(f"     h={r['h']:6.2f} T={r['T']:2d} {r['tier']:4s} {STRUCT_LABEL[r['k']]:5s}: "
              f"floor {r['gumbel']:.2f} -> {r['floor']:.2f} ({r['floor'] / r['gumbel']:.2f}x), "
              f"count {r['old_corr']:.2f} -> {r['corr']:.2f} [lo {r['corr_lo']:.2f}] "
              f"(was {r['old_verdict']})")
    print(f"\n  all verdict flips (any direction): {len(flips)}")

    # ---- the two at-risk cells named in the brief ----
    hr("STEP 4 — THE TWO AT-RISK ONSET CELLS (FLOOR_FIX_provisional §2), settled")
    AT_RISK = [(-13.25, 40, "lit", 3), (-13.25, 50, "lit", 2)]
    for c in AT_RISK:
        r = next(x for x in rows if (x["h"], x["T"], x["tier"], x["k"]) == c)
        st = "SURVIVES" if r["verdict"] == "ONSET" else "RETRACTED"
        print(f"  h={r['h']:.2f} T={r['T']} {r['tier']} {STRUCT_LABEL[r['k']]}: "
              f"floor {r['gumbel']:.2f} -> {r['floor']:.2f} ± {r['err']:.2f} "
              f"({r['floor'] / r['gumbel']:.2f}x), zf {r['zf']:.2f}")
        print(f"      count {r['old_corr']:.2f} -> {r['corr']:.2f}  "
              f"[at floor+bs: {r['corr_lo']:.2f}]   provisional bound said <= "
              f"{r['old_corr_lo']:.2f}")
        print(f"      ==> {r['verdict']}  ({st})")

    # ---- the faint frontier / h* claim ----
    hr("STEP 5 — THE FAINT FRONTIER: is h* still unbounded below? (SURFACE §4 ‡, suspended)")
    HMIN = min(H_GRID)
    print(f"  A column = (tier, structure, T). h* is UNBOUNDED BELOW in that column when the "
          f"faintest\n  grid strain h = {HMIN} is ONSET there (the frontier sits on the grid edge).\n")
    print(f"  {'tier':>5} {'struct':>6} {'T':>3} | {'old h* range':>16} {'edge?':>6} | "
          f"{'NEW h* range':>16} {'edge?':>6}")
    cols_old = cols_new = 0
    colrows = []
    for tier in TIERS:
        for k in STRUCTS:
            for T in T_GRID:
                oo = [r for r in rows if (r["tier"], r["k"], r["T"]) == (tier, k, T)
                      and r["old_verdict"] == "ONSET"]
                nn = [r for r in rows if (r["tier"], r["k"], r["T"]) == (tier, k, T)
                      and r["verdict"] == "ONSET"]
                eo = any(abs(r["h"] - HMIN) < 1e-9 for r in oo)
                en = any(abs(r["h"] - HMIN) < 1e-9 for r in nn)
                cols_old += eo; cols_new += en
                so = (f"{max(r['h'] for r in oo):.2f}..{min(r['h'] for r in oo):.2f}"
                      if oo else "none")
                sn = (f"{max(r['h'] for r in nn):.2f}..{min(r['h'] for r in nn):.2f}"
                      if nn else "none")
                print(f"  {tier:>5} {STRUCT_LABEL[k]:>6} {T:3d} | {so:>16} "
                      f"{'YES' if eo else '-':>6} | {sn:>16} {'YES' if en else '-':>6}")
                colrows.append((tier, k, T, eo, en))
    print(f"\n  columns with h* on the faint edge (unbounded below):  "
          f"OLD {cols_old} / 18   ->   RE-CUT {cols_new} / 18")
    if cols_new == 0:
        print("  ==> the ‡ claim is RETRACTED: with corrected floors no column's onset reaches the\n"
              "      faint grid edge. h* is BOUNDED BELOW everywhere in this box.")
    elif cols_new < cols_old:
        print(f"  ==> the ‡ claim is REINSTATED but WEAKENED: {cols_new} of 18 columns "
              f"(was {cols_old}).")
    else:
        print(f"  ==> the ‡ claim is REINSTATED as published: {cols_new} of 18 columns.")

    # ================= bank =================
    keys = ["h", "T", "k", "n", "zf", "gumbel", "gumbel_sd", "emp", "emp_sd", "floor", "err",
            "corr", "corr_lo", "corr_hi", "wrong", "corr_se", "old_corr", "old_corr_lo"]
    np.savez(
        os.path.join(REPORTS, "recut_surface.npz"),
        alpha=ALPHA, zf_gate=ZF_GATE, boot=BOOT, seed=SEED, qbar=QBAR,
        cols=np.array(keys),
        table=np.array([[r[c] for c in keys] for r in rows], float),
        tiers=np.array([r["tier"] for r in rows]),
        struct=np.array([STRUCT_LABEL[r["k"]] for r in rows]),
        estimator=np.array([r["est"] for r in rows]),
        touched=np.array([r["touched"] for r in rows]),
        verdict=np.array([r["verdict"] for r in rows]),
        old_verdict=np.array([r["old_verdict"] for r in rows]),
        n_onset=len(on_new), n_onset_banked=len(on_old),
        n_lost=len(lost), n_gained=len(gained),
        hstar_edge_cols_old=cols_old, hstar_edge_cols_new=cols_new,
        gateA_maxdev=float(max(dev.values())),
        gateB_maxdev=float(max(dcorr, dlo, dhi)),
        note=("counts RE-CUT from per-realisation sf_sig banks against the ADOPTED floor; "
              "floor = Gumbel where zero-fraction <= 0.20, else empirical q95 with bootstrap "
              "sd; ONSET iff count at (floor + err) > 1"))
    print(f"\nbanked -> reports/recut_surface.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
