"""
KINDLE g2 / D-7 — the fALL column, re-cut against the criterion-v2.2 validity gate.

THE LAST UNFINISHED PIECE OF THE FLOOR FIX (project_progress.md §10.16(g), §10.17 open-2).

The claim under test (SURFACE §10.13(f)): "21 cells clear onset on fALL, best 2.57/real,
and all 21 are 5+11." It stands on the pre-fix Gumbel estimator that criterion-v2.2 bounds,
because `surface_floors.npz` banks NO `fALL_zerofrac` -- so the validity gate was never
applied to it. RECUT banked the nullN offenders only.

RECOVERABILITY (the finding that makes this closable):
  surface.py:23    fALL := Gumbel over  nullN + nullA + nullL,  N >= 200 per cell
  surface.py:76    offender(dlnL, lnK, qmax) -- a PURE FUNCTION of three columns
  every sf_null{N,A,L}_*.npz banks dlnL_det, lnK, qmax
=> the fALL offender vectors were never banked but are FULLY RECOMPUTABLE from disk.
D-7 therefore closes as RECOVERED, not as UNRECOVERABLE-from-disk.

Gates, mirroring CW_transition/recut_surface.py exactly. No corrected number is emitted
until an uncorrected one is reproduced from the same raw columns:
  GATE A-ALL  recomputed oALL reproduces banked fALL, fALL_mu, fALL_beta, fALL_sd,
              fALL_n, fALL_emp                                   -> must be 0.0e+00
  GATE B-ALL  this scorer at the OLD fALL floors reproduces banked corr_fALL /
              corr_fALL_lo                                        -> must be 0.0e+00

CPU only. No GPU, no jax, no new realisations.
"""
import os, glob
import numpy as np
from scipy.stats import gumbel_r

# this file lives at <repo>/hpc_harbor/kindle/, so the repo root is THREE levels up
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "SURFACE_results")
REPORTS = os.path.join(ROOT, "reports")
KOUT = os.path.join(ROOT, "KINDLE_results")
os.makedirs(KOUT, exist_ok=True)

# ---- criterion-v2.1/2.2 constants, verbatim from surface.py / recut_surface.py ----
ALPHA = 0.05
Z_ALPHA = -np.log(-np.log(1.0 - ALPHA))
C_SD = 2.80
QBAR = 0.9
H_GRID = [-13.25, -13.0, -12.75, -12.5, -12.25, -12.0]
T_GRID = [30, 40, 50]
TIERS = ["lit", "vlbi"]
STRUCTS = [3, 5, 2]
STRUCT_LABEL = {3: "3+13", 5: "5+11", 2: "2+14"}

ZF_GATE = 0.20          # the criterion-v2.2 floor-validity gate
BOOT = 4000
SEED = 20260713         # same fixed seed as recut_surface.py: the bootstrap is reproducible


def hkey(h):
    return f"{abs(h):.2f}".replace(".", "")


def cell_tag(h, T, tier, k):
    return f"h{hkey(h)}_T{T}_{tier}_k{k}"


def cells():
    return [(h, T, tier, k) for T in T_GRID for k in STRUCTS
            for tier in TIERS for h in H_GRID]


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
    rng = np.random.default_rng(seed)
    off = np.asarray(off, float)
    n = len(off)
    q = [np.quantile(rng.choice(off, n, replace=True), 1.0 - alpha) for _ in range(B)]
    return float(np.std(q))


def adopt(gu, gu_sd, off, zf):
    """(floor, err, estimator) under the criterion-v2.2 gate."""
    if zf <= ZF_GATE:
        return float(gu), float(gu_sd), "gumbel"
    return emp_quantile(off), boot_sd(off), "emp_q95"


def load_signal(h, T, tier, k):
    fs = sorted(glob.glob(f"{OUT}/sf_sig_{cell_tag(h, T, tier, k)}_g*_n*.npz"))
    if not fs:
        return None
    D, K, Q, OT = [], [], [], []
    for f in fs:
        d = np.load(f, allow_pickle=True)
        D.append(np.nan_to_num(d["dlnL_det"], posinf=1e30))
        K.append(d["lnK"]); Q.append(d["qmax"]); OT.append(d["on_true"])
    return dict(dlnL=np.array(D), lnK=np.array(K), qmax=np.array(Q),
                on_true=np.array(OT), n=len(fs))


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
    # ============ STEP 1 — recompute the fALL offender vectors ============
    hr("STEP 1 — recompute the fALL offender vectors (nullN + nullA + nullL) from the raw banks")
    OFF, NPARTS = {}, {}
    for c in cells():
        tag = cell_tag(*c)
        parts = {}
        for var in ("nullN", "nullA", "nullL"):
            o = []
            for f in sorted(glob.glob(f"{OUT}/sf_{var}_{tag}_g*_n*.npz")):
                d = np.load(f, allow_pickle=True)
                o.append(offender(np.nan_to_num(d["dlnL_det"], posinf=1e30),
                                  d["lnK"], d["qmax"]))
            parts[var] = np.array(o, float)
        OFF[c] = np.concatenate([parts["nullN"], parts["nullA"], parts["nullL"]])
        NPARTS[c] = {k: len(v) for k, v in parts.items()}
    ns = np.array([len(OFF[c]) for c in cells()])
    ex = NPARTS[cells()[0]]
    print(f"  {len(OFF)} cells; per cell nullN/nullA/nullL = "
          f"{ex['nullN']}/{ex['nullA']}/{ex['nullL']}")
    print(f"  oALL per cell: min {ns.min()}  max {ns.max()}   (surface.py requires N >= 200)")
    if ns.min() < 200:
        print("  *** STOP: a cell has fewer than 200 fALL nulls.")
        return 1

    # ============ GATE A-ALL — floors reproduce ============
    hr("GATE A-ALL — recomputed oALL must reproduce the BANKED fALL floors (surface_floors.npz)")
    FB = np.load(os.path.join(REPORTS, "surface_floors.npz"), allow_pickle=True)
    fcols = [str(c) for c in FB["tab_cols"]]
    ftab, ftiers = FB["table"], [str(t) for t in FB["tiers"]]
    bank = {}
    for i in range(len(ftab)):
        r = {c: ftab[i, j] for j, c in enumerate(fcols)}
        bank[(float(r["h"]), int(r["T"]), ftiers[i], int(r["k"]))] = r

    dev = {"fALL": 0.0, "fALL_mu": 0.0, "fALL_beta": 0.0,
           "fALL_sd": 0.0, "fALL_n": 0.0, "fALL_emp": 0.0}
    for c in cells():
        o = OFF[c]
        fA, mu, beta, sd, n = gumbel_floor(o)
        got = dict(fALL=fA, fALL_mu=mu, fALL_beta=beta, fALL_sd=sd,
                   fALL_n=float(n), fALL_emp=emp_quantile(o))
        for kk in dev:
            dev[kk] = max(dev[kk], abs(got[kk] - float(bank[c][kk])))
    for kk, v in dev.items():
        print(f"  max |recomputed - banked|  {kk:12s} = {v:.3e}")
    if max(dev.values()) > 1e-9:
        print("  *** STOP: recomputed fALL offenders do not reproduce the bank. The offender "
              "statistic or the null family here is NOT the one that produced the surface.")
        return 1
    print("  GATE A-ALL PASSED — these ARE the vectors the published fALL floors were fitted to.")

    # ============ GATE B-ALL — counts reproduce under the OLD fALL floors ============
    hr("GATE B-ALL — this scorer at the OLD fALL floors must reproduce the BANKED corr_fALL")
    AB = np.load(os.path.join(REPORTS, "surface_analysis.npz"), allow_pickle=True)
    acols = [str(c) for c in AB["cols"]]
    atab = AB["table"]
    atiers = np.array([str(t) for t in AB["tiers"]])
    key2row = {(float(atab[i, acols.index("h")]), int(atab[i, acols.index("T")]),
                atiers[i], int(atab[i, acols.index("k")])): i for i in range(len(atab))}
    g = lambda c: atab[:, acols.index(c)]

    SIG = {}
    dA = dAlo = 0.0
    for c in cells():
        S = load_signal(*c)
        if S is None:
            print(f"  *** STOP: no signal bank for {cell_tag(*c)}")
            return 1
        SIG[c] = S
        b, i = bank[c], key2row[c]
        fA, sd = float(b["fALL"]), float(b["fALL_sd"])
        _, _, corr, _ = score(S, fA)
        _, _, corr_lo, _ = score(S, fA + sd)
        dA = max(dA, abs(per_real(corr, S["n"]) - g("corr_fALL")[i]))
        dAlo = max(dAlo, abs(per_real(corr_lo, S["n"]) - g("corr_fALL_lo")[i]))
    print(f"  max |recomputed - banked|  corr_fALL    = {dA:.3e}")
    print(f"  max |recomputed - banked|  corr_fALL_lo = {dAlo:.3e}")
    if max(dA, dAlo) > 1e-9:
        print("  *** STOP: the scorer here is not the banked scorer. No corrected fALL count "
              "may be emitted.")
        return 1
    print("  GATE B-ALL PASSED — the scorer here IS the one that produced the published fALL "
          "column.")

    # ============ STEP 2 — THE RE-CUT ============
    hr("STEP 2 — THE fALL RE-CUT: zero-fraction computed, validity gate applied, counts re-scored")
    rows = []
    for c in cells():
        h, T, tier, k = c
        o, S, b, i = OFF[c], SIG[c], bank[c], key2row[c]
        zf = float((o == 0.0).mean())
        gu, gu_sd = float(b["fALL"]), float(b["fALL_sd"])
        fl, err, est = adopt(gu, gu_sd, o, zf)
        _, _, corr, wrong = score(S, fl)
        _, _, corr_lo, _ = score(S, fl + err)
        cA, cAlo = per_real(corr, S["n"]), per_real(corr_lo, S["n"])
        wA = per_real(wrong, S["n"])
        old, old_lo = float(g("corr_fALL")[i]), float(g("corr_fALL_lo")[i])
        rows.append(dict(h=h, T=T, tier=tier, k=k, zf=zf, gumbel=gu, gumbel_sd=gu_sd,
                         emp=emp_quantile(o), floor=fl, err=err, est=est,
                         old=old, old_lo=old_lo, corr=cA, corr_lo=cAlo, wrong=wA,
                         old_verdict=verdict_of(old, old_lo),
                         new_verdict=verdict_of(cA, cAlo)))

    zfs = np.array([r["zf"] for r in rows])
    touched = [r for r in rows if r["zf"] > ZF_GATE]
    print(f"  fALL zero-fraction across 108 cells: min {zfs.min():.2f}  median "
          f"{np.median(zfs):.2f}  max {zfs.max():.2f}")
    print(f"  cells FAILING the validity gate (zf > {ZF_GATE:.0%}), Gumbel INVALID: "
          f"{len(touched)} / 108")
    print(f"  cells passing (Gumbel valid, floor untouched)                      : "
          f"{108 - len(touched)} / 108")

    # ---- the claim under test ----
    hr("STEP 3 — ADJUDICATION: 'the fALL column ignites: 21 cells, all of them 5+11'")
    old_on = [r for r in rows if r["old_verdict"] == "ONSET"]
    new_on = [r for r in rows if r["new_verdict"] == "ONSET"]
    old_marg = [r for r in rows if r["old_verdict"] == "MARGINAL"]
    new_marg = [r for r in rows if r["new_verdict"] == "MARGINAL"]
    print(f"  PUBLISHED (pre-fix Gumbel fALL floors) : ONSET {len(old_on):3d}   "
          f"MARGINAL {len(old_marg):3d}   below {108 - len(old_on) - len(old_marg):3d}")
    print(f"  RE-CUT    (criterion-v2.2 fALL floors) : ONSET {len(new_on):3d}   "
          f"MARGINAL {len(new_marg):3d}   below {108 - len(new_on) - len(new_marg):3d}")
    print()
    for lab, grp in (("PUBLISHED", old_on), ("RE-CUT", new_on)):
        bs = {}
        for r in grp:
            bs[STRUCT_LABEL[r["k"]]] = bs.get(STRUCT_LABEL[r["k"]], 0) + 1
        print(f"  {lab:10s} fALL-ONSET cells by structure: {bs if bs else '{}'}")

    flipped = [r for r in rows if r["old_verdict"] != r["new_verdict"]]
    print(f"\n  cells whose fALL VERDICT changed: {len(flipped)}")
    if flipped:
        print(f"\n  {'h':>7s} {'T':>3s} {'tier':>5s} {'struct':>7s} {'zf':>5s} "
              f"{'Gumbel':>8s} {'ADOPTED':>16s} {'old':>6s} {'new':>6s}  was -> now")
        for r in sorted(flipped, key=lambda r: (r["k"], r["h"])):
            print(f"  {r['h']:7.2f} {r['T']:3d} {r['tier']:>5s} {STRUCT_LABEL[r['k']]:>7s} "
                  f"{r['zf']:5.2f} {r['gumbel']:8.2f} {r['floor']:9.2f} +-{r['err']:4.2f} "
                  f"{r['old']:6.2f} {r['corr']:6.2f}  {r['old_verdict']} -> {r['new_verdict']}")

    if new_on:
        print(f"\n  the surviving fALL-ONSET cells (best first):")
        print(f"  {'h':>7s} {'T':>3s} {'tier':>5s} {'struct':>7s} {'zf':>5s} "
              f"{'ADOPTED floor':>16s} {'est':>7s} {'corr':>6s} {'@fl+e':>6s} {'wrong':>6s}")
        for r in sorted(new_on, key=lambda r: -r["corr"]):
            print(f"  {r['h']:7.2f} {r['T']:3d} {r['tier']:>5s} {STRUCT_LABEL[r['k']]:>7s} "
                  f"{r['zf']:5.2f} {r['floor']:9.2f} +-{r['err']:4.2f} {r['est']:>7s} "
                  f"{r['corr']:6.2f} {r['corr_lo']:6.2f} {r['wrong']:6.2f}")

    # ============ BANK ============
    hr("STEP 4 — BANK the fALL offenders (orientation DECLARED, per §10.16(d)) + the re-cut")
    cl = cells()
    off_mat = np.array([OFF[c] for c in cl], float)          # (108, 200)
    index = np.array([f"{c[0]:.2f}|{c[1]}|{c[2]}|{STRUCT_LABEL[c[3]]}" for c in cl])
    orientation = ("off_i <-> index[i] <-> meta row i. NO transpose, NO implied permutation. "
                   "Row i of `off` is cell i of every metadata column (h, T, tier, k, "
                   "zero_frac, emp_q95). Columns 0:100 are nullN, 100:150 nullA, 150:200 nullL, "
                   "concatenated exactly as surface.py builds oALL.")
    np.savez_compressed(
        os.path.join(REPORTS, "surface_fALL_offenders.npz"),
        off=off_mat, index=index, orientation=orientation,
        h=np.array([c[0] for c in cl]), T=np.array([c[1] for c in cl]),
        tier=np.array([c[2] for c in cl]), k=np.array([c[3] for c in cl]),
        n=np.array([len(OFF[c]) for c in cl]),
        n_nullN=np.array([NPARTS[c]["nullN"] for c in cl]),
        n_nullA=np.array([NPARTS[c]["nullA"] for c in cl]),
        n_nullL=np.array([NPARTS[c]["nullL"] for c in cl]),
        zero_frac=np.array([float((OFF[c] == 0.0).mean()) for c in cl]),
        emp_q95=np.array([emp_quantile(OFF[c]) for c in cl]),
    )
    print(f"  reports/surface_fALL_offenders.npz   ({off_mat.shape[0]} cells x "
          f"{off_mat.shape[1]} nulls, orientation declared)")

    keys = ["h", "T", "k", "zf", "gumbel", "gumbel_sd", "emp", "floor", "err",
            "old", "old_lo", "corr", "corr_lo", "wrong"]
    tab = np.array([[float(r[k]) for k in keys] for r in rows])
    np.savez_compressed(
        os.path.join(KOUT, "kindle_recut_fALL.npz"),
        table=tab, cols=np.array(keys),
        tiers=np.array([r["tier"] for r in rows]),
        est=np.array([r["est"] for r in rows]),
        old_verdict=np.array([r["old_verdict"] for r in rows]),
        new_verdict=np.array([r["new_verdict"] for r in rows]),
        zf_gate=ZF_GATE, alpha=ALPHA, boot=BOOT, seed=SEED,
    )
    print(f"  KINDLE_results/kindle_recut_fALL.npz  (108 cells, raw floors + counts + verdicts)")

    # readback gate
    Z = np.load(os.path.join(REPORTS, "surface_fALL_offenders.npz"), allow_pickle=True)
    rb = max(float(np.abs(Z["off"][i] - OFF[cl[i]]).max()) for i in range(len(cl)))
    rbz = max(abs(float(Z["zero_frac"][i]) - float((OFF[cl[i]] == 0.0).mean()))
              for i in range(len(cl)))
    print(f"\n  READBACK GATE  max|off - recomputed| = {rb:.3e}   "
          f"max|zero_frac - recomputed| = {rbz:.3e}")
    if max(rb, rbz) > 0.0:
        print("  *** STOP: the banked offenders do not read back.")
        return 1
    print("  READBACK GATE PASSED (0.0e+00) — the banked vectors ARE the ones re-cut here.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
