"""Bank SURFACE's raw nullN offender vectors — the one campaign that never did.

WHY. The adopted floor convention (ANCHOR §4) says: where the nullN zero-fraction exceeds
20%, the floor is the empirical (1-alpha) quantile WITH A BOOTSTRAP ERROR. A bootstrap error
cannot be computed from a summary — it needs the raw offender vector. CHORUS, IGNITE-2 and
ANCHOR all banked theirs; SURFACE banked only `fN_emp` and `fN_zerofrac`, so its corrected
floors had no error bar and its onset test could not be run at floor + err.

This file supplies them, for all 108 cells, recomputed from SURFACE_results/sf_nullN_*.

ORIENTATION IS DECLARED, NOT IMPLIED. The ANCHOR ladder bank taught this lesson the hard
way: a raw offender array whose row order does not match its sibling metadata columns is a
silent trap. So this bank ships an explicit `index` column — `index[i]` is the label of
`off_i`, in the SAME order as every metadata column here — and `orientation` states the
convention in words. A re-cut that trusts row i marries the right cell's label by
construction.

Lean: one 1-D float64 vector of 100 offenders per cell (108 x 100 doubles ~ 86 kB).

Run: python CW_transition/bank_surface_offenders.py      (CPU, seconds)
"""
import os, sys, glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "SURFACE_results")
REPORTS = os.path.join(ROOT, "reports")

ALPHA = 0.05
QBAR = 0.9
H_GRID = [-13.25, -13.0, -12.75, -12.5, -12.25, -12.0]
T_GRID = [30, 40, 50]
TIERS = ["lit", "vlbi"]
STRUCTS = [3, 5, 2]
STRUCT_LABEL = {3: "3+13", 5: "5+11", 2: "2+14"}


def hkey(h):
    return f"{abs(h):.2f}".replace(".", "")


def cell_tag(h, T, tier, k):
    return f"h{hkey(h)}_T{T}_{tier}_k{k}"


def cells():
    return [(h, T, tier, k) for T in T_GRID for k in STRUCTS
            for tier in TIERS for h in H_GRID]


def offender(dlnL, lnK, qmax):
    """surface.py, verbatim: largest dlnL among pulsars passing layers 1+3, else 0.0."""
    m = (dlnL > lnK) & (qmax > QBAR)
    return float(dlnL[m].max()) if m.any() else 0.0


def main():
    print("=" * 96)
    print("BANKING SURFACE's raw nullN offender vectors (108 cells)")
    print("=" * 96)

    out, index, meta = {}, [], []
    for i, (h, T, tier, k) in enumerate(cells()):
        tag = cell_tag(h, T, tier, k)
        fs = sorted(glob.glob(f"{OUT}/sf_nullN_{tag}_g*_n*.npz"))
        o = []
        for f in fs:
            d = np.load(f, allow_pickle=True)
            o.append(offender(np.nan_to_num(d["dlnL_det"], posinf=1e30), d["lnK"], d["qmax"]))
        o = np.array(o, float)
        out[f"off_{i}"] = o
        index.append(tag)
        meta.append((h, T, k, len(o), float((o == 0.0).mean()),
                     float(np.quantile(o, 1 - ALPHA))))
        if len(o) < 100:
            print(f"  *** STOP: {tag} has only {len(o)} nullN realisations (D2 needs >= 100)")
            return 1

    meta = np.array(meta, float)
    tiers = np.array([c[2] for c in cells()])

    # ---- the readback gate: these offenders must reproduce the BANKED floor summaries ----
    FB = np.load(os.path.join(REPORTS, "surface_floors.npz"), allow_pickle=True)
    fcols = [str(c) for c in FB["tab_cols"]]
    ftab, ftiers = FB["table"], [str(t) for t in FB["tiers"]]
    bank = {}
    for i in range(len(ftab)):
        r = {c: ftab[i, j] for j, c in enumerate(fcols)}
        bank[(float(r["h"]), int(r["T"]), ftiers[i], int(r["k"]))] = r
    dzf = demp = 0.0
    for i, c in enumerate(cells()):
        dzf = max(dzf, abs(meta[i, 4] - float(bank[c]["fN_zerofrac"])))
        demp = max(demp, abs(meta[i, 5] - float(bank[c]["fN_emp"])))
    print(f"  readback gate vs surface_floors.npz:  max|dev| zero-frac = {dzf:.3e}   "
          f"max|dev| emp_q95 = {demp:.3e}")
    if max(dzf, demp) > 1e-9:
        print("  *** STOP: recomputed offenders do not reproduce the banked floor summaries.")
        return 1
    print("  PASSED — these ARE the vectors the published floors were fitted to.")

    np.savez(
        os.path.join(REPORTS, "surface_nullN_offenders.npz"),
        alpha=ALPHA, qbar=QBAR, n_cells=len(index),
        index=np.array(index),
        orientation=("row i of every column here (index, h, T, k, tier, n, zero_frac, emp_q95) "
                     "labels the vector stored under key 'off_i'. One vector per cell, in "
                     "cells() order: T outer, then struct k, then tier, then h. There is NO "
                     "transpose and no implied permutation: off_i <-> index[i] <-> meta row i."),
        h=meta[:, 0], T=meta[:, 1].astype(int), k=meta[:, 2].astype(int), tier=tiers,
        struct=np.array([STRUCT_LABEL[int(x)] for x in meta[:, 2]]),
        n=meta[:, 3].astype(int), zero_frac=meta[:, 4], emp_q95=meta[:, 5],
        statistic=("offender = max over pulsars of dlnL_det among those passing "
                   "layer1 (dlnL > lnK) AND layer3 (qmax > 0.9); 0.0 if none pass"),
        source="SURFACE_results/sf_nullN_<cell>_g*_n*.npz",
        **out)

    print(f"\n  banked -> reports/surface_nullN_offenders.npz  "
          f"({len(index)} cells x {int(meta[0,3])} offenders)")
    print(f"  zero-fraction range: {meta[:,4].min():.2f} .. {meta[:,4].max():.2f}   "
          f"({int((meta[:,4] > 0.20).sum())} cells above the 20% validity gate)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
