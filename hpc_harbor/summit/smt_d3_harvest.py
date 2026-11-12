#!/usr/bin/env python
"""PHOENIX -- SUMMIT S1.4: the D3 (N_psr) single-dial ladder readout.

The D3 fan (12769215 warm + 12769216 main) ran to completion on 2026-07-25/26 and was
never harvested -- the session died with the banks on disk. This is that harvest.

RUNGS:  116 (INHERITED from the banked GLACIER ladder, never re-run)
        +30  -> 146   SUMMIT_results/d3_ext30/
        +100 -> 216   SUMMIT_results/d3_ext100/
Additions are clones of the top-N COMPASS-quality pulsars at seeded golden-angle
positions, THE BUNDLE INTACT (TOAs, design matrix, EM distance, K-table row, frozen
reference RN all the donor's; only name+position differ). SEED_D3 = 9300.

Scored on the SAME three contours as the plane, with SUMMIT's pre-registered HOLD bar:
>=2 consecutive refit intervals in which every commonly-fed member moves <=0.05 dex in
BOTH (log10 fgw, log10 mc).

CPU only.
"""
import glob, re, os
import numpy as np

NP_SRC, I_MC, I_FGW = 8, 3, 4
HOLD_DEX = 0.05
BAR = 0.5

SRC = [("116 (inherited)", "GLACIER_results", "gl2_r13p9_cell_none_s{s}_T30_lit"),
       ("146 (+30)",  "SUMMIT_results/d3_ext30",  "gl2_r13p9_cell_none_s{s}_T30_lit"),
       ("216 (+100)", "SUMMIT_results/d3_ext100", "gl2_r13p9_cell_none_s{s}_T30_lit")]
NULLS = [("216 (+100)", "SUMMIT_results/d3_ext100",
          ["gl2_r13p9_null0_none_s0_T30_lit", "gl2_r13p9_null0_none_s1_T30_lit",
           "gl2_r13p9_null1_none_s0_T30_lit"])]


def load_cell(d, stem):
    fs = sorted(glob.glob(f"{d}/{stem}_i*__*.npz"),
                key=lambda f: int(re.search(r"_i(\d+)__", f).group(1)))
    fs = [f for f in fs if "STOPANAT" not in f]
    if not fs:
        return None
    its = [int(re.search(r"_i(\d+)__", f).group(1)) for f in fs]
    return its, [np.load(f, allow_pickle=True) for f in fs]


def score(its, Z):
    nres = np.array([int(z["n_resolved"]) for z in Z])
    ncert = np.array([int(z["n_cert"]) for z in Z])
    wrong = np.array([int(z["wrong_cert"]) for z in Z])
    abg = np.array([float(z["a_bg"]) for z in Z])
    aeff = float(Z[0]["a_eff_drawn"])
    nd = [len(np.asarray(z["theta_rec"], float)) - NP_SRC * int(z["n_slot"]) for z in Z]
    fedth = []
    for z, n in zip(Z, nd):
        th = np.asarray(z["theta_rec"], float)
        fedth.append({int(k): (th[n + NP_SRC*k + I_FGW], th[n + NP_SRC*k + I_MC])
                      for k in np.flatnonzero(z["fed_mask"])})
    wmax, ncom = [], []
    for a, b in zip(fedth[:-1], fedth[1:]):
        com = set(a) & set(b)
        ncom.append(len(com))
        wmax.append(np.nan if not com else
                    max(max(abs(b[k][0]-a[k][0]), abs(b[k][1]-a[k][1])) for k in com))
    wmax = np.array(wmax, float)
    st = np.array([(not np.isnan(x)) and x <= HOLD_DEX and n > 0 for x, n in zip(wmax, ncom)], bool)
    hold = any(st[j] and st[j+1] for j in range(len(st)-1)) if len(st) > 1 else False
    return dict(n_it=len(its), npsr=int(nd[0]), feed=bool((nres > 0).any()),
                cert=bool((ncert > 0).any()), hold=bool(hold),
                nres_max=int(nres.max()), ncert_max=int(ncert.max()),
                wrong=int(wrong.sum()), bite=abg[0]-aeff, regr=abg[-1]-aeff,
                wander=float(np.nanmax(wmax)) if not np.all(np.isnan(wmax)) else np.nan,
                floor0=float(Z[0]["floor"]), zf0=float(Z[0]["zero_fraction"]))


print(f"{'rung':16s}{'sky':>4s}{'npsr':>6s}{'it':>4s}{'FEED':>6s}{'HOLD':>6s}{'CERT':>6s}"
      f"{'nres':>6s}{'wrong':>6s}{'bite_i0':>9s}{'regr_end':>9s}{'wanderMax':>10s}")
summary = {}
for label, d, tmpl in SRC:
    rows = []
    for s in range(4):
        r = load_cell(d, tmpl.format(s=s))
        if r is None:
            print(f"{label:16s}{s:>4d}{'':>6s}   -- MISSING --")
            continue
        v = score(*r)
        rows.append(v)
        print(f"{label:16s}{s:>4d}{v['npsr']:>6d}{v['n_it']:>4d}{int(v['feed']):>6d}"
              f"{int(v['hold']):>6d}{int(v['cert']):>6d}{v['nres_max']:>6d}{v['wrong']:>6d}"
              f"{v['bite']:>9.3f}{v['regr']:>9.3f}{v['wander']:>10.3f}")
    if rows:
        summary[label] = rows
    print()

print("NULL ARM (scrambled counterpart, the manufacturing control):")
for label, d, stems in NULLS:
    for st in stems:
        r = load_cell(d, st)
        if r is None:
            print(f"  {st}  -- MISSING (check for a STOPANAT bank)")
            continue
        v = score(*r)
        print(f"  {st:42s} it={v['n_it']} feed={int(v['feed'])} cert={int(v['cert'])} "
              f"wrong={v['wrong']} nres={v['nres_max']}")
    sa = sorted(glob.glob(f"{d}/*STOPANAT*"))
    for f in sa:
        print(f"  STOPANAT bank present: {os.path.basename(f)}")

print("\n" + "=" * 92)
print(f"{'rung':16s}{'feed/4':>8s}{'hold/4':>8s}{'cert/4':>8s}{'wrong':>7s}"
      f"{'mean bite':>11s}{'mean regr':>11s}{'mean wander':>13s}")
for label, rows in summary.items():
    print(f"{label:16s}{sum(r['feed'] for r in rows):>8d}{sum(r['hold'] for r in rows):>8d}"
          f"{sum(r['cert'] for r in rows):>8d}{sum(r['wrong'] for r in rows):>7d}"
          f"{np.mean([r['bite'] for r in rows]):>11.3f}"
          f"{np.mean([r['regr'] for r in rows]):>11.3f}"
          f"{np.nanmean([r['wander'] for r in rows]):>13.3f}")
print("=" * 92)

np.savez("reports/smt_d3_ladder.npz",
         rung=np.array([l for l, rs in summary.items() for _ in rs]),
         sky=np.array([i for l, rs in summary.items() for i, _ in enumerate(rs)]),
         **{k: np.array([r[k] for l, rs in summary.items() for r in rs])
            for k in ("npsr", "feed", "hold", "cert", "nres_max", "ncert_max",
                      "wrong", "bite", "regr", "wander", "floor0", "zf0")},
         note="SUMMIT D3 (N_psr) single-dial ladder at the feasible rung r13p9, circular "
              "arm, T=30, w=1. 116 inherited from the GLACIER ladder; +30/+100 are "
              "COMPASS-rank clones at golden-angle positions, SEED_D3=9300, bundle "
              "intact. HOLD = >=2 consecutive intervals, all commonly-fed members "
              "<=0.05 dex in (fgw,mc).")
print("banked reports/smt_d3_ladder.npz")
