#!/usr/bin/env python
"""SUMMIT Phase-0 harvest -- the Stage-2a (loudness x capability) plane from the
GLACIER ladder banks. CPU-only, reads GLACIER_results/gl2_* (+ gl1_* Stage-1 rerun
for the (1,30) Stage-1 context) and writes reports/summit_s0_harvest.npz.

Three per-cell states (the three contours):
  FIRST-BITE  : any iteration with n_resolved > 0
  CERT        : any iteration with n_cert > 0
  HOLD        : >=1 fed member AND >=2 consecutive refit intervals where every
                member fed across the interval moves <= 0.05 dex in BOTH
                (log10_fgw, log10_mc)   [the pre-registered M-step trust bar]
Safety ledger: wrong_cert summed over iterations; null-arm columns kept separate
(control first-bite baseline under the S4.15.1 amended semantics).
"""
import re, sys, glob
import numpy as np

RESULTS = "GLACIER_results"
NP_SRC, I_MC, I_FGW, ND = 8, 3, 4, 116
HOLD_DEX = 0.05

pat = re.compile(
    r"gl2_(?P<rung>r13p9|r13p5|r13p25)(?P<w>_w0p\d+)?_(?P<kind>cell|null\d+)_(?P<arm>e07|none)"
    r"_s(?P<sky>\d+)_T(?P<T>\d+)_(?P<tier>lit|vlbi|subpc)_i(?P<it>\d+)__(?P<lane>.+)\.npz")

cells = {}
for f in sorted(glob.glob(f"{RESULTS}/gl2_*_i*.npz") +
                glob.glob("SUMMIT_results/gl2_*_i*.npz")):
    m = pat.search(f)
    if not m:
        continue
    w = m.group("w") or "_w1"
    key = (m.group("rung"), m.group("kind"), m.group("arm"), int(m.group("sky")),
           w.lstrip("_"), int(m.group("T")), m.group("tier"))
    cells.setdefault(key, {})[int(m.group("it"))] = f

def fed_theta(z):
    fed = np.flatnonzero(z["fed_mask"])
    th = np.asarray(z["theta_rec"], float)
    nd = len(th) - NP_SRC * int(z["n_slot"])   # D3 cells have nd = 116 + n_add
    out = {}
    for k in fed:
        out[int(k)] = (th[nd + NP_SRC*k + I_FGW], th[nd + NP_SRC*k + I_MC])
    return out

rows = []
for key in sorted(cells):
    its = cells[key]
    iters = sorted(its)
    Z = {i: np.load(its[i], allow_pickle=True) for i in iters}
    nres = np.array([int(Z[i]["n_resolved"]) for i in iters])
    ncert = np.array([int(Z[i]["n_cert"]) for i in iters])
    wrong = np.array([int(Z[i]["wrong_cert"]) for i in iters])
    abg = np.array([float(Z[i]["a_bg"]) for i in iters])
    aeff = float(Z[iters[0]]["a_eff_drawn"])
    floor0 = float(Z[iters[0]]["floor"]); zf0 = float(Z[iters[0]]["zero_fraction"])
    tension = float(Z[iters[0]]["cond_tension_sigma"]) if "cond_tension_sigma" in Z[iters[0]] else np.nan
    hbr = float(Z[iters[0]]["cond_h_brightest"]) if "cond_h_brightest" in Z[iters[0]] else np.nan

    # per-interval max wander over members fed at BOTH ends
    fedth = {i: fed_theta(Z[i]) for i in iters}
    wmax, nfed_common = [], []
    for a, b in zip(iters[:-1], iters[1:]):
        common = set(fedth[a]) & set(fedth[b])
        nfed_common.append(len(common))
        if not common:
            wmax.append(np.nan); continue
        wmax.append(max(max(abs(fedth[b][k][0]-fedth[a][k][0]),
                            abs(fedth[b][k][1]-fedth[a][k][1])) for k in common))
    wmax = np.array(wmax, float)
    stable = np.array([(not np.isnan(x)) and x <= HOLD_DEX and n > 0
                       for x, n in zip(wmax, nfed_common)], bool)
    # >=2 CONSECUTIVE stable intervals
    hold = any(stable[j] and stable[j+1] for j in range(len(stable)-1)) if len(stable) > 1 else False

    rung, kind, arm, sky, w, T, tier = key
    rows.append(dict(rung=rung, kind=kind, arm=arm, sky=sky, w=w, T=T, tier=tier,
                     n_iters=len(iters),
                     first_bite=bool((nres > 0).any()), cert=bool((ncert > 0).any()),
                     hold=bool(hold),
                     nres_max=int(nres.max()), nres_i0=int(nres[0]),
                     ncert_max=int(ncert.max()), wrong_tot=int(wrong.sum()),
                     bite_i0=abg[0]-aeff, regress_end=abg[-1]-aeff,
                     abg_first=abg[0], abg_last=abg[-1], a_eff=aeff,
                     wander_max=float(np.nanmax(wmax)) if len(wmax) and not np.all(np.isnan(wmax)) else np.nan,
                     wander_min=float(np.nanmin(wmax)) if len(wmax) and not np.all(np.isnan(wmax)) else np.nan,
                     floor0=floor0, zf0=zf0, tension=tension, h_brightest=hbr))

hdr = ("rung  kind   arm  sky w     T  it FB CT HOLD nresM nc wr  bite_i0 regr_end "
       "wandMax wandMin  floor0  zf   tens")
print(hdr)
for r in rows:
    print(f"{r['rung']:5s} {r['kind']:6s} {r['arm']:4s} {r['sky']:2d} {r['w']:5s} "
          f"{r['T']:2d} {r['tier'][:5]:5s} {r['n_iters']:2d} "
          f"{int(r['first_bite']):2d} {int(r['cert']):2d} "
          f"{int(r['hold']):4d} {r['nres_max']:5d} {r['ncert_max']:2d} {r['wrong_tot']:2d} "
          f"{r['bite_i0']:+8.3f} {r['regress_end']:+8.3f} "
          f"{r['wander_max']:7.3f} {r['wander_min']:7.3f} {r['floor0']:7.2f} "
          f"{r['zf0']:4.2f} {r['tension']:6.1f}")

# lean bank
keys = list(rows[0].keys())
np.savez("reports/summit_s0_harvest.npz",
         **{k: np.array([r[k] for r in rows]) for k in keys},
         note="SUMMIT Phase-0 harvest of gl2 ladder banks; HOLD = >=2 consecutive "
              "intervals, all common-fed members <=0.05 dex in (fgw,mc); "
              "bite/regress = a_bg - a_eff_drawn at i0 / last iter.")
print(f"\nbanked reports/summit_s0_harvest.npz  ({len(rows)} cells)")
