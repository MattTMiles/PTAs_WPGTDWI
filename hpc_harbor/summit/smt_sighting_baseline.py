#!/usr/bin/env python
"""PHOENIX -- SUMMIT S0.2a: THE SIGHTING CALL, against the same-rung null baseline.

SUMMIT S0.2a posted two STRONG candidate closure events at r13p25 in the circular arm
and named the decisive test:

  "The completed baseline test is: do the rerunning null cells (12766445) late-feed
   their real members at the same rate/shape? A matching rate kills the attribution;
   a null arm that only ever i0-feeds (or never feeds) leaves the two candidates
   standing."

12766445 has since landed: the r13p25 none-arm nulls are banked 6/6 across
null0/null1/null2 x s0/s1. This script runs that comparison. CPU only.

SHAPE TAXONOMY (as SUMMIT measured it, applied identically to both arms):
  GRADUAL   -- >=2 consecutive strictly-decreasing refits before the crossing
               (the monotone multi-refit descent; the only class SUMMIT left standing
               as candidate closure evidence after STOP #S2)
  SINGLE    -- flat at/near 1.0, then one step across the bar (frontier jitter; the
               null demonstrated this onto a WRONG-SKY template at realistic SNR)
  HOVER     -- multiple iterations within +-0.05 of the bar before crossing

LANE CAVEAT, stated up front: the r13p25 SIGNAL cells were cut on dgx03/A100-80GB, the
r13p25 NULL reruns on hgx03/H200. conc_ratio is Fisher-derived and its measured
cross-lane wobble is 1.66e-02 (SG-D2a). That is 14x below the feed-bar margin and
1-2 orders below the descents being classified, so the SHAPE taxonomy is readable
cross-lane; a rate comparison at the bar would not be, and is not made here.
"""
import glob, re, sys
import numpy as np

RESULTS = "GLACIER_results"
BAR = 0.5
pat = re.compile(r"(?P<stem>gl2_r13p25_(?P<kind>cell|null\d+)_none_s(?P<sky>\d+)_T30_lit)"
                 r"_i(?P<it>\d+)__(?P<lane>.+)\.npz")

cells = {}
for f in sorted(glob.glob(f"{RESULTS}/gl2_r13p25_*_none_s*_T30_lit_i*.npz")):
    m = pat.search(f)
    if not m or "STOPANAT" in f:
        continue
    cells.setdefault((m.group("kind"), int(m.group("sky")), m.group("lane")), {})[int(m.group("it"))] = f


def classify(traj, it_fed):
    """traj = conc_ratio for one slot over iterations 0..it_fed (inclusive)."""
    pre = traj[:it_fed + 1]
    if len(pre) < 2:
        return "SINGLE", 0.0
    steps = np.diff(pre)
    # longest run of strictly-decreasing refits ending at the crossing
    run = 0
    for s in steps[::-1]:
        if s < -1e-9:
            run += 1
        else:
            break
    drop = float(pre[0] - pre[-1])
    near = np.sum(np.abs(pre[:-1] - BAR) <= 0.05)
    if run >= 2:
        return "GRADUAL", drop
    if near >= 2:
        return "HOVER", drop
    return "SINGLE", drop


rows = []
for key in sorted(cells):
    kind, sky, lane = key
    its = sorted(cells[key])
    Z = {i: np.load(cells[key][i], allow_pickle=True) for i in its}
    n_src = int(Z[its[0]]["n_src"])
    R = np.vstack([np.asarray(Z[i]["conc_ratio"], float) for i in its])   # (n_it, n_slot)
    fed = np.vstack([np.asarray(Z[i]["fed_mask"], bool) for i in its])
    # first iteration at which each slot is fed
    first_fed = {}
    for s in range(fed.shape[1]):
        w = np.flatnonzero(fed[:, s])
        if w.size:
            first_fed[s] = int(its[w[0]])
    i0_feeds = [s for s, i in first_fed.items() if i == its[0]]
    late = {s: i for s, i in first_fed.items() if i > its[0]}
    late_census = {s: i for s, i in late.items() if s < n_src}
    for s, i in sorted(late_census.items()):
        cls, drop = classify(R[:, s], its.index(i))
        rows.append(dict(kind=kind, sky=sky, lane=lane, slot=s, it=i, cls=cls, drop=drop,
                         traj=R[:its.index(i) + 1, s]))
    print(f"{kind:7s} s{sky} [{lane[:5]}] n_src={n_src:3d}  i0-feeds={len(i0_feeds):3d}  "
          f"late-feeds={len(late):2d}  late CENSUS feeds={len(late_census)}")

print("\n" + "=" * 96)
print("LATE CENSUS FEEDS -- shape taxonomy, both arms, same rung (r13p25, circular)")
print("=" * 96)
print(f"{'arm':8s}{'cell':10s}{'slot':>5s}{'fed@i':>6s}{'class':>9s}{'drop':>7s}   trajectory to feed")
for r in sorted(rows, key=lambda r: (r["kind"] != "cell", r["kind"], r["sky"])):
    arm = "SIGNAL" if r["kind"] == "cell" else "NULL"
    t = " -> ".join(f"{v:.3f}" for v in r["traj"])
    print(f"{arm:8s}{r['kind']+'/s'+str(r['sky']):10s}{r['slot']:>5d}{r['it']:>6d}"
          f"{r['cls']:>9s}{r['drop']:>7.3f}   {t}")

sig = [r for r in rows if r["kind"] == "cell"]
nul = [r for r in rows if r["kind"] != "cell"]
n_sig_cells = len({(r["kind"], r["sky"]) for r in rows if r["kind"] == "cell"})
n_nul_cells = len([k for k in cells if k[0] != "cell"])
print("\n" + "-" * 96)
print(f"SIGNAL arm: {len(sig)} late census feeds; classes "
      f"{ {c: sum(1 for r in sig if r['cls'] == c) for c in ('GRADUAL','SINGLE','HOVER')} }")
print(f"NULL   arm: {len(nul)} late census feeds over {n_nul_cells} cells; classes "
      f"{ {c: sum(1 for r in nul if r['cls'] == c) for c in ('GRADUAL','SINGLE','HOVER')} }")
ng = sum(1 for r in nul if r["cls"] == "GRADUAL")
sg = sum(1 for r in sig if r["cls"] == "GRADUAL")
print("-" * 96)
print(f"\nDECISIVE COMPARISON (SUMMIT S0.2a's stated test):")
print(f"  GRADUAL (monotone multi-refit descent) in SIGNAL arm : {sg}")
print(f"  GRADUAL (monotone multi-refit descent) in NULL   arm : {ng}")
if ng == 0 and sg > 0:
    print("  -> the null arm produces NO analog of the gradual class at the same rung.")
    print("     The attribution is NOT killed by the control. Candidates STAND.")
elif ng >= sg:
    print("  -> the control reproduces the gradual class at >= the signal rate.")
    print("     ATTRIBUTION KILLED -- the shape is not evidence of real structure.")
else:
    print("  -> the control produces the gradual class at a LOWER but nonzero rate.")
    print("     Attribution WEAKENED, not killed; quote the rate ratio, claim nothing.")
print("\nNOTE: this is the BASELINE result. Declaring a SIGHTING is charter STOP #5 -- Matt's call.")
