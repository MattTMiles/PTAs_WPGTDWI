#!/usr/bin/env python3
"""
VERNIER-28 -- the PAIRED comparison, recomputed from the banked run log.

The three ladders see bit-identical skies, arrays, noise draws and data (the
V-series ablation rule), so the honest comparison is paired: the difference in
conversion efficiency universe by universe, not the difference of two cell
means with independent error bars. This script reads the per-universe lines
`vernier28.py` prints and bootstraps those differences.

    python3 vernier28_paired.py [reports/vernier28.log]
"""
import re
import sys
import collections
import numpy as np

PAT = re.compile(
    r"\[\s*(\w+) r(\d+)\] K_av\s+(\d+)\s+gen R1/R2/R3\s+(\d+)/\s*(\d+)/\s*(\d+)"
    r"\s+phantom\s+(\d+)/\s*(\d+)/\s*(\d+)\s+S_end ([\d.]+)/([\d.]+)/([\d.]+)")
CELLS = ("C5", "C10", "C20", "C10b")


def boot(x, n=20000, seed=7):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    bs = np.array([x[rng.integers(0, x.size, x.size)].mean() for _ in range(n)])
    return x.mean(), float(bs.std(ddof=1)), np.percentile(bs, [2.5, 97.5])


def main(path):
    d = collections.defaultdict(list)
    for line in open(path):
        m = PAT.search(line)
        if m:
            g = m.groups()
            d[g[0]].append([float(g[2])] + [float(v) for v in g[3:]])
    n = sum(len(v) for v in d.values())
    print(f"  parsed {n} universes from {path}")
    print()
    A = np.array([v for c in CELLS for v in d[c]], float)
    kav, G, P = A[:, 0], A[:, 1:4], A[:, 4:7]
    eff = G / kav[:, None]
    print("  PAIRED over all universes -- same sky, same array, same noise draws,")
    print("  same Earth-block data; the ladders differ only in the fit.")
    print()
    for j, lab in enumerate(("R1", "R2", "R3")):
        m, se, ci = boot(eff[:, j])
        print(f"    conversion {lab} = {m:.3f} +- {se:.3f}  [{ci[0]:.3f}, {ci[1]:.3f}]")
    print()
    for j, lab in ((1, "R2 - R1"), (2, "R3 - R1")):
        m, se, ci = boot(eff[:, j] - eff[:, 0])
        verdict = ("consistent with zero" if ci[0] < 0 < ci[1] else "EXCLUDES zero")
        print(f"    DELTA {lab} = {m:+.3f} +- {se:.3f}  "
              f"[{ci[0]:+.3f}, {ci[1]:+.3f}]   {verdict}")
    m, se, ci = boot(eff[:, 1] - eff[:, 2])
    print(f"    DELTA R2 - R3 = {m:+.3f} +- {se:.3f}  [{ci[0]:+.3f}, {ci[1]:+.3f}]")
    print()
    print("  per cell:")
    for c in CELLS:
        a = np.array(d[c], float)
        e = a[:, 1:4] / a[:, 0:1]
        m21, s21, _ = boot(e[:, 1] - e[:, 0])
        m31, s31, _ = boot(e[:, 2] - e[:, 0])
        print(f"    {c:>5} n={len(a):2d}   R2-R1 {m21:+.3f}+-{s21:.3f}   "
              f"R3-R1 {m31:+.3f}+-{s31:.3f}")
    print()
    print(f"  phantoms per run: R1 {P[:,0].mean():.2f}  R2 {P[:,1].mean():.2f}  "
          f"R3 {P[:,2].mean():.2f}")
    print(f"  head-to-head on genuine: R2 > R1 in {(G[:,1]>G[:,0]).sum()}/{len(G)}, "
          f"tied {(G[:,1]==G[:,0]).sum()}/{len(G)}, "
          f"R2 < R1 {(G[:,1]<G[:,0]).sum()}/{len(G)}; "
          f"R3 > R1 in {(G[:,2]>G[:,0]).sum()}/{len(G)}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else "/home/mattm/projects/HSYMT/reports/vernier28.log")
