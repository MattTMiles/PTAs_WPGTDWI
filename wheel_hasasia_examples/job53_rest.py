# ================= JOB 53 -- analysis, Part 4 (observing span) =================================
import os, json
import numpy as np
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
from job53_lib import *
OUT = {}
print("== Part 4: observing span, natural skies, census priors, reference")
for T, mn, sn in ((15, "span15_nat", "stat15_nat"), (25, "span25_nat", "stat25_nat"), (45, "main_nat", "stat_nat")):
    M = ld(mn); S = ld(sn); rec = M["rec"]; info = S["info"]
    nsk = len(np.unique([r["sky"] for r in rec])); Ts = T*365.25*86400
    nb = np.array([(d["arr_bg"] >= 3).sum() for d in info]); no = np.array([(d["arr_opt"] >= 3).sum() for d in info]); nc = np.array([d["ncand"] for d in info])
    bu = np.concatenate([r["d_bun_on"] for r in rec]); mv = np.concatenate([r["d_mov"] for r in rec])
    cl = classes(rec)
    print(f"\n-- T = {T} yr: {nsk} skies of 200 requested, {len(rec)} pulsars ({len(rec)/nsk:.0f} per sky); bins to 33 nHz: {int(3.3e-8*Ts)+1}; candidates per sky {nc.mean():.0f} (min {nc.min()}); booked list {np.mean([d['n12'] for d in info]):.1f}")
    print(f"   detectable per sky: arr_bg>=3 {nb.mean():.2f}, arr_bg_opt>=3 {no.mean():.2f}; beta_unmod N=12 pairs 10/50/90% {np.percentile(bu, [10, 50, 90]).round(4)}; bins moved median {np.median(mv):.2f}; share of band power top 12: {np.mean([r['sh|12'] for r in rec]):.3f}")
    print(f"   classes: K>=60 {int(cl['K60'].sum())}, K>=600 {int(cl['K600'].sum())}; prior only {fmt3(prior_only(rec, 'natural'))}; median Kc {np.median([r['Kc'] for r in rec]):.0f}")
    OUT[T] = dict(nsk=nsk, nb=int(3.3e-8*Ts)+1, arr_bg=float(nb.mean()), arr_opt=float(no.mean()), bun=[float(x) for x in np.percentile(bu, [10, 50, 90])], ncand=float(nc.mean()))
    for n in (7, 12, "all"):
        row, rps = table_row(rec, f"w|ref|{n}|t1", "natural"); print(f"   reference N={str(n):>3s}: {fmt3(row)} | {rps:.2f} per sky"); OUT[T][str(n)] = row
        # same known-set tooth classes: K from the known set's coarsest tooth
    Kk = np.array([r["Kco|12"] for r in rec]); ok = ~np.array([r["w|ref|12|t1"] for r in rec])
    print(f"   N=12 by K of the known set's coarsest tooth: K<10 {ok[Kk < 10].mean() if (Kk < 10).any() else float('nan'):.3f} (n={int((Kk < 10).sum())}); 10-100 {ok[(Kk >= 10) & (Kk < 100)].mean():.3f} (n={int(((Kk >= 10) & (Kk < 100)).sum())}); >=100 {ok[Kk >= 100].mean():.3f} (n={int((Kk >= 100).sum())})")
json.dump(OUT, open("job53_data/rest.json", "w"), indent=1)
