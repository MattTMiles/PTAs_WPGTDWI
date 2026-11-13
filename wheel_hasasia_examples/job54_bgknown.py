# JOB 54 side test: known set chosen from ALL in-band sources detectable against the background, not only
# those with single-pulsar SNR >= 3 (collect50(cand_bg=300)).  Natural skies 100-139, census priors.
import os, sys, pickle, json
import numpy as np
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples"); sys.path.insert(0, ".")
from job53_lib import frac, classes
new = pickle.load(open("job54_data/run_mainbg.pkl", "rb")); st = pickle.load(open("job54_data/run_statbg.pkl", "rb"))
old = [r for r in pickle.load(open("job53_data/run_main_nat.pkl", "rb"))["rec"] if 100 <= r["sky"] < 140]
so = [d for d in pickle.load(open("job53_data/run_stat_nat.pkl", "rb"))["info"] if 100 <= d["sky"] < 140]
rn = new["rec"]; OUT = {}
print(f"natural skies 100-139: {len(rn)} pulsars (enlarged pool), {len(old)} (Job 53 pool)")
for lab, S in (("Job 53 pool (single-pulsar SNR >= 3, 1,300 loudest by it)", so), ("enlarged pool (+ 300 loudest against the background)", st["info"])):
    nb = np.array([(d["arr_bg"] >= 3).sum() for d in S]); no = np.array([(d["arr_opt"] >= 3).sum() for d in S])
    f12 = np.concatenate([d["fc"][np.argsort(d["arr_opt"])[::-1][:12]] for d in S])*1e9
    print(f"   {lab}: candidates {np.mean([d['ncand'] for d in S]):.0f}; detectable per sky arr_bg >= 3: {nb.mean():.2f}; arr_bg_opt >= 3: {no.mean():.2f}; Earth frequency of the 12 top-ranked 10/50/90%: {np.percentile(f12, [10, 50, 90]).round(1)} nHz; above 15 nHz: {np.mean(f12 > 15)*12:.1f} of 12")
    OUT[lab[:8]] = dict(arr_bg=float(nb.mean()), arr_opt=float(no.mean()))
def row(rec, key):
    ok = ~np.array([r[key] for r in rec]); s = np.array([r["sky"] for r in rec]); K = np.array([r["Kc"] for r in rec])
    a = frac(ok, np.ones(len(rec), bool), s); b = frac(ok, K >= 60, s); return f"{a[0]:.3f} [{a[1]:.3f}, {a[2]:.3f}] / {b[0]:.3f} [{b[1]:.3f}, {b[2]:.3f}]"
print("   fraction right, all pulsars / K >= 60 (booked list)      Job 53 pool            | enlarged pool")
for v in ("ref", "A", "off"):
    for n in (7, 12, 30):
        k = f"w|{v}|{n}|t1"; print(f"   {v:>3s} N={n:2d}: {row(old, k) if k in old[0] else '   -   '}   |   {row(rn, k)}")
bu = np.concatenate([r["d_bun_on"] for r in rn]); bo = np.concatenate([r["d_bun_on"] for r in old]); fP = np.concatenate([r["d_fP_on"] for r in rn])*1e9; mv = np.concatenate([r["d_mov"] for r in rn])
print(f"   beta_unmod of N = 12 pairs 10/50/90%: Job 53 pool {np.percentile(bo, [10, 50, 90]).round(4)}; enlarged pool {np.percentile(bu, [10, 50, 90]).round(4)}")
kf = np.array([r["Kfin|12"] for r in rn]); kfo = np.array([r["Kfin|12"] for r in old]); kc = np.array([r["Kco|12"] for r in rn]); kco = np.array([r["Kco|12"] for r in old])
print(f"   K of the N = 12 known set, median: finest tooth {np.median(kfo):.0f} -> {np.median(kf):.0f}; coarsest tooth {np.median(kco):.1f} -> {np.median(kc):.1f}")
fE = np.concatenate([r["d_fP_off"] for r in rn])*1e9; hi = fE > 15
print(f"   pairs whose source has Earth frequency above 15 nHz: {hi.mean():.3f}; their beta_unmod 10/50/90% {np.percentile(bu[hi], [10, 50, 90]).round(4) if hi.any() else '-'}; bins moved median {np.median(mv[hi]) if hi.any() else float('nan'):.1f}; pulsar-term frequency median {np.median(fP[hi]) if hi.any() else float('nan'):.1f} nHz")
json.dump(OUT, open("job54_data/bgknown.json", "w"))
