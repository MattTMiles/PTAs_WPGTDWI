# JOB 53 Part 0.4: natural against constructed skies (census priors, skies 100-119)
import os, json
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
from job53_lib import *
OUT = {}
R = {"constructed": (ld("main_con"), ld("stat_con")), "natural": (ld("main_nat20"), ld("stat_nat20"))}
lk = ld("mainlk_nat20")
print("== reproduction of Job 52 (constructed, census skies 100-119): all / K>=60 / K>=600")
for key, want in (("w|ref|12|t1", "0.053 / 0.021 / 0.000"), ("w|ref|all|t1", "0.271 / 0.188 / 0.080"), ("w|A|12|t1", "0.098 / 0.054 / 0.011")):
    row, _ = table_row(R["constructed"][0]["rec"], key, "constructed", ci=False); print(f"   {key}: {fmt3(row)}   (Job 52: {want})")
    assert fmt3(row) == want
for st, (M, S) in R.items():
    rec, info = M["rec"], S["info"]; nsk = len(info)
    print(f"\n== {st} skies: {nsk} skies, {len(rec)} pulsars, {len(np.unique(popid(rec, st)))} distinct populations")
    cl = classes(rec); print(f"   pulsars per class: all {len(rec)}, K>=60 {int(cl['K60'].sum())}, K>=600 {int(cl['K600'].sum())}")
    nb = np.array([(d["arr_bg"] >= 3).sum() for d in info]); no = np.array([(d["arr_opt"] >= 3).sum() for d in info]); nc = np.array([d["ncand"] for d in info])
    print(f"   sources per sky detectable against the background: arr_bg >= 3: mean {nb.mean():.1f} (min {nb.min()}, median {np.median(nb):.0f}, max {nb.max()}); arr_bg_opt >= 3: mean {no.mean():.1f} (min {no.min()}, median {np.median(no):.0f}, max {no.max()}); candidates {nc.mean():.0f}; sources {np.mean([d['nsrc'] for d in info]):.0f}")
    if st == "natural":
        nr = np.array([d["nrej"] for d in M["info"]]); print(f"   sources rejected by the CW upper limit in the draw, per sky: {nr.tolist()} (mean {nr.mean():.2f})")
    OUT[st] = dict(arr_bg=float(nb.mean()), arr_opt=float(no.mean()), ncand=float(nc.mean()))
    print("   share of the 2-30 nHz power (Earth + pulsar terms) held by the top N sources: " + ", ".join(f"N={n}: {np.mean([r[f'sh|{n}'] for r in rec]):.3f}" for n in (1, 3, 7, 12, 30, "all")))
    OUT[st]["share"] = {str(n): float(np.mean([r[f"sh|{n}"] for r in rec])) for n in (1, 3, 7, 12, 30, "all")}
    bm = np.array([r["d_bbin_max"] for r in rec]); bu = np.concatenate([r["d_bun_on"] for r in rec]); b12 = np.concatenate([r["d_bbin12"] for r in rec])
    print(f"   largest beta_bin per pulsar (over all candidates) 10/50/90%: {np.percentile(bm, [10, 50, 90]).round(3)}; max {bm.max():.3f}")
    print(f"   beta_unmod, N = 12 pairs ({len(bu)} pairs) 10/50/90%: {np.percentile(bu, [10, 50, 90]).round(4)}; 99%: {np.percentile(bu, 99):.3f}; max {bu.max():.3f}; fraction above 0.3: {np.mean(bu > 0.3):.4f}")
    print(f"   beta_bin, N = 12 pairs 10/50/90%: {np.percentile(b12, [10, 50, 90]).round(4)}")
    OUT[st]["bun"] = [float(x) for x in np.percentile(bu, [10, 50, 90])]; OUT[st]["bbin_max"] = [float(x) for x in np.percentile(bm, [10, 50, 90])]
    print(f"   prior only: {fmt3(prior_only(rec, st))}")
    print("   decode (all / K>=60 / K>=600) [95% interval, populations resampled] | right per sky")
    for v in ("ref", "A"):
        for n in (7, 12, "all"):
            row, rps = table_row(rec, f"w|{v}|{n}|t1", st); print(f"      {v:>3s} N={str(n):>3s}: {fmt3(row)} | {rps:.1f}")
            OUT[st][f"{v}|{n}"] = row
print("\n== natural skies 100-119, chirp mass from the amplitude lookup (as Jobs 50-52) instead of the drawn entry")
for n in (7, 12, "all"):
    a, _ = table_row(lk["rec"], f"w|ref|{n}|t1", "natural", ci=False); b, _ = table_row(R["natural"][0]["rec"], f"w|ref|{n}|t1", "natural", ci=False)
    print(f"   N={n}: lookup {fmt3(a)}   exact {fmt3(b)}")
json.dump(OUT, open("job53_data/skies.json", "w"), indent=1)
