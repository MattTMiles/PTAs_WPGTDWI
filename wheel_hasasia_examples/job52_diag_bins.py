# JOB 52 diagnostic: why do arms B and D get WORSE when timing noise is removed?
# Test: same extended configuration, zero timing noise, but bins above 33 nHz ignored.
import os, sys, pickle
os.environ.setdefault("MPLBACKEND", "Agg"); os.environ.pop("DISPLAY", None)
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job52_common as H, job52_run as R
os.chdir(C.J.BASE); C.J.load(verbose=False)
SK = [100, 101, 102, 103]
_, kw, hk, rk = R.RUNS["mainx_h"]
src = pickle.load(open("job52_data/parts/statx_h_100.pkl", "rb"))[1]
rank = {d["sky"]: (d["cand"], d["arr_opt"]) for d in src}
arms = [(v, n, "tinf") for v in ("B", "D") for n in (12, "all")]
old = [r for r in pickle.load(open("job52_data/run_mainx_h.pkl", "rb"))["rec"] if r["sky"] in SK]
rec, _ = C.collect50(SK, hook=H.Joint52(rank, ext=True, binmax=3.3e-8, arms=arms), **kw)
print(f"DIAG harness skies 100-103, {len(rec)} pulsars, zero timing noise, extended pool")
for a in arms:
    k = f"w|{a[0]}|{a[1]}|tinf"; k1 = f"w|{a[0]}|{a[1]}|t1"
    print(f"DIAG  {a[0]} N={a[1]}: bins to 100 nHz right {1-np.mean([r[k] for r in old]):.3f} | bins above 33 nHz ignored {1-np.mean([r[k] for r in rec]):.3f} | (x1 timing, bins to 100 nHz: {1-np.mean([r[k1] for r in old]):.3f})")
print(f"DIAG  |mean normalised product of neighbouring bins| of the unmodelled data: bins > 40 nHz {np.nanmean([r['corr_hi'] for r in rec]):.3f}; bins 2-10 nHz {np.mean([r['corr_lo'] for r in rec]):.3f} (independent bins: ~{1/np.sqrt(80):.2f} and ~{1/np.sqrt(11):.2f})")
