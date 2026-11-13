# JOB 52 -- rule 2 on one sky per prior set: new one-FFT joint fit vs Job 51's saved flags; timing
import os, sys, time, pickle
os.environ.setdefault("MPLBACKEND", "Agg"); os.environ.pop("DISPLAY", None)
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job52_common as H, job52_run as R
os.chdir(C.J.BASE); C.J.load(verbose=False)
for p, j51 in (("h", "joint_h_s1"), ("c", "joint_c_s1")):
    old = [r for r in pickle.load(open(f"job51_data/run_{j51}.pkl", "rb"))["rec"] if r["sky"] == 100]
    for nm in (f"main_{p}", f"mainx_{p}"):
        _, kw, hk, rk = R.RUNS[nm]
        src = pickle.load(open(f"job52_data/parts/{rk}_100.pkl", "rb"))[1]
        rank = {d["sky"]: (d["cand"], d["arr_opt"]) for d in src}
        t0 = time.time(); rec, info = C.collect50([100], hook=H.Joint52(rank, ext=hk[1]), **kw); dt = time.time()-t0
        if not hk[1]:
            a = np.array([r["w|ref|all|t1"] for r in rec]); b = np.array([r["wj"] for r in old])
            print(f"TEST {nm}: ref|all|t1 flags equal to Job 51 wj: {int((a==b).sum())}/{len(a)}; right {1-a.mean():.3f} (Job 51 {1-b.mean():.3f})")
            for n in (1, 7, 12, 30, 300):
                print(f"TEST   N={n}: right {1-np.mean([r[f'w|ref|{n}|t1'] for r in rec]):.3f}  share {np.median([r[f'sh|{n}'] for r in rec]):.3f}")
            for v in ("A", "C"):
                print(f"TEST   {v} all t1: right {1-np.mean([r[f'w|{v}|all|t1'] for r in rec]):.3f}; N=12: {1-np.mean([r[f'w|{v}|12|t1'] for r in rec]):.3f}")
            print(f"TEST   C: pairs moving >0.1 bin (all cand) median {np.median([r['c_mov01_all'] for r in rec]):.3f}; nseg max {max(r['c_nseg'] for r in rec)}; segcap {sum(r['c_segcap'] for r in rec)}")
        else:
            for v in ("B", "D"):
                print(f"TEST   {v} all t1: right {1-np.mean([r[f'w|{v}|all|t1'] for r in rec]):.3f}; N=12: {1-np.mean([r[f'w|{v}|12|t1'] for r in rec]):.3f}; pool {np.median([r['nC'] for r in rec]):.0f}")
        print(f"TEST {nm}: {dt:.0f} s for 1 sky, {len(rec)} pulsars; M max {max(r['j_M'] for r in rec)}", flush=True)
