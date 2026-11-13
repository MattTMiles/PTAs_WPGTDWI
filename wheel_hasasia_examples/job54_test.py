# JOB 54 test: Joint54 against Job 53's saved flags (natural sky 100), and timing of one sky
import os, sys, pickle, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import numpy as np
import job54_run as R
C = R.C; C.J.load(verbose=False); j = 100; t0 = time.time()
kw, hk, rk = R.RUNS["main"]; src = pickle.load(open(f"{rk}_{j:03d}.pkl", "rb"))[1]
arms = hk[1]+[("A", 12, "t1"), ("kn", 12, "t1"), ("ref", "all", "tinf")]
rec, info = C.collect50([j], hook=R.H.Joint54({d["sky"]: (d["cand"], d["arr_opt"]) for d in src}, arms, diag=True), **kw)
old = pickle.load(open(f"job53_data/parts/main_nat_{j:03d}.pkl", "rb"))[0]
assert len(rec) == len(old) == 79
for a in arms:
    k = f"w|{a[0]}|{a[1]}|{a[2]}"
    if k not in old[0]: continue
    same = sum(r[k] == o[k] for r, o in zip(rec, old)); dm = max(abs(r["m"+k[1:]]-o["m"+k[1:]]) for r, o in zip(rec, old))
    print(f"natural sky {j}, arm {a}: flags equal Job 53 for {same}/79 pulsars; max |margin diff| {dm:.1e}"); assert same == 79
db = max(float(np.max(np.abs(r["d_bun_on"]-o["d_bun_on"]))) for r, o in zip(rec, old)); print(f"beta_unmod (on) max diff vs Job 53: {db:.1e}")
print("new arms, fraction right (79 pulsars): " + ", ".join(f"{a[0]}|{a[1]}: {np.mean([not r[f'w|{a[0]}|{a[1]}|{a[2]}'] for r in rec]):.3f}" for a in arms if a[0] in ("q1", "q2", "q2s", "fo", "ao")))
print(f"PASS; one sky {time.time()-t0:.0f} s")
