# JOB 51 -- tests: (1) hooks leave Job 50 records bit-identical; (2) noiseless joint fit decodes
import os, sys, time, pickle
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job51_common as H
os.chdir(C.J.BASE); C.J.load(verbose=False)
A2 = [("on", "all", "oracle"), ("off", "all", "oracle")]
old = [r for r in pickle.load(open("job50_data/run_harness.pkl", "rb"))["rec"] if r["sky"] == 100]
t0 = time.time(); rec, info = C.collect50([100], arms=A2); t1 = time.time()-t0
keys = ["wf|on|all|oracle", "wb|on|all|oracle", "wf|off|all|oracle", "Ltrue", "Kc", "pri"]
same = all(a[k] == b[k] for a, b in zip(rec, old) for k in keys) and len(rec) == len(old)
print(f"(1) sky 100, marg=True, hook=None: {len(rec)} pulsars, all {len(keys)} keys identical to Job 50 run: {same}  [{t1:.0f} s]")
for nl in (True, False):
    t0 = time.time(); rec, info = C.collect50([100], marg=False, rate=1.0, arms=A2, hook=H.JointHook(noiseless=nl))
    K = np.array([r["Kc"] for r in rec]); wj = np.array([r["wj"] for r in rec]); ws = np.array([r["wsh"] for r in rec])
    wo = np.array([r["wf|on|all|oracle"] for r in rec])
    print(f"(2) noiseless={nl}: joint right {1-wj.mean():.3f} (K>=60: {1-wj[K>=60].mean():.3f}, n={int((K>=60).sum())}); "
          f"shared-data stack right {1-ws.mean():.3f} (K>=60 {1-ws[K>=60].mean():.3f}); Job50-style stack {1-wo.mean():.3f} "
          f"(K>=60 {1-wo[K>=60].mean():.3f}); cand share of GW power/bin median {np.median([r['j_share'] for r in rec]):.3f}; "
          f"N max {max(r['j_N'] for r in rec)}; cap {sum(r['j_cap'] for r in rec)}  [{time.time()-t0:.0f} s]")
t0 = time.time(); rec, info = C.collect50([100], arms=[("on", "all", "full"), ("off", "all", "full")], fine=False, hook=H.StatHook(store=True))
d = info[0]
print(f"(3) StatHook: arr_noise median {np.median(d['arr_noise']):.1f}, arr_bg {np.median(d['arr_bg']):.2f}, arr_opt {np.median(d['arr_opt']):.2f}; "
      f"n candidates with >=4 pulsars share>0.5: {(d['n_gw_any'][0]>=4).sum()} of {len(d['arr_bg'])}  [{time.time()-t0:.0f} s]")
