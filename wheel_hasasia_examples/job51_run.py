# ================= JOB 51 -- production runs (fork after load; one file per 4-sky chunk) ========
import os, sys, time, pickle
import multiprocessing as mp
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job51_common as H
os.chdir(C.J.BASE); os.makedirs("job51_data/parts", exist_ok=True)
SUF = "" if C.MASS == "z" else "_" + C.MASS
ORA = [("on", "all", "oracle"), ("off", "all", "oracle")]
FUL = [("on", "all", "full"), ("off", "all", "full")]
HELD = range(100, 200)
SN = {"1": 1.0, "3": 3.0, "10": 10.0, "100": 100.0, "inf": 1e6}       # "inf": timing noise / 1e6
RUNS = {}
for s in ("1", "3", "10"):      # Part 2: joint fit, fixed amplitude (rate 1), same skies as Job 50 held-out
    RUNS[f"joint_h_s{s}"] = (HELD, dict(marg=False, rate=1.0, sens_n=SN[s], arms=ORA), "joint", None)
    RUNS[f"joint_c_s{s}"] = (range(100, 120), dict(marg=False, rate=1.0, sens_n=SN[s], arms=ORA, priors="census"), "joint", None)
for s in ("1", "10", "inf"):    # Parts 3/4: statistics (no decode), Job 50 Part D setting (marginalised amplitude)
    RUNS[f"stat_h_s{s}"] = (HELD, dict(sens_n=SN[s], arms=FUL, fine=False), "stat", None)
    for v in ("bg", "opt"):     # Part 3: Part D row, full template term, dkappa from arr_bg / arr_opt
        RUNS[f"dec_{v}_s{s}"] = (HELD, dict(sens_n=SN[s], arms=FUL), None, (f"stat_h_s{s}", "arr_"+v))
RUNS["stat_c_s1"] = (HELD, dict(sens_n=1.0, arms=FUL, fine=False, priors="census"), "stat", None)

def work(task):
    name, skies = task
    _, kw, hk, ovr = RUNS[name]; kw = dict(kw)
    fn = f"job51_data/parts/{name}{SUF}_{skies[0]:03d}.pkl"
    if not os.path.exists(fn):
        t0 = time.time()
        if hk == "joint": kw["hook"] = H.JointHook()
        if hk == "stat": kw["hook"] = H.StatHook(store=name.endswith("_s1"))
        if ovr is not None:
            src = pickle.load(open(f"job51_data/parts/{ovr[0]}{SUF}_{skies[0]:03d}.pkl", "rb"))[1]
            kw["arr_ovr"] = {d["sky"]: d[ovr[1]] for d in src}
        res = C.collect50(list(skies), **kw)
        pickle.dump(res, open(fn + ".tmp", "wb")); os.replace(fn + ".tmp", fn)
        print(f"  [{name} skies {skies[0]}-{skies[-1]}] {time.time()-t0:.0f} s", flush=True)
    return name

if __name__ == "__main__":
    nw = int(sys.argv[1]); names = sys.argv[2:]
    t0 = time.time(); C.J.load(verbose=False); print(f"[load] {time.time()-t0:.0f} s", flush=True)
    tasks = []
    for n in names:
        sk = list(RUNS[n][0]); tasks += [(n, sk[a:a+4]) for a in range(0, len(sk), 4)]
    with mp.get_context("fork").Pool(nw) as pool:
        for k, name in enumerate(pool.imap_unordered(work, tasks)):
            if k % 25 == 0: print(f"[progress] {k+1}/{len(tasks)} chunks, {time.time()-t0:.0f} s", flush=True)
    for n in names:                                   # merge: one file per finished arm
        sk = list(RUNS[n][0]); rec = []; info = []
        for a in range(0, len(sk), 4):
            r_, i_ = pickle.load(open(f"job51_data/parts/{n}{SUF}_{sk[a]:03d}.pkl", "rb")); rec += r_; info += i_
        rec.sort(key=lambda r: (r["sky"], r["psr"])); info.sort(key=lambda d: d["sky"])
        pickle.dump(dict(rec=rec, info=info), open(f"job51_data/run_{n}{SUF}.pkl", "wb"))
        print(f"[{n}] records {len(rec)} skies {len(info)}", flush=True)
    print(f"TOTAL {time.time()-t0:.0f} s")
