# ================= JOB 52 -- production runs (fork after load; one file per chunk) ============
import os, sys, time, pickle
os.environ.setdefault("MPLBACKEND", "Agg"); os.environ.pop("DISPLAY", None)
import multiprocessing as mp
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job52_common as H
os.chdir(C.J.BASE); os.makedirs("job52_data/parts", exist_ok=True)
ORA = [("on", "all", "oracle"), ("off", "all", "oracle")]
FUL = [("on", "all", "full"), ("off", "all", "full")]
BASE = dict(marg=False, rate=1.0, sens_n=1.0, fine=False)
EXT = dict(fmax_cand=1.0e-7, cand_cap=5000, fgrid_max=1.05e-7)
SK = {"h": range(100, 200), "c": range(100, 120)}
PRI = {"h": {}, "c": {"priors": "census"}}
RUNS = {}
for p in "hc":
    RUNS[f"stat_{p}"] = (SK[p], dict(BASE, arms=FUL, **PRI[p]), ("stat",), None)
    RUNS[f"statx_{p}"] = (SK[p], dict(BASE, arms=FUL, **PRI[p], **EXT), ("stat",), None)
    RUNS[f"main_{p}"] = (SK[p], dict(BASE, arms=ORA, **PRI[p]), ("joint", False, "main", True, False), f"stat_{p}")
    RUNS[f"mainx_{p}"] = (SK[p], dict(BASE, arms=ORA, **PRI[p], **EXT), ("joint", True, "main", True, False), f"statx_{p}")
    RUNS[f"p5_{p}"] = (SK[p], dict(BASE, arms=ORA, **PRI[p]), ("joint", False, "p5", True, True), f"stat_{p}")
    RUNS[f"p5off_{p}"] = (SK[p], dict(BASE, arms=ORA, **PRI[p]), ("joint", False, "p5", False, True), f"stat_{p}")
CHUNK = {"h": 4, "c": 1}

def work(task):
    name, skies = task
    _, kw, hk, rk = RUNS[name]; kw = dict(kw)
    fn = f"job52_data/parts/{name}_{skies[0]:03d}.pkl"
    if not os.path.exists(fn):
        t0 = time.time()
        if hk[0] == "stat": kw["hook"] = H.Stat52(store=False)
        else:
            src = pickle.load(open(f"job52_data/parts/{rk}_{skies[0]:03d}.pkl", "rb"))[1]
            rank = {d["sky"]: (d["cand"], d["arr_opt"]) for d in src}
            kw["hook"] = H.Joint52(rank, ext=hk[1], mode=hk[2], ret=hk[3], part4=hk[4])
        res = C.collect50(list(skies), **kw)
        pickle.dump(res, open(fn + ".tmp", "wb")); os.replace(fn + ".tmp", fn)
        print(f"  [{name} skies {skies[0]}-{skies[-1]}] {time.time()-t0:.0f} s", flush=True)
    return name

if __name__ == "__main__":
    nw = int(sys.argv[1]); names = sys.argv[2:]
    t0 = time.time(); C.J.load(verbose=False); print(f"[load] {time.time()-t0:.0f} s", flush=True)
    tasks = []
    for n in names:
        sk = list(RUNS[n][0]); ch = CHUNK[n[-1]]
        tasks += [(n, sk[a:a+ch]) for a in range(0, len(sk), ch)]
    with mp.get_context("fork").Pool(nw) as pool:
        for k, name in enumerate(pool.imap_unordered(work, tasks)):
            if k % 10 == 0: print(f"[progress] {k+1}/{len(tasks)} chunks, {time.time()-t0:.0f} s", flush=True)
    for n in names:
        sk = list(RUNS[n][0]); ch = CHUNK[n[-1]]; rec = []; info = []
        for a in range(0, len(sk), ch):
            r_, i_ = pickle.load(open(f"job52_data/parts/{n}_{sk[a]:03d}.pkl", "rb")); rec += r_; info += i_
        rec.sort(key=lambda r: (r["sky"], r["psr"])); info.sort(key=lambda d: d["sky"])
        pickle.dump(dict(rec=rec, info=info), open(f"job52_data/run_{n}.pkl", "wb"))
        print(f"[{n}] records {len(rec)} skies {len(info)}", flush=True)
    print(f"TOTAL {time.time()-t0:.0f} s")
