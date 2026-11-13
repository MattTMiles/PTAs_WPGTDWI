# ================= JOB 54 -- production runs (natural skies; fork after load; one file per sky) ====
import os, sys, time, pickle
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import multiprocessing as mp
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job52_common as H52, job54_common as H
os.chdir(C.J.BASE); os.makedirs("job54_data/parts", exist_ok=True)
ORA = [("on", "all", "oracle"), ("off", "all", "oracle")]; FUL = [("on", "all", "full"), ("off", "all", "full")]
BASE = dict(marg=False, rate=1.0, sens_n=1.0, fine=False, priors="census", sky="natural")
SK = range(100, 300)
RUNS = {"main": (dict(BASE, arms=ORA), ("joint", H.arms_main(), True), "job53_data/parts/stat_nat")}
for T in (15, 25):
    RUNS[f"stat{T}"] = (dict(BASE, arms=FUL, T=T, noskip=True), ("stat",), None)
    RUNS[f"span{T}"] = (dict(BASE, arms=ORA, T=T, noskip=True), ("joint", H.arms_span(), True), f"job54_data/parts/stat{T}")


BGARMS = [(v, n, "t1") for v in ("ref", "A", "off") for n in (7, 12, 30)]
RUNS["statbg"] = (dict(BASE, arms=FUL, cand_bg=300), ("stat",), None)
RUNS["mainbg"] = (dict(BASE, arms=ORA, cand_bg=300), ("joint", BGARMS, True), "job54_data/parts/statbg")
RUNS["inf"] = (dict(BASE, arms=ORA), ("joint", [("ref", 12, "t1")], True), "job53_data/parts/stat_nat")
N4 = (7, 12, 30, "all"); SIGS = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0)
RUNS["sigr"] = (dict(BASE, arms=ORA), ("joint", [(v, n, "t1") for v in ("ref", "A") for n in N4], True, [(v, n, "t1") for v in ("ref", "A") for n in N4], SIGS, True), "job53_data/parts/stat_nat")
NLIST = (1, 2, 3, 5, 7, 12, 20, 30, 60, 100, 300, "all")
CURVE = [(v, n, "t1") for v in ("ref", "A") for n in NLIST] + [("ref", n, t) for n in N4 for t in ("t10", "tinf")]
RUNS["mainr"] = (dict(BASE, arms=ORA), ("joint", list(dict.fromkeys(H.arms_main()+CURVE)), True, (), (), True), "job53_data/parts/stat_nat")
for T in (15, 25): RUNS[f"span{T}r"] = (dict(BASE, arms=ORA, T=T, noskip=True), ("joint", H.arms_span(), True, (), (), True), f"job54_data/parts/stat{T}")
if os.environ.get("JOB54_SK"): SK = range(*[int(x) for x in os.environ["JOB54_SK"].split(",")])

def work(task):
    name, j = task; kw, hk, rk = RUNS[name]; kw = dict(kw)
    fn = f"job54_data/parts/{name}_{j:03d}.pkl"
    if not os.path.exists(fn):
        t0 = time.time()
        if hk[0] == "stat": kw["hook"] = H52.Stat52(store=False)
        else:
            src = pickle.load(open(f"{rk}_{j:03d}.pkl", "rb"))[1]
            kw["hook"] = H.Joint54({d["sky"]: (d["cand"], d["arr_opt"]) for d in src}, hk[1], diag=hk[2], sig_arms=hk[3] if len(hk) > 3 else (), sigs=hk[4] if len(hk) > 3 else (), randph=hk[5] if len(hk) > 3 else False)
        res = C.collect50([j], **kw)
        pickle.dump(res, open(fn+".tmp", "wb")); os.replace(fn+".tmp", fn)
        print(f"  [{name} sky {j}] {time.time()-t0:.0f} s", flush=True)
    return name


if __name__ == "__main__":
    nw = int(sys.argv[1]); names = sys.argv[2:]
    t0 = time.time(); C.J.load(verbose=False); print(f"[load] {time.time()-t0:.0f} s", flush=True)
    tasks = [(n, j) for n in names for j in SK]
    with mp.get_context("fork").Pool(nw) as pool:
        for k, name in enumerate(pool.imap_unordered(work, tasks)):
            if k % 20 == 0: print(f"[progress] {k+1}/{len(tasks)} sky runs, {time.time()-t0:.0f} s", flush=True)
    for n in names:
        rec = []; info = []
        for j in SK:
            r_, i_ = pickle.load(open(f"job54_data/parts/{n}_{j:03d}.pkl", "rb")); rec += r_; info += i_
        rec.sort(key=lambda r: (r["sky"], r["psr"])); info.sort(key=lambda d: d["sky"])
        pickle.dump(dict(rec=rec, info=info), open(f"job54_data/run_{n}.pkl", "wb"))
        print(f"[{n}] records {len(rec)} skies {len(info)}", flush=True)
    print(f"TOTAL {time.time()-t0:.0f} s")
