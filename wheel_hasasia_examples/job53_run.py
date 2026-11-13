# ================= JOB 53 -- production runs (fork after load; one file per sky) ===================
import os, sys, time, pickle
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
os.environ["MPLBACKEND"] = "Agg"; os.environ.pop("DISPLAY", None)
import multiprocessing as mp
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C, job52_common as H52, job53_common as H
os.chdir(C.J.BASE); os.makedirs("job53_data/parts", exist_ok=True)
ORA = [("on", "all", "oracle"), ("off", "all", "oracle")]
FUL = [("on", "all", "full"), ("off", "all", "full")]
BASE = dict(marg=False, rate=1.0, sens_n=1.0, fine=False, priors="census")
NNAT = int(os.environ.get("JOB53_NNAT", "100"))
SK = {"con": range(100, 120), "nat": range(100, 100+NNAT), "nat20": range(100, 120)}
RUNS = {}       # name: (skies, collector kwargs, hook spec, ranking run)
LK = [("ref", n, "t1") for n in H.N3]
RUNS["p1x_nat20"] = (SK["nat20"], dict(BASE, arms=ORA, sky="natural"), ("joint", [("kn", 12, "t1"), ("ref", 12, "t1"), ("off", 12, "t1")], (), (), True), "stat_nat")
RUNS["statlk_nat20"] = (SK["nat20"], dict(BASE, arms=FUL, sky="natural", mass="lookup"), ("stat",), None)
RUNS["mainlk_nat20"] = (SK["nat20"], dict(BASE, arms=ORA, sky="natural", mass="lookup"), ("joint", LK, (), (), False), "statlk_nat20")
for k, sky in (("con", "constructed"), ("nat", "natural")):
    RUNS[f"stat_{k}"] = (SK[k], dict(BASE, arms=FUL, sky=sky), ("stat",), None)
    RUNS[f"main_{k}"] = (SK[k], dict(BASE, arms=ORA, sky=sky), ("joint", H.arms_main(), (), (), True), f"stat_{k}")
    RUNS[f"sig_{k}"] = (SK[k], dict(BASE, arms=ORA, sky=sky), ("joint", [], H.arms_sig(), H.SIGS, False), f"stat_{k}")
    for T in (15, 25):
        RUNS[f"stat{T}_{k}"] = (SK[k], dict(BASE, arms=FUL, sky=sky, T=T), ("stat",), None)
        RUNS[f"span{T}_{k}"] = (SK[k], dict(BASE, arms=ORA, sky=sky, T=T), ("joint", H.arms_span(), (), (), True), f"stat{T}_{k}")


def work(task):
    name, j = task
    _, kw, hk, rk = RUNS[name]; kw = dict(kw)
    fn = f"job53_data/parts/{name}_{j:03d}.pkl"
    if not os.path.exists(fn):
        t0 = time.time()
        if hk[0] == "stat": kw["hook"] = H52.Stat52(store=False)
        else:
            src = pickle.load(open(f"job53_data/parts/{rk}_{j:03d}.pkl", "rb"))[1]
            rank = {d["sky"]: (d["cand"], d["arr_opt"]) for d in src}
            kw["hook"] = H.Joint53(rank, hk[1], hk[2], hk[3], diag=hk[4])
        rec, info = C.collect50([j], **kw)
        for d in info: d["nrej"] = C.NAT_REJ.get(j)
        pickle.dump((rec, info), open(fn + ".tmp", "wb")); os.replace(fn + ".tmp", fn)
        print(f"  [{name} sky {j}] {time.time()-t0:.0f} s", flush=True)
    return name


if __name__ == "__main__":
    nw = int(sys.argv[1]); names = sys.argv[2:]
    t0 = time.time(); C.J.load(verbose=False); print(f"[load] {time.time()-t0:.0f} s", flush=True)
    tasks = [(n, j) for n in names for j in RUNS[n][0]]
    with mp.get_context("fork").Pool(nw) as pool:
        for k, name in enumerate(pool.imap_unordered(work, tasks)):
            if k % 10 == 0: print(f"[progress] {k+1}/{len(tasks)} sky runs, {time.time()-t0:.0f} s", flush=True)
    for n in names:
        rec = []; info = []
        for j in RUNS[n][0]:
            r_, i_ = pickle.load(open(f"job53_data/parts/{n}_{j:03d}.pkl", "rb")); rec += r_; info += i_
        rec.sort(key=lambda r: (r["sky"], r["psr"])); info.sort(key=lambda d: d["sky"])
        pickle.dump(dict(rec=rec, info=info), open(f"job53_data/run_{n}.pkl", "wb"))
        print(f"[{n}] records {len(rec)} skies {len(info)}", flush=True)
    print(f"TOTAL {time.time()-t0:.0f} s")
