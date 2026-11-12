# ================= JOB 50 -- production runs (one pool, fork after load) ==============
import os, sys, time, pickle
import multiprocessing as mp
import numpy as np
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C
os.chdir(C.J.BASE); os.makedirs("job50_data", exist_ok=True)
SUF = "" if C.MASS == "z" else "_" + C.MASS      # JOB50_MASS=src -> run_*_src.pkl
ON2 = [("on", "top12", "machine"), ("on", "all", "machine")]
ALLOM = [(r, "all", t) for r in ("off", "on") for t in ("oracle", "machine")] + [("on", "top12", "machine"), ("off", "top12", "machine")]
RUNS = {
    "tau0":    (range(0, 200), dict(tau0=True, fine=False, arms=[("on", "top12", "machine"), ("off", "top12", "machine")])),
    "harness": (range(0, 200), dict()),
    "census":  (range(100, 200), dict(priors="census")),
    "sens3":   (range(100, 200), dict(sens_n=3.0, arms=ALLOM)),
    "sens10":  (range(100, 200), dict(sens_n=10.0, arms=ALLOM)),
}
for r_ in (0.33, 0.60, 1.0, 1.7, 2.7):
    RUNS[f"rate{r_:.2f}"] = (range(0, 200), dict(marg=False, rate=r_, arms=ON2 + [("off", "top12", "machine")]))

def work(task):
    name, skies, kw = task
    fn = f"job50_data/parts/{name}{SUF}_{skies[0]:03d}.pkl"
    if not os.path.exists(fn):
        t0 = time.time()
        res = C.collect50(list(skies), **kw)
        pickle.dump(res, open(fn + ".tmp", "wb")); os.replace(fn + ".tmp", fn)
        print(f"  [{name} skies {skies[0]}-{skies[-1]}] {time.time()-t0:.0f} s", flush=True)
    return name

if __name__ == "__main__":
    names = sys.argv[1:] or list(RUNS)
    os.makedirs("job50_data/parts", exist_ok=True)
    t0 = time.time(); C.J.load(verbose=False); print(f"[load] {time.time()-t0:.0f} s", flush=True)
    tasks = []
    for n in names:
        sk, kw = RUNS[n]; sk = list(sk)
        tasks += [(n, sk[a:a+4], kw) for a in range(0, len(sk), 4)]
    with mp.get_context("fork").Pool(22) as pool:
        for k, name in enumerate(pool.imap_unordered(work, tasks)):
            if k % 20 == 0: print(f"[progress] {k+1}/{len(tasks)} chunks, {time.time()-t0:.0f} s", flush=True)
    for n in names:
        sk, kw = RUNS[n]; sk = list(sk); rec = []; info = []
        for a in range(0, len(sk), 4):
            r_, i_ = pickle.load(open(f"job50_data/parts/{n}{SUF}_{sk[a]:03d}.pkl", "rb")); rec += r_; info += i_
        rec.sort(key=lambda r: (r["sky"], r["psr"])); info.sort(key=lambda d: d["sky"])
        pickle.dump(dict(rec=rec, info=info), open(f"job50_data/run_{n}{SUF}.pkl", "wb"))
        print(f"[{n}] records {len(rec)} skies {len(info)}", flush=True)
    print(f"TOTAL {time.time()-t0:.0f} s")
