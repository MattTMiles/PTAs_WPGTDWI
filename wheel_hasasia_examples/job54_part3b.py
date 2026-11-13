# JOB 54 Part 3, closing the gap: per-bin model with the share replaced by the inverse-variance combination over the
# bins the pulsar term leaks into (natural skies 100-159, N = 12)
import os, sys, pickle, json
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
import numpy as np, multiprocessing as mp
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples"); sys.path.insert(0, ".")
from job53_common import SIGS
import job54_parts34 as P
rec = pickle.load(open("job54_data/run_inf.pkl", "rb"))["rec"]
SG = {(r["sky"], r["psr"]): r for r in pickle.load(open("job53_data/run_sig_nat.pkl", "rb"))["rec"] if r["sky"] < 160}
def one(q):
    r = rec[q]; out = {}
    for var, P_ in (("mean", (1-np.clip(r["d_bun_on"].astype(float), 1e-9, 1-1e-9))/np.clip(r["d_bun_on"].astype(float), 1e-9, 1)), ("inv", 1.0/np.maximum(r["d_inf_on"].astype(float), 1e-12))):
        rg = np.random.default_rng([54, 33, q]); kap = r["d_kap_on"]; n = len(kap); C = (rg.normal(0, 1, n)+1j*rg.normal(0, 1, n))*np.sqrt(P_/2)
        for tag, sd in [("census", float(r["pri"]))]+[(f"{s:g}", float(s)) for s in SIGS]:
            Lt = r["Ltrue"]/r["pri"]*sd; out[f"{var}|{tag}"] = P.dec(np.exp(1j*kap*Lt)*(1+C), kap, 2.0/P_, sd, Lt)
    out["ratio"] = float(np.median(r["d_inf_on"]/(r["d_own_on"]/r["d_unm_on"])))
    return out
if __name__ == "__main__":
    with mp.get_context("fork").Pool(6) as pool: R = pool.map(one, range(len(rec)), chunksize=20)
    print(f"{len(rec)} pulsars (natural skies 100-159). Inverse-variance signal-to-unmodelled ratio / leakage-weighted-mean ratio, median over pairs: {np.median(np.concatenate([r['d_inf_on']/(r['d_own_on']/r['d_unm_on']) for r in rec])):.2f} (10-90%: {np.percentile(np.concatenate([r['d_inf_on']/(r['d_own_on']/r['d_unm_on']) for r in rec]), [10, 90]).round(2)})")
    print("   prior | sky result (same pulsars) | per-bin model, share from the weighted mean | difference | per-bin model, inverse-variance over bins | difference"); OUT = {}
    for tag in ["census"]+[f"{s:g}" for s in SIGS]:
        sky = np.mean([not r["w|ref|12|t1"] for r in rec]) if tag == "census" else np.mean([not SG[(r["sky"], r["psr"])][f"w|ref|12|t1|s{tag}"] for r in rec])
        a = np.mean([x["mean|"+tag] for x in R]); b = np.mean([x["inv|"+tag] for x in R])
        print(f"   {tag:>7s} | {sky:.3f} | {a:.3f} | {a-sky:+.3f} | {b:.3f} | {b-sky:+.3f}"); OUT[tag] = dict(sky=float(sky), mean=float(a), inv=float(b))
    json.dump(OUT, open("job54_data/part3b.json", "w"), indent=1)
