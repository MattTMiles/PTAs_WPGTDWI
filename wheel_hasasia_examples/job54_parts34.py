# ================= JOB 54 Parts 3 and 4 ===============================================================
# Part 3: does the one-source-per-bin model, fed the real shares, reproduce the sky result?
# Part 4: the loss from retardation split into "fainter" and "moved"; Earth against pulsar terms; is the
#         destination bin louder because of the migration?
import os, sys, json, pickle
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
import numpy as np
import multiprocessing as mp
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples"); sys.path.insert(0, ".")
from job53_lib import frac, classes, TW
from job53_common import SIGS
OUT = {}
rec = pickle.load(open("job54_data/run_mainr.pkl", "rb"))["rec"]; sky = np.array([r["sky"] for r in rec]); cl = classes(rec); one = np.ones(len(rec), bool)
fr = lambda key, m=None: frac(~np.array([r[key] for r in rec]), one if m is None else m, sky)
T45 = 45*365.25*86400

# ---------------- Part 3 ---------------------------------------------------------------------------------
def dec(z, kap, w, sd, Lt):
    dL = min(0.05, (2*np.pi/float(kap.max()))/6.0); M = int(np.floor(6*sd/dL))+1; best = -np.inf; Lb = 0.0
    for c0 in range(0, M, 16384):
        L = -3*sd+dL*np.arange(c0, min(c0+16384, M))
        S = (w[:, None]*np.real(z[:, None]*np.exp(-1j*kap[:, None]*L[None, :]))).sum(0)-L**2/(2*sd**2)
        k = int(np.argmax(S))
        if S[k] > best: best = S[k]; Lb = L[k]
    return abs(Lb-Lt) <= TW

def one_psr(a):
    q, seed = a; r = rec[q]; rg = np.random.default_rng([54, 3, seed, q]); kap = r["d_kap_on"]; b = np.clip(r["d_bun_on"].astype(float), 1e-9, 1-1e-9); n = len(kap)
    P = (1-b)/b; C = (rg.normal(0, 1, n)+1j*rg.normal(0, 1, n))*np.sqrt(P/2)
    own = r["d_own_on"].astype(float); unm = r["d_unm_on"].astype(float); fP = r["d_fP_on"]
    # variant "shared": known sources landing in the same bin share one draw of the unmodelled power there
    bi = np.rint(fP*T45).astype(int); ub, inv = np.unique(bi, return_inverse=True)
    Pb = np.array([unm[inv == g].mean() for g in range(len(ub))]); nbz = (rg.normal(0, 1, len(ub))+1j*rg.normal(0, 1, len(ub)))*np.sqrt(Pb/2)
    A = np.zeros((len(ub), n)); A[inv, np.arange(n)] = np.sqrt(own)
    res = {}
    for tag, sd in [("census", float(r["pri"]))]+[(f"{s:g}", float(s)) for s in SIGS]:
        Lt = r["Ltrue"]/r["pri"]*sd
        res["eq|"+tag] = dec(np.exp(1j*kap*Lt)*(1+C), kap, 2.0/P, sd, Lt)
        # shared-bin joint fit: S = sum_b [2 Re(d_b conj t_b) - |t_b|^2]/P_b
        d = A @ np.exp(1j*kap*Lt)+nbz
        dL = min(0.05, (2*np.pi/float(kap.max()))/6.0); M = int(np.floor(6*sd/dL))+1; best = -np.inf; Lb = 0.0
        for c0 in range(0, M, 16384):
            L = -3*sd+dL*np.arange(c0, min(c0+16384, M)); Tm = A @ np.exp(1j*kap[:, None]*L[None, :])
            S = ((2*np.real(np.conj(Tm)*d[:, None])-np.abs(Tm)**2)/Pb[:, None]).sum(0)-L**2/(2*sd**2); k = int(np.argmax(S))
            if S[k] > best: best = S[k]; Lb = L[k]
        res["sh|"+tag] = abs(Lb-Lt) <= TW
    res["nshared"] = n-len(ub)
    return res

if __name__ == "__main__":
    print("================ Part 3: the per-bin model fed the real shares (natural skies, N = 12 known sets) ================")
    rgs = np.random.default_rng(543); pick = np.sort(rgs.choice(len(rec), 6000, replace=False))
    with mp.get_context("fork").Pool(int(os.environ.get("JOB54_NW", "4"))) as pool: R = pool.map(one_psr, [(int(q), 0) for q in pick], chunksize=20)
    SG = pickle.load(open("job53_data/run_sig_nat.pkl", "rb"))["rec"]; rn = pickle.load(open("job53_data/run_main_nat.pkl", "rb"))["rec"]
    print(f"   {len(pick)} (pulsar, known set) pairs drawn from {len(rec)}; known sources sharing a bin with another known source: {np.mean([x['nshared'] for x in R]):.2f} of 12 per pulsar")
    print("   prior | sky result, all pulsars | sky result, same 6000 | per-bin model, real shares | difference | per-bin model with shared bins | difference")
    for tag in ["census"]+[f"{s:g}" for s in SIGS]:
        if tag == "census": skyv = float(np.mean([not r["w|ref|12|t1"] for r in rn])); sub = float(np.mean([not rn[q]["w|ref|12|t1"] for q in pick]))
        else: skyv = float(np.mean([not r[f"w|ref|12|t1|s{tag}"] for r in SG])); sub = float(np.mean([not SG[q][f"w|ref|12|t1|s{tag}"] for q in pick]))
        e = float(np.mean([x["eq|"+tag] for x in R])); s_ = float(np.mean([x["sh|"+tag] for x in R])); se = np.sqrt(max(e*(1-e), 1e-9)/len(pick))
        print(f"   {tag:>7s} | {skyv:.3f} | {sub:.3f} | {e:.3f} +- {se:.3f} | {e-sub:+.3f} | {s_:.3f} | {s_-sub:+.3f}")
        OUT.setdefault("part3", {})[tag] = dict(sky=skyv, sky_sub=sub, model=e, model_shared=s_, se=float(se))
    assert all(rec[q]["sky"] == rn[q]["sky"] and rec[q]["psr"] == rn[q]["psr"] for q in pick)

    # ---------------- Part 4.1 -------------------------------------------------------------------------------
    print("\n================ Part 4.1: frequency only / amplitude only (others retarded as in the reference) ================")
    for n in (7, 12, "all"):
        v = {a: fr(f"w|{a}|{n}|t1") for a in ("off", "ot", "fo", "ao", "ref")}; k = {a: fr(f"w|{a}|{n}|t1", cl["K60"]) for a in v}
        lo = v["ot"][0]-v["ref"][0]; lo2 = v["off"][0]-v["ref"][0]
        print(f"   N={n}: all pulsars: " + "; ".join(f"{a} {v[a][0]:.3f} [{v[a][1]:.3f}, {v[a][2]:.3f}]" for a in v))
        print(f"          K>=60:       " + "; ".join(f"{a} {k[a][0]:.3f}" for a in k))
        print(f"          loss (others-retarded baseline 'ot' -> on) {lo:.4f}: frequency only accounts for {(v['ot'][0]-v['fo'][0])/lo:.2f}, amplitude only {(v['ot'][0]-v['ao'][0])/lo:.2f}; against 'off' ({lo2:.4f}): {(v['off'][0]-v['fo'][0])/lo2:.2f}, {(v['off'][0]-v['ao'][0])/lo2:.2f}")
        OUT.setdefault("part41", {})[str(n)] = dict(all={a: v[a][:3] for a in v}, K60={a: k[a][:3] for a in k})
    pf = np.concatenate([r["d_pf"] for r in rec]); mv = np.concatenate([r["d_mov"] for r in rec])
    print(f"   N = 12 pairs: retarded / unretarded pulsar-term power (amplitude factor^2) 10/50/90%: {np.percentile(pf, [10, 50, 90]).round(3)}; bins moved 10/50/90%: {np.percentile(mv, [10, 50, 90]).round(2)}")

    # ---------------- Part 4.2 -------------------------------------------------------------------------------
    print("\n================ Part 4.2: Earth terms against pulsar terms in what a known pulsar term competes with (N = 12) ================")
    for st, lab in (("on", "retardation on (reference)"), ("off", "retardation off")):
        c = np.concatenate([r[f"d_cmp_{st}"] for r in rec]).astype(float); f = np.concatenate([r[f"d_fP_{st}"] for r in rec])*1e9; tot = c.sum(1)
        pul = (c[:, 1]+c[:, 2]+c[:, 3]/2)/tot; ear = (c[:, 0]+c[:, 3]/2)/tot; tim = c[:, 4]/tot
        for cls, m in (("all", np.ones(len(f), bool)), ("below 5 nHz", f < 5), ("5-10 nHz", (f >= 5) & (f < 10)), ("10-20 nHz", (f >= 10) & (f < 20)), ("above 20 nHz", f >= 20)):
            if m.sum() < 1: continue
            q = lambda x: np.percentile(x[m], [10, 50, 90]).round(3)
            print(f"   [{lab}] {cls:<13s} ({int(m.sum()):6d} pairs): pulsar terms (discrete + half the sea) 10/50/90% {q(pul)}; Earth terms (discrete + half the sea) {q(ear)}; timing {q(tim)}; of which sea halves: median {np.median((c[:, 3]/2/tot)[m]):.3f} each")
            OUT.setdefault("part42", {})[f"{st}|{cls}"] = dict(n=int(m.sum()), pulsar=q(pul).tolist(), earth=q(ear).tolist(), timing=q(tim).tolist())

    # ---------------- Part 4.3 -------------------------------------------------------------------------------
    print("\n================ Part 4.3: is the destination bin louder because of the migration? ================")
    Pon = np.array([r["b_PP_on"] for r in rec], float); Pof = np.array([r["b_PP_off"] for r in rec], float); Con = np.array([r["b_cnt_on"] for r in rec], float); Cof = np.array([r["b_cnt_off"] for r in rec], float)
    PE = np.array([r["b_PE"] for r in rec], float); sea = np.array([r["b_sea"] for r in rec], float)
    fb = np.arange(1, Pon.shape[1]+1)/T45*1e9
    rb = np.median(Pon/Pof, 0); print("   per bin, pulsar-term power of all discrete sources, on / off, median over pulsars: " + " ".join(f"{f:.1f}:{x:.2f}" for f, x in zip(fb, rb)))
    for cls, m in (("below 5 nHz", fb < 5), ("5-10 nHz", (fb >= 5) & (fb < 10)), ("10-20 nHz", (fb >= 10) & (fb < 20)), ("above 20 nHz", fb >= 20)):
        rp = Pon[:, m].sum(1)/Pof[:, m].sum(1); rc = Con[:, m].sum(1)/np.maximum(Cof[:, m].sum(1), 1)
        rt = (Pon[:, m].sum(1)+PE[:, m].sum(1)+sea[:, m].sum(1))/(Pof[:, m].sum(1)+PE[:, m].sum(1)+sea[:, m].sum(1))
        print(f"   {cls:<13s} ({int(m.sum())} bins): pulsar-term power on/off: median over pulsars {np.median(rp):.3f} (10-90% {np.percentile(rp, 10):.3f}-{np.percentile(rp, 90):.3f}); per-bin medians {rb[m].min():.2f}-{rb[m].max():.2f}; "
              f"count of candidate pulsar terms on/off: {np.median(rc):.3f} (mean count per pulsar on {Con[:, m].sum(1).mean():.0f}, off {Cof[:, m].sum(1).mean():.0f}); all gravitational-wave power on/off {np.median(rt):.3f}")
        OUT.setdefault("part43", {})[cls] = dict(power=float(np.median(rp)), count=float(np.median(rc)), total=float(np.median(rt)), bin_min=float(rb[m].min()), bin_max=float(rb[m].max()))
    json.dump(OUT, open("job54_data/parts34.json", "w"), indent=1)
