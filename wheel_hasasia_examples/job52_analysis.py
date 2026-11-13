# ================= JOB 52 -- analysis, Parts 1 and 2 (+ Figure 1) ==========================
import os, sys, json, pickle
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
OUT = {}
ld = lambda n: pickle.load(open(f"job52_data/run_{n}.pkl", "rb"))
NL = (1, 2, 3, 5, 7, 12, 20, 30, 60, 100, 300, "all"); N4 = (7, 12, 30, "all")
PL = {"c": "census priors, skies 100-119", "h": "harness priors, skies 100-199"}

def frac(rec, wrong, mask=None, nb=1000, seed=3):
    sk = np.array([r["sky"] for r in rec]); w = np.asarray(wrong, bool)
    m = np.ones(len(rec), bool) if mask is None else mask
    if m.sum() == 0: return (np.nan, np.nan, np.nan, 0)
    us = np.unique(sk); rg = np.random.default_rng(seed)
    num = np.array([np.sum(~w[m & (sk == s)]) for s in us]); den = np.array([np.sum(m & (sk == s)) for s in us])
    bs = [num[ii].sum()/max(den[ii].sum(), 1) for ii in rg.integers(0, len(us), (nb, len(us)))]
    return (float(num.sum()/den.sum()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), int(m.sum()))
f3 = lambda t: f"{t[0]:.3f} [{t[1]:.3f},{t[2]:.3f}]"

def row(rec, key):
    K = np.array([r["Kc"] for r in rec]); w = np.array([r[key] for r in rec])
    nsky = len({r["sky"] for r in rec}); ok = ~w
    return dict(all=frac(rec, w), k60=frac(rec, w, K >= 60), k600=frac(rec, w, K >= 600),
                per_sky=float(ok.sum()/nsky), per_sky_K3=float((ok & (K >= 3)).sum()/nsky), per_sky_K10=float((ok & (K >= 10)).sum()/nsky))

for p in "ch":
    R = ld(f"main_{p}")["rec"]; X = ld(f"mainx_{p}")["rec"]; S = ld(f"stat_{p}")["info"]; SX = ld(f"statx_{p}")["info"]
    assert [(r["sky"], r["psr"]) for r in R] == [(r["sky"], r["psr"]) for r in X]
    K = np.array([r["Kc"] for r in R]); Lt = np.array([r["Ltrue"] for r in R]); nsky = len(S)
    pw = np.abs(Lt) > 0.25
    po = dict(all=frac(R, pw), k60=frac(R, pw, K >= 60), k600=frac(R, pw, K >= 600))
    print(f"\n================ {PL[p]}: {len(R)} pulsars ({len(R)/nsky:.1f}/sky); K>=60: {int((K>=60).sum())}, K>=600: {int((K>=600).sum())} ================")
    print(f"   prior only (decode = prior centre): all {po['all'][0]:.3f} | K>=60 {po['k60'][0]:.3f} | K>=600 {po['k600'][0]:.3f}")
    # ranking facts
    in12 = []; nbg = []; nopt = []
    for d in S:
        rk = np.argsort(np.argsort(-d["arr_opt"])); in12.append(int((rk[d["is12"]] < 12).sum()))
        nbg.append(int((d["arr_bg"] >= 3).sum())); nopt.append(int((d["arr_opt"] >= 3).sum()))
    n12 = np.mean([d["is12"].sum() for d in S])
    print(f"   ranking by arr_bg_opt (fixed amplitude): booked sources in the top 12: mean {np.mean(in12):.1f} of {n12:.1f} booked per sky")
    print(f"   sources per sky with arr_bg >= 3: {np.mean(nbg):.1f}; with arr_bg_opt >= 3: {np.mean(nopt):.1f}; candidates per sky {np.mean([len(d['arr_opt']) for d in S]):.0f}")
    OUT[p] = dict(prior=po, in12=float(np.mean(in12)), nbooked=float(n12), n_bg3=float(np.mean(nbg)), n_opt3=float(np.mean(nopt)), part1={}, part2={})
    print(f"\n-- Part 1 [{PL[p]}; reference configuration; redshifted masses; fixed amplitude]")
    print("   N known | timing | all                 | K>=60               | K>=600              | right/sky (K>=3, K>=10) | share of band GW power")
    for tl in ("t1", "t10", "tinf"):
        for n in (NL if tl == "t1" else N4):
            key = f"w|ref|{n}|{tl}"; r = row(R, key); sh = float(np.median([x[f"sh|{n}"] for x in R]))
            r["share"] = sh; OUT[p]["part1"][f"{n}|{tl}"] = r
            print(f"   {str(n):>7s} | {tl:>5s}  | {f3(r['all'])} | {f3(r['k60'])} | {f3(r['k600'])} | {r['per_sky']:5.1f} ({r['per_sky_K3']:.1f}, {r['per_sky_K10']:.1f}) | {sh:.3f}")
    print(f"\n-- Part 2 [{PL[p]}]")
    print("   arm | N known | timing | all                 | K>=60               | K>=600")
    for n in N4:
        for tl in ("t1", "tinf"):
            for v, src in (("ref", R), ("A", R), ("B", X), ("C", R), ("D", X)):
                r = row(src, f"w|{v}|{n}|{tl}"); OUT[p]["part2"][f"{v}|{n}|{tl}"] = r
                print(f"   {v:>3s} | {str(n):>7s} | {tl:>5s}  | {f3(r['all'])} | {f3(r['k60'])} | {f3(r['k600'])}")
    g = lambda rec, k: float(np.nanmedian([x[k] for x in rec]))
    print("   Arm A: share of unmodelled GW power that was Earth-term (other sources' Earth terms + half the sea), median over bins and pulsars:")
    for n in N4:
        print(f"      N={n}: bins <= 10 nHz {g(R, f'fe_lo|{n}'):.3f}; bins 15-33 nHz {g(R, f'fe_hi|{n}'):.3f}; extended bins 15-40 nHz {g(X, f'fe_hi|{n}'):.3f}; all bins 2-30 nHz {g(R, f'fe_all|{n}'):.3f}")
    OUT[p]["armA_fe"] = {str(n): [g(R, f"fe_lo|{n}"), g(R, f"fe_hi|{n}"), g(X, f"fe_hi|{n}"), g(R, f"fe_all|{n}")] for n in N4}
    # Arm B counts
    hi_i = []; hi_bg = []; hi_opt = []; fhi = []; fphi = []; extra = []
    for d, d0 in zip(SX, S):
        m = d["fc"] > 3e-8; hi_i.append(int(m.sum())); hi_bg.append(int((d["arr_bg"][m] >= 3).sum())); hi_opt.append(int((d["arr_opt"][m] >= 3).sum()))
        fhi.append(d["fc"][m]); fphi.append(d["fP0"][m]); extra.append(len(d["fc"])-len(d0["fc"]))
    fhi = np.concatenate(fhi) if hi_i and sum(hi_i) else np.zeros(0); fphi = np.concatenate(fphi) if sum(hi_i) else np.zeros(0)
    print(f"   Arm B: sources per sky with Earth frequency > 30 nHz passing (i) booked cut vs timing noise: mean {np.mean(hi_i):.2f} (max {np.max(hi_i)}); (ii) arr_bg >= 3: {np.mean(hi_bg):.2f}; (iii) arr_bg_opt >= 3: {np.mean(hi_opt):.2f}; pool grows by {np.mean(extra):.2f} sources/sky")
    if len(fhi):
        print(f"          their Earth frequencies [nHz] 10/50/90%: {np.percentile(fhi*1e9,[10,50,90]).round(1)}; pulsar-term frequencies for the first pulsar of each sky: {np.percentile(fphi*1e9,[10,50,90]).round(1)}")
    for n in N4:
        print(f"          share of each bin's GW power held by the known set, bins 15-40 nHz (median over bins and pulsars): N={n}: {g(X, f'shhi|{n}'):.4f} (reference bins 15-33 nHz: {g(R, f'shhi|{n}'):.4f})")
    OUT[p]["armB"] = dict(n_i=float(np.mean(hi_i)), n_bg=float(np.mean(hi_bg)), n_opt=float(np.mean(hi_opt)),
                          shhi={str(n): [g(X, f"shhi|{n}"), g(R, f"shhi|{n}")] for n in N4})
    mv = lambda k: float(np.mean([x[k] for x in R]))
    print(f"   Arm C: (pulsar, source) pairs whose retarded frequency moves > 0.1 bin across +-3 sigma_d: every candidate {mv('c_mov01_all'):.3f}; top 30 {mv('c_mov01_30'):.3f}; top 12 {mv('c_mov01_12'):.3f}; "
          f"largest move median {np.median([x['c_dfmax'] for x in R]):.2f} bins; segment cap ({4000}) binds for {sum(x['c_segcap'] for x in R)} of {len(R)} pulsars (reg), {sum(x['c_segcap'] for x in X)} (ext)")
    OUT[p]["armC"] = dict(all=mv("c_mov01_all"), top30=mv("c_mov01_30"), top12=mv("c_mov01_12"))
    print(f"   grid: size cap binds for {sum(x['j_cap'] for x in R)} of {len(R)} pulsars; max step {max(x['j_dL'] for x in R):.3f} pc")
    # class table, finest tooth of the booked list
    Kf = np.array([r["Kf"] for r in R]); edges = [0, 1, 3, 10, 60, 600, np.inf]; labs = ["K<1", "1-3", "3-10", "10-60", "60-600", ">=600"]
    print(f"   by class of K (finest tooth of the booked source list), fraction right [95%], n pulsars:")
    OUT[p]["class"] = {}
    for v, src, n, tl in (("ref", R, "all", "t1"), ("D", X, "all", "t1"), ("D", X, "all", "tinf"), ("ref", R, 12, "t1"), ("D", X, 12, "t1"), ("ref", R, 7, "t1"), ("D", X, 7, "t1")):
        w = np.array([r[f"w|{v}|{n}|{tl}"] for r in src]); cells = []
        for a, b, lb in zip(edges[:-1], edges[1:], labs):
            m = (Kf >= a) & (Kf < b); fr = frac(src, w, m); cells.append(fr)
        OUT[p]["class"][f"{v}|{n}|{tl}"] = cells
        print(f"      {v:>3s} N={str(n):>3s} {tl:>4s}: " + " | ".join(f"{lb} {c[0]:.3f} [{c[1]:.2f},{c[2]:.2f}] n={c[3]}" for lb, c in zip(labs, cells)))
    pc = [frac(R, pw, (Kf >= a) & (Kf < b)) for a, b in zip(edges[:-1], edges[1:])]
    print("      prior only     : " + " | ".join(f"{lb} {c[0]:.3f}" for lb, c in zip(labs, pc)))
    OUT[p]["class"]["prior"] = pc

# ---------------- Figure 1 -----------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
p = "c"; P1 = OUT[p]["part1"]; nall = 1160
xs = [n if n != "all" else nall for n in NL]
for a_, key, ttl in ((ax[0], "all", "all pulsars"), (ax[1], "k60", "K ≥ 60")):
    for tl, col, lab in (("t1", "#1f77b4", "timing ×1"), ("t10", "#ff7f0e", "timing ×10"), ("tinf", "#2ca02c", "timing noise ÷ 10⁶")):
        nn = NL if tl == "t1" else N4
        x = [n if n != "all" else nall for n in nn]; y = [P1[f"{n}|{tl}"][key] for n in nn]
        a_.errorbar(x, [v[0] for v in y], yerr=[[v[0]-v[1] for v in y], [v[2]-v[0] for v in y]], marker="o", ms=4, lw=1.3, color=col, label=lab, capsize=2)
    a_.axhline(OUT[p]["prior"][key][0], ls=":", color="k", lw=1, label="prior only")
    a_.axvspan(OUT[p]["n_bg3"], OUT[p]["n_opt3"], color="0.85", label="detectable against background")
    a_.set_xscale("log"); a_.set_xlabel("number of sources granted as known (ranked by array S/N against background)")
    a_.set_ylabel("fraction of pulsar distances recovered"); a_.set_title(f"census priors, {ttl}"); a_.set_ylim(0, None)
ax[0].legend(fontsize=8, frameon=False)
plt.tight_layout(); plt.savefig("job52_fig1_fraction_vs_N.png", dpi=150); plt.close()
json.dump(OUT, open("job52_data/analysis_main.json", "w"), indent=1, default=float)
print("\nwrote job52_fig1_fraction_vs_N.png, job52_data/analysis_main.json")
