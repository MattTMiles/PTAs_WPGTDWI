# ================= JOB 52 -- analysis, Parts 3 to 6 (+ Figures 2, 3) =========================
import os, sys, json, pickle
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import i0e, i1e
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
OUT = {}
ld = lambda n: pickle.load(open(f"job52_data/run_{n}.pkl", "rb"))
have = lambda n: os.path.exists(f"job52_data/run_{n}.pkl")
N4 = (7, 12, 30, "all"); NL = (1, 2, 3, 5, 7, 12, 20, 30, 60, 100, 300, "all")

def bfrac(sk, ok, m, nb=1000, seed=3):
    us = np.unique(sk); rg = np.random.default_rng(seed)
    num = np.array([np.sum(ok[m & (sk == s)]) for s in us]); den = np.array([np.sum(m & (sk == s)) for s in us])
    if den.sum() == 0: return (np.nan, np.nan, np.nan)
    bs = [num[ii].sum()/max(den[ii].sum(), 1) for ii in rg.integers(0, len(us), (nb, len(us)))]
    return (float(num.sum()/den.sum()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)))

def smooth7(y):          # 7-point mean, window shrinks at the edges (no zero padding)
    return np.array([y[max(i-3, 0):i+4].mean() for i in range(len(y))])

# ---------------- Part 3 ------------------------------------------------------------------------
print("== Part 3: gate on the joint-fit margin (best peak minus the best value more than 0.25 pc away, prior included)")
print("   train = even-numbered skies, held-out = odd; 7-point edge-safe mean; scan 0.5-100% in 0.5% steps; LAST train crossing of 0.90; applied blind")
G = np.arange(0.5, 100.01, 0.5)
for p, lab in (("c", "census, skies 100-119"), ("h", "harness, skies 100-199")):
    R = ld(f"main_{p}")["rec"]; X = ld(f"mainx_{p}")["rec"]
    sk = np.array([r["sky"] for r in R]); K = np.array([r["Kc"] for r in R]); tr = sk % 2 == 0; nho = len(np.unique(sk[~tr]))
    for v, src in (("ref", R), ("D", X)):
        for n in (7, 12, "all"):
            m = np.array([r[f"m|{v}|{n}|t1"] for r in src]); ok = ~np.array([r[f"w|{v}|{n}|t1"] for r in src])
            o = np.argsort(-m[tr]); okt = ok[tr][o]; mt = m[tr][o]
            pur = np.array([okt[:max(int(round(g/100*len(okt))), 1)].mean() for g in G]); ps = smooth7(pur)
            key = f"{p}|{v}|{n}"
            if not (ps >= 0.90).any():
                print(f"   {lab} | {v:>3s} N={str(n):>3s}: no 0.90 catalogue (max smoothed train purity {ps.max():.3f}; ungated train purity {okt.mean():.3f})")
                OUT.setdefault("part3", {})[key] = dict(none=True, max_train=float(ps.max())); continue
            gi = np.where(ps >= 0.90)[0][-1]; thr = mt[max(int(round(G[gi]/100*len(okt))), 1)-1]
            sel = (~tr) & (m >= thr)
            pu = bfrac(sk, ok, sel); cnt = sel.sum()/nho
            us = np.unique(sk[~tr]); rg = np.random.default_rng(5); c = np.array([sel[sk == s].sum() for s in us])
            cb = [c[ii].mean() for ii in rg.integers(0, len(us), (1000, len(us)))]
            comp = [float(np.mean(K[sel] < 1)), float(np.mean(K[sel] < 3)), float(np.mean(K[sel] >= 10)), float(np.mean(K[sel] >= 60))] if sel.sum() else [np.nan]*4
            print(f"   {lab} | {v:>3s} N={str(n):>3s}: gate keeps {G[gi]:.1f}% of train; held-out catalogue {cnt:.2f} [{np.percentile(cb,2.5):.2f},{np.percentile(cb,97.5):.2f}] pulsars/sky at purity {pu[0]:.3f} [{pu[1]:.3f},{pu[2]:.3f}]; "
                  f"catalogue K<1 {comp[0]:.2f}, K<3 {comp[1]:.2f}, K>=10 {comp[2]:.2f}, K>=60 {comp[3]:.2f}  (ungated: {ok[~tr].sum()/nho:.1f} right/sky of 79)")
            OUT.setdefault("part3", {})[key] = dict(gate=float(G[gi]), per_sky=float(cnt), purity=pu, comp=comp)

# ---------------- Part 6 ------------------------------------------------------------------------
print("\n== Part 6: C = sum over bins ln(1 + P_known/P_unmodelled)")
rows = []     # (arm label, prior, C, Kfin, Kco, ok)
for p in "ch":
    R = ld(f"main_{p}")["rec"]; X = ld(f"mainx_{p}")["rec"]
    arms = [("ref", n, "t1", R) for n in NL]+[("ref", n, t, R) for n in N4 for t in ("t10", "tinf")] \
        + [(v, n, t, R) for v in ("A", "C") for n in N4 for t in ("t1", "tinf")] \
        + [(v, n, t, X) for v in ("B", "D") for n in N4 for t in ("t1", "tinf")]
    print(f"   [{'census' if p=='c' else 'harness'}] arm: C 10/50/90%")
    for v, n, t, src in arms:
        C = np.array([r[f"C|{v}|{n}|{t}"] for r in src]); ok = ~np.array([r[f"w|{v}|{n}|{t}"] for r in src])
        kf = np.array([r[f"Kfin|{n}"] for r in src]); kc = np.array([r[f"Kco|{n}"] for r in src])
        q = np.percentile(C, [10, 50, 90]); rows.append((f"{v}|{n}|{t}", p, C, kf, kc, ok))
        print(f"      {v:>3s} N={str(n):>3s} {t:>4s}: {q[0]:7.2f} / {q[1]:7.2f} / {q[2]:7.2f}")
        OUT.setdefault("part6_C", {})[f"{p}|{v}|{n}|{t}"] = [float(x) for x in q]
logi = lambda x, x0, s: 1.0/(1.0+np.exp(-(x-x0)/s))
edges = np.arange(-12, 30.01, 1.5); cen = 0.5*(edges[1:]+edges[:-1])
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6)); cmap = plt.get_cmap("turbo")
for a_, (kn, ki) in zip(axs, (("finest", 3), ("coarsest", 4))):
    xall = np.concatenate([r[2]-np.log(np.maximum(r[ki], 1e-9)) for r in rows]); oall = np.concatenate([r[5] for r in rows]).astype(float)
    (x0, s), _ = curve_fit(logi, xall, oall, p0=[3.0, 2.0], maxfev=20000)
    print(f"   K from the {kn} tooth of the known set, all arms pooled ({len(xall)} decodes): logistic 50% at C - ln K = {x0:.2f} nats, 90% at {x0+s*np.log(9):.2f} (width {s:.2f})")
    idx = np.digitize(xall, edges)-1
    print("      binned: " + "  ".join(f"[{edges[b]:.1f},{edges[b+1]:.1f}): {oall[idx==b].mean():.2f} (n={int((idx==b).sum())})" for b in range(len(cen)) if (idx == b).sum() >= 200))
    per = {}
    for ai, r in enumerate(rows):
        x = r[2]-np.log(np.maximum(r[ki], 1e-9)); o = r[5].astype(float)
        try:
            (a0,), _ = curve_fit(lambda xx, a: logi(xx, a, s), x, o, p0=[x0], maxfev=5000); per[(r[0], r[1])] = float(a0)
        except Exception: per[(r[0], r[1])] = np.nan
        ii = np.digitize(x, edges)-1; yy = [o[ii == b].mean() if (ii == b).sum() >= 40 else np.nan for b in range(len(cen))]
        a_.plot(cen, yy, "-", lw=0.7, color=cmap(ai/len(rows)), alpha=0.8)
    v = np.array(list(per.values())); nm = list(per.keys())
    grp = lambda f: np.nanmedian([per[k] for k in nm if f(k)])
    print(f"      per-arm 50% point (width fixed): min {np.nanmin(v):.2f}, median {np.nanmedian(v):.2f}, max {np.nanmax(v):.2f}, spread (max-min) {np.nanmax(v)-np.nanmin(v):.2f} nats")
    print(f"      medians: census {grp(lambda k: k[1]=='c'):.2f}, harness {grp(lambda k: k[1]=='h'):.2f} | ref {grp(lambda k: k[0].startswith('ref')):.2f}, A {grp(lambda k: k[0].startswith('A')):.2f}, "
          f"B {grp(lambda k: k[0].startswith('B')):.2f}, C {grp(lambda k: k[0].startswith('C|')):.2f}, D {grp(lambda k: k[0].startswith('D')):.2f} | "
          f"N<=12 {grp(lambda k: k[0].split('|')[1] in ('1','2','3','5','7','12')):.2f}, N=20-300 {grp(lambda k: k[0].split('|')[1] in ('20','30','60','100','300')):.2f}, every candidate {grp(lambda k: k[0].split('|')[1]=='all'):.2f} | "
          f"x1 {grp(lambda k: k[0].endswith('t1')):.2f}, zero noise {grp(lambda k: k[0].endswith('tinf')):.2f}")
    OUT.setdefault("part6_fit", {})[kn] = dict(x50=float(x0), x90=float(x0+s*np.log(9)), width=float(s), per_arm={f"{k[1]}|{k[0]}": val for k, val in per.items()})
    xx = np.linspace(-12, 30, 200); a_.plot(xx, logi(xx, x0, s), "k--", lw=1.5, label=f"pooled logistic: 50% at {x0:.1f}, 90% at {x0+s*np.log(9):.1f}")
    a_.set_xlabel(f"C − ln K  (K from the {kn} tooth of the known set) [nats]"); a_.set_ylabel("fraction decoded correctly"); a_.legend(fontsize=8, frameon=False)
    a_.set_title(f"{len(rows)} arms (colour = arm; census and harness)")
plt.tight_layout(); plt.savefig("job52_fig3_information.png", dpi=150); plt.close()

# ---------------- Parts 4 and 5 (need the p5 runs) ------------------------------------------------
if have("p5_c") and have("p5off_c"):
    Pn = ld("p5_c")["rec"]; Po = ld("p5off_c")["rec"]
    def rice(beta):      # exact <cos dphi> and rms dphi for phasor 1 + n, E|n|^2 = (1-beta)/beta
        s2 = (1-beta)/beta; k = 1.0/s2                       # k = A^2/sigma^2 with sigma^2 = E|n|^2
        c = np.sqrt(np.pi*k)/2*(i0e(k/2)+i1e(k/2))           # <cos> = sqrt(pi k)/2 e^{-k/2}[I0(k/2)+I1(k/2)]
        rg = np.random.default_rng(1); z = 1+np.sqrt(s2/2)*(rg.normal(size=200000)+1j*rg.normal(size=200000))
        return c, float(np.sqrt(np.mean(np.angle(z)**2)))
    print("\n== Part 4: phase error measured in the bins, N = 12, census priors, reference configuration")
    print("   (pulsar term spread over bins: amplitude-weighted sum over bins of datum x conj(true phasor))")
    P4 = {}
    for nm, rec in (("retardation on", Pn), ("retardation off", Po)):
        for tl in ("t1", "tinf"):
            bu = np.concatenate([r["p4_bun_"+tl] for r in rec]); bb = np.concatenate([r["p4_bbin_"+tl] for r in rec])
            print(f"   {nm}, timing {'x1' if tl=='t1' else 'zero'}: {len(bu)} pairs ({len(bu)/len(rec):.1f}/pulsar); beta_bin 10/50/90%: {np.percentile(bb,[10,50,90]).round(4)}; beta_unmod 10/50/90%: {np.percentile(bu,[10,50,90]).round(4)}")
            P4[f"{nm}|{tl}|pct"] = dict(bbin=[float(x) for x in np.percentile(bb, [10, 50, 90])], bun=[float(x) for x in np.percentile(bu, [10, 50, 90])])
    be = np.geomspace(0.003, 0.99, 9)
    figd = {}
    for tl in ("t1", "tinf"):
        bu = np.concatenate([r["p4_bun_"+tl] for r in Pn]).astype(float); dp = np.concatenate([r["p4_dphi_"+tl] for r in Pn]).astype(float)
        print(f"   retardation on, timing {'x1' if tl=='t1' else 'zero'}: bins of beta_unmod | n pairs | measured rms dphi | exact rms | small-angle sqrt((1-b)/(2b)) | measured <cos> +- se | exact <cos>")
        tab = []
        for a, b in zip(be[:-1], be[1:]):
            m = (bu >= a) & (bu < b)
            if m.sum() < 5: print(f"      [{a:.3f},{b:.3f}): n={int(m.sum())}"); continue
            bm = float(np.median(bu[m])); c_ex, r_ex = rice(bm)
            cm = float(np.mean(np.cos(dp[m]))); se = float(np.std(np.cos(dp[m]))/np.sqrt(m.sum())); rm = float(np.sqrt(np.mean(dp[m]**2)))
            # exact value averaged over the betas in the bin
            cex = float(np.mean([rice(x)[0] for x in np.quantile(bu[m], np.linspace(0.05, 0.95, 10))]))
            tab.append((a, b, int(m.sum()), bm, rm, r_ex, float(np.sqrt((1-bm)/(2*bm))), cm, se, cex))
            print(f"      [{a:.3f},{b:.3f}) | {int(m.sum()):6d} | {rm:.3f} | {r_ex:.3f} | {np.sqrt((1-bm)/(2*bm)):.3f} | {cm:.3f} +- {se:.3f} | {cex:.3f}")
        P4[tl] = tab; figd[tl] = tab
    print(f"   below 0.003: {np.mean(np.concatenate([r['p4_bun_t1'] for r in Pn]) < 0.003):.3f} of pairs; above 0.99: {np.mean(np.concatenate([r['p4_bun_t1'] for r in Pn]) > 0.99):.3f}")
    OUT["part4"] = P4
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
    bx = np.geomspace(0.003, 0.99, 100); a1.plot(bx, [rice(x)[0] for x in bx], "k-", label="exact (Rice)")
    for tl, col, lb in (("t1", "#1f77b4", "measured, timing ×1"), ("tinf", "#2ca02c", "measured, timing noise ÷ 10⁶")):
        t = figd[tl]; a1.errorbar([x[3] for x in t], [x[7] for x in t], yerr=[2*x[8] for x in t], fmt="o", color=col, ms=4, capsize=2, label=lb)
    a1.set_xscale("log"); a1.set_ylabel("⟨cos δφ⟩"); a1.legend(frameon=False, fontsize=8); a1.set_title("N = 12 known, census priors: phase error in the bins")
    hb = np.geomspace(1e-4, 1, 50)
    a2.hist(np.concatenate([r["p4_bun_t1"] for r in Pn]), bins=hb, histtype="step", color="#d62728", label="retardation on")
    a2.hist(np.concatenate([r["p4_bun_t1"] for r in Po]), bins=hb, histtype="step", color="#7f7f7f", label="retardation off")
    a2.set_xlabel("β_unmod (pair's pulsar-term power ÷ (itself + unmodelled power))"); a2.set_ylabel("pairs"); a2.legend(frameon=False, fontsize=8); a2.set_xlim(1e-4, 1)
    plt.tight_layout(); plt.savefig("job52_fig2_phase_error.png", dpi=150); plt.close()

    print("\n== Part 5: one change at a time, census priors, skies 100-119, x1 (fractions right: all / K>=60 / K>=600)")
    sk = np.array([r["sky"] for r in Pn]); K = np.array([r["Kc"] for r in Pn]); P5 = {}
    def trip(rec, key):
        ok = ~np.array([r[key] for r in rec]); return [bfrac(sk, ok, np.ones(len(K), bool)), bfrac(sk, ok, K >= 60), bfrac(sk, ok, K >= 600)]
    f3 = lambda t: f"{t[0]:.3f} [{t[1]:.3f},{t[2]:.3f}]"
    for n in (12, "all"):
        for lab, rec, key in (("reference", Pn, f"w|ref|{n}|t1"), ("retardation off", Po, f"w|ref|{n}|t1"), ("sea removed", Pn, f"w|nosea|{n}|t1"),
                              ("discrete sources outside the known set removed", Pn, f"w|nodisc|{n}|t1"), ("timing noise removed", Pn, f"w|ref|{n}|tinf"),
                              ("amplitude 3.7e-15", Pn, f"w|ref|{n}|a37"), ("amplitude 10.6e-15", Pn, f"w|ref|{n}|a106"),
                              ("known Earth terms left in", Pn, f"w|earthin|{n}|t1")):
            t = trip(rec, key); P5[f"{n}|{lab}"] = t
            print(f"   N={str(n):>3s} | {lab:48s} | {f3(t[0])} | {f3(t[1])} | {f3(t[2])}")
    OUT["part5"] = P5
json.dump(OUT, open("job52_data/analysis_rest.json", "w"), indent=1, default=float)
print("\nwrote job52_data/analysis_rest.json, job52_fig3_information.png" + (", job52_fig2_phase_error.png" if have("p5_c") else ""))
