# ================= JOB 53 -- analysis, Parts 1 and 3 (+ Figures 3, 4) ============================
import os, json
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erf
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
from job53_lib import *
from job53_common import NLIST, N4, N3, SIGS
OUT = {}
NAT = ld("main_nat"); CON = ld("main_con"); SNAT = ld("stat_nat")
rn, rc = NAT["rec"], CON["rec"]
nsk = len(np.unique([r["sky"] for r in rn]))
print(f"natural skies: {nsk} distinct populations (seeds {5000+min(r['sky'] for r in rn)}-{5000+max(r['sky'] for r in rn)}), {len(rn)} pulsars; constructed: {len(np.unique([r['sky'] for r in rc]))} skies, {len(rc)} pulsars")
cl = classes(rn); print(f"natural: pulsars per class all {len(rn)}, K>=60 {int(cl['K60'].sum())}, K>=600 {int(cl['K600'].sum())}; prior only {fmt3(prior_only(rn, 'natural'))}")
nb_ = np.array([(d["arr_bg"] >= 3).sum() for d in SNAT["info"]]); no_ = np.array([(d["arr_opt"] >= 3).sum() for d in SNAT["info"]])
nr = np.array([d["nrej"] for d in NAT["info"]])
print(f"natural: detectable per sky arr_bg>=3 mean {nb_.mean():.2f} (10/50/90% {np.percentile(nb_, [10, 50, 90])}); arr_bg_opt>=3 mean {no_.mean():.2f} (10/50/90% {np.percentile(no_, [10, 50, 90])}); candidates {np.mean([d['ncand'] for d in SNAT['info']]):.0f}; upper-limit rejections per sky mean {nr.mean():.2f}, max {nr.max()}")
OUT["det"] = dict(arr_bg=float(nb_.mean()), arr_opt=float(no_.mean()), nrej=float(nr.mean()))
bu = np.concatenate([r["d_bun_on"] for r in rn]); OUT["bun"] = [float(x) for x in np.percentile(bu, [10, 50, 90])]
print(f"natural: beta_unmod N=12 pairs ({len(bu)}) 10/50/90%: {np.percentile(bu, [10, 50, 90]).round(4)}; share of band power top N: " + ", ".join(f"{n}: {np.mean([r[f'sh|{n}'] for r in rn]):.3f}" for n in NLIST))


def tab(rec, st, keys, title):
    print(f"\n-- {title}  (all / K>=60 / K>=600 [95% interval, populations resampled] | right per sky)")
    for lab, key in keys:
        row, rps = table_row(rec, key, st); print(f"   {lab:<34s}: {fmt3(row)} | {rps:.2f}"); OUT.setdefault("tab", {})[f"{st}|{key}"] = dict(row=row, per_sky=rps)


# ---------------- Part 1 ---------------------------------------------------------------------------
print("\n================ Part 1: what retardation costs ================")
AR = (("off", "off"), ("kn", "known only"), ("ot", "others only"), ("ref", "on (reference)"))
for n in N3:
    tab(rn, "natural", [(f"N={n} {lab}", f"w|{v}|{n}|t1") for v, lab in AR], f"natural skies, N = {n}")
tab(rc, "constructed", [(f"N=12 {lab}", f"w|{v}|12|t1") for v, lab in AR], "constructed skies (100-119), N = 12")
print("\n-- attribution (all pulsars, natural): fraction right")
for n in N3:
    f = {v: float(np.mean(~np.array([r[f"w|{v}|{n}|t1"] for r in rn]))) for v, _ in AR}
    loss = f["off"]-f["ref"]
    print(f"   N={n}: off {f['off']:.4f}, known only {f['kn']:.4f}, others only {f['ot']:.4f}, on {f['ref']:.4f}; total loss {loss:.4f}; "
          f"own term alone removes {(f['off']-f['kn'])/loss:.2f} of it, others alone {(f['off']-f['ot'])/loss:.2f}; own term given others {(f['ot']-f['ref'])/loss:.2f}, others given own {(f['kn']-f['ref'])/loss:.2f}; factor off/on {f['off']/f['ref']:.2f}")
    OUT.setdefault("attr", {})[str(n)] = f
print("\n-- per-pair bookkeeping, N = 12 pairs (pulsar, known source), on / off")
for st, rec in (("natural", rn), ("constructed", rc)):
    own = np.concatenate([r["d_own_on"]/r["d_own_off"] for r in rec]); unm = np.concatenate([r["d_unm_on"]/r["d_unm_off"] for r in rec])
    bb = np.concatenate([r["d_bun_on"]/r["d_bun_off"] for r in rec]); mv = np.concatenate([r["d_mov"] for r in rec])
    bon = np.concatenate([r["d_bun_on"] for r in rec]); bof = np.concatenate([r["d_bun_off"] for r in rec])
    q = lambda x: np.percentile(x, [10, 50, 90]).round(3)
    print(f"   [{st}] {len(own)} pairs: own pulsar-term power on/off 10/50/90% {q(own)}; unmodelled power it competes with on/off {q(unm)}; beta_unmod on/off {q(bb)}; bins moved {q(mv)}")
    print(f"      beta_unmod off 10/50/90% {q(bof)}; on {q(bon)}; pairs moved > 1 bin {np.mean(mv > 1):.3f}, > 0.1 bin {np.mean(mv > 0.1):.3f}")
    print(f"      log decomposition of the median drop in beta/(1-beta): own power x{np.median(own):.3f}, competing power x{np.median(unm):.3f}; geometric means: own x{np.exp(np.mean(np.log(own))):.3f}, competing x{np.exp(np.mean(np.log(unm))):.3f}")
    OUT.setdefault("pair", {})[st] = dict(own=q(own).tolist(), unm=q(unm).tolist(), bun=q(bb).tolist(), mov=q(mv).tolist(), gm_own=float(np.exp(np.mean(np.log(own)))), gm_unm=float(np.exp(np.mean(np.log(unm)))))
print("\n-- composition of the unmodelled power where a known pulsar term lands (retardation on, N = 12; share of the total, median over pairs / mean over pairs)")
NM = ("Earth terms of outside sources", "outside pulsar terms from > 1 bin higher", "outside pulsar terms moved <= 1 bin", "sea", "timing noise")
for st, rec in (("natural", rn), ("constructed", rc)):
    cmp = np.concatenate([r["d_cmp"] for r in rec]).astype(float); fP = np.concatenate([r["d_fP"] for r in rec])*1e9
    sh = cmp/cmp.sum(1, keepdims=True)
    for lab, m in (("all pairs", np.ones(len(fP), bool)), ("below 5 nHz", fP < 5), ("5-10 nHz", (fP >= 5) & (fP < 10)), ("10-20 nHz", (fP >= 10) & (fP < 20)), ("above 20 nHz", fP >= 20)):
        if m.sum() == 0: print(f"   [{st}] {lab}: no pairs"); continue
        print(f"   [{st}] {lab:<13s} ({int(m.sum()):6d} pairs): " + "; ".join(f"{nm} {np.median(sh[m, k]):.3f}/{sh[m, k].mean():.3f}" for k, nm in enumerate(NM)))
        OUT.setdefault("cmp", {})[f"{st}|{lab}"] = dict(median=[float(np.median(sh[m, k])) for k in range(5)], mean=[float(sh[m, k].mean()) for k in range(5)], n=int(m.sum()))

# ---------------- Part 3a --------------------------------------------------------------------------
print("\n================ Part 3a: the curve, natural skies ================")
tab(rn, "natural", [(f"reference N={n}", f"w|ref|{n}|t1") for n in NLIST], "natural, reference, timing x1")
tab(rn, "natural", [(f"Arm A N={n}", f"w|A|{n}|t1") for n in NLIST], "natural, Arm A, timing x1")
tab(rn, "natural", [(f"reference N={n} {t}", f"w|ref|{n}|{t}") for n in N4 for t in ("t1", "t10", "tinf")], "natural, reference, timing levels")
tab(rc, "constructed", [(f"reference N={n}", f"w|ref|{n}|t1") for n in NLIST], "constructed (skies 100-119), reference, timing x1")
print("\n-- Arm A / reference (all pulsars; K>=60), natural: " + "; ".join(f"N={n}: {np.mean(~np.array([r[f'w|A|{n}|t1'] for r in rn]))/np.mean(~np.array([r[f'w|ref|{n}|t1'] for r in rn])):.2f}" for n in NLIST))
print("-- largest change with timing (all pulsars), natural: " + "; ".join(f"N={n}: x10 {np.mean(~np.array([r[f'w|ref|{n}|t10'] for r in rn]))-np.mean(~np.array([r[f'w|ref|{n}|t1'] for r in rn])):+.4f}, zero {np.mean(~np.array([r[f'w|ref|{n}|tinf'] for r in rn]))-np.mean(~np.array([r[f'w|ref|{n}|t1'] for r in rn])):+.4f}" for n in N4))
# Figure 3
ncm = np.mean([r["nC"] for r in rn]); xN = [n if n != "all" else ncm for n in NLIST]
def curve(rec, st, v, t, m=None):
    pop = popid(rec, st); out = []
    for n in NLIST:
        k = f"w|{v}|{n}|{t}"
        if k not in rec[0]: out.append((np.nan,)*3); continue
        ok = ~np.array([r[k] for r in rec]); out.append(frac(ok, np.ones(len(rec), bool) if m is None else m, pop)[:3])
    return np.array(out)
fig, axs = plt.subplots(2, 1, figsize=(7.2, 8.0), sharex=True)
for ax, cn, ttl in ((axs[0], None, "all pulsars"), (axs[1], "K60", "pulsars with K ≥ 60")):
    mn = None if cn is None else classes(rn)[cn]; mc_ = None if cn is None else classes(rc)[cn]
    R_ = curve(rn, "natural", "ref", "t1", mn); A_ = curve(rn, "natural", "A", "t1", mn); Cc = curve(rc, "constructed", "ref", "t1", mc_)
    ax.axvspan(nb_.mean(), no_.mean(), color="0.85", lw=0, label="detectable against the background (natural)")
    ax.fill_between(xN, R_[:, 0], A_[:, 0], color="#1f77b4", alpha=0.15, lw=0)
    ax.errorbar(xN, R_[:, 0], yerr=[R_[:, 0]-R_[:, 1], R_[:, 2]-R_[:, 0]], fmt="o-", color="#1f77b4", ms=4, capsize=2, label="reference (timing ×1)")
    ax.errorbar(xN, A_[:, 0], yerr=[A_[:, 0]-A_[:, 1], A_[:, 2]-A_[:, 0]], fmt="s-", color="#2ca02c", ms=4, capsize=2, label="unknown Earth terms removed perfectly (Arm A)")
    for t, mk, lb in (("t10", "^", "reference, timing ×10 better"), ("tinf", "v", "reference, no timing noise")):
        T_ = curve(rn, "natural", "ref", t, mn); ax.plot(xN, T_[:, 0], mk, color="#ff7f0e" if t == "t10" else "#d62728", ms=6, mfc="none", label=lb)
    ax.plot(xN, Cc[:, 0], "--", color="0.3", lw=0.9, label="constructed skies (Jobs 49–52), reference")
    po = float(prior_only(rn, "natural")[0 if cn is None else 1]); ax.axhline(po, color="k", ls=":", lw=1, label="prior only")
    ax.set_xscale("log"); ax.set_ylabel(f"fraction right, {ttl}"); ax.set_ylim(0, None)
    if cn is None:
        ax2 = ax.secondary_yaxis("right", functions=(lambda y: y*79, lambda y: y/79)); ax2.set_ylabel("pulsars right per sky (of 79)")
        ax.legend(fontsize=7, frameon=False, loc="upper left")
axs[1].set_xlabel("number of sources whose parameters are granted, N (last point: every candidate)")
plt.tight_layout(); plt.savefig("job53_fig3_curve.png", dpi=160); plt.close()

# ---------------- Part 3b --------------------------------------------------------------------------
print("\n================ Part 3b: fraction right against prior width, natural skies ================")
SG = ld("sig_nat")["rec"]; pop = popid(SG, "natural"); one = np.ones(len(SG), bool)
pri = np.array([r["pri"] for r in SG]); Ltu = np.array([r["Ltrue"] for r in SG])/pri
def crossx(y, level):
    x = np.array(SIGS)
    for a in range(len(x)-1):
        if y[a] >= level >= y[a+1]:
            f = (y[a]-level)/max(y[a]-y[a+1], 1e-12); return float(10**(np.log10(x[a])+f*(np.log10(x[a+1])-np.log10(x[a]))))
    return np.nan
po = [float(np.mean(np.abs(Ltu*s) <= TW)) for s in SIGS]
print("   prior only (measured; analytic erf(0.25/(sqrt2 sigma))): " + "  ".join(f"{s:g} pc: {p:.3f} ({erf(0.25/(np.sqrt(2)*s)):.3f})" for s, p in zip(SIGS, po)))
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True); cols = {7: "#1f77b4", 12: "#ff7f0e", 30: "#2ca02c", "all": "#d62728"}
for ax, v, ttl in ((axs[0], "ref", "reference"), (axs[1], "A", "unknown Earth terms removed perfectly (Arm A)")):
    for n in N4:
        F = np.array([frac(~np.array([r[f"w|{v}|{n}|t1|s{s:g}"] for r in SG]), one, pop)[:3] for s in SIGS])
        K = [np.median([r[f"Kco|{n}|s{s:g}"] for r in SG]) for s in SIGS]
        print(f"   {v:>3s} N={str(n):>3s}: " + "  ".join(f"{s:g}: {f[0]:.3f} [{f[1]:.3f},{f[2]:.3f}]" for s, f in zip(SIGS, F)) + f" | crosses 0.9 at {crossx(F[:, 0], 0.9):.2f} pc, 0.5 at {crossx(F[:, 0], 0.5):.2f} pc, 0.1 at {crossx(F[:, 0], 0.1):.2f} pc")
        OUT.setdefault("sig", {})[f"{v}|{n}"] = dict(frac=F.tolist(), c90=crossx(F[:, 0], 0.9), c50=crossx(F[:, 0], 0.5), c10=crossx(F[:, 0], 0.1), Kco=[float(k) for k in K])
        ax.errorbar(SIGS, F[:, 0], yerr=[F[:, 0]-F[:, 1], F[:, 2]-F[:, 0]], fmt="o-", color=cols[n], ms=4, capsize=2, label=f"N = {n}" if n != "all" else "every candidate")
    ax.plot(SIGS, po, "k:", label="prior only"); ax.set_xscale("log"); ax.set_title(ttl, fontsize=10); ax.set_xlabel("prior width σ_d of every pulsar [pc]")
    for yy in (0.9, 0.5, 0.1): ax.axhline(yy, color="0.8", lw=0.6, zorder=0)
axs[0].set_ylabel("fraction right (natural skies)"); axs[0].legend(fontsize=8, frameon=False)
print("   median K (coarsest tooth of the known set) per sigma, N=12: " + "  ".join(f"{s:g}: {k:.1f}" for s, k in zip(SIGS, OUT["sig"]["ref|12"]["Kco"])))
OUT["sig_prior_only"] = po
plt.tight_layout(); plt.savefig("job53_fig4_prior_width.png", dpi=160); plt.close()

# ---------------- Part 3c --------------------------------------------------------------------------
print("\n================ Part 3c: the Job 52 Part 3 gate, unchanged (train = even skies, held out = odd) ================")
def smooth7(y): return np.array([y[max(i-3, 0):i+4].mean() for i in range(len(y))])
G = np.arange(0.5, 100.01, 0.5)
def gate(rec, mk, wk, lab, Kkey=None):
    sk = np.array([r["sky"] for r in rec]); tr = sk % 2 == 0; nho = len(np.unique(sk[~tr]))
    m = np.array([r[mk] for r in rec]); ok = ~np.array([r[wk] for r in rec])
    o = np.argsort(-m[tr]); okt = ok[tr][o]; mt = m[tr][o]
    pur = np.array([okt[:max(int(round(g/100*len(okt))), 1)].mean() for g in G]); ps = smooth7(pur)
    if not (ps >= 0.90).any():
        print(f"   {lab}: no 0.90 catalogue (best smoothed train purity {ps.max():.3f}; ungated {ok[~tr].sum()/nho:.1f} right per sky)"); OUT.setdefault("gate", {})[lab] = dict(none=True, max_train=float(ps.max())); return
    gi = np.where(ps >= 0.90)[0][-1]; thr = mt[max(int(round(G[gi]/100*len(okt))), 1)-1]
    sel = (~tr) & (m >= thr); pu = frac(ok, sel, sk); us = np.unique(sk[~tr]); c = np.array([sel[sk == s].sum() for s in us])
    rg = np.random.default_rng(5); cb = [c[ii].mean() for ii in rg.integers(0, len(us), (1000, len(us)))]
    extra = ""
    if Kkey is not None and sel.sum():
        K = np.array([r[Kkey] for r in rec]); extra = f"; catalogue K<3 {np.mean(K[sel] < 3):.2f}, K>=10 {np.mean(K[sel] >= 10):.2f}, K>=60 {np.mean(K[sel] >= 60):.2f}"
    print(f"   {lab}: gate keeps {G[gi]:.1f}% of train; held-out catalogue {sel.sum()/nho:.2f} [{np.percentile(cb, 2.5):.2f}, {np.percentile(cb, 97.5):.2f}] pulsars per sky at purity {pu[0]:.3f} [{pu[1]:.3f}, {pu[2]:.3f}] (ungated {ok[~tr].sum()/nho:.1f} right per sky){extra}")
    OUT.setdefault("gate", {})[lab] = dict(gate=float(G[gi]), per_sky=float(sel.sum()/nho), purity=pu[:3])
for v in ("ref", "A"):
    for n in N3:
        gate(rn, f"m|{v}|{n}|t1", f"w|{v}|{n}|t1", f"census priors, {v}, N={n}", "Kc")
for s in (1, 3):
    for v in ("ref", "A"):
        for n in N3:
            gate(SG, f"m|{v}|{n}|t1|s{s}", f"w|{v}|{n}|t1|s{s}", f"sigma_d = {s} pc, {v}, N={n}", f"Kco|{n}|s{s}")
json.dump(OUT, open("job53_data/main.json", "w"), indent=1)
print("\nwrote job53_data/main.json, job53_fig3_curve.png, job53_fig4_prior_width.png")
