# ================= JOB 50 -- tables and figure from saved runs ========================
import os, sys, json, pickle, csv
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job50_common as C
J = C.J; J.NS.setdefault("alias_score24", None)
os.chdir(J.BASE)
OUT = {}
TW = 0.25
SUF = "" if C.MASS == "z" else "_" + C.MASS    # JOB50_MASS=src: source-frame-mass retardation runs
ld = lambda n: pickle.load(open(f"job50_data/run_{n}{'' if n == 'tau0' else SUF}.pkl", "rb"))
A12 = [C.ARMNAME(a) for a in C.ARMS12]


def frac(rec, key, mask, nb=1000, seed=3):
    """fraction correct with a sky-cluster bootstrap (per-sky counts), 5/95 %."""
    sk = np.array([r["sky"] for r in rec])[mask]
    ok = ~np.array([r[key] for r in rec], bool)[mask]
    if len(ok) == 0: return (np.nan, np.nan, np.nan, 0)
    u = np.unique(sk); per = np.array([[ok[sk == j].sum(), (sk == j).sum()] for j in u])
    rg = np.random.default_rng(seed); bs = rg.integers(0, len(u), (nb, len(u)))
    t = per[bs].sum(1); b = t[:, 0]/np.maximum(t[:, 1], 1)
    return (float(ok.mean()), float(np.quantile(b, .05)), float(np.quantile(b, .95)), int(len(ok)))

f3 = lambda t: f"{t[0]:.3f} [{t[1]:.3f}, {t[2]:.3f}]"

# ---------------- B2: tau = 0 reproduces p2_mf ---------------------------------------
R49 = pickle.load(open("job49_data/rec_marg.pkl", "rb"))
t0 = ld("tau0")["rec"]
same_all = np.mean([r["wb|on|top12|machine"] == r["wb|off|top12|machine"] for r in t0])
B0 = [r for r in t0 if r["booked"]]
same49 = sum(a["wb|on|top12|machine"] == b["wrong_p2_mf"] for a, b in zip(B0, R49["rec"]))
recb = [dict(sky=r["sky"], wrong=r["wb|on|top12|machine"]) for r in B0]
bk = J.book_sw(recb, alias=R49["alias"])
print(f"== B2 tau=0: p2r(tau=0) == p2 flags on all {len(t0)} pulsars: {same_all:.5f};"
      f" vs Job 49 p2_mf on booked records: {same49}/{len(R49['rec'])};"
      f" booked gate: {bk['y']:.2f} @ {bk['pur']:.4f} (Job 49: 36.56 @ 0.9083)")
OUT["B2"] = dict(same_all=float(same_all), same49=int(same49), n49=len(R49["rec"]), y=float(bk["y"]), pur=float(bk["pur"]))

# ---------------- Part C ---------------------------------------------------------------
def tableC(run, arms, label):
    D = ld(run); rec = [r for r in D["rec"] if r["sky"] >= 100]
    Kc = np.array([r["Kc"] for r in rec]); Kf = np.array([r["Kf"] for r in rec])
    Lt = np.array([r["Ltrue"] for r in rec])
    allm = np.ones(len(rec), bool)
    prior = {n: float((np.abs(Lt[m]) <= TW).mean()) if m.any() else np.nan
             for n, m in (("K60", Kc >= 60), ("K600", Kc >= 600), ("all", allm), ("K60f", Kf >= 60), ("K600f", Kf >= 600))}
    print(f"\n== Part C [{label}] held-out skies 100-199: {len(rec)} pulsars ({len(rec)/100:.1f}/sky);"
          f" n(Kc>=60)={int((Kc>=60).sum())}, n(Kc>=600)={int((Kc>=600).sum())}, n(Kf>=60)={int((Kf>=60).sum())}, n(Kf>=600)={int((Kf>=600).sum())}")
    print(f"   prior-only (decode = prior centre): K>=60 {prior['K60']:.3f}, K>=600 {prior['K600']:.3f}, all {prior['all']:.3f};"
          f" finest-K: {prior['K60f']:.3f}, {prior['K600f']:.3f}")
    rows = {}
    for a in arms:
        r_ = dict(n_med=float(np.median([r["n|"+a] for r in rec])), info_med=float(np.median([r["i|"+a] for r in rec])))
        for g in ("f", "b"):
            for n, m in (("K60", Kc >= 60), ("K600", Kc >= 600), ("all", allm), ("K60f", Kf >= 60), ("K600f", Kf >= 600)):
                r_[f"{g}_{n}"] = frac(rec, f"w{g}|{a}", m)
        rows[a] = r_
        print(f"  {a:>20}: FINE K>=60 {f3(r_['f_K60'])}  K>=600 {f3(r_['f_K600'])}  all {f3(r_['f_all'])}"
              f" | finest-K {r_['f_K60f'][0]:.3f}/{r_['f_K600f'][0]:.3f} | BOOKED-GRID {r_['b_K60'][0]:.3f}/{r_['b_K600'][0]:.3f}/{r_['b_all'][0]:.3f}"
              f" | n_src {r_['n_med']:.0f} info {r_['info_med']:.2f}")
    return dict(rows=rows, prior=prior, n=len(rec)), rec

OUT["C_harness"], recH = tableC("harness", A12, "harness priors")
OUT["C_census"], recC = tableC("census", A12, "census priors")

# grid-step facts
dl = np.array([r["dL_off"] for r in recH]); pri = np.array([r["pri"] for r in recH])
Kc = np.array([r["Kc"] for r in recH])
print(f"\n  booked grid step = 6 sigma_d/600: median over K>=60 pulsars {np.median(6*pri[Kc>=60]/600):.3f} pc,"
      f" over K>=600 {np.median(6*pri[Kc>=600]/600):.2f} pc (window 0.25 pc); fine grid step max {dl.max():.3f} pc;"
      f" fine grid relaxed (FFT cap) for {np.mean(dl > 0.0501):.4f} of pulsars")
OUT["grid"] = dict(step60=float(np.median(6*pri[Kc >= 60]/600)), step600=float(np.median(6*pri[Kc >= 600]/600)))

# ---------------- C2 anchors -------------------------------------------------------------
inf = ld("census")["info"]
a12 = np.concatenate([d["anchors"][d["is12"]] for d in inf]); aall = np.concatenate([d["anchors"] for d in inf])
infH = [d for d in ld("harness")["info"] if d["sky"] >= 100]
h12 = np.concatenate([d["anchors"][d["is12"]] for d in infH]); hall = np.concatenate([d["anchors"] for d in infH])
print(f"\n== C2 anchors: pulsars (of 79 red-noise-clean per sky) with K < 1 for a given source")
for lab, v in (("census priors, booked top-12 sources", a12), ("census priors, every candidate", aall),
               ("harness priors, booked top-12 sources", h12), ("harness priors, every candidate", hall)):
    print(f"   {lab:>40}: per (sky, source) median {np.median(v):.0f}, 10/90 {np.quantile(v,[.1,.9]).tolist()}, share with >= 3: {np.mean(v >= 3):.3f}")
OUT["C2"] = dict(census12=[float(np.median(a12)), float(np.mean(a12 >= 3))], census_all=[float(np.median(aall)), float(np.mean(aall >= 3))],
                 harness12=[float(np.median(h12)), float(np.mean(h12 >= 3))], harness_all=[float(np.median(hall)), float(np.mean(hall >= 3))])

# ---------------- Part D sensitivity ---------------------------------------------------
print("\n== Part D: timing noise / n (strains untouched), arm p2r every candidate, held-out skies")
OUT["D"] = {}
for n, run in ((1, "harness"), (3, "sens3"), (10, "sens10")):
    rec = [r for r in ld(run)["rec"] if r["sky"] >= 100]
    Kc_ = np.array([r["Kc"] for r in rec]); al = np.ones(len(rec), bool)
    for t in ("oracle", "machine"):
        a = f"on|all|{t}"
        v = dict(ncand=float(np.mean([r["ncand"] for r in rec])), info=float(np.median([r["i|"+a] for r in rec])),
                 K60=frac(rec, "wf|"+a, Kc_ >= 60), K600=frac(rec, "wf|"+a, Kc_ >= 600), all=frac(rec, "wf|"+a, al))
        OUT["D"][f"{n}|{t}"] = v
        print(f"  n={n:2d} {t:>8}: candidates/pulsar {v['ncand']:.0f} | K>=60 {f3(v['K60'])} | K>=600 {f3(v['K600'])} | all {f3(v['all'])} | info {v['info']:.2f}")

# ---------------- Part E gate ------------------------------------------------------------
def gateE(sk, wr, al, ntr=100, Kc=None):
    cal = sk < ntr; tst = sk >= ntr; u = np.unique(sk[tst]); nte = len(u)
    qs = np.linspace(0.5, 100, 200); k = np.ones(7)
    th = np.array([np.quantile(al[cal], q/100) for q in qs])
    pu = np.array([1-wr[cal & (al <= t)].mean() if (cal & (al <= t)).sum() >= 10 else np.nan for t in th])
    v = np.isfinite(pu)
    ps = np.convolve(np.where(v, pu, 0), k, "same")/np.maximum(np.convolve(v*1.0, k, "same"), 1)
    cr = [a for a in range(len(ps)-1) if ps[a] >= 0.90 >= ps[a+1]]
    if ps[-1] >= 0.90:
        thr, gate = th[-1], 100.0
    elif cr:
        a = cr[-1]; f = (ps[a]-0.90)/max(ps[a]-ps[a+1], 1e-9)
        thr = th[a]+f*(th[a+1]-th[a]); gate = qs[a]+f*(qs[a+1]-qs[a])
    else:
        return dict(ok=False, maxpur=float(np.nanmax(ps)))
    s = tst & (al <= thr)
    per = np.array([[(s & (sk == j) & ~wr).sum(), (s & (sk == j)).sum()] for j in u])
    rg = np.random.default_rng(23); bs = rg.integers(0, nte, (1000, nte)); t = per[bs].sum(1)
    out = dict(ok=True, gate=float(gate), pur=float(1-wr[s].mean()), y=float((~wr[s]).sum()/nte),
               yb=np.quantile(t[:, 0]/nte, [.05, .95]).tolist(), pb=np.quantile(t[:, 0]/np.maximum(t[:, 1], 1), [.05, .95]).tolist(),
               n=int(s.sum()))
    if Kc is not None:
        out["fracK3"] = float((Kc[s] < 3).mean()); out["fracK1"] = float((Kc[s] < 1).mean())
        out["nK10_right"] = float((s & (Kc >= 10) & ~wr).sum()/nte)
    return out

def fmtE(g):
    if not g["ok"]: return f"no 0.90 catalogue (max train purity {g['maxpur']:.3f})"
    s = f"{g['y']:.2f} [{g['yb'][0]:.2f}, {g['yb'][1]:.2f}] @ {g['pur']:.4f} [{g['pb'][0]:.3f}, {g['pb'][1]:.3f}], gate {g['gate']:.1f}%"
    if "fracK3" in g: s += f"; catalogue K<1 {g['fracK1']:.2f}, K<3 {g['fracK3']:.2f}; right with K>=10: {g['nK10_right']:.2f}/sky"
    return s

print("\n== Part E: gate = edge-safe 7-pt mean, scan 0.5-100 %, LAST train crossing of 0.90, blind on held-out skies")
H = ld("harness"); BH = [r for r in H["rec"] if r["booked"]]
assert len(BH) == len(R49["rec"])
sk = np.array([r["sky"] for r in BH]); al = R49["alias"]; KcB = np.array([r["Kc"] for r in BH])
assert np.allclose([r["alias"] for r in BH], al)
E = {}
E["incumbent"] = gateE(sk, np.array([r["wrong_inc"] for r in R49["rec"]], bool), al, Kc=KcB)
E["p2_mf (booked grid)"] = gateE(sk, np.array([r["wrong_p2_mf"] for r in R49["rec"]], bool), al, Kc=KcB)
for lab, key in (("p2_mf (fine grid)", "wf|off|top12|machine"),
                 ("p2r_mf top 12, machine template (booked grid)", "wb|on|top12|machine"),
                 ("p2r_mf top 12, machine template (fine grid)", "wf|on|top12|machine"),
                 ("p2r_mf every candidate, machine template (booked grid)", "wb|on|all|machine"),
                 ("p2r_mf every candidate, machine template (fine grid)", "wf|on|all|machine")):
    E[lab] = gateE(sk, np.array([r[key] for r in BH], bool), al, Kc=KcB)
for lab, g in E.items():
    ung = None
    print(f"  {lab:>56}: {fmtE(g)}")
OUT["E"] = E
print("  amplitude table, arm p2r_mf every candidate, machine template:")
EA = {}
for rate, A in ((0.33, 3.7), (0.60, 5.0), (1.0, 6.4), (1.7, 8.3), (2.7, 10.6)):
    Rr = ld(f"rate{rate:.2f}"); Bq = [r for r in Rr["rec"] if r["booked"]]
    sk_ = np.array([r["sky"] for r in Bq]); al_ = np.array([r["alias"] for r in Bq]); Kq = np.array([r["Kc"] for r in Bq])
    for g_, key in (("fine", "wf|on|all|machine"), ("booked grid", "wb|on|all|machine"), ("top12 fine", "wf|on|top12|machine")):
        w_ = np.array([r[key] for r in Bq], bool)
        ge = gateE(sk_, w_, al_, Kc=Kq); EA[f"{A}|{g_}"] = ge
        print(f"    A={A:4.1f}e-15 [{g_:>11}] ungated {1-w_.mean():.3f}: {fmtE(ge)}")
OUT["E_amp"] = EA

# ---------------- figure: fraction correct vs K -----------------------------------------
def kbin(rec, key, edges):
    Kc_ = np.array([r["Kc"] for r in rec]); ok = ~np.array([r[key] for r in rec], bool)
    return [ok[(Kc_ >= a) & (Kc_ < b)].mean() if ((Kc_ >= a) & (Kc_ < b)).sum() >= 15 else np.nan for a, b in zip(edges[:-1], edges[1:])]
edges = np.array([0.01, 1, 3, 10, 30, 100, 300, 1000, 1e5]); xc = np.array([0.5, 1.7, 5.5, 17, 55, 170, 550, 3000])
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
col = {"oracle": "C0", "machine": "C1", "full": "0.45"}
for a_, rec, ttl in ((ax[0], recH, "harness priors"), (ax[1], recC, "census priors")):
    for (r_, s_, t_) in C.ARMS12:
        a = C.ARMNAME((r_, s_, t_))
        a_.plot(xc, kbin(rec, "wf|"+a, edges), color=col[t_], ls="-" if r_ == "on" else "--",
                marker="o" if s_ == "all" else "x", ms=5, lw=1.3)
    Lt = np.array([r["Ltrue"] for r in rec]); Kc_ = np.array([r["Kc"] for r in rec])
    a_.plot(xc, [(np.abs(Lt[(Kc_ >= a) & (Kc_ < b)]) <= TW).mean() if ((Kc_ >= a) & (Kc_ < b)).sum() >= 15 else np.nan
                 for a, b in zip(edges[:-1], edges[1:])], color="k", ls=":", lw=1.2)
    a_.set_xscale("log"); a_.set_xlabel("K = 6 sigma_d / coarsest tooth spacing"); a_.set_title(ttl); a_.grid(alpha=.25, lw=.5)
ax[0].set_ylabel("fraction decoded correctly (fine grid, ungated)")
from matplotlib.lines import Line2D
ax[1].legend(handles=[Line2D([], [], color="C0", label="oracle templates"), Line2D([], [], color="C1", label="machine template term"),
                      Line2D([], [], color="0.45", label="full template term (dk * L)"),
                      Line2D([], [], color="k", ls="-", label="retardation on"), Line2D([], [], color="k", ls="--", label="retardation off"),
                      Line2D([], [], color="k", marker="o", ls="", label="every candidate"), Line2D([], [], color="k", marker="x", ls="", label="top 12"),
                      Line2D([], [], color="k", ls=":", label="prior only")], fontsize=8, loc="upper right")
plt.tight_layout(); plt.savefig(f"job50_fraction_vs_K{SUF}.png", dpi=140); plt.close()
json.dump(OUT, open(f"job50_data/analysis{SUF}.json", "w"), indent=1, default=float)
print(f"\nwrote job50_fraction_vs_K{SUF}.png, job50_data/analysis{SUF}.json")
