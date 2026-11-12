# ================= JOB 49 -- tables from saved files: D2, D3, D4, D5, E, Part C in-situ =
# Reads job49_data/*.json|pkl written by job49_D2D3.py.  Every printed number is
# read from those files.  D4 recomputes inventory counts (needs the machine
# namespace).  E fits purity against (N_eff, beta, K) in two forms.
import os, sys, json, pickle, csv, glob
import numpy as np
from scipy.special import i0e, i1e
from scipy.optimize import minimize
from scipy import stats
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job49_common as J
os.chdir(J.BASE)
OUT = {}
ARMS = J.ARMS
LAB = {"inc": "incumbent", "p1_inc": "P1 D-inc", "p1_mf": "P1 D-mf",
       "p2_inc": "P2 D-inc", "p2_mf": "P2 D-mf", "p2g_inc": "P2g D-inc", "p2g_mf": "P2g D-mf"}

# ---------------- D2 / D3 tables -----------------------------------------------------
def fmt(b):
    if b is None: return "--"
    if not b["ok"]: return f"no 0.90 crossing (max train {b['maxpur']:.3f})"
    return (f"{b['y']:.2f} [{b['yb'][0]:.2f}, {b['yb'][1]:.2f}] @ {b['pur']:.4f} "
            f"[{b['pb'][0]:.3f}, {b['pb'][1]:.3f}], gate {b['gate']:.1f}%")
rows = []
for fn in ["job49_data/D2D3_marg.json"] + sorted(glob.glob("job49_data/D2D3_rate*.json")):
    d = json.load(open(fn))
    print(f"\n== {d['cfg']}: repro {'OK' if d['repro_ok'] else 'FAIL'} ({d['repro_txt']}); "
          f"switch identical {d['identical']}; n={d['n']}")
    for arm in ARMS:
        a = d["arms"][arm]
        print(f"  {LAB[arm]:>10}: ungated {a['ungated_purity']:.4f} | re-frozen: {fmt(a['refrozen'])}"
              f" | OLD threshold: {fmt(a['oldgate'])}")
        rows.append([d["cfg"], arm, a["ungated_purity"],
                     *(a["refrozen"].get(k) if a["refrozen"]["ok"] else None for k in ("y",)),
                     *(a["refrozen"]["yb"] if a["refrozen"]["ok"] else [None, None]),
                     a["refrozen"]["pur"] if a["refrozen"]["ok"] else None,
                     *(a["refrozen"]["pb"] if a["refrozen"]["ok"] else [None, None]),
                     a["refrozen"]["gate"] if a["refrozen"]["ok"] else None,
                     a["refrozen"].get("maxpur"),
                     a["oldgate"]["y"] if a["oldgate"] else None,
                     a["oldgate"]["pur"] if a["oldgate"] else None])
with open("job49_D2D3_table.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["config", "arm", "ungated_purity", "yield_refrozen", "y_q10", "y_q90",
                 "purity_refrozen", "p_q10", "p_q90", "gate_pct", "max_train_purity",
                 "yield_oldthr", "purity_oldthr"])
    for r in rows: cw.writerow(["" if v is None else (f"{v:.5g}" if isinstance(v, float) else v) for v in r])
print("wrote job49_D2D3_table.csv")

# ---------------- coherence helpers ---------------------------------------------------
def rice_cos_g(g):
    g = np.maximum(g, 1e-12)
    return np.sqrt(np.pi)/2*np.sqrt(g)*(i0e(g/2)+i1e(g/2))

def comb_terms(r, dataset):
    """per-comb (beta, f_old, f_rice) for a record and dataset.
    dataset 'inc': the generator's implied share beta = 1/(1+sigma^2) (sigma incl.
    meas + fuzzy).  'p1'/'p2': the TRUE share of the datum,
    beta = 1/(1 + m Xc + 2/rho^2), fuzzy as coherence factor e^{-s_fz^2/2}."""
    sg = r["sig"]; sfz = r["sfz"]; Xc = r["Xc"]; rho = r["rho"]
    if dataset == "inc":
        beta = 1/(1+sg**2)
        f_old = np.exp(-sg**2)
        f_rice = rice_cos_g(beta/(1-beta))**2
    else:
        m = 1 if dataset == "p1" else 2
        P = m*Xc + 2/np.maximum(rho, 1e-9)**2
        beta = 1/(1+P)
        f_old = np.exp(-(1-beta)/beta - sfz**2)
        f_rice = rice_cos_g(1/P)**2*np.exp(-sfz**2)
    return beta, f_old, f_rice

# ---------------- D5 ------------------------------------------------------------------
R = pickle.load(open("job49_data/rec_marg.pkl", "rb"))
rec, al = R["rec"], R["alias"]
K = np.array([r["K"] for r in rec]); sk = np.array([r["sky"] for r in rec])
print("\n== D5 (headline config, harness priors; ungated = argmax correct among certifiable records)")
d5 = {}
for lo, lab in ((0, "all"), (60, "K>=60"), (600, "K>=600")):
    m = K >= lo
    d5[lab] = {"n": int(m.sum())}
    s = f"  {lab:>7} n={m.sum():5d}:"
    for arm in ARMS:
        key = "wrong_inc" if arm == "inc" else f"wrong_{arm}"
        w = np.array([r[key] for r in rec], bool)
        # sky-cluster bootstrap of the correct fraction (per-sky counts)
        u = np.unique(sk[m]); per = {j: (int((~w[m & (sk == j)]).sum()), int((m & (sk == j)).sum())) for j in u}
        rg = np.random.default_rng(5); bt = []
        for b in range(1000):
            bs = rg.choice(u, len(u)); n_ = sum(per[j][1] for j in bs)
            if n_: bt.append(sum(per[j][0] for j in bs)/n_)
        fr = float((~w[m]).mean()) if m.any() else np.nan
        ci = np.quantile(bt, [.05, .95]).tolist() if bt else [np.nan, np.nan]
        d5[lab][arm] = dict(frac=fr, ci90=ci)
        s += f" {LAB[arm]} {fr:.3f}[{ci[0]:.3f},{ci[1]:.3f}]"
    print(s)
print(f"  harness prior widths: median {np.median([r['pri'] for r in rec]):.2f} pc; "
      f"K quantiles 50/90/99% = {np.quantile(K, [.5, .9, .99]).round(1).tolist()}")
OUT["D5"] = d5

# ---------------- Part C in-situ: standardised residuals of the booked generator -----
z = np.concatenate([r["eps_res"]/r["sig"] for r in rec])
ks = stats.kstest(z, "norm")
print(f"\n== Part C in-situ (headline run, {len(z)} comb draws): (eps - kappa L)/sigma mean {z.mean():+.4f}"
      f" sd {z.std():.4f} excess kurt {stats.kurtosis(z):+.3f}; KS vs N(0,1) p = {ks.pvalue:.3f}")
OUT["C_insitu"] = dict(n=int(len(z)), mean=float(z.mean()), sd=float(z.std()),
                       kurt=float(stats.kurtosis(z)), ks_p=float(ks.pvalue))
# physical arms' phase-error coherence vs Rice prediction, by comb share
for tag, m in (("p1", 1), ("p2", 2)):
    dph = np.concatenate([r[f"dphi_{tag}"] for r in rec])
    pred = np.concatenate([rice_cos_g(1/(m*r["Xc"]+2/np.maximum(r["rho"], 1e-9)**2))*np.exp(-r["sfz"]**2/2) for r in rec])
    print(f"   {tag}: mean cos(dphi) {np.cos(dph).mean():.4f} vs Rice x fuzzy prediction {pred.mean():.4f}")
    OUT["C_insitu"][tag] = dict(cos_meas=float(np.cos(dph).mean()), cos_pred=float(pred.mean()))
inc_coh = np.concatenate([np.exp(-r["sig"]**2/2) for r in rec])
print(f"   incumbent: mean cos(eps-kL) {np.cos(np.concatenate([r['eps_res'] for r in rec])).mean():.4f}"
      f" vs e^(-sig^2/2) {inc_coh.mean():.4f}")

# ---------------- E: criterion refit ---------------------------------------------------
def fit_form(x, ok):
    """logistic P(correct) = 1/(1+exp(-(a + b ln x))); returns a, b."""
    lx = np.log(np.maximum(x, 1e-9))
    nll = lambda p: -np.sum(ok*(-np.logaddexp(0, -(p[0]+p[1]*lx))) + (~ok)*(-np.logaddexp(0, p[0]+p[1]*lx)))
    r = minimize(nll, [0.0, 1.0], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-9, maxiter=5000))
    return r.x

def collapse(x, ok, ab, nb=12):
    lx = np.log(np.maximum(x, 1e-9))
    e = np.quantile(lx, np.linspace(0, 1, nb+1)); e[-1] += 1e-9
    bi = np.clip(np.digitize(lx, e)-1, 0, nb-1)
    obs = np.array([ok[bi == b].mean() for b in range(nb)])
    fit = np.array([np.mean(1/(1+np.exp(-(ab[0]+ab[1]*lx[bi == b])))) for b in range(nb)])
    r2 = 1-np.sum((obs-fit)**2)/np.sum((obs-obs.mean())**2)
    return float(r2), obs, fit

print("\n== E: purity vs criterion ratio r = sum_k f(beta_k) / (2 ln K), K >= 2 records")
E = {}
for ds, arms in (("inc", ("inc",)), ("p1", ("p1_inc", "p1_mf")), ("p2", ("p2_inc", "p2_mf"))):
    mK = K >= 2
    TB = [comb_terms(r, ds) for r in rec]
    for arm in arms:
        key = "wrong_inc" if arm == "inc" else f"wrong_{arm}"
        ok = ~np.array([r[key] for r in rec], bool)
        for form, idx in (("E-old", 1), ("E-rice", 2)):
            x = np.array([TB[n][idx].sum()/2 for n in range(len(rec))])/np.log(np.maximum(K, 1.0001))
            x, okm = x[mK], ok[mK]
            ab = fit_form(x, okm)
            r2, obs, fit = collapse(x, okm, ab)
            cp = {p: float(np.exp((np.log(p/(1-p))-ab[0])/ab[1])) for p in (0.5, 0.7, 0.8, 0.9, 0.95)}
            # residual vs beta: per-record median comb beta, 6 quantile bins
            bmed = np.array([np.median(TB[n][0]) for n in range(len(rec))])[mK]
            pr = 1/(1+np.exp(-(ab[0]+ab[1]*np.log(np.maximum(x, 1e-9)))))
            e = np.quantile(bmed, np.linspace(0, 1, 7)); e[-1] += 1e-9
            bi = np.clip(np.digitize(bmed, e)-1, 0, 5)
            resid = [(float(np.median(bmed[bi == b])), float((okm[bi == b]-pr[bi == b]).mean()),
                      float(np.sqrt(pr[bi == b].mean()*(1-pr[bi == b].mean())/max((bi == b).sum(), 1))))
                     for b in range(6)]
            ll = float(np.sum(okm*np.log(pr+1e-12)+(~okm)*np.log(1-pr+1e-12)))
            E[f"{arm}|{form}"] = dict(n=int(mK.sum()), a=float(ab[0]), b=float(ab[1]), c=cp, R2=r2,
                                      loglik=ll, resid_vs_beta=resid)
            print(f"  {LAB[arm]:>10} {form:>6}: n={mK.sum()} c(p) " +
                  " ".join(f"{p}:{v:.3g}" for p, v in cp.items()) +
                  f" | R2 {r2:.3f} | lnL {ll:.1f} | resid vs beta: " +
                  " ".join(f"{b_:.2f}:{r_:+.3f}" for b_, r_, s_ in resid))
OUT["E"] = E

# ---------------- DIAGNOSTIC gate (NOT the booked procedure): edge-safe, 0.5-100 % ----
# book()'s 7-point mean uses np.convolve(mode="same"), which ZERO-PADS the last 3
# points of the 0.5-60 % search range; when train purity is still >= 0.90 at 60 %,
# the "crossing" is manufactured at ~58.6 %.  Diagnostic: same smoothing kernel,
# normalised by the number of available points (no zero padding), range to 100 %.
def book_diag(rec, al, key, ntr=100):
    sk = np.array([r["sky"] for r in rec]); wr = np.array([r[key] for r in rec], bool)
    cal = sk < ntr; tst = sk >= ntr; u = np.unique(sk[tst]); nte = len(u)
    qs = np.linspace(0.5, 100, 200); k = np.ones(7)
    th = [np.quantile(al[cal], q/100) for q in qs]
    pu = np.array([1-wr[cal & (al <= t)].mean() if (cal & (al <= t)).sum() >= 10 else np.nan for t in th])
    v = np.isfinite(pu)
    ps = np.convolve(np.where(v, pu, 0), k, "same")/np.maximum(np.convolve(v*1.0, k, "same"), 1)
    thr = None
    for a in range(len(ps)-1):
        if ps[a] >= 0.90 >= ps[a+1]:
            f = (ps[a]-0.90)/max(ps[a]-ps[a+1], 1e-9); thr = th[a]+f*(th[a+1]-th[a])
            gate = qs[a]+f*(qs[a+1]-qs[a]); break
    if thr is None:
        if ps[-1] >= 0.90: thr, gate = th[-1], 100.0       # never drops: everything passes
        else: return dict(ok=False, maxpur=float(np.nanmax(ps)))
    s = tst & (al <= thr)
    per = {j: (int((s & (sk == j) & ~wr).sum()), int((s & (sk == j)).sum())) for j in u}
    rg = np.random.default_rng(23); by = []; bp = []
    for b in range(400):
        bs = rg.choice(u, size=nte, replace=True)
        r_ = sum(per[j][0] for j in bs); n_ = sum(per[j][1] for j in bs)
        by.append(r_/nte); bp.append(r_/max(n_, 1))
    return dict(ok=True, gate=float(gate), pur=float(1-wr[s].mean()), y=float((~wr[s]).sum()/nte),
                yb=np.quantile(by, [.1, .9]).tolist(), pb=np.quantile(bp, [.1, .9]).tolist())

print("\n== DIAGNOSTIC gate (edge-safe smoothing, search 0.5-100 %) -- NOT the booked procedure")
DG = {}
for fn in ["job49_data/rec_marg.pkl"] + sorted(glob.glob("job49_data/rec_rate*.pkl")):
    tag = os.path.basename(fn)[4:-4]
    Rr = pickle.load(open(fn, "rb")); rr, aa = Rr["rec"], Rr["alias"]
    DG[tag] = {}
    for arm in ARMS:
        key = "wrong_inc" if arm == "inc" else f"wrong_{arm}"
        b = book_diag(rr, aa, key); DG[tag][arm] = b
        print(f"  {tag:>9} {LAB[arm]:>10}: " + (f"{b['y']:.2f} [{b['yb'][0]:.2f},{b['yb'][1]:.2f}] @ {b['pur']:.4f}"
              f" [{b['pb'][0]:.3f},{b['pb'][1]:.3f}] gate {b['gate']:.1f}%" if b["ok"] else
              f"no 0.90 crossing (max train {b['maxpur']:.3f})"))
OUT["diag_gate"] = DG

# ---------------- D4: usable sources per sky --------------------------------------
NS = J.load(verbose=False)
from scipy.optimize import brentq
# physical cut: the SAME minimum coherence the old cut admitted.  Old cut sigma < pi
# admits combs whose assumed coherence e^{-sigma^2/2} >= e^{-pi^2/2} = 0.0072.  Under
# the physical (Rice) law with the two-term ledger, coherence = rice(gamma), gamma =
# 1/(2 X + 2/rho^2), X = res/Pt (pre-depro, like the cut itself).
c_floor = np.exp(-np.pi**2/2)
g_floor = brentq(lambda g: rice_cos_g(g)-c_floor, 1e-12, 1.0)
print(f"\n== D4: physical cut = Rice coherence >= e^(-pi^2/2) = {c_floor:.4f} <=> gamma >= {g_floor:.3g}")
T_OBS = NS["T_OBS"]; YR_S = NS["YR_S"]
cnt = {k: [] for k in ("cand", "old_cut", "old_nmax", "phys1_cut", "phys2_cut", "phys2_nmax")}
rgA = np.random.default_rng(555)
for j in range(200):
    A = 6.4e-15*np.exp(rgA.normal(0, 0.32)); rt = float(np.clip((A/6.4e-15)**2, 0.2, 3.5))
    sd, f2, h2, P0, sea = NS["_pop18"](j % 50, 0.35)
    T = 45.0*YR_S; nt = 45.0*365.25/14.0
    rp = (h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100.0*1e-9); arr = rp*np.sqrt(0.2*120)
    cand = np.where((arr >= 3.0) & (f2 >= 2e-9) & (f2 <= 3e-8) & (rp >= 3.0))[0]   # inv_rate l.493
    if len(cand) > 1300: cand = cand[np.argsort(rp[cand])[::-1][:1300]]           # inv_rate l.494
    o = np.argsort(f2); fs = f2[o]; Pb = P0[o]; half = 100.0/T
    be, bc = NS["bin_edges"], NS["bin_cent"]
    X = np.empty(len(cand))
    for a, k in enumerate(cand):                                                   # inv_rate l.498-505
        lo = np.searchsorted(fs, f2[k]-half); hi = np.searchsorted(fs, f2[k]+half)
        c = np.sinc((fs[lo:hi]-f2[k])*T)**2*Pb[lo:hi]*(T/T_OBS)
        b = int(np.clip(np.digitize(f2[k], be)-1, 0, len(bc)-1))
        X[a] = max((c.sum()-P0[k]*(T/T_OBS)+sea[b]*(T_OBS/T))*rt, 0)/(P0[k]*(T/T_OBS))
    sg = np.sqrt(X+1/rp[cand]**2)
    g1 = 1/(X+2/rp[cand]**2); g2 = 1/(2*X+2/rp[cand]**2)
    cnt["cand"].append(len(cand)); cnt["old_cut"].append(int((sg < np.pi).sum()))
    cnt["old_nmax"].append(min(12, int((sg < np.pi).sum())))
    cnt["phys1_cut"].append(int((g1 >= g_floor).sum())); cnt["phys2_cut"].append(int((g2 >= g_floor).sum()))
    cnt["phys2_nmax"].append(min(12, int((g2 >= g_floor).sum())))
for k, v in cnt.items():
    v = np.array(v); print(f"  {k:>11}: mean {v.mean():7.2f} per sky, median {np.median(v):.0f}, 10/90 {np.quantile(v,[.1,.9]).tolist()}")
# per-pulsar usable combs in the decode (N_eff as booked: records' comb counts)
ne = np.array([len(r["sig"]) for r in rec])
b_inc = np.concatenate([1/(1+r["sig"]**2) for r in rec])
b_tr2 = np.concatenate([1/(1+2*r["Xc"]+2/np.maximum(r["rho"], 1e-9)**2) for r in rec])
print(f"  per certifiable pulsar: usable combs (sigma_inc<pi after depro+fuzzy) mean {ne.mean():.2f}, median {np.median(ne):.0f}")
print(f"  comb share: machine s=1/(1+sigma^2) median {np.median(b_inc):.3f}; true two-term share median {np.median(b_tr2):.3f}")
print(f"  sigma_inc median over decode combs {np.median(np.concatenate([r['sig'] for r in rec])):.2f} rad")
OUT["D4"] = dict(c_floor=float(c_floor), g_floor=float(g_floor),
                 counts={k: dict(mean=float(np.mean(v)), median=float(np.median(v))) for k, v in cnt.items()},
                 neff_mean=float(ne.mean()), share_machine_median=float(np.median(b_inc)),
                 share_true2_median=float(np.median(b_tr2)))
json.dump(OUT, open("job49_data/analysis.json", "w"), indent=1, default=float)
print("\nwrote job49_data/analysis.json")
