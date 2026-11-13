# ================= JOB 54: prior-width slice with random pulsar-term phases for every candidate ===================
# Job 53 Part 3b/3c rerun (run "sigr": same skies, pulsars, known sets and draws; one change: job54_common.phasors B["phi"]).
import os, sys, csv, json, pickle, shutil
import numpy as np
from scipy.special import erf
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples"); sys.path.insert(0, ".")
from job53_lib import frac, classes, TW
N4 = (7, 12, 30, "all"); SIGS = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0); OUT = {}
new = pickle.load(open("job54_data/run_sigr.pkl", "rb"))["rec"]; old = pickle.load(open("job53_data/run_sig_nat.pkl", "rb"))["rec"]; o53 = pickle.load(open("job53_data/run_main_nat.pkl", "rb"))["rec"]
_have = {r["sky"] for r in new}; old = [r for r in old if r["sky"] in _have]; o53 = [r for r in o53 if r["sky"] in _have]; assert len(old) == len(new) == len(o53)
print(f"corrected run: {len(_have)} natural populations (skies {min(_have)}-{max(_have)}), {len(new)} pulsars; Job 53 numbers below are for the same pulsars")
sky = np.array([r["sky"] for r in new]); one = np.ones(len(new), bool); cl = classes(new)
fr = lambda rec, k, m=None: frac(~np.array([r[k] for r in rec]), one if m is None else m, sky)
print("== census priors: random phases against Job 53 (all pulsars / K>=60 / K>=600)")
for v in ("ref", "A"):
    for n in N4:
        k = f"w|{v}|{n}|t1"; a = [fr(new, k, m) for m in cl.values()]; b = [fr(o53, k, m) for m in cl.values()]
        print(f"   {v:>3s} N={str(n):>3s}: random phases {a[0][0]:.3f} [{a[0][1]:.3f}, {a[0][2]:.3f}] / {a[1][0]:.3f} / {a[2][0]:.3f}   Job 53 {b[0][0]:.3f} / {b[1][0]:.3f} / {b[2][0]:.3f}   pulsars that changed {sum(x[k] != y[k] for x, y in zip(new, o53))}")
        OUT[f"census|{v}|{n}"] = dict(new=[a[q][:3] for q in range(3)], old=[b[q][0] for q in range(3)])
if not os.path.exists("talk_data/prior_slice_job53.csv"): shutil.copy("talk_data/prior_slice.csv", "talk_data/prior_slice_job53.csv")
pri = np.array([r["pri"] for r in new]); Ltu = np.array([r["Ltrue"] for r in new])/pri; rows = []
def crossx(y, level):
    x = np.array(SIGS)
    for a in range(len(x)-1):
        if y[a] >= level >= y[a+1]:
            f = (y[a]-level)/max(y[a]-y[a+1], 1e-12); return float(10**(np.log10(x[a])+f*(np.log10(x[a+1])-np.log10(x[a]))))
    return np.nan
print("\n== prior-width slice (fraction right; every pulsar the same prior width)")
T = {}
for v in ("ref", "A"):
    for n in N4:
        F = np.array([fr(new, f"w|{v}|{n}|t1|s{s:g}")[:3] for s in SIGS]); O = np.array([fr(old, f"w|{v}|{n}|t1|s{s:g}")[0] for s in SIGS]); T[(v, n)] = F
        print(f"   {v:>3s} N={str(n):>3s} random phases: " + "  ".join(f"{s:g}: {f[0]:.3f}" for s, f in zip(SIGS, F)) + f" | crosses 0.9 / 0.5 / 0.1 at {crossx(F[:, 0], 0.9):.2f} / {crossx(F[:, 0], 0.5):.2f} / {crossx(F[:, 0], 0.1):.1f} pc")
        print(f"                Job 53:        " + "  ".join(f"{s:g}: {f:.3f}" for s, f in zip(SIGS, O)))
        OUT[f"sig|{v}|{n}"] = dict(new=F.tolist(), old=O.tolist(), cross=[crossx(F[:, 0], x) for x in (0.9, 0.5, 0.1)])
for qi, s in enumerate(SIGS):
    row = [s, f"{np.mean(np.abs(Ltu*s) <= TW):.4f}", f"{erf(0.25/(np.sqrt(2)*s)):.4f}"]
    for v in ("ref", "A"):
        for n in N4: row += [f"{x:.4f}" for x in T[(v, n)][qi]]
    rows.append(row)
with open("talk_data/prior_slice.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["prior_width_pc", "prior_only_measured", "prior_only_analytic"] + [f"{'reference' if v == 'ref' else 'earth_terms_removed'}_N{n if n != 'all' else '_every_candidate'}{sfx}" for v in ("ref", "A") for n in N4 for sfx in ("", "_lo95", "_hi95")]); w.writerows(rows)
print("wrote talk_data/prior_slice.csv (corrected); Job 53 version kept as talk_data/prior_slice_job53.csv")
# per-bin model (job54_parts34.py, Part 3) against the corrected skies
p3 = json.load(open("job54_data/parts34.json"))["part3"]
print("\n== Part 3 again: per-bin model with real shares against the corrected sky result (N = 12, reference, all pulsars)")
for tag in ["census"]+[f"{s:g}" for s in SIGS]:
    sk = fr(new, "w|ref|12|t1")[0] if tag == "census" else fr(new, f"w|ref|12|t1|s{tag}")[0]
    print(f"   {tag:>7s}: sky, random phases {sk:.3f} | per-bin model {p3[tag]['model']:.3f} ({p3[tag]['model']-sk:+.3f}) | shared bins {p3[tag]['model_shared']:.3f} ({p3[tag]['model_shared']-sk:+.3f}) | Job 53 sky {p3[tag]['sky']:.3f}")
    OUT[f"p3|{tag}"] = dict(sky=float(sk), model=p3[tag]["model"], shared=p3[tag]["model_shared"], sky53=p3[tag]["sky"])
# gate (Job 52 Part 3, unchanged)
print("\n== certification gate (train even skies, held out odd), corrected")
def smooth7(y): return np.array([y[max(i-3, 0):i+4].mean() for i in range(len(y))])
G = np.arange(0.5, 100.01, 0.5)
def gate(rec, mk, wk, lab):
    tr = sky % 2 == 0; nho = len(np.unique(sky[~tr])); m = np.array([r[mk] for r in rec]); ok = ~np.array([r[wk] for r in rec])
    o = np.argsort(-m[tr]); okt = ok[tr][o]; mt = m[tr][o]; ps = smooth7(np.array([okt[:max(int(round(g/100*len(okt))), 1)].mean() for g in G]))
    if not (ps >= 0.90).any(): print(f"   {lab}: no 0.90 catalogue (best {ps.max():.3f})"); OUT["gate|"+lab] = None; return
    gi = np.where(ps >= 0.90)[0][-1]; thr = mt[max(int(round(G[gi]/100*len(okt))), 1)-1]; sel = (~tr) & (m >= thr); pu = frac(ok, sel, sky)
    print(f"   {lab}: {sel.sum()/nho:.2f} pulsars per sky at purity {pu[0]:.3f} [{pu[1]:.3f}, {pu[2]:.3f}] (ungated {ok[~tr].sum()/nho:.1f} right per sky)"); OUT["gate|"+lab] = dict(per_sky=float(sel.sum()/nho), purity=pu[:3])
for s in (1, 3):
    for v in ("ref", "A"):
        for n in (7, 12, "all"): gate(new, f"m|{v}|{n}|t1|s{s}", f"w|{v}|{n}|t1|s{s}", f"sigma_d = {s} pc, {v}, N={n}")
for v in ("ref", "A"):
    for n in (7, 12, "all"): gate(new, f"m|{v}|{n}|t1", f"w|{v}|{n}|t1", f"census, {v}, N={n}")
json.dump(OUT, open("job54_data/tables2.json", "w"), indent=1)
