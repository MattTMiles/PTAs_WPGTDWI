# ================= JOB 51 -- analysis: Part 2, 3, 4 tables ================================
import os, sys, json, pickle
import numpy as np
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
OUT = {}
ld = lambda n: pickle.load(open(f"job51_data/run_{n}.pkl", "rb"))
ld50 = lambda n: pickle.load(open(f"job50_data/run_{n}.pkl", "rb"))
have = lambda n: os.path.exists(f"job51_data/run_{n}.pkl")

def frac(rec, key, mask=None, nb=1000, seed=3):
    """fraction right (per pulsar), 95% interval by resampling skies"""
    sk = np.array([r["sky"] for r in rec]); w = np.array([bool(r[key]) for r in rec])
    m = np.ones(len(rec), bool) if mask is None else mask
    if m.sum() == 0: return (np.nan, np.nan, np.nan, 0)
    us = np.unique(sk); rg = np.random.default_rng(seed)
    num = np.array([np.sum(~w[m & (sk == s)]) for s in us]); den = np.array([np.sum(m & (sk == s)) for s in us])
    bs = []
    for _ in range(nb):
        ii = rg.integers(0, len(us), len(us)); bs.append(num[ii].sum()/max(den[ii].sum(), 1))
    return (float(num.sum()/den.sum()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), int(m.sum()))
f3 = lambda t: f"{t[0]:.3f} [{t[1]:.3f}, {t[2]:.3f}]"

def three(rec, key):
    K = np.array([r["Kc"] for r in rec])
    return frac(rec, key, K >= 60), frac(rec, key, K >= 600), frac(rec, key)

# ---------------- rule 3: Job 50 reference values from its saved run -------------------------
print("== Rule 3: Job 50 references (job50_data/run_harness.pkl, run_census.pkl; collect50, uncommitted 2026-11-11)")
for nm, sel in (("harness", lambda r: r["sky"] >= 100), ("census", lambda r: r["sky"] >= 100)):
    rec = [r for r in ld50(nm)["rec"] if sel(r)]
    a, b, c = three(rec, "wf|on|all|oracle")
    print(f"   {nm:8s} p2r every candidate oracle fine grid, {len(rec)} pulsars: K>=60 {a[0]:.3f}  K>=600 {b[0]:.3f}  all {c[0]:.3f}")
    OUT[f"ref_{nm}"] = [a, b, c]
    if nm == "census":
        rec20 = [r for r in rec if r["sky"] < 120]; a, b, c = three(rec20, "wf|on|all|oracle")
        print(f"   census, skies 100-119 only ({len(rec20)} pulsars): K>=60 {a[0]:.3f}  K>=600 {b[0]:.3f}  all {c[0]:.3f}")
        OUT["ref_census20"] = [a, b, c]

# ---------------- Part 2 ----------------------------------------------------------------------
print("\n== Part 2: joint oracle fit, retardation on, fixed amplitude (rate 1), redshifted masses")
for pr, lab in (("h", "harness priors, skies 100-199"), ("c", "census priors, skies 100-119")):
    for s in ("1", "3", "10"):
        n = f"joint_{pr}_s{s}"
        if not have(n): print(f"   {n}: not run"); continue
        rec = ld(n)["rec"]; K = np.array([r["Kc"] for r in rec])
        row = {}
        for key, nm in (("wf|on|all|oracle", "stack, own noise per candidate (Job 50 method)"),
                        ("wsh", "stack on the shared bins"), ("wj", "JOINT FIT on the shared bins")):
            a, b, c = three(rec, key); row[key] = [a, b, c]
            print(f"   {lab} | timing x{s:>3s} | {nm:46s} | K>=60 {f3(a)} (n={a[3]}) | K>=600 {f3(b)} (n={b[3]}) | all {f3(c)}")
        g = lambda k: float(np.median([r[k] for r in rec]))
        row["share"] = dict(median_bins=g("j_share"), power_weighted=g("j_share_sum"), bins_le_10nHz=g("j_share_lo"),
                            pulsar_terms_only=g("j_shareP"), timing_over_unmodelled=g("j_tim_over_un"),
                            ncand=g("ncand"), cap=int(sum(r["j_cap"] for r in rec)), dLmax=float(max(r["j_dL"] for r in rec)),
                            npsr=len(rec))
        sh = row["share"]
        print(f"      candidates' share of each bin's GW power, median over bins (2-30 nHz) and pulsars: {sh['median_bins']:.3f}; "
              f"bins <= 10 nHz: {sh['bins_le_10nHz']:.3f}; summed over the band: {sh['power_weighted']:.3f}; "
              f"candidates/pulsar {sh['ncand']:.0f}; timing share of unmodelled power (median bin) {sh['timing_over_unmodelled']:.3f}; "
              f"grid cap binds for {sh['cap']} of {sh['npsr']} pulsars (max step {sh['dLmax']:.3f} pc)")
        OUT[n] = row

# ---------------- Part 3 ----------------------------------------------------------------------
if have("stat_h_s1"):
    print("\n== Part 3: array SNR of each candidate's Earth term (x1 candidate list), harness priors, marginalised amplitude, skies 100-199")
    info = ld("stat_h_s1")["info"]
    LV = (("x1", 1.0), ("x10", 10.0), ("x100", 100.0), ("zero", np.inf))
    P3 = {}
    for lv, s in LV:
        acc = {v: {"a12": [], "aall": [], "d12": [], "dall": []} for v in ("noise", "bg", "opt")}
        need = {"12": [], "all": []}
        for d in info:
            is12 = d["is12"]; Nt = d["Nt"]/s**2; Ps = d["Psig"]
            arrs = {"noise": np.sqrt(24*2*Ps/Nt) if np.isfinite(s) else np.full(len(Ps), np.inf),
                    "bg": np.sqrt(24*np.mean(2*Ps[None, :]/(Nt[None, :]+d["oth_bg"].astype(float)), axis=0)),
                    "opt": np.sqrt(24*np.mean(2*Ps[None, :]/(Nt[None, :]+d["oth_opt"].astype(float)), axis=0))}
            for v, a in arrs.items():
                acc[v]["a12"].append(a[is12]); acc[v]["aall"].append(a)
                dl = d["cpair"]/a[None, :]
                acc[v]["d12"].append(dl[:, is12].ravel().astype(np.float32)); acc[v]["dall"].append(dl.ravel().astype(np.float32))
            need["12"].append(d["cpair"][:, is12].ravel()); need["all"].append(d["cpair"].ravel())
        P3[lv] = {}
        for v in ("noise", "bg", "opt"):
            for grp, ka, kd in (("booked top 12", "a12", "d12"), ("every candidate", "aall", "dall")):
                a = np.concatenate(acc[v][ka]); dl = np.concatenate(acc[v][kd])
                q = np.percentile(dl, [10, 50, 90])
                r = dict(arr_med=float(np.median(a)), arr_p90=float(np.percentile(a, 90)), dkL=[float(x) for x in q],
                         lt1=float(np.mean(dl < 1)), ltpi=float(np.mean(dl < np.pi)), npairs=int(len(dl)), nsrc=int(len(a)))
                P3[lv][f"{v}|{grp}"] = r
                print(f"   timing {lv:>4s} | arr_{v:5s} | {grp:15s} | array SNR median {r['arr_med']:9.3g}, 90th {r['arr_p90']:9.3g} | "
                      f"dkappa*L 10/50/90%: {q[0]:8.3g} / {q[1]:8.3g} / {q[2]:8.3g} rad | pairs < 1 rad {r['lt1']:.4f}, < pi {r['ltpi']:.4f}")
    for grp in ("12", "all"):
        c = np.concatenate(need[grp])
        P3[f"arr_needed_median_pair_{grp}"] = float(np.median(c))
        print(f"   array SNR needed for the median (pulsar, source) pair to reach 1 rad, {'booked top 12' if grp=='12' else 'every candidate'}: {np.median(c):.0f}"
              f"  ({len(c)} pairs: {len(c)/len(info):.0f} per sky)")
    OUT["part3_stats"] = P3
    print("\n   Part D row rerun: every candidate, FULL template term (dkappa*L_true), retardation on, fine grid, marginalised amplitude")
    r50 = [r for r in ld50("harness")["rec"] if r["sky"] >= 100]; a, b, c = three(r50, "wf|on|all|full")
    print(f"   Job 50, arr as coded, timing x1                : K>=60 {f3(a)} | K>=600 {f3(b)} | all {f3(c)}")
    prc = float(np.mean([abs(r['Ltrue']) <= 0.25 for r in r50]))
    print(f"   prior only (decode = prior centre)             : all {prc:.3f}")
    OUT["part3_dec"] = {"job50_noise_x1": [a, b, c], "prior_only_all": prc}
    for s in ("1", "10", "inf"):
        for v in ("bg", "opt"):
            n = f"dec_{v}_s{s}"
            if not have(n): print(f"   {n}: not run"); continue
            rec = ld(n)["rec"]; a, b, c = three(rec, "wf|on|all|full")
            print(f"   dkappa from arr_{v:3s}, timing {'zero' if s=='inf' else 'x'+s:>4s}                  : K>=60 {f3(a)} | K>=600 {f3(b)} | all {f3(c)} | candidates/pulsar {np.median([r['ncand'] for r in rec]):.0f}")
            OUT["part3_dec"][n] = [a, b, c]

# ---------------- Part 4 ----------------------------------------------------------------------
THR = (0.5, 0.2, 0.1)
for nm, lab in (("stat_h_s1", "harness priors"), ("stat_c_s1", "census priors")):
    if not have(nm): continue
    info = ld(nm)["info"]; P4 = {}
    print(f"\n== Part 4 [{lab}]: pulsars (of {np.mean([d['nclean'] for d in info]):.1f} red-noise-clean per sky) in which a source's pulsar term holds a share > threshold of its bin; retardation on; {len(info)} skies")
    for sdef, sl in (("gw", "share of GW power (machine's s, no depro factor)"), ("tot", "share of GW + timing power")):
        for kn, kl in (("any", "all pulsars"), ("K1", "pulsars with K < 1 for that source"), ("K10", "pulsars with K < 10 for that source")):
            for grp in ("booked top 12", "every candidate"):
                for ti, t in enumerate(THR):
                    per = []; n4 = []; n10 = []; n30 = []
                    for d in info:
                        c = d[f"n_{sdef}_{kn}"][ti]; c = c[d["is12"]] if grp.startswith("booked") else c
                        per.append(c); n4.append((c >= 4).sum()); n10.append((c >= 10).sum()); n30.append((c >= 30).sum())
                    per = np.concatenate(per)
                    r = dict(med=float(np.median(per)), p90=float(np.percentile(per, 90)), mx=int(per.max()),
                             n4=float(np.mean(n4)), n10=float(np.mean(n10)), n30=float(np.mean(n30)),
                             nsrc=float(len(per)/len(info)))
                    P4[f"{sdef}|{kn}|{grp}|{t}"] = r
                    if sdef == "gw" or kn == "any":
                        print(f"   {sl[:22]:22s} | {kl:34s} | {grp:15s} ({r['nsrc']:6.1f} sources/sky) | share > {t}: pulsars per source median {r['med']:.0f}, 90th {r['p90']:.0f}, max {r['mx']} | "
                              f"sources per sky with >=4: {r['n4']:.2f}, >=10: {r['n10']:.2f}, >=30: {r['n30']:.2f}")
    OUT["part4_"+nm] = P4
json.dump(OUT, open("job51_data/analysis.json", "w"), indent=1, default=float)
print("\nwrote job51_data/analysis.json")
