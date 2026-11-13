# ================= JOB 54 Part 1D-1G: tables behind the talk figures ==============================
# Exports saved Job 53 results (job53_data/run_*.pkl, heat_*.npz) and the Job 54 runs (job54_data/run_*.pkl).
import os, sys, csv, json, pickle
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
from scipy.special import i0e, i1e, erf
from scipy.ndimage import uniform_filter
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples"); sys.path.insert(0, ".")
from job53_lib import frac, popid, classes, TW
from job53_common import NLIST, N4, SIGS
TD = "talk_data/"; os.makedirs(TD, exist_ok=True)
def W(name, header, rows):
    with open(TD+name, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)
    print(f"wrote {name}: {len(rows)} rows x {len(header)} columns")
L53 = lambda n: pickle.load(open(f"job53_data/run_{n}.pkl", "rb")); L54 = lambda n: pickle.load(open(f"job54_data/run_{n}.pkl", "rb"))
f3 = lambda x: f"{x:.4f}"
def fr(rec, key, m=None, pop=None):
    ok = ~np.array([r[key] for r in rec]); pop = np.array([r["sky"] for r in rec]) if pop is None else pop
    return frac(ok, np.ones(len(rec), bool) if m is None else m, pop)

# ---------------- curve.csv (Job 53 saved per-pulsar results) -------------------------------------
rn53 = L53("main_nat")["rec"]; rn = L54("mainr")["rec"]; S53 = L53("stat_nat"); cl = classes(rn); nsky = len(S53["info"])
nbg = float(np.mean([(d["arr_bg"] >= 3).sum() for d in S53["info"]])); nop = float(np.mean([(d["arr_opt"] >= 3).sum() for d in S53["info"]]))
ncm = float(np.mean([r["nC"] for r in rn])); po = [float(np.mean(np.abs(np.array([r["Ltrue"] for r in rn])[m]) <= TW)) for m in cl.values()]
rows = []
for arm, lab in (("ref", "reference"), ("A", "earth_terms_removed")):
    for tl, tlab in (("t1", "x1"), ("t10", "x10"), ("tinf", "zero")):
        for n in NLIST:
            k = f"w|{arm}|{n}|{tl}"
            if k not in rn[0]: continue
            a = [fr(rn, k, m) for m in cl.values()]
            rows.append([lab, tlab, n if n != "all" else "every_candidate", f"{n if n != 'all' else ncm:.0f}" if n == "all" else n] + [f3(v) for t in a for v in t[:3]] + [f"{a[0][0]*79:.3f}", f3(po[0]), f3(po[1]), f3(po[2]), f"{nbg:.3f}", f"{nop:.3f}"])
W("curve.csv", ["arm", "timing", "N_known", "N_known_number", "frac_all", "frac_all_lo95", "frac_all_hi95", "frac_K60", "frac_K60_lo95", "frac_K60_hi95", "frac_K600", "frac_K600_lo95", "frac_K600_hi95",
                "pulsars_right_per_sky_of_79", "prior_only_all", "prior_only_K60", "prior_only_K600", "detectable_per_sky_arr_bg_ge3", "detectable_per_sky_arr_bg_opt_ge3"], rows)

# ---------------- prior_slice.csv -------------------------------------------------------------------
SG = L53("sig_nat")["rec"]; pri = np.array([r["pri"] for r in SG]); Ltu = np.array([r["Ltrue"] for r in SG])/pri; rows = []
for s in SIGS:
    row = [s, f3(float(np.mean(np.abs(Ltu*s) <= TW))), f3(float(erf(0.25/(np.sqrt(2)*s))))]
    for arm in ("ref", "A"):
        for n in N4:
            a = fr(SG, f"w|{arm}|{n}|t1|s{s:g}"); row += [f3(a[0]), f3(a[1]), f3(a[2])]
    rows.append(row)
W("prior_slice_job53.csv", ["prior_width_pc", "prior_only_measured", "prior_only_analytic"] + [f"{'reference' if arm == 'ref' else 'earth_terms_removed'}_N{n if n != 'all' else '_every_candidate'}{sfx}" for arm in ("ref", "A") for n in N4 for sfx in ("", "_lo95", "_hi95")], rows)

# ---------------- heatmaps + contours ------------------------------------------------------------------
Z = np.load("job53_data/heat_grid.npz"); XS, BS = Z["XS"], Z["BS"]; KK = json.load(open("job53_data/heat_K.json")); crow = []
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
for Ns in (2, 4, 8, 12):
    G = Z[f"G{Ns}"]; SM = uniform_filter(G, size=3, mode="nearest")
    W(f"heatmap_Ns{Ns}.csv", ["prior_width_pc", "beta", "fraction_right", "fraction_right_3x3_mean", "K_coarsest_tooth", "K_finest_tooth"],
      [[f"{x:.5g}", f"{b:.5f}", f3(G[iy, ix]), f3(SM[iy, ix]), f"{x*KK[str(Ns)][0]:.4g}", f"{x*KK[str(Ns)][1]:.4g}"] for iy, b in enumerate(BS) for ix, x in enumerate(XS)])
    fig, ax = plt.subplots(); cs = ax.contour(np.log10(XS), np.log(BS/(1-BS)), SM, levels=[0.5, 0.9])
    for lev, segs in zip(cs.levels, cs.allsegs):
        for si, seg in enumerate(segs):
            crow += [[Ns, lev, si, f"{10**p[0]:.5g}", f"{1/(1+np.exp(-p[1])):.5f}"] for p in seg]
    plt.close(fig)
W("heatmap_contours.csv", ["N_s", "level_fraction_right", "segment", "prior_width_pc", "beta"], crow)

# ---------------- jitter_points.csv, beta_hist.csv (Job 54 run, natural skies, N = 12 pairs) -----------------------
rec = rn
def rice(beta):
    s2 = (1-beta)/beta; k = 1.0/s2; c = np.sqrt(np.pi*k)/2*(i0e(k/2)+i1e(k/2))
    rg = np.random.default_rng(1); z = 1+np.sqrt(s2/2)*(rg.normal(size=200000)+1j*rg.normal(size=200000))
    return float(c), float(np.sqrt(np.mean(np.angle(z)**2)))
be = np.geomspace(0.003, 0.99, 13); rows = []; hb = np.concatenate([[0], np.geomspace(1e-3, 1.0, 61)]); hrows = {}
for st, lab in (("on", "on"), ("off", "off")):
    bu = np.concatenate([r[f"d_bun_{st}"] for r in rec]).astype(float); dp = np.concatenate([r[f"d_dphi_{st}"] for r in rec]).astype(float)
    for a, b in zip(be[:-1], be[1:]):
        m = (bu >= a) & (bu < b)
        if m.sum() < 20: continue
        bm = float(np.median(bu[m])); qs = np.quantile(bu[m], np.linspace(0.05, 0.95, 10)); ex = [rice(x) for x in qs]
        c = np.cos(dp[m]); rows.append([lab, f"{a:.5f}", f"{b:.5f}", f"{bm:.5f}", int(m.sum()), f3(c.mean()), f3(c.std()/np.sqrt(m.sum())), f3(np.sqrt(np.mean(dp[m]**2))),
                                        f3(np.mean([e[0] for e in ex])), f3(np.mean([e[1] for e in ex])), f3(np.sqrt((1-bm)/(2*bm))), f3(np.sqrt((1-bm)/bm))])
    hrows[lab] = np.histogram(bu, bins=hb)[0]; print(f"   retardation {lab}: {len(bu)} pairs; beta_unmod 10/50/90% {np.percentile(bu, [10, 50, 90]).round(4)}; below 0.003: {np.mean(bu < 0.003):.3f}")
W("jitter_points.csv", ["retardation", "beta_lo", "beta_hi", "beta_median", "n_pairs", "cos_dphi_measured", "cos_dphi_se", "rms_dphi_measured_rad", "cos_dphi_exact_rice", "rms_dphi_exact_rice_rad",
                        "rms_small_angle_sqrt((1-b)/(2b))_rad", "rms_old_law_sqrt((1-b)/b)_rad"], rows)
W("beta_hist.csv", ["beta_lo", "beta_hi", "n_pairs_retardation_on", "n_pairs_retardation_off"], [[f"{hb[k]:.5g}", f"{hb[k+1]:.5g}", int(hrows["on"][k]), int(hrows["off"][k])] for k in range(len(hb)-1)])

# ---------------- ladder.csv --------------------------------------------------------------------------------------
sky = np.array([r["sky"] for r in rec]); Kb = np.array([r["Kc"] for r in rec]); rows = []
J53 = {15: {d["sky"] for d in L53("span15_nat")["info"]}, 25: {d["sky"] for d in L53("span25_nat")["info"]}}
SP = {T: L54(f"span{T}r")["rec"] for T in (15, 25)}
def lrow(order, step, label, n, recs, key, skies="all 200", mask=None):
    s_ = np.array([r["sky"] for r in recs]); m0 = np.ones(len(recs), bool) if mask is None else mask
    kb = np.array([r["Kc"] for r in recs]); kk = np.array([r[f"Kco|{n}"] for r in recs])
    a = fr(recs, key, m0, s_); b = fr(recs, key, m0 & (kb >= 60), s_); c = fr(recs, key, m0 & (kk >= 60), s_)
    rows.append([n, order, step, label, skies, len(np.unique(s_[m0])), int(m0.sum())] + [f"{v:.3f}" for v in a[:3]] + [f"{a[0]*80:.2f}"] + [f"{v:.3f}" for v in b[:3]] + [b[3]] + [f"{v:.3f}" for v in c[:3]] + [c[3]])
    return a[0]
for n in (12, 7):
    for order, s2, l2, l3 in (("sea_first", "q2", "add the sea", "add every other binary (no retardation anywhere)"), ("binaries_first", "q2s", "add every other binary (no sea, no retardation)", "add the sea")):
        lrow(order, "1", "known binaries + timing noise only (quiet sky, no retardation)", n, rec, f"w|q1|{n}|t1")
        lrow(order, "2", l2, n, rec, f"w|{s2}|{n}|t1"); lrow(order, "3", l3, n, rec, f"w|off|{n}|t1")
        lrow(order, "4", "add retardation (the reference)", n, rec, f"w|ref|{n}|t1")
        lrow(order, "5", "step 4 with timing x10 better", n, rec, f"w|ref|{n}|t10"); lrow(order, "6", "step 4 with timing noise / 1e6", n, rec, f"w|ref|{n}|tinf")
        for T in (25, 15):
            s_ = np.array([r["sky"] for r in SP[T]]); inj = np.array([x in J53[T] for x in s_])
            lrow(order, f"7_{T}yr", f"step 4 with a {T}-year span", n, SP[T], f"w|ref|{n}|t1")
            lrow(order, f"7_{T}yr", f"step 4 with a {T}-year span", n, SP[T], f"w|ref|{n}|t1", skies=f"the {len(J53[T])} skies Job 53 used", mask=inj)
            lrow(order, f"7_{T}yr", f"step 4 with a {T}-year span", n, SP[T], f"w|ref|{n}|t1", skies=f"the {200-len(J53[T])} skies Job 53 skipped", mask=~inj)
W("ladder.csv", ["N_known", "order", "step", "what_is_added", "skies", "n_skies", "n_pulsars", "frac_all", "frac_all_lo95", "frac_all_hi95", "count_of_80_pulsars",
                 "frac_K60_booked_list", "frac_K60_booked_lo95", "frac_K60_booked_hi95", "n_pulsars_K60_booked", "frac_K60_known_set", "frac_K60_known_lo95", "frac_K60_known_hi95", "n_pulsars_K60_known"], rows)
for r in rows:
    if r[1] == "sea_first" or r[2] in ("2", "3"): print("   ", r[0], r[1], r[2], r[4], "| all", r[7], f"[{r[8]}, {r[9]}]", "| K>=60 booked", r[11], "| K>=60 known set", r[15], "|", r[3])
# check against Job 53: skies used -> same numbers
for T, old in ((15, "span15_nat"), (25, "span25_nat")):
    o = {(r["sky"], r["psr"]): r["w|ref|12|t1"] for r in L53(old)["rec"]}; nw = {(r["sky"], r["psr"]): r["w|ref|12|t1"] for r in SP[T]}
    print(f"   T = {T} yr: on the skies Job 53 used, flags equal Job 53 (aligned phases) for {sum(nw[k] == v for k, v in o.items())}/{len(o)} pulsars; fraction right {1-np.mean(list(o.values())):.4f} (Job 53) -> {1-np.mean([nw[k] for k in o]):.4f}")
for n_ in (7, 12, "all"):
    k_ = f"w|ref|{n_}|t1"; print(f"   T = 45 yr reference N = {n_}: flags equal Job 53 (aligned phases) for {sum(a[k_] == b[k_] for a, b in zip(rec, rn53))}/{len(rec)} pulsars; fraction right {1-np.mean([b[k_] for b in rn53]):.4f} (Job 53) -> {1-np.mean([a[k_] for a in rec]):.4f}")

# ---------------- answers_example.csv (natural sky 100, N = 12, reference) -------------------------------------------
r100 = [r for r in rec if r["sky"] == 100]; rows = []
for q, r in enumerate(r100):
    cen = r["L_pc"]-r["Ltrue"]; Lh = r["Lh|ref|12|t1"]
    rows.append([r["psr"], q, f"{r['L_pc']:.4f}", f"{cen:.4f}", f"{r['pri']:.4f}", f"{r['Kc']:.3f}", f"{r['Kco|12']:.3f}", f"{r['Kfin|12']:.3f}", f"{cen+Lh:.4f}", f"{Lh-r['Ltrue']:.4f}",
                 int(not r["w|ref|12|t1"]), f"{r['m|ref|12|t1']:.4f}", int(q % 4 == 0)])
W("answers_example.csv", ["pulsar_index", "row", "true_distance_pc", "prior_centre_pc", "prior_width_pc", "K_booked_coarsest_tooth", "K_known_set_coarsest_tooth", "K_known_set_finest_tooth",
                          "recovered_distance_pc", "recovered_minus_true_pc", "right_within_0.25pc", "joint_fit_margin_lnL", "marked_for_animation_every_4th_row"], rows)
print(f"   sky 100: {sum(r[10] for r in rows)} of {len(rows)} right; marked {sum(r[12] for r in rows)}")

# ---------------- fisher_curves.csv -------------------------------------------------------------------------------------
P2 = np.load("/home/mattm/projects/HSYMT/CW_transition/prong2_results.npz"); rows = []
for pre, lab in (("ps", "fixed_per_source_strain"), ("tot", "fixed_total_power")):
    for k, n in enumerate(P2[pre+"_n"]):
        rows.append([lab, "8 log-uniform bins 3-20 nHz (prong2 main sweep)", int(n), int(P2[pre+"_n_seeds"])] + [f"{P2[pre+'_'+q][k]:.6e}" for q in ("cond", "cond_lo", "cond_hi", "marg", "marg_lo", "marg_hi")] + [f"{P2[pre+'_marg'][k]/P2[pre+'_cond'][k]:.5f}"])
for r in csv.DictReader(open("/home/mattm/projects/HSYMT/PAPER_figures/job48b_turnover_per_bin.csv")):
    lab = {"B8": "confined to 8 bins", "B4": "confined to 4 bins", "B2": "confined to 2 bins", "B1": "confined to 1 bin", "SF": "one frequency"}[r["arm"]]
    rows.append(["fixed_per_source_strain_job48b", lab+f" (B = {r['B']})", int(r["N"]), int(r["n_skies"]), r["I_cond_median"], r["I_cond_p16"], r["I_cond_p84"], r["I_marg_median"], r["I_marg_p16"], r["I_marg_p84"], r["ratio_median_of_medians"]])
W("fisher_curves.csv", ["mode", "band", "n_sources", "n_skies", "info_conditional_median_pc^-2", "info_conditional_p16", "info_conditional_p84", "info_marginal_median_pc^-2", "info_marginal_p16", "info_marginal_p84", "marginal_over_conditional"], rows)
