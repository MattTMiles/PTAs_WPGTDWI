"""JOB 54 post-processing: tidy table, figures, numbers for the draft, skymap file.
usage: job54_loc_post.py [tables|maps]"""
import numpy as np, os, sys, glob, json, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import job54_loc_core as C

HERE = C.HERE; SKY = os.path.join(C.DATA, "sky")
NPH = C.NPH_LIST; SIGLS = (0.1, 0.3, 1.0, 3.0, 10.0); NPART = (4, 10, 116)
FA = np.load(os.path.join(C.DATA, "fisher_areas.npz"))
FB = np.load(os.path.join(C.DATA, "fisher_blocks.npz"), allow_pickle=True)


def q(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return (np.nan,) * 3 if len(x) == 0 else (np.median(x), np.percentile(x, 16), np.percentile(x, 84))


def load_b(f, snr, kind="main", tag=""):
    """dict of arrays [n_src, K] plus list of source indices."""
    fs = sorted(glob.glob(os.path.join(SKY, f"{kind}_s???_f{f}_snr{snr}{tag}.npz")))
    if not fs:
        return None
    out = {"isrc": []}
    for fn in fs:
        d = np.load(fn, allow_pickle=True); out["isrc"].append(int(d["isrc"]))
        for k in ("A90", "A99", "n_isl", "A_truth_isl", "truth_in", "zoom_frac", "n_base_members_out", "pmin_over_h", "max_rel"):
            out.setdefault(k, []).append(d[k])
        for k in ("h_deg", "n_cells", "a_mask_deg2"):
            out.setdefault(k, []).append(float(d[k]))
    out = {k: np.array(v) for k, v in out.items()}
    # validity rules (see draft): area valid if the base grid resolves the phased set (>= 4 cells per finest fringe),
    # or >= 90% of the probability lies in the zoom window (64x finer), or the 90% set holds >= 1000 base cells
    # (broad set; an under-sampled fringe is then sampled at quasi-random phase).
    res = out["pmin_over_h"] >= 4
    out["area_ok"] = res | (out["zoom_frac"] >= 0.9) | (out["n_base_members_out"] >= 1000)
    out["isl_ok"] = res | (out["n_base_members_out"] == 0)
    out["A90v"] = np.where(out["area_ok"], out["A90"], np.nan)
    out["n_islv"] = np.where(out["isl_ok"], out["n_isl"], np.nan)
    out["A_truthv"] = np.where(out["isl_ok"], out["A_truth_isl"], np.nan)
    return out


def tables():
    rows = []; J = {}
    hdr = ["f_gw_nHz", "snr", "ordering", "n_ph", "sigma_L_pc", "measure", "n_sources", "area90_deg2_median",
           "area90_deg2_p16", "area90_deg2_p84", "n_islands_median", "truth_island_area_deg2_median",
           "n_sources_unresolved", "frac_truth_in_90pct_set"]
    # ---- measure (a)
    for f in (10, 20):
        for snr in (10, 20):
            for o in ("near", "rand"):
                for sl in (0.01,) + SIGLS:
                    for n in sorted(set(NPH + [30])):
                        for mode, cond, name in (("free", "", "a"), ("free", "_cond", "a_cond"), ("census", "", "a_census")):
                            A = FA[f"f{f}_snr{snr}_{o}_sl{sl}_n{n}_{mode}{cond}"]
                            m, lo, hi = q(A)
                            rows.append([f, snr, {"near": "nearest_first", "rand": "random"}[o], n, sl, name, len(A), m, lo, hi, "", "", 0, ""])
    # ---- measure (b)
    B = {}
    for f in (10, 20):
        for snr in (10, 20):
            b = load_b(f, snr)
            if b is None:
                continue
            B[(f, snr)] = b
            for io, o in enumerate(("nearest_first", "random")):
                for j, n in enumerate(NPH):
                    k = io * 11 + j
                    m, lo, hi = q(b["A90v"][:, k])
                    rows.append([f, snr, o, n, 0.01, "b", int(np.isfinite(b["A90v"][:, k]).sum()), m, lo, hi,
                                 q(b["n_islv"][:, k])[0], q(b["A_truthv"][:, k])[0], int((~b["area_ok"][:, k]).sum()),
                                 float(np.mean(b["truth_in"][:, k]))])
    bp = load_b(10, 10, "partial")
    if bp is not None:
        k = 1
        for sl in SIGLS:
            for n in NPART:
                m, lo, hi = q(bp["A90v"][:, k])
                rows.append([10, 10, "nearest_first", n, sl, "b", int(np.isfinite(bp["A90v"][:, k]).sum()), m, lo, hi,
                             q(bp["n_islv"][:, k])[0], q(bp["A_truthv"][:, k])[0], int((~bp["area_ok"][:, k]).sum()),
                             float(np.mean(bp["truth_in"][:, k]))])
                k += 1
    os.makedirs(os.path.join(HERE, "talk_data"), exist_ok=True)
    with open(os.path.join(HERE, "talk_data", "localisation.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(hdr)
        for r in rows:
            w.writerow([f"{x:.6g}" if isinstance(x, (float, np.floating)) else x for x in r])
    R = {(r[0], r[1], r[2], r[3], r[4], r[5]): r for r in rows}

    def cell(key, isl=False):
        if key not in R or not np.isfinite(R[key][7]):
            return "–"
        r = R[key]; s = f"{r[7]:.3g} [{r[8]:.3g}, {r[9]:.3g}]"
        if isl and r[10] != "":
            s += f" / {r[10]:.0f}" if np.isfinite(r[10]) else " / n.r."
        return s
    md = []
    for f in (10, 20):
        for snr in (10, 20):
            if (f, snr) not in B:
                continue
            nb = len(B[(f, snr)]["isrc"])
            md.append(f"\n**{f} nHz, SNR {snr}.** 90% area in deg^2, median [16th, 84th percentile]. (a): 200 sources. (b): {nb} sources; "
                      "after the slash, the median number of islands (n.r. = not resolved).\n")
            md.append("| n_ph | (a) nearest | (a) random | (a_cond) nearest | (b) nearest / islands | (b) random / islands | (b) truth-island area, nearest (deg^2) | (b) cases unresolved near/rand |")
            md.append("|---|---|---|---|---|---|---|---|")
            for n in NPH:
                kb = (f, snr, "nearest_first", n, 0.01, "b"); kr = (f, snr, "random", n, 0.01, "b")
                ti = R[kb][11]
                md.append(f"| {n if n < 116 else 'all (116)'} | {cell((f, snr, 'nearest_first', n, 0.01, 'a'))} | {cell((f, snr, 'random', n, 0.01, 'a'))} | "
                          f"{cell((f, snr, 'nearest_first', n, 0.01, 'a_cond'))} | {cell(kb, True)} | {cell(kr, True)} | "
                          f"{ti:.3g} | {R[kb][12]}/{R[kr][12]} |")
    # partial table
    md.append("\n**Partial knowledge, 10 nHz, SNR 10, nearest-first.** 90% area in deg^2, median [16th, 84th]. "
              "(a_lin) = linearised Fisher with a Gaussian distance prior sigma_L on the phased set (200 sources). "
              f"(b) = distance marginalised over the prior ({0 if bp is None else len(bp['isrc'])} sources) / median islands.\n")
    md.append("| sigma_L (pc) | n_ph=4 (a_lin) | n_ph=4 (b) | n_ph=10 (a_lin) | n_ph=10 (b) | all (a_lin) | all (b) |")
    md.append("|---|---|---|---|---|---|---|")
    for sl in (0.01,) + SIGLS:
        md.append(f"| {sl} | " + " | ".join(cell((10, 10, "nearest_first", n, sl, m), m == "b") for n in NPART for m in ("a", "b")) + " |")
    open(os.path.join(C.DATA, "tables.md"), "w").write("\n".join(md))

    # ---- numbers for the draft
    def med(f, snr, o, n, m, sl=0.01):
        k = (f, snr, o, n, sl, m)
        return R[k][7] if k in R else np.nan
    for f in (10, 20):
        for snr in (10, 20):
            for o in ("nearest_first", "random"):
                for m in ("a", "b"):
                    A0, Aall = med(f, snr, o, 0, m), med(f, snr, o, 116, m)
                    if not np.isfinite(A0):
                        continue
                    e = {"A0": A0, "Aall": Aall, "full_gain_dex": float(np.log10(A0 / Aall))}
                    for n in (3, 10, 30):
                        if m == "a":
                            An = med(f, snr, o, n, m)
                        elif n == 30:
                            An = float(np.exp(np.interp(np.log(30), np.log([20, 40]), np.log([med(f, snr, o, 20, m), med(f, snr, o, 40, m)]))))
                        else:
                            An = med(f, snr, o, n, m)
                        e[f"n{n}"] = dict(area=An, factor=A0 / An, frac_of_log_gain=float(np.log10(A0 / An) / np.log10(A0 / Aall)))
                    # first n_ph within 10x of all-phased
                    e["first_n_within_10x_of_all"] = next((n for n in NPH if med(f, snr, o, n, m) <= 10 * Aall), None)
                    J[f"gain_f{f}_snr{snr}_{o}_{m}"] = e
    # per-source checks
    for (f, snr), b in B.items():
        src = b["isrc"]
        fa = lambda n, suf: FA[f"f{f}_snr{snr}_near_sl0.01_n{n}_free{suf}"][src]
        e = {"n_src": len(src), "sources": src.tolist()}
        for n, j in ((0, 0), (116, 10), (80, 9)):
            e[f"b_over_a_marg_n{n}"] = q(b["A90"][:, j] / fa(n, ""))
            e[f"b_over_a_cond_n{n}"] = q(b["A90"][:, j] / fa(n, "_cond"))
        e["b_n0_A99_median_deg2"] = float(np.median(b["A99"][:, 0])); e["mask_deg2_median"] = float(np.median(b["a_mask_deg2"]))
        e["h_deg_range"] = [float(b["h_deg"].min()), float(b["h_deg"].max())]
        e["max_rel_max"] = float(b["max_rel"].max())
        # W2: 3 random phased
        r3 = b["A90v"][:, 0] / b["A90v"][:, 11 + 3]; e["W2_ratio_n0_over_n3_rand"] = q(r3)
        e["W2_ratio_of_medians_rand"] = float(np.nanmedian(b["A90v"][:, 0]) / np.nanmedian(b["A90v"][:, 14]))
        vi = np.isfinite(b["n_islv"][:, 14]); e["W2_n_src_islands_resolved_rand_n3"] = int(vi.sum())
        e["W2_frac_multi_island_rand_n3"] = float(np.mean(b["n_islv"][vi, 14] > 1)) if vi.any() else float("nan"); e["W2_islands_median_rand_n3"] = q(b["n_islv"][:, 14])[0]
        e["W2_ratio_n0_over_n3_near"] = q(b["A90v"][:, 0] / b["A90v"][:, 3]); e["W2_frac_multi_island_near_n3"] = float(np.nanmean(b["n_islv"][:, 3] > 1))
        # W3: b / a per n_ph (same sources), nearest & random
        for io, o in enumerate(("near", "rand")):
            e[f"W3_b_over_a_{o}"] = {str(n): q(b["A90v"][:, io * 11 + j] / FA[f"f{f}_snr{snr}_{o}_sl0.01_n{n}_free"][src])[0] for j, n in enumerate(NPH)}
            e[f"W3_b_over_acond_{o}"] = {str(n): q(b["A90v"][:, io * 11 + j] / FA[f"f{f}_snr{snr}_{o}_sl0.01_n{n}_free_cond"][src])[0] for j, n in enumerate(NPH)}
        e["frac_truth_in_near"] = b["truth_in"][:, :11].mean(0).tolist()
        e["zoom_frac_median_near"] = np.median(b["zoom_frac"][:, :11], 0).tolist()
        e["pmin_over_h_median_near"] = np.median(b["pmin_over_h"][:, :11], 0).tolist(); e["pmin_over_h_median_rand"] = np.median(b["pmin_over_h"][:, 11:], 0).tolist()
        J[f"check_f{f}_snr{snr}"] = e
    if bp is not None:
        J["partial"] = {"sources": bp["isrc"].tolist(), "h_deg": bp["h_deg"].tolist(),
                        "W5_gain_n10_sl3_b": float(np.nanmedian(bp["A90v"][:, 0]) / np.nanmedian(bp["A90v"][:, 1 + 3 * 3 + 1])),
                        "W5_gain_n10_sl3_b_persrc": q(bp["A90v"][:, 0] / bp["A90v"][:, 1 + 3 * 3 + 1]),
                        "gain_b_persrc": {f"sl{sl}_n{n}": q(bp["A90v"][:, 0] / bp["A90v"][:, 1 + i * 3 + j])[0] for i, sl in enumerate(SIGLS) for j, n in enumerate(NPART)},
                        "max_rel_max": float(bp["max_rel"].max())}
    J["W5_gain_n10_sl3_a_lin"] = float(med(10, 10, "nearest_first", 0, "a") / med(10, 10, "nearest_first", 10, "a", 3.0))
    # half-step test
    fh = os.path.join(SKY, "main_s002_f10_snr10_half.npz"); f0 = os.path.join(SKY, "main_s002_f10_snr10.npz")
    if os.path.exists(fh) and os.path.exists(f0):
        a, b_ = np.load(f0), np.load(fh)
        J["half_step"] = dict(h_deg=[float(a["h_deg"]), float(b_["h_deg"])], n_cells=[int(a["n_cells"]), int(b_["n_cells"])],
                              A90_base=a["A90"].tolist(), A90_half=b_["A90"].tolist(), ratio=(b_["A90"] / a["A90"]).tolist(),
                              isl_base=a["n_isl"].tolist(), isl_half=b_["n_isl"].tolist(), pmin_over_h_base=a["pmin_over_h"].tolist(),
                              truth_isl_base=a["A_truth_isl"].tolist(), truth_isl_half=b_["A_truth_isl"].tolist())
    # census counts for the drawn array
    sg = FB["sigL_pc"]; L = FB["L_pc"]
    J["census"] = {"n_lt_0.4pc": int((sg < 0.4).sum()), "n_lt_1pc": int((sg < 1).sum()), "n_le_1pc": int((sg <= 1).sum()), "n_lt_3pc": int((sg < 3).sum()),
                   "L_pc_min_med_max": [float(L.min()), float(np.median(L)), float(L.max())], "sig_pc_min_med_max": [float(sg.min()), float(np.median(sg)), float(sg.max())],
                   "L_sorted_pc": np.sort(L).tolist(), "n_unique_census_rows": int(len(np.unique(FB["census_idx"])))}
    # spread over which pulsars are phased: 30 random nested permutations, Fisher (a), 10 nHz SNR 10
    import job54_loc_fisher as Fm
    Ap, bb, cc = FB["f10_Aperp"], FB["f10_b"], FB["f10_c"]; th = FB["src_theta"]; sp = {}
    for n in (3, 4, 10, 30):
        meds = []
        for p in range(30):
            order = np.random.default_rng(900 + p).permutation(C.N_PSR)
            sig = np.full(C.N_PSR, np.inf); sig[order[:n]] = 0.01
            meds.append(np.median([Fm.area_from_F(Fm.F8(Ap[i], bb[i], cc[i], sig, 10.0), th[i]) for i in range(0, 200, 4)]))
        sp[str(n)] = dict(min=float(np.min(meds)), p16=float(np.percentile(meds, 16)), med=float(np.median(meds)), p84=float(np.percentile(meds, 84)), max=float(np.max(meds)),
                          nearest=float(np.median(FA[f"f10_snr10_near_sl0.01_n{n}_free"][::4])))
    J["spread_which_pulsars_fisher_f10_snr10"] = sp
    json.dump(J, open(os.path.join(C.DATA, "numbers.json"), "w"), indent=1, default=float)

    # ---- figures
    cA, cB = "#2a78d6", "#eb6834"; X0 = 0.45
    xs = lambda ns: np.array([X0 if n == 0 else n for n in ns], float)
    sc = C.N_PSR / 79.0
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.5), sharex=True, sharey=True)
    for ia, snr in enumerate((10, 20)):
        for ib, f in enumerate((10, 20)):
            ax = axs[ia, ib]
            rowsm = [(2.5e6, "real priors", (0 * sc, 0.34 * sc), (2.3 * sc, 3.6 * sc)),
                     (4e5, "Earth terms removed", (0.39 * sc, 1.5 * sc), (3.9 * sc, 6.6 * sc)),
                     (6e4, "all sources known, Earth terms removed", (26.8 * sc, 26.8 * sc), (39.6 * sc, 39.6 * sc))]
            for yl, lab, cr, rr in rowsm:
                for (lo, hi), col, tg in ((cr, "#8f8e84", "C"), (rr, "#3d3d3a", "R")):
                    lo2 = max(lo, X0)
                    ax.plot([lo2, max(hi, lo2)], [yl, yl], color=col, lw=5, solid_capstyle="butt", marker="|" if lo == hi else None, ms=10, mew=3)
                    ax.plot([np.sqrt(lo2 * max(hi, lo2))] * 2, [2e-7, yl], color=col, lw=0.7, ls=":")
                    ax.text(np.sqrt(lo2 * max(hi, lo2)), yl * 1.5, tg, fontsize=7, ha="center", va="bottom", color="#3d3d3a")
                ax.text(140, yl, lab, fontsize=7, ha="right", va="center", color="#3d3d3a") if lo != hi else ax.text(22, yl, lab, fontsize=7, ha="right", va="center", color="#3d3d3a")
            for o, ls, mk in (("nearest_first", "-", "o"), ("random", "--", "s")):
                na = sorted(set(NPH + [30]))
                v = np.array([R[(f, snr, o, n, 0.01, "a")][7:10] for n in na], float)
                ax.fill_between(xs(na), v[:, 1], v[:, 2], color=cA, alpha=0.12, lw=0)
                ax.plot(xs(na), v[:, 0], ls, color=cA, lw=2, marker=mk, ms=5, label=f"(a) Fisher, {o.replace('_', '-')}")
                if (f, snr, o, 0, 0.01, "b") in R:
                    v = np.array([R[(f, snr, o, n, 0.01, "b")][7:10] for n in NPH], float)
                    ax.fill_between(xs(NPH), v[:, 1], v[:, 2], color=cB, alpha=0.15, lw=0)
                    ax.plot(xs(NPH), v[:, 0], ls, color=cB, lw=2, marker=mk, ms=5, label=f"(b) sky likelihood, {o.replace('_', '-')}")
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(X0 * 0.8, 150); ax.set_ylim(2e-7, 1.5e7)
            ax.set_xticks([X0, 1, 2, 3, 4, 6, 10, 20, 40, 80, 116]); ax.set_xticklabels(["0", "1", "2", "3", "4", "6", "10", "20", "40", "80", "all"])
            ax.axvline(0.7, color="#c3c2b7", lw=0.8); ax.minorticks_off()
            ax.grid(True, which="major", color="#e4e3dc", lw=0.6); ax.set_axisbelow(True)
            for sname in ("top", "right"):
                ax.spines[sname].set_visible(False)
            ax.set_title(f"{f} nHz, SNR {snr}" + ("" if (f, snr) in B else "  [(b) not run]"), fontsize=11, loc="left")
            if ia == 1:
                ax.set_xlabel("number of phased pulsars n_ph (of 116; sigma_L = 0.01 pc)")
            if ib == 0:
                ax.set_ylabel("90% sky area (deg$^2$)")
            if ia == 0 and ib == 0:
                ax.legend(fontsize=8, loc="lower left", frameon=False, bbox_to_anchor=(0.0, 0.05))
                ax.text(0.5, 3e-7, "top bars: supply of distances per sky, scaled x116/79 (C = certified at 90% purity, R = recovered)", fontsize=6.5, color="#3d3d3a", va="bottom")
    fig.suptitle("90% sky area vs number of phased pulsars. Lines: median over sources; bands: 16–84%. (a) 200 sources, (b) 20 sources.", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "job54_fig1_area_vs_nph.png"), dpi=140); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(13, 4.6), sharey=True)
    sls = (0.01,) + SIGLS
    for ax, n in zip(axs, NPART):
        v = np.array([R[(10, 10, "nearest_first", n, sl, "a")][7:10] for sl in sls], float)
        ax.fill_between(sls, v[:, 1], v[:, 2], color=cA, alpha=0.12, lw=0)
        ax.plot(sls, v[:, 0], "-o", color=cA, lw=2, ms=5, label="(a) linearised Fisher, Gaussian prior")
        if (10, 10, "nearest_first", n, 0.1, "b") in R:
            v = np.array([R[(10, 10, "nearest_first", n, sl, "b")][7:10] for sl in sls], float)
            ax.fill_between(sls, v[:, 1], v[:, 2], color=cB, alpha=0.15, lw=0)
            ax.plot(sls, v[:, 0], "-s", color=cB, lw=2, ms=5, label="(b) sky likelihood, distance marginalised")
        a0 = R[(10, 10, "nearest_first", 0, 0.01, "b")][7] if (10, 10, "nearest_first", 0, 0.01, "b") in R else np.nan
        ax.axhline(a0, color="#6b6a63", lw=1, ls=":"); ax.text(0.011, a0 * 1.3, "n_ph = 0, measure (b)", fontsize=8, color="#3d3d3a")
        ax.axvline(299792458.0 / 1e-8 / 3.0856775814913673e16 / (2 * np.pi), color="#c3c2b7", lw=1)
        ax.text(0.16, 3e-6, "lambda/2pi\n= 0.155 pc", fontsize=8, color="#3d3d3a")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(f"n_ph = {n if n < 116 else 'all (116)'}, nearest-first", fontsize=11, loc="left")
        ax.set_xlabel("distance error of the phased set, sigma_L (pc)"); ax.grid(True, color="#e4e3dc", lw=0.6); ax.set_axisbelow(True)
        for sname in ("top", "right"):
            ax.spines[sname].set_visible(False)
    axs[0].set_ylabel("90% sky area (deg$^2$)"); axs[0].legend(fontsize=8, frameon=False, loc="center left")
    fig.suptitle("10 nHz, SNR 10: 90% sky area vs distance error of the phased set (unphased pulsars: phase unknown). Median and 16–84% band.", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "job54_fig2_area_vs_sigmaL.png"), dpi=140); plt.close(fig)
    # typical source rule
    if (10, 10) in B:
        b = B[(10, 10)]
        la = np.log10(b["A90"][:, :11]); dist = np.abs(la - np.median(la, 0)).sum(1)
        J["typical"] = dict(isrc=int(b["isrc"][np.argmin(dist)]), score_dex_sum=float(dist.min()), scores=dist.tolist())
        json.dump(J, open(os.path.join(C.DATA, "numbers.json"), "w"), indent=1, default=float)
        print("typical source", J["typical"]["isrc"])
    print(open(os.path.join(C.DATA, "tables.md")).read())


def maps():
    J = json.load(open(os.path.join(C.DATA, "numbers.json"))); i = J["typical"]["isrc"]
    d = np.load(os.path.join(SKY, f"main_s{i:03d}_f10_snr10_maps.npz"), allow_pickle=True)
    d0 = np.load(os.path.join(SKY, f"main_s{i:03d}_f10_snr10.npz"), allow_pickle=True)
    out = dict(n_ph=np.array(NPH), lnL_rel_max=d["disp"], u_deg=d["disp_u_deg"].astype(np.float32), v_deg=d["disp_v_deg"].astype(np.float32),
               ra_deg=d["disp_ra_deg"], dec_deg=d["disp_dec_deg"], true_ra_deg=d["true_ra_deg"], true_dec_deg=d["true_dec_deg"],
               area90_deg2=d0["A90"][:11], n_islands=d0["n_isl"][:11], area90_truth_island_deg2=d0["A_truth_isl"][:11],
               zoom_n_ph=np.array(NPH[6:]), zoom_lnL_rel_max=d["zoom"][6:], zoom_u_deg=d["zoom_u_deg"].astype(np.float32), zoom_v_deg=d["zoom_v_deg"].astype(np.float32),
               zoom_ra_deg=d["zoom_ra_deg"], zoom_dec_deg=d["zoom_dec_deg"], source_index=i, f_gw_nHz=10.0, snr=10.0,
               fine_step_deg=float(d0["h_deg"]), zoom_fine_step_deg=float(d["zoom_fine_step_deg"]), phased_order=FB["order_near"],
               readme=np.array("lnL_rel_max[k]: expected (noiseless) log-likelihood minus its value at the true position, for n_ph[k] nearest-first (nested) phased pulsars; "
                               "each display cell holds the MAXIMUM of the fine evaluation over that cell (float16). Axes: v_deg (rows) x u_deg (cols) = local lat/lon offsets "
                               "from the true position in a frame rotated so the truth is at (0,0); u east, v north; ra_deg/dec_deg give RA/Dec of each cell centre. "
                               "area90_deg2 etc. come from the fine computation (step fine_step_deg, plus 64x zoom at the truth), not from these display grids. "
                               "zoom_*: same for n_ph >= 10 on a grid 30x finer (linear) centred on the truth; block max of an 8x finer evaluation."))
    fn = os.path.join(HERE, "talk_data", "skymaps.npz"); np.savez_compressed(fn, **out)
    print("skymaps", os.path.getsize(fn) / 1e6, "MB; peak of each map", d["disp"].astype(float).max((1, 2)), d["zoom"].astype(float).max((1, 2)))


if __name__ == "__main__":
    (maps if len(sys.argv) > 1 and sys.argv[1] == "maps" else tables)()
