"""EASEL-4 — six explainer figures, banked npz only.

READ-ONLY. No likelihood evaluation, no new realisations. Every number is either read
straight out of a banked npz or re-cut from banked raw statistics with the campaign's own
definitions (criterion-v2.1: layer 1 dlnL > lnK, layer 2 dlnL > floor, layer 3 q_max > 0.9;
offender = max dlnL among pulsars passing layers 1+3).

FLOOR CONVENTION (RECUT, adopted — supersedes the Gumbel-everywhere convention EASEL-2 used):
the Gumbel-MLE alpha=0.05 quantile is valid only where the null's zero-fraction is <= 20 %
(ANCHOR §4). Above that it describes a point mass at zero, not a tail, and the floor is the
EMPIRICAL 95th percentile with a bootstrap error. F3/F4/F5 therefore read `floor_adopted`
(NOT `fN`, which still carries the Gumbel) and the already-re-cut `corr`/`corr_lo` columns
from the _recut banks. `corr_lo` = count at floor + its error = the onset/CONFIRMED test.

Run:  python explainer_figures.py
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle
from scipy.stats import gumbel_r

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"
OUT = f"{REPO}/EXPLAINER_results"

# criterion-v2.1 constants, verbatim from hpc_harbor/surface/surface.py
QBAR = 0.9
ALPHA = 0.05
Z_ALPHA = -np.log(-np.log(1.0 - ALPHA))   # 2.9702
C_SD = 2.80                               # sd(floor) = C_SD * beta / sqrt(N)

plt.rcParams.update({
    "font.size": 15, "axes.labelsize": 17, "axes.titlesize": 18,
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 130, "savefig.bbox": "tight", "axes.linewidth": 1.2,
})
INK = "#1b1b1b"
BLUE = "#2166ac"
RED = "#b2182b"
GREEN = "#1a7d3f"
AMBER = "#d98c00"
GREY = "#8c8c8c"


def offender(dlnL, lnK, qmax):
    """Campaign's offender statistic: largest dlnL among pulsars passing layers 1+3."""
    m = (dlnL > lnK) & (qmax > QBAR)
    return float(dlnL[m].max()) if m.any() else 0.0


def load_cell(pat):
    rows = []
    for f in sorted(glob.glob(pat)):
        d = np.load(f, allow_pickle=True)
        rows.append((np.nan_to_num(d["dlnL_det"], posinf=1e30),
                     d["lnK"], d["qmax"], d["on_true"]))
    return rows


# ======================================================================
# F1 — THE PROBLEM: one pulsar's fringe comb
# ======================================================================
def f1():
    scan = np.load(f"{REPO}/lnLs_GWAmp_phase_connected/runD_3CW_test/J0437-4715_cw_p_dist/"
                   f"psrTerm/runD_3CW_test_J0437-4715_cw_p_dist_psrTerm_0_0.npz",
                   allow_pickle=True)
    D = scan["scan_values"]                 # kpc
    y = scan["logls"][0]
    y = y - y.max()                         # relative log-likelihood

    a0 = np.load(f"{REPO}/CW_transition/anchor_a0_priors.npz", allow_pickle=True)
    i = list(a0["fname"]).index("J0437-4715")
    D_em = float(a0["dist_feather_kpc"][i])      # 0.1549 kpc
    s_em = float(a0["sig_feather_kpc"][i])       # 0.001 kpc = 1 pc

    geo = np.load(f"{REPO}/reports/geo_dlnl_bank.npz", allow_pickle=True)
    names = list(geo["names"])
    Kmed = np.median(geo["K_counted"], axis=0)   # per-pulsar median over 40 sky draws
    K_j0437 = Kmed[names.index("J0437-4715")]

    fig, ax = plt.subplots(1, 3, figsize=(21, 6.4))
    fig.subplots_adjust(wspace=0.42)
    box = dict(boxstyle="round,pad=0.45", fc="white", ec=GREY, alpha=0.93)

    # (a) the whole comb
    ax[0].plot(D, y, lw=0.4, color=BLUE, alpha=0.85)
    ax[0].axvline(D_em, color=RED, lw=2.2, zorder=5)
    ax[0].set_xlabel("Distance to the pulsar  (kpc)")
    ax[0].set_ylabel("How well the model fits\n(relative log-likelihood)")
    ax[0].set_title("(a) The likelihood is a comb", loc="left", weight="bold")
    ax[0].set_xlim(0.1, 5.0)
    ax[0].set_ylim(-2350, 620)
    ax[0].text(0.50, 0.055, "about 810 near-identical peaks —\nevery one fits almost as well",
               transform=ax[0].transAxes, ha="center", fontsize=13.5, color=INK, bbox=box)
    ax[0].annotate("true distance", xy=(D_em, 180), xytext=(1.05, 430),
                   fontsize=14, color=RED, weight="bold",
                   arrowprops=dict(arrowstyle="->", color=RED, lw=1.8))

    # (b) zoom onto the EM prior window
    lo, hi = 0.100, 0.225
    m = (D >= lo) & (D <= hi)
    Dz, yz = D[m], y[m]
    ax[1].plot(Dz, yz, lw=1.6, color=BLUE, marker="o", ms=2.6, mfc=BLUE, mec="none")
    pk = np.where((yz[1:-1] > yz[:-2]) & (yz[1:-1] > yz[2:]))[0] + 1
    ax[1].plot(Dz[pk], yz[pk], "v", ms=10, color=AMBER, mec=INK, mew=0.6,
               label="candidate slots (comb teeth)", ls="none", zorder=5)
    ax[1].axvspan(D_em - s_em, D_em + s_em, color=RED, alpha=0.30, zorder=0,
                  label="what telescopes already know\n(distance prior, ±1 pc)")
    ax[1].axvline(D_em, color=RED, lw=2.0)
    ax[1].set_xlabel("Distance to the pulsar  (kpc)")
    ax[1].set_ylabel("How well the model fits\n(relative log-likelihood)")
    ax[1].set_title("(b) Which tooth is the pulsar on?", loc="left", weight="bold")
    ax[1].set_xlim(lo, hi)
    ax[1].set_ylim(-2250, 1150)
    ax[1].legend(loc="lower center", framealpha=0.96, fontsize=12)
    ax[1].text(0.50, 0.965, "teeth 5.9 pc apart here, so the prior\n"
                            "nearly isolates one — the EASY case",
               transform=ax[1].transAxes, ha="center", va="top",
               fontsize=13, color=INK, bbox=box)
    ax[1].text(0.50, -0.30, "NOTE: this is a banked legacy single-source scan — its comb is\n"
                            "COARSER than the census. Panel (c) is the real distribution.",
               transform=ax[1].transAxes, ha="center", va="top", fontsize=11.5,
               color=RED, style="italic")

    # (c) how many teeth actually sit inside the prior window, across the real array
    ax[2].hist(np.log10(np.maximum(Kmed, 1)), bins=26, color=BLUE, alpha=0.80,
               edgecolor="white", lw=0.7)
    ax[2].axvline(np.log10(K_j0437), color=GREEN, lw=3.2, zorder=5)
    ax[2].axvline(np.log10(np.median(Kmed)), color=RED, lw=3.2, ls="--", zorder=5)
    ax[2].set_xlabel("Candidate slots inside the distance prior")
    ax[2].set_ylabel("Number of pulsars")
    ax[2].set_title("(c) The real array is not the easy case", loc="left", weight="bold")
    ax[2].set_xticks([0, 1, 2, 3, 4, 5])
    ax[2].set_xticklabels(["1", "10", "100", "1,000", "10,000", "100,000"])
    ax[2].set_xlim(-0.4, 5.4)
    ax[2].set_ylim(0, 30)
    ax[2].annotate(f"J0437−4715:\n{K_j0437:.0f} slots\n(the anchor)",
                   xy=(np.log10(K_j0437), 19.5), xytext=(0.15, 24.0),
                   fontsize=13, color=GREEN, weight="bold", bbox=box,
                   arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.8))
    ax[2].annotate(f"typical pulsar:\n{np.median(Kmed):,.0f} slots",
                   xy=(np.log10(np.median(Kmed)), 19.5), xytext=(2.75, 25.4),
                   fontsize=13, color=RED, weight="bold", bbox=box,
                   arrowprops=dict(arrowstyle="->", color=RED, lw=1.8))

    fig.suptitle("F1 · THE PROBLEM — the signal pins the pulsar's distance, "
                 "but only to within one wavelength",
                 fontsize=20, weight="bold", y=1.02)
    fig.savefig(f"{OUT}/F1_the_problem.png")
    plt.close(fig)
    print(f"F1 ok  (K: J0437 {K_j0437:.0f}, census median {np.median(Kmed):,.0f})")


# ======================================================================
# F2 — THE WALL: the registration-tolerance ladder
# ======================================================================
def f2():
    F = np.load(f"{REPO}/CW_transition/trackB_F2_ladder.npz", allow_pickle=True)
    tol = np.asarray(F["tol_sky"], float)
    tol = tol[np.isfinite(tol)]

    LOOSEST = tol.max()          # 1.852e-3 — the most forgiving pulsar-source pair
    BLIND_LO, BLIND_HI = 0.05, 0.21   # F-stat seed floor -> soft-EM endpoint (banked)
    PULLIN = 1e-4                # fixed-integer pull-in basin (L2b/L2c, banked)

    fig, ax = plt.subplots(figsize=(15.5, 8.0))
    box = dict(boxstyle="round,pad=0.4", fc="white", ec=GREY, alpha=0.93)

    ax.axvspan(np.log10(LOOSEST), np.log10(BLIND_LO), color=GREY, alpha=0.22, zorder=0)
    ax.axvspan(np.log10(BLIND_LO), np.log10(BLIND_HI), color=RED, alpha=0.16, zorder=0)
    ax.axvspan(-6.8, np.log10(LOOSEST), color=GREEN, alpha=0.13, zorder=0)

    ax.hist(np.log10(tol), bins=34, color=BLUE, alpha=0.92, edgecolor="white", lw=0.7,
            zorder=3)
    ax.axvline(np.log10(LOOSEST), color=GREEN, lw=2.8, zorder=4)
    ax.axvline(np.log10(BLIND_LO), color=RED, lw=2.8, zorder=4)
    ax.axvline(np.log10(PULLIN), color=INK, lw=1.8, ls=":", zorder=4)

    ax.set_xlim(-6.5, -0.35)
    ax.set_ylim(0, 66)
    ax.set_xlabel("How accurately you must already know where the source sits on the sky\n"
                  "(further left = you must know it more precisely)")
    ax.set_ylabel("Number of pulsar–source pairs")
    ax.set_xticks(np.log10([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels(["$10^{-6}$", "$10^{-5}$", "$10^{-4}$", "$10^{-3}$",
                        "$10^{-2}$", "$10^{-1}$"])

    ax.text(-4.55, 61, "WHERE THE ANSWER LIVES\nall 348 pulsar–source pairs",
            ha="center", va="center", fontsize=15, weight="bold", color=GREEN, zorder=6)
    ax.text(np.log10(np.sqrt(LOOSEST * BLIND_LO)), 61,
            "NOTHING\nIN BETWEEN",
            ha="center", va="center", fontsize=15, weight="bold", color=INK, zorder=6)
    ax.text(np.log10(0.102), 61, "WHERE A BLIND\nSEARCH GETS YOU",
            ha="center", va="center", fontsize=15, weight="bold", color=RED, zorder=6)

    ax.annotate("zero of the 348 pairs\nland in this gap —\nit is 27–112× wide",
                xy=(np.log10(8e-3), 20), xytext=(np.log10(6.5e-3), 44),
                ha="center", fontsize=13.5, weight="bold", color=INK, bbox=box, zorder=6)
    ax.annotate("the single most forgiving pair\n(the nearest pulsar, J0711)",
                xy=(np.log10(LOOSEST), 11), xytext=(-4.30, 33),
                fontsize=12.5, color=GREEN, weight="bold", bbox=box, zorder=6,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.8))
    ax.annotate("locking onto the right tooth\nneeds this side of the line",
                xy=(np.log10(PULLIN), 6), xytext=(-6.35, 22),
                fontsize=12.5, color=INK, bbox=box, zorder=6,
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.5))

    ax.set_title("F2 · THE WALL — the accuracy you need and the accuracy you can get\n"
                 "do not overlap, and nothing bridges them",
                 loc="left", weight="bold", fontsize=19)
    fig.savefig(f"{OUT}/F2_the_wall.png")
    plt.close(fig)
    print(f"F2 ok  (loosest {LOOSEST:.2e}, blind {BLIND_LO}–{BLIND_HI}, "
          f"pairs in gap: {int(((tol > LOOSEST) & (tol < BLIND_LO)).sum())})")


# ======================================================================
# F3 — THE HONEST CRITERION: null vs signal at one cell
# ======================================================================
def f3():
    tag = "h1300_T30_lit_k3"          # h = -13.00, T = 30 yr, lit priors, census structure
    N = load_cell(f"{REPO}/SURFACE_results/sf_nullN_{tag}_g*_n*.npz")
    S = load_cell(f"{REPO}/SURFACE_results/sf_sig_{tag}_g*_n*.npz")

    off_null = np.array([offender(d, k, q) for d, k, q, _ in N])
    off_sig = np.array([offender(d, k, q) for d, k, q, _ in S])

    # The ADOPTED floor for this cell, straight from the re-cut bank. This cell's null has a
    # 27 % zero-fraction, so its Gumbel fails ANCHOR §4's validity gate and the adopted floor
    # is the empirical q95 +/- bootstrap. The Gumbel is kept only to draw the irony.
    R = np.load(f"{REPO}/reports/surface_analysis_recut.npz", allow_pickle=True)
    ci = {c: i for i, c in enumerate(str(c) for c in R["cols"])}
    T = R["table"]
    m = ((np.round(T[:, ci["h"]], 2) == -13.00) & (T[:, ci["T"]] == 30)
         & (R["tiers"] == "lit") & (R["struct"] == "3+13"))
    r = T[int(np.flatnonzero(m)[0])]
    floor = r[ci["floor_adopted"]]        # 16.60  (empirical q95)
    sd = r[ci["floor_adopted_sd"]]        # +/- 1.60 (bootstrap)
    gum, gum_sd = r[ci["fN"]], r[ci["fN_sd"]]   # 19.46 +/- 1.40 — INVALID here
    zf = r[ci["fN_zerofrac"]]             # 0.27 > 0.20 gate
    c_new = r[ci["corr"]]                 # 0.60, already re-cut against floor_adopted
    OLD = 13.48                       # IGNITE's max-of-10 floor for this cell (banked, SURFACE D-4)

    def counts(fl):
        c = w = 0
        for d, k, q, ot in S:
            cert = (d > np.maximum(k, fl)) & (q > QBAR)
            c += int((cert & ot).sum())
            w += int((cert & ~ot).sum())
        return c / len(S), w / len(S)

    c_old, _ = counts(OLD)
    # gate: this scorer, at the adopted floor, must reproduce the bank's re-cut count
    assert abs(counts(floor)[0] - c_new) < 1e-9, "F3 scorer disagrees with recut bank"

    fig, ax = plt.subplots(figsize=(14.5, 7.6))
    bins = np.linspace(0, 46, 32)
    ax.hist(off_null, bins=bins, color=GREY, alpha=0.85, edgecolor="white", lw=0.8,
            label=f"PURE NOISE, no source present  (N = {len(off_null)})")
    ax.hist(off_sig, bins=bins, color=BLUE, alpha=0.80, edgecolor="white", lw=0.8,
            label=f"REAL SOURCE injected  (N = {len(off_sig)})")

    ax.axvspan(floor - sd, floor + sd, color=GREEN, alpha=0.20, zorder=1)
    ax.axvline(floor, color=GREEN, lw=3.0, zorder=4)
    ax.axvline(OLD, color=RED, lw=3.0, ls="--", zorder=4)
    # the bar this figure USED to draw — a Gumbel this cell's own zero-fraction forbids
    ax.axvline(gum, color=GREY, lw=1.8, ls=":", zorder=3)

    ax.set_xlabel("Best evidence for a certified pulsar in one dataset\n"
                  "(how much better the best candidate fits than no source, in nats)")
    ax.set_ylabel("Number of datasets")
    ax.set_xlim(0, 46)
    ax.set_ylim(0, 34)

    ax.annotate(f"THE OLD BAR ({OLD:.1f})\nset by the loudest of just 10\nnoise-only trials\n"
                f"→ {c_old:.2f} 'certifications' per dataset\n→ looked like a detection",
                xy=(OLD, 22), xytext=(0.6, 25.4), fontsize=13.5, color=RED, weight="bold",
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.8))
    ax.annotate(f"THE HONEST BAR ({floor:.1f} ± {sd:.1f})\nthe 95th percentile pure noise\n"
                f"actually reaches (N = 100)\n→ {c_new:.2f} per dataset → RETRACTED",
                xy=(floor, 13), xytext=(24.5, 20.5), fontsize=13.5, color=GREEN, weight="bold",
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.8))
    ax.annotate("noise alone routinely\nreaches this far",
                xy=(17.6, 2.9), xytext=(22.0, 7.6), fontsize=13, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.4))
    ax.annotate(f"the fitted bar this figure\nused to draw ({gum:.1f}) — see below",
                xy=(gum, 25.5), xytext=(22.5, 25.8), fontsize=12, color=GREY,
                arrowprops=dict(arrowstyle="->", color=GREY, lw=1.4))

    # the retraction is FIRMER, and this is why: the corrected bar is LOWER — more permissive,
    # the direction that would have rescued the detection — and the count still does not reach 1.
    ax.text(0.5, -0.27,
            f"The corrected bar is LOWER than the fitted one it replaces "
            f"({floor:.1f} < {gum:.1f}) — the change runs in the direction that\n"
            f"would have RESCUED the detection, and the count still fails to reach 1. "
            f"The retraction no longer\nleans on the stricter bar: it is FIRMER than when this "
            f"figure was first drawn.",
            transform=ax.transAxes, ha="center", va="top", fontsize=12.5, color=INK,
            bbox=dict(boxstyle="round,pad=0.45", fc="#eef5ee", ec=GREEN, alpha=0.95))
    ax.text(0.5, -0.54,
            f"Footnote, and the figure owes it: the very cell chosen to teach honest "
            f"floor-fitting had a dishonest floor.\n{zf:.0%} of its noise trials certify nothing "
            f"at all, so the Gumbel fit (grey dotted) describes a point mass at zero\n"
            f"rather than a tail — ANCHOR §4 forbids it here. The green bar is the empirical "
            f"percentile instead.",
            transform=ax.transAxes, ha="center", va="top", fontsize=11.5, color=INK,
            style="italic",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f4f4f4", ec=GREY, alpha=0.95))

    ax.legend(loc="upper right", framealpha=0.96)
    ax.set_title("F3 · THE HONEST CRITERION — the same evidence, two different bars\n"
                 "(one cell: source strength h = $10^{-13}$, 30-year array)",
                 loc="left", weight="bold", fontsize=19)
    fig.savefig(f"{OUT}/F3_the_honest_criterion.png")
    plt.close(fig)
    print(f"F3 ok  (adopted floor {floor:.2f}+/-{sd:.2f} [emp q95, zf={zf:.2f}] vs "
          f"INVALID Gumbel {gum:.2f}+/-{gum_sd:.2f}; counts {c_old:.2f} -> {c_new:.2f})")


# ======================================================================
# F4 — THE SWITCH: count vs eccentricity (CHORUS)
# ======================================================================
def f4():
    A = np.load(f"{REPO}/reports/ch_analysis_recut.npz", allow_pickle=True)
    tags = [str(t) for t in A["surface_tags"]]
    corr = A["surface_corr"]            # re-cut count at the ADOPTED floor
    corr_lo = A["surface_corr_lo"]      # count at floor + bootstrap error -> the CONFIRMED test
    fl = A["surface_floor_adopted"]
    fl_sd = A["surface_floor_adopted_sd"]

    def g(a, tag):
        return float(a[tags.index(tag)])

    ecc = ["m0e00", "m1e03", "m1e05", "m1e07"]           # ONE eccentric member
    lit = [g(corr, f"{t}_lit") for t in ecc]
    vlb = [g(corr, f"{t}_vlbi") for t in ecc]
    lit_lo = [g(corr_lo, f"{t}_lit") for t in ecc]
    vlb_lo = [g(corr_lo, f"{t}_vlbi") for t in ecc]
    flit = [g(fl, f"{t}_lit") for t in ecc]
    flit_sd = [g(fl_sd, f"{t}_lit") for t in ecc]
    fvlb = [g(fl, f"{t}_vlbi") for t in ecc]
    fvlb_sd = [g(fl_sd, f"{t}_vlbi") for t in ecc]

    lever = lit[3] / lit[0]                              # 14.8x, circular -> e=0.7 (lit)

    x = np.arange(4)
    w = 0.36
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(17.5, 7.6),
                                  gridspec_kw={"width_ratios": [1.62, 1]})

    b1 = ax.bar(x - w / 2, lit, w, color=BLUE, label="today's distance knowledge")
    b2 = ax.bar(x + w / 2, vlb, w, color=AMBER, label="sharper (VLBI) distances")
    # whisker down to the count at floor + error: a bar CONFIRMS only if this stays above 1
    for xs, c, c_lo in ((x - w / 2, lit, lit_lo), (x + w / 2, vlb, vlb_lo)):
        ax.vlines(xs, c_lo, c, color=INK, lw=1.6, zorder=4)
        ax.hlines(c_lo, np.asarray(xs) - 0.10, np.asarray(xs) + 0.10, color=INK, lw=1.6,
                  zorder=4)
    ax.axhline(1.0, color=RED, lw=2.6, ls="--", zorder=3)
    ax.text(3.46, 1.18, "the bar:  1 certified pulsar per dataset", color=RED,
            fontsize=13, weight="bold", ha="right")
    for bars in (b1, b2):
        for b in bars:
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.14,
                    f"{b.get_height():.2f}", ha="center", fontsize=13, weight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(["circular\n(e = 0)", "slightly oval\n(e = 0.3)",
                        "oval\n(e = 0.5)", "very oval\n(e = 0.7)"])
    ax.set_xlabel("Shape of ONE of the three loud orbits")
    ax.set_ylabel("Pulsars certified per dataset")
    ax.set_ylim(0, 8.3)
    ax.legend(loc="upper right", framealpha=0.95, fontsize=12)
    ax.set_title("(a) ONE oval orbit turns it on — but only from e = 0.5",
                 loc="left", weight="bold")

    ax.annotate(f"NOT ENOUGH\n{lit[1]:.2f} lit · {vlb[1]:.2f} vlbi\n"
                f"neither holds at the\nbar's own error",
                xy=(1.20, 1.06), xytext=(0.50, 2.35),
                fontsize=12.5, color=RED, weight="bold",
                bbox=dict(boxstyle="round,pad=0.35", fc="#fdeeee", ec=RED, alpha=0.95),
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.8))
    ax.annotate(f"THE SWITCH-ON\n{lit[2]:.2f} — confirmed", xy=(2 - w / 2 + 0.19, lit[2] + 0.10),
                xytext=(2.13, 4.35), fontsize=13, color=GREEN, weight="bold",
                bbox=dict(boxstyle="round,pad=0.35", fc="#eef5ee", ec=GREEN, alpha=0.95),
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.8))
    ax.text(3.0, 6.35, f"{lever:.1f}× the circular count", ha="center",
            fontsize=13, color=GREEN, weight="bold")
    ax.text(0.5, -0.235, "Whisker = the count at the bar's own error bar; a cell CONFIRMS only "
                         "if it stays above 1.\nThat is why e = 0.3 fails in BOTH tiers, though "
                         "VLBI's 1.03 grazes the bar.",
            transform=ax.transAxes, ha="center", va="top", fontsize=11, color=INK,
            style="italic")

    # ---- inset: the SAME e = 0.3 that fails above, with TWO eccentric members ----
    axi = ax.inset_axes([0.035, 0.545, 0.275, 0.375])
    xi = np.arange(2)
    n1 = [g(corr, "m1e03_lit"), g(corr, "m1e03_vlbi")]
    n2 = [g(corr, "m2e03_lit"), g(corr, "m2e03_vlbi")]
    n1_lo = [g(corr_lo, "m1e03_lit"), g(corr_lo, "m1e03_vlbi")]
    n2_lo = [g(corr_lo, "m2e03_lit"), g(corr_lo, "m2e03_vlbi")]
    axi.bar(xi - 0.20, n1, 0.38, color=GREY, label="ONE member")
    axi.bar(xi + 0.20, n2, 0.38, color=GREEN, label="TWO members")
    for xs, c, c_lo in ((xi - 0.20, n1, n1_lo), (xi + 0.20, n2, n2_lo)):
        axi.vlines(xs, c_lo, c, color=INK, lw=1.3, zorder=4)
        axi.hlines(c_lo, np.asarray(xs) - 0.07, np.asarray(xs) + 0.07, color=INK, lw=1.3,
                   zorder=4)
    for xs, c in ((xi - 0.20, n1), (xi + 0.20, n2)):
        for xx, cc in zip(xs, c):
            axi.text(xx, cc + 0.09, f"{cc:.2f}", ha="center", fontsize=10, weight="bold")
    axi.axhline(1.0, color=RED, lw=1.8, ls="--", zorder=3)
    axi.set_xticks(xi)
    axi.set_xticklabels(["today's\ndistances", "VLBI"], fontsize=10)
    axi.set_ylim(0, 3.9)
    axi.set_yticks([0, 1, 2, 3])
    axi.tick_params(labelsize=9)
    axi.set_ylabel("certified", fontsize=10)
    axi.legend(fontsize=9, loc="upper left", framealpha=0.95, handlelength=1.1,
               borderpad=0.3, labelspacing=0.25)
    axi.set_title("at e = 0.3, TWO members\nDO turn it on (both tiers)",
                  fontsize=11, weight="bold", loc="left", color=GREEN)
    for s in ("top", "right"):
        axi.spines[s].set_visible(True)
    axi.patch.set_alpha(0.97)

    # ---- (b) the honesty check, on the CORRECTED floors ----
    ax2.errorbar(x, flit, yerr=flit_sd, marker="o", ms=11, lw=2.6, color=BLUE,
                 capsize=6, label="today's distances")
    ax2.errorbar(x, fvlb, yerr=fvlb_sd, marker="s", ms=11, lw=2.6, color=AMBER,
                 capsize=6, label="sharper distances")
    ax2.set_xticks(x)
    ax2.set_xticklabels(["0", "0.3", "0.5", "0.7"])
    ax2.set_xlabel("Orbit shape (eccentricity)")
    ax2.set_ylabel("Height of the bar\n(evidence noise alone can fake, nats)")
    ax2.set_ylim(0, 16)
    ax2.legend(loc="lower right", framealpha=0.95)
    ax2.set_title("(b) ...and the bar does NOT drop to let it in", loc="left", weight="bold")
    ax2.annotate(f"at the switch-on the bar is\nHIGHER than circular "
                 f"({flit[0]:.1f} → {flit[2]:.1f})\nand the count still rises "
                 f"{lit[2] / lit[0]:.1f}×",
                 xy=(2, flit[2] + flit_sd[2] + 0.25), xytext=(0.06, 13.4), fontsize=12,
                 color=INK,
                 bbox=dict(boxstyle="round,pad=0.35", fc="#f4f4f4", ec=GREY, alpha=0.95),
                 arrowprops=dict(arrowstyle="->", color=INK, lw=1.4))
    ax2.text(0.5, -0.235, "The bar is not monotone in e: it peaks at e = 0.3 — the one\n"
                          "shape that fails. No eccentric bar sits below the circular one.",
             transform=ax2.transAxes, ha="center", va="top", fontsize=10.5, color=INK,
             style="italic")

    fig.suptitle("F4 · THE SWITCH — at exactly the same loudness, orbit shape decides whether "
                 "you measure anything —\nand how MANY orbits carry the shape decides how oval "
                 "they must be",
                 fontsize=19, weight="bold", y=1.05)
    fig.savefig(f"{OUT}/F4_the_switch.png")
    plt.close(fig)
    print(f"F4 ok  (1-member lit {np.round(lit, 2)}, lo {np.round(lit_lo, 2)}; "
          f"2-member e=0.3 lit {n2[0]:.2f} vlbi {n2[1]:.2f}; "
          f"floors lit {np.round(flit, 2)}; lever {lever:.1f}x)")


# ======================================================================
# F5 — THE MAP: SURFACE's onset surface
# ======================================================================
def f5():
    A = np.load(f"{REPO}/reports/surface_analysis_recut.npz", allow_pickle=True)
    cols = [str(c) for c in A["cols"]]
    T = A["table"]
    tiers = np.array([str(t) for t in A["tiers"]])
    ci = {c: i for i, c in enumerate(cols)}
    REFIT = A["estimator"] == "emp_q95"     # the 15 cells RECUT re-fitted off the Gumbel

    H = np.array([-13.25, -13.0, -12.75, -12.5, -12.25, -12.0])
    Ts = np.array([30, 40, 50])

    def grid(k, key, arr=None):
        m = (T[:, ci["k"]] == k) & (tiers == "lit")
        G = np.full((len(Ts), len(H)), np.nan)
        for j in np.flatnonzero(m):
            r = T[j]
            a = list(Ts).index(int(r[ci["T"]]))
            b = list(H).index(round(r[ci["h"]], 2))
            G[a, b] = r[ci[key]] if arr is None else arr[j]
        return G

    fig, axes = plt.subplots(1, 2, figsize=(18, 6.6))
    panels = [(3, "(a) The population we think we have\n(3 loud sources of 16)"),
              (5, "(b) Two more loud sources\n(5 of 16) — the same array"),]
    vmax = 8.0
    n_refit = 0
    for ax, (k, title) in zip(axes, panels):
        C = grid(k, "corr")
        CL = grid(k, "corr_lo")
        RF = grid(k, "corr", arr=REFIT.astype(float))
        im = ax.imshow(C, origin="lower", aspect="auto", cmap="viridis",
                       vmin=0, vmax=vmax,
                       extent=[-0.5, len(H) - 0.5, -0.5, len(Ts) - 0.5])
        # bold onset contour: the confirmed-onset boundary (count survives floor + error)
        ax.contour(np.arange(len(H)), np.arange(len(Ts)), CL, levels=[1.0],
                   colors="white", linewidths=4.5)
        ax.contour(np.arange(len(H)), np.arange(len(Ts)), CL, levels=[1.0],
                   colors=RED, linewidths=2.4)
        # RECUT re-fitted these floors off the Gumbel and onto the empirical q95. They are
        # settled, not pending: a thin dotted border states the provenance, nothing more.
        # (The former "////" hatch meant "floor pending re-fit, count is an UPPER BOUND".
        #  Nothing on this map is an upper bound any more — every floor is resolved.)
        for a in range(len(Ts)):
            for b in range(len(H)):
                if RF[a, b] > 0.5:
                    ax.add_patch(Rectangle((b - 0.5, a - 0.5), 1, 1, fill=False,
                                           edgecolor="white", lw=1.6, ls=(0, (2, 2)),
                                           zorder=4))
                    n_refit += 1
        for a in range(len(Ts)):
            for b in range(len(H)):
                v = C[a, b]
                ax.text(b, a, f"{v:.1f}", ha="center", va="center", fontsize=14.5,
                        weight="bold", color="white" if v < 0.55 * vmax else "black",
                        zorder=6)
        ax.set_xticks(range(len(H)))
        ax.set_xticklabels(["faintest\n$10^{-13.25}$", "$10^{-13}$", "$10^{-12.75}$",
                            "$10^{-12.5}$", "$10^{-12.25}$", "loudest\n$10^{-12}$"],
                           fontsize=12.5)
        ax.set_yticks(range(len(Ts)))
        ax.set_yticklabels([f"{t} yr" for t in Ts])
        ax.set_xlabel("How loud the sources are")
        ax.set_ylabel("How long you watch")
        ax.set_title(title, loc="left", weight="bold", fontsize=17)

    kbox = dict(boxstyle="round,pad=0.35", fc="black", ec="white", alpha=0.55)
    axes[0].annotate("ON — right of the red line",
                     xy=(4.0, 0.30), xytext=(2.72, 0.66), fontsize=13,
                     color="white", weight="bold", bbox=kbox,
                     arrowprops=dict(arrowstyle="->", color="white", lw=2.0))
    axes[0].annotate("OFF", xy=(1.0, 0.30), xytext=(0.30, 0.62), fontsize=13.5,
                     color="white", weight="bold", bbox=kbox,
                     arrowprops=dict(arrowstyle="->", color="white", lw=2.0))
    axes[1].text(2.9, 0.30, "the red line has moved all the way\n"
                            "to the faint edge — almost everything is ON",
                 ha="center", va="center", fontsize=13, color="white",
                 weight="bold", bbox=kbox)

    cb = fig.colorbar(im, ax=axes, fraction=0.032, pad=0.02)
    cb.set_label("Pulsars certified per dataset", fontsize=15)
    cb.ax.tick_params(labelsize=13)

    refit_proxy = Patch(facecolor="#4c4c4c", edgecolor="white", ls=(0, (2, 2)), lw=1.6,
                        label="floor re-fit on the empirical 95th percentile (settled)")
    axes[0].legend(handles=[refit_proxy], loc="upper center",
                   bbox_to_anchor=(1.10, -0.30), fontsize=13, framealpha=0.96)
    fig.text(0.5, -0.30,
             "Dotted border: too many of this cell's noise-only trials certify nothing at all "
             "for a Gumbel tail-fit to be valid (ANCHOR §4),\nso its floor is the empirical 95th "
             "percentile with a bootstrap error. Those floors are now RE-FIT and the counts are "
             "final —\nno cell on this map is an upper bound. The re-fit moved four cells across "
             "the boundary: two onsets died at the faint edge and two were born.",
             ha="center", fontsize=11.5, color=INK,
             bbox=dict(boxstyle="round,pad=0.4", fc="#f4f4f4", ec=GREY, alpha=0.95))

    fig.suptitle("F5 · THE MAP — what it takes to switch the measurement on\n"
                 "red line = the on/off boundary, drawn only where it survives its own error bar",
                 fontsize=19, weight="bold", y=1.06)
    fig.savefig(f"{OUT}/F5_the_map.png")
    plt.close(fig)
    print(f"F5 ok  ({n_refit} lit cells re-fit on emp q95; N_onset(all 108) = "
          f"{int(A['n_onset'])}, lost {int(A['n_lost'])}, gained {int(A['n_gained'])})")


# ======================================================================
# F6 — THE PAYOFF: distance precision vs certified pulsars
# ======================================================================
def f6():
    S = np.load(f"{REPO}/SIREN_results/siren_summary.npz", allow_pickle=True)
    frac = S["B_0p1pc__frac_DL"]        # (f, mc, seedset), fraction
    gf = S["grid_f"]
    gmc = S["grid_mc"]
    DL = S["DL_Mpc"]
    nseed = [0, 1, 2, 3, 5]             # tags 0..4 are the frequency-ranked nested sequence

    M4 = np.load(f"{REPO}/reports/atlas_m2m4_summary.npz")["m4"]
    # cols: forb, mc, estar, sig_mc_comb, fDL_comb, fDL_1, fDL_3
    row = M4[(M4[:, 0] == -8.0) & (M4[:, 1] == 9.5)][0]
    m4_comb = row[4] * 100

    fig, ax = plt.subplots(figsize=(14.5, 8.6))
    fig.subplots_adjust(bottom=0.34)
    box = dict(boxstyle="round,pad=0.4", fc="white", ec=GREY, alpha=0.93)

    ax.axhspan(10, 30, color=GREEN, alpha=0.16, zorder=0)
    ax.text(5.28, 17.0, "USEFUL FOR COSMOLOGY\n(the 10–30% class)", ha="right",
            fontsize=14.5, color=GREEN, weight="bold", zorder=5)

    cmap = plt.get_cmap("plasma")
    for i in range(3):
        for j in range(3):
            y = frac[i, j, :5] * 100
            c = cmap(0.08 + 0.78 * (j * 3 + i) / 8)
            ax.plot(nseed, y, marker="o", ms=8, lw=2.2, color=c, alpha=0.92,
                    label=f"$M_c=10^{{{gmc[j]:.1f}}}M_\\odot$ at {DL[i, j]:.1f} Mpc")

    ax.plot([0], [m4_comb], marker="*", ms=32, color=RED, mec=INK, mew=1.2, ls="none",
            zorder=6, label="an OVAL orbit, zero certified pulsars (it clocks itself)")

    ax.set_yscale("log")
    ax.set_xticks(nseed)
    ax.set_xlabel("Number of pulsars you have certified")
    ax.set_ylabel("Distance uncertainty to the\nblack-hole pair  (%)")
    ax.set_xlim(-0.35, 5.35)
    ax.set_ylim(4, 900)
    ax.set_yticks([5, 10, 30, 100, 300])
    ax.set_yticklabels(["5%", "10%", "30%", "100%", "300%"])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              framealpha=0.96, fontsize=12.5)

    ax.annotate("with none certified, the distance\nis essentially unmeasured (332%)",
                xy=(0.05, 332), xytext=(0.75, 640), fontsize=13.5, color=INK, bbox=box,
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.5))
    ax.annotate("for the heavier pairs, 3–5 certified\npulsars reach the useful band",
                xy=(3.05, 11.5), xytext=(1.30, 42), fontsize=13.5, color=INK,
                weight="bold", bbox=box,
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.5))
    ax.annotate(f"{m4_comb:.0f}% — for free", xy=(0.03, m4_comb), xytext=(0.45, 5.6),
                fontsize=14.5, color=RED, weight="bold", bbox=box,
                arrowprops=dict(arrowstyle="->", color=RED, lw=1.8))

    ax.set_title("F6 · THE PAYOFF — what one certified pulsar buys\n"
                 "(each line is one black-hole pair; a single loud source at $h=10^{-13.25}$)",
                 loc="left", weight="bold", fontsize=19)
    fig.savefig(f"{OUT}/F6_the_payoff.png")
    plt.close(fig)
    print(f"F6 ok  (N0 {frac[1,1,0]*100:.0f}%, N3 {frac[1,1,3]*100:.1f}%, "
          f"ATLAS comb-only {m4_comb:.1f}%)")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    f1(); f2(); f3(); f4(); f5(); f6()
    print("\nall six written to", OUT)
