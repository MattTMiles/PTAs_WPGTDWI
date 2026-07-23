"""P3 — THE HONEST CRITERION.

(a) null vs signal offender distributions at a representative cell
    (h=-13.00, T=30, lit, census structure), old max-of-10 bar vs the honest
    adopted bar, the retraction shown;
(b) floor validity — per-cell null zero-fraction with the 20% Gumbel gate, and
    the Gumbel-vs-empirical divergence where the gate fails (errs PERMISSIVE);
(c) ANCHOR's realism ladder — per-cell floor shifts consistent with zero while
    the realised noise rms rises 32%.

Sources: SURFACE_results/sf_{nullN,sig}_h1300_T30_lit_k3_*.npz (re-cut, CPU);
reports/surface_analysis_recut.npz; reports/anchor_ladder.npz. The re-cut
scorer is gated against the banked re-cut count before drawing. The realised
per-rung rms values (1.993->2.634 us) are quoted from the tracked ANCHOR
report table (ANCHOR_realdata_null.md); they are not npz-banked.
"""
import glob

import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"
QBAR = 0.9
OLD_FLOOR = 13.48   # IGNITE-family max-of-10 floor for this cell (banked, SURFACE D-4 context)
RMS_US = {"R0": 1.993, "R1": 2.030, "R2": 2.207, "R3": 2.634}  # ANCHOR report table (prose)


def offender(dlnL, lnK, qmax):
    m = (dlnL > lnK) & (qmax > QBAR)
    return float(dlnL[m].max()) if m.any() else 0.0


def load_cell(pat):
    rows = []
    for f in sorted(glob.glob(pat)):
        d = np.load(f, allow_pickle=True)
        rows.append((np.nan_to_num(d["dlnL_det"], posinf=1e30),
                     d["lnK"], d["qmax"], d["on_true"]))
    return rows


def main():
    # ---------------- panel (a): the representative cell ----------------
    tag = "h1300_T30_lit_k3"
    N = load_cell(f"{REPO}/SURFACE_results/sf_nullN_{tag}_g*_n*.npz")
    S = load_cell(f"{REPO}/SURFACE_results/sf_sig_{tag}_g*_n*.npz")
    off_null = np.array([offender(*r[:3]) for r in N])
    off_sig = np.array([offender(*r[:3]) for r in S])

    R = np.load(f"{REPO}/reports/surface_analysis_recut.npz", allow_pickle=True)
    ci = {c: i for i, c in enumerate(str(c) for c in R["cols"])}
    T = R["table"]
    tiers = np.array([str(t) for t in R["tiers"]])
    struct = np.array([str(s) for s in R["struct"]])
    m = ((np.round(T[:, ci["h"]], 2) == -13.00) & (T[:, ci["T"]] == 30)
         & (tiers == "lit") & (struct == "3+13"))
    r = T[int(np.flatnonzero(m)[0])]
    floor, sd = r[ci["floor_adopted"]], r[ci["floor_adopted_sd"]]
    gum, zf = r[ci["fN"]], r[ci["fN_zerofrac"]]
    c_new = r[ci["corr"]]

    def count(fl):
        c = 0
        for d, k, q, ot in S:
            c += int(((d > np.maximum(k, fl)) & (q > QBAR) & ot).sum())
        return c / len(S)

    c_old = count(OLD_FLOOR)
    assert abs(count(floor) - c_new) < 1e-9, "P3a scorer disagrees with recut bank"

    fig, (ax, axv, axl) = plt.subplots(
        1, 3, figsize=(st.DOUBLE, 2.25), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    fig.subplots_adjust(wspace=0.40)

    bins = np.linspace(0, 46, 32)
    ax.hist(off_null, bins=bins, color=st.GREY, alpha=0.85, edgecolor="white",
            lw=0.3, label=f"pure noise ($N={len(off_null)}$)")
    ax.hist(off_sig, bins=bins, color=st.BLUE, alpha=0.75, edgecolor="white",
            lw=0.3, label=f"signal ($N={len(off_sig)}$)")
    ax.axvspan(floor - sd, floor + sd, color=st.GREEN, alpha=0.20, zorder=1)
    ax.axvline(floor, color=st.GREEN, lw=1.2, zorder=4)
    ax.axvline(OLD_FLOOR, color=st.VERMIL, lw=1.2, ls="--", zorder=4)
    ax.axvline(gum, color=st.GREY, lw=0.8, ls=":", zorder=3)
    ax.set_xlim(0, 46)
    ax.set_xlabel("offender statistic (nat)")
    ax.set_ylabel("realisations")
    st.panel_label(ax, "(a)")
    ax.legend(loc="upper right", fontsize=5.6)
    ax.annotate(f"old bar (max-of-10):\n{OLD_FLOOR:.1f} $\\to$ {c_old:.2f}/real",
                xy=(OLD_FLOOR, 16), xytext=(0.02, 0.55),
                textcoords="axes fraction", fontsize=5.6, color=st.VERMIL,
                weight="bold",
                arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.6))
    ax.annotate(f"honest bar (emp. q95):\n{floor:.1f}$\\pm${sd:.1f} $\\to$ "
                f"{c_new:.2f}/real\nRETRACTED",
                xy=(floor, 8), xytext=(0.50, 0.44),
                textcoords="axes fraction", fontsize=5.6, color=st.GREEN,
                weight="bold",
                arrowprops=dict(arrowstyle="->", color=st.GREEN, lw=0.6))
    ax.annotate(f"invalid Gumbel ({gum:.1f});\nnull {zf:.0%} silent",
                xy=(gum, 21), xytext=(0.55, 0.74), textcoords="axes fraction",
                fontsize=5.2, color=st.GREY,
                arrowprops=dict(arrowstyle="->", color=st.GREY, lw=0.5))

    # ---------------- panel (b): floor validity across the box ----------------
    zfs = T[:, ci["fN_zerofrac"]]
    rat = T[:, ci["fN"]] / T[:, ci["fN_emp"]]      # Gumbel / empirical q95
    for tr, col, mk in (("lit", st.BLUE, "o"), ("vlbi", st.ORANGE, "s")):
        mm = tiers == tr
        axv.plot(zfs[mm], rat[mm], mk, ms=2.4, mfc=col, mec="none", alpha=0.75,
                 label=f"SURFACE, {tr}")
    An = np.load(f"{REPO}/reports/anchor_ladder.npz", allow_pickle=True)
    amask = An["emp_q95"] > 0
    axv.plot(An["zero_frac"][amask], An["gumbel"][amask] / An["emp_q95"][amask],
             "^", ms=2.6, mfc=st.PURPLE, mec="none", alpha=0.75,
             label="ANCHOR ladder (T=15)")
    imin = np.argmin(An["gumbel"][amask] / An["emp_q95"][amask])
    rmin = (An["gumbel"][amask] / An["emp_q95"][amask])[imin]
    axv.annotate(f"$\\times${1/rmin:.1f} too permissive\n(zero-fraction "
                 f"{An['zero_frac'][amask][imin]:.0%})",
                 xy=(An["zero_frac"][amask][imin], rmin),
                 xytext=(0.45, 0.30), textcoords="axes fraction", fontsize=5.4,
                 color=st.VERMIL, weight="bold",
                 arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.6))
    axv.axvspan(0.20, 1.0, color=st.VERMIL, alpha=0.08, zorder=0)
    axv.axvline(0.20, color=st.VERMIL, lw=0.9)
    axv.axhline(1.0, color=st.INK, lw=0.6)
    axv.set_yscale("log")
    axv.set_ylim(0.28, 2.0)
    axv.set_yticks([0.33, 0.5, 1.0, 2.0])
    axv.set_yticklabels(["1/3", "1/2", "1", "2"])
    axv.set_xlabel("null zero-fraction")
    axv.set_ylabel("Gumbel floor / empirical q95")
    st.panel_label(axv, "(b)")
    axv.legend(loc="upper right", fontsize=4.8)
    axv.text(0.02, 0.04, "Gumbel valid only left of the gate;\nabove it the fit "
             "errs PERMISSIVE\n(criterion-v2.2 $\\to$ empirical q95)",
             transform=axv.transAxes, fontsize=5.2, va="bottom")
    n_refit = int((np.array([str(e) for e in R["estimator"]]) == "emp_q95").sum())

    # ---------------- panel (c): the ANCHOR realism ladder ----------------
    A = np.load(f"{REPO}/reports/anchor_ladder.npz", allow_pickle=True)
    rung = np.array([str(x) for x in A["rung"]])
    hh = A["h"]; tier = np.array([str(x) for x in A["tier"]])
    demp = A["d_emp"]          # empirical-q95 shift vs the R0 control, per cell
    mwp = A["mw_p"]
    order = ["R1", "R2", "R3"]
    xpos = {k: i for i, k in enumerate(order)}
    cells = sorted({(h, t) for h, t in zip(hh, tier)})
    rng = np.random.default_rng(5)
    nsig = 0
    for (h, t) in cells:
        xs, ys, ps = [], [], []
        for rg in order:
            mm = (rung == rg) & (hh == h) & (tier == t)
            if mm.any():
                xs.append(xpos[rg]); ys.append(float(demp[mm][0]))
                ps.append(float(mwp[mm][0]))
        col = st.BLUE if t == "lit" else st.ORANGE
        j = rng.normal(0, 0.05)
        axl.plot(np.array(xs) + j, ys, "-o", ms=2.0, lw=0.6, color=col, alpha=0.6)
        for x0, y0, p0 in zip(xs, ys, ps):
            if p0 < 0.05 / 40:
                axl.plot([x0 + j], [y0], "o", ms=5, mfc="none", mec=st.VERMIL,
                         mew=0.9, zorder=6)
                nsig += 1
    axl.axhline(0, color=st.INK, lw=0.6)
    axl.set_xticks(range(3))
    axl.set_xticklabels(["R1\nRN/GWB\nmis-spec", "R2\n+unmod.\nDM GP",
                         "R3\n+non-Gauss.\ntails"], fontsize=5.4)
    axl.set_ylabel("$\\Delta$q95 floor vs control (nat)")
    st.panel_label(axl, "(c)")
    axl.plot([], [], "-o", ms=2, lw=0.6, color=st.BLUE, label="lit cells")
    axl.plot([], [], "-o", ms=2, lw=0.6, color=st.ORANGE, label="vlbi cells")
    axl.plot([], [], "o", ms=5, mfc="none", mec=st.VERMIL, mew=0.9,
             label="signif. (Bonferroni)")
    axl.legend(loc="upper left", fontsize=5.2)
    axl.text(0.97, 0.03,
             "realised noise rms rises\n$1.99\\to2.63\\,\\mu$s ($\\times$1.32)$^{*}$\n"
             "floors: flat (median $-0.18$ nat)",
             transform=axl.transAxes, ha="right", fontsize=5.6, weight="bold")

    st.savefig(fig, "P3_the_honest_criterion")
    print(f"P3: cell floor {floor:.2f}+/-{sd:.2f} (zf {zf:.2f}), counts "
          f"{c_old:.2f}->{c_new:.2f}; {n_refit} cells emp_q95; "
          f"anchor signif cells circled: {nsig}")


if __name__ == "__main__":
    main()
