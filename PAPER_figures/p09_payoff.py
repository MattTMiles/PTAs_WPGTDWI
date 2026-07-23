"""P9 — THE PAYOFF.

(a) sigma(log10 Mc) and sigma(DL)/DL vs N_cert (SIREN arm B, 9 source cells),
    dark-siren band 10-30%, the eccentric self-siren star (ATLAS M4);
(b) sky area vs N_cert with the host-galaxy-unique line;
(c) the growth-path ladder with the science thresholds (3 useful, 5
    saturated) — every rung read back from a bank.

Sources: SIREN_results/siren_summary.npz (arm B_0p1pc: Cramer-Rao bounds on
zero-noise Asimov data with fringe integers GIVEN — every sigma is a LOWER
bound); reports/atlas_m2m4_summary.npz (M4 self-siren);
reports/ch_analysis_recut.npz (mixture counts, criterion-v2.2 floors);
reports/surface_analysis_recut.npz (best cell in the box);
reports/geo_summary.npz (GEO strict counts + union);
CW_transition/stagec_p2b_results.npz (readable sub-array: SNR_pterm > 1).
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"
NSEED = [0, 1, 2, 3, 5]


def main():
    S = np.load(f"{REPO}/SIREN_results/siren_summary.npz", allow_pickle=True)
    frac = S["B_0p1pc__frac_DL"]
    smc = S["B_0p1pc__sig_log10_mc"]
    sky = S["B_0p1pc__sky_deg2"]
    gmc = S["grid_mc"]
    DL = S["DL_Mpc"]

    M4 = np.load(f"{REPO}/reports/atlas_m2m4_summary.npz")["m4"]
    row = M4[(M4[:, 0] == -8.0) & (M4[:, 1] == 9.5)][0]
    m4_pct = row[4] * 100

    fig = plt.figure(figsize=(st.DOUBLE, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 0.95, 1.35], wspace=0.45)
    ax = fig.add_subplot(gs[0])
    axs = fig.add_subplot(gs[1])
    axg = fig.add_subplot(gs[2])

    cmap = plt.get_cmap("viridis")
    # ---------- (a) sigma(DL)/DL, with sigma(Mc) inset ----------
    ax.axhspan(10, 30, color=st.GREEN, alpha=0.15, zorder=0)
    ax.text(5.15, 16.5, "dark-siren\nuseful", ha="right", fontsize=5.4,
            color=st.GREEN, weight="bold")
    for i in range(3):
        for j in range(3):
            c = cmap(0.05 + 0.85 * (j * 3 + i) / 8)
            ax.plot(NSEED, frac[i, j, :5] * 100, "-o", ms=2.0, lw=0.8, color=c,
                    alpha=0.9)
    ax.plot([0], [m4_pct], "*", ms=11, color=st.VERMIL, mec=st.INK, mew=0.5,
            zorder=6)
    ax.annotate(f"eccentric self-siren:\n{m4_pct:.0f}% with 0 seeds",
                xy=(0.06, m4_pct), xytext=(0.30, 0.10),
                textcoords="axes fraction", fontsize=5.4, color=st.VERMIL,
                weight="bold",
                arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.6))
    ax.set_yscale("log")
    ax.set_ylim(4, 900)
    ax.set_yticks([5, 10, 30, 100, 300])
    ax.set_yticklabels(["5", "10", "30", "100", "300"])
    ax.set_xticks(NSEED)
    ax.set_xlabel("certified pulsar terms $N_{\\rm cert}$")
    ax.set_ylabel("$\\sigma(D_L)/D_L$ (%)")
    st.panel_label(ax, "(a)")
    axi = ax.inset_axes([0.58, 0.60, 0.40, 0.37])
    axi.patch.set_alpha(1.0)
    axi.set_zorder(6)
    for i in range(3):
        for j in range(3):
            c = cmap(0.05 + 0.85 * (j * 3 + i) / 8)
            axi.plot(NSEED, smc[i, j, :5], "-", lw=0.6, color=c, alpha=0.9)
    axi.set_yscale("log")
    axi.tick_params(labelsize=4.4, length=2)
    axi.set_ylabel("$\\sigma(\\log M_c)$ (dex)", fontsize=4.8, labelpad=1)
    axi.set_xticks(NSEED)
    for sp in ("top", "right"):
        axi.spines[sp].set_visible(True)

    # ---------- (b) sky area ----------
    for i in range(3):
        for j in range(3):
            c = cmap(0.05 + 0.85 * (j * 3 + i) / 8)
            axs.plot(NSEED, sky[i, j, :5], "-o", ms=2.0, lw=0.8, color=c,
                     alpha=0.9,
                     label=(f"$M_c=10^{{{gmc[j]:.1f}}}$, {DL[i, j]:.0f} Mpc"
                            if i == 1 else None))
    axs.axhline(1e-3, color=st.GREEN, lw=1.0, ls="--")
    axs.text(0.06, 1.35e-3, "host-galaxy-unique", ha="left", fontsize=5.4,
             color=st.GREEN, weight="bold",
             transform=axs.get_yaxis_transform())
    axs.set_yscale("log")
    axs.set_xticks(NSEED)
    axs.set_ylim(2e-4, 12)
    axs.set_xlabel("certified pulsar terms $N_{\\rm cert}$")
    axs.set_ylabel("90% sky area (deg$^2$)")
    st.panel_label(axs, "(b)")
    axs.legend(loc="upper right", fontsize=4.4, handlelength=1.4)
    axs.text(0.5, -0.34,
             "Cramér–Rao bounds, zero-noise Asimov, fringe integers GIVEN — "
             "every $\\sigma$ is a LOWER bound (SIREN convention)",
             transform=axs.transAxes, ha="center", va="top", fontsize=4.8,
             style="italic")

    # ---------- (c) the growth path ----------
    C = np.load(f"{REPO}/reports/ch_analysis_recut.npz", allow_pickle=True)
    tags = [str(t) for t in C["surface_tags"]]
    corr = C["surface_corr"]
    A = np.load(f"{REPO}/reports/surface_analysis_recut.npz", allow_pickle=True)
    aci = {str(c): i for i, c in enumerate(A["cols"])}
    best = float(A["table"][:, aci["corr"]].max())
    G = np.load(f"{REPO}/reports/geo_summary.npz", allow_pickle=True)
    strict_mean = float(np.mean(G["counts_strict"]))
    union = int((G["freq"] > 0).sum())
    B = np.load(f"{REPO}/CW_transition/stagec_p2b_results.npz")
    readable = int((B["snr_pterm"] > 1).sum())

    marks = [
        (corr[tags.index("m1e03_lit")], "single $e$=0.3 member\n(below the bar)",
         st.VERMIL, -1),
        (strict_mean, "GEO strict count\n(mean of 40 skies)", st.GREY, +1),
        (best, "best cell in the\n108-cell box", st.BLUE, +2),
        (union, "certify in $\\geq$1 sky\n(GEO union)", st.INK, -2),
        (readable, "readable sub-array\n(SNR$_{\\rm pterm}$>1)", st.PURPLE, +1),
    ]
    axg.set_xscale("log")
    axg.set_xlim(0.5, 300)
    axg.set_ylim(-1.8, 1.9)
    axg.plot([0.55, 116], [0, 0], lw=3.2, color=st.LIGHT,
             solid_capstyle="round", zorder=1)
    axg.axvspan(3, 5, color=st.GREEN, alpha=0.15, zorder=0)
    axg.text(3.9, 1.62, "3 useful / 5 saturated", ha="center", fontsize=5.4,
             color=st.GREEN, weight="bold")
    for v, lab, c, lvl in marks:
        y = 0.42 * lvl
        num = f"{v:.2f}" if v < 10 else f"{v:.0f}"
        axg.plot([v, v], [0, y * 0.72], lw=0.8, color=c, zorder=3)
        axg.plot([v], [0], "o", ms=5, color=c, mec="white", mew=0.8, zorder=5)
        axg.text(v, y * 0.85, f"{num}\n{lab}", ha="center",
                 va=("bottom" if lvl > 0 else "top"), fontsize=4.9, color=c,
                 weight="bold")
    axg.plot([116], [0], "o", ms=5, color=st.GREY, mec="white", mew=0.8,
             zorder=5)
    axg.text(116, -0.30, "116\nfull array\n(never reached)", ha="center",
             va="top", fontsize=4.9, color=st.GREY, weight="bold")
    axg.set_yticks([])
    axg.spines["left"].set_visible(False)
    axg.set_xticks([1, 3, 5, 10, 18, 30, 70, 116])
    axg.set_xticklabels(["1", "3", "5", "10", "18", "30", "70", "116"])
    axg.set_xlabel("certified pulsars (counts / pool sizes)")
    st.panel_label(axg, "(c)")

    st.savefig(fig, "P9_the_payoff")
    print(f"P9: m4 {m4_pct:.1f}%, ladder = 0.70:{corr[tags.index('m1e03_lit')]:.2f} "
          f"strict:{strict_mean:.2f} best:{best:.2f} union:{union} "
          f"readable:{readable}")


if __name__ == "__main__":
    main()
