"""P1 — THE PROBLEM.

(a) one pulsar's fringe comb (likelihood vs distance) with the EM prior window
    shaded;
(b) the census K distribution — all 116 pulsars' wrong-fringe counts inside
    their real priors, J0437 marked, "no free anchors" stated.

Sources: legacy banked lnL(D) scan (the only banked likelihood-vs-distance
curve in the repo — its comb is COARSER than the census's, stated on panel);
CW_transition/anchor_a0_priors.npz; reports/geo_dlnl_bank.npz (K_counted,
median over 40 sky draws).
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"


def main():
    scan = np.load(f"{REPO}/lnLs_GWAmp_phase_connected/runD_3CW_test/"
                   f"J0437-4715_cw_p_dist/psrTerm/"
                   f"runD_3CW_test_J0437-4715_cw_p_dist_psrTerm_0_0.npz",
                   allow_pickle=True)
    D = scan["scan_values"]
    y = scan["logls"][0]
    y = y - y.max()

    a0 = np.load(f"{REPO}/CW_transition/anchor_a0_priors.npz", allow_pickle=True)
    i = list(a0["fname"]).index("J0437-4715")
    D_em = float(a0["dist_feather_kpc"][i])
    s_em = float(a0["sig_feather_kpc"][i])

    geo = np.load(f"{REPO}/reports/geo_dlnl_bank.npz", allow_pickle=True)
    names = list(geo["names"])
    Kmed = np.median(geo["K_counted"], axis=0)
    K_j0437 = Kmed[names.index("J0437-4715")]
    n_anchor = int((Kmed <= 1).sum())      # strict anchors: K <= 1

    fig, (ax, axz, axk) = plt.subplots(
        1, 3, figsize=(st.DOUBLE, 2.15), gridspec_kw={"width_ratios": [1.15, 1, 1]})
    fig.subplots_adjust(wspace=0.42)

    # (a) the comb
    ax.plot(D, y, lw=0.25, color=st.BLUE, alpha=0.9)
    ax.axvline(D_em, color=st.VERMIL, lw=1.0, zorder=5)
    ax.set_xlim(0.1, 5.0)
    ax.set_ylim(-2350, 700)
    ax.set_xlabel("pulsar distance $L$ (kpc)")
    ax.set_ylabel("$\\Delta\\ln\\mathcal{L}$")
    st.panel_label(ax, "(a)")
    ax.annotate("true $L$", xy=(D_em, 300), xytext=(1.3, 420), fontsize=6.5,
                color=st.VERMIL,
                arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.7))
    ax.text(0.97, 0.05, "legacy banked scan;\ncomb coarser than census",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=5.2,
            style="italic", color=st.GREY)

    # (a, zoom) the EM prior window
    lo, hi = 0.100, 0.225
    m = (D >= lo) & (D <= hi)
    axz.plot(D[m], y[m], lw=0.7, color=st.BLUE)
    axz.axvspan(D_em - s_em, D_em + s_em, color=st.VERMIL, alpha=0.28, lw=0,
                label="EM prior ($\\pm 1\\sigma$)")
    axz.axvline(D_em, color=st.VERMIL, lw=0.8)
    axz.set_xlim(lo, hi)
    axz.set_xlabel("pulsar distance $L$ (kpc)")
    axz.set_ylabel("$\\Delta\\ln\\mathcal{L}$")
    st.panel_label(axz, "(b)")
    axz.legend(loc="lower left", fontsize=5.8)
    axz.set_title("which fringe?", fontsize=7, style="italic", pad=2)

    # (c) census K distribution
    axk.hist(np.log10(np.maximum(Kmed, 1)), bins=26, color=st.BLUE, alpha=0.85,
             edgecolor="white", lw=0.4)
    axk.axvline(np.log10(K_j0437), color=st.GREEN, lw=1.2, zorder=5)
    axk.axvline(np.log10(np.median(Kmed)), color=st.VERMIL, lw=1.2, ls="--",
                zorder=5)
    axk.set_xticks([0, 1, 2, 3, 4, 5])
    axk.set_xticklabels(["$1$", "$10$", "$10^2$", "$10^3$", "$10^4$", "$10^5$"])
    axk.set_xlim(-0.3, 5.3)
    axk.set_xlabel("fringes inside EM prior, $K$")
    axk.set_ylabel("pulsars")
    st.panel_label(axk, "(c)")
    axk.annotate(f"J0437$-$4715\n$K={K_j0437:.0f}$",
                 xy=(np.log10(K_j0437), 12), xytext=(0.04, 0.62),
                 textcoords="axes fraction", fontsize=6, color=st.GREEN,
                 arrowprops=dict(arrowstyle="->", color=st.GREEN, lw=0.7))
    axk.annotate(f"median\n$K={np.median(Kmed):,.0f}$",
                 xy=(np.log10(np.median(Kmed)) - 0.06, 17), xytext=(0.42, 0.86),
                 textcoords="axes fraction", fontsize=6, color=st.VERMIL,
                 ha="right",
                 arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.7))
    axk.text(0.03, 0.04, f"$K\\leq 1$ (free anchors): {n_anchor}/116",
             transform=axk.transAxes, fontsize=6, weight="bold")

    st.savefig(fig, "P1_the_problem")
    print(f"P1: K_J0437={K_j0437:.0f}, median K={np.median(Kmed):.0f}, "
          f"anchors K<=1: {n_anchor}")


if __name__ == "__main__":
    main()
