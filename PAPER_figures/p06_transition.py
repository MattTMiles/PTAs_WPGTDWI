"""P6 — THE TRANSITION (prong 2).

(a) marginal/conditional distance information vs N_CW with the confusion knee,
    and the array-scaling inset (knee/N* = 0.40 N_psr^1.03 vs Boyle & Pen 2N/7);
(b) the coherence law — marg/cond collapse under phase wander for both banked
    functional forms; the 0.5-crossing location is the invariant (0.30x shift
    between forms; absolute abscissa values are not physical);
(c) the compounding map — confusion x coherence residual vs the factorised
    product (commuting reduction), two-sided structure, peak |residual| 0.202.

Sources: CW_transition/prong2_results.npz (tot_* curve + knee);
CW_transition/stagec_p2c_results.npz (A, b fit, 5 array sizes);
CW_transition/prong2_coherence.npz (exp1_linear/saturating, x_half_*);
CW_transition/prong2_Ntc_map.npz (resid_commuting map, x_grid=N/N*, Y_grid).
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"


def main():
    P = np.load(f"{REPO}/CW_transition/prong2_results.npz", allow_pickle=True)
    n = P["tot_n"]
    ratio = P["tot_marg"] / P["tot_cond"]
    rlo = P["tot_marg_lo"] / P["tot_cond"]
    rhi = P["tot_marg_hi"] / P["tot_cond"]
    knee = float(P["diag_T15_b3-20_knee"])

    S = np.load(f"{REPO}/CW_transition/stagec_p2c_results.npz")
    A_, b_ = float(S["A"]), float(S["b"])
    npsr, kk = S["n_psr"], S["knee_over_Nstar"]

    C = np.load(f"{REPO}/CW_transition/prong2_coherence.npz", allow_pickle=True)
    xg = C["x_grid"]
    lin, sat = C["exp1_linear"], C["exp1_saturating"]
    xh_l, xh_s = float(C["x_half_linear"]), float(C["x_half_saturating"])

    M = np.load(f"{REPO}/CW_transition/prong2_Ntc_map.npz", allow_pickle=True)
    xN, Y = M["x_grid"], M["Y_grid"]
    RES = M["resid_commuting"]
    peak = float(M["max_factor_dev_commuting"])

    fig, (ax, axc, axr) = plt.subplots(
        1, 3, figsize=(st.DOUBLE, 2.35), gridspec_kw={"width_ratios": [1.2, 1, 1.15]})
    fig.subplots_adjust(wspace=0.42)

    # ---------- (a) confusion ----------
    ax.fill_between(n, rlo, rhi, color=st.BLUE, alpha=0.18, lw=0)
    ax.plot(n, ratio, "-o", ms=2.4, lw=1.0, color=st.BLUE)
    ax.axvline(knee, color=st.VERMIL, lw=0.9, ls="--")
    ax.axhline(0.5, color=st.GREY, lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(1, 1100)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("number of CW sources $N$")
    ax.set_ylabel("marginal / conditional information")
    st.panel_label(ax, "(a)")
    ax.annotate(f"confusion knee\n$N\\approx{knee:.0f}$", xy=(knee, 0.60),
                xytext=(0.60, 0.80), textcoords="axes fraction", fontsize=5.8,
                color=st.VERMIL,
                arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.6))
    axi = ax.inset_axes([0.14, 0.07, 0.42, 0.36])
    xx = np.linspace(12, 230, 100)
    axi.plot(xx, A_ * xx ** b_, lw=0.9, color=st.BLUE)
    axi.plot(xx, (2 / 7) * xx, lw=0.8, ls="--", color=st.GREY)
    axi.plot(npsr, kk, "o", ms=2.2, color=st.VERMIL)
    axi.set_xscale("log"); axi.set_yscale("log")
    axi.set_xticks([15, 60, 200]); axi.set_xticklabels(["15", "60", "200"],
                                                       fontsize=4.5)
    axi.set_yticks([10, 100]); axi.set_yticklabels(["10", "100"], fontsize=4.5)
    axi.tick_params(length=2)
    axi.set_xlabel("$N_{\\rm psr}$", fontsize=5, labelpad=0.0)
    axi.set_ylabel("knee/$N^*$", fontsize=5, labelpad=0.0)
    axi.text(0.97, 0.10, f"${A_:.2f}\\,N^{{{b_:.2f}}}$ ($\\approx 2N/5$)\n"
                         "vs Boyle–Pen $2N/7$ (dashed)",
             transform=axi.transAxes, fontsize=4.6, va="bottom", ha="right")
    for s in ("top", "right"):
        axi.spines[s].set_visible(True)

    # ---------- (b) coherence ----------
    axc.plot(xg, lin, lw=1.0, color=st.BLUE, label="linear-wander form")
    axc.plot(xg, sat, lw=1.0, ls="--", color=st.ORANGE, label="saturating form")
    axc.axhline(0.5, color=st.GREY, lw=0.6, ls=":")
    for xh, col in ((xh_l, st.BLUE), (xh_s, st.ORANGE)):
        axc.axvline(xh, color=col, lw=0.6, ls=":", alpha=0.8)
    axc.set_xscale("log")
    axc.set_xlim(xg.min(), xg.max())
    axc.set_ylim(0, 1.0)
    axc.set_xlabel("phase-wander abscissa (scaled; see note)")
    axc.set_ylabel("marginal / conditional information")
    st.panel_label(axc, "(b)")
    axc.legend(loc="upper right", fontsize=5.4)
    axc.annotate(f"0.5-crossings differ by\n$0.30\\times$ between forms —\n"
                 "the LOCATION ($\\mathrm{SNR}^2\\sigma_\\phi^2\\approx1$)\n"
                 "is the invariant",
                 xy=(xh_l, 0.5), xytext=(0.30, 0.62), textcoords="axes fraction",
                 fontsize=5.4,
                 arrowprops=dict(arrowstyle="->", color=st.INK, lw=0.5))
    axc.text(0.02, 0.03, "beyond the knee: $\\sigma_L\\to d_L\\sigma_\\phi/2\\pi$,\n"
             "an SNR-INDEPENDENT wander floor;\nabsolute abscissa not physical",
             transform=axc.transAxes, fontsize=5.2, style="italic")

    # ---------- (c) compounding ----------
    vmax = 0.25
    im = axr.pcolormesh(xN, Y, RES.T, cmap=st.DIV_CMAP, vmin=-vmax, vmax=vmax,
                        shading="nearest")
    axr.set_xscale("log"); axr.set_yscale("log")
    axr.contour(xN, Y, M["R"].T, levels=[0.5], colors=[st.INK], linewidths=0.9)
    axr.set_xlabel("$N/N^*$ (confusion axis)")
    axr.set_ylabel("coherence axis $\\mathrm{SNR}^2\\sigma_\\phi^2$")
    st.panel_label(axr, "(c)")
    cb = fig.colorbar(im, ax=axr, fraction=0.05, pad=0.02)
    cb.set_label("$R - R_{\\rm conf}R_{\\rm coh}$", fontsize=6)
    cb.ax.tick_params(labelsize=5)
    axr.text(0.03, 0.05,
             f"peak $|$resid$|$ = {abs(peak):.3f}\nline: $R=0.5$ (curved —\n"
             "the transitions COMPOUND)",
             transform=axr.transAxes, fontsize=5.4, color=st.INK)

    st.savefig(fig, "P6_the_transition")
    print(f"P6: ratio at N=1/100/398/1000 = "
          f"{[round(float(ratio[list(n).index(v)]), 3) for v in (1, 100, 398, 1000)]}, "
          f"knee {knee:.1f}, A={A_:.4f}, b={b_:.4f}, peak resid {peak:.4f}")


if __name__ == "__main__":
    main()
