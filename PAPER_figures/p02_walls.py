"""P2 — THE WALLS.

(a) the F2 ladder: all 348 pulsar-source sky registration tolerances vs the
    measured blind-float band, the 27-113x gap shaded, zero rungs inside;
(b) the clocks: baseline T for a blind search to reach the first rung, under
    BOTH banked SNR scalings (p=0.5 baseline ~11,000 yr; p=1.5 optimistic);
(c) the tier referendum + the >20x chirp-mass-wall bound.

Sources: CW_transition/trackB_F2_ladder.npz (tol_sky, 348 pairs);
CW_transition/PENCIL_t_crossing.npz (all crossing scalars, both scalings);
(c) STORY.md S4.2.8/S4.2.10 canonical values — their npz (b1_step2_table.npz,
b1_breakeven_curve.npz) were produced on cronus and are NOT in this checkout
(repo-wide find is empty); the on-disk b1_referendum_tier{A,C}.npz are the
SUPERSEDED 2-seed runs and are not plotted. Stated on the panel.
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"

BLIND_LO, BLIND_HI = 0.05, 0.21   # blind-float band: seed floor -> endpoint (banked Track-B close-out)
PULLIN = 1e-4                     # fixed-integer pull-in basin (L2b/L2c, banked)

# STORY S4.2.8 (frozen 4-seed protocol PRIMARY per adjudicated D-2) + S4.2.10 curve
REF_F = [0.0847, 0.0481, 0.0323]
REF_ERR = [0.131, 0.0227, 0.0134]
REF_TRAIL_C = (0.0431, 0.0185)          # 5-seed auto-ingest, retained in the trail
BE_LAM = np.array([1.0, 0.3, 0.12, 0.05])
BE_F = np.array([0.032, 0.107, 0.289, 0.673])


def main():
    F = np.load(f"{REPO}/CW_transition/trackB_F2_ladder.npz", allow_pickle=True)
    tol = np.asarray(F["tol_sky"], float)
    tol = tol[np.isfinite(tol)]
    LOOSEST = tol.max()
    n_gap = int(((tol > LOOSEST) & (tol < BLIND_LO)).sum())
    assert len(tol) == 348 and n_gap == 0

    P = np.load(f"{REPO}/CW_transition/PENCIL_t_crossing.npz", allow_pickle=True)
    Tl05, Tl15 = float(P["Tcross_loose_p05"]), float(P["Tcross_loose_p15"])
    Tc05, Tc15 = float(P["Tcross_l2c_p05"]), float(P["Tcross_l2c_p15"])
    T0 = float(P["T0"])

    fig, (ax, axb, axc) = plt.subplots(
        1, 3, figsize=(st.DOUBLE, 2.25), gridspec_kw={"width_ratios": [1.35, 0.85, 1.0]})
    fig.subplots_adjust(wspace=0.40)

    # ---- (a) the ladder
    ax.axvspan(np.log10(LOOSEST), np.log10(BLIND_LO), color=st.LIGHT, alpha=0.55,
               zorder=0)
    ax.axvspan(np.log10(BLIND_LO), np.log10(BLIND_HI), color=st.VERMIL, alpha=0.18,
               zorder=0)
    ax.hist(np.log10(tol), bins=34, color=st.BLUE, alpha=0.95, edgecolor="white",
            lw=0.3, zorder=3)
    ax.axvline(np.log10(LOOSEST), color=st.GREEN, lw=1.0, zorder=4)
    ax.axvline(np.log10(PULLIN), color=st.INK, lw=0.8, ls=":", zorder=4)
    ax.set_xlim(-6.4, -0.4)
    ax.set_ylim(0, 66)
    ax.set_xticks(np.log10([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]))
    ax.set_xticklabels(["$10^{-6}$", "$10^{-5}$", "$10^{-4}$", "$10^{-3}$",
                        "$10^{-2}$", "$10^{-1}$"])
    ax.set_xlabel("sky registration tolerance (scaled)")
    ax.set_ylabel("pulsar–source pairs")
    st.panel_label(ax, "(a)")
    ax.text(np.log10(np.sqrt(LOOSEST * BLIND_LO)), 64,
            f"27–112$\\times$ gap\n$N_{{\\rm pairs}}=0$", ha="center", va="top",
            fontsize=6.2, weight="bold")
    ax.text(np.log10(0.102), 45, "blind\nfloat", ha="center", va="top",
            fontsize=6.2, weight="bold", color=st.VERMIL)
    ax.text(-5.5, 64, "all 348 rungs", ha="center", va="top", fontsize=6.2,
            weight="bold", color=st.BLUE)
    ax.annotate("loosest rung\n(J0711, nearest)", xy=(np.log10(LOOSEST), 12),
                xytext=(-2.75, 30), ha="center", fontsize=5.6, color=st.GREEN,
                arrowprops=dict(arrowstyle="->", color=st.GREEN, lw=0.6))
    ax.annotate("pull-in basin", xy=(np.log10(PULLIN), 5), xytext=(-6.2, 18),
                fontsize=5.6, color=st.INK,
                arrowprops=dict(arrowstyle="->", color=st.INK, lw=0.6))

    # ---- (b) the clocks
    y = np.arange(4)[::-1]
    vals = [Tl05, Tl15, Tc05, Tc15]
    labs = [f"loosest rung\n$\\sigma_{{\\rm sky}}\\propto T^{{-1/2}}$",
            "loosest rung\n$T^{-3/2}$ (optimistic)",
            "pull-in basin\n$T^{-1/2}$", "pull-in basin\n$T^{-3/2}$"]
    cols = [st.BLUE, st.SKY, st.VERMIL, st.ORANGE]
    axb.barh(y, vals, color=cols, edgecolor=st.INK, lw=0.4, height=0.62)
    axb.axvline(T0, color=st.GREEN, lw=1.0, ls="--")
    axb.text(T0 * 1.35, 2.5, "today\n(15 yr)", fontsize=5.6, color=st.GREEN,
             weight="bold", va="center")
    axb.set_xscale("log")
    axb.set_xlim(8, 3e7)
    axb.set_yticks(y)
    axb.set_yticklabels(labs, fontsize=5.4)
    axb.set_xlabel("baseline $T$ to cross (yr)")
    st.panel_label(axb, "(b)")
    for yy, v in zip(y, vals):
        axb.text(v * 1.5, yy, f"{v:,.0f}", va="center", fontsize=5.6,
                 weight="bold")

    # ---- (c) the tier referendum + break-even bound
    x = np.arange(3)
    axc.bar(x, REF_F, 0.55, yerr=REF_ERR, capsize=2.5,
            color=[st.GREY, st.ORANGE, st.BLUE], edgecolor=st.INK, lw=0.4,
            error_kw=dict(lw=0.8))
    axc.plot([2], [REF_TRAIL_C[0]], "o", ms=3, mfc="none", mec=st.INK, mew=0.7,
             zorder=6)
    axc.axhline(0.95, color=st.GREEN, lw=1.0, ls="--")
    axc.text(-0.42, 0.885, "required: 0.95", fontsize=6, color=st.GREEN,
             weight="bold")
    axc.set_xticks(x)
    axc.set_xticklabels(["sky\nonly", "+ EM\nperiod", "+ host\n$D_L$ (1%)"],
                        fontsize=6)
    axc.set_ylim(0, 1.05)
    axc.set_ylabel("posterior mass on truth, $f$")
    st.panel_label(axc, "(c)")

    axi = axc.inset_axes([0.35, 0.30, 0.62, 0.50])
    axi.plot(1.0 / BE_LAM, BE_F, marker="o", ms=2.5, lw=1.0, color=st.BLUE)
    axi.axhline(0.95, color=st.GREEN, lw=0.7, ls="--")
    axi.axvline(20, color=st.VERMIL, lw=0.9)
    axi.set_xscale("log")
    axi.set_xlim(0.8, 40)
    axi.set_ylim(0, 1.02)
    axi.set_xticks([1, 3, 8, 20])
    axi.set_xticklabels(["1", "3", "8", "20"], fontsize=5)
    axi.tick_params(labelsize=5, length=2)
    axi.set_xlabel("$\\sigma(\\log M_c)$ improvement", fontsize=5.2, labelpad=1)
    axi.set_ylabel("$f$", fontsize=5.4, labelpad=1)
    axi.text(19, 0.10, "$>20\\times$\nbound", fontsize=5.2, color=st.VERMIL,
             ha="right", weight="bold")
    for s in ("top", "right"):
        axi.spines[s].set_visible(True)

    axc.text(0.5, -0.30,
             "canonical STORY values (frozen 4-seed protocol; D-2);\n"
             "source npz cronus-only — not in this checkout",
             transform=axc.transAxes, ha="center", va="top", fontsize=5.0,
             style="italic", color=st.GREY)

    st.savefig(fig, "P2_the_walls")
    print(f"P2: loosest {LOOSEST:.3e}, gap pairs {n_gap}, Tcross loose "
          f"{Tl05:.0f}/{Tl15:.0f} yr, l2c {Tc05:.0f}/{Tc15:.0f} yr")


if __name__ == "__main__":
    main()
