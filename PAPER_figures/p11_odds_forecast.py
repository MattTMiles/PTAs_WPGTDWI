"""P11 — ODDS + FORECAST.

(a,b,c) LOTTERY's observer panel standalone: P(switch-on) under the
    channel-budget criterion, under the conservative threshold rule, and their
    one-sided disagreement, over the (f_ecc, e_char) mixture grid;
(d) a reserved half-panel for GLACIER's drain curve — layout only, marked
    PENDING: the Stage-1/2b fan is queued (held by dependency) and the drain
    bank A_bg(iter) does not yet exist; only the zero-noise Stage-0 precursor
    (reports/glacier_g2_population.npz) is banked, and it is not the
    measurement.

Sources: reports/lottery_tier1.npz (P_chan, P_thr, P_dis over f_grid x
e_grid; N_DRAW_MIX = 50,000 per cell); GLACIER status read, not touched.
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"


def main():
    L = np.load(f"{REPO}/reports/lottery_tier1.npz", allow_pickle=True)
    fg, eg = L["f_grid"], L["e_grid"]
    panels = [("P_chan", "channel-budget criterion\n($n_{\\rm active}\\geq 30$)"),
              ("P_thr", "conservative threshold rule\n($e\\geq0.5$ single / "
                        "$e\\geq0.3$ pair)"),
              ("P_dis", "disagreement (one-sided)\nchannel $-$ threshold")]

    fig = plt.figure(figsize=(st.DOUBLE, 2.35))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.1], wspace=0.30)
    ims = []
    for k, (key, title) in enumerate(panels):
        ax = fig.add_subplot(gs[k])
        im = ax.pcolormesh(fg, eg, L[key], cmap=st.SEQ_CMAP, vmin=0, vmax=1,
                           shading="nearest")
        ims.append(im)
        ax.set_title(title, fontsize=5.6, pad=4)
        ax.set_xlabel("$f_{\\rm ecc}$", labelpad=1)
        if k == 0:
            ax.set_ylabel("$e_{\\rm char}$")
        else:
            ax.set_yticklabels([])
        ax.text(0.04, 0.975, f"({'abc'[k]})", transform=ax.transAxes,
                color="white", weight="bold", fontsize=8, va="top")
        ax.tick_params(labelsize=5)
        heat_axes = heat_axes + [ax] if k else [ax]
    cb = fig.colorbar(ims[0], ax=heat_axes, fraction=0.025, pad=0.015)
    cb.set_label("$P$(switch-on)", fontsize=6)
    cb.ax.tick_params(labelsize=5)

    # (d) the reserved GLACIER half-panel
    axg = fig.add_subplot(gs[3])
    axg.set_xlim(0, 16)
    axg.set_ylim(0, 1)
    axg.set_xlabel("resolution iteration", labelpad=1)
    axg.set_ylabel("$A_{\\rm bg}$ (background amplitude)")
    axg.set_yticks([])
    axg.tick_params(labelsize=5)
    st.panel_label(axg, "(d)", x=0.03, y=0.955)
    axg.text(0.5, 0.60, "PENDING", ha="center", va="center", fontsize=14,
             weight="bold", color=st.ORANGE, alpha=0.85,
             transform=axg.transAxes)
    axg.text(0.5, 0.32,
             "GLACIER drain curve $A_{\\rm bg}$(iter):\nfan queued; no drain "
             "bank exists yet.\nReserved layout — only the zero-noise\nStage-0 "
             "precursor is banked\n(glacier_g2_population.npz).",
             ha="center", va="center", fontsize=5.2, style="italic",
             transform=axg.transAxes)
    axg.set_title("GLACIER: eating the\nbackground (reserved)", fontsize=5.8,
                  pad=3)

    st.savefig(fig, "P11_odds_forecast")
    print(f"P11: grid {L['P_chan'].shape}, P_chan range "
          f"{L['P_chan'].min():.2f}-{L['P_chan'].max():.2f}; GLACIER PENDING")


if __name__ == "__main__":
    main()
