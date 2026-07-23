"""P5 — THE ONSET MAP.

The full corrected surface: correct-certification-count heatmaps over (h, T)
for the COMPLETE structure x prior-tier grid — all six (structure, tier)
panels, 108 cells. The onset contour is drawn on the count at floor+error
(corr_lo > 1, the CONFIRMED test); its calibration-error band is the region
between the corr>1 and corr_lo>1 contours. MARGINAL cells (clear at the floor,
die at floor+error) are hatched per the floor-provenance convention; cells
whose floor was re-fit on the empirical q95 (Gumbel invalid, criterion-v2.2)
carry a dotted border.

Source: reports/surface_analysis_recut.npz (108 cells x 30 signal
realisations, per-cell floors N=100, alpha=0.05, criterion-v2.2).
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"

H = [-13.25, -13.0, -12.75, -12.5, -12.25, -12.0]
TS = [30, 40, 50]
STRUCTS = ["2+14", "3+13", "5+11"]
TIERS = ["lit", "vlbi"]


def main():
    A = np.load(f"{REPO}/reports/surface_analysis_recut.npz", allow_pickle=True)
    ci = {str(c): i for i, c in enumerate(A["cols"])}
    T = A["table"]
    tiers = np.array([str(t) for t in A["tiers"]])
    struct = np.array([str(s) for s in A["struct"]])
    verd = np.array([str(v) for v in A["verdicts"]])
    est = np.array([str(e) for e in A["estimator"]])

    def grid(k, tier, key=None, flag=None):
        mm = (struct == k) & (tiers == tier)
        G = np.full((len(TS), len(H)), np.nan)
        for j in np.flatnonzero(mm):
            r = T[j]
            a = TS.index(int(r[ci["T"]]))
            b = H.index(round(r[ci["h"]], 2))
            if key is not None:
                G[a, b] = r[ci[key]]
            else:
                G[a, b] = flag[j]
        return G

    fig, axes = plt.subplots(2, 3, figsize=(st.DOUBLE, 3.6), sharex=True,
                             sharey=True)
    fig.subplots_adjust(wspace=0.06, hspace=0.16, right=0.90)
    vmax = 8.0
    n_onset_check = 0
    for irow, tier in enumerate(TIERS):
        for icol, k in enumerate(STRUCTS):
            ax = axes[irow, icol]
            Cg = grid(k, tier, "corr")
            CL = grid(k, tier, "corr_lo")
            MG = grid(k, tier, flag=(verd == "MARGINAL").astype(float))
            RF = grid(k, tier, flag=(est == "emp_q95").astype(float))
            n_onset_check += int((grid(k, tier,
                                       flag=(verd == "ONSET").astype(float)) > 0.5).sum())
            im = ax.imshow(Cg, origin="lower", aspect="auto", cmap=st.SEQ_CMAP,
                           vmin=0, vmax=vmax,
                           extent=[-0.5, len(H) - 0.5, -0.5, len(TS) - 0.5])
            xg, yg = np.arange(len(H)), np.arange(len(TS))
            # calibration-error band: thin dashed contour = bar at the nominal
            # floor; bold contour = bar at floor+error (the onset test). The
            # region between the two is the floor-calibration error band.
            ax.contour(xg, yg, Cg, levels=[1.0], colors="white",
                       linewidths=0.8, linestyles="--")
            ax.contour(xg, yg, CL, levels=[1.0], colors="white", linewidths=2.0)
            ax.contour(xg, yg, CL, levels=[1.0], colors=[st.VERMIL],
                       linewidths=0.9)
            for a in range(len(TS)):
                for b in range(len(H)):
                    if MG[a, b] > 0.5:
                        st.marginal_hatch(ax, b - 0.5, a - 0.5)
                    if RF[a, b] > 0.5:
                        ax.add_patch(Rectangle((b - 0.5, a - 0.5), 1, 1,
                                               fill=False, edgecolor="white",
                                               lw=0.8, ls=(0, (1.5, 1.5)),
                                               zorder=4))
                    v = Cg[a, b]
                    ax.text(b, a, f"{v:.1f}", ha="center", va="center",
                            fontsize=5.4, weight="bold",
                            color="white" if v < 0.55 * vmax else "black",
                            zorder=6)
            ax.set_xticks(range(len(H)))
            ax.set_xticklabels(["$-13.25$", "$-13$", "$-12.75$", "$-12.5$",
                                "$-12.25$", "$-12$"], fontsize=5.4, rotation=45)
            ax.set_yticks(range(len(TS)))
            ax.set_yticklabels([f"{t}" for t in TS], fontsize=6)
            if irow == 1:
                ax.set_xlabel("$\\log_{10} h$")
            if icol == 0:
                ax.set_ylabel(f"{tier} tier\n$T$ (yr)")
            if irow == 0:
                nl = {"2+14": "2 loud + 14 faint", "3+13": "3 loud + 13 faint (census)",
                      "5+11": "5 loud + 11 faint"}[k]
                ax.set_title(nl, fontsize=7, pad=3)

    assert n_onset_check == int(A["n_onset"]) == 59, n_onset_check
    cax = fig.add_axes([0.92, 0.13, 0.015, 0.74])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("correct certifications / realisation", fontsize=7)
    cb.ax.tick_params(labelsize=5.5)

    fig.subplots_adjust(bottom=0.17)
    fig.text(0.5, 0.012,
             "solid contour: onset boundary on the count at floor$+$error; thin dashed contour: "
             "bar at the nominal floor — the region between is the calibration-error band. "
             "Hatched: MARGINAL cells (3). Dotted border: floor re-fit on the empirical q95 "
             "(null zero-fraction $>20\\%$; 15 cells). $N_{\\rm onset}=59$ of 108.",
             ha="center", fontsize=5.6, style="italic")
    st.savefig(fig, "P5_the_onset_map")
    print(f"P5: onset cells {n_onset_check} (bank {int(A['n_onset'])}), "
          f"marginal {int(A['n_marginal'])}, refit {int((est=='emp_q95').sum())}")


if __name__ == "__main__":
    main()
