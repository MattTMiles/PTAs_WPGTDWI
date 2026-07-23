"""P15 — CALIBRATION + THE SADDLE.

(a) S-4 coverage, exactly as banked: per-realisation credible level of truth
    under the 8-parameter loud-source sampler — scenario B (control,
    white+RN; N=2, both inside their 90% box) vs scenario C (full GWB
    marginalisation; N=1, truth at the 94.3rd percentile, outside). N is
    stated on the panel; the brief's "N=5" does not match the bank, and the
    SAMPLER_dev.md prose table quotes C values (0.903, 0.813) that are not in
    the npz — per the ARTIFACT READBACK convention the npz is plotted.
(b) truth-is-a-saddle: the banked Hessian condition numbers (3.2e8 scenario
    C); the eigenvalue spectrum ("3/8 negative") was read off a log
    (sampler_s4b.log) that is NOT in this checkout and is quoted as
    PENDING-provenance text, not drawn as data.

Sources: SAMPLER_results/s4_scenB.npz, s4_scenC.npz (keys c, inside90, cond,
seeds). See GAP_REPORT.md.
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"


def main():
    B = np.load(f"{REPO}/SAMPLER_results/s4_scenB.npz", allow_pickle=True)
    C = np.load(f"{REPO}/SAMPLER_results/s4_scenC.npz", allow_pickle=True)
    cB = np.asarray(B["c"], float).ravel()
    cC = np.asarray(C["c"], float).ravel()
    inB = np.asarray(B["inside90"], bool).ravel()
    inC = np.asarray(C["inside90"], bool).ravel()
    condB = np.asarray(B["cond"], float).ravel()
    condC = np.asarray(C["cond"], float).ravel()

    fig, (ax, axs) = plt.subplots(1, 2, figsize=(st.SINGLE * 2.05, 2.3),
                                  gridspec_kw={"width_ratios": [1.15, 1]})
    fig.subplots_adjust(wspace=0.35)

    # (a) coverage
    for x0, cc, ii, col, lab in ((0, cB, inB, st.BLUE,
                                  f"B: control, white+RN ($N$={len(cB)})"),
                                 (1, cC, inC, st.ORANGE,
                                  f"C: + full GWB marg. ($N$={len(cC)})")):
        for k, (v, ins) in enumerate(zip(cc, ii)):
            ax.plot([x0 + (k - (len(cc) - 1) / 2) * 0.14], [v],
                    "o" if ins else "X", ms=6, mfc=col, mec=st.INK, mew=0.5)
        ax.plot([], [], "o", ms=5, mfc=col, mec=st.INK, mew=0.5, ls="none",
                label=lab)
    ax.axhline(0.90, color=st.VERMIL, lw=0.9, ls="--")
    ax.text(-0.44, 0.915, "90% credible level", fontsize=5.4, color=st.VERMIL)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["scenario B\n(control)", "scenario C\n(full GWB marg.)"],
                       fontsize=6)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("credible level of truth")
    st.panel_label(ax, "(a)")
    ax.legend(loc="center left", fontsize=5.0)
    ax.text(0.5, -0.30,
            "circles: truth inside its 90% box; X: outside. Banked N = 2 (B) "
            "+ 1 (C) realisations —\nthe run was killed before the planned "
            "ensemble; N stated, nothing extrapolated.",
            transform=ax.transAxes, ha="center", va="top", fontsize=4.8,
            style="italic")

    # (b) the saddle mechanism
    xs = [0, 0.35, 1.0]
    vals = list(condB) + list(condC)
    cols = [st.BLUE, st.BLUE, st.ORANGE]
    labs = ["B seed\n770501", "B seed\n771501", "C seed\n770001"]
    b = axs.bar([0, 0.6, 1.5], vals, 0.45, color=cols, edgecolor=st.INK, lw=0.4)
    for r, v in zip(b, vals):
        axs.text(r.get_x() + r.get_width() / 2, v * 1.3, f"{v:.1e}",
                 ha="center", fontsize=5.2, weight="bold")
    axs.set_yscale("log")
    axs.set_ylim(1e7, 6e10)
    axs.set_xticks([0, 0.6, 1.5])
    axs.set_xticklabels(labs, fontsize=5.4)
    axs.set_ylabel("Hessian condition number at truth")
    st.panel_label(axs, "(b)")
    axs.text(0.5, 0.94,
             "truth is a SADDLE: quoted spectrum '3/8 negative\neigenvalues' "
             "was read off sampler_s4b.log —\nlog absent from this checkout; "
             "only cond is banked.\nMode-finding fails at truth even standing "
             "on it.",
             transform=axs.transAxes, ha="center", va="top", fontsize=5.2,
             bbox=st.NOTE_BOX)

    st.savefig(fig, "P15_calibration_saddle")
    print(f"P15: cB={cB}, cC={cC}, condC={condC}")


if __name__ == "__main__":
    main()
