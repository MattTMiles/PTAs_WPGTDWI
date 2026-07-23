"""P14 — THE ARCHAEOLOGY: why previous samplers' confidence was not evidence.

(a) q_max vs loudness with snap-purity overlaid: median q_max rises
    0.16 -> 0.94 across the loudness ladder while the purity of confident
    snaps falls 0.79 -> 0.37 — confidence rises exactly where correctness
    collapses;
(b) P(confident AND true) vs loudness — the small, flat joint quantity — with
    the MK9 corner cited: the regime where delta-on-truth held (5 pulsars,
    1 source, h = 10^-12, 1 us white noise, coverage 13/13) is a corner the
    realistic bank never revisits. The MK9 corner is notebook-only
    (CW_node_sampling/data_likelihood_sandbox_MK9_executed.ipynb, quoted via
    SAMPLER_dev.md §8.2); it is cited as an annotation, not plotted as data.

Sources: SAMPLER_results/archaeology_snap.npz (IGNITE signal bank, lit tier,
T=30); raw re-derivation bank reports/ignite_bank.npz.
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"


def main():
    A = np.load(f"{REPO}/SAMPLER_results/archaeology_snap.npz",
                allow_pickle=True)
    h = A["h"]
    q = A["med_qmax"]
    pur = A["purity_q90"]
    fq90 = A["frac_q90"]
    pct = A["p_confident_and_true"]

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(st.SINGLE * 2.05, 2.3))
    fig.subplots_adjust(wspace=0.35)

    # (a) confidence vs purity
    ax.plot(h, q, "-o", ms=3, lw=1.1, color=st.BLUE,
            label="median $q_{\\rm max}$ (confidence)")
    ax.plot(h, pur, "-s", ms=3, lw=1.1, color=st.VERMIL,
            label="purity of $q>0.9$ snaps (correctness)")
    ax.plot(h, fq90, "-^", ms=2.6, lw=0.8, color=st.GREY, alpha=0.8,
            label="fraction with $q>0.9$")
    for x0, y0, t, va in ((h[0], q[0] + 0.045, f"{q[0]:.2f}", "bottom"),
                          (h[-1], q[-1] - 0.06, f"{q[-1]:.2f}", "top")):
        ax.text(x0, y0, t, fontsize=5.4, color=st.BLUE, ha="center", va=va)
    for x0, y0, t in ((h[0], pur[0], f"{pur[0]:.2f}"),
                      (h[-1], pur[-1], f"{pur[-1]:.2f}")):
        ax.text(x0, y0 + 0.045, t, fontsize=5.4, color=st.VERMIL, ha="center")
    ax.set_xlabel("$\\log_{10} h$")
    ax.set_ylabel("fraction / probability")
    ax.set_ylim(0, 1.05)
    st.panel_label(ax, "(a)")
    ax.legend(loc="center left", fontsize=5.2)
    ax.text(0.03, 0.04, "confidence rises exactly where\ncorrectness collapses",
            transform=ax.transAxes, fontsize=5.6, weight="bold")

    # (b) the joint quantity + the MK9 corner citation
    axb.plot(h, pct, "-o", ms=3, lw=1.1, color=st.GREEN,
             label="P(confident AND true)")
    axb.set_xlabel("$\\log_{10} h$")
    axb.set_ylabel("P($q>0.9$ AND true fringe)")
    axb.set_ylim(0, 0.55)
    st.panel_label(axb, "(b)")
    axb.legend(loc="upper left", fontsize=5.4)
    axb.text(-13.24, 0.34,
             "the MK9 corner (where $\\delta$-on-truth held):\n"
             "5 pulsars, 1 source, $h=10^{-12}$, $1\\,\\mu$s white,\n"
             "coverage 13/13 — a regime $\\sim$56$\\times$ louder\nthan realistic; "
             "notebook-only artifact$^{*}$,\ncited not plotted",
             fontsize=5.2, va="top", bbox=st.NOTE_BOX)
    axb.text(0.97, 0.03, "$^{*}$data_likelihood_sandbox_MK9_executed.ipynb\n"
             "via SAMPLER_dev.md §8.2 — no npz bank",
             transform=axb.transAxes, ha="right", fontsize=4.6,
             style="italic", color=st.GREY)

    st.savefig(fig, "P14_the_archaeology")
    print(f"P14: q {q[0]:.2f}->{q[-1]:.2f}, purity {pur[0]:.2f}->{pur[-1]:.2f}, "
          f"joint {pct.min():.2f}-{pct.max():.2f}")


if __name__ == "__main__":
    main()
