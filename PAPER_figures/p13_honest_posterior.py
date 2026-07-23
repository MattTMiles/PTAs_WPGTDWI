"""P13 — THE HONEST POSTERIOR (sampler, two-regime by design).

STATUS: LAYOUT ONLY — ALL THREE PANELS PENDING. The brief's instruction is
"source everything, improvise NOTHING", and none of the three panels can be
sourced to a bank on this machine:

(a) above-onset fringe-marginalised distance posterior — the SAMPLER g3
    product; SAMPLER_dev.md §3: "the posterior gates did not run to
    completion on cronus" (g1/g2/g3/SBC all NOT RUN; ~7 GPU-hr/chain, killed
    by host-RAM OOM). No chain or fringe-posterior npz exists in this repo.
(b) the same machinery at census loudness — same gate; the diffuse-at-census
    statement is a stated prediction in SAMPLER_dev.md §4, not a banked
    posterior.
(c) sigma_MLE vs sigma_sampler ("4.5-17x") — no such comparison exists in the
    SAMPLER record; a sampler sigma would require the chains that were never
    produced. (The repo's 4-17x is a DIFFERENT quantity: ATLAS sky-area
    shrink per SNR doubling.)

This figure reserves the paper's layout and states the missing artifact per
panel, exactly as P11 does for GLACIER. See GAP_REPORT.md.
"""
import matplotlib.pyplot as plt

import easel5_style as st

MSG = [
    ("(a) above-onset realisation", "fringe-marginalised $p(L\\,|\\,d)$\n"
     "concentrated on the correct fringe",
     "SAMPLER g3: NOT RUN\n(no chain banked; cronus\nOOM, ~7 GPU-hr/chain)"),
    ("(b) census loudness", "same machinery — honestly\nmultimodal/diffuse\n"
     "(concentration tracks evidence)",
     "SAMPLER g3: NOT RUN\n(stated prediction only,\nSAMPLER_dev.md §4)"),
    ("(c) $\\sigma_{\\rm MLE}$ vs $\\sigma_{\\rm sampler}$",
     "linear-parameter covariance\nunderestimate, per parameter",
     "no such comparison in\nthe SAMPLER record;\nrequires the unrun chains"),
]


def main():
    fig, axes = plt.subplots(1, 3, figsize=(st.DOUBLE, 2.1))
    fig.subplots_adjust(wspace=0.30)
    for ax, (tag, what, why) in zip(axes, MSG):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(tag, fontsize=6.5, loc="left", pad=3)
        ax.text(0.5, 0.72, "PENDING", ha="center", va="center", fontsize=13,
                weight="bold", color=st.ORANGE, alpha=0.85,
                transform=ax.transAxes)
        ax.text(0.5, 0.46, what, ha="center", va="center", fontsize=5.6,
                transform=ax.transAxes)
        ax.text(0.5, 0.16, why, ha="center", va="center", fontsize=5.2,
                style="italic", color=st.GREY, transform=ax.transAxes)
        for s in ("top", "right"):
            ax.spines[s].set_visible(True)
    fig.suptitle("P13 — reserved layout: the sampler's posterior gates have "
                 "not run; no panel is sourceable to a bank on this machine",
                 fontsize=7, style="italic", y=1.00)
    st.savefig(fig, "P13_honest_posterior_PENDING")
    print("P13: layout reserved; all three panels PENDING (see GAP_REPORT.md)")


if __name__ == "__main__":
    main()
