"""P10 — GEOMETRY + WARNING.

(a) per-pulsar certification frequency across the 40-sky GEO ensemble (the
    geometry lottery), census names coloured; 98 of 116 never certify;
(b) RING — sky-localisation bias vs SNR with coverage collapse: the case for
    pc-class distance campaigns as bias PREVENTION.

Sources: reports/geo_summary.npz (40 isotropic sky draws, census population,
zero noise, honest gate); RING_results/npz/cell_fgw-8_A_snr{5,10,20}.npz
(scenario A, 10 realisations/cell; 'feather' = today's real distance priors,
'tier3' = near-exact pc-class distances).
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]
SNR = [5, 10, 20]


def main():
    g = np.load(f"{REPO}/reports/geo_summary.npz", allow_pickle=True)
    nm = [str(x) for x in g["names"]]
    fr = g["freq"]
    keep = np.where(fr > 0)[0]
    order = keep[np.argsort(-fr[keep])]

    D = {}
    for t in ("feather", "tier3"):
        bias, ins = [], []
        for s in SNR:
            d = np.load(f"{REPO}/RING_results/npz/cell_fgw-8_A_snr{s}.npz",
                        allow_pickle=True)
            bias.append(np.median(d[f"{t}__map_offset_local_deg"]))
            ins.append(np.mean(d[f"{t}__inside90_local"]) * 100)
        D[t] = (np.array(bias), np.array(ins))

    fig = plt.figure(figsize=(st.DOUBLE, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 0.9, 0.9], wspace=0.45)
    ax = fig.add_subplot(gs[0])
    axb = fig.add_subplot(gs[1])
    axc = fig.add_subplot(gs[2])

    # (a) the lottery
    y = np.arange(len(order))[::-1]
    cols = [st.VERMIL if nm[i] in CENSUS else st.BLUE for i in order]
    ax.barh(y, fr[order] * 100, color=cols, edgecolor="white", lw=0.3,
            height=0.75)
    ax.set_yticks(y)
    ax.set_yticklabels([nm[i] for i in order], fontsize=4.8)
    ax.set_xlabel("skies in which the pulsar certifies (%)")
    ax.set_xlim(0, 108)
    st.panel_label(ax, "(a)")
    for yy, i in zip(y, order):
        ax.text(fr[i] * 100 + 1.5, yy, f"{fr[i]*100:.0f}", va="center",
                fontsize=4.4)
    ax.text(0.97, 0.06,
            f"only {len(order)} of 116 ever certify;\n"
            f"the other {116-len(order)} never do, in any sky.\n"
            "red: the census names — reproduced\ntogether in 0 of 40 skies",
            transform=ax.transAxes, ha="right", fontsize=5.4, bbox=st.NOTE_BOX)

    # (b) bias does not shrink
    for t, lab, col in (("feather", "today's distances", st.VERMIL),
                        ("tier3", "pc-class distances", st.GREEN)):
        axb.plot(SNR, D[t][0], "-o", ms=3, lw=1.0, color=col, label=lab)
    axb.set_xticks(SNR)
    axb.set_xlim(3.5, 21.5)
    axb.set_ylim(0, 7.2)
    axb.set_xlabel("source SNR")
    axb.set_ylabel("sky MAP offset from truth (deg)")
    st.panel_label(axb, "(b)")
    axb.legend(loc="center right", fontsize=5.4)
    axb.text(0.5, 0.86, "bias is SNR-independent —\nwrong distances MOVE the\n"
             "answer, not blur it", transform=axb.transAxes, ha="center",
             fontsize=5.4, color=st.VERMIL, weight="bold")

    # (c) coverage collapse
    for t, lab, col in (("feather", "today's distances", st.VERMIL),
                        ("tier3", "pc-class distances", st.GREEN)):
        axc.plot(SNR, D[t][1], "-o", ms=3, lw=1.0, color=col, label=lab)
    for s, v in zip(SNR, D["feather"][1]):
        axc.text(s, v + 5, f"{v:.0f}%", ha="center", fontsize=5.6,
                 weight="bold", color=st.VERMIL)
    axc.axhline(90, color=st.INK, lw=0.7, ls="--")
    axc.text(20.8, 92.5, "nominal 90%", ha="right", fontsize=5.2)
    axc.set_xticks(SNR)
    axc.set_xlim(3.5, 21.5)
    axc.set_ylim(-6, 112)
    axc.set_xlabel("source SNR")
    axc.set_ylabel("truth inside the 90% box (%)")
    st.panel_label(axc, "(c)")
    axc.text(0.05, 0.10, "louder $\\Rightarrow$ MORE confidently\nwrong: distance "
             "campaigns are\nbias PREVENTION", transform=axc.transAxes,
             fontsize=5.4, color=st.INK, weight="bold")

    st.savefig(fig, "P10_geometry_warning")
    print(f"P10: union {len(order)}; feather inside90 {D['feather'][1]}, "
          f"tier3 {D['tier3'][1]}")


if __name__ == "__main__":
    main()
