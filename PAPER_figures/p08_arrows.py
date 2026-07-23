"""P8 — THE ARROWS.

(a) SPARK: the reservoir detection statistic AND the null floor per certified
    state — coherence raises the signal and LOWERS the floor (selectivity);
    the negative-signal-gain recruit (sC_m1e05) marked. N stated honestly:
    floors rest on 1300 nulls each; the signal side is ONE noise realisation;
    tabled N_cert are the MAX over 30 banked realisations (means 5.5/3.2/2.8).
(b) SPARK-3: the crossing ledger at the true marginal width (4/5 units
    survive, scrambled-clean) with both Fisher bounds, plus the
    model-quality-law panel (floor 118 unmodelled -> 276 conditional-width ->
    744 prior-width; truth arm 148 recovers).

Sources: SPARK_results/spark_s2c.npz; SPARK3_results/spark3_final_verdict.npz,
spark3_ledger.npz, s3r_A_g3_r2_k8.npz.
"""
import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"

STATES = [("s0", 0), ("sC_m2e03", 10), ("sC_m1e05", 10), ("sC_m1e07", 14),
          ("sMAX", 116)]


def main():
    S = np.load(f"{REPO}/SPARK_results/spark_s2c.npz", allow_pickle=True)

    fig = plt.figure(figsize=(st.DOUBLE, 2.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 1.0], wspace=0.45)
    ax = fig.add_subplot(gs[0])
    axl = fig.add_subplot(gs[1])
    axq = fig.add_subplot(gs[2])

    # ---------- (a) SPARK selectivity ----------
    rng = np.random.default_rng(3)
    xs = np.arange(len(STATES))
    for i, (stt, ncert) in enumerate(STATES):
        tf = np.asarray(S[f"twoF_{stt}"], float)
        fl = float(S[f"floor_{stt}_floor"])
        clr = int((tf > fl).sum())
        ax.plot(np.full(13, i) + rng.normal(0, 0.05, 13), tf, "o", ms=2.2,
                mfc=st.BLUE, mec="none", alpha=0.65, zorder=3)
        ax.plot([i - 0.28, i + 0.28], [fl, fl], color=st.VERMIL, lw=1.4,
                zorder=4)
        ax.plot([i - 0.28, i + 0.28],
                [np.median(tf)] * 2, color=st.INK, lw=0.9, ls=":", zorder=4)
        ax.text(i, 1.2, f"{clr}/13", ha="center", fontsize=6, weight="bold",
                color=st.GREEN)
    ax.plot([], [], "o", ms=2.4, mfc=st.BLUE, mec="none",
            label="13 reservoir sources (2F)")
    ax.plot([], [], color=st.VERMIL, lw=1.4, label="null floor (1300 nulls)")
    ax.plot([], [], color=st.INK, lw=0.9, ls=":", label="signal median")
    ax.annotate("recruits on the FLOOR DROP alone\n(signal gain $-0.20$)",
                xy=(2, 16.2), xytext=(0.24, 0.86), textcoords="axes fraction",
                fontsize=5.4, color=st.VERMIL, weight="bold",
                arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.6))
    ax.set_xticks(xs)
    ax.set_xticklabels([f"$N_{{\\rm cert}}$={n}\n{s.replace('_', ' ')}"
                        for s, n in STATES], fontsize=5.2)
    ax.set_ylabel("detection statistic $2F$ (nat)")
    ax.set_ylim(0, 80)
    st.panel_label(ax, "(a)")
    ax.legend(loc="upper left", fontsize=5.2)
    ax.text(0.5, -0.34,
            "oracle-anchored, no trials factor; signal side is ONE noise draw; "
            "$N_{\\rm cert}$ = max over 30 banked realisations (means 5.5/3.2/2.8)",
            transform=ax.transAxes, ha="center", va="top", fontsize=4.8,
            style="italic")

    # ---------- (b) SPARK-3 crossing ledger ----------
    V = np.load(f"{REPO}/SPARK3_results/spark3_final_verdict.npz",
                allow_pickle=True)
    units = [str(u) for u in V["unit"]]
    c_opt = np.asarray(V["cross_opt"], float)
    c_marg = np.asarray(V["cross_marg"], float)
    surv = np.asarray(V["survives"], bool)
    x = np.arange(len(units))
    w = 0.36
    axl.bar(x - w / 2, c_opt, w, color=st.SKY,
            label="optimistic Fisher bound")
    axl.bar(x + w / 2, c_marg, w, color=st.BLUE,
            label="TRUE marginal width")
    for xi, (sv, cm) in enumerate(zip(surv, c_marg)):
        axl.text(xi, max(c_opt[xi], cm) + 0.25,
                 "$\\checkmark$" if sv else "$\\times$", ha="center", fontsize=8,
                 color=(st.GREEN if sv else st.VERMIL), weight="bold")
    axl.set_xticks(x)
    axl.set_xticklabels(units, fontsize=6)
    axl.set_xlabel("crossing-eligible unit")
    axl.set_ylabel("fringe crossings")
    axl.set_ylim(0, 8.2)
    st.panel_label(axl, "(b)")
    axl.legend(loc="upper right", fontsize=4.8)
    axl.text(0.02, 0.98,
             f"at the true marginal width:\n{int(V['n_survive'])}/"
             f"{int(V['n_ready'])} units survive,\nscrambled-clean "
             f"({'yes' if bool(V['scrambled_clean']) else 'NO'}) — "
             f"{str(V['verdict'])}",
             transform=axl.transAxes, va="top", fontsize=5.4, weight="bold",
             color=st.GREEN)
    # pessimistic bound: zero crossings everywhere (ledger cross_pes = 0)
    L = np.load(f"{REPO}/SPARK3_results/spark3_ledger.npz", allow_pickle=True)
    axl.text(0.98, 0.55, f"pessimistic (prior-width)\nbound: "
             f"{int(L['cross_pes'])} crossings anywhere", fontsize=5.0,
             transform=axl.transAxes, va="top", ha="right", color=st.GREY)

    # ---------- (c) the model-quality law ----------
    U = np.load(f"{REPO}/SPARK3_results/s3r_A_g3_r2_k8.npz", allow_pickle=True)
    labs = ["unmodelled", "well-modelled\n(truth)", "conditional-\nwidth model",
            "prior-width\nmodel"]
    floors = [float(U["floor_un"]), float(U["floor_tr"]),
              float(U["floor_opt"]), float(U["floor_pes"])]
    cols = [st.GREY, st.GREEN, st.ORANGE, st.VERMIL]
    b = axq.bar(np.arange(4), floors, 0.6, color=cols, edgecolor=st.INK, lw=0.4)
    for r, v in zip(b, floors):
        axq.text(r.get_x() + r.get_width() / 2, v * 1.06, f"{v:.0f}",
                 ha="center", fontsize=6, weight="bold")
    axq.set_xticks(np.arange(4))
    axq.set_xticklabels(labs, fontsize=4.4)
    axq.set_ylabel("certification floor (nat)")
    axq.set_yscale("log")
    axq.set_ylim(60, 1600)
    st.panel_label(axq, "(c)")
    axq.annotate("a prior-width model\nINFLATES the null\n($118\\to744$): worse\nthan no model at all",
                 xy=(2.8, 700), xytext=(0.03, 0.72), textcoords="axes fraction",
                 fontsize=5.2, color=st.VERMIL, weight="bold",
                 arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.6))
    axq.text(0.5, -0.34, "unit A, realisation r2, coherence rung k8 "
             "(SPARK-3 §5.0); the law GLACIER now carries",
             transform=axq.transAxes, ha="center", va="top", fontsize=4.8,
             style="italic")

    st.savefig(fig, "P8_the_arrows")
    print(f"P8: floors {[round(float(S[f'floor_{s}_floor']), 3) for s, _ in STATES]}; "
          f"ledger {units} survive {surv.tolist()}; quality floors {floors}")


if __name__ == "__main__":
    main()
