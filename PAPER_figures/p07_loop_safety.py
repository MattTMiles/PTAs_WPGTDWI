"""P7 — THE LOOP: SAFETY.

(a) ALL banked soft-loop trajectories, rendered as the three-claim
    decomposition (S8.5.0, no borrowed adjectives):
      1. seeded at/near truth (IGNITE-2 + CHORUS + EMBER truth arm): HOLDS,
         occasionally +1, and SELF-CLEANS scrambled false alarms;
      2. honest cold starts (EMBER map arm): INERT — DeltaN = 0.000 in all 9
         cells;
      3. danger mode (scrambled counterpart): manufacturing, whose boundary is
         MOTION — with the manufacture-anatomy inset (accepted M-step motion
         vs DeltaN across the 112-realisation scrambled census).
(b) the EMBER predictor comparison: motion sensitivity 1.00 vs engagement
    0.60; the joint (engaged AND mobile) form refuted.

Sources: IGNITE2_results/ig_sloop*.npz (40) + CHORUS_results/ch_sloop*.npz
(30) trajectory banks; EMBER_results/ember_analysis.npz (396 realisation
rows, per-cell verdicts); EMBER_results/ember_predictors.npz (112-realisation
scrambled census, contingency tables, Fisher p, logistic betas).
"""
import glob

import numpy as np
import matplotlib.pyplot as plt

import easel5_style as st

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"


def main():
    # ---------------- trajectories ----------------
    fs = sorted(glob.glob(f"{REPO}/IGNITE2_results/ig_sloop*_*.npz")) + \
         sorted(glob.glob(f"{REPO}/CHORUS_results/ch_sloop*_*.npz"))
    T = []
    for f in fs:
        d = np.load(f, allow_pickle=True)
        if "traj_n_cert_true" not in d.files:
            continue
        T.append(dict(cert=np.asarray(d["traj_n_cert_true"], float),
                      wrong=np.asarray(d["traj_wrong"], float),
                      scram=bool(d["scrambled"])))
    grew = [t for t in T if t["cert"][-1] > t["cert"][0]]
    shed = [t for t in T if t["wrong"][-1] < t["wrong"][0]]
    special = {id(t) for t in grew} | {id(t) for t in shed}
    flat = [t for t in T if id(t) not in special]

    E = np.load(f"{REPO}/EMBER_results/ember_analysis.npz", allow_pickle=True)
    r_arm = np.array([str(a) for a in E["r_arm"]])
    r_scr = np.asarray(E["r_scr"], bool)
    r_dN = np.asarray(E["r_dN"], float)
    map_dN = r_dN[(r_arm == "map") & ~r_scr]
    truth_dN = r_dN[(r_arm == "truth") & ~r_scr]
    verdicts = [str(v) for v in E["verdict"]]
    assert all(v == "HOLDS" for v in verdicts) and len(verdicts) == 9
    assert float(np.abs(map_dN).max()) == 0.0        # the INERT claim, re-read

    P = np.load(f"{REPO}/EMBER_results/ember_predictors.npz", allow_pickle=True)
    scr_dN = np.asarray(P["scr_dN"], float)
    scr_mstep = np.asarray(P["scr_mstep"], float)
    scr_man = np.asarray(P["scr_manufacture"], bool)
    mo_thr = float(P["motion_thr"])
    fp = np.asarray(P["fisher_p"], float)
    fpn = [str(x) for x in P["fisher_p_names"]]
    beta = np.asarray(P["logit_beta"], float)

    fig = plt.figure(figsize=(st.DOUBLE, 2.65))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.42)
    ax = fig.add_subplot(gs[0])
    axm = fig.add_subplot(gs[1])
    axp = fig.add_subplot(gs[2])

    # (a) trajectories: claim 1 (seeded) + shed (self-clean)
    rng = np.random.default_rng(7)
    for t in flat:
        nn = len(t["cert"])
        ax.plot(np.arange(nn), t["cert"] + rng.normal(0, 0.045, nn), lw=0.5,
                color=st.GREY, alpha=0.4, zorder=2)
    for t in shed:
        ax.plot(np.arange(len(t["wrong"])), t["wrong"], lw=1.1, color=st.VERMIL,
                ls="--", marker="o", ms=2.2, alpha=0.9, zorder=4)
    for t in grew:
        ax.plot(np.arange(len(t["cert"])), t["cert"], lw=1.3, color=st.GREEN,
                marker="o", ms=2.4, zorder=5)
    ax.plot([], [], color=st.GREY, lw=1.0,
            label=f"held flat ({len(flat)}/{len(T)})")
    ax.plot([], [], color=st.GREEN, lw=1.2, marker="o", ms=2.4,
            label=f"gained a true cert (+1) ({len(grew)})")
    ax.plot([], [], color=st.VERMIL, lw=1.2, ls="--", marker="o", ms=2.2,
            label=f"false alarm SHED by loop ({len(shed)})")
    ax.set_xlabel("loop iteration")
    ax.set_ylabel("certified pulsars\n(solid true / dashed false)")
    ax.set_xlim(-0.2, 9.4)
    ax.set_ylim(-0.35, 4.4)
    st.panel_label(ax, "(a)", y=1.10)
    ax.legend(loc="upper right", fontsize=5.2)
    ax.set_title("claim 1 — seeded: HOLDS, self-cleans (70 loops)",
                 fontsize=6.2, loc="left", pad=8)
    ax.text(0.02, 0.03,
            f"claim 2 — honest cold starts (EMBER map arm, "
            f"$N={len(map_dN)}$):\nINERT, $\\Delta N=0.000$ in all 9 cells "
            f"(max$|\\Delta N|$={np.abs(map_dN).max():.0f});\n"
            f"EMBER truth arm: repair $\\Delta N$ up to "
            f"+{truth_dN.max():.0f} (claim 1)",
            transform=ax.transAxes, fontsize=5.2, va="bottom",
            bbox=st.NOTE_BOX)

    # (b) manufacture anatomy: motion vs dN over the scrambled census
    j = rng.normal(0, 0.03, len(scr_dN))
    axm.plot(np.maximum(scr_mstep[~scr_man], 2e-3), scr_dN[~scr_man] + j[~scr_man],
             "o", ms=2.2, mfc=st.GREY, mec="none", alpha=0.6,
             label=f"scrambled, safe ({int((~scr_man).sum())})")
    axm.plot(np.maximum(scr_mstep[scr_man], 2e-3), scr_dN[scr_man],
             "D", ms=4.0, mfc=st.VERMIL, mec=st.INK, mew=0.5,
             label=f"MANUFACTURED ({int(scr_man.sum())})")
    axm.axvline(mo_thr, color=st.INK, lw=0.7, ls=":")
    axm.text(mo_thr * 0.85, -3.4, "motion\nthreshold", fontsize=5.2, ha="right")
    axm.set_xscale("log")
    axm.set_xlabel("accepted M-step motion (scaled; floor at $2\\times10^{-3}$)")
    axm.set_ylabel("$\\Delta N$ (false certs grown)")
    st.panel_label(axm, "(b)", y=1.10)
    axm.legend(loc="lower left", fontsize=5.2)
    axm.set_title(f"claim 3 — motion under a wrong cpart ($N={len(scr_dN)}$)",
                  fontsize=6.2, loc="left", pad=8)
    axm.text(0.97, 0.05, "every manufacturer\nMOVED (5/5)",
             transform=axm.transAxes, ha="right", fontsize=5.6,
             color=st.VERMIL, weight="bold")

    # (c) predictor comparison
    names = ["engagement\n($C_0\\geq1$)", "MOTION\n(step$>0$)",
             "engaged\nAND mobile", "engaged\nOR mobile"]
    # sens/spec derived from the banked contingency tables (TP,FP,FN,TN);
    # the report rounds motion specificity to 0.72, the bank gives 80/107=0.75
    sens, spec = [], []
    for key in ["cont_engagement", "cont_motion", "cont_and", "cont_or"]:
        tp, fp_, fn, tn = np.asarray(P[key], float)
        sens.append(tp / (tp + fn))
        spec.append(tn / (tn + fp_))
    assert np.allclose(sens, [0.60, 1.00, 0.60, 1.00], atol=0.005)
    x = np.arange(4)
    w = 0.36
    axp.bar(x - w / 2, sens, w, color=st.BLUE, label="sensitivity")
    axp.bar(x + w / 2, spec, w, color=st.SKY, label="specificity")
    for xi, (nm, pv) in enumerate(zip(fpn, fp)):
        axp.text(xi, 1.045, f"p={pv:.3g}", ha="center", fontsize=5.0,
                 color=(st.GREEN if pv < 0.01 else st.INK))
    axp.set_xticks(x)
    axp.set_xticklabels(names, fontsize=4.5)
    axp.set_ylim(0, 1.18)
    axp.set_ylabel("fraction")
    st.panel_label(axp, "(c)")
    axp.legend(loc="upper left", fontsize=5.0, borderaxespad=0.2)
    axp.annotate("pre-registered joint\nboundary REFUTED\n(misses 2/5)",
                 xy=(2, 0.62), xytext=(1.25, 0.30), fontsize=5.2,
                 color=st.VERMIL, weight="bold",
                 arrowprops=dict(arrowstyle="->", color=st.VERMIL, lw=0.6))
    axp.text(0.5, -0.30,
             f"logistic: $\\beta_{{\\rm motion}}=+{beta[2]:.3f}$ > "
             f"$\\beta_{{\\rm eng}}=+{beta[1]:.3f}$; positive control: 7/7 "
             "legitimate ignitions moved",
             transform=axp.transAxes, fontsize=4.8, va="top", ha="center",
             style="italic")

    st.savefig(fig, "P7_the_loop_safety")
    print(f"P7: {len(T)} seeded loops (flat {len(flat)}, grew {len(grew)}, "
          f"shed {len(shed)}); map-arm max|dN|=0; manufacturers "
          f"{int(scr_man.sum())}/112; fisher {dict(zip(fpn, np.round(fp, 4)))}")


if __name__ == "__main__":
    main()
