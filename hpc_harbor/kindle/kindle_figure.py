"""KINDLE figure — the F15/F11-KINDLE panel. House style = ignite2_figures.py.
Panel A: the ratio-gain degeneracy (transition census: NaN at ignitions, 1.000 tautologies).
Panel B: additive dN vs ratio gain, across the 70-trajectory banked corpus.
Panel C: the off-truth pilot — arm(a) truth vs arm(b-FIX) 1sig(mc), dN and margin flow
         [drawn only if pilot banks exist; otherwise annotated 'pending GPU'].
CPU-only, matplotlib Agg."""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/data/taylor_group/matt_miles/PTAs_WPGTDWI"
KRES = f"{REPO}/KINDLE_results"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C_TRUE, C_ALL, C_WRONG, C_FIX = "#1baf7a", "#2a78d6", "#d0362f", "#8047c9"

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "font.size": 9,
    "axes.titlesize": 9.5, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.8, "legend.frameon": False,
})


def load_corpus():
    rows = []
    for pat, camp in [(f"{REPO}/reports/ig_sloop*.npz", "IG2"),
                      (f"{REPO}/CHORUS_results/ch_sloop*.npz", "CHO")]:
        for f in sorted(glob.glob(pat)):
            d = np.load(f, allow_pickle=True)
            if bool(d["scrambled"]):
                continue
            nc = d["traj_n_cert"].astype(int)
            rows.append((camp, nc, np.asarray(d["traj_gain"], float)))
    return rows


def main():
    corpus = load_corpus()
    pilot = sorted(glob.glob(f"{KRES}/k_pilot*.npz"))
    have_pilot = len(pilot) > 0
    ncol = 3 if have_pilot else 2
    fig, axs = plt.subplots(1, ncol, figsize=(4.6 * ncol, 4.3))

    # ---- Panel A: transition census ----
    ax = axs[0]
    trans = {}
    for _, nc, g in corpus:
        for i in range(1, len(nc)):
            a, b = int(nc[i - 1]), int(nc[i])
            lab = "NaN" if not np.isfinite(g[i]) else f"{g[i]:.1f}"
            key = (a, b)
            trans.setdefault(key, {"NaN": 0, "fin": 0})
            trans[key]["NaN" if lab == "NaN" else "fin"] += 1
    labels, nan_c, fin_c, colors = [], [], [], []
    for (a, b), c in sorted(trans.items()):
        labels.append(f"{a}→{b}")
        nan_c.append(c["NaN"]); fin_c.append(c["fin"])
        colors.append(C_TRUE if (a == 0 and b > 0) else INK2)
    y = np.arange(len(labels))
    ax.barh(y, nan_c, color="#c9c7c1", label="ratio gain = NaN (|C|=0 guard)")
    ax.barh(y, fin_c, left=nan_c, color=C_ALL, label="ratio gain finite")
    for i, (a, b) in enumerate(sorted(trans.keys())):
        if a == 0 and b > 0:
            ax.text(nan_c[i] + fin_c[i] + 0.4, i, "IGNITION\n(NaN'd)", va="center",
                    fontsize=7, color=C_TRUE, fontweight="bold")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("count of (iteration transitions), 50 honest runs")
    ax.set_ylabel("|C| transition (prev → now)")
    ax.set_title("(a) every ignition is NaN'd; every finite\ngain is an n→n tautology",
                 loc="left")
    ax.legend(fontsize=7, loc="lower right")

    # ---- Panel B: additive dN vs ratio gain support ----
    ax = axs[1]
    dN = np.array([int(nc[-1] - nc[0]) for _, nc, _ in corpus])
    allg = np.concatenate([g for _, _, g in corpus])
    fin = np.isfinite(allg)
    # additive dN histogram
    vals, cnts = np.unique(dN, return_counts=True)
    ax.bar(vals - 0.15, cnts, width=0.3, color=C_TRUE, label="additive ΔN (per run)")
    gv, gc = np.unique(np.round(allg[fin], 1), return_counts=True)
    ax.bar(gv + 0.15, gc, width=0.3, color=C_ALL, label="ratio gain (finite entries)")
    ax.axvline(0, color=INK2, lw=0.8, ls=":")
    ax.set_xlabel("value")
    ax.set_ylabel("count")
    ax.set_title("(b) additive ΔN has support {−1,0,+1};\nratio gain only {0, 1.000}",
                 loc="left")
    ax.set_xticks([-1, 0, 1])
    ax.legend(fontsize=7.5, loc="upper right")

    # ---- Panel C: the off-truth pilot ----
    if have_pilot:
        ax = axs[2]
        arms = {"truth": [], "fix_mc": []}
        for f in pilot:
            d = np.load(f, allow_pickle=True)
            if bool(d["scrambled"]):
                continue
            sm = str(d["start_mode"])
            arms.setdefault(sm, []).append(
                (int(d["C0_start"]), int(d["C_fix"]), float(d["dN_final"]),
                 float(np.asarray(d["traj_margin_flow"]).sum())))
        for k, (sm, col, lab) in enumerate([("truth", INK2, "arm (a) truth"),
                                            ("fix_mc", C_FIX, "arm (b-FIX) 1σ(mc)")]):
            vals = arms.get(sm, [])
            if not vals:
                continue
            dNs = [v[2] for v in vals]
            ax.scatter(np.full(len(dNs), k) + np.linspace(-0.12, 0.12, len(dNs)),
                       dNs, color=col, s=42, label=lab, zorder=3)
        ax.axhline(0, color=INK2, lw=0.8, ls=":")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["truth", "b-FIX 1σ(mc)"])
        ax.set_ylabel("ΔN = |C_fix| − |C_0(start)|")
        ax.set_title("(c) pilot: does a known mc displacement\nget repaired?", loc="left")
        ax.legend(fontsize=7.5, loc="best")
    else:
        # annotate the placeholder into panel B region — nothing to draw
        pass

    fig.suptitle("KINDLE — the loop-gain statistic was degenerate; the ignition line is "
                 "measured off-truth (F15/F11-KINDLE)", fontsize=10.5, x=0.5, y=1.02)
    fig.tight_layout()
    out = f"{REPO}/reports/KINDLE_gain_contour.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}  (pilot panel: {'YES' if have_pilot else 'pending GPU'})")


if __name__ == "__main__":
    main()
