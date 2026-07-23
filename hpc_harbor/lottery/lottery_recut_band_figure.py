"""LOTTERY re-cut band figure: the m1 e=0.3 grade as a band over the GWB-basis ensemble,
with m2 e=0.3 as the far-from-bar contrast that proves the rule is not just noise-widening."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
BANK = os.path.abspath(os.path.join(HERE, "..", "..", "reports", "lottery_recut_band.npz"))
FIG = os.path.abspath(os.path.join(HERE, "..", "..", "LOTTERY_recut_band.png"))

d = np.load(BANK, allow_pickle=True)
rows = json.loads(str(d["rows_json"]))
summary = {s["cell"]: s for s in json.loads(str(d["summary_json"]))}

GC = {"CONFIRMED": "#1b7837", "MARGINAL": "#e08214", "below": "#b2182b"}
CELLS = ["m1e03_lit", "m2e03_lit"]
TITLES = {"m1e03_lit": "m1 e=0.30  —  THE RECUT HEADLINE CELL\n(\"a single e=0.3 member does not switch the count on\")",
          "m2e03_lit": "m2 e=0.30  —  far-from-bar CONTRAST\n(the control that shows the rule is not just noise-widening)"}

fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.4), constrained_layout=True,
                         gridspec_kw={"width_ratios": [1, 1.25]})

for ci, cell in enumerate(CELLS):
    rs = [r for r in rows if r["cell"] == cell]
    s = summary[cell]
    ys = np.arange(len(rs))[::-1]
    labs = [f"{r['basis']}  {r['host'].split()[0]}/{r['cpus']}thr" for r in rs]

    # ---- LEFT: floor per basis state, with its own bootstrap error ----------
    ax = axes[ci][0]
    for y, r in zip(ys, rs):
        ax.errorbar(r["floor"], y, xerr=r["err"], fmt="o", ms=8, color=GC[r["grade_band"]],
                    mec="k", mew=0.7, ecolor="0.45", capsize=4, zorder=3)
    ax.axvspan(s["floor_min"], s["floor_max"], color="0.85", zorder=0,
               label=f"across-basis spread = {s['floor_spread']:.2f} nat")
    ax.set_yticks(ys); ax.set_yticklabels(labs, fontsize=8)
    ax.set_xlabel("refit null floor (nat)")
    ax.set_title(f"floor vs GWB basis   —   spread / mean error = {s['spread_over_err']:.2f}",
                 fontsize=9.5)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.grid(axis="x", alpha=0.25)

    # ---- RIGHT: certified count band vs the >1 bar --------------------------
    ax = axes[ci][1]
    for y, r in zip(ys, rs):
        lo, hi = r["count_lo"], r["count_hi"]
        ax.plot([lo, hi], [y, y], "-", lw=6, color=GC[r["grade_band"]], alpha=0.35,
                solid_capstyle="butt", zorder=2)
        ax.plot(r["count"], y, "o", ms=8, color=GC[r["grade_band"]], mec="k", mew=0.7, zorder=3)
        ax.annotate(r["grade_band"], (max(hi, r["count"]), y), textcoords="offset points",
                    xytext=(9, -3), fontsize=8, color=GC[r["grade_band"]], fontweight="bold")
    ax.axvline(1.0, color="k", ls="--", lw=1.4, zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels([])
    ax.set_xlabel("certified count per realisation   (bar spans floor∓err)")
    straddle = ("count band STRADDLES the bar" if s["straddles_bar"]
                else "count band CLEARS the bar in every basis state")
    ax.set_title(f"grade band   —   ENSEMBLE: {s['ensemble_grade']}\n{straddle}", fontsize=9.5,
                 color=("#b2182b" if s["straddles_bar"] else "#1b7837"))
    ax.set_xlim(0, max(3.9, max(r["count_hi"] for r in rs) * 1.35))
    ax.set_ylim(ys.min() - 0.7, ys.max() + 0.7)
    ax.grid(axis="x", alpha=0.25)
    # inside the axes, hugging the dashed bar -- keeps clear of the panel title
    ax.text(0.97, ys.min() - 0.5, "the >1 switch-on bar", fontsize=7.5, ha="right",
            va="center", color="0.25")

    axes[ci][0].text(-0.34, 0.5, TITLES[cell], transform=axes[ci][0].transAxes,
                     rotation=90, va="center", ha="center", fontsize=8.5, fontweight="bold")

leg = [Line2D([], [], marker="o", ls="", mfc=GC[g], mec="k", label=g)
       for g in ("CONFIRMED", "MARGINAL", "below")]
fig.legend(handles=leg, fontsize=8.5, ncol=3, loc="lower center", frameon=False,
           bbox_to_anchor=(0.5, -0.012))
fig.suptitle("LOTTERY re-cut — the m1 e=0.3 (lit) grade quoted as a BAND over a GWB-basis ensemble\n"
             "pre-registered two-sided rule: CONFIRMED iff count(floor+err)>1 · below iff count(floor−err)≤1 · "
             "else MARGINAL (the bar lies inside 1σ)", fontsize=11)
fig.savefig(FIG, dpi=130, bbox_inches="tight")
print(f"figure -> {FIG}")
