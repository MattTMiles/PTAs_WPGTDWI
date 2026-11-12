#!/usr/bin/env python
"""PHOENIX -- THE STAGE-2a PLANE: the paper's closing figure.

Four layers over (sky loudness x array capability), one panel per igniter arm:

  L1 FEED     (first-bite)    -- sequential single-hue fill, 0-4 skies
  L2 HOLD     (M-step trust)  -- open circles, one per holding sky   [shape, not colour]
  L3 CERTIFY  (N_cert > 0)    -- star marker                          [shape, not colour]
  L4 SAFETY   (manufacture)   -- reserved status red, hatch + label; the wrong-cert
                                 boundary, plus the structure-limited shading that
                                 says WHY certification never arrives at the feasible
                                 rung (it is not a sensitivity limit)

  + the SUMMIT S2 open-question box over the unmeasured (loudness x capability) region.

Identity is never colour-alone: every contour after the first is a distinct SHAPE, and
the safety layer carries hatch + text as well as its status hue.

Reads reports/summit_s0_harvest.npz (96 cells, re-harvested 2026-07-29 over the fully
drained array ladder). CPU only.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

d = np.load("reports/summit_s0_harvest.npz", allow_pickle=True)

RUNGS = ["r13p9", "r13p5", "r13p25"]                     # y, bottom -> top (loudness up)
RUNG_LBL = ["r13p9\nNG15-feasible\n(0$\\sigma$)",
            "r13p5\n+7.9$\\sigma$",
            "r13p25\n+12.5$\\sigma$"]
CAPS = [("w1", 30), ("w1", 40), ("w0p5", 30), ("w0p5", 40), ("w0p25", 30), ("w0p25", 40)]
CAP_LBL = ["rms\nT=30\n(today)", "rms\nT=40", "rms/2\nT=30", "rms/2\nT=40",
           "rms/4\nT=30\n(DSA)", "rms/4\nT=40\n(DSA)"]
ARMS = [("none", "CIRCULAR igniter (e = 0)"), ("e07", "ECCENTRIC igniter (e = 0.7)")]

SEQ = plt.get_cmap("Blues")          # sequential, one hue, light -> dark (magnitude)
STATUS_RED = "#b2182b"               # reserved status: manufacture / wrong certs
INK, MUTED = "#1a1a1a", "#6b6b6b"


def cellstats(arm, rung, w, T):
    m = ((d["arm"] == arm) & (d["rung"] == rung) & (d["kind"] == "cell") &
         (d["tier"] == "lit") & (d["w"] == w) & (d["T"] == T))
    if not m.any():
        return None
    return dict(n=int(m.sum()),
                feed=int(d["first_bite"][m].sum()),
                hold=int(d["hold"][m].sum()),
                cert=int(d["cert"][m].sum()),
                wrong=int(d["wrong_tot"][m].sum()))


fig, axes = plt.subplots(1, 2, figsize=(15.2, 6.3), sharey=True)

for ax, (arm, arm_title) in zip(axes, ARMS):
    ax.set_title(arm_title, fontsize=12, fontweight="bold", color=INK, pad=10)
    for yi, rung in enumerate(RUNGS):
        for xi, (w, T) in enumerate(CAPS):
            s = cellstats(arm, rung, w, T)
            if s is None:                                   # unmeasured -> SUMMIT S2
                ax.add_patch(Rectangle((xi - .5, yi - .5), 1, 1, facecolor="#f2f2f2",
                                       edgecolor="#cfcfcf", hatch="//", lw=.8, zorder=1))
                continue
            frac = s["feed"] / s["n"]
            ax.add_patch(Rectangle((xi - .5, yi - .5), 1, 1,
                                   facecolor=SEQ(0.08 + 0.80 * frac),
                                   edgecolor="white", lw=2, zorder=2))
            # L1 the count, in ink (text never wears the series colour)
            ax.text(xi, yi + .30, f"{s['feed']}/{s['n']}", ha="center", va="center",
                    fontsize=11, fontweight="bold",
                    color="white" if frac > .55 else INK, zorder=5)
            # L2 HOLD -- one open circle per holding sky
            for j in range(s["hold"]):
                ax.plot(xi - .225 + .15 * j, yi - .02, "o", ms=7.5, mfc="none",
                        mec="white" if frac > .55 else INK, mew=1.7, zorder=5)
            # L3 CERTIFY -- star per certifying sky
            for j in range(s["cert"]):
                ax.plot(xi - .18 + .16 * j, yi - .30, "*", ms=13,
                        color="#f0c000", mec=INK, mew=.7, zorder=6)
            # L4 SAFETY -- wrong certs: status hue + hatch + count
            if s["wrong"]:
                ax.add_patch(Rectangle((xi - .5, yi - .5), 1, 1, facecolor="none",
                                       edgecolor=STATUS_RED, hatch="xx", lw=3, zorder=7))
                ax.text(xi + .34, yi + .34, f"{s['wrong']}", ha="center", va="center",
                        fontsize=9.5, fontweight="bold", color="white", zorder=9,
                        bbox=dict(boxstyle="circle,pad=.22", fc=STATUS_RED, ec="none"))

    ax.set_xlim(-.5, len(CAPS) - .5); ax.set_ylim(-.62, len(RUNGS) + .45)
    ax.set_xticks(range(len(CAPS))); ax.set_xticklabels(CAP_LBL, fontsize=8.5, color=MUTED)
    ax.set_xlabel("ARRAY CAPABILITY  $\\longrightarrow$", fontsize=10.5,
                  color=INK, fontweight="bold", labelpad=8)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)

axes[0].set_yticks(range(len(RUNGS)))
axes[0].set_yticklabels(RUNG_LBL, fontsize=9, color=MUTED)
axes[0].set_ylabel("SKY LOUDNESS  $\\longrightarrow$", fontsize=10.5,
                   color=INK, fontweight="bold", labelpad=8)

# ---- the two annotations that make it the closing figure --------------------
ax0, ax1 = axes
# (iv) the open question -- drawn FIRST so annotations sit on top of it
for ax in axes:
    ax.add_patch(Rectangle((0.5, 0.5), len(CAPS) - 1, 2, facecolor="none",
                           edgecolor="#9a9a9a", lw=1.6, ls=(0, (5, 3)), zorder=3))
ax1.text(3.55, 1.72, "SUMMIT §2  OPEN QUESTION\ncomposed cells\n"
                     "(loudness × capability never jointly run)",
         ha="center", va="center", fontsize=8.8, color="#5a5a5a", style="italic",
         fontweight="bold", zorder=9,
         bbox=dict(boxstyle="round,pad=.35", fc="white", ec="#c8c8c8", lw=.9, alpha=.93))

# (i) the manufacture boundary -- eccentric arm, above +8 sigma (in the headroom)
ax1.annotate("MANUFACTURE BOUNDARY — wrong certs, $q$=1.000, on_true=FALSE;\n"
             "eccentric arm ONLY, tension $\\geq$ +7.9$\\sigma$. 20/20 killed by the D2 rigidity gate.",
             xy=(0.0, 2.52), xytext=(0.55, 3.16), fontsize=8.6, color=STATUS_RED,
             fontweight="bold", ha="left", va="center",
             arrowprops=dict(arrowstyle="->", color=STATUS_RED, lw=1.8,
                             connectionstyle="arc3,rad=-.25"))

# (ii) certification is structure/venue-limited, not sensitivity-limited (headroom)
ax0.annotate("CERTIFICATION NEVER ARRIVES ALONG THE FEASIBLE ROW —\n"
             "and capability is not what stops it. 0/24 certs at every (rms, T);\n"
             "the brightest conditioned member lands at ~2.3 nHz, fringe-starved.\n"
             "STRUCTURE / VENUE-LIMITED, not sensitivity-limited.",
             xy=(2.6, 0.52), xytext=(-0.38, 3.05), fontsize=8.3, color=INK, ha="left",
             va="center",
             bbox=dict(boxstyle="round,pad=.40", fc="#fffbe6", ec="#e0c86a", lw=1.1),
             arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.4,
                             connectionstyle="arc3,rad=-.18"))

# (iii) the NEW crossing this re-harvest exposes (below the plane, in the margin)
ax0.annotate("NEW (2026-07-29 re-harvest): sky 3 crosses FEED at rms/2 and rms/4, T=30\n"
             "(3/4 $\\to$ 4/4) — the ONE sensitivity-driven first-bite crossing on the plane.",
             xy=(4.0, 0.34), xytext=(1.62, 1.42), fontsize=8.3, color="#0b5394",
             fontweight="bold", ha="left", va="center", annotation_clip=False,
             bbox=dict(boxstyle="round,pad=.34", fc="#eaf2fb", ec="#7ba7d7", lw=1.0),
             arrowprops=dict(arrowstyle="->", color="#0b5394", lw=1.6,
                             connectionstyle="arc3,rad=.22"))

handles = [
    Patch(fc=SEQ(0.88), ec="white", label="FEED — first-bite ON (fill = skies fed / 4)"),
    Line2D([], [], marker="o", ls="none", mfc="none", mec=INK, mew=1.7, ms=7.5,
           label="HOLD — M-step trust (one ring per holding sky)"),
    Line2D([], [], marker="*", ls="none", color="#f0c000", mec=INK, mew=.7, ms=13,
           label="CERTIFY — N_cert > 0 (one star per certifying sky)"),
    Patch(fc="none", ec=STATUS_RED, hatch="xx", lw=2, label="MANUFACTURE — wrong certs (count in badge)"),
    Patch(fc="#f2f2f2", ec="#cfcfcf", hatch="//", label="not measured — SUMMIT §2 territory"),
]
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=9,
           bbox_to_anchor=(.5, -.075), labelcolor=INK)

fig.suptitle("THE STAGE-2a PLANE — feed, hold, certify, and the safety boundary "
             "over sky loudness × array capability",
             fontsize=13.5, fontweight="bold", color=INK, y=1.005)
fig.text(.5, .945, "GLACIER Stage-2a dual ladder, 96 harvested cells, r13p9 array ladder "
                   "fully drained. HOLD = ≥2 consecutive refit intervals with every "
                   "commonly-fed member ≤0.05 dex in (log$_{10}f_{gw}$, log$_{10}m_c$).",
         ha="center", fontsize=8.6, color=MUTED)
fig.text(.5, .012,
         "Certification is scored under the ADOPTED criterion (v2.2). The single true-cert "
         "cell (r13p5/none, ★) certifies at i0 and i4 — non-consecutive — and survives a "
         "full-null floor re-cut\n(floor 0.000±0.03, 97% of nulls silent; the binding bar is "
         "the trials term lnK+0.578=1.27). A persistence rule (proposed v2.3, NOT adopted) "
         "would delete it. No drain layer is\nshown: the feed-free reference is offset and "
         "two-sided (LEDGER B1), so bite/regress are not quotable per-cell — every eccentric "
         "cell here fed nothing.",
         ha="center", fontsize=7.4, color=MUTED, style="italic")

fig.tight_layout(rect=[0, .155, 1, .90])
fig.savefig("reports/STAGE2A_PLANE.png", dpi=200, bbox_inches="tight",
            facecolor="white")
print("wrote reports/STAGE2A_PLANE.png")

# ---- the table view (accessibility: identity never colour-alone) ------------
print(f"\n{'arm':6s}{'rung':8s}{'w':7s}{'T':>4s}{'feed':>6s}{'hold':>6s}{'cert':>6s}{'wrong':>7s}")
for arm, _ in ARMS:
    for rung in RUNGS:
        for w, T in CAPS:
            s = cellstats(arm, rung, w, T)
            if s is None:
                print(f"{arm:6s}{rung:8s}{w:7s}{T:>4d}{'—':>6s}{'—':>6s}{'—':>6s}{'—':>7s}")
            else:
                print(f"{arm:6s}{rung:8s}{w:7s}{T:>4d}{s['feed']:>3d}/{s['n']:<2d}"
                      f"{s['hold']:>3d}/{s['n']:<2d}{s['cert']:>3d}/{s['n']:<2d}{s['wrong']:>7d}")
