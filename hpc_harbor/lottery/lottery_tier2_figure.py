"""LOTTERY tier-2 calibration figure: certified count vs n_active (the proxy's axis),
floor overlaid, grades marked, basis-control grade-flip highlighted."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BANK = os.path.abspath(os.path.join(HERE, "..", "..", "reports", "lottery_tier2.npz"))
FIG  = os.path.abspath(os.path.join(HERE, "..", "..", "LOTTERY_tier2_calibration.png"))

d = np.load(BANK, allow_pickle=True)
rows = json.loads(str(d["rows_json"]))
basis = json.loads(str(d["basis_json"]))

GC = {"CONFIRMED": "#1b7837", "MARGINAL": "#e08214", "below": "#b2182b"}
def is_dgx(r): return r.get("carried")
def is_hgx_ctrl(r): return r.get("is_control")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.5, 5.6), constrained_layout=True,
                               gridspec_kw={"width_ratios": [1.7, 1]})

# ---- LEFT: count vs n_active, common-basis (hgx) cells + carried dgx anchors ----
for r in rows:
    x = r["n_active"]; y = r["corr"]; g = r["grade"]
    mk = "s" if is_hgx_ctrl(r) else "o"
    axL.errorbar(x, y, yerr=abs(r["corr"] - r["corr_lo"]), fmt=mk, ms=11, color=GC[g],
                 mec="k", mew=0.8, ecolor="0.5", capsize=3, zorder=3)
    lab = r["tag"].replace("_lit", "").replace(" (hgx)", "*")
    axL.annotate(lab, (x, y), textcoords="offset points", xytext=(7, 5), fontsize=7.5)
# carried dgx anchors (open diamonds)
dgx = [dict(tag="m1e03", n_active=23, corr=0.70, corr_lo=0.43, grade="below"),
       dict(tag="m2e03", n_active=30, corr=2.77, corr_lo=2.07, grade="CONFIRMED"),
       dict(tag="m3e07", n_active=109, corr=4.07, corr_lo=3.60, grade="CONFIRMED")]
for r in dgx:
    axL.plot(r["n_active"], r["corr"], "D", ms=12, mfc="none",
             mec=GC[r["grade"]], mew=1.8, zorder=2)

axL.axhline(1.0, color="k", ls="--", lw=1.2, label="the >1 switch-on bar")
axL.axvline(30, color="#3b3bcf", ls=":", lw=1.6, label=r"proxy switch $n_{\rm active}=30$")
axL.set_xlabel(r"$n_{\rm active}$  (channel budget — the analytic proxy's axis)")
axL.set_ylabel("certified count per realisation  (at refit floor)")
axL.set_title("(A) certified count vs channel budget — the proxy is NOT a monotonic switch\n"
              "filled = re-derived here (hgx03, one basis); open ◇ = CHORUS banked (dgx); "
              "* = basis-control", fontsize=9.5)
axL.set_xlim(19, 37)
from matplotlib.lines import Line2D
leg = [Line2D([], [], marker="o", ls="", mfc=GC["CONFIRMED"], mec="k", label="CONFIRMED"),
       Line2D([], [], marker="o", ls="", mfc=GC["MARGINAL"], mec="k", label="MARGINAL"),
       Line2D([], [], marker="o", ls="", mfc=GC["below"], mec="k", label="below"),
       Line2D([], [], ls="--", color="k", label="the >1 bar"),
       Line2D([], [], ls=":", color="#3b3bcf", label=r"proxy switch ($n_{\rm active}{=}30$)")]
axL.legend(handles=leg, fontsize=8, loc="upper left", framealpha=0.9)
axL.annotate("m3 e=0.25: MOST channels (34)\nbut HIGHEST floor (13.6) → only MARGINAL",
             (34, 1.50), xytext=(28.6, 3.4), fontsize=7.5,
             arrowprops=dict(arrowstyle="->", color="0.4"))

# ---- RIGHT: the basis control — same cell, two hosts ----
labels = [b["tag"].replace("_lit", "") for b in basis]
xs = np.arange(len(basis))
w = 0.36
axR.bar(xs - w/2, [b["bank_floor"] for b in basis], w, label="floor dgx (banked)",
        color="#8073ac", edgecolor="k")
axR.bar(xs + w/2, [b["mine_floor"] for b in basis], w, label="floor hgx (here)",
        color="#b2abd2", edgecolor="k", hatch="//")
for i, b in enumerate(basis):
    axR.plot(xs[i]-w/2, b["bank_count"], "kD", ms=8)
    axR.plot(xs[i]+w/2, b["mine_count"], "kd", ms=8, mfc="w")
    axR.annotate((f"dgx:{b['bank_grade']}\nhgx:{b['mine_grade']}") if b['bank_grade']!=b['mine_grade']
                 else b['bank_grade'], (xs[i], 5.4),
                 ha="center", fontsize=7.5,
                 color=("#b2182b" if b['bank_grade']!=b['mine_grade'] else "0.3"),
                 fontweight=("bold" if b['bank_grade']!=b['mine_grade'] else "normal"))
axR.set_ylim(0, 13)
axR.axhline(1.0, color="k", ls="--", lw=1, label="the >1 bar (count)")
axR.set_xticks(xs); axR.set_xticklabels(labels)
axR.set_ylabel("floor (bars, nat)  ·  count (◇ dgx / ◇ hgx)")
axR.set_title("(B) BASIS CONTROL: same cell, two hosts\n"
              "m2e03 robust; m1e03 floor −1.4 nat → grade FLIPS below→CONFIRMED", fontsize=9.5)
axR.legend(fontsize=7.5, loc="upper left")

fig.suptitle("LOTTERY tier-2 — validating the channel-budget switch-on proxy against certified-count scoring "
             "(criterion-v2.2, floors refit per cell)", fontsize=11.5)
fig.savefig(FIG, dpi=130)
print(f"figure -> {FIG}")
