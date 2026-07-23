#!/usr/bin/env python3
"""GLACIER g2a-ii forensic figure -- WHY the fit gate pegged, and what un-pegs it.

Three banks (GLACIER_results/, lane hgx03_NVIDIAH200, zero-noise gate-class):
  g2forensic_T15_n32   -- the smoke venue that failed (incumbent band 10-32 nHz)
  g2forensic_T30_n256  -- the main-gate venue: same reading (H-FAINT)
  g2forensic_T30_n256_flom8p7 -- the REMEDY-A study (band extended to ~2 nHz)

L: the profiles. Population AND the true-A0 control are monotone from the grid floor at
   both incumbent venues -- the array carries no information on the in-band amplitude.
M: the decomposition at T=30: the quadratic (signal) term is FLAT against the Occam term
   at every amplitude -- there is nothing to fit, not a mis-set tolerance.
R: remedy A -- the same fit at the band the array actually hears.
CPU, no GPU. House palette (glacier_g2_figure.py).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("GLACIER_OUT",
                     os.path.abspath(os.path.join(_HERE, "..", "..", "GLACIER_results")))
OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "GLACIER_g2aii_forensic.png"))
LANE = "hgx03_NVIDIAH200"

BLUE, RED, INK, GREY = "#4a7ab5", "#c1432e", "#2a2a2a", "#8a8a8a"


def bank(stem):
    return np.load(f"{RES}/{stem}__{LANE}.npz", allow_pickle=True)


z15 = bank("g2forensic_T15_n32")
z30 = bank("g2forensic_T30_n256")
zA = bank("g2forensic_T30_n256_flom8p7")

fig, (axL, axM, axR) = plt.subplots(1, 3, figsize=(15.5, 4.5))

# --- L: profiles at the incumbent band -- everything monotone from the floor ---
for z, ls, tag in ((z15, "--", "T=15, n=32"), (z30, "-", "T=30, n=256")):
    g = z["grid"]
    axL.plot(g, z["lnl_pop"] - z["lnl_pop"][0], ls, color=BLUE, lw=3.2, alpha=0.55,
             label=f"population ({tag})")
    axL.plot(g, z["lnl_ctrl"] - z["lnl_ctrl"][0], ls, color=RED, lw=1.3,
             label=f"TRUE-A0 control ({tag})")
axL.axvline(float(z30["a_eff_log10"]), color=INK, lw=0.9, ls=":")
axL.text(float(z30["a_eff_log10"]) + 0.04, -300, r"$A_{\rm eff}$ (drawn)", rotation=90,
         fontsize=8, color=INK, va="bottom")
axL.set_xlabel(r"$\log_{10} A_{\rm background}$")
axL.set_ylabel(r"$\Delta\ln L$ from the grid floor")
axL.set_title("Incumbent band (10-32 nHz):\neven a TRUE NG15 background pegs at the floor")
axL.legend(frameon=False, fontsize=8, loc="lower left")
axL.grid(alpha=0.25)

# --- M: the decomposition at the main-gate venue ---
g = z30["grid"]
axM.semilogy(g, np.maximum(z30["quad_pop"] - z30["quad_pop"][0], 1e-3), "-",
             color=BLUE, lw=1.6, label="signal term rise (population)")
axM.semilogy(g, np.maximum(z30["quad_ctrl"] - z30["quad_ctrl"][0], 1e-3), "-",
             color=RED, lw=1.6, label="signal term rise (control)")
axM.semilogy(g, np.maximum(z30["occam_pop"] - z30["occam_pop"][0], 1e-3), "-",
             color=GREY, lw=1.6, label="Occam term rise")
axM.set_xlabel(r"$\log_{10} A_{\rm background}$")
axM.set_ylabel("rise from the grid floor (nats)")
axM.set_title("T=30, n=256, decomposed:\nno signal term to fit -- the venue is DEAF in-band")
axM.legend(frameon=False, fontsize=8, loc="upper left")
axM.grid(alpha=0.25, which="both")

# --- R: remedy A -- the band the array hears ---
gA = zA["grid"]
axR.plot(gA, zA["lnl_pop"] - zA["lnl_pop"].max(), "-", color=BLUE, lw=1.6,
         label="population")
axR.plot(gA, zA["lnl_ctrl"] - zA["lnl_ctrl"].max(), "-", color=RED, lw=1.6,
         label="TRUE-A0 control")
aeffA, ahatA, sigA = float(zA["a_eff_log10"]), float(zA["ahat_pop"]), float(zA["sig_pop"])
axR.axvline(aeffA, color=INK, lw=0.9, ls=":")
if np.isfinite(sigA):
    axR.axvspan(ahatA - sigA, ahatA + sigA, color=BLUE, alpha=0.12)
axR.set_xlabel(r"$\log_{10} A_{\rm background}$")
axR.set_ylabel(r"$\Delta\ln L$ from the peak")
axR.set_ylim(-40, 3)
edge = bool(zA["edge_pop"])
axR.set_title(f"REMEDY A (band down to 2 nHz), T=30, n=256:\n"
              f"$\\hat A$ = {ahatA:.2f} $\\pm$ {sigA:.3f} vs $A_{{\\rm eff}}$ = "
              f"{aeffA:.2f}{'  [STILL EDGE]' if edge else ''}")
axR.legend(frameon=False, fontsize=8, loc="lower left")
axR.grid(alpha=0.25)

fig.suptitle("GLACIER g2a-ii forensic: the fit gate pegged because the generative band is "
             "outside the array's hearing -- not a bug, and not fixable by tolerance",
             fontsize=11, y=1.03)
fig.tight_layout()
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print(f"wrote {OUT}")
