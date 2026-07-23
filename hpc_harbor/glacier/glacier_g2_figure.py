#!/usr/bin/env python3
"""GLACIER Stage-0 figure — the census power hierarchy and the arithmetic drain curve.

Reads reports/glacier_g2_population.npz (banked by glacier_population.py). CPU, no GPU.
This is the ZERO-NOISE drain that EXISTS to be measured before the loop runs — NOT a loop result.
The campaign's headline drain-vs-knee figure is DEFERRED with the wide launch (GLACIER_capstone.md).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.abspath(os.path.join(_HERE, "..", "..", "reports", "glacier_g2_population.npz"))
OUT = os.path.abspath(os.path.join(_HERE, "..", "..", "GLACIER_g2_drain.png"))

d = np.load(NPZ, allow_pickle=True)
P = d["P16_residual"]; is_loud = d["is_loud16"]; Abg = d["drain_Abg"]
A_total = float(d["A_total"]); loud_frac = float(d["loud_drain_frac"])

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.5, 4.4))

# --- left: per-source power hierarchy (loud vs faint reservoir) ---
order = np.argsort(-P)
Ps = P[order]; loud_s = is_loud[order]
x = np.arange(len(Ps))
axL.bar(x[loud_s], Ps[loud_s] / Ps.sum(), color="#c1432e", label="loud (3)")
axL.bar(x[~loud_s], Ps[~loud_s] / Ps.sum(), color="#4a7ab5", label="faint reservoir (13)")
axL.set_yscale("log")
axL.set_xlabel("source (power-ranked)")
axL.set_ylabel("fraction of total residual power")
axL.set_title(f"Census power hierarchy\nloud carry {loud_frac:.1%}; reservoir = {1-loud_frac:.1%}")
axL.legend(frameon=False, fontsize=9)

# --- right: the arithmetic drain curve A_background as sources resolve ---
n = np.arange(len(Abg))
axR.plot(n, Abg / A_total, "-o", color="#2a2a2a", ms=3.5, lw=1.4)
axR.axvspan(0, 3, color="#c1432e", alpha=0.10)
axR.text(1.5, 0.55, "loud\nresolve", ha="center", va="center", fontsize=8, color="#c1432e")
axR.set_xlabel("N_resolved (sources fed, loudest-power first)")
axR.set_ylabel(r"$A_{\rm background}\,/\,A_{\rm total}$")
axR.set_title("The arithmetic drain (frozen census, zero-noise)\nconservation exact; drain -> 0 when all fed")
axR.set_ylim(-0.03, 1.03)
axR.grid(alpha=0.25)
# annotate the residual floor after the loud tier
after_loud = Abg[3] / A_total
axR.annotate(f"after loud tier: {after_loud:.1%} of A_total\nremains in the reservoir",
             xy=(3, after_loud), xytext=(6.5, 0.35), fontsize=8,
             arrowprops=dict(arrowstyle="->", color="#555"))

fig.suptitle("GLACIER Stage-0 (g2): the drain that EXISTS to be measured — before the loop",
             fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(OUT, dpi=130, bbox_inches="tight")
print(f"wrote {OUT}")
