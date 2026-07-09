"""LANES figures: (1) lane ladder vs e over the F2 gap; (2) K_eff contour over (e,Mc,f)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lanes_eccentric_ladder import (lane_ladder, harmonic_spectrum, TOL_LOOSEST,
                                     FLOAT_CEIL, fundamental_match_decay, yr)

ECC = [0.1, 0.3, 0.5, 0.7, 0.9]

# ================= FIG 1: lane ladder vs e over F2 gap =================
fig, ax = plt.subplots(figsize=(9, 6))
# F2 gap bands
ax.axhspan(TOL_LOOSEST, FLOAT_CEIL, color="crimson", alpha=0.08, zorder=0)
ax.axhline(FLOAT_CEIL, color="crimson", lw=1.6, ls="--")
ax.axhline(TOL_LOOSEST, color="darkorange", lw=1.6, ls="--")
ax.text(0.955, FLOAT_CEIL*1.05, "float ceiling 0.05\n(blind-search floor)",
        color="crimson", ha="right", va="bottom", fontsize=9)
ax.text(0.955, TOL_LOOSEST*0.72, "F2 loosest rung 1.85e-3", color="darkorange",
        ha="right", va="top", fontsize=9)
ax.text(0.12, np.sqrt(TOL_LOOSEST*FLOAT_CEIL), "THE GAP\n(27x, 0 rungs)",
        color="crimson", ha="center", va="center", fontsize=11, weight="bold", alpha=0.75)

for i, e in enumerate(ECC):
    Lp = lane_ladder(e, "power")   # power-weighted usable band
    Ls = lane_ladder(e, "snr")     # residual-SNR-weighted usable band
    x = e
    # power-weighted rungs (open markers)
    ax.scatter([x]*len(Lp["rung"]), Lp["rung"], s=22, facecolors="none",
               edgecolors="steelblue", zorder=5,
               label="power-weighted usable rungs" if i == 0 else None)
    ax.plot([x, x], [Lp["rung"].min(), Lp["rung"].max()], color="steelblue", lw=1, alpha=0.5)
    # snr-weighted widest lane (filled star = the physically-detected widest lane)
    ax.scatter([x], [Ls["widest"]], marker="*", s=160, color="green", zorder=6,
               label="residual-SNR widest lane (physical)" if i == 0 else None)
    ax.annotate(f"n_pk={Lp['n_peak']}", (x, Lp["rung"].max()), textcoords="offset points",
                xytext=(6, 3), fontsize=8, color="steelblue")

ax.set_yscale("log")
ax.set_xlabel("eccentricity e")
ax.set_ylabel("registration lane tolerance (scaled source-sky units)")
ax.set_title("Eccentric harmonic lane ladder vs the F2 gap\n"
             "(anchor A2: power-dominant harmonic at f_gw = F2's 1.85e-3 rung)")
ax.set_ylim(3e-4, 8e-2)
ax.set_xlim(0.0, 1.0)
ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
ax.grid(True, which="both", alpha=0.2)
fig.tight_layout()
fig.savefig("lanes_ladder.png", dpi=130)
print("saved lanes_ladder.png")

# ================= FIG 2: K_eff contour over (e, Mc, f) =================
# K_eff ~ pi / dphi_cycle  (resolvable fundamental fringes before pulsar-lag decoherence)
Mc_panels = [1e8, 1e9, 5e9]
e_ax = np.linspace(0.05, 0.9, 40)
f_ax = np.logspace(-9.2, -8.0, 40)   # f_orb Hz
tau = 3e3 * yr

fig2, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
for k, Mc in enumerate(Mc_panels):
    Keff = np.zeros((len(f_ax), len(e_ax)))
    for i, f in enumerate(f_ax):
        for j, e in enumerate(e_ax):
            d = fundamental_match_decay(Mc, f, e, tau)
            dphi = max(d["dphi_cycle"], 1e-12)
            Keff[i, j] = np.pi / dphi
    Keff = np.clip(Keff, 0.3, 1e4)
    ax = axes[k]
    cs = ax.contourf(e_ax, f_ax, np.log10(Keff), levels=20, cmap="viridis")
    # K_eff = 1 and K_eff = 3 boundaries
    for lev, col, lab in [(1.0, "red", "K_eff=1"), (3.0, "white", "K_eff=3")]:
        c = ax.contour(e_ax, f_ax, Keff, levels=[lev], colors=col, linewidths=2)
        ax.clabel(c, fmt=lab, fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("eccentricity e")
    ax.set_title(f"Mc = {Mc:.0e} Msun  (tau_p = 3 kyr)")
    if k == 0:
        ax.set_ylabel("orbital fundamental f_orb (Hz)")
    cb = fig2.colorbar(cs, ax=ax)
    cb.set_label("log10 K_eff (fundamental)")
fig2.suptitle("Fundamental-lane coherence K_eff = pi/dphi_cycle over (e, Mc, f_orb): "
              "evolution kills the fundamental lane only in the high-Mc/high-f/high-e corner",
              fontsize=11)
fig2.tight_layout()
fig2.savefig("lanes_kcontour.png", dpi=130)
print("saved lanes_kcontour.png")
