"""LOTTERY tier-1 figure + validation-set selection (reads reports/lottery_tier1.npz)."""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
BANK = os.path.abspath(os.path.join(HERE, "..", "..", "reports", "lottery_tier1.npz"))
FIG  = os.path.abspath(os.path.join(HERE, "..", "..", "LOTTERY_tier1_mix.png"))

d = np.load(BANK, allow_pickle=True)
f_grid = d["f_grid"]; e_grid = d["e_grid"]
P_chan = d["P_chan"]; P_thr = d["P_thr"]; P_dis = d["P_dis"]
E_chan = d["E_chan"]; E_thr = d["E_thr"]
named = json.loads(str(d["named_json"]))

# extent for imshow: x=f_ecc, y=e_char ; origin lower
ext = [f_grid[0]-0.025, f_grid[-1]+0.025, e_grid[0]-0.025, e_grid[-1]+0.025]

fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), constrained_layout=True)

def panel(ax, Z, title, cmap, norm=None, vmin=None, vmax=None):
    im = ax.imshow(Z, origin="lower", aspect="auto", extent=ext, cmap=cmap,
                   norm=norm, vmin=vmin, vmax=vmax)
    ax.set_xlabel(r"$f_{\rm ecc}$  (fraction of loud members eccentric)")
    ax.set_ylabel(r"$e_{\rm char}$")
    ax.set_title(title, fontsize=11)
    return im

# A: channel-budget P
imA = panel(axes[0], P_chan, r"(A) $P$(switch-on) $-$ CHANNEL BUDGET  ($n_{\rm active}\geq30$)",
            "viridis", vmin=0, vmax=1)
cs = axes[0].contour(f_grid, e_grid, P_chan, levels=[0.5], colors="w", linewidths=2)
axes[0].clabel(cs, fmt="P=0.5", fontsize=8)
fig.colorbar(imA, ax=axes[0], shrink=0.9, label="P")

# B: threshold-form P
imB = panel(axes[1], P_thr, r"(B) $P$(switch-on) $-$ THRESHOLD  (any $e\geq0.5$ or pair $\geq0.3$)",
            "viridis", vmin=0, vmax=1)
cs = axes[1].contour(f_grid, e_grid, P_thr, levels=[0.5], colors="w", linewidths=2)
axes[1].clabel(cs, fmt="P=0.5", fontsize=8)
fig.colorbar(imB, ax=axes[1], shrink=0.9, label="P")

# C: disagreement (channel - threshold); always >=0 (threshold-ON subset of channel-ON)
imC = panel(axes[2], P_dis, r"(C) DISAGREEMENT  $P_{\rm chan}-P_{\rm thr}$  (proxy over-permissive)",
            "magma", vmin=0, vmax=max(0.5, float(P_dis.max())))
fig.colorbar(imC, ax=axes[2], shrink=0.9, label=r"$\Delta P$")

fig.suptitle("LOTTERY tier-1  --  mixed-population switch-on: channel-budget mechanism vs quoted threshold rule\n"
             "3+13 census, each loud member eccentric ($e{=}e_{\\rm char}$) w.p. $f_{\\rm ecc}$; "
             f"N={int(json.loads(str(d['conventions']))['N_DRAW_MIX'])} draws/cell",
             fontsize=12)
fig.savefig(FIG, dpi=130)
print(f"figure -> {FIG}")

# ---- th_not_ch structural check across mix -----------------------------------
print(f"\nmix: P_dis min={P_dis.min():.4f} max={P_dis.max():.4f}  "
      f"(min>=0 confirms threshold-ON subset of channel-ON)")

# ============================================================================
# VALIDATION-SET SELECTION
# ============================================================================
print("\n" + "=" * 78)
print("VALIDATION-SET CANDIDATES for tier-2")
print("=" * 78)

def mix_cell(fe, ec):
    j = int(np.argmin(np.abs(f_grid - fe))); i = int(np.argmin(np.abs(e_grid - ec)))
    return dict(kind="mix", f_ecc=float(f_grid[j]), e_char=float(e_grid[i]),
                P_chan=float(P_chan[i, j]), P_thr=float(P_thr[i, j]),
                P_dis=float(P_dis[i, j]))

# (i) disagreement band: largest P_dis, spread across the grid + the lnN peak=0.3 families
flat = [(P_dis[i, j], i, j) for i in range(len(e_grid)) for j in range(len(f_grid))]
flat.sort(reverse=True)
print("\n(i) DISAGREEMENT BAND -- top mix cells by (P_chan - P_thr):")
seen = []
for val, i, j in flat[:40]:
    fe, ec = float(f_grid[j]), float(e_grid[i])
    # spread: avoid near-duplicates
    if any(abs(fe-s[0])<0.14 and abs(ec-s[1])<0.14 for s in seen):
        continue
    seen.append((fe, ec))
    print(f"    f_ecc={fe:.2f} e_char={ec:.2f}  P_chan={P_chan[i,j]:.3f} "
          f"P_thr={P_thr[i,j]:.3f}  dP={val:.3f}")
    if len(seen) >= 6:
        break

# (ii) P ~ 0.5 (steep part) -- pick where P_chan closest to 0.5 (the proxy being validated)
print("\n(ii) STEEP PART -- mix cells with P_chan ~ 0.5:")
half = [(abs(P_chan[i, j]-0.5), i, j) for i in range(len(e_grid)) for j in range(len(f_grid))]
half.sort()
seenh = []
for val, i, j in half[:40]:
    fe, ec = float(f_grid[j]), float(e_grid[i])
    if any(abs(fe-s[0])<0.2 and abs(ec-s[1])<0.2 for s in seenh):
        continue
    seenh.append((fe, ec))
    print(f"    f_ecc={fe:.2f} e_char={ec:.2f}  P_chan={P_chan[i,j]:.3f} "
          f"P_thr={P_thr[i,j]:.3f}")
    if len(seenh) >= 4:
        break

# (iii) concordant anchors
print("\n(iii) CONCORDANT ANCHORS:")
# high: both ~1
hi = None; lo = None
for i in range(len(e_grid)):
    for j in range(len(f_grid)):
        if P_chan[i, j] > 0.98 and P_thr[i, j] > 0.98 and hi is None and f_grid[j] > 0.5:
            hi = (i, j)
# low: both ~0 but not the trivial all-circular edge; pick small f_ecc, low e
for i in range(len(e_grid)):
    for j in range(len(f_grid)):
        if P_chan[i, j] < 0.02 and P_thr[i, j] < 0.02 and f_grid[j] > 0.0 and lo is None:
            lo = (i, j)
for label, cell in [("concordant-HIGH", hi), ("concordant-LOW", lo)]:
    if cell:
        i, j = cell
        print(f"    {label}: f_ecc={f_grid[j]:.2f} e_char={e_grid[i]:.2f}  "
              f"P_chan={P_chan[i,j]:.3f} P_thr={P_thr[i,j]:.3f}")

# named-distribution disagreement leaders (for reference)
print("\nNAMED distributions ranked by disagreement:")
for r in sorted(named, key=lambda r: -r["P_disagree"])[:5]:
    print(f"    {r['name']:<22} P_chan={r['P_chan']:.3f} P_thr={r['P_thr']:.3f} "
          f"dis={r['P_disagree']:.3f}")
