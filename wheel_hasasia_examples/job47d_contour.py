# ================= JOB 47d -- REFINED 0.9 CONTOUR ON THE 47c MAP ===============
# The map's CELLS STAY RAW (no smoothing of displayed values). Changes:
#   1. Every cell within ONE CELL of the current (47c) 0.9 contour is re-measured
#      at NREAL_HI >= 1000 realisations (same protocol, same seed formula, more
#      trials); all other cells keep their 47c values (250 realisations).
#   2. The 0.9 CONTOUR ONLY is traced through the field smoothed with a 3x3 MEAN
#      (uniform filter, nearest-edge padding) -- stated here and in the report.
#      Displayed cells are the raw (unsmoothed) per-cell fractions.
# Caption addition: the sigma_d below which the prior window (6 sigma_d) is
# narrower than the coarsest tooth DeltaL, i.e. K = 6 sigma_d/DeltaL < 1 -- the
# left edge of the map is prior-only there.
import os, json, time, csv
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter

BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
os.chdir(BASE)
t_wall = time.time()

# ---- machine namespace: Job 47b/c fast loader, verbatim -------------------------
NS = {"__name__": "__main__"}
nb = json.load(open("active_wheel_demo_v4.ipynb"))
cells = {i: "".join(c["source"]) for i, c in enumerate(nb["cells"])
         if c["cell_type"] == "code"}

def ex(i, lo=None, hi=None, tag=""):
    src = cells[i]
    if lo is not None:
        src = "\n".join(src.splitlines()[lo-1:hi])
    exec(compile(src, f"<v4-cell-{i}{':'+tag if tag else ''}>", "exec"), NS)

for i in (1, 3, 5, 7, 9, 11):
    ex(i)
ex(14, 116, 116); ex(14, 168, 168); ex(14, 780, 786); ex(14, 1586, 1591)
ex(14, 1072, 1085); ex(14, 3023, 3029); ex(14, 4043, 4049); ex(14, 4154, 4201)
ex(15, 13, 17); ex(17, 1, 127); ex(18, 12, 17); ex(27, 67, 78)
print("  machine namespace loaded (sliced fast path)")
_pop18 = NS["_pop18"]
C_PC_S = NS["C_PC_S"]
TOOTH_W = 0.25

lines43 = open("job43c_rebooking.py").read().splitlines()
def run_span(lo, hi):
    exec(compile("\n".join(lines43[lo-1:hi]), f"<43c:{lo}-{hi}>", "exec"), NS)
run_span(485, 509); run_span(661, 689); run_span(712, 756)
inv_rate = NS["inv_rate"]

r = NS["book"](NS["collect36m"](fuzzy=True, full_dmu=True, marg=True))
print(f"  ANCHOR CHECK: {r['y']:.2f} @ {r['pur']:.4f} -- booked 36.27 @ 0.9009")
assert abs(r["y"] - 36.27) < 0.5 and abs(r["pur"] - 0.9009) < 0.005, "anchor moved -- STOP"

SETS = []
for j in range(50):
    sd, f2, h2, P0, sea = _pop18(j, 0.35)
    rres = inv_rate(f2, h2, P0, sea, 45.0)
    SETS.append(None if rres is None else rres[2])   # fk sorted by nats
print(f"  comb sets: {sum(s is not None for s in SETS)}/50 skies")

NCOMB = 4
def purity_cell(sigd, sigphi, nreal, seed0=470_000_000):
    right = 0; dl = []
    for t in range(nreal):
        fk = SETS[t % 50]
        if fk is None:
            fk = SETS[(t+1) % 50]
        fk = fk[:NCOMB]
        rg = np.random.default_rng(seed0 + t*7919 + int(1e6*sigd) % 99991
                                   + int(1e6*sigphi) % 9973)
        mu = np.arccos(rg.uniform(-0.9, 0.9, len(fk)))
        kap = 2*np.pi*fk*(1-np.cos(mu))/C_PC_S
        dl.append(2*np.pi/kap.min())
        Lt = float(np.clip(rg.normal(0, sigd), -2.9*sigd, 2.9*sigd))
        Lg = np.linspace(-3*sigd, 3*sigd, 601)
        eps = kap*Lt + rg.normal(0, sigphi, len(kap))
        per = (kap[:, None]*Lg[None, :] - eps[:, None] + np.pi) % (2*np.pi) - np.pi
        S = -(per**2).sum(0)/(2*sigphi**2) - Lg**2/(2*sigd**2)
        right += int(abs(float(Lg[int(np.argmax(S))]) - Lt) <= TOOTH_W)
    return right/nreal, np.array(dl)

# ---- load the raw 47c grid ------------------------------------------------------
rows = list(csv.DictReader(open("job47c_grid.csv")))
XS = np.array(sorted({float(r["sigma_d_pc"]) for r in rows}))
YS = np.array(sorted({float(r["sigma_phi_rad"]) for r in rows}))
PUR = np.zeros((len(YS), len(XS)))
NRE = np.full(PUR.shape, 250, int)
for r_ in rows:
    PUR[np.argmin(np.abs(YS-float(r_["sigma_phi_rad"]))),
        np.argmin(np.abs(XS-float(r_["sigma_d_pc"])))] = float(r_["purity"])
print(f"  47c grid loaded: {PUR.shape[0]} x {PUR.shape[1]}")

# ---- cells within one cell of the current 0.9 contour ---------------------------
above = PUR >= 0.9
straddle = np.zeros_like(above)
straddle[:, :-1] |= above[:, :-1] != above[:, 1:]
straddle[:, 1:] |= above[:, :-1] != above[:, 1:]
straddle[:-1, :] |= above[:-1, :] != above[1:, :]
straddle[1:, :] |= above[:-1, :] != above[1:, :]
from scipy.ndimage import binary_dilation
band = binary_dilation(straddle, iterations=1)        # within one cell
print(f"  contour band: {int(band.sum())} of {band.size} cells -> re-measured"
      f" at >= 1000 realisations")

# ---- refine ---------------------------------------------------------------------
NREAL_HI = 1000
t0 = time.time()
dl_all = []
for iy, ix in zip(*np.where(band)):
    PUR[iy, ix], dl = purity_cell(XS[ix], YS[iy], NREAL_HI)
    NRE[iy, ix] = NREAL_HI
    dl_all.append(dl)
dl_all = np.concatenate(dl_all)
print(f"  refinement done ({time.time()-t0:.0f} s)")

# ---- contour from the 3x3-mean-smoothed field (contour ONLY; cells stay raw) ----
SM = uniform_filter(PUR, size=3, mode="nearest")

def cross_x(row, xs=XS, level=0.9):
    for a in range(len(xs)-1):
        if row[a] >= level >= row[a+1]:
            f = (row[a]-level)/max(row[a]-row[a+1], 1e-9)
            return float(10**(np.log10(xs[a]) + f*(np.log10(xs[a+1])-np.log10(xs[a]))))
    return np.nan

# ---- caption numbers ------------------------------------------------------------
dl_med = float(np.median(dl_all))
k1 = np.quantile(dl_all/6.0, [0.1, 0.5, 0.9])
print("\n  CAPTION NUMBERS:")
print(f"    coarsest-comb DeltaL median {dl_med:.2f} pc; K < 1 (prior window"
      f" narrower than the coarsest tooth) below sigma_d = DeltaL/6 ="
      f" {k1[1]:.2f} pc [10/90 over realisations: {k1[0]:.2f}, {k1[2]:.2f}]"
      f" -- the map's left edge is prior-only there")
for sp in (0.25, 0.5, 1.0, 1.5, 2.0):
    iy = np.argmin(np.abs(YS-sp))
    print(f"    smoothed-field 0.9 crossing near sigma_phi = {sp:4.2f}"
          f" (grid row {YS[iy]:.2f}): sigma_d = {cross_x(SM[iy]):.2f} pc")

# ---- CSV (per-cell n_real now varies) -------------------------------------------
with open("job47d_grid.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["sigma_d_pc", "sigma_phi_rad", "purity", "n_real"])
    for iy in range(len(YS)):
        for ix in range(len(XS)):
            cw.writerow([f"{XS[ix]:.4f}", f"{YS[iy]:.4f}",
                         f"{PUR[iy,ix]:.4f}", NRE[iy, ix]])
print("  wrote job47d_grid.csv")

# ---- figure (identical presentation to 47c; contour from SM) --------------------
fig, ax = plt.subplots(figsize=(6.8, 4.8))
xe = np.sqrt(XS[:-1]*XS[1:]); xe = np.concatenate([[XS[0]**2/xe[0]], xe, [XS[-1]**2/xe[-1]]])
ye = 0.5*(YS[:-1]+YS[1:]); ye = np.concatenate([[2*YS[0]-ye[0]], ye, [2*YS[-1]-ye[-1]]])
pc = ax.pcolormesh(xe, ye, PUR, cmap="Blues", vmin=0.0, vmax=1.0, shading="flat")
ax.set_xscale("log")
ax.set_xlim(xe[0], xe[-1]); ax.set_ylim(ye[0], ye[-1])
ax.axhline(2.0, color="0.25", lw=1.0, ls="--")
cs = ax.contour(XS, YS, SM, levels=[0.9], colors="0.15", linewidths=0.9)
edge_hit = any(len(p.vertices) and ((p.vertices[:, 0] <= XS[0]*1.0001).any()
               or (p.vertices[:, 0] >= XS[-1]*0.9999).any())
               for p in cs.get_paths())
if edge_hit:
    cs.remove()
    print("  0.9 contour runs off a side edge -> OMITTED")
else:
    print("  0.9 contour (3x3-mean-smoothed field) closes cleanly -> DRAWN")
ax.set_xlabel(r"uncertainty on the pulsar distance prior $\sigma_d$ [pc]")
ax.set_ylabel(r"GW phase jitter $\sigma_\phi$ [rad]")
sax = ax.secondary_yaxis("right", functions=(lambda sp: sp, lambda sp: sp))
sticks = np.array([0.9, 0.7, 0.5, 0.3, 0.1])
sax.set_yticks(np.sqrt((1-sticks)/sticks))
sax.set_yticklabels([f"{s:g}" for s in sticks])
sax.set_ylabel("s = source's share of its frequency bin")
cb = fig.colorbar(pc, ax=ax, pad=0.13)
cb.set_label("fraction of trials where the correct distance is found")
plt.tight_layout()
plt.savefig("fig47d_criterion_heatmap.png", dpi=150, bbox_inches="tight")
print("  wrote fig47d_criterion_heatmap.png")
print(f"\n  TOTAL RUNTIME: {(time.time()-t_wall)/60:.1f} min (CPU-class)")
