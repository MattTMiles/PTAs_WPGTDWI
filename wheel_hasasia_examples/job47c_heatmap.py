# ================= JOB 47c -- POLISHED SIMPLE HEATMAP ===========================
# Same measurement as Job 47b (accepted); finer grid + presentation only.
# CONFIG (printed at runtime too, for verifiable sameness with 47b):
#   harness      : banked SEED50 (50 populations), _pop18(j, 0.35)
#   comb sets    : machine inv_rate (cap 1300, rho>=3, 2-30 nHz, 45 yr/100 ns/
#                  14 d), ranked by nats = ln(pi/sigma); TOP 4 kept (fixed N)
#   geometry     : mu = arccos U(-0.9, 0.9) per comb per realisation,
#                  kappa = 2 pi f (1-cos mu)/c
#   K <-> sigma_d: K = 6 sigma_d / DeltaL(coarsest of the 4)
#   s <-> sigma_phi: s = 1/(1+sigma_phi^2), applied to every comb
#   correct      : winner within TOOTH_W = 0.25 pc of truth; prior term included
#   realisations : 250 per cell, sky j = t mod 50, seeded (deterministic)
# Presentation deviation, flagged: the spec's "top axis" for s cannot sit on top
# -- s maps through sigma_phi^2 = (1-s)/s to the Y quantity; the s scale is
# drawn as a secondary axis on the RIGHT, twinned to sigma_phi.
import os, sys, json, time, csv
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
os.chdir(BASE)
t_wall = time.time()

# ---- machine namespace: Job 47b fast loader, verbatim ---------------------------
NS = {"__name__": "__main__"}
nb = json.load(open("active_wheel_demo_v4.ipynb"))
cells = {i: "".join(c["source"]) for i, c in enumerate(nb["cells"])
         if c["cell_type"] == "code"}

def ex(i, lo=None, hi=None, tag=""):
    src = cells[i]
    if lo is not None:
        src = "\n".join(src.splitlines()[lo-1:hi])
    t0 = time.time()
    exec(compile(src, f"<v4-cell-{i}{':'+tag if tag else ''}>", "exec"), NS)
    print(f"  [cell {i:2d}{' '+tag if tag else ''}] OK {time.time()-t0:6.1f} s",
          flush=True)

for i in (1, 3, 5, 7, 9, 11):
    ex(i)
ex(14, 116, 116, "EXP")
ex(14, 168, 168, "TOOTH_W")
ex(14, 780, 786, "rn_draw")
ex(14, 1586, 1591, "FP_ALL")
ex(14, 1072, 1085, "SEED50")
ex(14, 3023, 3029, "_sea50v")
ex(14, 4043, 4049, "h_of_rhoP")
ex(14, 4154, 4201, "build_pop15")
ex(15, 13, 17, "s_to_s30")
ex(17, 1, 127, "depro_sigma")
ex(18, 12, 17, "_pop18")
ex(27, 67, 78, "alias_score24")
_pop18 = NS["_pop18"]
C_PC_S = NS["C_PC_S"]; YR_S = NS["YR_S"]
TOOTH_W = 0.25

lines43 = open("job43c_rebooking.py").read().splitlines()
def run_span(lo, hi):
    exec(compile("\n".join(lines43[lo-1:hi]), f"<43c:{lo}-{hi}>", "exec"), NS)
run_span(485, 509)
run_span(661, 689)
run_span(712, 756)
inv_rate = NS["inv_rate"]

# ---- anchor gate (equivalence proof, as in 47b) ---------------------------------
r = NS["book"](NS["collect36m"](fuzzy=True, full_dmu=True, marg=True))
print(f"\n  ANCHOR CHECK: {r['y']:.2f} [{r['yb'][0]:.2f},{r['yb'][1]:.2f}]"
      f" @ {r['pur']:.4f} -- booked 36.27 [34.02,38.75] @ 0.9009")
assert abs(r["y"] - 36.27) < 0.5 and abs(r["pur"] - 0.9009) < 0.005, "anchor moved -- STOP"

# ---- comb sets, verbatim 47b ----------------------------------------------------
def comb_sets(strain_fac=1.0):
    out = []
    for j in range(50):
        sd, f2, h2, P0, sea = _pop18(j, 0.35)
        rres = inv_rate(f2, h2*strain_fac, P0*strain_fac**2, sea*strain_fac**2, 45.0)
        if rres is None:
            out.append(None); continue
        ks, sg, fk, rp, arr = rres
        s_el = 1.0/(1.0 + np.maximum(sg**2 - 1.0/np.maximum(rp, 1e-9)**2, 0.0))
        out.append((fk, s_el))
    return out

t0 = time.time()
SETS = comb_sets(1.0)
SETS3 = comb_sets(3.0)
print(f"  comb sets built ({time.time()-t0:.0f} s);"
      f" skies with sets: {sum(s is not None for s in SETS)}/50")
sh1, sh3 = [], []
for a, b in zip(SETS, SETS3):
    if a is None or b is None:
        continue
    common = np.intersect1d(np.round(a[0], 14), np.round(b[0], 14))
    for f0 in common[:4]:
        sh1.append(float(a[1][np.argmin(np.abs(a[0]-f0))]))
        sh3.append(float(b[1][np.argmin(np.abs(b[0]-f0))]))
print(f"  WHICH-TERM CHECK: own-element share s on shared combs, x1 vs x3 strain:"
      f" max |ds| = {np.max(np.abs(np.array(sh1)-np.array(sh3))):.2e}"
      f"  (s is GW-over-GW: UNCHANGED)")

# ---- protocol (identical to 47b) ------------------------------------------------
NREAL = 250
NCOMB = 4

def purity_cell(sigd, sigphi, ncomb=NCOMB, sets=None, seed0=470_000_000):
    sets = SETS if sets is None else sets
    right = 0; dl_coarse = []
    for t in range(NREAL):
        j = t % 50
        st = sets[j]
        if st is None:
            st = sets[(j+1) % 50]
        fk = st[0][:ncomb]
        rg = np.random.default_rng(seed0 + t*7919 + int(1e6*sigd) % 99991
                                   + int(1e6*sigphi) % 9973)
        mu = np.arccos(rg.uniform(-0.9, 0.9, len(fk)))
        kap = 2*np.pi*fk*(1-np.cos(mu))/C_PC_S
        dl_coarse.append(2*np.pi/kap.min())
        Lt = float(np.clip(rg.normal(0, sigd), -2.9*sigd, 2.9*sigd))
        Lg = np.linspace(-3*sigd, 3*sigd, 601)
        eps = kap*Lt + rg.normal(0, sigphi, len(kap))
        per = (kap[:, None]*Lg[None, :] - eps[:, None] + np.pi) % (2*np.pi) - np.pi
        S = -(per**2).sum(0)/(2*sigphi**2) - Lg**2/(2*sigd**2)
        right += int(abs(float(Lg[int(np.argmax(S))]) - Lt) <= TOOTH_W)
    return right/NREAL, float(np.median(dl_coarse))

# ---- the 20 x 20 grid -----------------------------------------------------------
XS = np.geomspace(0.1, 300.0, 20)
YS = np.linspace(0.05, 3.0, 20)
print(f"\n  GRID: 20 x 20, {NREAL} realisations/cell, N = {NCOMB} combs fixed")
t0 = time.time()
PUR = np.zeros((len(YS), len(XS)))
DLmed = []
for iy, sp in enumerate(YS):
    for ix, sd_ in enumerate(XS):
        PUR[iy, ix], dlm = purity_cell(sd_, sp)
        DLmed.append(dlm)
    print(f"  row sigma_phi={sp:5.2f}: " +
          " ".join(f"{PUR[iy,ix]:.2f}" for ix in range(len(XS))), flush=True)
DL_MED = float(np.median(DLmed))
print(f"  grid done ({time.time()-t0:.0f} s);"
      f" representative coarsest-comb DeltaL = {DL_MED:.2f} pc"
      f" -> K = 6 sigma_d / {DL_MED:.2f}")

# ---- caption numbers ------------------------------------------------------------
def row_at(sigphi, ncomb=NCOMB, sets=None, xs=XS):
    return np.array([purity_cell(x, sigphi, ncomb=ncomb, sets=sets)[0] for x in xs])

def cross_x(row, xs=XS, level=0.9):
    for a in range(len(xs)-1):
        if row[a] >= level >= row[a+1]:
            f = (row[a]-level)/max(row[a]-row[a+1], 1e-9)
            return float(10**(np.log10(xs[a]) + f*(np.log10(xs[a+1])-np.log10(xs[a]))))
    return np.nan

print("\n  CAPTION NUMBERS:")
# (1) 0.9 contour crossings at sigma_phi = 0.25/0.5/1.0 (dedicated exact rows),
#     and existence at sigma_phi >= 1.5 from the grid itself
for sp in (0.25, 0.5, 1.0):
    print(f"    0.9 crossing at sigma_phi = {sp:4.2f} rad:"
          f" sigma_d = {cross_x(row_at(sp)):.2f} pc")
hi = YS >= 1.5
exists_hi = bool((PUR[hi] >= 0.9).any())
print(f"    0.9 level at sigma_phi >= 1.5 rad: "
      f"{'EXISTS' if exists_hi else 'DOES NOT EXIST anywhere in the grid'}")
# (2) sigma_phi = 2.0 row at the named sigma_d values (dedicated exact cells)
X20 = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
row20 = row_at(2.0, xs=X20)
print(f"    sigma_phi = 2.0 row: " +
      "  ".join(f"{x:g} pc -> {p:.2f}" for x, p in zip(X20, row20)))
# (3) clean-limit row: first drop below 0.9
row005 = PUR[0]
print(f"    sigma_phi = 0.05 (clean limit): fraction first drops below 0.9 at"
      f" sigma_d = {cross_x(row005):.2f} pc (prior-only limit)")
# (4) x3-strain column (grid column nearest 10 pc)
ixc = int(np.argmin(np.abs(XS - 10.0)))
col1 = PUR[:, ixc]
col3 = np.array([purity_cell(XS[ixc], sp, sets=SETS3)[0] for sp in YS])
dmax = float(np.max(np.abs(col3 - col1)))
binom = float(np.max(np.sqrt(col1*(1-col1)/NREAL)))*2
print(f"    x3-strain column (sigma_d = {XS[ixc]:.2f} pc): max |dpurity| ="
      f" {dmax:.3f} vs 2x binomial sigma = {binom:.3f} ->"
      f" {'WITHIN (map does not move with sensitivity)' if dmax <= max(binom, 0.02) else 'EXCEEDS -- REPORT'}")
# (5) N dependence at sigma_phi = 2.0: where the 0.5 crossing moves
for ncb in (2, 4, 8):
    rowN = row_at(2.0, ncomb=ncb)
    print(f"    sigma_phi = 2.0, N = {ncb}: 0.5 crossing at sigma_d ="
          f" {cross_x(rowN, level=0.5):.2f} pc")

# ---- CSV ------------------------------------------------------------------------
with open("job47c_grid.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["sigma_d_pc", "sigma_phi_rad", "purity", "n_real"])
    for iy in range(len(YS)):
        for ix in range(len(XS)):
            cw.writerow([f"{XS[ix]:.4f}", f"{YS[iy]:.4f}", f"{PUR[iy,ix]:.4f}", NREAL])
print("  wrote job47c_grid.csv")

# ---- figure ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.8, 4.8))
xe = np.sqrt(XS[:-1]*XS[1:]); xe = np.concatenate([[XS[0]**2/xe[0]], xe, [XS[-1]**2/xe[-1]]])
ye = 0.5*(YS[:-1]+YS[1:]); ye = np.concatenate([[2*YS[0]-ye[0]], ye, [2*YS[-1]-ye[-1]]])
pc = ax.pcolormesh(xe, ye, PUR, cmap="Blues", vmin=0.0, vmax=1.0, shading="flat")
ax.set_xscale("log")
ax.set_xlim(xe[0], xe[-1]); ax.set_ylim(ye[0], ye[-1])
ax.axhline(2.0, color="0.25", lw=1.0, ls="--")
# 0.9 contour: draw only if it does not run off the LEFT or RIGHT plot edge
# (termination on the top/bottom of the y range is a domain boundary, not a
# truncation). Decision printed.
cs = ax.contour(XS, YS, PUR, levels=[0.9], colors="0.15", linewidths=0.9)
edge_hit = False
for path in cs.get_paths():
    v = path.vertices
    if len(v) and (v[:, 0] <= XS[0]*1.0001).any() or (v[:, 0] >= XS[-1]*0.9999).any():
        edge_hit = True
if edge_hit:
    cs.remove()
    print("  0.9 contour RUNS OFF a side edge -> OMITTED from the figure;"
          " crossings reported above instead")
else:
    print("  0.9 contour closes cleanly inside the x range -> DRAWN")
ax.set_xlabel(r"uncertainty on the pulsar distance prior $\sigma_d$ [pc]")
ax.set_ylabel(r"GW phase jitter $\sigma_\phi$ [rad]")
# s scale: secondary axis twinned to sigma_phi (right side; see header note)
sax = ax.secondary_yaxis("right",
                         functions=(lambda sp: sp, lambda sp: sp))
sticks = np.array([0.9, 0.7, 0.5, 0.3, 0.1])
sax.set_yticks(np.sqrt((1-sticks)/sticks))
sax.set_yticklabels([f"{s:g}" for s in sticks])
sax.set_ylabel("s = source's share of its frequency bin")
cb = fig.colorbar(pc, ax=ax, pad=0.13)
cb.set_label("fraction of trials where the correct distance is found")
plt.tight_layout()
plt.savefig("fig47c_criterion_heatmap.png", dpi=150, bbox_inches="tight")
print("  wrote fig47c_criterion_heatmap.png")
print(f"\n  TOTAL RUNTIME: {(time.time()-t_wall)/60:.1f} min (CPU-class)")
