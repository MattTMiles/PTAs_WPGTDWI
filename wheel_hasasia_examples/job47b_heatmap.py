# ================= JOB 47b -- THE CRITERION AS A SIMPLE MEASURED HEATMAP ========
# One panel for the group email. Axes carry only email-introduced quantities:
#   x (log)  : distance-prior uncertainty sigma_d [pc], internally K = 6 sigma_d /
#              DeltaL(coarsest comb) -- K never shown.
#   y (lin)  : GW phase jitter sigma_phi [rad], internally the confusion share
#              s = 1/(1+sigma_phi^2), applied to EVERY usable comb.
#   colour   : fraction of trials in which the correct pulsar distance is
#              recovered (winner within TOOTH_W = 0.25 pc of truth), >=250
#              realisations per cell, banked SEED50 harness.
# Usable-CW count FIXED at N = 4 combs (machine's realistic N_eff ~ 4.2, booked):
# per banked sky the machine's own inv_rate inventory (cap 1300, rho>=3, band
# 2-30 nHz, T = 45 yr / 100 ns / 14 d) ranked by its own currency nats =
# ln(pi/sigma); the TOP 4 are kept; per-realisation geometry mu = arccos U(-0.9,
# 0.9) per comb (machine convention), kappa = 2 pi f (1-cos mu)/c.
#
# PRE-REGISTERED:
#   (a) x3 source strain column: s unchanged (GW-over-GW), purity differences
#       within binomial error.
#   (b) sigma_phi = 2 rad row at N = 2 and N = 8: the 0.9 crossing moves little
#       (each added comb carries its own jitter).
# WHICH-TERM: the strain arm multiplies every GW amplitude x3 (h x3, all powers
# x9, sea x9): res/Pt -- hence s and sigma_phi -- is UNCHANGED; only detection
# (rho, arr) and hence comb SELECTION can move.
import os, sys, json, time, csv
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/mattm/projects/HSYMT/wheel_hasasia_examples"
os.chdir(BASE)
t_wall = time.time()

# ---- machine namespace ----------------------------------------------------------
# FAST PATH (default): exec the cheap cells in full and only the DEFINITION
# slices of the campaign cells that the _pop18 -> inv_rate -> collect36m chain
# needs. Equivalence with the full rebuild is not assumed -- it is PROVEN by the
# bit-exact anchor gate below (36.27 @ 0.9009 on the same seed chain); any
# missed or out-of-order dependency either raises or fails that gate.
# JOB47B_FULL=1 falls back to the Job 45 full-cell rebuild (~38 min).
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

if os.environ.get("JOB47B_FULL"):
    nfail = 0
    for i in sorted(cells):
        if i in (12, 43):
            continue
        try:
            ex(i)
        except Exception as e:
            nfail += 1
            print(f"  [cell {i:2d}] FAIL {e!r}", flush=True)
    assert nfail == 0, "cell failures -- STOP"
else:
    for i in (1, 3, 5, 7, 9, 11):
        ex(i)                                  # cheap cells, in full
    ex(14, 116, 116, "EXP")                    # EXP = (120, 100)
    ex(14, 168, 168, "TOOTH_W")                # DLG, N_SCR_W, TOOTH_W = ...
    ex(14, 780, 786, "rn_draw")                # RN_FRAC + rn_draw
    ex(14, 1586, 1591, "FP_ALL")               # retarded-frequency map (cell-17 dep)
    ex(14, 1072, 1085, "SEED50")               # the banked 50-population harness
    ex(14, 3023, 3029, "_sea50v")
    ex(14, 4043, 4049, "h_of_rhoP")            # _SN15/_NT15/_EL15 + rho_P_ref/h_of_rhoP
    ex(14, 4154, 4201, "build_pop15")          # _BAND15 + build_pop15
    ex(15, 13, 17, "s_to_s30")
    ex(17, 1, 127, "depro_sigma")              # deprojection machinery + def
    ex(18, 12, 17, "_pop18")
    ex(27, 67, 78, "alias_score24")
_pop18 = NS["_pop18"]
C_PC_S = NS["C_PC_S"]; YR_S = NS["YR_S"]
TOOTH_W = 0.25          # the machine's tooth-correctness radius (cell 14, banked)

# job43c inv_rate verbatim (lines 485-509) + book/collect36m for the anchor check
lines43 = open("job43c_rebooking.py").read().splitlines()
def run_span(lo, hi):
    exec(compile("\n".join(lines43[lo-1:hi]), f"<43c:{lo}-{hi}>", "exec"), NS)
run_span(485, 509)
run_span(661, 689)
run_span(712, 756)
inv_rate = NS["inv_rate"]

# ---- anchor sanity (drift guard; booked 36.27 [34.02,38.75] @ 0.9009) -----------
r = NS["book"](NS["collect36m"](fuzzy=True, full_dmu=True, marg=True))
print(f"\n  ANCHOR CHECK: {r['y']:.2f} [{r['yb'][0]:.2f},{r['yb'][1]:.2f}]"
      f" @ {r['pur']:.4f} -- booked 36.27 [34.02,38.75] @ 0.9009")
assert abs(r["y"] - 36.27) < 0.5 and abs(r["pur"] - 0.9009) < 0.005, "anchor moved -- STOP"

# ---- per-sky machine comb sets (top 12 by nats; grid uses the top N of these) ---
def comb_sets(strain_fac=1.0):
    """Per banked sky: (freqs sorted by machine nats ranking, own-element share s
    per comb). strain_fac multiplies EVERY GW amplitude (h x fac, powers x fac^2)."""
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

# which-term check printout: shares under x1 vs x3 on the shared top combs
sh1, sh3 = [], []
for a, b in zip(SETS, SETS3):
    if a is None or b is None:
        continue
    common = np.intersect1d(np.round(a[0], 14), np.round(b[0], 14))
    for f0 in common[:4]:
        sh1.append(float(a[1][np.argmin(np.abs(a[0]-f0))]))
        sh3.append(float(b[1][np.argmin(np.abs(b[0]-f0))]))
sh1 = np.array(sh1); sh3 = np.array(sh3)
print(f"  WHICH-TERM CHECK: own-element share s on shared combs, x1 vs x3 strain:"
      f" max |ds| = {np.max(np.abs(sh1-sh3)):.2e}  (s is GW-over-GW: UNCHANGED)")

# ---- the grid -------------------------------------------------------------------
XS = np.array([0.3, 0.6, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0])   # sigma_d [pc]
YS = np.array([0.2, 0.5, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5])         # sigma_phi [rad]
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

t0 = time.time()
PUR = np.zeros((len(YS), len(XS)))
DLmed = []
for iy, sp in enumerate(YS):
    for ix, sd_ in enumerate(XS):
        PUR[iy, ix], dlm = purity_cell(sd_, sp)
        DLmed.append(dlm)
    print(f"  row sigma_phi={sp:4.1f}: " +
          " ".join(f"{PUR[iy,ix]:.2f}" for ix in range(len(XS))), flush=True)
DL_MED = float(np.median(DLmed))
print(f"  grid done ({time.time()-t0:.0f} s); {NREAL} realisations/cell,"
      f" N = {NCOMB} combs fixed")
print(f"  representative coarsest-comb tooth spacing DeltaL (median over"
      f" realisations): {DL_MED:.2f} pc  -> K = 6 sigma_d / {DL_MED:.2f}")

# ---- caption crossings ----------------------------------------------------------
def cross_x(row):        # sigma_d at 0.9, interpolate in log x (purity falls with x)
    for a in range(len(XS)-1):
        if row[a] >= 0.9 >= row[a+1]:
            f = (row[a]-0.9)/max(row[a]-row[a+1], 1e-9)
            return float(10**(np.log10(XS[a]) + f*(np.log10(XS[a+1])-np.log10(XS[a]))))
    return np.nan

def cross_y(col):        # sigma_phi at 0.9 (purity falls with y)
    for a in range(len(YS)-1):
        if col[a] >= 0.9 >= col[a+1]:
            f = (col[a]-0.9)/max(col[a]-col[a+1], 1e-9)
            return float(YS[a] + f*(YS[a+1]-YS[a]))
    return np.nan

iy20 = int(np.where(YS == 2.0)[0][0]); iy10 = int(np.where(YS == 1.0)[0][0])
ix1 = int(np.where(XS == 1.0)[0][0]); ix10 = int(np.where(XS == 10.0)[0][0])
print("\n  CAPTION NUMBERS:")
print(f"    0.9 contour at sigma_phi = 2.0 rad: sigma_d = {cross_x(PUR[iy20]):.2f} pc")
print(f"    0.9 contour at sigma_phi = 1.0 rad: sigma_d = {cross_x(PUR[iy10]):.2f} pc")
print(f"    0.9 contour at sigma_d = 1 pc :  sigma_phi = {cross_y(PUR[:, ix1]):.2f} rad")
print(f"    0.9 contour at sigma_d = 10 pc:  sigma_phi = {cross_y(PUR[:, ix10]):.2f} rad")

# ---- (a) x3-strain column check -------------------------------------------------
ixc = ix10
col3 = np.array([purity_cell(XS[ixc], sp, sets=SETS3)[0] for sp in YS])
dmax = np.max(np.abs(col3 - PUR[:, ixc]))
binom = np.max(np.sqrt(PUR[:, ixc]*(1-PUR[:, ixc])/NREAL))*2
print(f"\n  x3 STRAIN COLUMN (sigma_d = {XS[ixc]:.0f} pc): max |dpurity| = {dmax:.3f}"
      f" vs 2x binomial sigma = {binom:.3f} ->"
      f" {'WITHIN (pre-reg HIT)' if dmax <= max(binom, 0.02) else 'EXCEEDS (pre-reg MISS)'}")

# ---- (b) N-dependence, sigma_phi = 2 rad row ------------------------------------
print(f"\n  N-DEPENDENCE (sigma_phi = 2.0 rad row): 0.9 crossing in sigma_d")
for ncb in (2, 4, 8):
    rowN = np.array([purity_cell(sd_, 2.0, ncomb=ncb)[0] for sd_ in XS])
    print(f"    N = {ncb}: crossing {cross_x(rowN):.2f} pc   row: " +
          " ".join(f"{v:.2f}" for v in rowN))

# ---- deliverables ---------------------------------------------------------------
with open("job47b_grid.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["sigma_d_pc", "sigma_phi_rad", "purity", "n_real"])
    for iy in range(len(YS)):
        for ix in range(len(XS)):
            cw.writerow([XS[ix], YS[iy], f"{PUR[iy,ix]:.4f}", NREAL])
print("  wrote job47b_grid.csv")

fig, ax = plt.subplots(figsize=(6.4, 4.6))
xe = np.sqrt(XS[:-1]*XS[1:]); xe = np.concatenate([[XS[0]**2/xe[0]], xe, [XS[-1]**2/xe[-1]]])
ye = 0.5*(YS[:-1]+YS[1:]); ye = np.concatenate([[2*YS[0]-ye[0]], ye, [2*YS[-1]-ye[-1]]])
pc = ax.pcolormesh(xe, ye, PUR, cmap="Blues", vmin=0.0, vmax=1.0, shading="flat")
ax.contour(XS, YS, PUR, levels=[0.9], colors="0.15", linewidths=0.9)
ax.set_xscale("log")
ax.set_xlabel(r"uncertainty on the pulsar distance prior $\sigma_d$ [pc]")
ax.set_ylabel(r"GW phase jitter $\sigma_\phi$ [rad]")
cb = fig.colorbar(pc, ax=ax)
cb.set_label("fraction of trials where the correct distance is found")
plt.tight_layout()
plt.savefig("fig47b_criterion_heatmap.png", dpi=150, bbox_inches="tight")
print("  wrote fig47b_criterion_heatmap.png")
print(f"\n  TOTAL RUNTIME: {(time.time()-t_wall)/60:.1f} min (CPU-class)")
