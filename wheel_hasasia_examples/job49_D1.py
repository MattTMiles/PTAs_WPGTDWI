# ================= JOB 49 D1 -- the 47c/47d heatmap, incumbent vs physical =========
# (1) Reproduce job47c_grid.csv (20x20, 250 real./cell) with the booked purity_cell
#     (job47c_heatmap.py:108-128, exec'd from the file) and job47d_grid.csv (47d band
#     refinement at 1000 real.).  Exact equality required.
# (2) purity_cell_sw: same seed formula, same skies/combs/geometry/L_true, and the
#     incumbent eps draw is consumed exactly as booked (bit-identical incumbent).
#     Physical datum per comb, from a SEPARATE stream:
#        z = e^{i kappa L_true} (1 + C),  E|C|^2 = P = (1-s)/s = sigma_phi^2
#     where s is the cell's share (47c: s = 1/(1+sigma_phi^2)).  C is the comb's
#     ACTUAL contributor mix (its 64 loudest in-window sources as explicit phasors,
#     random phases, + the rest and the sea as an exact-power complex Gaussian),
#     rescaled to total power P.  47c has no measurement noise / fuzzy term; none added.
#     D-inc: booked wrapped-quadratic decoder with sigma_phi.  D-mf: 2 Re(z e^{-i k L})/P.
# (3) Each arm gets the 47d treatment: band cells within one cell of ITS 0.9 contour
#     re-measured at 1000 realisations; contour from the 3x3-mean field; cells raw.
import os, sys, json, time, csv
import numpy as np
from scipy.ndimage import uniform_filter, binary_dilation
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/mattm/projects/HSYMT/wheel_hasasia_examples")
import job49_common as J

os.chdir(J.BASE); os.makedirs("job49_data", exist_ok=True)
t_all = time.time()
NS = J.load(verbose=False)
_pop18 = NS["_pop18"]; inv_rate = NS["inv_rate"]; C_PC_S = NS["C_PC_S"]
TOOTH_W = 0.25
print(f"[load] {time.time()-t_all:.0f} s", flush=True)

# ---- comb sets exactly as 47d (fk sorted by nats) + contributor shapes ----------
SETS, SHAPE = [], []
for j in range(50):
    sd, f2, h2, P0, sea = _pop18(j, 0.35)
    rr = inv_rate(f2, h2, P0, sea, 45.0)
    if rr is None:
        SETS.append(None); SHAPE.append(None); continue
    SETS.append(rr[2])
    cb = J._sky_contrib(f2, P0, sea, rr[0], 45.0, 1.0)
    # normalise each comb's mix to unit total power
    SHAPE.append([(top/cf, rem/cf) for top, rem, cf in cb])
print(f"  comb sets: {sum(s is not None for s in SETS)}/50 skies", flush=True)

# ---- booked purity_cell (47c:108-128), exec'd from the file ---------------------
L47 = open("job47c_heatmap.py").read().splitlines()
NS47 = dict(np=np, SETS=[None if s is None else (s, None) for s in SETS],
            NREAL=250, NCOMB=4, C_PC_S=C_PC_S, TOOTH_W=TOOTH_W)
exec(compile("\n".join(L47[108-1:128]), "<47c:108-128>", "exec"), NS47)
purity_cell_booked = NS47["purity_cell"]
NCOMB = 4


def purity_cell_sw(sigd, sigphi, nreal, seed0=470_000_000):
    """47c/47d purity_cell + physical arms.  Returns fractions (inc, phys_inc, phys_mf)."""
    ok = np.zeros(3, int)
    P = sigphi**2
    for t in range(nreal):
        jj = t % 50
        fk = SETS[jj]; shp = SHAPE[jj]
        if fk is None:
            fk = SETS[(jj+1) % 50]; shp = SHAPE[(jj+1) % 50]
        fk = fk[:NCOMB]
        rg = np.random.default_rng(seed0 + t*7919 + int(1e6*sigd) % 99991
                                   + int(1e6*sigphi) % 9973)
        mu = np.arccos(rg.uniform(-0.9, 0.9, len(fk)))
        kap = 2*np.pi*fk*(1-np.cos(mu))/C_PC_S
        Lt = float(np.clip(rg.normal(0, sigd), -2.9*sigd, 2.9*sigd))
        Lg = np.linspace(-3*sigd, 3*sigd, 601)
        eps = kap*Lt + rg.normal(0, sigphi, len(kap))               # booked
        per = (kap[:, None]*Lg[None, :] - eps[:, None] + np.pi) % (2*np.pi) - np.pi
        S = -(per**2).sum(0)/(2*sigphi**2) - Lg**2/(2*sigd**2)
        ok[0] += int(abs(float(Lg[int(np.argmax(S))]) - Lt) <= TOOTH_W)
        # physical datum
        rph = np.random.default_rng([49, 47, t, int(1e6*sigd), int(1e6*sigphi)])
        u = np.ones(len(fk), complex)
        for c in range(len(fk)):
            top, rem = shp[c]
            th = rph.uniform(0, 2*np.pi, len(top))
            u[c] += np.sum(np.sqrt(P*top)*np.exp(1j*th))
            g = rph.normal(0, 1, 2)*np.sqrt(P*rem/2)
            u[c] += g[0] + 1j*g[1]
        epp = kap*Lt + np.angle(u)
        per = (kap[:, None]*Lg[None, :] - epp[:, None] + np.pi) % (2*np.pi) - np.pi
        S = -(per**2).sum(0)/(2*sigphi**2) - Lg**2/(2*sigd**2)
        ok[1] += int(abs(float(Lg[int(np.argmax(S))]) - Lt) <= TOOTH_W)
        z = np.exp(1j*kap*Lt)*u
        S = (2.0/P)*np.real(z[:, None]*np.exp(-1j*kap[:, None]*Lg[None, :])).sum(0) \
            - Lg**2/(2*sigd**2)
        ok[2] += int(abs(float(Lg[int(np.argmax(S))]) - Lt) <= TOOTH_W)
    return ok/nreal


XS = np.geomspace(0.1, 300.0, 20)
YS = np.linspace(0.05, 3.0, 20)
ARM = ("incumbent", "physical D-inc", "physical D-mf")

# ---- (1)+(2) grid ---------------------------------------------------------------
t0 = time.time()
G = np.zeros((3, 20, 20)); Gb = np.zeros((20, 20))
for iy, sp in enumerate(YS):
    for ix, x in enumerate(XS):
        G[:, iy, ix] = purity_cell_sw(x, sp, 250)
        Gb[iy, ix] = purity_cell_booked(x, sp)[0]
print(f"  grid done {time.time()-t0:.0f} s", flush=True)
ref = {}
for r in csv.DictReader(open("job47c_grid.csv")):
    ref[(round(float(r["sigma_d_pc"]), 4), round(float(r["sigma_phi_rad"]), 4))] = float(r["purity"])
d47c = max(abs(round(G[0, iy, ix], 4) - ref[(round(XS[ix], 4), round(YS[iy], 4))])
           for iy in range(20) for ix in range(20))
dbk = float(np.max(np.abs(G[0]-Gb)))
print(f"  REPRO 47c: max |switch - job47c_grid.csv| = {d47c:.2e}; max |switch - booked fn| = {dbk:.2e}")

# ---- (3) 47d band refinement, per arm --------------------------------------------
NHI = 1000
NRE = np.full((3, 20, 20), 250)
G_ref = G.copy()
def band_of(P):
    above = P >= 0.9
    st = np.zeros_like(above)
    st[:, :-1] |= above[:, :-1] != above[:, 1:]; st[:, 1:] |= above[:, :-1] != above[:, 1:]
    st[:-1, :] |= above[:-1, :] != above[1:, :]; st[1:, :] |= above[:-1, :] != above[1:, :]
    return binary_dilation(st, iterations=1)
bands = [band_of(G[a]) for a in range(3)]
union = bands[0] | bands[1] | bands[2]
t0 = time.time()
# 47d reads the grid back from job47c_grid.csv (4-decimal strings), so its seeds use
# the ROUNDED sigma_d / sigma_phi; reproduce that exactly (same cells, same realisations)
XSr = np.round(XS, 4); YSr = np.round(YS, 4)
for iy, ix in zip(*np.where(union)):
    v = purity_cell_sw(XSr[ix], YSr[iy], NHI)
    for a in range(3):
        if bands[a][iy, ix]:
            G_ref[a, iy, ix] = v[a]; NRE[a, iy, ix] = NHI
print(f"  refinement: {int(union.sum())} cells at {NHI} ({time.time()-t0:.0f} s)", flush=True)
ref47d = {}
for r in csv.DictReader(open("job47d_grid.csv")):
    ref47d[(round(float(r["sigma_d_pc"]), 4), round(float(r["sigma_phi_rad"]), 4))] = (float(r["purity"]), int(r["n_real"]))
d47d = max(abs(round(G_ref[0, iy, ix], 4) - ref47d[(round(XS[ix], 4), round(YS[iy], 4))][0])
           for iy in range(20) for ix in range(20))
nmis = sum(NRE[0, iy, ix] != ref47d[(round(XS[ix], 4), round(YS[iy], 4))][1]
           for iy in range(20) for ix in range(20))
print(f"  REPRO 47d: max |incumbent refined - job47d_grid.csv| = {d47d:.2e}; n_real mismatches {nmis}")
SM = np.array([uniform_filter(G_ref[a], size=3, mode="nearest") for a in range(3)])


def cross_x(row, xs, level=0.9):
    for a in range(len(xs)-1):
        if row[a] >= level >= row[a+1]:
            f = (row[a]-level)/max(row[a]-row[a+1], 1e-9)
            return float(10**(np.log10(xs[a]) + f*(np.log10(xs[a+1])-np.log10(xs[a]))))
    return np.nan


# ---- dedicated rows: crossings at beta = 0.8/0.5/0.3/0.2, and sqrt2-rescaled ----
out = dict(repro_47c=d47c, repro_47d=d47d, repro_47d_nreal_mismatch=int(nmis), rows={})
XR = np.geomspace(0.1, 300.0, 40)
print("\n  0.9-CROSSING sigma_d [pc], dedicated rows, 1000 real./cell (40 log cells):")
print(f"  {'beta':>5} {'sigma_phi':>9} | {'incumbent':>10} {'phys D-inc':>11} {'phys D-mf':>10}"
      f" | incumbent at sigma_phi/sqrt2")
for b in (0.8, 0.5, 0.3, 0.2):
    sp = float(np.sqrt((1-b)/b))
    R = np.array([purity_cell_sw(x, sp, NHI) for x in XR]).T
    Rh = np.array([purity_cell_sw(x, sp/np.sqrt(2), NHI)[0] for x in XR])
    cx = [cross_x(R[a], XR) for a in range(3)]; ch = cross_x(Rh, XR)
    out["rows"][str(b)] = dict(sigma_phi=sp, xs=XR.tolist(), frac=R.tolist(), cross=cx,
                               inc_sqrt2=Rh.tolist(), cross_inc_sqrt2=ch)
    print(f"  {b:5.2f} {sp:9.3f} | {cx[0]:10.3f} {cx[1]:11.3f} {cx[2]:10.3f} | {ch:.3f}", flush=True)

print("\n  RECOVERED FRACTION along beta = 0.2 (sigma_phi = 2.0), 4000 real./cell:")
X02 = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0])
R02 = np.array([purity_cell_sw(x, 2.0, 4000) for x in X02])
out["beta02"] = dict(xs=X02.tolist(), frac=R02.tolist(), n=4000)
for x, v in zip(X02, R02):
    se = np.sqrt(v*(1-v)/4000)
    print(f"    sigma_d {x:6.1f} pc: " + "  ".join(f"{ARM[a]} {v[a]:.3f}+-{se[a]:.3f}" for a in range(3)))

# ---- CSV + JSON -----------------------------------------------------------------
with open("job49_D1_grid.csv", "w", newline="") as fh:
    cw = csv.writer(fh)
    cw.writerow(["sigma_d_pc", "sigma_phi_rad", "beta", "arm", "frac_raw", "n_real", "frac_3x3_contour_field"])
    for a in range(3):
        for iy in range(20):
            for ix in range(20):
                cw.writerow([f"{XS[ix]:.4f}", f"{YS[iy]:.4f}", f"{1/(1+YS[iy]**2):.5f}", ARM[a],
                             f"{G_ref[a,iy,ix]:.4f}", NRE[a, iy, ix], f"{SM[a,iy,ix]:.4f}"])
json.dump(out, open("job49_data/D1_rows.json", "w"), indent=1)
np.savez("job49_data/D1_grid.npz", XS=XS, YS=YS, G250=G, G=G_ref, NRE=NRE, SM=SM)
print("  wrote job49_D1_grid.csv, job49_data/D1_rows.json, job49_data/D1_grid.npz")

# ---- figure: beta on the vertical axis -------------------------------------------
BE = 1/(1+YS**2)
ye_s = 0.5*(YS[:-1]+YS[1:]); ye_s = np.concatenate([[max(2*YS[0]-ye_s[0], 0.0)], ye_s, [2*YS[-1]-ye_s[-1]]])
ye_b = 1/(1+ye_s**2)                              # beta edges (decreasing)
xe = np.sqrt(XS[:-1]*XS[1:]); xe = np.concatenate([[XS[0]**2/xe[0]], xe, [XS[-1]**2/xe[-1]]])
fig, ax = plt.subplots(1, 4, figsize=(18, 4.6), gridspec_kw=dict(width_ratios=[1, 1, 1, 1.05]))
for a in range(3):
    pc = ax[a].pcolormesh(xe, ye_b, G_ref[a], cmap="Blues", vmin=0, vmax=1, shading="flat")
    ax[a].contour(XS, BE, SM[a], levels=[0.9], colors="0.1", linewidths=1.0)
    ax[a].set_xscale("log"); ax[a].set_title(ARM[a])
    ax[a].set_xlabel(r"$\sigma_d$ [pc]")
    sx = ax[a].secondary_yaxis("right", functions=(lambda b: np.sqrt(np.clip((1-b)/np.clip(b, 1e-6, 1), 0, None)),
                                                    lambda s: 1/(1+s**2)))
    if a == 2: sx.set_ylabel(r"old axis $\sigma_\phi=\sqrt{(1-\beta)/\beta}$ [rad]")
ax[0].set_ylabel(r"$\beta$ = source's share of its bin")
fig.colorbar(pc, ax=ax[2], pad=0.18).set_label("fraction correct")
cols = ("0.15", "C0", "C3")
for a in range(3):
    ax[3].contour(XS, BE, SM[a], levels=[0.9], colors=cols[a], linewidths=1.4)
    ax[3].plot([], [], color=cols[a], label=ARM[a])
ax[3].set_xscale("log"); ax[3].set_xlim(xe[0], xe[-1]); ax[3].set_ylim(ye_b[-1], ye_b[0])
ax[3].set_title("0.9 contours (3x3-mean field)"); ax[3].legend(fontsize=8, loc="lower left")
ax[3].set_xlabel(r"$\sigma_d$ [pc]")
for a_ in ax: a_.grid(alpha=.2, lw=.5)
plt.tight_layout(); plt.savefig("job49_D1_heatmaps.png", dpi=140); plt.close()
print("  wrote job49_D1_heatmaps.png")
print(f"TOTAL {time.time()-t_all:.0f} s")
