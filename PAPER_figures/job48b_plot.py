"""JOB 48b -- figure in Figure 1's style + knees, overlay check, collapse panel."""
import json, csv, time
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
t0 = time.time()
OUT = "/home/mattm/projects/HSYMT/PAPER_figures/"
z = np.load(OUT + "job48_data/job48b_results.npz")
d = np.load("/home/mattm/projects/HSYMT/CW_transition/prong2_results.npz", allow_pickle=True)
grid = z["grid"]; NMAX = 1000; NPSR = 116
FIG1 = [1, 2, 3, 4, 6, 10, 16, 25, 40, 63, 100, 158, 251, 398, 631, 1000]
ARMS = [("B8", "8 frequency bins (same as prev figure)", 8.047), ("B4", "4 frequency bins", 4.0),
        ("B2", "2 frequency bins", 2.0), ("B1", "1 frequency bin", 1.0), ("SF", "all at one frequency", 1.0)]
rng = np.random.default_rng(4803)

# ---- Figure 1's exact axes (rebuild its plot to read the autoscaled y-limits) ----
BLUE, RED = "#1f5fa8", "#c8102e"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "legend.frameon": False})
f0, a0 = plt.subplots(figsize=(6.4, 5.2)); n1 = d["ps_n"]
a0.fill_between(n1, d["ps_cond_lo"], d["ps_cond_hi"]); a0.plot(n1, d["ps_cond"], "o-")
a0.fill_between(n1, d["ps_marg_lo"], d["ps_marg_hi"]); a0.plot(n1, d["ps_marg"], "s--")
a0.set_xscale("log"); a0.set_yscale("log"); YL = a0.get_ylim(); XL0 = a0.get_xlim(); plt.close(f0)
print(f"Figure 1 axes: x {XL0[0]:.3g}-{XL0[1]:.3g}, y {YL[0]:.3g}-{YL[1]:.3g}")

S = {}
for k, lab, B in ARMS:
    CM = z[f"{k}_CM"]; MM = z[f"{k}_MM"].copy(); nok = z[f"{k}_nok"]
    last_ok = []
    for s in range(MM.shape[0]):
        bad = ~np.isfinite(MM[s])
        if bad.any():
            j = np.where(~bad)[0].max(); last_ok.append(MM[s, j] / CM[s, grid[j]-1])
            MM[s, bad] = 0.0      # past Cholesky breakdown: see report (exact zero for SF; below floor for B1)
    Cg = CM[:, grid - 1]
    S[k] = dict(CM=CM, Cg=Cg, MM=MM, B=B, lab=lab, nok=nok, last_ok=np.array(last_ok))

def curve(Cg, MM):
    return np.median(MM, 0) / np.median(Cg, 0)

def cross(x, y, lv=0.5):
    for i in range(len(x) - 1):
        if y[i] >= lv > y[i+1]:
            t = (y[i] - lv) / (y[i] - y[i+1]); return float(np.exp(np.log(x[i]) + t*np.log(x[i+1]/x[i])))
    return np.nan

J48 = {"B8": 415.7, "B4": 274.0, "B2": 158.4, "B1": 102.4, "SF": 57.7}
J48I = {"B8": (411.4, 429.0), "B4": (269.3, 274.4), "B2": (156.2, 159.4), "B1": (101.2, 105.1), "SF": (55.4, 58.6)}
print("\nKNEES (median-over-skies ratio crosses 0.5; bootstrap over skies, 1000 reps)")
print(f"  {'arm':4s} {'N_knee':>8s} {'boot 16-84':>18s} {'boot 2.5-97.5':>18s} | Job 48 [boot 16-84]      (N_psr/2)(B+0.74)")
K = {}
for k, lab, B in ARMS:
    s = S[k]; kn = cross(grid, curve(s["Cg"], s["MM"]))
    bs = []
    for _ in range(1000):
        ii = rng.integers(0, s["MM"].shape[0], s["MM"].shape[0])
        bs.append(cross(grid, curve(s["Cg"][ii], s["MM"][ii])))
    bs = np.array(bs); K[k] = (kn, bs)
    pred = NPSR/2*(1.0 if k == "SF" else B + 0.74)
    print(f"  {k:4s} {kn:8.1f} [{np.percentile(bs,16):7.1f},{np.percentile(bs,84):7.1f}]"
          f" [{np.percentile(bs,2.5):7.1f},{np.percentile(bs,97.5):7.1f}] | {J48[k]:6.1f} [{J48I[k][0]:.1f},{J48I[k][1]:.1f}]"
          f"   {pred:6.1f}")

print("\nCHOLESKY BREAKDOWN (sources factored before the first non-positive pivot)")
for k, lab, B in ARMS:
    s = S[k]; nok = s["nok"]
    msg = f"  {k}: all 1000 in {np.sum(nok == NMAX)}/{len(nok)} skies"
    if (nok < NMAX).any():
        msg += (f"; breakdown at N = {nok.min()}-{nok.max()} (median {int(np.median(nok))});"
                f" ratio at the last factored grid point: median {np.median(s['last_ok']):.1e}, max {np.max(s['last_ok']):.1e}")
    print(msg)
print(f"  plot floor as a ratio at N=800: {YL[0]/np.median(S['B1']['CM'][:,799]):.1e}")

# ---- overlay check vs Figure 1's red curve ----
idx = [list(grid).index(n) for n in FIG1]
m8 = np.median(S["B8"]["MM"], 0)[idx]; mf = d["ps_marg"]; lo, hi = d["ps_marg_lo"], d["ps_marg_hi"]
rel = (m8 - mf)/mf
c8 = np.median(S["B8"]["Cg"], 0)[idx]; relc = (c8 - d["ps_cond"])/d["ps_cond"]
inside = np.sum((m8 >= lo) & (m8 <= hi))
print(f"\nOVERLAY vs Figure 1 (16 grid points): marginal max |rel diff| = {np.max(np.abs(rel)):.3f}"
      f" (at N={FIG1[int(np.argmax(np.abs(rel)))]}), median |rel diff| {np.median(np.abs(rel)):.3f};"
      f" {inside}/16 new medians inside Figure 1's 16-84 band")
print("   N   : " + " ".join(f"{n:6d}" for n in FIG1))
print("   rel : " + " ".join(f"{v:+6.3f}" for v in rel))
print(f"   conditional max |rel diff| = {np.max(np.abs(relc)):.3f}")
spr = max(np.ptp([np.median(S[k]["CM"][:, n-1]) for k, _, _ in ARMS]) / np.median(S["B8"]["CM"][:, n-1]) for n in FIG1)
print(f"   conditional: max spread of the arm medians across all 5 arms at Figure 1's N = {spr:.3f}")

# ---- the figure ----
fig, a = plt.subplots(figsize=(6.4, 5.2))
Nall = np.arange(1, NMAX + 1)
C = S["B8"]["CM"]; cmed = np.median(C, 0)
a.fill_between(Nall, np.percentile(C, 16, 0), np.percentile(C, 84, 0), color=BLUE, alpha=0.2, lw=0)
a.plot(Nall, cmed, "-", color=BLUE, lw=2, marker="o", ms=4.5, markevery=[n-1 for n in FIG1],
       label="conditional information")
cmap = plt.get_cmap("Reds"); shade = {"B8": 0.45, "B4": 0.62, "B2": 0.80, "B1": 0.97}
for k, lab, B in ARMS:
    s = S[k]; col = "k" if k == "SF" else cmap(shade[k])
    med = np.median(s["MM"], 0); plo = np.percentile(s["MM"], 16, 0); phi_ = np.percentile(s["MM"], 84, 0)
    fl = 1e-30                                    # zeros -> tiny so lines run off the bottom (clipped, not truncated)
    a.fill_between(grid, np.maximum(plo, fl), np.maximum(phi_, fl), color=col, alpha=0.15, lw=0)
    a.plot(grid, np.maximum(med, fl), "--", color=col, lw=2, marker="s", ms=4.5, markevery=idx, label=lab)
a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(XL0); a.set_ylim(YL)
a.set_xlabel("number of CW sources on the sky, N  (array fixed at N_psr = 116)")
a.set_ylabel("distance information per pulsar, $I_{LL} = 1/\\sigma_L^2$  [kpc$^{-2}$]")
sec = a.secondary_yaxis("right", functions=(lambda I: 1e3/np.sqrt(I), lambda s: (1e3/s)**2))
sec.set_ylabel("distance precision $\\sigma_L$  [pc]")
a.grid(alpha=0.25, which="both"); a.legend(fontsize=8.5, loc="upper left")
fig.savefig(OUT + "job48b_turnover_per_bin_fig1style.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---- Step 2: collapse on N/(B+0.74) ----
fig, a = plt.subplots(figsize=(6.4, 5.2))
dev = {}
for k, lab, B in ARMS:
    s = S[k]; Be = 1.0 if k == "SF" else B + 0.74; col = "k" if k == "SF" else cmap(shade[k])
    x = grid/Be; r = s["MM"]/s["Cg"]
    med = curve(s["Cg"], s["MM"])
    a.fill_between(x, np.percentile(r, 16, 0), np.percentile(r, 84, 0), color=col, alpha=0.15, lw=0)
    a.plot(x, med, "-", color=col, lw=1.8, label=lab)
    kk = K[k][0]/Be; m = x < kk
    dd = (med - (1 - x/NPSR))[m]; dev[k] = dd[np.argmax(np.abs(dd))]
xx = np.geomspace(0.1, 1000, 400)
a.plot(xx, np.clip(1 - xx/NPSR, 0, None), "--", color="0.5", lw=1, label="1 − x/N$_{psr}$")
a.axhline(0.5, color="k", lw=0.8, alpha=0.5)
a.set_xscale("log"); a.set_ylim(0, 1.05); a.set_xlim(0.8, 700)
a.set_xlabel("effective sources per frequency bin, N/(B + 0.74)")
a.set_ylabel("surviving fraction of distance information,  $I_{\\rm marg}/I_{\\rm cond}$")
a.grid(alpha=0.25, which="both"); a.legend(fontsize=8.5, loc="lower left")
fig.savefig(OUT + "job48b_collapse_Beff.png", dpi=200, bbox_inches="tight"); plt.close(fig)
print("\nCOLLAPSE on N/(B+0.74): largest deviation from 1 - x/N_psr below the knee: " +
      ", ".join(f"{k} {v:+.3f}" for k, v in dev.items()))

# ---- CSV ----
with open(OUT + "job48b_turnover_per_bin.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["arm", "B", "N", "n_skies", "I_cond_median", "I_cond_p16", "I_cond_p84",
                "I_marg_median", "I_marg_p16", "I_marg_p84", "ratio_median_of_medians"])
    for k, lab, B in ARMS:
        s = S[k]
        for j, N in enumerate(grid):
            c = s["Cg"][:, j]; m = s["MM"][:, j]
            w.writerow([k, B, int(N), len(c), f"{np.median(c):.5e}", f"{np.percentile(c,16):.5e}",
                        f"{np.percentile(c,84):.5e}", f"{np.median(m):.5e}", f"{np.percentile(m,16):.5e}",
                        f"{np.percentile(m,84):.5e}", f"{np.median(m)/np.median(c):.6f}"])
print(f"\nplotting + analysis wall: {time.time()-t0:.1f} s")
