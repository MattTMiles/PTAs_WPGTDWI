# ================= JOB 54 Part 2 (detection reach): strain at array SNR 10 against number of phased pulsars ===
# Array: the localisation array (job54_loc_core.build_all: 116 pulsars, census distances, nested orderings).
# Sky- and orientation-averaged matched-filter SNR^2 for unit strain, per pulsar:
#   Earth term   e(f)   = 0.2 (n_toa/2) [1/(2 pi f sigma)]^2 / R(f)
#   pulsar term  p_a(f) = < w(c) (n_toa/2) [1/(2 pi f sigma)]^2 (f_P/f)^(-2/3) / R(f_P) >_c ,  c = cos(mu) uniform in (-1, 1),
#                         w(c) = 0.6 (1+c)^2/4  (mean 0.2: the response is strongest where the lag L(1-c)/c is shortest),
#                         f_P = f (1 + tau/t_c)^(-3/8),  tau = L_a (1-c)/c
#   0.2 per pulsar is the machine's convention (job43c_rebooking.py:492  arr=rp*np.sqrt(0.2*120)).
# Filters:  (i)  rho^2 = h^2 N_psr e                 Earth terms only
#           (ii) rho^2 = h^2 [N_psr e + sum_phased p_a]      + phased pulsar terms, coherent; the rest ignored
#           (iii) rho^2 = (ii) + sum_unphased x_a^2/(4 + 4 x_a), x_a = h^2 p_a   unphased pulsar terms with their phase
#                 marginalised = an energy statistic (2 degrees of freedom per pulsar, optimal weights); this is its
#                 deflection squared.  Cross terms between Earth and pulsar terms average to zero over phase and are dropped.
# Noise:   (A) timing noise only: R = 1.
#          (B) timing noise + every other gravitational-wave contribution in the bin, per pulsar, as Job 51's arr_bg
#              (job51_common.py:112  self.oth_bg.append(np.maximum(loc["gEs"](fc)-Pc*sc+gP+2*bs, 0)*rt)
#               job51_common.py:125  arr = lambda oth: np.sqrt(0.2*120*np.mean(2.0*Psig[None, :]/(Nt[None, :]+oth), axis=0)) ):
#              R(f) = 1 + [gE(f) + <gP(f)> + 2 sea(f)]/N_t(f) on the first natural sky (sky 100), 45-yr bins.
# Timing set-ups: "loc" = the localisation array (15 yr, 7 d, 100 ns); "machine" = Job 51's (45 yr, 14 d, 100 ns).  (B) needs
# the machine's 45-yr bins, so it is given for "machine" timing only.
import os, sys, csv, json, pickle
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
os.chdir("/home/mattm/projects/HSYMT/wheel_hasasia_examples")
YR = 365.25*86400; PC_S = 1.0/9.715611890180196e-09; TSUN = 4.925490947e-6; SIG = 1e-7
NPH = [0, 1, 2, 3, 4, 6, 10, 20, 40, 80, 116]
# ---- array (same draws as job54_loc_core.build_all) ----------------------------------------------------
rows = list(csv.DictReader(open("job45_t2_parallax_census.csv")))
D = np.array([float(r["D_kpc"]) for r in rows])*1000.0
L = D[np.random.default_rng(54).integers(0, len(D), 116)]; NP = len(L)
ORD = {"nearest_first": np.argsort(L, kind="stable"), "random": np.random.default_rng(5401).permutation(NP)}
if os.path.exists("job54_data/array_check.json"):
    chk = json.load(open("job54_data/array_check.json")); assert np.allclose(chk["L_pc"], L) and chk["order_near"] == ORD["nearest_first"].tolist() and chk["order_rand"] == ORD["random"].tolist()
    print("array equal to job54_loc_core.build_all: distances and both orderings identical")
CAP = pickle.load(open("job54_data/sens_capture_sky100.pkl", "rb"))
NT = {"loc": 15*365.25/7.0, "machine": 45*365.25/14.0}
cn, cw = np.polynomial.legendre.leggauss(96); wc = 0.6*(1+cn)**2/4*cw/2          # nodes in c, weights incl. response
Rbg = 1.0+(CAP["gE"]+CAP["gP_mean"]+2*CAP["bsea"])*CAP["rt"]/(CAP["Ncal"]*CAP["fg"]**3)
from scipy.ndimage import median_filter
# for the CURVES the background is smoothed (running median over 41 grid points = 3.6 nHz), so that a loud source does not raise
# the curve at its own frequency; the listed binaries use the unsmoothed field with their own power removed (Job 51's rule)
Rbg_s = 1.0+10**median_filter(np.log10((CAP["gE"]+CAP["gP_mean"]+2*CAP["bsea"])*CAP["rt"]), size=41, mode="nearest")/(CAP["Ncal"]*CAP["fg"]**3)
SMOOTH = [True]
def Rf(f, ver): return np.ones_like(f) if ver == "A" else np.interp(f, CAP["fg"], Rbg_s if SMOOTH[0] else Rbg)
def tchirp(f, mc): return (5.0/256.0)*(TSUN*mc)**(-5.0/3.0)*(np.pi*f)**(-8.0/3.0)
def parts(f, tc, timing, ver, own=None):
    """e(f) per pulsar [scalar per f] and p_a(f) [n_f, N_psr], unit strain.  own: the source's own Earth-term power, removed
    from the noise at its Earth frequency (as Job 51 does); used for the listed binaries, not for the curves."""
    f = np.atleast_1d(f); base = (NT[timing]/2)/(2*np.pi*f*SIG)**2
    Re = Rf(f, ver) if (own is None or ver == "A") else 1.0+np.maximum((Rf(f, ver)-1.0)*CAP["Ncal"]*f**3-own*CAP["rt"], 0)/(CAP["Ncal"]*f**3)
    e = 0.2*base/Re
    tau = L[None, :, None]*PC_S*(1-cn[None, None, :])
    fr = np.ones((len(f), NP, len(cn))) if tc is None else (1+tau/np.atleast_1d(tc)[:, None, None])**(-3.0/8.0)
    p = base[:, None]*np.sum(wc[None, None, :]*fr**(-2.0/3.0)/Rf(f[:, None, None]*fr, ver), axis=2)
    return e, p
def h_at(e, p, ph, filt, snr=10.0):
    """strain giving SNR = snr; e [n_f], p [n_f, N_psr], ph = boolean phased mask"""
    if filt == "i": return snr/np.sqrt(NP*e)
    coh = NP*e+p[:, ph].sum(1)
    if filt == "ii" or ph.all(): return snr/np.sqrt(coh)
    lo = np.zeros(len(e)); hi = snr/np.sqrt(coh)
    for _ in range(60):
        h = 0.5*(lo+hi); x = h[:, None]**2*p[:, ~ph]; r2 = h**2*coh+np.sum(x**2/(4+4*x), 1)
        lo = np.where(r2 < snr**2, h, lo); hi = np.where(r2 < snr**2, hi, h)
    return 0.5*(lo+hi)

FG = np.geomspace(2e-9, 1e-7, 120); rowsC = []; CUR = {}
for chirp, mc in (("Mc_1e9_Msun", 1e9), ("non_evolving", None)):
    tc = None if mc is None else tchirp(FG, mc)
    for timing, ver in (("loc", "A"), ("machine", "A"), ("machine", "B")):
        e, p = parts(FG, tc, timing, ver)
        for on, o in ORD.items():
            for n in NPH:
                ph = np.zeros(NP, bool); ph[o[:n]] = True
                for filt in ("i", "ii", "iii"):
                    h = h_at(e, p, ph, filt); CUR[(chirp, timing, ver, on, n, filt)] = h
                    rowsC += [[f"{f*1e9:.5f}", n, on, filt, f"{ver}_{timing}_timing", chirp, f"{x:.5e}"] for f, x in zip(FG, h)]
with open("talk_data/sensitivity_vs_nph.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["frequency_nHz", "n_phased", "ordering", "filter", "noise_version", "source", "strain_at_snr10"]); w.writerows(rowsC)
print(f"wrote talk_data/sensitivity_vs_nph.csv: {len(rowsC)} rows")
print("\n== ratio of the strain at SNR 10 to Earth-terms-only (filter i); median over 2-100 nHz [value at 5, 10, 30 nHz]")
iF = [int(np.argmin(np.abs(FG-x))) for x in (5e-9, 1e-8, 3e-8)]; RAT = {}
for chirp in ("Mc_1e9_Msun", "non_evolving"):
    for timing, ver in (("loc", "A"), ("machine", "A"), ("machine", "B")):
        for on in ORD:
            for filt in ("ii", "iii"):
                line = []
                for n in NPH:
                    r = CUR[(chirp, timing, ver, on, n, filt)]/CUR[(chirp, timing, ver, on, n, "i")]; RAT[(chirp, timing, ver, on, filt, n)] = r
                    line.append(f"{n}: {np.median(r):.3f}")
                r = RAT[(chirp, timing, ver, on, filt, NP)]
                print(f"   {chirp:<12s} {ver} {timing:<7s} {on:<13s} filter {filt:<3s}: " + "  ".join(line) + f"   | all phased at 5/10/30 nHz: {r[iF[0]]:.3f} / {r[iF[1]]:.3f} / {r[iF[2]]:.3f}")
json.dump({"|".join(map(str, k)): [float(np.median(v)), float(v.min()), float(v.max())] for k, v in RAT.items()}, open("job54_data/sens_ratio.json", "w"))

# ---- the binaries of the first natural sky ----------------------------------------------------------------
SMOOTH[0] = False
f2, h2, tcs = CAP["f2"], CAP["h2"], CAP["tc"]; o = ORD["nearest_first"]; SN = {}
for ver in ("A", "B"):
    snr = np.zeros((len(NPH), len(f2)))
    for c0 in range(0, len(f2), 1500):
        s = slice(c0, c0+1500); e, p = parts(f2[s], tcs[s], "machine", ver, own=CAP["P0"][s]*CAP["sc"])
        for k, n in enumerate(NPH): snr[k, s] = h2[s]*np.sqrt(NP*e+p[:, o[:n]].sum(1))
    SN[ver] = snr
inb = (f2 >= 2e-9) & (f2 <= 1e-7); inb30 = (f2 >= 2e-9) & (f2 <= 3e-8)
print(f"\n== binaries of natural sky 100 above the curve (machine timing, 116 pulsars, filter ii, nearest-first), of {len(f2)} drawn sources; counts for Earth frequency 2-30 nHz (Job 51's band), then 2-100 nHz")
print("   n_phased | version A: SNR >= 10, SNR >= 3 | version B: SNR >= 10, SNR >= 3")
CNT = {}
for k, n in enumerate(NPH):
    CNT[n] = [int(((SN[v][k] >= t) & inb30).sum()) for v in ("A", "B") for t in (10, 3)] + [int(((SN[v][k] >= t) & inb).sum()) for v in ("A", "B") for t in (10, 3)]; print(f"   {n:4d} | {CNT[n][0]:5d}, {CNT[n][1]:5d} | {CNT[n][2]:4d}, {CNT[n][3]:4d}   (2-100 nHz: {CNT[n][4]}, {CNT[n][5]} | {CNT[n][6]}, {CNT[n][7]})")
ab = CAP["arr_bg"]; print(f"   sources below 2 nHz with version-B SNR >= 3 at n_ph = 0 (not counted; below the band): {int(((SN['B'][0] >= 3) & (f2 < 2e-9)).sum())}")
print(f"   check against Job 51's own number on this sky: arr_bg >= 3 for {int((ab >= 3).sum())} sources (120 pulsars, own power removed from the noise); here, version B, n_ph = 0, SNR >= 3: {CNT[0][3]}")
cidx = CAP["cand"]; mine = SN["B"][0][cidx]*np.sqrt(120/NP); rel = mine/ab; big = ab >= 1
print(f"   per-source ratio (this script, rescaled to 120 pulsars) / (Job 51 arr_bg), candidates with arr_bg >= 1: 10/50/90% {np.percentile(rel[big], [10, 50, 90]).round(3)}")
mineA = SN["A"][0][cidx]*np.sqrt(120/NP)/CAP["arr_noise"]; print(f"   same for version A against Job 51's timing-only arr: min {mineA.min():.4f}, max {mineA.max():.4f}")
# binaries above the (smoothed) version-B curve at SNR 3, Earth terms only
hc3 = np.interp(f2, FG, CUR[("Mc_1e9_Msun", "machine", "B", "nearest_first", 0, "i")]*0.3)
print(f"   binaries lying above the drawn version-B curve (SNR 3, none phased), 2-30 nHz: {int(((h2 > hc3) & inb30).sum())}; of the {int(((SN['B'][0] >= 3) & inb30).sum())} with SNR >= 3, above the drawn curve: {int(((h2 > hc3) & inb30 & (SN['B'][0] >= 3)).sum())}")
d3 = np.where((SN["B"][0] >= 3) & inb30)[0]; print(f"   detectable binaries (B, SNR >= 3): frequency {f2[d3].min()*1e9:.1f}-{f2[d3].max()*1e9:.1f} nHz, strain {h2[d3].min():.1e}-{h2[d3].max():.1e}")
top = np.argsort(-h2/(2*np.pi*f2))[:3000]; rg = np.random.default_rng(5403)
ra = rg.uniform(0, 360, len(f2)); dec = np.degrees(np.arcsin(rg.uniform(-1, 1, len(f2))))
with open("talk_data/sources_strain.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["f_earth_nHz", "strain_amplitude", "ra_deg_assigned", "dec_deg_assigned"] + [f"snr_A_nph{n}" for n in NPH] + [f"snr_B_nph{n}" for n in NPH])
    for q in top: w.writerow([f"{f2[q]*1e9:.5f}", f"{h2[q]:.5e}", f"{ra[q]:.3f}", f"{dec[q]:.3f}"] + [f"{SN['A'][k][q]:.4f}" for k in range(len(NPH))] + [f"{SN['B'][k][q]:.4f}" for k in range(len(NPH))])
print("wrote talk_data/sources_strain.csv: 3000 loudest by timing-residual amplitude h/(2 pi f); sky positions are ASSIGNED at random (the population draw has none)")
json.dump(dict(counts={str(n): v for n, v in CNT.items()}, columns=["A_snr10_2-30nHz", "A_snr3_2-30nHz", "B_snr10_2-30nHz", "B_snr3_2-30nHz", "A_snr10_2-100nHz", "A_snr3_2-100nHz", "B_snr10_2-100nHz", "B_snr3_2-100nHz"]), open("job54_data/sens_counts.json", "w"))

# ---- Figure 3 ---------------------------------------------------------------------------------------------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.8))
a1.scatter(f2[inb]*1e9, h2[inb], s=2, color="0.6", alpha=0.5, lw=0, label="binaries of the first natural sky")
dd = (SN["B"][0] >= 3) & inb30; a1.scatter(f2[dd]*1e9, h2[dd], s=22, facecolor="none", edgecolor="#d62728", lw=1, label="signal-to-noise ≥ 3 with none phased")
cols = plt.get_cmap("viridis")(np.linspace(0, 0.9, 5))
for c, n in zip(cols, (0, 1, 3, 10, NP)):
    a1.plot(FG*1e9, CUR[("Mc_1e9_Msun", "machine", "B", "nearest_first", n, "ii")]*0.3, color=c, lw=1.6, label=f"{n if n < NP else 'all 116'} phased")
a1.plot(FG*1e9, CUR[("Mc_1e9_Msun", "machine", "A", "nearest_first", 0, "i")]*0.3, "k:", lw=1, label="timing noise only, none phased")
a1.set_xscale("log"); a1.set_yscale("log"); a1.set_xlim(3, 35); a1.set_ylim(1e-16, 6e-15); a1.set_xlabel("gravitational-wave frequency [nHz]"); a1.set_ylabel("strain amplitude")
a1.set_title("strain at array signal-to-noise 3, against timing noise + background (B)", fontsize=9); a1.legend(fontsize=7, frameon=False)
for (ver, timing, ls) in (("A", "loc", ":"), ("A", "machine", "--"), ("B", "machine", "-")):
    for on, c in (("nearest_first", "#1f77b4"), ("random", "#d62728")):
        a2.plot(range(len(NPH)), [np.median(RAT[("Mc_1e9_Msun", timing, ver, on, "ii", n)]) for n in NPH], ls, color=c, marker="o", ms=3, label=f"{ver}, {timing} timing, {on.replace('_', ' ')}")
a2.set_xticks(range(len(NPH))); a2.set_xticklabels([str(n) if n < NP else "all" for n in NPH]); a2.set_xlabel("number of phased pulsars"); a2.set_ylabel("strain at SNR 10 ÷ Earth-terms-only (median over 2–100 nHz)")
a2.legend(fontsize=7, frameon=False); a2.set_ylim(0.6, 1.02)
plt.tight_layout(); plt.savefig("job54_fig3_reach.png", dpi=160); print("wrote job54_fig3_reach.png")
