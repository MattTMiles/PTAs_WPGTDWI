"""JOB 46c -- Panel A: strain vs frequency with arcs for the 10 loudest; Panel B: in-band ruler ladder, retardation ON vs OFF."""
import sys, numpy as np, matplotlib
ARCS = sys.argv[1] if len(sys.argv) > 1 else "all"
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
d = np.load("partB_sky.npz"); fE, hE, fP, sea = d["f_src"], d["h_src"], d["fp41"], d["S_faint"]
T = float(d["T_OBS"]); YR = 3.15576e7; edges = d["bin_edges"]; nb = len(edges)-1
hP = hE*(fP/fE)**(2/3)
band = (15e-9, 40e-9); inb = lambda f: (f >= band[0]) & (f < band[1])
S_PI = 1/(1+np.pi**2)
def shares(f, h):
    P = h**2*f*T; b = np.digitize(f, edges)-1
    tot = np.bincount(b, weights=P, minlength=nb) + sea
    return P/tot[b]
sON = shares(fP, hP)            # retardation ON: pulsar terms at retarded f
sOFF = shares(fE, hE)           # retardation OFF: tau = 0, pulsar terms at Earth f (same as Earth-term shares)
iON = np.where(inb(fP))[0]; iOFF = np.where(inb(fE))[0]
rON = iON[np.argsort(sON[iON])[::-1]]; rOFF = iOFF[np.argsort(sOFF[iOFF])[::-1]]
ladON, ladOFF = sON[rON], sOFF[rOFF]
nON, nOFF = int((ladON > S_PI).sum()), int((ladOFF > S_PI).sum())
sig = lambda s: np.sqrt((1-s)/s)
medON, medOFF = np.median(sig(ladON[:20])), np.median(sig(ladOFF[:20]))
rankON = np.argmax(ladON < S_PI)+1 if (ladON < S_PI).any() else None
rankOFF = np.argmax(ladOFF < S_PI)+1 if (ladOFF < S_PI).any() else None
loudrank = np.empty(len(hE), int); loudrank[np.argsort(hE)[::-1]] = np.arange(1, len(hE)+1)   # 1 = loudest Earth strain
top10 = np.argsort(hE)[::-1][:10] if ARCS == "all" else np.where(fE >= 15e-9)[0][np.argsort(hE[fE >= 15e-9])[::-1][:10]]
print(f"in-band pulsar terms: ON {len(iON)}  OFF {len(iOFF)}")
print(f"usable rulers (s > {S_PI:.3f}) in band: ON = {nON}, OFF = {nOFF}; first rank below line: ON {rankON}, OFF {rankOFF}")
print(f"median sigma_phi over top-20: ON {medON:.2f} rad, OFF {medOFF:.2f} rad;  s>=0.5 in band: ON {(ladON>=.5).sum()} OFF {(ladOFF>=.5).sum()}")
print("ON ladder top 10  (loudness rank by h_E, f_E->f_P nHz, s):")
for k, j in enumerate(rON[:10]): print(f"   rank {k+1:2d}: loudness #{loudrank[j]:5d}  f_E {fE[j]*1e9:5.1f} -> f_P {fP[j]*1e9:5.1f}  s={ladON[k]:.3f}")
print("OFF ladder top 10 (loudness rank by h_E, f_E nHz, s):")
for k, j in enumerate(rOFF[:10]): print(f"   rank {k+1:2d}: loudness #{loudrank[j]:5d}  f_E {fE[j]*1e9:5.1f}  s={ladOFF[k]:.3f}")
print(f"Panel A arcs ({ARCS}): loudness rank, f_E -> f_P, s_ON, s_OFF, in band ON?")
if False: print(": f_E -> f_P nHz, s_ON, s_OFF:")
for j in top10: print(f"   loudness #{loudrank[j]:2d}: {fE[j]*1e9:5.1f} -> {fP[j]*1e9:5.1f}   s_ON={sON[j]:.3f}  s_OFF={sOFF[j]:.3f}   in band ON? {bool(inb(fP[j]))}")
print(f"above 50 nHz: Earth {(fE>=50e-9).sum()} pulsar {(fP>=50e-9).sum()}")
print(f"ON ladder top-10 members among the arc set: {len(set(rON[:10]) & set(top10))}; OFF top-10 among arc set: {len(set(rOFF[:10]) & set(top10))}; among the 100 loudest: {int((loudrank[rON[:10]]<=100).sum())}; loudness ranks of ON top-10: {sorted(loudrank[rON[:10]].tolist())}")

GREY, BLUE, BAND = "0.55", "#1f5fa8", "#f3e3b5"
plt.rcParams.update({"font.size": 9.5, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw=dict(wspace=0.25))
# Panel A
a.axvspan(band[0]*1e9, band[1]*1e9, color=BAND, alpha=0.6, lw=0, zorder=0)
a.loglog(fE*1e9, hE, "o", mfc="none", mec=GREY, ms=2.6, mew=0.4, label="Earth term", zorder=2, rasterized=True)
a.loglog(fP*1e9, hP, "o", color=BLUE, ms=2.4, mec="none", alpha=0.7, label="pulsar term", zorder=3, rasterized=True)
for j in top10:
    a.add_patch(FancyArrowPatch((fE[j]*1e9, hE[j]), (fP[j]*1e9, hP[j]), connectionstyle="arc3,rad=0.3", arrowstyle="-|>,head_length=4,head_width=2.5", color="#c8102e", lw=1.0, zorder=4, shrinkA=2, shrinkB=2))
a.plot([], [], color="#c8102e", lw=1.0, label="10 loudest sources: Earth term → pulsar term" if ARCS == "all" else "10 loudest sources with Earth term ≥ 15 nHz: Earth term → pulsar term")
a.set_xlim(0.9, 110); a.set_xlabel("gravitational-wave frequency [nHz]"); a.set_ylabel("strain amplitude h")
a.text(np.sqrt(band[0]*band[1])*1e9, hE.min()*1.0, "PTA sensitivity /\ncomb band\n15–40 nHz", ha="center", va="bottom", fontsize=8)
a.set_ylim(hE.min()*0.7, hE.max()*2.2)
a.text(0.02, 0.03, f"above 50 nHz: Earth {(fE>=50e-9).sum()}, pulsar {(fP>=50e-9).sum()}", transform=a.transAxes, fontsize=8.5)
from matplotlib.lines import Line2D
handles = [Line2D([], [], marker="o", ls="none", mfc="none", mec=GREY, ms=5, label="Earth term"), Line2D([], [], marker="o", ls="none", color=BLUE, ms=5, label="pulsar term"), Line2D([], [], color="#c8102e", lw=1.0, label=("10 loudest sources: Earth term → pulsar term" if ARCS == "all" else "10 loudest sources with Earth term ≥ 15 nHz:\nEarth term → pulsar term"))]
a.legend(handles=handles, loc="upper right", fontsize=8)
a.set_title("where the loud sources go", fontsize=9.5, loc="left")
# Panel B
R = 50; x = np.arange(1, R+1)
b.semilogy(x, ladOFF[:R], "o", mfc="none", mec=GREY, ms=5, mew=1.0, label="retardation OFF (τ = 0: pulsar terms at their Earth frequencies)")
b.semilogy(x, ladON[:R], "o", color=BLUE, ms=5, label="retardation ON (the real sky)")
b.axhline(0.5, color="0.3", ls="--", lw=0.9); b.text(R, 0.5*1.1, "s = 0.5:  σ_φ = 1 rad", ha="right", fontsize=8, color="0.3")
b.axhline(S_PI, color="0.3", ls=":", lw=0.9); b.text(R/2, S_PI/1.7, "s = 1/(1+π²) ≈ 0.09:  σ_φ = π, comb unreadable", ha="center", va="top", fontsize=8, color="0.3")
b.set_xlim(0, R+1); b.set_ylim(1.5e-2, 1.3); b.set_xlabel("rank of the pulsar term within 15–40 nHz, sorted by its share s (1 = largest)")
b.set_ylabel("s = the pulsar term's share of the total GW power in its 1/T bin")
b.text(0.02, 0.03, f"usable rulers (s > 0.09) in band:  ON = {nON},  OFF = {nOFF}\nmedian σ_φ over the top 20:  ON = {medON:.2f} rad,  OFF = {medOFF:.2f} rad", transform=b.transAxes, fontsize=8.5, va="bottom")
b.legend(loc="upper right", fontsize=8)
b.set_title("the pulsar terms that reach the band, ranked by how much of their bin they own —\nwith and without the pile-up", fontsize=9.5, loc="left")
b.grid(alpha=0.2, lw=0.4, axis="y")
fig.text(0.5, -0.04, f"frozen simulation-machine sky (active_wheel_demo_v4, draw_population seed 3); pulsar 41 at L = 2105 pc; bins = 1/T at T = {T/YR:.0f} yr; s includes the faint sea in every bin's total.\nPanel A shows the {len(fE):,} sources above the hybrid split (its per-bin strain floor is the sawtooth lower envelope); the faint sea is not drawn.", ha="center", fontsize=7.5, color="0.35")
fig.savefig("/home/mattm/projects/HSYMT/PAPER_figures/job46_pileup_histogram.png" if ARCS == "all" else "/home/mattm/projects/HSYMT/PAPER_figures/job46_pileup_histogram_arcs15.png", dpi=200, bbox_inches="tight"); print("saved")
