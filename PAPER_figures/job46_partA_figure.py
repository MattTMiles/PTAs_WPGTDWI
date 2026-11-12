"""JOB 46 Part A v2 -- three separate bare plots from prong2_results.npz (no annotations)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = np.load("/home/mattm/projects/HSYMT/CW_transition/prong2_results.npz", allow_pickle=True)
n = d["ps_n"]; Nstar = float(d["N_star"]); knee = float(d["knee_N"]); npsr = int(d["ps_n_psr"]); nsky = int(d["ps_n_seeds"])
BLUE, RED = "#1f5fa8", "#c8102e"
OUT = "/home/mattm/projects/HSYMT/PAPER_figures/"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "legend.frameon": False})
XL = f"number of CW sources on the sky, N  (array fixed at N_psr = {npsr})"
def verticals(a):
    a.axvline(Nstar, color="gray", ls=":", lw=1.3, label="one source per frequency bin, N = TΔf")
    a.axvline(knee, color=RED, ls=":", lw=1.3, label=f"N_knee = {knee/(npsr*Nstar):.2f} N_psr TΔf")
def info_plot(pfx, fname):
    fig, a = plt.subplots(figsize=(6.4, 5.2))
    c, cl, ch = d[pfx+"_cond"], d[pfx+"_cond_lo"], d[pfx+"_cond_hi"]
    m, ml, mh = d[pfx+"_marg"], d[pfx+"_marg_lo"], d[pfx+"_marg_hi"]
    a.fill_between(n, cl, ch, color=BLUE, alpha=0.2, lw=0)
    a.plot(n, c, "o-", color=BLUE, lw=2, ms=4.5, label="conditional information")
    a.fill_between(n, ml, mh, color=RED, alpha=0.2, lw=0)
    a.plot(n, m, "s--", color=RED, lw=2, ms=4.5, label="marginal information")
    verticals(a)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(XL)
    a.set_ylabel("distance information per pulsar, $I_{LL} = 1/\\sigma_L^2$  [kpc$^{-2}$]")
    sec = a.secondary_yaxis("right", functions=(lambda I: 1e3/np.sqrt(I), lambda s: (1e3/s)**2))
    sec.set_ylabel("distance precision $\\sigma_L$  [pc]")
    a.grid(alpha=0.25, which="both")
    a.legend(fontsize=8.5, loc="lower right" if pfx == "ps" else "lower left")
    fig.savefig(OUT+fname, dpi=200, bbox_inches="tight"); print("saved", fname)
info_plot("ps", "job46_info_fixed_per_source_strain.png")
info_plot("tot", "job46_info_fixed_total_power.png")
fig, a = plt.subplots(figsize=(6.4, 5.2))
a.plot(n, d["ps_marg"]/d["ps_cond"], "o-", color=BLUE, lw=2, ms=5, label="fixed per-source strain")
a.plot(n, d["tot_marg"]/d["tot_cond"], "s--", color="#7a3fa0", lw=1.6, ms=6, mfc="none", label="fixed total power")
a.axhline(0.5, color="k", lw=0.8, alpha=0.5)
verticals(a)
a.set_xscale("log"); a.set_ylim(0, 1.05); a.set_xlabel(XL)
a.set_ylabel("surviving fraction of distance information,  $I_{\\rm marg}/I_{\\rm cond}$")
a.grid(alpha=0.25, which="both"); a.legend(fontsize=9, loc="lower left")
a.text(0.03, 0.55, "the two normalisations coincide:\nthe surviving fraction depends only on\nhow many sources share the band,\nnot on their loudness", transform=a.transAxes, fontsize=8.5, va="bottom")
fig.savefig(OUT+"job46_info_surviving_fraction.png", dpi=200, bbox_inches="tight"); print("saved surviving")
