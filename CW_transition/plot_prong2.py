"""Plot the prong-2 transition results (scaled run, v2).

Usage:
    python plot_prong2.py [path/to/prong2_results.npz]

Six panels (2x3):
  A. fixed per-source : conditional vs marginal info vs N, with population overlay (+∝N)
  B. fixed total power : conditional vs marginal info vs N
  C. marg/cond vs N    : uniform ps & tot overlap (mode-independent) + population departs
  D. knee diagnostic   : marg/cond vs N/N* across (T,band) -> knee tracks N* (offset ~50x)
  E. sigma_L (pc)      : J0437-like pulsar, fixed_total (degrades) + fixed_persource;
                         sub-parsec and EM-prior lines
  F. rcond sweep       : marg/cond for rcond in {1e-8,1e-10,1e-12} -> knee stability
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

npz_path = sys.argv[1] if len(sys.argv) > 1 else "prong2_results.npz"
out_png = os.path.join(os.path.dirname(os.path.abspath(npz_path)),
                       "prong2_distance_information_transition.png")

d = np.load(npz_path, allow_pickle=True)
def g(pfx, k): return d[f"{pfx}_{k}"]
N_star = float(d["ps_N_star"])
knee = float(d["knee_N"]) if "knee_N" in d.files else np.nan

fig, ax = plt.subplots(2, 3, figsize=(17, 9))

# ---- A: fixed per-source cond/marg + population overlay ----
a = ax[0, 0]; n = g("ps", "n")
a.fill_between(n, g("ps", "cond_lo"), g("ps", "cond_hi"), color="tab:blue", alpha=.15)
a.plot(n, g("ps", "cond"), "o-", color="tab:blue", lw=2, ms=4, label="conditional")
a.fill_between(n, g("ps", "marg_lo"), g("ps", "marg_hi"), color="tab:red", alpha=.15)
a.plot(n, g("ps", "marg"), "s-", color="tab:red", lw=2, ms=4, label="marginal")
if "pop_n" in d.files:
    a.plot(g("pop", "n"), g("pop", "marg"), "^--", color="tab:orange", lw=1.6, ms=4,
           label="marginal (population)")
a.plot(n, g("ps", "cond")[0] * n, "k--", lw=1, alpha=.5, label=r"$\propto N$")
a.axvline(N_star, color="gray", ls=":", lw=1.2)
if np.isfinite(knee):
    a.axvline(knee, color="crimson", ls=":", lw=1.2)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("N"); a.set_ylabel(r"distance info $I_{LL}$ [kpc$^{-2}$]")
a.set_title("Fixed per-source strain", fontsize=10); a.grid(alpha=.25, which="both")
a.legend(fontsize=7.5)

# ---- B: fixed total cond/marg ----
a = ax[0, 1]; n = g("tot", "n")
a.fill_between(n, g("tot", "cond_lo"), g("tot", "cond_hi"), color="tab:blue", alpha=.15)
a.plot(n, g("tot", "cond"), "o-", color="tab:blue", lw=2, ms=4, label="conditional")
a.fill_between(n, g("tot", "marg_lo"), g("tot", "marg_hi"), color="tab:red", alpha=.15)
a.plot(n, g("tot", "marg"), "s-", color="tab:red", lw=2, ms=4, label="marginal")
a.axvline(N_star, color="gray", ls=":", lw=1.2)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("N"); a.set_ylabel(r"distance info $I_{LL}$ [kpc$^{-2}$]")
a.set_title("Fixed total power", fontsize=10); a.grid(alpha=.25, which="both")
a.legend(fontsize=7.5)

# ---- C: marg/cond, mode-independence + population departure ----
a = ax[0, 2]
a.plot(g("ps", "n"), g("ps", "marg") / g("ps", "cond"), "o-", color="tab:blue", lw=2,
       ms=4, label="fixed per-source")
a.plot(g("tot", "n"), g("tot", "marg") / g("tot", "cond"), "s--", color="tab:gray",
       lw=1.6, ms=4, label="fixed total (overlaps)")
if "pop_n" in d.files:
    a.plot(g("pop", "n"), g("pop", "marg") / g("pop", "cond"), "^-", color="tab:orange",
           lw=2, ms=4, label="population (k=3 loud)")
    dn = int(d["pop_depart_N"]) if "pop_depart_N" in d.files else None
    if dn and dn > 0:
        a.axvline(dn, color="tab:orange", ls=":", lw=1.2)
        a.text(dn * 1.05, 0.1, f"departs N={dn}", color="tab:orange", fontsize=7.5,
               rotation=90, va="bottom")
mdev = float(d["mode_indep_maxdev"])
a.axvline(N_star, color="gray", ls=":", lw=1.2); a.axhline(0.5, color="k", lw=.6, alpha=.4)
a.set_xscale("log"); a.set_xlabel("N"); a.set_ylabel("marginal / conditional")
a.set_title(f"marg/cond: modes overlap (Δ={mdev:.3f}); population breaks it", fontsize=9.5)
a.set_ylim(0, 1.05); a.grid(alpha=.25, which="both"); a.legend(fontsize=7.5)

# ---- D: knee diagnostic collapse ----
a = ax[1, 0]
labels = [str(x) for x in d["diag_labels"]] if "diag_labels" in d.files else []
cmap = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(labels), 1)))
for c, lab in zip(cmap, labels):
    nn = d[f"diag_{lab}_n"]; r = d[f"diag_{lab}_ratio"]; ns = float(d[f"diag_{lab}_Nstar"])
    a.plot(nn / ns, r, "-o", color=c, lw=1.5, ms=3, alpha=.9,
           label=fr"{lab} ($N^*$={ns:.0f})")
kr = np.nan
if labels:
    kr = float(np.median([float(d[f"diag_{l}_knee"]) / float(d[f"diag_{l}_Nstar"])
                          for l in labels]))
a.axvline(1.0, color="gray", ls=":", lw=1.2)
if np.isfinite(kr):
    a.axvline(kr, color="crimson", ls="--", lw=1.3)
    a.text(kr * 1.05, 0.9, fr"knee$\approx{kr:.0f}N^*$", color="crimson", fontsize=8,
           rotation=90, va="top")
a.axhline(0.5, color="k", lw=.6, alpha=.4)
a.set_xscale("log"); a.set_xlabel(r"$N/N^*$ ($N^*=T\Delta f$)")
a.set_ylabel("marginal / conditional")
a.set_title(r"Knee $\propto N^*$ (T-scaling); array offset $\sim50\times$", fontsize=9.5)
a.set_ylim(0, 1.05); a.grid(alpha=.25, which="both"); a.legend(fontsize=6.5, loc="lower left")

# ---- E: sigma_L in pc for J0437-like ----
a = ax[1, 1]
if "j_n" in d.files:
    jn = d["j_n"]
    a.plot(jn, d["j_sigma_pc"], "o-", color="tab:red", lw=2, ms=4,
           label="fixed total (degrades)")
    if "j_sigma_pc_ps" in d.files:
        a.plot(jn, d["j_sigma_pc_ps"], "s--", color="tab:blue", lw=1.6, ms=4,
               label="fixed per-source")
    em = float(d["em_sigma_pc"])
    a.axhline(1.0, color="k", ls=":", lw=1.2); a.text(jn[0], 1.05, "sub-parsec", fontsize=7.5)
    a.axhline(em, color="tab:green", ls="--", lw=1.2)
    a.text(jn[0], em * 1.1, f"EM prior {em} pc", color="tab:green", fontsize=7.5)
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("N"); a.set_ylabel(r"$\sigma_L = 1/\sqrt{I_{LL}^{\rm marg}}$  [pc]")
a.set_title("J0437-like distance precision\n(absolute scale schematic)", fontsize=9.5)
a.grid(alpha=.25, which="both"); a.legend(fontsize=7.5)

# ---- F: rcond sweep ----
a = ax[1, 2]
rcv = d["rcond_values"] if "rcond_values" in d.files else []
cmap2 = ["tab:purple", "tab:green", "tab:brown"]
for c, rc in zip(cmap2, rcv):
    key = f"ps_marg_rcond_{rc:.0e}"
    if key in d.files:
        r = d[key] / d["ps_cond"]
        a.plot(d["ps_n"], r, "-o", color=c, lw=1.6, ms=3, label=f"rcond={rc:.0e}")
a.axhline(0.5, color="k", lw=.6, alpha=.4)
spread = float(d["rcond_knee_spread"]) * 100 if "rcond_knee_spread" in d.files else np.nan
a.set_xscale("log"); a.set_xlabel("N"); a.set_ylabel("marginal / conditional")
a.set_title(f"rcond sweep: knee stable (spread {spread:.1f}%)", fontsize=9.5)
a.set_ylim(0, 1.05); a.grid(alpha=.25, which="both"); a.legend(fontsize=7.5)

fig.suptitle("Prong 2 (info-only): CW->GWB distance-information transition  "
             "(116 psr, 15 yr weekly)", fontsize=12, y=1.0)
fig.tight_layout()
fig.savefig(out_png, dpi=125, bbox_inches="tight")
print("saved", out_png)
