"""Plot the prong-2 transition results (scaled run).

Usage:
    python plot_prong2.py [path/to/prong2_results.npz]

Defaults to ./prong2_results.npz and writes the PNG next to the npz.

Three panels:
  0. fixed per-source strain : conditional vs marginal distance info vs N (+ ∝N ref)
  1. fixed total power       : conditional vs marginal distance info vs N
  2. diagnostics             : marg/cond vs N/N* -- the (T, band) configs collapse onto
     one knee at N/N*~1 (knee tracks N* = T·Δf), and the two power modes overlap
     (marg/cond is mode-independent until the large-N break).
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

fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))

# ---- panels 0,1: conditional vs marginal, per mode ----
for col, (pfx, title) in enumerate([
        ("ps",  "Fixed per-source strain\n(Mihir independence line)"),
        ("tot", "Fixed total power\n(CW -> background fragmentation)")]):
    a = ax[col]; n = g(pfx, "n")
    a.fill_between(n, g(pfx, "cond_lo"), g(pfx, "cond_hi"), color="tab:blue", alpha=.18)
    a.plot(n, g(pfx, "cond"), "o-", color="tab:blue", lw=2, ms=4,
           label="conditional (phases known)")
    a.fill_between(n, g(pfx, "marg_lo"), g(pfx, "marg_hi"), color="tab:red", alpha=.18)
    a.plot(n, g(pfx, "marg"), "s-", color="tab:red", lw=2, ms=4,
           label="marginal (phases marginalised)")
    if pfx == "ps":
        a.plot(n, g(pfx, "cond")[0] * n, "k--", lw=1, alpha=.6, label=r"$\propto N$ reference")
    a.axvline(N_star, color="gray", ls=":", lw=1.4)
    a.annotate(r"$N^\ast=T\Delta f$", xy=(N_star, a.get_ylim()[0]),
               xytext=(N_star * 1.15, None), color="gray", fontsize=8, rotation=90,
               va="bottom")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("number of CW sources  $N$")
    a.set_ylabel(r"distance Fisher info  $I_{LL}$  [kpc$^{-2}$]")
    a.set_title(title, fontsize=10); a.grid(alpha=.25, which="both"); a.legend(fontsize=8)

# ---- panel 2: diagnostics (knee-tracking + mode-independence) ----
a = ax[2]
labels = [str(x) for x in d["diag_labels"]] if "diag_labels" in d.files else []
cmap = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(labels), 1)))
for c, lab in zip(cmap, labels):
    n = d[f"diag_{lab}_n"]; r = d[f"diag_{lab}_ratio"]; ns = float(d[f"diag_{lab}_Nstar"])
    knee = float(d[f"diag_{lab}_knee"])
    a.plot(n / ns, r, "-", color=c, lw=1.6, ms=3, marker="o", alpha=.9,
           label=fr"{lab}  ($N^*$={ns:.0f}, knee={knee:.0f})")

# overlay the two main power modes (same T, band -> same N*) to show mode-independence
for pfx, c, lab in [("ps", "k", "main: fixed per-source"),
                    ("tot", "tab:gray", "main: fixed total")]:
    n = g(pfx, "n"); r = g(pfx, "marg") / g(pfx, "cond")
    a.plot(n / N_star, r, "--", color=c, lw=1.8, alpha=.8, label=lab)

mdev = float(d["mode_indep_maxdev"]) if "mode_indep_maxdev" in d.files else float("nan")
# empirical knee location in units of N*: median over the diagnostic configs
kr = np.nan
if labels:
    kr = float(np.median([float(d[f"diag_{l}_knee"]) / float(d[f"diag_{l}_Nstar"])
                          for l in labels]))
a.axvline(1.0, color="gray", ls=":", lw=1.2)
a.text(1.05, 0.06, r"naive $N^\ast=T\Delta f$", color="gray", fontsize=7.5,
       rotation=90, va="bottom")
if np.isfinite(kr):
    a.axvline(kr, color="crimson", ls="--", lw=1.4)
    a.text(kr * 1.05, 0.92, fr"empirical knee $\approx{kr:.0f}\,N^\ast$",
           color="crimson", fontsize=8, rotation=90, va="top")
a.axhline(0.5, color="k", lw=.6, alpha=.5)
a.set_xscale("log")
a.set_xlabel(r"$N / N^\ast$   ($N^\ast = T\,\Delta f$)")
a.set_ylabel("marginal / conditional")
a.set_title(r"Knee $\propto N^\ast$ (const knee/$N^\ast$ vs $T$); "
            "modes overlap exactly\n"
            fr"(max$|\Delta$(marg/cond)$|$={mdev:.3f}; array offset $\approx{kr:.0f}\times$)",
            fontsize=9.5)
a.set_ylim(0, 1.05); a.grid(alpha=.25, which="both"); a.legend(fontsize=7, loc="lower left")

fig.suptitle("Prong 2 (info-only): distance information across the resolvable->confusion "
             "transition  (116 psr, 15 yr weekly)", fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(out_png, dpi=130, bbox_inches="tight")
print("saved", out_png)
