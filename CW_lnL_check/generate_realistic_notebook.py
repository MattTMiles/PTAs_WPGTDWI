"""Generate notebook 05: realistic first-detection CW parameters.

- log10_h drawn per source from [-14.3, -13.8] (GWB floor to loud)
- log10_fgw drawn per source from [-8.3, -7.7] (5-20 nHz)
- log10_mc drawn per source from [8.5, 9.5]
- GWB injected at NG15 level: log10_A=-14.6, gamma=13/3
- White noise equad = -6 (1 us, realistic)
- Fourier components reduced to 14 (major Cholesky speedup vs 30)
- N_CW in [1, 2, 3, 5, 10] (realistic early-detection numbers)
"""
import json
import uuid


def md(source):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [source],
        "id": uuid.uuid4().hex[:8],
    }


def code(source):
    return {
        "cell_type": "code",
        "metadata": {},
        "source": [source],
        "outputs": [],
        "execution_count": None,
        "id": uuid.uuid4().hex[:8],
    }


def notebook(cells):
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3 (discotech)",
                "language": "python",
                "name": "discotech",
            },
            "language_info": {"name": "python", "version": "3.12.0"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def write_nb(path, cells):
    with open(path, "w") as f:
        json.dump(notebook(cells), f, indent=1)
    print(f"Wrote {path}")


cells = [
    md(r"""# Realistic First-Detection CWs: Distance Landscape vs N_CW

**Setup (tuned to what PTAs actually expect to detect first):**

- **CW strain per source:** `log10_h ~ Uniform(-14.0, -13.5)`
  - Lower bound ≈ single-CW detection threshold (NG 12.5-yr h95 ~ 9×10⁻¹⁵ at 10 nHz).
    Unlike the GWB, a monochromatic CW has no frequency integration, so its
    detectable floor sits ~3-4× above the GWB characteristic strain.
  - Upper bound = moderately loud (SNR ~10-20 in realistic noise)
- **CW frequency per source:** `log10_fgw ~ Uniform(-8.0, -7.5)` (10-32 nHz, PTA CW sensitivity bucket)
- **CW chirp mass per source:** `log10_mc ~ Uniform(8.5, 9.5)` (massive SMBHBs)
- **Sky position, inclination, phase, psi:** uniform on sphere / uniform
- **GWB injected:** `log10_A = -14.6, gamma = 13/3` (NG15 central value, SMBHB spectrum)
- **White noise:** `log10_equad = -6` (1 μs, typical IPTA precision)
- **Fourier components:** 14 (still covers up to ~30 nHz in 15-yr data)

**Question:** Under *realistic* first-detection conditions, does adding more CW sources
sharpen the pulsar-distance constraint? This mirrors notebook 02 but with physically
motivated source parameters instead of fixed strain and generic scenarios.
"""),
    code(r"""import sys, importlib
sys.path.insert(0, ".")
import cw_helpers
importlib.reload(cw_helpers)
from cw_helpers import *
import matplotlib.pyplot as plt
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

%matplotlib inline
plt.rcParams.update({"figure.dpi": 120, "font.size": 11})
"""),
    md("## 1. Load pulsars"),
    code(r"""Npulsars = 116  # all available pulsars
ent_psrs, disco_psrs = load_pulsars(Npulsars)
print(f"Loaded {len(ent_psrs)} pulsars")
"""),
    md(r"""## 2. Configuration
"""),
    code(r"""seed = 1234
ncw_values = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

LOG10_H_RANGE   = (-14.0, -13.5)   # single-CW detectable (above NG 12.5yr h95 ~ 9e-15)
LOG10_FGW_RANGE = (-8.0, -7.5)     # 10-32 nHz, peak PTA CW sensitivity bucket
LOG10_MC_RANGE  = (8.5, 9.5)

GWB_LOG10_A = -14.6
GWB_GAMMA   = 13.0 / 3.0
LOG10_EQUAD = -6.0    # realistic 1 us white noise
N_COMPONENTS = 14     # Fourier modes for GP marginalisation
"""),
    md(r"""## 3. Run N_CW scan with realistic draws

For each N_CW:
1. Build enterprise PTA at reduced components and realistic equad
2. Draw CW params per source from the realistic ranges
3. Simulate (includes realistic GWB)
4. Build discovery likelihood with matching noise/components
5. Compute joint Δln L (truth vs best wrong mode in EM prior)
6. Also scan the first pulsar's distance for landscape plotting
"""),
    code(r"""all_results = {}

for ncw in ncw_values:
    print(f"\n{'='*60}\nN_CW = {ncw}\n{'='*60}")

    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    # Build PTA with realistic equad and reduced components
    pta, cw_block_names, Tspan = build_enterprise_pta(
        ent_psrs, ncw,
        components=N_COMPONENTS,
        log10_equad=LOG10_EQUAD,
    )

    # Realistic per-source draws + NG15 GWB
    inj_params = generate_injection_params(
        pta, ent_psrs, ncw, cw_block_names,
        log10_h=None,  # unused for "realistic"
        scenario="realistic",
        rng=rng,
        log10_h_range=LOG10_H_RANGE,
        log10_fgw_range=LOG10_FGW_RANGE,
        log10_mc_range=LOG10_MC_RANGE,
        gwb_log10_A=GWB_LOG10_A,
        gwb_gamma=GWB_GAMMA,
    )

    # Report realized source parameters
    for name in cw_block_names:
        h_s = inj_params[f"{name}_log10_h"]
        f_s = inj_params[f"{name}_log10_fgw"]
        m_s = inj_params[f"{name}_log10_mc"]
        print(f"  {name}: log10_h={h_s:.2f}, log10_fgw={f_s:.2f}, log10_mc={m_s:.2f}")

    # Simulate (with GWB at realistic level)
    sim_resids = simulate(pta, inj_params)
    resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}

    # Build discovery likelihood with cached-Cholesky fast path (GP hyperparams
    # frozen to truth; GWB+red noise fixed, CW params varied). This is required
    # at N_psr=116 because the slow path hits a jax tracer issue in jnp.block
    # at large array sizes, and it is also ~10-30x faster per evaluation.
    logl_fn, param_keys, base_vals = build_fast_scan_likelihood(
        disco_psrs, resid_map, ncw, inj_params, cw_block_names,
        components=N_COMPONENTS,
        log10_equad=LOG10_EQUAD,
    )

    ll_truth = float(logl_fn(base_vals))
    print(f"  ln L(truth) = {ll_truth:.2f}")

    # Scan first pulsar's distance (for landscape)
    psr_e = ent_psrs[0]
    psr_d = disco_psrs[0]
    true_dist = psr_e.pdist[0]
    dist_key = f"{psr_e.name}_cw_p_dist"
    disco_suffixes = ["" if i == 0 else f"_{i+1}" for i in range(ncw)]

    all_dL = []
    for i, suffix in enumerate(disco_suffixes):
        ent_name = cw_block_names[i]
        dL_i = compute_mode_spacing(
            inj_params[f"{ent_name}_cos_gwtheta"],
            inj_params[f"{ent_name}_gwphi"],
            inj_params[f"{ent_name}_log10_fgw"],
            psr_d.pos,
        )
        all_dL.append(dL_i)
    min_dL = min(all_dL)
    max_dL = max(all_dL)

    n_modes_scan = 10
    scan_half = n_modes_scan * max_dL
    scan_min = max(0.01, true_dist - scan_half)
    scan_max = true_dist + scan_half

    scan_d, scan_ll = scan_pulsar_distance(
        logl_fn, base_vals, param_keys, dist_key,
        scan_min, scan_max, n_points=3000,
        n_components=N_COMPONENTS,
    )
    landscape = analyze_peaks(scan_d, scan_ll, true_dist, mode_spacing=min_dL)
    landscape["scan_d"] = scan_d
    landscape["scan_ll"] = scan_ll
    landscape["true_dist"] = true_dist
    landscape["min_dL"] = min_dL
    landscape["max_dL"] = max_dL

    # Joint best wrong mode within EM prior (all pulsars)
    joint = compute_joint_best_wrong_in_prior(
        logl_fn, base_vals, param_keys,
        ent_psrs, disco_psrs, inj_params,
        cw_block_names=cw_block_names, n_sigma=3.0,
        n_components=N_COMPONENTS,
    )

    all_results[ncw] = {
        "landscape_first_psr": landscape,
        "ll_truth": ll_truth,
        "ll_best_wrong": joint["ll_best_wrong"],
        "delta_lnL_joint": joint["delta_lnL"],
        "joint_per_pulsar": joint["per_pulsar"],
        "inj_params": inj_params,
    }
    print(f"  Joint Δln L (truth vs best wrong in EM prior): {joint['delta_lnL']:.2f}")

    jax.clear_caches()

print("\nAll runs complete.")
"""),
    md(r"""## 4. First-pulsar distance landscape vs N_CW
"""),
    code(r"""first_psr_name = ent_psrs[0].name
fig, axes = plt.subplots(len(ncw_values), 1, figsize=(14, 3 * len(ncw_values)), sharex=False)

for ax_idx, ncw in enumerate(ncw_values):
    ax = axes[ax_idx]
    r = all_results[ncw]["landscape_first_psr"]
    scan_d, scan_ll = r["scan_d"], r["scan_ll"]
    baseline = r["truth_logl"]

    ax.plot(scan_d, scan_ll - baseline, "b-", lw=0.5, alpha=0.8)

    if len(r["peak_distances"]) > 0:
        truth_mask = np.abs(r["peak_distances"] - r["true_dist"]) < 0.3 * r["min_dL"]
        wrong_mask = ~truth_mask
        if np.any(wrong_mask):
            ax.plot(r["peak_distances"][wrong_mask],
                    r["peak_heights"][wrong_mask] - baseline,
                    "o", color="orange", ms=3, alpha=0.6)
        if np.any(truth_mask):
            ax.plot(r["peak_distances"][truth_mask],
                    r["peak_heights"][truth_mask] - baseline,
                    "o", color="green", ms=5)

    ax.axvline(r["true_dist"], color="red", ls="--", lw=0.8, alpha=0.6)
    delta = r["delta_lnL_truth_vs_best_wrong"]
    delta_str = f"{delta:.1f}" if not np.isnan(delta) else "N/A"
    ax.set_ylabel(r"$\Delta \ln L$")
    ax.set_title(
        f"N_CW = {ncw}  |  peaks = {r['n_peaks']}  |  "
        f"dL range [{r['min_dL']:.3f}, {r['max_dL']:.3f}] kpc  |  "
        r"$\Delta \ln L$(truth $-$ best wrong) = " + delta_str,
        fontsize=10,
    )
    ax.grid(True, alpha=0.2)

axes[-1].set_xlabel("p_dist [kpc]")
fig.suptitle(f"Realistic CWs: {first_psr_name} distance landscape", fontsize=13, y=1.01)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 5. Joint Δln L vs N_CW
"""),
    code(r"""joint_deltas = [all_results[n]["delta_lnL_joint"] for n in ncw_values]

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(ncw_values, joint_deltas, "s-", lw=2, ms=10, color="steelblue")
ax.set_xlabel("Number of CW sources (realistic draws)")
ax.set_ylabel(r"Joint $\Delta \ln L$ (truth $-$ best wrong in EM prior)")
ax.set_title("Distance importance vs N_CW under realistic first-detection conditions")
ax.axhline(0, color="k", ls="--", lw=0.5, alpha=0.5)
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 6a. Diagnostic: why is Δln L what it is?

For each N_CW, break down per-pulsar behaviour without needing the full distance
landscapes. We look at:

- **distribution of per-pulsar Δln L** — if most pulsars sit at Δln L > 0, truth is the
  per-pulsar ll max and stacking should be constructive. If many are ≤ 0, the truth is
  not the local max (noise-absorption regime).
- **best-wrong offset in units of mode-spacing** — are alternate peaks at expected
  phase-wrap distances (integer multiples of dL), or clustered near truth (noise)?
- **peak height ratio** — how much do wrong peaks beat truth by, on average?
- **n_peaks_in_prior** — how many modes live in ±3σ. Growing with N_CW means more
  chances to noise-fit.
"""),
    code(r"""fig, axes = plt.subplots(2, 2, figsize=(14, 10))

cmap = plt.cm.viridis
colors = [cmap(i / max(1, len(ncw_values)-1)) for i in range(len(ncw_values))]

# (a) Histogram of per-pulsar Δln L across pulsars for each N_CW
ax = axes[0, 0]
for i, ncw in enumerate(ncw_values):
    pp = all_results[ncw]["joint_per_pulsar"]
    deltas = np.array([v["delta_lnL"] for v in pp.values() if not np.isnan(v["delta_lnL"])])
    ax.hist(deltas, bins=25, alpha=0.4, color=colors[i], label=f"N_CW={ncw}",
            histtype="step", lw=2)
ax.axvline(0, color="k", ls="--", lw=0.8)
ax.set_xlabel(r"Per-pulsar $\Delta \ln L$ (truth $-$ best wrong)")
ax.set_ylabel("# pulsars")
ax.set_title("Per-pulsar Δln L: is truth the local max?")
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)

# (b) Best-wrong offset in mode-spacing units
ax = axes[0, 1]
for i, ncw in enumerate(ncw_values):
    pp = all_results[ncw]["joint_per_pulsar"]
    offs = np.array([v["best_wrong_offset_modes"] for v in pp.values()
                     if not np.isnan(v["best_wrong_offset_modes"])])
    if len(offs) > 0:
        ax.hist(offs, bins=np.arange(-10, 11, 0.5), alpha=0.4, color=colors[i],
                label=f"N_CW={ncw}", histtype="step", lw=2)
ax.axvline(0, color="k", ls="--", lw=0.8)
ax.set_xlabel("Best-wrong offset / mode-spacing dL")
ax.set_ylabel("# pulsars")
ax.set_title("Where does the alternate peak sit? (integer = phase-wrap mode)")
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)

# (c) n_peaks_in_prior vs N_CW (box/violin across pulsars)
ax = axes[1, 0]
n_peaks_data = []
for ncw in ncw_values:
    pp = all_results[ncw]["joint_per_pulsar"]
    n_peaks_data.append([v["n_peaks_in_prior"] for v in pp.values()])
parts = ax.boxplot(n_peaks_data, positions=ncw_values, widths=0.5,
                   patch_artist=True, showmeans=True)
for patch, c in zip(parts["boxes"], colors):
    patch.set_facecolor(c); patch.set_alpha(0.5)
ax.set_xlabel("N_CW")
ax.set_ylabel("# peaks in ±3σ EM prior")
ax.set_title("Mode density vs N_CW (more modes = more chances to noise-fit)")
ax.grid(True, alpha=0.3)

# (d) Fraction of pulsars with Δln L ≤ 0 vs N_CW
ax = axes[1, 1]
frac_neg = []
median_delta = []
for ncw in ncw_values:
    pp = all_results[ncw]["joint_per_pulsar"]
    deltas = np.array([v["delta_lnL"] for v in pp.values() if not np.isnan(v["delta_lnL"])])
    frac_neg.append(np.mean(deltas <= 0) if len(deltas) else np.nan)
    median_delta.append(np.median(deltas) if len(deltas) else np.nan)
ax.plot(ncw_values, frac_neg, "o-", lw=2, ms=8, color="firebrick",
        label="Fraction with Δln L ≤ 0")
ax.set_xlabel("N_CW")
ax.set_ylabel("Fraction of pulsars", color="firebrick")
ax.tick_params(axis="y", labelcolor="firebrick")
ax.set_ylim(0, 1.05)
ax.grid(True, alpha=0.3)
ax2 = ax.twinx()
ax2.plot(ncw_values, median_delta, "s--", lw=2, ms=7, color="steelblue",
         label="Median per-pulsar Δln L")
ax2.set_ylabel("Median per-pulsar Δln L", color="steelblue")
ax2.tick_params(axis="y", labelcolor="steelblue")
ax2.axhline(0, color="k", ls=":", lw=0.6)
ax.set_title("Is truth the per-pulsar max? (want frac=0, median > 0)")

fig.suptitle("Per-pulsar diagnostics for the joint Δln L metric", fontsize=13, y=1.00)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 6b. Per-pulsar breakdown at largest N_CW
"""),
    code(r"""ncw_show = ncw_values[-1]
joint_pp = all_results[ncw_show]["joint_per_pulsar"]

names = list(joint_pp.keys())
deltas = [joint_pp[n]["delta_lnL"] for n in names]

fig, ax = plt.subplots(figsize=(10, 5))
ax.bar(range(len(names)), deltas, color="steelblue", alpha=0.7)
ax.set_xticks(range(len(names)))
ax.set_xticklabels([n.split("-")[0] if "-" in n else n for n in names], rotation=45, ha="right")
ax.set_ylabel(r"$\Delta \ln L$ (truth $-$ best wrong within 3σ EM prior)")
ax.set_title(f"Per-pulsar distance constraint at N_CW = {ncw_show}")
ax.axhline(0, color="k", ls="--", lw=0.5)
ax.grid(True, alpha=0.2, axis="y")
fig.tight_layout()
plt.show()
"""),
    md(r"""## 7. Summary table
"""),
    code(r"""print("=" * 70)
print("REALISTIC CW FIRST-DETECTION SCAN")
print(f"  log10_h     ~ U{LOG10_H_RANGE}")
print(f"  log10_fgw   ~ U{LOG10_FGW_RANGE}")
print(f"  log10_mc    ~ U{LOG10_MC_RANGE}")
print(f"  GWB         log10_A = {GWB_LOG10_A}, gamma = {GWB_GAMMA:.3f}")
print(f"  white noise log10_equad = {LOG10_EQUAD}  ({10**LOG10_EQUAD*1e6:.2f} us)")
print(f"  Fourier components = {N_COMPONENTS}")
print("=" * 70)
print()
print(f"{'N_CW':>5s} | {'ln L(truth)':>13s} | {'ln L(wrong)':>13s} | {'joint Δln L':>12s}")
print("-" * 55)
for n in ncw_values:
    r = all_results[n]
    print(f"{n:5d} | {r['ll_truth']:13.2f} | {r['ll_best_wrong']:13.2f} | {r['delta_lnL_joint']:12.2f}")
"""),
    md(r"""## Interpretation

Compared to notebook 02 (fixed log10_h = -13.5, no realistic GWB, equad = -8):

- **GWB confusion**: with HD-correlated GWB at realistic amplitude, per-pulsar residuals
  now contain a common stochastic component that competes with the CW signal at the
  weakest strains. Sources near log10_h = -14.3 are at the GWB floor → near-zero
  per-source contribution to the distance constraint.
- **Source mix**: each realization has a mix of loud (log10_h ~ -13.9) and marginal
  (log10_h ~ -14.3) sources. Expect the loud sources to dominate the distance
  constraint; stacking marginal ones should help less than in the noiseless/clean case.
- **Realistic noise**: log10_equad = -6 reduces per-pulsar SNR compared to 10 ns. This is
  what actual PTAs operate at.
- **N_CW range**: first detections realistically feature 1-3 resolvable sources, not
  hundreds, so this scan is scientifically meaningful.
"""),
]

write_nb("05_realistic_first_detection.ipynb", cells)
print("Done!")
