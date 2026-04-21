"""Generate the noise-free CW distance landscape notebook."""
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


cells = [
    md(r"""# Noise-Free CW: Distance Likelihood Landscape

**Purpose:** Demonstrate that, in the absence of noise, the phase-connected CW likelihood
is essentially degenerate in pulsar distance for a single CW source, and that adding
multiple CW sources breaks this degeneracy.

By injecting pure CW signals (no white noise, no red noise, no GWB) and evaluating
the likelihood with white-noise-only noise model, we isolate the effect of pulsar
distance on the CW model fit without any noise-driven scatter.

**Key prediction:**
- Single CW: all distance modes have essentially the same likelihood (Δln L ≈ 0)
- Multiple CWs: only the true distance satisfies all CW phase constraints simultaneously
"""),
    code(r"""import sys, importlib
sys.path.insert(0, ".")
import cw_helpers
importlib.reload(cw_helpers)
from cw_helpers import *
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import warnings
warnings.filterwarnings("ignore")

%matplotlib inline
plt.rcParams.update({"figure.dpi": 120, "font.size": 11})
"""),
    md("## 1. Load pulsars"),
    code(r"""Npulsars = 10
_, disco_psrs = load_pulsars(Npulsars)
print(f"Loaded {Npulsars} pulsars (discovery only, no enterprise needed for noise-free injection)")
for p in disco_psrs:
    print(f"  {p.name:16s}  pdist = {p.pdist[0]:.4f} kpc")
"""),
    md(r"""## 2. Single CW: noise-free injection

Inject a single CW source and scan each pulsar's distance. We expect the likelihood
to be periodic in distance with nearly uniform peak heights — each mode corresponds
to a valid phase assignment for the pulsar term and is equally good.
"""),
    code(r"""cw_single = [{
    "cos_gwtheta": 0.3,
    "gwphi": 2.5,
    "cos_inc": -0.2,
    "log10_mc": 9.0,
    "log10_fgw": -8.0,
    "log10_h": -13.5,
    "phase0": 1.0,
    "psi": 0.7,
}]

rmap_single, true_dists = inject_noisefree_cw(disco_psrs, cw_single)
logl_fn_1, pkeys_1, bvals_1 = build_noisefree_likelihood(disco_psrs, rmap_single, 1, cw_single)

ll_truth_1 = float(logl_fn_1(bvals_1))
print(f"ln L(truth) = {ll_truth_1:.4f}")
print()

# Show CW signal RMS per pulsar
for psr in disco_psrs:
    rms = np.std(rmap_single[psr.name]) * 1e9
    dL = compute_mode_spacing(0.3, 2.5, -8.0, psr.pos)
    print(f"  {psr.name:16s}: CW rms = {rms:8.2f} ns, dL = {dL:.6f} kpc")
"""),
    code(r"""# Scan each pulsar's distance
single_cw_scans = {}

for psr in disco_psrs:
    dkey = f"{psr.name}_cw_p_dist"
    if dkey not in pkeys_1:
        continue

    true_d = psr.pdist[0]
    dL = compute_mode_spacing(0.3, 2.5, -8.0, psr.pos)

    # Scan over 10 mode spacings
    scan_half = 10 * dL
    scan_min = max(0.01, true_d - scan_half)
    scan_max = true_d + scan_half

    scan_d, scan_ll = scan_pulsar_distance(
        logl_fn_1, bvals_1, pkeys_1, dkey,
        scan_min, scan_max, n_points=4000,
    )

    single_cw_scans[psr.name] = {
        "scan_d": scan_d,
        "scan_ll": scan_ll,
        "true_dist": true_d,
        "dL": dL,
    }

print("Distance scans complete for all pulsars")
"""),
    md(r"""### Single CW distance landscape

Each panel shows $\Delta \ln L$ = ln L(d) - ln L(truth) for one pulsar.
The red dashed line marks the injected (true) distance.

If all peaks are at the same height (Δln L ≈ 0), the distance is uninformative
for a single CW — any mode is as good as any other.
"""),
    code(r"""ncols = 2
nrows = (len(single_cw_scans) + 1) // 2
fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3 * nrows))
axes = axes.flatten()

for i, (pname, sc) in enumerate(sorted(single_cw_scans.items())):
    ax = axes[i]
    rel_ll = sc["scan_ll"] - ll_truth_1

    ax.plot(sc["scan_d"], rel_ll, "b-", lw=0.5, alpha=0.8)
    ax.axvline(sc["true_dist"], color="red", ls="--", lw=0.8, alpha=0.6, label="injected dist")

    # Find and mark peaks
    prom = max(0.01, 0.01 * (np.max(rel_ll) - np.min(rel_ll)))
    peaks, _ = find_peaks(sc["scan_ll"], prominence=prom, distance=3)
    if len(peaks) > 0:
        peak_heights = rel_ll[peaks]
        ax.plot(sc["scan_d"][peaks], peak_heights, "o", color="orange", ms=2, alpha=0.6)

        # Annotate peak height range
        ax.set_title(
            f"{pname}  |  dL={sc['dL']:.5f} kpc  |  "
            f"peak range: [{np.min(peak_heights):.4f}, {np.max(peak_heights):.4f}]",
            fontsize=9,
        )
    else:
        ax.set_title(f"{pname}  |  dL={sc['dL']:.5f} kpc", fontsize=9)

    ax.set_xlabel("p_dist [kpc]")
    ax.set_ylabel(r"$\Delta \ln L$")
    ax.grid(True, alpha=0.2)
    if i == 0:
        ax.legend(fontsize=8)

for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle("Single CW (noise-free): Distance landscape — all modes equally good?",
             fontsize=13, y=1.01)
fig.tight_layout()
plt.show()
"""),
    md(r"""### Peak height statistics (single CW)

For each pulsar, collect all peak heights and show their range.
If the range is ≈ 0, the distance is uninformative.
"""),
    code(r"""print(f"{'Pulsar':16s} | {'N peaks':>8s} | {'Peak max':>10s} | {'Peak min':>10s} | {'Range':>10s} | {'dL (kpc)':>10s}")
print("-" * 75)

for pname in sorted(single_cw_scans.keys()):
    sc = single_cw_scans[pname]
    rel_ll = sc["scan_ll"] - ll_truth_1
    prom = max(0.01, 0.01 * (np.max(rel_ll) - np.min(rel_ll)))
    peaks, _ = find_peaks(sc["scan_ll"], prominence=prom, distance=3)

    if len(peaks) > 0:
        ph = rel_ll[peaks]
        print(f"{pname:16s} | {len(peaks):8d} | {np.max(ph):10.4f} | {np.min(ph):10.4f} | "
              f"{np.max(ph)-np.min(ph):10.4f} | {sc['dL']:10.6f}")
    else:
        print(f"{pname:16s} | {'0':>8s} | {'N/A':>10s} | {'N/A':>10s} | {'N/A':>10s} | {sc['dL']:10.6f}")
"""),
    md(r"""## 3. Multiple CWs: noise-free injection

Now inject multiple CW sources (N=2, 4, 7, 10) with a single shared distance per pulsar.
Each CW creates its own periodic distance modulation with a different spacing $\Delta L$.
The true distance is the only point where ALL CW phases are simultaneously correct.
"""),
    code(r"""rng = np.random.default_rng(42)

def make_cw_params(n_sources, rng, scenario="well_separated"):
    # Generate n_sources CW parameter dicts.
    sources = []
    for i in range(n_sources):
        if scenario == "close_freq":
            log10_fgw = -8.0 + rng.uniform(-0.1, 0.1)
            cos_gwtheta = rng.uniform(-1, 1)
            gwphi = rng.uniform(0, 2 * np.pi)
        elif scenario == "close_sky":
            log10_fgw = rng.uniform(-8.5, -7.5)
            cos_gwtheta = np.clip(0.3 + rng.uniform(-0.05, 0.05), -1, 1)
            gwphi = 2.5 + rng.uniform(-0.055, 0.055)
        else:  # well_separated
            log10_fgw = rng.uniform(-8.5, -7.5)
            cos_gwtheta = rng.uniform(-1, 1)
            gwphi = rng.uniform(0, 2 * np.pi)

        sources.append({
            "cos_gwtheta": cos_gwtheta,
            "gwphi": gwphi,
            "cos_inc": rng.uniform(-1, 1),
            "log10_mc": 9.0,
            "log10_fgw": log10_fgw,
            "log10_h": -13.5,
            "phase0": rng.uniform(0, 2 * np.pi),
            "psi": rng.uniform(0, np.pi),
        })
    return sources
"""),
    code(r"""ncw_values = [1, 2, 4, 7, 10, 25, 50, 100, 250, 500]
test_psr = disco_psrs[0]  # Use first pulsar as representative
test_psr_name = test_psr.name
scenario = "well_separated"

multi_cw_results = {}

for ncw in ncw_values:
    print(f"\n--- N_CW = {ncw} ---")

    rng_local = np.random.default_rng(42)
    cw_params = make_cw_params(ncw, rng_local, scenario=scenario)

    # Inject noise-free
    rmap, _ = inject_noisefree_cw(disco_psrs, cw_params)
    logl_fn, pkeys, bvals = build_noisefree_likelihood(disco_psrs, rmap, ncw, cw_params)

    ll_truth = float(logl_fn(bvals))
    print(f"  ln L(truth) = {ll_truth:.2f}")

    # Compute mode spacings for all CWs at test pulsar
    all_dL = []
    for cw_p in cw_params:
        dL = compute_mode_spacing(cw_p["cos_gwtheta"], cw_p["gwphi"],
                                  cw_p["log10_fgw"], test_psr.pos)
        all_dL.append(dL)
        print(f"  CW: f_gw=10^{cw_p['log10_fgw']:.2f} Hz, dL={dL:.6f} kpc")

    min_dL = min(all_dL)
    max_dL = max(all_dL)

    # Scan test pulsar distance
    dkey = f"{test_psr_name}_cw_p_dist"
    true_d = test_psr.pdist[0]
    scan_half = 10 * max_dL
    scan_min = max(0.01, true_d - scan_half)
    scan_max = true_d + scan_half

    scan_d, scan_ll = scan_pulsar_distance(
        logl_fn, bvals, pkeys, dkey,
        scan_min, scan_max, n_points=5000,
    )

    # Peak analysis
    analysis = analyze_peaks(scan_d, scan_ll, true_d, mode_spacing=min_dL)

    multi_cw_results[ncw] = {
        "scan_d": scan_d,
        "scan_ll": scan_ll,
        "true_dist": true_d,
        "all_dL": all_dL,
        "min_dL": min_dL,
        "max_dL": max_dL,
        "ll_truth": ll_truth,
        "analysis": analysis,
    }

    delta = analysis["delta_lnL_truth_vs_best_wrong"]
    print(f"  Peaks found: {analysis['n_peaks']}")
    print(f"  Delta ln L (truth - best wrong): {delta:.4f}" if not np.isnan(delta) else
          f"  Delta ln L: N/A")

    # Joint: best wrong mode within EM prior for each pulsar
    joint = compute_joint_best_wrong_in_prior(
        logl_fn, bvals, pkeys,
        disco_psrs, disco_psrs, cw_params, n_sigma=3.0,
    )
    multi_cw_results[ncw]["ll_wrong_joint"] = joint["ll_best_wrong"]
    multi_cw_results[ncw]["delta_joint"] = joint["delta_lnL"]
    print(f"  Joint Delta ln L (best wrong in EM prior): {joint['delta_lnL']:.4f}")

    jax.clear_caches()

print("\nDone!")
"""),
    md(r"""### Distance landscape comparison

As N_CW increases, the truth peak should progressively dominate the landscape.
The wrong-mode peaks should drop in height because the single distance cannot
simultaneously satisfy all CW phase constraints.
"""),
    code(r"""fig, axes = plt.subplots(len(ncw_values), 1, figsize=(14, 3.5 * len(ncw_values)))

for ax_idx, ncw in enumerate(ncw_values):
    ax = axes[ax_idx]
    r = multi_cw_results[ncw]
    scan_d, scan_ll = r["scan_d"], r["scan_ll"]
    baseline = r["ll_truth"]
    rel_ll = scan_ll - baseline

    ax.plot(scan_d, rel_ll, "b-", lw=0.4, alpha=0.8)
    ax.axvline(r["true_dist"], color="red", ls="--", lw=1, alpha=0.6, label="injected dist")

    # Mark peaks
    a = r["analysis"]
    if len(a["peak_distances"]) > 0:
        truth_mask = np.abs(a["peak_distances"] - r["true_dist"]) < 0.3 * r["min_dL"]
        wrong_mask = ~truth_mask

        if np.any(wrong_mask):
            wrong_rel = a["peak_heights"][wrong_mask] - baseline
            ax.plot(a["peak_distances"][wrong_mask], wrong_rel,
                    "o", color="orange", ms=2, alpha=0.5, label="wrong peaks")
        if np.any(truth_mask):
            truth_rel = a["peak_heights"][truth_mask] - baseline
            ax.plot(a["peak_distances"][truth_mask], truth_rel,
                    "o", color="green", ms=5, label="truth peak")

    delta = a["delta_lnL_truth_vs_best_wrong"]
    delta_str = f"{delta:.4f}" if not np.isnan(delta) else "N/A"
    ax.set_ylabel(r"$\Delta \ln L$")
    ax.set_title(f"N_CW = {ncw}  |  "
                 r"$\Delta \ln L$(truth $-$ best wrong) = " + delta_str +
                 f"  |  joint shift: {r['delta_joint']:.4f}",
                 fontsize=10)
    ax.grid(True, alpha=0.2)
    if ax_idx == 0:
        ax.legend(fontsize=8, loc="lower right")

axes[-1].set_xlabel(f"p_dist [kpc] — {test_psr_name}")
fig.suptitle(f"Noise-free CW: {test_psr_name} distance landscape vs N_CW ({scenario})",
             fontsize=13, y=1.01)
fig.tight_layout()
plt.show()
"""),
    md("### Summary plot: Δln L vs N_CW"),
    code(r"""fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Per-pulsar best wrong peak
ax = axes[0]
per_psr_delta = [multi_cw_results[n]["analysis"]["delta_lnL_truth_vs_best_wrong"]
                 for n in ncw_values]
ax.plot(ncw_values, per_psr_delta, "o-", color="steelblue", lw=2, ms=8)
ax.set_xlabel("Number of CW sources")
ax.set_ylabel(r"$\Delta \ln L$ (truth $-$ best wrong peak)")
ax.set_title(f"Per-pulsar ({test_psr_name}): truth peak advantage")
ax.grid(True, alpha=0.3)

# Joint shift
ax = axes[1]
joint_delta = [multi_cw_results[n]["delta_joint"] for n in ncw_values]
ax.plot(ncw_values, joint_delta, "s-", color="firebrick", lw=2, ms=8)
ax.set_xlabel("Number of CW sources")
ax.set_ylabel(r"$\Delta \ln L$ (truth $-$ best wrong in EM prior)")
ax.set_title("Joint: best wrong mode within EM prior")
ax.grid(True, alpha=0.3)

fig.suptitle(f"Noise-free: distance importance vs N_CW ({scenario})", fontsize=13, y=1.02)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 4. Scenario comparison (noise-free)

Repeat the multi-CW analysis for all three scenarios: close in frequency,
close in sky position, and well-separated.
"""),
    code(r"""all_scenarios = ["close_freq", "close_sky", "well_separated"]
scenario_labels = {
    "close_freq": "Close in frequency",
    "close_sky": "Close in sky position",
    "well_separated": "Well-separated",
}

scenario_results = {}

for scenario in all_scenarios:
    print(f"\n{'='*50}")
    print(f"Scenario: {scenario_labels[scenario]}")
    print(f"{'='*50}")

    for ncw in ncw_values:
        rng_local = np.random.default_rng(42)
        cw_params = make_cw_params(ncw, rng_local, scenario=scenario)

        rmap, _ = inject_noisefree_cw(disco_psrs, cw_params)
        logl_fn, pkeys, bvals = build_noisefree_likelihood(disco_psrs, rmap, ncw, cw_params)

        ll_truth = float(logl_fn(bvals))

        # Per-pulsar scan
        dkey = f"{test_psr_name}_cw_p_dist"
        true_d = test_psr.pdist[0]
        all_dL = [compute_mode_spacing(cp["cos_gwtheta"], cp["gwphi"],
                                       cp["log10_fgw"], test_psr.pos) for cp in cw_params]
        min_dL = min(all_dL)
        max_dL = max(all_dL)
        scan_half = 10 * max_dL
        scan_d, scan_ll = scan_pulsar_distance(
            logl_fn, bvals, pkeys, dkey,
            max(0.01, true_d - scan_half), true_d + scan_half, n_points=5000,
        )
        analysis = analyze_peaks(scan_d, scan_ll, true_d, mode_spacing=min_dL)

        # Joint: best wrong mode within EM prior
        joint = compute_joint_best_wrong_in_prior(
            logl_fn, bvals, pkeys,
            disco_psrs, disco_psrs, cw_params, n_sigma=3.0,
        )

        key = (scenario, ncw)
        scenario_results[key] = {
            "per_psr_delta": analysis["delta_lnL_truth_vs_best_wrong"],
            "joint_delta": joint["delta_lnL"],
            "n_peaks": analysis["n_peaks"],
        }

        print(f"  N_CW={ncw:2d}: per-psr Δ={analysis['delta_lnL_truth_vs_best_wrong']:.4f}, "
              f"joint Δ={joint['delta_lnL']:.4f}")

        jax.clear_caches()

print("\nAll scenarios complete!")
"""),
    code(r"""fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Per-pulsar
ax = axes[0]
for scenario in all_scenarios:
    deltas = [scenario_results[(scenario, n)]["per_psr_delta"] for n in ncw_values]
    ax.plot(ncw_values, deltas, "o-", label=scenario_labels[scenario], lw=2, ms=8)
ax.set_xlabel("Number of CW sources")
ax.set_ylabel(r"$\Delta \ln L$ (truth $-$ best wrong peak)")
ax.set_title(f"Per-pulsar ({test_psr_name})")
ax.legend()
ax.grid(True, alpha=0.3)

# Joint
ax = axes[1]
for scenario in all_scenarios:
    deltas = [scenario_results[(scenario, n)]["joint_delta"] for n in ncw_values]
    ax.plot(ncw_values, deltas, "s-", label=scenario_labels[scenario], lw=2, ms=8)
ax.set_xlabel("Number of CW sources")
ax.set_ylabel(r"$\Delta \ln L$ (truth $-$ best wrong in EM prior)")
ax.set_title("Joint: best wrong in EM prior")
ax.legend()
ax.grid(True, alpha=0.3)

fig.suptitle("Noise-free: Δln L vs N_CW by scenario", fontsize=13, y=1.02)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 5. Summary
"""),
    code(r"""print("=" * 80)
print("NOISE-FREE SUMMARY")
print("=" * 80)

print(f"\n{'Scenario':25s} | " + " | ".join(f"N_CW={n}" for n in ncw_values))
print("-" * (25 + 3 + len(ncw_values) * 14))

print("\nPer-pulsar Delta ln L (truth peak - best wrong peak):")
for scenario in all_scenarios:
    row = f"{scenario_labels[scenario]:25s} | "
    for ncw in ncw_values:
        d = scenario_results[(scenario, ncw)]["per_psr_delta"]
        row += f" {d:10.4f}   | "
    print(row)

print(f"\nJoint Delta ln L (truth vs best wrong mode within EM prior):")
for scenario in all_scenarios:
    row = f"{scenario_labels[scenario]:25s} | "
    for ncw in ncw_values:
        d = scenario_results[(scenario, ncw)]["joint_delta"]
        row += f" {d:10.4f}   | "
    print(row)
"""),
    md(r"""## Interpretation

**Noise-free results isolate the pure geometric/phase effect:**

1. **Single CW (N_CW=1):** The truth peak and all wrong-mode peaks should have
   essentially the same likelihood. This confirms that for a single source,
   the pulsar distance is uninformative — each mode is a valid phase assignment.

2. **Multiple CWs:** The truth peak should progressively dominate as N_CW grows.
   Each additional CW source imposes an additional phase constraint that only
   the true distance can satisfy simultaneously.

3. **Scenario comparison:** If the effect is similar across scenarios (close freq,
   close sky, well-separated), the distance constraint is a generic multi-source
   effect, not specific to "confused" sources.

The noise-free case provides the cleanest test because:
- No noise scatter to confuse peak heights
- No red noise GP degrees of freedom to absorb signal
- The only thing that matters is the CW phase matching
"""),
]

with open("00_noisefree_distance_landscape.ipynb", "w") as f:
    json.dump(notebook(cells), f, indent=1)
print("Wrote 00_noisefree_distance_landscape.ipynb")
