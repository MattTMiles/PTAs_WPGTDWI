"""Generate the analysis notebooks as .ipynb files."""
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


# ====================================================================
# Notebook 1: Single CW distance landscape
# ====================================================================

nb1_cells = [
    md(r"""# Single CW: Pulsar Distance Likelihood Landscape

**Question:** For a single continuous wave source with correctly identified global parameters
and noise, how much does the log-likelihood change when pulsar distances are correct vs wrong?

We inject a single CW at three signal strengths (`log10_h = -14.0, -13.5, -13.0`),
fix all global CW + noise parameters at truth, and scan each pulsar's distance to map
the oscillatory likelihood landscape. The key metric is the height difference between
the truth-distance peak and the best alternative (wrong-mode) peak.

**Expectation:** For a single CW, the likelihood is periodic in each pulsar's distance
with spacing $\Delta L = 1 / (f_{\rm gw} \cdot \kappa_{\rm kpc/c} \cdot |1 - \cos\mu|)$.
Alternative mode peaks should be nearly as high as the truth peak, because each mode
corresponds to a valid phase assignment for the pulsar term.
"""),
    code(r"""import sys, importlib
sys.path.insert(0, ".")
import cw_helpers
importlib.reload(cw_helpers)
from cw_helpers import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

%matplotlib inline
plt.rcParams.update({"figure.dpi": 120, "font.size": 11})
"""),
    md("## 1. Load pulsars"),
    code(r"""Npulsars = 10
ent_psrs, disco_psrs = load_pulsars(Npulsars)
print(f"Loaded {len(ent_psrs)} pulsars:")
for p in ent_psrs:
    print(f"  {p.name:16s}  pdist = {p.pdist[0]:.3f} +/- {p.pdist[1]:.3f} kpc   Ntoas = {len(p.toas)}")
"""),
    md("## 2. Build enterprise PTA (1 CW source) for simulation"),
    code(r"""num_cw = 1
pta, cw_block_names, Tspan = build_enterprise_pta(ent_psrs, num_cw)
print(f"PTA built: {len(pta.param_names)} parameters, Tspan = {Tspan/3.156e7:.1f} yr")
"""),
    md(r"""## 3. Scan distance landscape at three signal strengths

For each `log10_h`, we:
1. Generate injection parameters and simulate data
2. Build the discovery likelihood (fixing all params at truth)
3. Scan each pulsar's distance over ~10 mode spacings
4. Find peaks and compare truth peak to best wrong-mode peak
"""),
    code(r"""log10_h_values = [-14.0, -13.5, -13.0]
seed = 1234

# Storage for all results
all_results = {}

for log10_h in log10_h_values:
    print(f"\n{'='*60}")
    print(f"log10_h = {log10_h}")
    print(f"{'='*60}")

    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    # Generate injection
    inj_params = generate_injection_params(
        pta, ent_psrs, num_cw, cw_block_names, log10_h, scenario="single", rng=rng,
    )

    # Simulate
    sim_resids = simulate(pta, inj_params)
    resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}

    # Build likelihood
    logl_fn, param_keys, base_vals = build_disco_likelihood(
        disco_psrs, resid_map, num_cw, inj_params, cw_block_names,
    )

    # Warm up JIT
    _ = logl_fn(base_vals)
    print(f"  ln L(truth) = {float(logl_fn(base_vals)):.2f}")

    # Get CW params for mode spacing calculation
    cw_cos_gwtheta = inj_params["cw_cos_gwtheta"]
    cw_gwphi = inj_params["cw_gwphi"]
    cw_log10_fgw = inj_params["cw_log10_fgw"]

    results = {}
    for psr_e, psr_d in zip(ent_psrs, disco_psrs):
        pname = psr_e.name
        true_dist = psr_e.pdist[0]
        dist_key = f"{pname}_cw_p_dist"

        if dist_key not in param_keys:
            continue

        # Mode spacing
        dL = compute_mode_spacing(cw_cos_gwtheta, cw_gwphi, cw_log10_fgw, psr_d.pos)

        # Scan range: truth +/- 5 mode spacings, but stay positive
        n_modes_scan = 10
        scan_half = n_modes_scan * dL
        scan_min = max(0.01, true_dist - scan_half)
        scan_max = true_dist + scan_half

        print(f"  {pname}: dL = {dL:.4f} kpc, scanning [{scan_min:.3f}, {scan_max:.3f}]")

        scan_d, scan_ll = scan_pulsar_distance(
            logl_fn, base_vals, param_keys, dist_key,
            scan_min, scan_max, n_points=3000,
        )

        analysis = analyze_peaks(scan_d, scan_ll, true_dist, mode_spacing=dL)
        analysis["scan_d"] = scan_d
        analysis["scan_ll"] = scan_ll
        analysis["true_dist"] = true_dist
        analysis["dL"] = dL
        results[pname] = analysis

        delta = analysis["delta_lnL_truth_vs_best_wrong"]
        print(f"    -> {analysis['n_peaks']} peaks, "
              f"Delta ln L (truth - best wrong) = {delta:.2f}" if not np.isnan(delta) else
              f"    -> {analysis['n_peaks']} peaks, no wrong peaks found")

    all_results[log10_h] = results
    jax.clear_caches()

print("\nDone!")
"""),
    md(r"""## 4. Distance landscape plots

Each panel shows the log-likelihood as a function of one pulsar's distance (all other
params fixed at truth). The vertical red dashed line marks the injected distance.
Peaks are marked with dots. The truth peak is green; wrong-mode peaks are orange.
"""),
    code(r"""for log10_h in log10_h_values:
    results = all_results[log10_h]
    psr_names = list(results.keys())
    n = len(psr_names)
    ncols = 2
    nrows = (n + 1) // 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.5 * nrows))
    axes = axes.flatten()

    for i, pname in enumerate(psr_names):
        ax = axes[i]
        r = results[pname]
        scan_d, scan_ll = r["scan_d"], r["scan_ll"]

        # Normalize to truth peak
        baseline = r["truth_logl"]
        ax.plot(scan_d, scan_ll - baseline, "b-", lw=0.5, alpha=0.8)

        # Mark peaks
        if len(r["peak_distances"]) > 0:
            truth_mask = np.abs(r["peak_distances"] - r["true_dist"]) < 0.3 * r["dL"]
            wrong_mask = ~truth_mask

            if np.any(wrong_mask):
                ax.plot(
                    r["peak_distances"][wrong_mask],
                    r["peak_heights"][wrong_mask] - baseline,
                    "o", color="orange", ms=3, alpha=0.7, label="wrong-mode peaks",
                )
            if np.any(truth_mask):
                ax.plot(
                    r["peak_distances"][truth_mask],
                    r["peak_heights"][truth_mask] - baseline,
                    "o", color="green", ms=5, label="truth peak",
                )

        ax.axvline(r["true_dist"], color="red", ls="--", lw=0.8, alpha=0.6, label="injected")
        ax.set_xlabel("p_dist [kpc]")
        ax.set_ylabel(r"$\Delta \ln L$")
        ax.set_title(f"{pname}  (dL={r['dL']:.4f} kpc)")
        if i == 0:
            ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(f"Single CW: Distance landscape  |  log10_h = {log10_h}", fontsize=14, y=1.01)
    fig.tight_layout()
    plt.show()
"""),
    md(r"""## 5. Peak height comparison: truth vs wrong-mode peaks

For each pulsar, we show the distribution of wrong-mode peak heights relative
to the truth peak. If Δln L ≈ 0 for most wrong-mode peaks, the distance
doesn't matter much for a single CW.
"""),
    code(r"""fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

for ax_idx, log10_h in enumerate(log10_h_values):
    ax = axes[ax_idx]
    results = all_results[log10_h]

    psr_names_sorted = sorted(results.keys())
    for j, pname in enumerate(psr_names_sorted):
        r = results[pname]
        if r["n_peaks"] < 2:
            continue

        # Δln L for each wrong-mode peak relative to truth
        truth_mask = np.abs(r["peak_distances"] - r["true_dist"]) < 0.3 * r["dL"]
        wrong_heights = r["peak_heights"][~truth_mask]
        if len(wrong_heights) == 0:
            continue

        delta = r["truth_logl"] - wrong_heights
        ax.scatter(
            np.full(len(delta), j), delta,
            s=8, alpha=0.5, color="steelblue",
        )

    ax.set_xticks(range(len(psr_names_sorted)))
    ax.set_xticklabels([n.split("-")[0] if "-" in n else n for n in psr_names_sorted],
                       rotation=45, ha="right", fontsize=8)
    ax.set_xlabel("Pulsar")
    ax.set_title(f"log10_h = {log10_h}")
    ax.axhline(0, color="k", ls="--", lw=0.5)
    ax.grid(True, alpha=0.2)

axes[0].set_ylabel(r"$\Delta \ln L$ (truth peak $-$ wrong peak)")
fig.suptitle("Single CW: Peak height differences (truth vs wrong-mode)", fontsize=13, y=1.02)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 6. Summary statistics

Table of Δln L between truth peak and best wrong-mode peak, per pulsar and signal strength.
"""),
    code(r"""print(f"{'Pulsar':16s} | " + " | ".join(f"log10_h={h:5.1f}" for h in log10_h_values))
print("-" * (16 + 3 + len(log10_h_values) * 18))

psr_names = sorted(all_results[log10_h_values[0]].keys())
for pname in psr_names:
    row = f"{pname:16s} | "
    for log10_h in log10_h_values:
        r = all_results[log10_h].get(pname, {})
        delta = r.get("delta_lnL_truth_vs_best_wrong", np.nan)
        n_peaks = r.get("n_peaks", 0)
        row += f"  {delta:7.2f} ({n_peaks:3d}p) | "
    print(row)

print()
print("(Np) = number of peaks found in scan range")
print("Delta ln L > 0 means truth peak is higher than best wrong-mode peak")
"""),
    md(r"""## 7. Joint distance effect: best wrong mode within EM prior

For each pulsar, find the best non-truth peak **within the EM distance prior**
(3-sigma range). Then evaluate the joint likelihood with ALL pulsars simultaneously
at their respective best wrong modes. This answers: "given realistic EM distance
constraints, how much worse is the best alternative distance assignment?"
"""),
    code(r"""print(f"{'Signal':12s} | {'ln L (truth)':>14s} | {'ln L (wrong)':>14s} | {'Delta ln L':>12s} | {'Per-pulsar details'}")
print("-" * 90)

for log10_h in log10_h_values:
    # Re-build the likelihood for this signal strength
    rng = np.random.default_rng(seed)
    np.random.seed(seed)
    inj_params = generate_injection_params(
        pta, ent_psrs, num_cw, cw_block_names, log10_h, scenario="single", rng=rng,
    )
    sim_resids = simulate(pta, inj_params)
    resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}
    logl_fn, param_keys, base_vals = build_disco_likelihood(
        disco_psrs, resid_map, num_cw, inj_params, cw_block_names,
    )

    joint = compute_joint_best_wrong_in_prior(
        logl_fn, base_vals, param_keys,
        ent_psrs, disco_psrs, inj_params,
        cw_block_names=cw_block_names, n_sigma=3.0,
    )

    print(f"h={log10_h:5.1f}     | {joint['ll_truth']:14.2f} | {joint['ll_best_wrong']:14.2f} | {joint['delta_lnL']:12.2f} |")
    for pname, r in joint["per_pulsar"].items():
        sig = [pe.pdist[1] for pe in ent_psrs if pe.name == pname][0]
        print(f"  {pname:16s}: Delta={r['delta_lnL']:8.2f}, peaks_in_prior={r['n_peaks_in_prior']:4d}, "
              f"EM prior width={sig:.4f} kpc")
    jax.clear_caches()
"""),
    md(r"""## Interpretation

- If Δln L between truth and best wrong-mode peak is **small** (≲ a few) for individual
  pulsars, then knowing the correct distance provides only a marginal likelihood improvement
  for a single CW.
- If the joint Δln L (all pulsars at their best wrong mode within EM prior) is small,
  the contention that "correct distances cause a dramatic likelihood jump" is not supported
  for the single-source case.

Compare these results with the multi-CW analysis in `02_multi_cw_distance_landscape.ipynb`.
"""),
]

# ====================================================================
# Notebook 2: Multi-CW distance landscape
# ====================================================================

nb2_cells = [
    md(r"""# Multi-CW: Pulsar Distance Likelihood Landscape

**Question:** How does the likelihood penalty for wrong pulsar distances change
when multiple CW sources share the same pulsar distance parameter?

We test N_CW = {1, 2, 4, 7, 10} under three scenarios:
- **Close in frequency:** all CWs within ±0.1 dex of $f_{\rm gw} = 10^{-8}$ Hz
- **Close in sky position:** all CWs within ~15° on sky
- **Well-separated:** random sky positions and frequencies

The physical argument: with multiple CWs, the single distance per pulsar must
simultaneously be at a peak of ALL CW pulsar-term oscillations. Since different
CWs produce oscillations with (generically) incommensurate spacings, only the
true distance satisfies all constraints simultaneously.
"""),
    code(r"""import sys, importlib
sys.path.insert(0, ".")
import cw_helpers
importlib.reload(cw_helpers)
from cw_helpers import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

%matplotlib inline
plt.rcParams.update({"figure.dpi": 120, "font.size": 11})
"""),
    md("## 1. Load pulsars"),
    code(r"""Npulsars = 10
ent_psrs, disco_psrs = load_pulsars(Npulsars)
print(f"Loaded {Npulsars} pulsars")
"""),
    md(r"""## 2. Configuration

We fix `log10_h = -13.5` (moderate-strong signal) and vary N_CW and scenario.
For each configuration we scan a representative pulsar's distance and compare
the truth peak to wrong-mode peaks.
"""),
    code(r"""log10_h = -13.5
seed = 1234
ncw_values = [1, 2, 4, 7, 10, 25, 50, 100, 250, 500]
scenarios = ["close_freq", "close_sky", "well_separated"]
scenario_labels = {
    "close_freq": "Close in frequency",
    "close_sky": "Close in sky position",
    "well_separated": "Well-separated",
}
"""),
    md(r"""## 3. Run distance scans

For each (N_CW, scenario) combination:
1. Build enterprise PTA with N_CW sources
2. Generate injection parameters
3. Simulate data
4. Build discovery likelihood
5. Scan first 3 pulsars' distances and analyze peaks

This gives us the oscillatory structure and peak-height comparisons.
"""),
    code(r"""all_results = {}

for scenario in scenarios:
    print(f"\n{'#'*70}")
    print(f"# Scenario: {scenario_labels[scenario]}")
    print(f"{'#'*70}")

    for ncw in ncw_values:
        key = (scenario, ncw)
        print(f"\n  N_CW = {ncw}")

        rng = np.random.default_rng(seed)
        np.random.seed(seed)

        # Build PTA
        pta, cw_block_names, Tspan = build_enterprise_pta(ent_psrs, ncw)

        # Generate injection
        inj_params = generate_injection_params(
            pta, ent_psrs, ncw, cw_block_names, log10_h, scenario=scenario, rng=rng,
        )

        # Simulate
        sim_resids = simulate(pta, inj_params)
        resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}

        # Build likelihood
        logl_fn, param_keys, base_vals = build_disco_likelihood(
            disco_psrs, resid_map, ncw, inj_params, cw_block_names,
        )

        _ = logl_fn(base_vals)
        ll_truth = float(logl_fn(base_vals))
        print(f"    ln L(truth) = {ll_truth:.2f}")

        # Scan first 3 pulsars
        results_per_psr = {}
        disco_suffixes = ["" if i == 0 else f"_{i+1}" for i in range(ncw)]

        for psr_e, psr_d in zip(ent_psrs[:3], disco_psrs[:3]):
            pname = psr_e.name
            true_dist = psr_e.pdist[0]
            dist_key = f"{pname}_cw_p_dist"

            if dist_key not in param_keys:
                continue

            # Mode spacing for the first CW (representative)
            cw1_cos_gwtheta = inj_params["cw_cos_gwtheta"]
            cw1_gwphi = inj_params["cw_gwphi"]
            cw1_log10_fgw = inj_params["cw_log10_fgw"]
            dL = compute_mode_spacing(cw1_cos_gwtheta, cw1_gwphi, cw1_log10_fgw, psr_d.pos)

            # Also compute mode spacings for ALL CWs
            all_dL = []
            for suffix in disco_suffixes:
                cgwt_key = f"cw_cos_gwtheta{suffix}"
                gwphi_key = f"cw_gwphi{suffix}"
                fgw_key = f"cw_log10_fgw{suffix}"
                # Map back to enterprise names
                ent_name = cw_block_names[disco_suffixes.index(suffix)]
                _dL = compute_mode_spacing(
                    inj_params[f"{ent_name}_cos_gwtheta"],
                    inj_params[f"{ent_name}_gwphi"],
                    inj_params[f"{ent_name}_log10_fgw"],
                    psr_d.pos,
                )
                all_dL.append(_dL)

            min_dL = min(all_dL)

            # Scan: use the smallest mode spacing to set resolution
            n_modes_scan = 10
            scan_half = n_modes_scan * max(all_dL)
            scan_min = max(0.01, true_dist - scan_half)
            scan_max = true_dist + scan_half

            print(f"    {pname}: dL range = [{min(all_dL):.4f}, {max(all_dL):.4f}] kpc")

            scan_d, scan_ll = scan_pulsar_distance(
                logl_fn, base_vals, param_keys, dist_key,
                scan_min, scan_max, n_points=4000,
            )

            analysis = analyze_peaks(scan_d, scan_ll, true_dist, mode_spacing=min_dL)
            analysis["scan_d"] = scan_d
            analysis["scan_ll"] = scan_ll
            analysis["true_dist"] = true_dist
            analysis["all_dL"] = all_dL
            analysis["min_dL"] = min_dL
            results_per_psr[pname] = analysis

            delta = analysis["delta_lnL_truth_vs_best_wrong"]
            print(f"      -> {analysis['n_peaks']} peaks, "
                  f"Delta ln L = {delta:.2f}" if not np.isnan(delta) else
                  f"      -> {analysis['n_peaks']} peaks, no wrong peaks")

        # Joint wrong-distance: best wrong mode within EM prior for each pulsar
        joint = compute_joint_best_wrong_in_prior(
            logl_fn, base_vals, param_keys,
            ent_psrs, disco_psrs, inj_params,
            cw_block_names=cw_block_names, n_sigma=3.0,
        )

        all_results[key] = {
            "psr_results": results_per_psr,
            "ll_truth": ll_truth,
            "ll_best_wrong": joint["ll_best_wrong"],
            "delta_lnL_joint": joint["delta_lnL"],
            "joint_per_pulsar": joint["per_pulsar"],
        }
        print(f"    Joint Delta ln L (best wrong in EM prior): {joint['delta_lnL']:.2f}")

        jax.clear_caches()

print("\nAll scans complete!")
"""),
    md(r"""## 4. Distance landscape plots: comparing N_CW

For the first pulsar, show how the landscape changes as N_CW increases.
"""),
    code(r"""first_psr = ent_psrs[0].name

for scenario in scenarios:
    fig, axes = plt.subplots(len(ncw_values), 1, figsize=(14, 3 * len(ncw_values)), sharex=False)

    for ax_idx, ncw in enumerate(ncw_values):
        ax = axes[ax_idx]
        key = (scenario, ncw)
        r = all_results[key]["psr_results"].get(first_psr)

        if r is None:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center")
            continue

        scan_d, scan_ll = r["scan_d"], r["scan_ll"]
        baseline = r["truth_logl"]

        ax.plot(scan_d, scan_ll - baseline, "b-", lw=0.5, alpha=0.8)

        # Mark peaks
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
        ax.set_title(f"N_CW = {ncw}  |  peaks = {r['n_peaks']}  |  "
                     r"$\Delta \ln L$(truth-best wrong) = " + delta_str, fontsize=10)
        ax.grid(True, alpha=0.2)

    axes[-1].set_xlabel("p_dist [kpc]")
    fig.suptitle(f"{scenario_labels[scenario]}: {first_psr} distance landscape  |  "
                 f"log10_h = {log10_h}", fontsize=13, y=1.01)
    fig.tight_layout()
    plt.show()
"""),
    md(r"""## 5. Summary: Δln L vs N_CW for each scenario

The key plot: how does the truth-vs-wrong-mode Δln L scale with the number of CW sources?
"""),
    code(r"""fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: Per-pulsar best wrong-mode Δln L (averaged over scanned pulsars)
ax = axes[0]
for scenario in scenarios:
    avg_deltas = []
    for ncw in ncw_values:
        key = (scenario, ncw)
        deltas = [
            r["delta_lnL_truth_vs_best_wrong"]
            for r in all_results[key]["psr_results"].values()
            if not np.isnan(r["delta_lnL_truth_vs_best_wrong"])
        ]
        avg_deltas.append(np.mean(deltas) if deltas else np.nan)
    ax.plot(ncw_values, avg_deltas, "o-", label=scenario_labels[scenario], lw=2, ms=8)

ax.set_xlabel("Number of CW sources")
ax.set_ylabel(r"Mean $\Delta \ln L$ (truth $-$ best wrong peak)")
ax.set_title("Per-pulsar: truth peak advantage")
ax.legend()
ax.grid(True, alpha=0.3)

# Panel 2: Joint Δln L (best wrong mode within EM prior)
ax = axes[1]
for scenario in scenarios:
    joint_deltas = [all_results[(scenario, ncw)]["delta_lnL_joint"] for ncw in ncw_values]
    ax.plot(ncw_values, joint_deltas, "s-", label=scenario_labels[scenario], lw=2, ms=8)

ax.set_xlabel("Number of CW sources")
ax.set_ylabel(r"$\Delta \ln L$ (truth $-$ best wrong in EM prior)")
ax.set_title("Joint: best wrong mode within EM prior")
ax.legend()
ax.grid(True, alpha=0.3)

fig.suptitle(f"Distance importance vs N_CW  |  log10_h = {log10_h}", fontsize=13, y=1.02)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 6. Peak height distributions

For each (scenario, N_CW), show the distribution of wrong-mode peak heights
relative to the truth peak, pooled across all scanned pulsars.
"""),
    code(r"""fig, axes = plt.subplots(1, 3, figsize=(16, 5))

for ax_idx, scenario in enumerate(scenarios):
    ax = axes[ax_idx]

    for ncw in ncw_values:
        key = (scenario, ncw)
        all_deltas = []
        for r in all_results[key]["psr_results"].values():
            if r["n_peaks"] < 2:
                continue
            truth_mask = np.abs(r["peak_distances"] - r["true_dist"]) < 0.3 * r["min_dL"]
            wrong_heights = r["peak_heights"][~truth_mask]
            if len(wrong_heights) > 0:
                all_deltas.extend(r["truth_logl"] - wrong_heights)

        if all_deltas:
            parts = ax.violinplot([all_deltas], positions=[ncw], widths=0.8, showmedians=True)
            for pc in parts["bodies"]:
                pc.set_alpha(0.4)

    ax.set_xlabel("N_CW")
    ax.set_ylabel(r"$\Delta \ln L$ (truth $-$ wrong peak)")
    ax.set_title(scenario_labels[scenario])
    ax.set_xticks(ncw_values)
    ax.grid(True, alpha=0.2)

fig.suptitle(f"Wrong-mode peak height deficit distribution  |  log10_h = {log10_h}",
             fontsize=13, y=1.02)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 7. Signal strength dependence (multi-CW)

Repeat the joint Δln L measurement at three signal strengths for each scenario
with N_CW = 4 (moderate number of sources).
"""),
    code(r"""log10_h_values = [-14.0, -13.5, -13.0]
ncw_test = 4

strength_results = {}

for log10_h_test in log10_h_values:
    for scenario in scenarios:
        key = (scenario, log10_h_test)
        print(f"  {scenario_labels[scenario]}, log10_h = {log10_h_test}")

        rng = np.random.default_rng(seed)
        np.random.seed(seed)

        pta_t, cw_names_t, _ = build_enterprise_pta(ent_psrs, ncw_test)
        inj_t = generate_injection_params(
            pta_t, ent_psrs, ncw_test, cw_names_t, log10_h_test, scenario=scenario, rng=rng,
        )
        sim_t = simulate(pta_t, inj_t)
        rmap_t = {getattr(p, "name", p): y for p, y in zip(pta_t.pulsars, sim_t)}
        logl_t, pkeys_t, bvals_t = build_disco_likelihood(
            disco_psrs, rmap_t, ncw_test, inj_t, cw_names_t,
        )

        joint_t = compute_joint_best_wrong_in_prior(
            logl_t, bvals_t, pkeys_t,
            ent_psrs, disco_psrs, inj_t,
            cw_block_names=cw_names_t, n_sigma=3.0,
        )
        strength_results[key] = joint_t["delta_lnL"]
        print(f"    Delta ln L = {joint_t['delta_lnL']:.2f}")

        jax.clear_caches()

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
for scenario in scenarios:
    deltas = [strength_results[(scenario, h)] for h in log10_h_values]
    ax.plot(log10_h_values, deltas, "o-", label=scenario_labels[scenario], lw=2, ms=8)

ax.set_xlabel("log10(h)")
ax.set_ylabel(r"$\Delta \ln L$ (correct $-$ wrong distances)")
ax.set_title(f"Joint distance effect vs signal strength  |  N_CW = {ncw_test}")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 8. Summary table
"""),
    code(r"""print("=" * 80)
print("SUMMARY: Joint Delta ln L (truth vs best wrong mode within EM prior)")
print("=" * 80)
print(f"\n{'Scenario':25s} | " + " | ".join(f"N_CW={n}" for n in ncw_values))
print("-" * (25 + 3 + len(ncw_values) * 12))
for scenario in scenarios:
    row = f"{scenario_labels[scenario]:25s} | "
    for ncw in ncw_values:
        delta = all_results[(scenario, ncw)]["delta_lnL_joint"]
        row += f" {delta:8.2f}  | "
    print(row)

print(f"\n\nSignal strength dependence (N_CW = {ncw_test}):")
print(f"{'Scenario':25s} | " + " | ".join(f"h={h:5.1f}" for h in log10_h_values))
print("-" * (25 + 3 + len(log10_h_values) * 12))
for scenario in scenarios:
    row = f"{scenario_labels[scenario]:25s} | "
    for h in log10_h_values:
        delta = strength_results[(scenario, h)]
        row += f" {delta:8.2f}  | "
    print(row)
"""),
    md(r"""## Interpretation

**Key findings to look for:**

1. **Single CW (N_CW=1):** If Δln L ≈ 0 between truth peak and best wrong-mode peak,
   then for a single source, knowing the correct distance does not dramatically improve
   the likelihood. Each wrong mode is an equally valid phase assignment.

2. **Multiple CWs:** If Δln L grows with N_CW, this confirms that multiple CW sources
   with shared pulsar distances create a *combinatorial constraint* — the correct distance
   is the only one that simultaneously satisfies all CW pulsar-term phases.

3. **Close vs separated:** If the effect is similar regardless of whether CWs are close
   in frequency/sky or well-separated, then the distance constraint is generic to
   multi-source analyses, not specific to "confused" sources.

4. **Signal strength:** Stronger signals should amplify the effect because the pulsar-term
   SNR is higher, making the oscillation amplitude larger.

Compare with the single-CW analysis in `01_single_cw_distance_landscape.ipynb`.
"""),
]

# ====================================================================
# Write notebooks
# ====================================================================

write_nb("01_single_cw_distance_landscape.ipynb", nb1_cells)
write_nb("02_multi_cw_distance_landscape.ipynb", nb2_cells)
print("Done!")
