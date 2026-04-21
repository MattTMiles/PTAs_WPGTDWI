"""Generate notebook varying number of pulsars (fixed N_CW=10, well-separated, h=-13)."""
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
    md(r"""# Distance Constraint vs Number of Pulsars

**Motivation:** Notebooks 00–03 showed that well-separated CW sources with strong
signals (log10_h = -13) produce the largest Δln L between truth and wrong pulsar
distances. This notebook fixes the CW configuration (N_CW = 10, well-separated,
log10_h = -13) and instead varies the **number of pulsars** in the array.

The joint Δln L metric evaluates: for each pulsar, find the best non-truth distance
mode within the EM prior (±3σ), then measure the total likelihood penalty when ALL
pulsars are simultaneously shifted to their respective best wrong modes.

More pulsars means more independent distance constraints — each pulsar provides its
own set of distance modes with different spacings (because each pulsar has a different
sky position relative to the CW sources). The question is: how does the joint penalty
scale with array size?

**Cases tested:**
- N_pulsars = 2, 4, 7, 10, 15, 20, 30, 50
- Both noise-free and noisy (realistic white + red noise + GWB)
"""),
    code(r"""import sys, importlib
sys.path.insert(0, ".")
import cw_helpers
importlib.reload(cw_helpers)
from cw_helpers import *
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

%matplotlib inline
plt.rcParams.update({"figure.dpi": 120, "font.size": 11})
"""),
    md("## 1. Configuration"),
    code(r"""npsr_values = [2, 4, 7, 10, 15, 20, 30, 50]
ncw = 10
log10_h = -13.0
seed = 1234
scenario = "well_separated"

# Pre-load the maximum number of pulsars we'll need
max_pulsars = max(npsr_values)
all_ent_psrs, all_disco_psrs = load_pulsars(max_pulsars)
print(f"Loaded {len(all_ent_psrs)} pulsars (max needed: {max_pulsars})")
print(f"Configuration: N_CW={ncw}, log10_h={log10_h}, scenario={scenario}")

# Show the pulsars we'll use
print(f"\nPulsar pool:")
for i, p in enumerate(all_ent_psrs):
    print(f"  {i+1:2d}. {p.name:16s}  pdist = {p.pdist[0]:.3f} +/- {p.pdist[1]:.3f} kpc")
"""),
    md(r"""## 2. Noise-free case

For each array size, inject N_CW=10 well-separated CW sources (using the same
random seed so the CW parameters are identical across runs), then measure the
joint Δln L. Only the number of pulsars changes.
"""),
    code(r"""noisefree_results = {}

for npsr in npsr_values:
    ent_psrs = all_ent_psrs[:npsr]
    disco_psrs = all_disco_psrs[:npsr]

    print(f"\n--- N_pulsars = {npsr} ---")

    rng_local = np.random.default_rng(42)

    # Generate CW parameters (same for all array sizes)
    def make_cw_params(n_sources, rng):
        sources = []
        for i in range(n_sources):
            sources.append({
                "cos_gwtheta": rng.uniform(-1, 1),
                "gwphi": rng.uniform(0, 2 * np.pi),
                "cos_inc": rng.uniform(-1, 1),
                "log10_mc": 9.0,
                "log10_fgw": rng.uniform(-8.5, -7.5),
                "log10_h": log10_h,
                "phase0": rng.uniform(0, 2 * np.pi),
                "psi": rng.uniform(0, np.pi),
            })
        return sources

    cw_params = make_cw_params(ncw, rng_local)

    # Inject noise-free CW
    rmap, _ = inject_noisefree_cw(disco_psrs, cw_params)
    logl_fn, pkeys, bvals = build_noisefree_likelihood(disco_psrs, rmap, ncw, cw_params, include_gwb=False)

    ll_truth = float(logl_fn(bvals))
    print(f"  ln L(truth) = {ll_truth:.2f}")

    # Joint best wrong mode within EM prior
    joint = compute_joint_best_wrong_in_prior(
        logl_fn, bvals, pkeys,
        disco_psrs, disco_psrs, cw_params, n_sigma=3.0,
    )

    noisefree_results[npsr] = {
        "ll_truth": joint["ll_truth"],
        "ll_best_wrong": joint["ll_best_wrong"],
        "delta_joint": joint["delta_lnL"],
        "per_pulsar": joint["per_pulsar"],
    }

    print(f"  ln L(best wrong) = {joint['ll_best_wrong']:.2f}")
    print(f"  Joint Delta ln L = {joint['delta_lnL']:.2f}")

    # Per-pulsar breakdown
    for pname, r in joint["per_pulsar"].items():
        print(f"    {pname:16s}: Delta={r['delta_lnL']:8.2f}, peaks_in_prior={r['n_peaks_in_prior']:4d}")

    jax.clear_caches()

print("\nNoise-free complete!")
"""),
    md(r"""## 3. Noisy case

Same CW sources but with realistic noise simulation via enterprise.
"""),
    code(r"""noisy_results = {}

for npsr in npsr_values:
    ent_psrs = all_ent_psrs[:npsr]
    disco_psrs = all_disco_psrs[:npsr]

    print(f"\n--- N_pulsars = {npsr} ---")

    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    # Build enterprise PTA
    pta, cw_names, Tspan = build_enterprise_pta(ent_psrs, ncw)

    # Generate injection (well-separated scenario)
    inj_params = generate_injection_params(
        pta, ent_psrs, ncw, cw_names, log10_h, scenario=scenario, rng=rng,
    )

    # Simulate
    sim_resids = simulate(pta, inj_params)
    resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}

    # Build discovery likelihood
    logl_fn, pkeys, bvals = build_disco_likelihood(
        disco_psrs, resid_map, ncw, inj_params, cw_names, include_gwb=False,
    )

    ll_truth = float(logl_fn(bvals))
    print(f"  ln L(truth) = {ll_truth:.2f}")

    # Joint best wrong mode within EM prior
    joint = compute_joint_best_wrong_in_prior(
        logl_fn, bvals, pkeys,
        ent_psrs, disco_psrs, inj_params,
        cw_block_names=cw_names, n_sigma=3.0,
    )

    noisy_results[npsr] = {
        "ll_truth": joint["ll_truth"],
        "ll_best_wrong": joint["ll_best_wrong"],
        "delta_joint": joint["delta_lnL"],
        "per_pulsar": joint["per_pulsar"],
    }

    print(f"  ln L(best wrong) = {joint['ll_best_wrong']:.2f}")
    print(f"  Joint Delta ln L = {joint['delta_lnL']:.2f}")

    for pname, r in joint["per_pulsar"].items():
        print(f"    {pname:16s}: Delta={r['delta_lnL']:8.2f}, peaks_in_prior={r['n_peaks_in_prior']:4d}")

    jax.clear_caches()

print("\nNoisy complete!")
"""),
    md(r"""## 4. Joint Δln L vs number of pulsars
"""),
    code(r"""fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: absolute Δln L
ax = axes[0]
nf_deltas = [noisefree_results[n]["delta_joint"] for n in npsr_values]
ny_deltas = [noisy_results[n]["delta_joint"] for n in npsr_values]

ax.plot(npsr_values, nf_deltas, "o-", label="Noise-free", lw=2, ms=8, color="steelblue")
ax.plot(npsr_values, ny_deltas, "s-", label="Noisy", lw=2, ms=8, color="firebrick")
ax.set_xlabel("Number of pulsars")
ax.set_ylabel(r"Joint $\Delta \ln L$ (truth $-$ best wrong in EM prior)")
ax.set_title(f"N_CW={ncw}, well-separated, log10_h={log10_h}")
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xticks(npsr_values)

# Right: Δln L per pulsar (normalized)
ax = axes[1]
nf_per_psr = [noisefree_results[n]["delta_joint"] / n for n in npsr_values]
ny_per_psr = [noisy_results[n]["delta_joint"] / n for n in npsr_values]

ax.plot(npsr_values, nf_per_psr, "o-", label="Noise-free", lw=2, ms=8, color="steelblue")
ax.plot(npsr_values, ny_per_psr, "s-", label="Noisy", lw=2, ms=8, color="firebrick")
ax.set_xlabel("Number of pulsars")
ax.set_ylabel(r"Joint $\Delta \ln L$ / N_pulsars")
ax.set_title(f"Per-pulsar contribution (normalized)")
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xticks(npsr_values)

fig.tight_layout()
plt.show()
"""),
    md(r"""## 5. Per-pulsar breakdown

Show each pulsar's individual contribution to the joint penalty as the array grows.
This reveals whether some pulsars are much more constraining than others (e.g. due
to tighter EM distance priors or more favourable sky positions relative to the CW
sources).
"""),
    code(r"""fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for col, (results, label) in enumerate([
    (noisefree_results, "Noise-free"),
    (noisy_results, "Noisy"),
]):
    ax = axes[col]

    # Collect per-pulsar Δln L for largest array to rank pulsars
    largest = max(npsr_values)
    psr_names = list(results[largest]["per_pulsar"].keys())
    psr_deltas = {pname: results[largest]["per_pulsar"][pname]["delta_lnL"] for pname in psr_names}
    sorted_psrs = sorted(psr_deltas.keys(), key=lambda x: psr_deltas[x], reverse=True)

    # For each N_pulsars, show cumulative per-pulsar contributions
    for npsr in npsr_values:
        per_psr = results[npsr]["per_pulsar"]
        names = list(per_psr.keys())
        deltas = [per_psr[n]["delta_lnL"] for n in names]
        ax.barh(
            [f"{n} (N={npsr})" for n in names],
            deltas,
            alpha=0.15,
        )

    # Show breakdown for largest array clearly
    per_psr = results[largest]["per_pulsar"]
    names = list(per_psr.keys())
    deltas = [per_psr[n]["delta_lnL"] for n in names]
    sort_idx = np.argsort(deltas)[::-1]
    sorted_names = [names[i] for i in sort_idx]
    sorted_deltas = [deltas[i] for i in sort_idx]

    ax.barh(sorted_names, sorted_deltas, color="steelblue" if col == 0 else "firebrick", alpha=0.7)
    ax.set_xlabel(r"Per-pulsar $\Delta \ln L$ (truth $-$ best wrong)")
    ax.set_title(f"{label} — N_pulsars={largest}")
    ax.grid(True, alpha=0.2, axis="x")

fig.suptitle(f"Per-pulsar distance constraints | N_CW={ncw}, well-separated, log10_h={log10_h}",
             fontsize=13, y=1.02)
fig.tight_layout()
plt.show()
"""),
    md(r"""## 6. Scaling analysis

Does Δln L scale linearly with N_pulsars (each pulsar contributes independently)
or is there a different scaling?
"""),
    code(r"""fig, ax = plt.subplots(1, 1, figsize=(8, 5))

nf_deltas = np.array([noisefree_results[n]["delta_joint"] for n in npsr_values])
ny_deltas = np.array([noisy_results[n]["delta_joint"] for n in npsr_values])
npsr_arr = np.array(npsr_values, dtype=float)

# Linear fit
from numpy.polynomial import polynomial as P
nf_coeffs = P.polyfit(npsr_arr, nf_deltas, 1)
ny_coeffs = P.polyfit(npsr_arr, ny_deltas, 1)

ax.plot(npsr_values, nf_deltas, "o", color="steelblue", ms=10, label="Noise-free")
ax.plot(npsr_values, ny_deltas, "s", color="firebrick", ms=10, label="Noisy")

x_fit = np.linspace(0, max(npsr_values) * 1.1, 100)
ax.plot(x_fit, P.polyval(x_fit, nf_coeffs), "--", color="steelblue", alpha=0.5,
        label=f"Linear fit (slope={nf_coeffs[1]:.1f}/pulsar)")
ax.plot(x_fit, P.polyval(x_fit, ny_coeffs), "--", color="firebrick", alpha=0.5,
        label=f"Linear fit (slope={ny_coeffs[1]:.1f}/pulsar)")

ax.set_xlabel("Number of pulsars")
ax.set_ylabel(r"Joint $\Delta \ln L$")
ax.set_title(f"Scaling: N_CW={ncw}, well-separated, log10_h={log10_h}")
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, max(npsr_values) * 1.1)
ax.set_ylim(0, None)

fig.tight_layout()
plt.show()

print(f"Noise-free: slope = {nf_coeffs[1]:.1f} Δln L per pulsar, intercept = {nf_coeffs[0]:.1f}")
print(f"Noisy:      slope = {ny_coeffs[1]:.1f} Δln L per pulsar, intercept = {ny_coeffs[0]:.1f}")
"""),
    md("## 7. Summary table"),
    code(r"""print("=" * 90)
print(f"SUMMARY: Joint Delta ln L vs N_pulsars")
print(f"N_CW={ncw}, well-separated, log10_h={log10_h}")
print("=" * 90)

print(f"\n{'N_pulsars':>10s} | {'Noise-free':>14s} | {'Noisy':>14s} | {'NF/pulsar':>12s} | {'Noisy/pulsar':>12s}")
print("-" * 72)
for npsr in npsr_values:
    nf = noisefree_results[npsr]["delta_joint"]
    ny = noisy_results[npsr]["delta_joint"]
    print(f"{npsr:10d} | {nf:14.2f} | {ny:14.2f} | {nf/npsr:12.2f} | {ny/npsr:12.2f}")
"""),
    md(r"""## Interpretation

- If Δln L scales **linearly** with N_pulsars, each pulsar contributes independently
  to the distance constraint. This would mean the per-pulsar Δln L is roughly constant.

- If Δln L scales **sub-linearly**, additional pulsars provide diminishing returns —
  perhaps because newly added pulsars have worse EM priors or less favourable geometry.

- If Δln L scales **super-linearly**, pulsars are reinforcing each other — the joint
  constraint is more powerful than the sum of individual constraints.

The per-pulsar breakdown (Section 5) reveals which pulsars dominate the constraint
and whether some are effectively free-riders with negligible individual Δln L.
"""),
]

with open("04_vary_pulsars.ipynb", "w") as f:
    json.dump(notebook(cells), f, indent=1)
print("Wrote 04_vary_pulsars.ipynb")
