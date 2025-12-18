# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is Q1 of the HSYMT (Hasasia-styled Synthetic MAP Tests) project, focusing on **MAP-based ring PTA experiments** for continuous gravitational wave (CW) detection and localization. The codebase simulates pulsar timing array (PTA) data with injected CW signals and performs maximum a posteriori (MAP) optimization using JAX-accelerated likelihoods.

**Primary Research Question (Q1):**

**How does detection & localization change with pulsar distance constraints?**

This question investigates the interplay between pulsar distance prior uncertainties and CW signal recovery. The analysis quantifies:
1. How CW phase evolution breaks distance-parameter degeneracies
2. The gains from improved distance measurements (e.g., GAIA will improve by factor of 1.6)
3. Detection and localization performance across varying distance precision

**Related Project Questions:**
- **Q2**: How do pulsar distance constraints improve as you add CW information? (multi-source, joint vs factorized models)
- **Q3**: What changes when phase-connected pulsar-term models cannot be used? (under consideration)
- **Q4**: How does precise distance knowledge for some pulsars help infer distances of others?

## Q1 Experimental Design

**Baseline Configuration:**
- **Ring PTA**: 30 pulsars arranged perpendicular to CW source direction
- **CW Signal**: Single source at moderate frequency (log10_fgw ~ -8) to allow evolution without extreme phase wrapping
- **SNR Targets**: 5, 10, 20 (SNR used as control parameter)

**Distance Prior Tiers:**
- **Real** (scale=1.0): Current IPTA DR2-like uncertainties
- **Tier 1 Loosened** (scale=3.0 or 10.0): Shows what CW can do alone, quantifies degeneracy breaking from phase evolution, relevant for future PTA pulsars
- **Tier 2 Tightened** (scale<1.0, e.g., 0.6): Simulates improved measurements; GAIA alone gives factor of 1.6 improvement
- **Tier 3 Exact** (scale=0.0): Dream scenario with near-perfect distance knowledge

**Noise Scenarios:**
- **Scenario A**: White noise only (WN) - establish CW SNR baseline
- **Scenario B**: White noise + red noise (WN+RN) - more realistic single-pulsar noise
- **Scenario C**: White noise + red noise + gravitational wave background (WN+RN+GWB) - full realistic PTA

**Metrics to Measure:**
1. **Detection Statistic**: Bayes factor or Δlog(L) to assess signal detectability
2. **Sky Localization Area**: 90% credible region; check if true position is recovered
3. **Distance Recovery**:
   - Posterior σ(dist) / injected_dist (absolute performance)
   - Prior σ(dist) / posterior σ(dist) (CW-driven improvement)
4. **Credible Level at Injection**: Assess calibration of posterior uncertainties

## Core Architecture

### Three-Module Design

1. **`sim.py`**: Simulation framework
   - Loads pulsars from `data_products/` (feather files)
   - Constructs ring geometries around CW source directions
   - Builds Enterprise PTAs with noise models (WN, RN, GWB)
   - Simulates residuals with injected CW signals
   - Mirrors data into Discovery likelihood objects (JAX-compatible)

2. **`optimize.py`**: MAP optimization engines
   - Simple momentum-based gradient ascent (`run_map_gradient_ascent`)
   - NumPyro SVI with AutoDelta guide (`run_svi_map`) - requires `pulsar_map_noise_estimates` package
   - Optional LBFGS refinement via `jaxopt` (`run_lbfgs`)

3. **`run_q1_map.py`**: Command-line entry point
   - Orchestrates full pipeline: draw parameters → simulate → optimize → report
   - See Quick Start below for usage

### Data Flow

```
Enterprise PTA (simulation) → residuals → Discovery likelihood (JAX) → MAP optimizer → results
```

**Enterprise vs Discovery:**
- **Enterprise**: Used for simulation only (defines noise models, timing models, CW signals)
- **Discovery**: Used for likelihood evaluation during optimization (JAX-compiled, GPU-compatible)
- The two frameworks must be kept in sync: ring positions, residuals, and parameter naming conventions

### Scenario System

Three noise scenarios controlled by `Q1Config.scenario` (see Q1 Experimental Design above):
- **A**: White noise only (WN)
- **B**: White noise + red noise (WN+RN)
- **C**: White noise + red noise + gravitational wave background (WN+RN+GWB)

These scenarios are applied sequentially in Q1 analysis: first find CW SNR in WN-only, then add RN, then add GWB.

### Distance Tiers

Control pulsar distance prior uncertainty via `DistanceTier(name, scale)`:
- `scale=1.0`: Real uncertainties (IPTA DR2-like)
- `scale>1.0`: Loosen priors (e.g., 3.0 or 10.0 for Tier 1)
- `scale<1.0`: Tighten priors (e.g., 0.6 for GAIA-improved Tier 2)
- `scale=0.0`: Near-exact distances (Tier 3, sigma → 1e-6)

**Implementation Note:** Distance tiers are applied by multiplying the pulsar's original distance sigma in `apply_distance_tier()` (sim.py:178). The tier scale directly affects the Gaussian prior width used in likelihood evaluation.

### SNR Targeting (Planned)

The Q1 experimental design calls for SNR control at specific values (5, 10, 20). This requires:
1. Computing SNR for a given log10_h injection
2. Iteratively scaling log10_h to hit target SNR
3. This functionality is listed as a "next step" in README_Q1.md:23

Currently, log10_h is set manually in `draw_cw_population()` (sim.py:131). To implement SNR targeting, you'll need to add an SNR calculation utility and a root-finding loop to adjust log10_h before finalizing the injection parameters.

## Quick Start

### Running a single MAP optimization

```bash
cd /home/mattm/projects/HSYMT/MAIN_PROJECT_QUESTIONS/Q1
python run_q1_map.py --scenario B --tier-scale 1.0 --tier-name real --opt-mode svi
```

**Key flags:**
- `--scenario {A,B,C}`: Noise scenario
- `--tier-scale FLOAT`: Distance sigma multiplier
- `--tier-name STR`: Label for bookkeeping
- `--n-pulsars INT`: Ring size (default 30)
- `--opt-mode {svi,map}`: Optimizer choice (svi recommended)
- `--steps INT`: Gradient ascent iterations (ignored for svi)
- `--lr FLOAT`: Learning rate for map mode
- `--seed INT`: Reproducibility seed

### Interactive exploration

Use Jupyter notebooks:
- `run_q1_map.ipynb`: Interactive MAP pipeline with parameter scans and LBFGS refinement
- `run_q1_sample.ipynb`: Sampling experiments (if applicable)
- `run_q1_det_stat.ipynb`: Detection statistic analysis

## Key Implementation Details

### Ring Geometry Construction

Ring pulsars are placed perpendicular to the mean CW direction(s):
- Radius auto-expands if multiple CW sources are separated
- Default radius: 20° (configurable via `Q1Config.ring_radius_deg`)
- Implementation: `build_ring_positions()` in sim.py:90

### Phase-Connected Models

By default, CW signals use **phase-connected pulsar terms**:
- `Q1Config.include_pterm=True`: Include pulsar term
- `Q1Config.inj_evolve=True`: Inject evolving signals
- `Q1Config.rec_evolve=True`: Recover with evolving model
- Discovery uses `make_phase_connected_binary()` from discovery package

### TOA Error Forcing

Before simulation, TOA errors are forced to 1e-6 via `set_toaerrs()` to avoid numerical outliers. This prevents Enterprise's sparse Cholesky from failing on ill-conditioned covariance matrices.

### JAX Configuration

`sim.py` sets critical JAX environment variables:
```python
XLA_PYTHON_CLIENT_PREALLOCATE=false  # Avoid large GPU memory grabs
XLA_PYTHON_CLIENT_ALLOCATOR=platform
JAX_DISABLE_MMAP_CACHE=1
XLA_FLAGS=--xla_gpu_autotune_level=2
jax.config.update("jax_enable_x64", True)  # Required for likelihood precision
```

Do not modify these without understanding implications for numerical stability.

### Distance Marginalization

`sim.py` provides `build_discovery_loglike_marginalized()` (experimental):
- Analytically marginalizes pulsar distances using Gaussian approximation
- Pulsar term evaluated at prior mean, attenuated by coherence factor κ
- κ = exp(-0.5 * [2ω₀(1-cos μ)(kpc/c)σ_d]²)
- Removes distance parameters from optimization (reduces dimensionality)

## Dependencies

**Required:**
- `enterprise` / `enterprise_extensions`: PTA simulation framework
- `discovery`: JAX-based likelihood library
- `jax` / `jax.numpy`: Automatic differentiation and GPU support
- `numpy`, `scipy`, `matplotlib`
- `sksparse.cholmod`: Sparse Cholesky for large covariance matrices

**Optional:**
- `numpyro`: For SVI-based MAP (`--opt-mode svi`)
- `pulsar_map_noise_estimates`: SVI helper functions (see fallback path in optimize.py:108)
- `jaxopt`: For LBFGS refinement

**Data Requirements:**
- Pulsar feather files in `/home/mattm/projects/HSYMT/data_products/`
- Default path: `Q1Config.feather_dir` (two levels up from Q1/)

## Common Development Patterns

### Adding a new optimizer

1. Implement function in `optimize.py` with signature:
   ```python
   def run_new_optimizer(
       logpost: Callable[[Array, Sequence[str]], float],
       param_keys: Sequence[str],
       init_params: Params,
       **kwargs
   ) -> MAPResult
   ```

2. Add option to `run_q1_map.py` argparse and conditional

3. Update notebooks to test interactively

### Adding a new scenario

1. Extend `Q1Config.scenario` choices in sim.py:64
2. Update `draw_noise_params()` to handle new scenario's noise parameters
3. Modify `build_enterprise_pta()` to add new signal blocks if needed
4. Mirror changes in `build_discovery_loglike()` for Discovery side

### Debugging likelihood mismatches

If Enterprise simulation and Discovery likelihood disagree:
1. Check ring positions are identical: `apply_ring_positions()` called on both
2. Verify residuals copied correctly: `residual_map` in `build_discovery_loglike()`
3. Ensure distance priors match: `apply_distance_tier()` vs Discovery's `psr.pdist`
4. Check parameter naming: Enterprise uses `{name}_{param}`, Discovery must match exactly

### Parameter naming conventions

- CW global params: `cw_{param}` (single source) or `cw_{param}_{idx+1}` (multi-source)
  - Examples: `cw_log10_h`, `cw_gwphi`, `cw2_log10_fgw`
- Pulsar distances: `{psr_name}_cw_p_dist`
- Noise params: `{psr_name}_rednoise_log10_A`, `gwb_log10_A`, etc.

## Project Context

**Q1's Role in Broader HSYMT Study:**

Q1 is the first of four targeted research questions investigating CW detection/localization and pulsar distance constraints. The broader project structure:

- **Q1** (this directory): Distance constraints → detection/localization performance
- **Q2** (planned): CW information → improved pulsar distance inference (multi-source, joint vs factorized)
- **Q3** (under consideration): Phase-unconnected vs phase-connected models
- **Q4** (planned): Distance constraint propagation across pulsar arrays

**Relationship to Parent Directory:**

The parent directory (`/home/mattm/projects/HSYMT/`) contains exploratory work that informed Q1's design:
- **Exploratory notebooks**: `data_likelihood_sandbox.ipynb`, `hessian_check_*.ipynb`, `GWB_recovery_tests.ipynb`
- **Terminal scripts**: `lnL_GW_recovery_phase_connected*.py` (ring PTA prototypes)
- **Result archives**: `lnLs_master_tests_ringPTA/`, `lnLs_GWAmp*/` (earlier experiments)

Q1 scaffolding refactors these exploratory scripts into reusable modules (`sim.py`, `optimize.py`) for systematic parameter sweeps.

**Important Design Philosophy:**

From `main_project_questions_readme.md:48-49`: The ring PTA should respect the real PTA's distance uncertainty distribution. When drawing pulsars from the larger sample, retain the fractional uncertainty distribution to ensure "fair" comparisons between tiers.

Git branch: `MM_playground` (see gitStatus above for current state)

## Planned Sweep Infrastructure (To Be Implemented)

The Q1 experimental design requires sweeps over SNR × distance_tier × scenario, with result logging for analysis. Key missing pieces:

1. **SNR targeting**: Automatically scale log10_h to hit target SNR (5, 10, 20)
2. **Sweep orchestration**: Loop over parameter combinations, running multiple realizations per configuration
3. **Result logging**: Save MAP results, detection statistics, and metrics to structured format (npz, pkl, or HDF5)
4. **Analysis pipeline**: Compute metrics from saved results:
   - Detection statistics (BF or Δlog(L))
   - Sky localization area (90% credible regions from Hessian or sampling)
   - Distance recovery metrics (posterior σ vs prior σ, bias)
5. **Plotting utilities**: Visualize results as functions of SNR and distance tier

**Suggested Implementation Approach:**

Create `run_q1_sweep.py` that:
- Wraps `run_q1_map.py` logic in nested loops (SNR, tier, scenario, realization)
- Saves results to `results/` subdirectory with systematic naming
- Logs runtime, convergence status, and basic metrics

Create `analyze_q1_results.py` that:
- Loads results from sweep
- Computes detection statistics and credible regions
- Generates plots for paper/presentations

For Hessian-based credible regions, see `run_q1_map.ipynb` cells with LBFGS-B refinement and parameter scans.

## Caveats and Known Issues

1. **SVI requires external package**: `pulsar_map_noise_estimates` must be installed or in path (see optimize.py:108 for fallback)
2. **LBFGS can drift**: Unbounded LBFGS may push distance parameters to extreme values; use bounded variants (LBFGSB) with prior-informed bounds
3. **Sparse Cholesky failures**: If simulation crashes, check for negative TOA errors or singular covariance matrices
4. **GPU memory**: XLA may preallocate large chunks; set environment variables as in sim.py:36-40
5. **Multi-source indexing**: When `n_cw_sources>1`, parameter suffixes start at 1 for second source (not 0)
