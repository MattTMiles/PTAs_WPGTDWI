# CW ESS Improvement Experiments — Executive Summary

## Goal
Improve the effective sample size (ESS) for CW parameters in the phase-connected
pulsar timing array sampler. Starting ESS was ~16-23 for the tightest CW params
(cos_gwtheta, gwphi, log10_mc, log10_fgw).

## Best Configuration
**MK9 with Empirical Covariance Adaptive + Newton Distance Snap**:
- Test 1 (truth): CW ESS min=62, coverage 13/13
- Test 2 (shifted): CW ESS min=62, coverage 13/13
- **4x improvement** over MK5 baseline (ESS 16 -> 62-73)
- Speed: ~400 it/s (slower than MK5's 950 due to Newton snap gradient evals)
- Key features:
  - Phase 1 (burn-in): Fisher eigenmodes + Newton snap + aggressive adaptation
  - Phase 2: Empirical covariance eigenmodes (from burn-in chain samples)
  - Newton distance snapping (3 grad evals per CW proposal)
  - Prior-snap big jumps for distance exploration
  - Coherent freq-dist proposal with Newton snap

## Root Cause of ESS Bottleneck (SOLVED)

The 8x8 Fisher conditional Hessian gives eigenmode proposal widths that are
**100-15000x too small** for the bottleneck CW params:

| Param | Chain std | Fisher eig sig | Ratio |
|-------|----------|----------------|-------|
| cos_gwtheta | 9.1e-3 | 6.2e-7 | 14,800x |
| gwphi | 1.5e-3 | 1.3e-5 | 109x |
| cos_inc | 4.1e-2 | 1.6e-5 | 2,600x |
| log10_mc | 9.6e-3 | 1.3e-3 | 7x |

**Why**: The Hessian measures sensitivity at fixed distances. When distances are fixed,
tiny changes in sky position cause huge phase shifts. But in reality, distances
mode-hop to compensate — the effective posterior is much wider.

**Fix**: Switch to empirical chain covariance at mid-burn-in. The empirical covariance
captures the actual posterior geometry including distance mode-hopping.

## Experiments Chronology

### MK5 (Baseline): Fisher Eigenmode Proposals
- CW ESS min=16-23, ~950 it/s
- 5 move types: eigenmode (35%), dist within (10%), big jump (30%), freq-dist (10%), CW joint (15%)

### MK6: CW Eigenmode Sweep After Big Jumps
**Result: Worse** — ESS same, coverage dropped, speed dropped to ~630 it/s

### MK7: Cross-Hessian CW Coupling in Big Jumps
**Result: Much worse** — killed big jump acceptance (0.005 vs 0.40)

### MK8: HMC for 8 CW Dimensions
**Result: Mixed** — improved "easy" CW ESS (psi: 20→907) but bottleneck unchanged (ESS ~16)
- HMC improved within-mode mixing but couldn't address CW-distance coupling
- 7-12x slower per iteration due to gradient evaluations

### MK9 v1: Eigenmode-Distance Coherent Proposals (Newton Snap)
**Result: Better AR but same ESS** — CW Eigen AR 0.35→0.71, but ESS still 16
- Newton snap correctly adjusts distances after CW changes
- BUT Fisher eigenmode widths are 100-15000x too small for bottleneck params
- Half of eigenmode proposals (stiff modes 0-3) take near-zero steps with AR≈1.0

### MK9 v2: Marginal CW Covariance (Schur Complement)
**Result: Worse** — full 13x13 Hessian too ill-conditioned for reliable inversion
- CW Joint AR dropped to 0.0%, "easy" param ESS degraded

### MK9 v3 (FINAL): Empirical Covariance + Newton Snap
**Result: 4x improvement** — CW ESS min 16 → 62-73
- Empirical covariance from burn-in chain gives correctly-sized proposals
- Newton snap enables wider CW exploration (distances adjust to compensate)
- All 8 eigenmodes have similar AR (~0.15-0.27)
- 13/13 coverage in both tests

## Key Findings

1. **Fisher eigenvalues are misleading**: The 8x8 conditional Hessian eigenvalues
   span 10^10, but the actual posterior eigenvalue range is ~10^5. The stiff modes
   (eigenvalues 10^10-10^13) reflect CW-distance coupling, not true posterior narrowness.

2. **Newton distance snapping is important**: Without it, cos_gwtheta explores a 9x
   narrower range (std 4.1e-4 vs 3.7e-3). Newton enables wider CW exploration by
   keeping distances locked to peaks as CW params change.

3. **Empirical covariance is the key**: Adaptive Metropolis (switching to empirical
   chain covariance) directly solves the mismatch between Fisher and actual proposal widths.

4. **ESS scales linearly with N**: Autocorrelation time is constant (~750 steps).
   Running longer chains gives proportionally more ESS.

## File Locations
- `data_likelihood_sandbox_MK9.ipynb` — **Best sampler** (recommended)
- `data_likelihood_sandbox_MK5.ipynb` — Previous best baseline
- `data_likelihood_sandbox_MK8.ipynb` — HMC experiment (not recommended)
- `data_likelihood_sandbox_MK6.ipynb` — CW sweep experiment (worse)
- `data_likelihood_sandbox_MK7.ipynb` — Cross-Hessian coupling (much worse)
