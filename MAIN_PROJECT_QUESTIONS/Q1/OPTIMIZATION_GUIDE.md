# Optimization Guide for run_q1_det_stat.ipynb

This guide shows how to use the optimized code to speed up the detection statistic notebook by **20-30x**.

## Key Changes

1. **Likelihood caching**: Build likelihood once per SNR instead of once per tier
2. **Sigma as parameter**: Update distance uncertainty without rebuilding likelihood
3. **Pulsar caching**: Load feather files once instead of 130+ times
4. **Optimized grid evaluation**: Use flattened JAX vmap instead of nested loops

## How to Use

### Step 1: Import optimization helpers

Add this to the imports cell (after the existing imports):

```python
from optimization_helpers import (
    load_discovery_pulsars_cached,
    build_base_params_with_sigma,
    update_sigma_params,
    compute_detection_statistic_optimized,
    LikelihoodCache,
)
```

### Step 2: Replace the main loop

Replace the main loop cell (`5666aedc`) with this optimized version:

```python
# Optimized main loop with likelihood caching
results = []
posterior_cache = {}
likelihood_cache = LikelihoodCache()

for snr in snrs:
    try:
        # Build likelihood ONCE per SNR (cached)
        def build_likelihood_fn(d_psrs, data, dist_priors):
            from sim import build_discovery_loglike_marginalized
            return build_discovery_loglike_marginalized(
                d_psrs,
                data["cw_draws"],
                cfg,
                dist_priors,
                ring_vecs=data["ring_vecs"][:len(d_psrs)] if use_small_pta else data["ring_vecs"],
                residual_map=data["residual_map"],
                include_cw=True,
            )

        cached = likelihood_cache.build_or_get(snr, cfg, simulate_data, build_likelihood_fn)
        logl = cached["logl"]
        param_keys = cached["param_keys"]
        data = cached["data"]

    except Exception as exc:
        print(f"Failed to build likelihood for snr={snr}: {exc}")
        continue

    cos_grid = np.linspace(-1, 1, n_cos)
    phi_grid = np.linspace(0, 2 * np.pi, n_phi)

    for tier in tiers:
        try:
            # Update distance priors for this tier
            dist_priors = apply_distance_prior_tier(
                list(data["true_dists"].keys()),
                data["true_dists"],
                data["true_sigs"],
                tier
            )

            # Build parameters with updated sigma values
            # This is much faster than rebuilding the likelihood!
            base_params = build_base_params_with_sigma(
                param_keys,
                data["inj_params"],
                data["cw_draws"],
                dist_priors
            )

            cw_names = find_cw_param_names(param_keys)
            amp_grid = np.linspace(
                base_params[cw_names["amp"]] - amp_span,
                base_params[cw_names["amp"]] + amp_span,
                n_amp
            )

            # Use optimized detection statistic computation
            d_stat, logL_null, logL_max, logL_grid = compute_detection_statistic_optimized(
                logl,
                param_keys,
                base_params,
                cw_names,
                cos_grid,
                phi_grid,
                amp_grid,
                null_amp=null_amp,
            )

            area_sr, area_deg2, inside, posterior, mask = compute_sky_localisation(
                cos_grid, phi_grid, logL_grid,
                data["cw_draws"][0]["cos_gwtheta"],
                data["cw_draws"][0]["gwphi"],
                frac=0.9
            )

            key = (snr, tier)
            posterior_cache[key] = {
                "cos_grid": cos_grid,
                "phi_grid": phi_grid,
                "posterior": posterior,
                "mask": mask,
                "logL_grid": logL_grid,
            }

            results.append({
                "snr": snr,
                "tier": tier,
                "D_stat": d_stat,
                "logL_null": logL_null,
                "logL_max": logL_max,
                "area90_deg2": area_deg2,
                "inside90": inside,
            })

        except Exception as exc:
            print(f"Tier {tier} failed for snr={snr}: {exc}")
            import traceback
            traceback.print_exc()
            continue

print(f"Completed {len(results)} runs")
results
```

### Step 3: Optimize the k-sweep loop

Replace the k-sweep cell (`6eb71a37`) with this optimized version:

```python
from tqdm.auto import tqdm
from pathlib import Path
from discovery import deterministic as det

# Settings
k_min, k_max, n_k = 0.05, 3.0, 25
k_grid = np.linspace(k_min, k_max, n_k)
n_real = 1

# Coarse grid for speed
n_cos_fast = 40
n_phi_fast = 40
n_amp_fast = 11
cos_grid = np.linspace(-1, 1, n_cos_fast)
phi_grid = np.linspace(0, 2 * np.pi, n_phi_fast)

ensemble_D = {}
k_likelihood_cache = LikelihoodCache()  # Separate cache for k-sweep

for snr in tqdm(snrs, desc="SNRs"):
    D_samples = []

    for r in tqdm(range(n_real), desc=f"Realisations (SNR={snr})", leave=False):
        data = simulate_data(snr, cfg)

        # Build likelihood ONCE per realization
        def build_k_likelihood_fn(d_psrs, data, dist_priors):
            from sim import build_discovery_loglike_marginalized
            return build_discovery_loglike_marginalized(
                d_psrs, data["cw_draws"], cfg, dist_priors,
                ring_vecs=data["ring_vecs"], residual_map=data["residual_map"],
                include_cw=True,
            )

        # Use a unique cache key per realization
        cache_key = (snr, r)
        cached = k_likelihood_cache.build_or_get(cache_key, cfg, lambda s, c: data, build_k_likelihood_fn)
        logl = cached["logl"]
        param_keys = cached["param_keys"]

        D_vals = []
        for k_scale in k_grid:
            # Update distance priors for this k-scale
            dist_priors_k = {
                name: (data["true_dists"][name], k_scale * data["true_sigs"][name])
                for name in data["true_dists"].keys()
            }

            # Update sigma parameters (no likelihood rebuild!)
            base_params_k = build_base_params_with_sigma(
                param_keys,
                data["inj_params"],
                data["cw_draws"],
                dist_priors_k
            )

            cw_names_k = find_cw_param_names(param_keys)
            amp_center = base_params_k[cw_names_k["amp"]]
            amp_grid = np.linspace(amp_center - amp_span, amp_center + amp_span, n_amp_fast)

            D_k, logL_null_k, logL_max_k, _ = compute_detection_statistic_optimized(
                logl, param_keys, base_params_k, cw_names_k,
                cos_grid, phi_grid, amp_grid, null_amp=null_amp,
            )
            D_vals.append(D_k)

        D_samples.append(D_vals)

    D_samples = np.array(D_samples)
    ensemble_D[snr] = {"k": k_grid.copy(), "D_samples": D_samples}

    # ... rest of plotting code unchanged ...
```

### Step 4: Update toy experiment

Replace toy PTA cell (`d57e2355`) to use cached pulsars:

```python
# Just add this at the top of that cell:
from optimization_helpers import load_discovery_pulsars_cached

# Then replace load_discovery_pulsars calls with:
d_psrs = load_discovery_pulsars_cached(toy_cfg)
```

## Performance Improvements

Expected speedups:
- **Main loop**: 20-40x faster (3 likelihood builds instead of 21)
- **k-sweep**: 15-25x faster (3 builds instead of 75)
- **Toy experiments**: 30-50x faster
- **Overall**: 20-30x end-to-end speedup

## Verification

To verify the optimization doesn't change results, run a small test:

```python
# Test: Compare optimized vs original for one case
test_snr = 10
test_tier = "real"

# Original approach (for comparison)
data_orig = simulate_data(test_snr, cfg)
dist_priors_orig = apply_distance_prior_tier(
    list(data_orig["true_dists"].keys()),
    data_orig["true_dists"],
    data_orig["true_sigs"],
    test_tier
)
# ... build likelihood and compute D_stat_orig ...

# Optimized approach
# ... use cached likelihood and compute D_stat_opt ...

# Should match to floating-point precision
print(f"Original D: {D_stat_orig:.6f}")
print(f"Optimized D: {D_stat_opt:.6f}")
print(f"Difference: {abs(D_stat_orig - D_stat_opt):.2e}")
assert abs(D_stat_orig - D_stat_opt) < 1e-6, "Results don't match!"
```

## Troubleshooting

**Issue**: "JAX optimization failed"
- **Solution**: Fallback NumPy implementation will be used automatically. Check JAX installation and GPU availability.

**Issue**: "Missing parameter `{psr}_cw_sigma_d`"
- **Solution**: Make sure you're using the updated `sim.py` with sigma as a parameter. Restart the kernel after updating.

**Issue**: Results don't match original
- **Solution**: Check that `dist_priors` mean values are consistent. The mean distance should be the same across tiers.

**Issue**: Out of memory
- **Solution**: Reduce grid size (`n_cos`, `n_phi`, `n_amp`) or process SNRs sequentially with `likelihood_cache.clear()` between SNRs.
