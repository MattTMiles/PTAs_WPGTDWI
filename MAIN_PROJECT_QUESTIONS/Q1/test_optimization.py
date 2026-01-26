"""Test script to verify optimization preserves correctness.

Run this after implementing the optimizations to ensure results match.
"""

import sys
from pathlib import Path
import numpy as np

# Add Q1 directory to path
sys.path.insert(0, str(Path(__file__).parent))

from sim import (
    Q1Config,
    DistanceTier,
    draw_cw_population,
    load_enterprise_pulsars,
    load_discovery_pulsars,
    build_enterprise_pta,
    build_discovery_loglike_marginalized,
    draw_noise_params,
    inject_cw_params,
    simulate_enterprise_residuals,
    set_toaerrs,
    make_rng,
)

from optimization_helpers import (
    load_discovery_pulsars_cached,
    build_base_params_with_sigma,
    update_sigma_params,
)


def test_sigma_parameter_equivalence():
    """Test that sigma-as-parameter gives same results as closure capture."""
    print("Testing sigma parameter equivalence...")

    cfg = Q1Config(
        n_pulsars=5,  # Small for fast test
        scenario="A",
        seed=42,
        distance_tier=DistanceTier("real", 1.0),
        include_pterm=True,
    )

    # Simulate data
    rng = make_rng(cfg.seed)
    cw_draws = draw_cw_population(cfg, rng)
    cw_draws[0]["log10_fgw"] = -8.5

    e_psrs = load_enterprise_pulsars(cfg)
    e_psrs = e_psrs[:cfg.n_pulsars]
    set_toaerrs(e_psrs, sigma=1e-6)

    true_dists = {p.name: float(p.pdist[0]) for p in e_psrs}
    true_sigs = {p.name: float(p.pdist[1]) if len(np.atleast_1d(p.pdist)) > 1 else 0.1 for p in e_psrs}

    pta, ring_vecs = build_enterprise_pta(e_psrs, cfg, cw_draws)
    pta.set_default_params({})
    noise_params = draw_noise_params(pta, cfg, rng, psrs=e_psrs)
    inj_params = inject_cw_params(noise_params, cw_draws, cfg)
    sim_resids = simulate_enterprise_residuals(pta, inj_params, sparse_cholesky=True)
    name_to_resid = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}

    d_psrs = load_discovery_pulsars(cfg)
    d_psrs = d_psrs[:cfg.n_pulsars]
    residual_map = {psr.name: name_to_resid[psr.name] for psr in d_psrs}

    # Test 1: Build likelihood with sigma=0.1 for all pulsars
    print("  Test 1: Building likelihood with initial sigma=0.1...")
    dist_priors_1 = {psr.name: (true_dists[psr.name], 0.1) for psr in d_psrs}
    logl_1, param_keys_1 = build_discovery_loglike_marginalized(
        d_psrs, cw_draws, cfg, dist_priors_1,
        ring_vecs=ring_vecs, residual_map=residual_map, include_cw=True,
    )

    base_params_1 = build_base_params_with_sigma(param_keys_1, inj_params, cw_draws, dist_priors_1)
    logL_1 = float(logl_1(base_params_1))
    print(f"    logL with sigma=0.1: {logL_1:.6f}")

    # Test 2: Update sigma to 0.5 using parameter approach
    print("  Test 2: Updating sigma to 0.5 via parameters...")
    dist_priors_2 = {psr.name: (true_dists[psr.name], 0.5) for psr in d_psrs}
    base_params_2 = update_sigma_params(base_params_1.copy(), param_keys_1, dist_priors_2)
    logL_2 = float(logl_1(base_params_2))  # Reuse same likelihood object!
    print(f"    logL with sigma=0.5: {logL_2:.6f}")

    # Test 3: Build new likelihood with sigma=0.5 to verify
    print("  Test 3: Building fresh likelihood with sigma=0.5 for verification...")
    logl_3, param_keys_3 = build_discovery_loglike_marginalized(
        d_psrs, cw_draws, cfg, dist_priors_2,
        ring_vecs=ring_vecs, residual_map=residual_map, include_cw=True,
    )
    base_params_3 = build_base_params_with_sigma(param_keys_3, inj_params, cw_draws, dist_priors_2)
    logL_3 = float(logl_3(base_params_3))
    print(f"    logL with sigma=0.5 (fresh build): {logL_3:.6f}")

    # Verify results match
    diff = abs(logL_2 - logL_3)
    print(f"\n  Difference between updated and fresh: {diff:.2e}")

    if diff < 1e-8:
        print("  ✓ PASS: Sigma update matches fresh build!")
        return True
    else:
        print(f"  ✗ FAIL: Difference {diff:.2e} exceeds tolerance 1e-8")
        return False


def test_pulsar_caching():
    """Test that cached pulsar loading works correctly."""
    print("\nTesting pulsar caching...")

    cfg = Q1Config(n_pulsars=10, seed=42)

    # First load
    print("  First load (should hit disk)...")
    psrs_1 = load_discovery_pulsars_cached(cfg)
    n_psrs_1 = len(psrs_1)

    # Second load (should hit cache)
    print("  Second load (should hit cache)...")
    psrs_2 = load_discovery_pulsars_cached(cfg)
    n_psrs_2 = len(psrs_2)

    # Verify counts match
    if n_psrs_1 == n_psrs_2 == cfg.n_pulsars:
        print(f"  ✓ PASS: Loaded {n_psrs_1} pulsars consistently")
        return True
    else:
        print(f"  ✗ FAIL: Inconsistent counts: {n_psrs_1} vs {n_psrs_2}")
        return False


def test_parameter_updates():
    """Test that parameter updates work correctly."""
    print("\nTesting parameter updates...")

    param_keys = ["cw_log10_h", "cw_gwphi", "J0030+0451_cw_sigma_d", "J0437-4715_cw_sigma_d"]
    inj_params = {"cw_log10_h": -12.0, "cw_gwphi": 1.5}
    cw_draws = [{"log10_h": -12.0, "gwphi": 1.5, "suffix": ""}]
    dist_priors = {"J0030+0451": (0.3, 0.05), "J0437-4715": (0.15, 0.01)}

    params = build_base_params_with_sigma(param_keys, inj_params, cw_draws, dist_priors)

    # Check sigma values
    sigma_1 = params.get("J0030+0451_cw_sigma_d")
    sigma_2 = params.get("J0437-4715_cw_sigma_d")

    print(f"  J0030+0451 sigma: {sigma_1}")
    print(f"  J0437-4715 sigma: {sigma_2}")

    if sigma_1 == 0.05 and sigma_2 == 0.01:
        print("  ✓ PASS: Sigma values set correctly")
        return True
    else:
        print(f"  ✗ FAIL: Expected (0.05, 0.01), got ({sigma_1}, {sigma_2})")
        return False


if __name__ == "__main__":
    print("=" * 60)
    print("Q1 Optimization Test Suite")
    print("=" * 60)

    results = []

    try:
        results.append(("Sigma parameter equivalence", test_sigma_parameter_equivalence()))
    except Exception as e:
        print(f"  ✗ FAIL: Exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Sigma parameter equivalence", False))

    try:
        results.append(("Pulsar caching", test_pulsar_caching()))
    except Exception as e:
        print(f"  ✗ FAIL: Exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Pulsar caching", False))

    try:
        results.append(("Parameter updates", test_parameter_updates()))
    except Exception as e:
        print(f"  ✗ FAIL: Exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Parameter updates", False))

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {name}")

    all_passed = all(passed for _, passed in results)
    print("=" * 60)
    if all_passed:
        print("All tests passed! ✓")
        sys.exit(0)
    else:
        print("Some tests failed! ✗")
        sys.exit(1)
