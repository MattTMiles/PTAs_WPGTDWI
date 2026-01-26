"""Optimization helpers for Q1 detection statistic notebook.

These functions enable caching the likelihood object and reusing it across
different distance prior tiers by only updating the sigma parameters.
"""

from pathlib import Path
import copy
import jax
import jax.numpy as jnp
import numpy as np
import discovery as ds


# Global cache for pulsars to avoid repeated disk I/O
_PULSAR_CACHE = {}


def load_discovery_pulsars_cached(cfg):
    """Load discovery pulsars with caching to avoid repeated disk reads.

    Returns deep copies to avoid mutation issues across different uses.
    """
    cache_key = (cfg.n_pulsars, str(cfg.feather_dir))

    if cache_key not in _PULSAR_CACHE:
        files = sorted(Path(cfg.feather_dir).glob("*.feather"))
        _PULSAR_CACHE[cache_key] = [ds.Pulsar.read_feather(f) for f in files[:cfg.n_pulsars]]

    # Return deep copies to avoid mutation issues
    return [copy.deepcopy(psr) for psr in _PULSAR_CACHE[cache_key]]


def build_base_params_with_sigma(param_keys, inj_params, cw_draws, dist_priors):
    """Build parameter dict including sigma values for each pulsar.

    This extends the standard build_base_params to include the distance
    uncertainty (sigma) as parameters, allowing tier changes without
    rebuilding the likelihood.

    Args:
        param_keys: List of parameter names
        inj_params: Injection parameters dict
        cw_draws: CW parameter draws
        dist_priors: Dict mapping pulsar name to (mean, sigma) tuple

    Returns:
        Complete parameter dict including sigma values
    """
    # Build base params (CW and noise parameters)
    params = {k: inj_params.get(k, 0.0) for k in param_keys}

    # Add CW draw parameters
    for draw in cw_draws:
        suffix = draw.get("suffix", "")
        for field in ["cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi"]:
            name = f"cw_{field}{suffix}"
            if name in params:
                params[name] = draw[field]

    # Add sigma values for each pulsar
    for key in param_keys:
        if key.endswith("_cw_sigma_d"):
            psr_name = key.replace("_cw_sigma_d", "")
            if psr_name in dist_priors:
                _, sigma = dist_priors[psr_name]
                params[key] = sigma

    return params


def update_sigma_params(base_params, param_keys, dist_priors):
    """Update only the sigma parameters in base_params dict.

    This is more efficient than rebuilding the entire params dict.

    Args:
        base_params: Existing parameter dict (will be modified in place)
        param_keys: List of parameter names
        dist_priors: Dict mapping pulsar name to (mean, sigma) tuple

    Returns:
        Updated base_params (same object, modified in place)
    """
    for key in param_keys:
        if key.endswith("_cw_sigma_d"):
            psr_name = key.replace("_cw_sigma_d", "")
            if psr_name in dist_priors:
                _, sigma = dist_priors[psr_name]
                base_params[key] = sigma
    return base_params


def compute_detection_statistic_optimized(logl, param_keys, base_params, cw_names,
                                         cos_grid, phi_grid, amp_grid, null_amp):
    """Optimized detection statistic computation using flattened JAX vmap.

    This version uses a single vmap over flattened grid points instead of
    nested vmaps, which is more efficient for JAX compilation and execution.

    Args:
        logl: Log-likelihood callable
        param_keys: List of parameter names
        base_params: Base parameters dict
        cw_names: Dict with keys 'cos', 'phi', 'amp' mapping to parameter names
        cos_grid: Array of cos(theta) values
        phi_grid: Array of phi values
        amp_grid: Array of log10_h values
        null_amp: Amplitude value for null hypothesis

    Returns:
        Tuple of (D_stat, logL_null, logL_max, logL_grid_2d)
    """
    try:
        base_vec = jnp.array([base_params[k] for k in param_keys])
        idx_cos = param_keys.index(cw_names["cos"])
        idx_phi = param_keys.index(cw_names["phi"])
        idx_amp = param_keys.index(cw_names["amp"])

        # Create flattened grid of (cos, phi, amp) points
        cos_mesh, phi_mesh, amp_mesh = jnp.meshgrid(
            jnp.array(cos_grid),
            jnp.array(phi_grid),
            jnp.array(amp_grid),
            indexing='ij'
        )
        grid_points = jnp.stack([
            cos_mesh.ravel(),
            phi_mesh.ravel(),
            amp_mesh.ravel()
        ], axis=1)

        @jax.jit
        def eval_point(point):
            """Evaluate logL at a single (cos, phi, amp) grid point."""
            x = base_vec.at[idx_cos].set(point[0])
            x = x.at[idx_phi].set(point[1])
            x = x.at[idx_amp].set(point[2])
            params = {k: x[i] for i, k in enumerate(param_keys)}
            return logl(params)

        # Single vmap over all grid points (more efficient than nested vmaps)
        eval_grid = jax.vmap(eval_point)
        logL_flat = eval_grid(grid_points)

        # Reshape and take max over amplitude dimension
        logL_grid_3d = logL_flat.reshape(len(cos_grid), len(phi_grid), len(amp_grid))
        logL_grid_2d = logL_grid_3d.max(axis=2)  # Max over amplitude
        logL_max = float(logL_grid_2d.max())

        # Null hypothesis evaluation
        x_null = base_vec.at[idx_amp].set(null_amp)
        params_null = {k: x_null[i] for i, k in enumerate(param_keys)}
        logL_null = float(logl(params_null))

        d_stat = 2 * (logL_max - logL_null)
        return d_stat, logL_null, logL_max, np.array(logL_grid_2d)

    except Exception as e:
        # Fallback to NumPy implementation if JAX fails
        print(f"JAX optimization failed ({e}), using fallback NumPy implementation")
        logL_grid = np.zeros((len(cos_grid), len(phi_grid)))
        for i, cosv in enumerate(cos_grid):
            for j, phiv in enumerate(phi_grid):
                best = -np.inf
                for ampv in amp_grid:
                    params = dict(base_params)
                    params[cw_names["cos"]] = float(cosv)
                    params[cw_names["phi"]] = float(phiv)
                    params[cw_names["amp"]] = float(ampv)
                    lp = float(logl(params))
                    if lp > best:
                        best = lp
                logL_grid[i, j] = best
        logL_max = float(logL_grid.max())
        params_null = dict(base_params)
        params_null[cw_names["amp"]] = null_amp
        logL_null = float(logl(params_null))
        d_stat = 2 * (logL_max - logL_null)
        return d_stat, logL_null, logL_max, logL_grid


class LikelihoodCache:
    """Cache for likelihood objects to avoid rebuilding for each tier.

    This class manages cached likelihoods per SNR, allowing sigma parameters
    to be updated without rebuilding the entire likelihood structure.
    """

    def __init__(self):
        self._cache = {}

    def build_or_get(self, snr, cfg, simulate_data_fn, build_likelihood_fn):
        """Get cached likelihood or build a new one.

        Args:
            snr: SNR value (cache key)
            cfg: Q1Config object
            simulate_data_fn: Function to simulate data (called only on cache miss)
            build_likelihood_fn: Function to build likelihood (called only on cache miss)

        Returns:
            Dict with keys: 'logl', 'param_keys', 'data', 'd_psrs'
        """
        if snr not in self._cache:
            # Cache miss - build new likelihood
            data = simulate_data_fn(snr, cfg)

            # Use cached pulsar loading
            d_psrs = load_discovery_pulsars_cached(cfg)
            if hasattr(cfg, 'use_small_pta') and cfg.use_small_pta:
                d_psrs = [psr for psr in d_psrs if psr.name in data["residual_map"]]

            # Build likelihood with initial dist_priors
            # The sigma values will be updated per-tier without rebuilding
            initial_dist_priors = {
                psr.name: (data["true_dists"][psr.name], data["true_sigs"].get(psr.name, 0.1))
                for psr in d_psrs
            }

            logl, param_keys = build_likelihood_fn(
                d_psrs,
                data,
                initial_dist_priors
            )

            self._cache[snr] = {
                "logl": logl,
                "param_keys": param_keys,
                "data": data,
                "d_psrs": d_psrs,
            }

        return self._cache[snr]

    def clear(self):
        """Clear the cache."""
        self._cache.clear()
