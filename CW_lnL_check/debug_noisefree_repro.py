"""Minimal noise-free CW distance likelihood reproduction.

Runs a B1953+29 distance scan with one injected CW and progressively adds
likelihood components to identify which one moves the maximum away from truth.
"""

import argparse
import time

import jax
import jax.numpy as jnp
import numpy as np

import discovery as ds

from cw_helpers import (
    MultiSourceDelay,
    build_noisefree_likelihood,
    inject_noisefree_cw,
    load_pulsars,
)

jax.config.update("jax_enable_x64", True)


CW_PARAMS = [{
    "cos_gwtheta": 0.3,
    "gwphi": 2.5,
    "cos_inc": -0.2,
    "log10_mc": 9.0,
    "log10_fgw": -8.0,
    "log10_h": -12.0,
    "phase0": 1.0,
    "psi": 0.7,
}]


def build_variant_likelihood(
    disco_psrs,
    residual_map,
    cw_params_list,
    *,
    include_timing=False,
    include_red=False,
    include_gwb=False,
    components=30,
    timing_variance=1e-14,
    log10_equad=-8.0,
):
    """Build discovery likelihood with toggled components."""
    num_cw = len(cw_params_list)
    noisedict = {}
    for psr in disco_psrs:
        noisedict[f"{psr.name}_KAT_MKBF_efac"] = 1.0
        noisedict[f"{psr.name}_KAT_MKBF_log10_ecorr"] = -8.0
        noisedict[f"{psr.name}_KAT_MKBF_log10_t2equad"] = log10_equad

    noise_terms = {
        psr.name: ds.makenoise_measurement(psr, noisedict=noisedict)
        for psr in disco_psrs
    }
    timing_terms = {
        psr.name: ds.makegp_timing(psr, variance=timing_variance)
        for psr in disco_psrs
    }

    cgp = None
    ggp = None
    if include_red or include_gwb:
        span = ds.getspan(disco_psrs)
        cgp = ds.makecommongp_fourier(
            disco_psrs, ds.powerlaw, components, span, name="rednoise"
        )
        if include_gwb:
            ggp = ds.makeglobalgp_fourier(
                disco_psrs, ds.powerlaw, ds.hd_orf, components, span, name="gwb"
            )

    pulsar_likes = []
    for psr in disco_psrs:
        args = [
            np.array(residual_map[psr.name], copy=True),
            noise_terms[psr.name],
        ]
        if include_timing:
            args.append(timing_terms[psr.name])
        args.append(MultiSourceDelay(psr, num_cw, include_pterm=True))
        pulsar_likes.append(ds.PulsarLikelihood(args))

    fml = ds.ArrayLikelihood(pulsar_likes, commongp=cgp, globalgp=ggp)
    logl = fml.logL

    suffixes = ["" if i == 0 else f"_{i+1}" for i in range(num_cw)]
    param_dict = {}
    for suffix, cw_p in zip(suffixes, cw_params_list):
        for key in MultiSourceDelay.global_params:
            param_dict[f"cw_{key}{suffix}"] = cw_p[key]
    for psr in disco_psrs:
        param_dict[f"{psr.name}_cw_p_dist"] = psr.pdist[0]
        param_dict[f"{psr.name}_rednoise_log10_A"] = -20.0
        param_dict[f"{psr.name}_rednoise_gamma"] = 4.0
    param_dict["gwb_log10_A"] = -20.0
    param_dict["gwb_gamma"] = 4.333

    order = {k: i for i, k in enumerate(logl.params)}
    keys = [
        key for key, _ in sorted(param_dict.items(), key=lambda kv: order.get(kv[0], 10**9))
        if key in logl.params
    ]
    base = jnp.array([param_dict[k] for k in keys], dtype=jnp.float64)

    def wrapped(x_array):
        return logl({k: v for k, v in zip(keys, x_array)})

    return wrapped, keys, base


def scan_distance(logl_fn, base_values, param_keys, dist_key, center, width, n_points):
    idx = param_keys.index(dist_key)
    grid = np.linspace(max(1e-6, center - width), center + width, n_points)
    grid = np.unique(np.concatenate([grid, [center]]))
    values = []
    for dist in grid:
        x = base_values.at[idx].set(dist)
        values.append(float(logl_fn(x)))
    return grid, np.array(values)


def summarize_variant(name, builder, disco_psrs, residual_map, true_dist, dist_key, args):
    t0 = time.time()
    logl_fn, param_keys, base_values = builder()
    truth_ll = float(logl_fn(base_values))
    grid, vals = scan_distance(
        logl_fn,
        base_values,
        param_keys,
        dist_key,
        true_dist,
        args.width,
        args.points,
    )
    imax = int(np.argmax(vals))
    truth_idx = int(np.where(grid == true_dist)[0][0])
    print(f"\n[{name}]")
    print(f"  build+scan_s: {time.time() - t0:.2f}")
    print(f"  params: {len(param_keys)}")
    print(f"  lnL_truth: {truth_ll:.12g}")
    print(f"  scan_truth_lnL: {vals[truth_idx]:.12g}")
    print(f"  scan_max_lnL: {vals[imax]:.12g}")
    print(f"  true_dist_kpc: {true_dist:.12g}")
    print(f"  max_dist_kpc: {grid[imax]:.12g}")
    print(f"  max_minus_truth_kpc: {grid[imax] - true_dist:.12g}")
    print(f"  lnL_truth_minus_max: {vals[truth_idx] - vals[imax]:.12g}")
    print(f"  truth_is_global_max: {truth_idx == imax or np.isclose(vals[truth_idx], vals[imax])}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=101)
    parser.add_argument("--width", type=float, default=1.0)
    parser.add_argument("--components", type=int, default=30)
    parser.add_argument("--psr", default="B1953+29")
    args = parser.parse_args()

    ent_psrs, disco_psrs = load_pulsars(5)
    names = [p.name for p in disco_psrs]
    if args.psr not in names:
        raise ValueError(f"{args.psr} not in first 5 discovery pulsars: {names}")

    residual_map, true_dists = inject_noisefree_cw(disco_psrs, CW_PARAMS)
    max_template_error = 0.0
    check_map, _ = inject_noisefree_cw(disco_psrs, CW_PARAMS)
    for psr in disco_psrs:
        max_template_error = max(
            max_template_error,
            float(np.max(np.abs(residual_map[psr.name] - check_map[psr.name]))),
        )

    target = next(p for p in disco_psrs if p.name == args.psr)
    dist_key = f"{target.name}_cw_p_dist"
    true_dist = true_dists[target.name]

    print(f"loaded_pulsars: {names}")
    print(f"target: {target.name}")
    print(f"max_template_minus_injection: {max_template_error:.12g}")
    print(f"scan_width_kpc: {args.width}")
    print(f"scan_points_plus_truth: {args.points}+truth")

    variants = [
        (
            "cw_only",
            lambda: build_variant_likelihood(
                disco_psrs, residual_map, CW_PARAMS, components=args.components
            ),
        ),
        (
            "cw_plus_timing",
            lambda: build_variant_likelihood(
                disco_psrs,
                residual_map,
                CW_PARAMS,
                include_timing=True,
                components=args.components,
            ),
        ),
        (
            "cw_plus_timing_plus_red",
            lambda: build_variant_likelihood(
                disco_psrs,
                residual_map,
                CW_PARAMS,
                include_timing=True,
                include_red=True,
                components=args.components,
            ),
        ),
        (
            "build_noisefree_likelihood",
            lambda: build_noisefree_likelihood(
                disco_psrs,
                residual_map,
                1,
                CW_PARAMS,
                include_gwb=True,
            ),
        ),
    ]
    for name, builder in variants:
        summarize_variant(name, builder, disco_psrs, residual_map, true_dist, dist_key, args)


if __name__ == "__main__":
    main()
