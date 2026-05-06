"""Deterministic multi-CW pulsar-distance recovery by pairwise coordinate ascent.

This is the production-facing counterpart to the 1CW evolution notebook:

- inject N_CW known continuous waves into noise-free PTA residuals
- keep CW parameters fixed at truth
- optimize only pulsar distances
- initialize chains at EM prior means or random prior draws
- hold non-scanned distances at the chain's current best guess
- update distances via local 2D brute-force scans over pulsar pairs

Truth is used only for scoring/diagnostics, not during updates.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "../CW_lnL_check")

from cw_helpers import (  # noqa: E402
    build_disco_likelihood,
    build_enterprise_pta,
    compute_mode_spacing,
    build_pure_cw_likelihood,
    generate_injection_params,
    inject_noisefree_cw,
    load_pulsars,
    scan_pair_window_2d,
    simulate,
)


OUTDIR = Path("02_multicw_pairwise_outputs")

CW_LIBRARY = [
    dict(
        cos_gwtheta=0.30,
        gwphi=2.50,
        cos_inc=-0.20,
        phase0=1.00,
        psi=0.70,
        log10_h=-12.0,
        log10_mc=9.00,
        log10_fgw=-8.00,
    ),
    dict(
        cos_gwtheta=-0.50,
        gwphi=0.80,
        cos_inc=0.40,
        phase0=2.10,
        psi=1.30,
        log10_h=-12.0,
        log10_mc=9.00,
        log10_fgw=-7.80,
    ),
    dict(
        cos_gwtheta=0.05,
        gwphi=4.20,
        cos_inc=-0.65,
        phase0=0.35,
        psi=2.30,
        log10_h=-12.0,
        log10_mc=9.00,
        log10_fgw=-8.20,
    ),
    dict(
        cos_gwtheta=0.70,
        gwphi=5.40,
        cos_inc=0.10,
        phase0=1.70,
        psi=0.20,
        log10_h=-12.0,
        log10_mc=9.00,
        log10_fgw=-7.90,
    ),
    dict(
        cos_gwtheta=-0.10,
        gwphi=3.30,
        cos_inc=0.80,
        phase0=2.90,
        psi=1.80,
        log10_h=-12.0,
        log10_mc=9.00,
        log10_fgw=-8.10,
    ),
]


def select_pulsars(n_psr: int):
    """Load all pulsars, sort by EM distance prior width, keep tightest n_psr."""
    ent_all, disco_all = load_pulsars(None)
    pairs = sorted(zip(ent_all, disco_all), key=lambda pair: pair[1].pdist[1])
    selected = pairs[:n_psr]
    ent_psrs = [pair[0] for pair in selected]
    disco_psrs = [pair[1] for pair in selected]
    return ent_psrs, disco_psrs, len(disco_all)


def make_pair_schedule(n_psr: int, mode: str):
    """Return ordered pairs of selected pulsar indices."""
    if mode == "all":
        return [(i, j) for i in range(n_psr) for j in range(i + 1, n_psr)]
    if mode == "anchor":
        return [(0, j) for j in range(1, n_psr)]
    if mode == "anchor_ladder":
        pairs = [(0, j) for j in range(1, n_psr)]
        pairs += [(j, j + 1) for j in range(1, n_psr - 1)]
        return list(dict.fromkeys(pairs))
    raise ValueError(f"unknown pair schedule mode: {mode}")


def set_distance_values(values, param_keys, disco_psrs, distances):
    """Return copy of values with all pulsar-distance parameters replaced."""
    out = np.array(values, dtype=float).copy()
    for psr, dist in zip(disco_psrs, distances):
        key = f"{psr.name}_cw_p_dist"
        if key in param_keys:
            out[param_keys.index(key)] = dist
    return out


def get_distance_values(values, param_keys, disco_psrs):
    vals = []
    for psr in disco_psrs:
        key = f"{psr.name}_cw_p_dist"
        vals.append(float(values[param_keys.index(key)]))
    return np.asarray(vals)


def draw_initial_distances(disco_psrs, rng, init_kind, prior_sigma_clip):
    means = np.array([psr.pdist[0] for psr in disco_psrs], dtype=float)
    sigmas = np.array([psr.pdist[1] for psr in disco_psrs], dtype=float)
    if init_kind == "mean":
        return means
    if init_kind == "random":
        lo = np.maximum(0.01, means - prior_sigma_clip * sigmas)
        hi = means + prior_sigma_clip * sigmas
        return rng.uniform(lo, hi)
    raise ValueError(f"unknown init kind: {init_kind}")


def min_mode_spacings(disco_psrs, cw_params_list):
    spacings = []
    for psr in disco_psrs:
        vals = [
            compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"], cw["log10_fgw"], psr.pos)
            for cw in cw_params_list
        ]
        spacings.append(float(np.nanmin(vals)))
    return np.asarray(spacings)


def apply_cw_overrides(cw_params_list, args):
    out = [dict(cw) for cw in cw_params_list]
    for cw in out:
        if args.log10_h is not None:
            cw["log10_h"] = args.log10_h
        if args.log10_mc is not None:
            cw["log10_mc"] = args.log10_mc
    return out


def cw_list_from_enterprise_params(enterprise_params, cw_block_names):
    keys = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi")
    out = []
    for name in cw_block_names:
        out.append({key: float(enterprise_params[f"{name}_{key}"]) for key in keys})
    return out


def add_phenomenological_noise(residual_map, disco_psrs, args):
    """Add simple residual-domain disturbances for robustness sweeps.

    This is not a full stochastic-process likelihood.  It is a controlled
    stress test: white residual noise, independent red-ish low-frequency
    sinusoids, and common red-ish sinusoids across all pulsars.
    """
    if not (args.extra_white_rms or args.extra_red_rms or args.extra_common_red_rms):
        return residual_map

    rng = np.random.default_rng(args.noise_seed)
    noisy = {name: np.array(vals, copy=True) for name, vals in residual_map.items()}

    common_cache = {}
    for psr in disco_psrs:
        vals = noisy[psr.name]
        toas = np.asarray(psr.toas, dtype=float)
        t = toas - np.mean(toas)
        span = np.ptp(toas)
        if args.extra_white_rms:
            vals = vals + rng.normal(0.0, args.extra_white_rms, size=vals.shape)
        if args.extra_red_rms:
            vals = vals + red_like_series(t, span, args.extra_red_rms, rng)
        if args.extra_common_red_rms:
            # Same random Fourier coefficients, evaluated at each pulsar's TOAs.
            if "coeffs" not in common_cache:
                common_cache["coeffs"] = make_red_coeffs(args.extra_common_red_rms, rng)
            vals = vals + red_like_series(t, span, args.extra_common_red_rms, rng, common_cache["coeffs"])
        noisy[psr.name] = vals
    return noisy


def make_red_coeffs(rms, rng, n_terms=6):
    amps = rng.normal(size=(2, n_terms))
    weights = np.arange(1, n_terms + 1, dtype=float) ** -2.0
    amps = amps * weights[None, :]
    norm = np.sqrt(np.mean(amps**2))
    return amps * (rms / norm if norm > 0 else 0.0)


def red_like_series(t, span, rms, rng, coeffs=None, n_terms=6):
    if span <= 0:
        return np.zeros_like(t)
    if coeffs is None:
        coeffs = make_red_coeffs(rms, rng, n_terms=n_terms)
    x = np.zeros_like(t, dtype=float)
    for idx in range(coeffs.shape[1]):
        freq = (idx + 1) / span
        phase = 2.0 * np.pi * freq * t
        x += coeffs[0, idx] * np.sin(phase) + coeffs[1, idx] * np.cos(phase)
    return x


def score(values, param_keys, disco_psrs, truth_dist, mode_spacings):
    recovered = get_distance_values(values, param_keys, disco_psrs)
    err = recovered - truth_dist
    sigmas = np.array([psr.pdist[1] for psr in disco_psrs], dtype=float)
    mode_err = err / mode_spacings
    return {
        "max_abs_kpc": float(np.nanmax(np.abs(err))),
        "median_abs_kpc": float(np.nanmedian(np.abs(err))),
        "max_abs_sigma": float(np.nanmax(np.abs(err / sigmas))),
        "median_abs_sigma": float(np.nanmedian(np.abs(err / sigmas))),
        "max_abs_modes": float(np.nanmax(np.abs(mode_err))),
        "median_abs_modes": float(np.nanmedian(np.abs(mode_err))),
        "frac_within_half_mode": float(np.nanmean(np.abs(mode_err) < 0.5)),
        "frac_within_one_mode": float(np.nanmean(np.abs(mode_err) < 1.0)),
    }


def run_chain(
    chain_id,
    init_kind,
    logl_fn,
    param_keys,
    base_values,
    disco_psrs,
    cw_params_list,
    truth_dist,
    mode_spacings,
    pair_schedule,
    rng,
    args,
):
    current_dist = draw_initial_distances(
        disco_psrs, rng, init_kind, args.prior_sigma_clip
    )
    current_values = set_distance_values(
        base_values, param_keys, disco_psrs, current_dist
    )
    history = []
    start_lnL = float(logl_fn(current_values))
    print(f"chain {chain_id:02d} init={init_kind:6s} lnL={start_lnL:.3f}")

    stale_sweeps = 0
    last_sweep_lnL = start_lnL
    for sweep in range(args.max_sweeps):
        sweep_start = time.time()
        before = float(logl_fn(current_values))
        n_updates = 0
        pairs_this_sweep = pair_schedule
        if args.max_pairs_per_sweep:
            pairs_this_sweep = pair_schedule[: args.max_pairs_per_sweep]

        for pair_idx, (i, j) in enumerate(pairs_this_sweep):
            psr_i = disco_psrs[i]
            psr_j = disco_psrs[j]
            current_dist = get_distance_values(current_values, param_keys, disco_psrs)
            res = scan_pair_window_2d(
                logl_fn,
                current_values,
                param_keys,
                psr_i,
                psr_j,
                cw_params_list,
                center_i=current_dist[i],
                center_j=current_dist[j],
                half_width_modes=args.half_width_modes,
                points_per_mode=args.points_per_mode,
                min_points=args.min_points,
                chunk_size=args.chunk_size,
                include_truth=False,
            )
            new_values = np.array(current_values, dtype=float).copy()
            new_values[param_keys.index(f"{psr_i.name}_cw_p_dist")] = res["peak_i"]
            new_values[param_keys.index(f"{psr_j.name}_cw_p_dist")] = res["peak_j"]
            new_lnL = float(res["peak_lnL"])
            old_lnL = float(logl_fn(current_values))

            # Coordinate ascent should be monotone. A tiny tolerance avoids
            # rejecting equivalent modes due to numerical roundoff.
            if new_lnL >= old_lnL - args.accept_tol:
                current_values = new_values
                n_updates += 1

            if args.progress_every and (pair_idx + 1) % args.progress_every == 0:
                now = float(logl_fn(current_values))
                print(
                    f"  chain {chain_id:02d} sweep {sweep:02d} "
                    f"pair {pair_idx + 1:04d}/{len(pairs_this_sweep):04d} "
                    f"lnL={now:.3f}"
                )

        after = float(logl_fn(current_values))
        delta = after - before
        summary = score(current_values, param_keys, disco_psrs, truth_dist, mode_spacings)
        row = {
            "sweep": sweep,
            "lnL": after,
            "delta_lnL": delta,
            "n_updates": n_updates,
            "seconds": time.time() - sweep_start,
            **summary,
        }
        history.append(row)
        print(
            f"chain {chain_id:02d} sweep {sweep:02d} "
            f"lnL={after:.3f} dlnL={delta:.3g} "
            f"within0.5mode={summary['frac_within_half_mode']:.2%} "
            f"time={row['seconds']:.1f}s"
        )

        if abs(after - last_sweep_lnL) < args.convergence_tol:
            stale_sweeps += 1
        else:
            stale_sweeps = 0
        last_sweep_lnL = after
        if stale_sweeps >= args.convergence_sweeps:
            break

    recovered = get_distance_values(current_values, param_keys, disco_psrs)
    return {
        "chain_id": chain_id,
        "init_kind": init_kind,
        "history": history,
        "final_lnL": float(logl_fn(current_values)),
        "final_score": score(current_values, param_keys, disco_psrs, truth_dist, mode_spacings),
        "recovered_distances": recovered.tolist(),
        "initial_lnL": start_lnL,
    }


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--n-psr", type=int, default=80)
    p.add_argument("--n-cw", type=int, default=2)
    p.add_argument("--n-chains", type=int, default=4)
    p.add_argument("--max-sweeps", type=int, default=2)
    p.add_argument("--data-mode", choices=["pure", "stochastic"], default="pure")
    p.add_argument("--stochastic-scenario", default="well_separated",
                   choices=["single", "close_freq", "close_sky", "well_separated", "realistic"])
    p.add_argument("--components", type=int, default=30)
    p.add_argument("--gwb-log10-a", type=float, default=-17.5)
    p.add_argument("--gwb-gamma", type=float, default=4.333)
    p.add_argument("--no-gwb", action="store_true")
    p.add_argument("--log10-h", type=float, default=None)
    p.add_argument("--log10-mc", type=float, default=None)
    p.add_argument("--extra-white-rms", type=float, default=0.0)
    p.add_argument("--extra-red-rms", type=float, default=0.0)
    p.add_argument("--extra-common-red-rms", type=float, default=0.0)
    p.add_argument("--noise-seed", type=int, default=24680)
    p.add_argument("--pair-mode", choices=["anchor", "anchor_ladder", "all"], default="anchor_ladder")
    p.add_argument("--max-pairs-per-sweep", type=int, default=0)
    p.add_argument("--half-width-modes", type=float, default=3.0)
    p.add_argument("--points-per-mode", type=float, default=6.0)
    p.add_argument("--min-points", type=int, default=41)
    p.add_argument("--chunk-size", type=int, default=0)
    p.add_argument("--prior-sigma-clip", type=float, default=3.0)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--convergence-tol", type=float, default=0.1)
    p.add_argument("--convergence-sweeps", type=int, default=2)
    p.add_argument("--accept-tol", type=float, default=1e-8)
    p.add_argument("--progress-every", type=int, default=25)
    return p.parse_args()


def main():
    args = parse_args()
    if args.chunk_size == 0:
        args.chunk_size = None
    if args.n_cw > len(CW_LIBRARY):
        raise ValueError(f"n_cw={args.n_cw} exceeds CW_LIBRARY size={len(CW_LIBRARY)}")

    OUTDIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(args.seed)

    ent_psrs, disco_psrs, total_loaded = select_pulsars(args.n_psr)
    t0 = time.time()
    if args.data_mode == "pure":
        cw_params_list = apply_cw_overrides(CW_LIBRARY[: args.n_cw], args)
        residual_map, _ = inject_noisefree_cw(disco_psrs, cw_params_list)
        residual_map = add_phenomenological_noise(residual_map, disco_psrs, args)
        logl_fn, param_keys, base_values = build_pure_cw_likelihood(
            disco_psrs,
            residual_map,
            cw_params_list,
        )
        enterprise_params = None
        cw_block_names = None
    else:
        # Proper stochastic injection path from data_likelihood_sandbox.ipynb:
        # enterprise PTA simulation draws white noise, intrinsic red noise, and
        # HD-correlated GWB from the injected hyperparameters. Discovery then
        # evaluates a matching red/GWB GP likelihood with those hyperparameters
        # fixed during distance optimization.
        np.random.seed(args.noise_seed)
        pta, cw_block_names, _ = build_enterprise_pta(
            ent_psrs, args.n_cw, components=args.components
        )
        enterprise_params = generate_injection_params(
            pta,
            ent_psrs,
            args.n_cw,
            cw_block_names,
            log10_h=-12.0 if args.log10_h is None else args.log10_h,
            scenario=args.stochastic_scenario,
            rng=rng,
            gwb_log10_A=args.gwb_log10_a,
            gwb_gamma=args.gwb_gamma,
        )
        if args.log10_mc is not None:
            for name in cw_block_names:
                enterprise_params[f"{name}_log10_mc"] = args.log10_mc
        cw_params_list = cw_list_from_enterprise_params(enterprise_params, cw_block_names)
        sim_resids = simulate(pta, enterprise_params, sparse_cholesky=True)
        residual_map = {
            getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)
        }
        logl_fn, param_keys, base_values = build_disco_likelihood(
            disco_psrs,
            residual_map,
            num_cw=args.n_cw,
            enterprise_params=enterprise_params,
            cw_block_names=cw_block_names,
            components=args.components,
            include_gwb=not args.no_gwb,
        )

    truth_dist = np.array([psr.pdist[0] for psr in disco_psrs], dtype=float)
    mode_spacings = min_mode_spacings(disco_psrs, cw_params_list)
    pair_schedule = make_pair_schedule(len(disco_psrs), args.pair_mode)
    if args.max_pairs_per_sweep:
        pair_schedule_used = pair_schedule[: args.max_pairs_per_sweep]
    else:
        pair_schedule_used = pair_schedule

    print(f"loaded pulsars total: {total_loaded}")
    print(f"selected pulsars: {len(disco_psrs)} tightest distance priors")
    print(f"data mode: {args.data_mode}")
    print(f"CWs: {len(cw_params_list)}")
    print(f"pair mode: {args.pair_mode}; pairs/sweep: {len(pair_schedule_used)}")
    print(f"mode spacing min/median/max [kpc]: "
          f"{np.nanmin(mode_spacings):.3g} / {np.nanmedian(mode_spacings):.3g} / {np.nanmax(mode_spacings):.3g}")

    bad_residuals = [
        name for name, vals in residual_map.items() if not np.all(np.isfinite(vals))
    ]
    if bad_residuals:
        raise RuntimeError(
            "CW injection produced non-finite residuals for pulsars: "
            + ", ".join(bad_residuals[:10])
            + (" ..." if len(bad_residuals) > 10 else "")
        )
    truth_lnL = float(logl_fn(base_values))
    if not np.isfinite(truth_lnL):
        raise RuntimeError("truth lnL is non-finite; check CW parameters and selected pulsars")
    print(f"build+JIT truth lnL={truth_lnL:.3f}; seconds={time.time() - t0:.1f}")

    init_kinds = ["mean"] + ["random"] * max(0, args.n_chains - 1)
    results = []
    for chain_id, init_kind in enumerate(init_kinds):
        result = run_chain(
            chain_id,
            init_kind,
            logl_fn,
            param_keys,
            base_values,
            disco_psrs,
            cw_params_list,
            truth_dist,
            mode_spacings,
            pair_schedule,
            rng,
            args,
        )
        results.append(result)

    payload = {
        "config": vars(args),
        "truth_lnL": truth_lnL,
        "pulsars": [
            {
                "name": psr.name,
                "dist_mean_kpc": float(psr.pdist[0]),
                "dist_sigma_kpc": float(psr.pdist[1]),
                "min_mode_spacing_kpc": float(mode_spacings[ii]),
            }
            for ii, psr in enumerate(disco_psrs)
        ],
        "cw_params_list": cw_params_list,
        "enterprise_params": enterprise_params,
        "cw_block_names": cw_block_names,
        "results": results,
    }

    stamp = time.strftime("%Y%m%d_%H%M%S")
    json_path = OUTDIR / f"run_ncw{args.n_cw}_npsr{args.n_psr}_{args.pair_mode}_{stamp}.json"
    json_path.write_text(json.dumps(payload, indent=2))

    npz_path = json_path.with_suffix(".npz")
    np.savez(
        npz_path,
        truth_distances=truth_dist,
        mode_spacings=mode_spacings,
        recovered_distances=np.array([r["recovered_distances"] for r in results]),
        final_lnL=np.array([r["final_lnL"] for r in results]),
        truth_lnL=truth_lnL,
        pulsar_names=np.array([psr.name for psr in disco_psrs]),
    )
    print(f"saved {json_path}")
    print(f"saved {npz_path}")


if __name__ == "__main__":
    main()
