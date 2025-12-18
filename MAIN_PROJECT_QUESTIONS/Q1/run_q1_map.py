"""Minimal end-to-end MAP run for one Q1 scenario."""

from __future__ import annotations

import argparse
from pathlib import Path

from sim import (
    DistanceTier,
    Q1Config,
    apply_distance_tier,
    build_discovery_loglike,
    build_enterprise_pta,
    draw_cw_population,
    draw_noise_params,
    inject_cw_params,
    load_discovery_pulsars,
    load_enterprise_pulsars,
    make_rng,
    set_toaerrs,
    simulate_enterprise_residuals,
)
from optimize import make_logposterior, run_map_gradient_ascent, run_svi_map


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a single Q1 MAP optimisation.")
    parser.add_argument("--scenario", choices=["A", "B", "C"], default="A", help="Noise scenario: A WN, B WN+RN, C WN+RN+GWB.")
    parser.add_argument("--tier-scale", type=float, default=1.0, help="Distance sigma multiplier (1=real, <1 tightens, >1 loosens).")
    parser.add_argument("--tier-name", default="real", help="Label for the tier.")
    parser.add_argument("--n-pulsars", type=int, default=30)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--steps", type=int, default=200, help="Gradient ascent iterations.")
    parser.add_argument("--lr", type=float, default=5e-3, help="Learning rate for the MAP ascent.")
    parser.add_argument("--opt-mode", choices=["map", "svi"], default="svi", help="Optimiser: simple MAP ascent or NumPyro SVI AutoDelta.")
    parser.add_argument("--grad-clip", type=float, default=None, help="Gradient clipping value for SVI (None to disable).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = Q1Config(
        scenario=args.scenario,
        distance_tier=DistanceTier(args.tier_name, args.tier_scale),
        seed=args.seed,
        n_pulsars=args.n_pulsars,
    )

    rng = make_rng(cfg.seed)
    cw_draws = draw_cw_population(cfg, rng)

    e_psrs = load_enterprise_pulsars(cfg)
    set_toaerrs(e_psrs, sigma=1e-6)
    true_distances = apply_distance_tier(e_psrs, cfg.distance_tier)
    pta, ring_vecs = build_enterprise_pta(e_psrs, cfg, cw_draws)
    pta.set_default_params({})

    base_params = draw_noise_params(pta, cfg, rng, psrs=e_psrs)
    inj_params = inject_cw_params(base_params, cw_draws, cfg)

    sim_resids = simulate_enterprise_residuals(pta, inj_params, sparse_cholesky=True)
    name_to_resid = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}

    d_psrs = load_discovery_pulsars(cfg)
    residual_map = {psr.name: name_to_resid[psr.name] for psr in d_psrs}

    logl, param_keys = build_discovery_loglike(d_psrs, cw_draws, cfg, ring_vecs=ring_vecs, residual_map=residual_map)
    logpost = make_logposterior(logl)

    init_params = {k: inj_params.get(k, 0.0) for k in param_keys}
    if args.opt_mode == "svi":
        result = run_svi_map(
            logpost,
            param_keys,
            init_params,
            batch_size=500,
            max_batches=20,
            patience=3,
            peak_lr=1e-2,
            rng_seed=cfg.seed or 0,
            gradient_clipping_val=args.grad_clip,
        )
    else:
        result = run_map_gradient_ascent(
            logpost,
            param_keys,
            init_params,
            learning_rate=args.lr,
            max_steps=args.steps,
        )

    print(f"Finished MAP run. converged={result.converged} grad_norm={result.grad_norm:.3e}")
    if result.trace:
        print(f"logpost trace (first/last): {result.trace[0]:.3f} -> {result.trace[-1]:.3f}")
    print("Selected parameters (subset):")
    for key in param_keys[: min(8, len(param_keys))]:
        val = result.params[key]
        inj = inj_params.get(key, None)
        delta = None if inj is None else val - inj
        if delta is None:
            print(f"  {key}: {val:.4f}")
        else:
            print(f"  {key}: {val:.4f} (inj {inj:.4f}, delta {delta:+.2e})")


if __name__ == "__main__":
    main()
