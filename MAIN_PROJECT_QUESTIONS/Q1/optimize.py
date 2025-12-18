"""Lightweight MAP helpers for discovery likelihoods."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Sequence

import jax
import jax.numpy as jnp


Array = jnp.ndarray
Params = Dict[str, float]


@dataclass
class MAPResult:
    params: Params
    trace: List[float]
    grad_norm: float
    converged: bool


def dict_to_array(param_keys: Sequence[str], params: Params) -> Array:
    return jnp.array([params[k] for k in param_keys], dtype=jnp.float64)


def array_to_dict(param_keys: Sequence[str], values: Array) -> Params:
    return {k: float(v) for k, v in zip(param_keys, values)}


def make_logposterior(logl_fn: Callable[[Params], float]) -> Callable[[Array, Sequence[str]], float]:
    """Wrap discovery logL(params) to work on flat arrays for JAX optimisation."""

    def _logpost(x: Array, keys: Sequence[str]) -> float:
        params = {k: v for k, v in zip(keys, x)}
        return logl_fn(params)

    return _logpost


def run_map_gradient_ascent(
    logpost: Callable[[Array, Sequence[str]], float],
    param_keys: Sequence[str],
    init_params: Params,
    learning_rate: float = 1e-2,
    max_steps: int = 400,
    beta: float = 0.9,
) -> MAPResult:
    """Simple momentum-based ascent for quick MAP prototyping."""
    x0 = dict_to_array(param_keys, init_params)

    def safe_logpost(z):
        val = logpost(z, param_keys)
        return val

    def step(state, _):
        x, v = state
        val, grad = jax.value_and_grad(safe_logpost)(x)
        v_new = beta * v + (1.0 - beta) * grad
        x_new = x + learning_rate * v_new
        return (x_new, v_new), (val, grad)

    (xf, vf), (vals, grads) = jax.lax.scan(step, (x0, jnp.zeros_like(x0)), None, length=max_steps)

    # pick last finite iterate; if none, fall back to start
    finite_mask = jnp.isfinite(vals)
    last_finite_idx = jnp.where(finite_mask, size=1, fill_value=0)[0][-1]
    xf = jax.lax.cond(
        finite_mask.any(),
        lambda _: xf,
        lambda _: x0,
        operand=None,
    )
    final_grad = grads[last_finite_idx] if finite_mask.any() else grads[0]
    final_val = vals[last_finite_idx] if finite_mask.any() else vals[0]

    grad_norm = float(jnp.linalg.norm(final_grad))
    out = MAPResult(
        params=array_to_dict(param_keys, xf),
        trace=[float(v) for v in vals],
        grad_norm=grad_norm,
        converged=grad_norm < 1e-3,
    )
    return out


def run_svi_map(
    logpost: Callable[[Array, Sequence[str]], float],
    param_keys: Sequence[str],
    init_params: Params,
    batch_size: int = 1000,
    max_batches: int = 50,
    patience: int = 3,
    peak_lr: float = 1e-2,
    rng_seed: int | None = None,
    gradient_clipping_val: float | None = None,
) -> MAPResult:
    """MAP via NumPyro AutoDelta SVI with early stopping."""
    try:
        import numpyro
        import numpyro.distributions as dist
        from numpyro.infer.autoguide import AutoDelta
        try:
            from pulsar_map_noise_estimates.map_noise_estimate import setup_svi, run_svi_early_stopping
        except ModuleNotFoundError:
            import os, sys
            default_path = "/home/mattm/soft/pulsar-map-noise-estimates/src"
            if os.path.isdir(default_path) and default_path not in sys.path:
                sys.path.insert(0, default_path)
            from pulsar_map_noise_estimates.map_noise_estimate import setup_svi, run_svi_early_stopping
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("NumPyro + pulsar_map_noise_estimates are required for SVI MAP") from exc

    def model():
        params = {}
        for k in param_keys:
            init = init_params.get(k, 0.0)
            if k.endswith("_cw_p_dist"):
                loc = jnp.log(jnp.clip(init if init > 0 else 1.0, 1e-3))
                z = numpyro.sample(k, dist.Normal(loc, 1.0))
                params[k] = jnp.exp(z)
            else:
                params[k] = numpyro.sample(k, dist.Normal(loc=init, scale=5.0))
        x = jnp.array([params[k] for k in param_keys])
        numpyro.factor("loglike", logpost(x, param_keys))

    guide = AutoDelta(model)
    svi = setup_svi(
        model,
        guide,
        max_epochs=max_batches * batch_size,
        peak_lr=peak_lr,
        num_warmup_steps=batch_size,
        gradient_clipping_val=gradient_clipping_val,
    )

    rng_key = jax.random.PRNGKey(0 if rng_seed is None else rng_seed)
    params = run_svi_early_stopping(
        rng_key,
        svi,
        batch_size=batch_size,
        patience=patience,
        max_num_batches=max_batches,
    )

    # strip any AutoDelta suffixes
    cleaned_raw = {key.removesuffix("_auto_loc"): val for key, val in params.items()}
    clean = {}
    for k in param_keys:
        v = cleaned_raw.get(k)
        if v is None:
            continue
        if k.endswith("_cw_p_dist"):
            clean[k] = float(jnp.exp(v))
        else:
            clean[k] = float(v)
    return MAPResult(params=clean, trace=[], grad_norm=float("nan"), converged=True)


def run_lbfgs(
    logpost: Callable[[Array, Sequence[str]], float],
    param_keys: Sequence[str],
    init_params: Params,
    max_iters: int = 200,
) -> MAPResult:
    """LBFGS using JAX's built-in optimizer to refine a MAP point."""
    try:
        from jaxopt import LBFGS
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("jaxopt is required for LBFGS refinement") from exc

    def fun(x):
        return -logpost(x, param_keys)  # jaxopt minimizes

    solver = LBFGS(fun=fun, maxiter=max_iters)
    x0 = dict_to_array(param_keys, init_params)
    result = solver.run(x0)
    xf = result.params
    grad = jax.grad(fun)(xf)
    grad_norm = float(jnp.linalg.norm(grad))
    return MAPResult(
        params=array_to_dict(param_keys, xf),
        trace=[float(-fun(xf))],
        grad_norm=grad_norm,
        converged=result.state.error is not None and result.state.error < 1e-4,
    )
