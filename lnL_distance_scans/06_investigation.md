# 06 - Realistic Noise Failure Investigation

Lean diagnostics for realistic-noise distance recovery. No beam optimizer, no scaling runs.


```python
from pathlib import Path
import sys
import time
import copy
import importlib
import itertools
import concurrent.futures as cf

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd
import scipy.optimize
from scipy.signal import find_peaks
from IPython.display import display

HERE = Path.cwd()
sys.path.insert(0, str((HERE / "../CW_lnL_check").resolve()))

import cw_helpers
importlib.reload(cw_helpers)

from cw_helpers import (
    build_disco_likelihood,
    build_fast_scan_likelihood,
    build_enterprise_pta,
    build_pure_cw_likelihood,
    compute_mode_spacing,
    generate_injection_params,
    inject_noisefree_cw,
    load_pulsars,
    make_distance_optimizer,
    scan_pulsar_distance,
    simulate,
)

# Control panel
DATA_MODE = "stochastic"
N_PSR = 40
N_CW = 4
LOG10_H = -12.0
LOG10_MC = None
RNG_SEED = 12345
NOISE_SEED = 24680

STOCHASTIC_SCENARIO = "well_separated"
COMPONENTS = 30
INCLUDE_GWB = True
GWB_LOG10_A = -14.5
GWB_GAMMA = 13 / 3
INCLUDE_RN = True
RN_COMPONENTS = 30

TRUTH_DISTANCE_MODE = "gaussian_prior_draw"
TRUTH_SIGMA_CLIP = 3.0
TRUTH_SIGMA_MULTIPLIER = 1.0
MODE_PRIOR_SIGMA_WIDTH = 3.0

# Timing knobs
FAST_GRADIENT_MAXITER = 100
INITIAL_GRADIENT_MAXITER = 100
FAST_MODE_STARTS = 20
FAST_PRIOR_STARTS = 5
OPT_FACTR = 1e7
DIAGNOSTIC_SCAN_POINTS = 2000
MAX_1D_SCAN_POINTS = 3000
N_CONTEXT_RANDOM = 2
MODE_SAMPLES_PER_MODE = 6
TOP_MODE_CANDIDATES = 200
PARALLEL_CANDIDATE_SCANS = True
MAX_CANDIDATE_WORKERS = 4
MAKE_PLOTS = False
FAILED_PULSAR_INDICES = None  # set to a list from notebook 05 to diagnose that exact failed set

# Experiment gates
RUN_SANITY_GRADIENT = True
RUN_EXP1 = True
RUN_EXP2 = True
RUN_EXP3 = True
RUN_EXP4 = True
RUN_EXP5 = True

REPRESENTATIVE_INDICES = [0, 4, 10, 20, 34]
PROGRESSIVE_N_HARD = 10
BLOCK_SIZE = 8
BLOCK_CYCLES = 1
EXP4_STRATEGIES = ["sigma_sorted_8"]  # expand to random_8/sky_8 after timing is acceptable

rng = np.random.default_rng(RNG_SEED)
np.random.seed(NOISE_SEED)

print(HERE)
```

    /home/mattm/miniforge3/envs/discotech/lib/python3.12/site-packages/enterprise/signals/utils.py:13: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
      from pkg_resources import Requirement, resource_filename


    /home/mattm/projects/HSYMT/lnL_distance_scans



```python
CW_LIBRARY = [
    dict(cos_gwtheta=0.30, gwphi=2.50, cos_inc=-0.20, phase0=1.00, psi=0.70, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.00),
    dict(cos_gwtheta=-0.50, gwphi=0.80, cos_inc=0.40, phase0=2.10, psi=1.30, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.80),
    dict(cos_gwtheta=0.05, gwphi=4.20, cos_inc=-0.65, phase0=0.35, psi=2.30, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.20),
    dict(cos_gwtheta=0.70, gwphi=5.40, cos_inc=0.10, phase0=1.70, psi=0.20, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.90),
    dict(cos_gwtheta=-0.10, gwphi=3.30, cos_inc=0.80, phase0=2.90, psi=1.80, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.10),
    dict(cos_gwtheta=0.45, gwphi=1.60, cos_inc=-0.35, phase0=0.80, psi=2.70, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.70),
    dict(cos_gwtheta=-0.75, gwphi=5.90, cos_inc=0.55, phase0=2.40, psi=0.45, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.35),
    dict(cos_gwtheta=0.18, gwphi=0.25, cos_inc=-0.85, phase0=1.35, psi=1.05, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.60),
    dict(cos_gwtheta=-0.32, gwphi=2.05, cos_inc=0.25, phase0=3.05, psi=2.05, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.45),
    dict(cos_gwtheta=0.88, gwphi=3.75, cos_inc=-0.05, phase0=0.15, psi=0.95, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.95),
    dict(cos_gwtheta=-0.62, gwphi=4.75, cos_inc=0.72, phase0=1.95, psi=2.85, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.25),
    dict(cos_gwtheta=0.58, gwphi=0.95, cos_inc=-0.48, phase0=2.75, psi=1.55, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.75),
    dict(cos_gwtheta=-0.18, gwphi=5.15, cos_inc=0.08, phase0=0.55, psi=0.10, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.05),
    dict(cos_gwtheta=0.02, gwphi=1.20, cos_inc=-0.70, phase0=2.25, psi=2.45, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.85),
    dict(cos_gwtheta=-0.88, gwphi=3.95, cos_inc=0.38, phase0=1.15, psi=1.75, log10_h=-12.0, log10_mc=9.00, log10_fgw=-8.30),
    dict(cos_gwtheta=0.36, gwphi=4.55, cos_inc=-0.18, phase0=2.55, psi=0.60, log10_h=-12.0, log10_mc=9.00, log10_fgw=-7.65),
]


def clipped_normal(mean, sigma, rng, clip):
    draw = rng.normal(mean, sigma)
    lo = np.maximum(0.01, mean - clip * sigma)
    hi = mean + clip * sigma
    return np.clip(draw, lo, hi)


def select_pulsars(n_psr):
    ent_all, disco_all = load_pulsars(None)
    pairs = sorted(zip(ent_all, disco_all), key=lambda pair: pair[1].pdist[1])
    selected = pairs[:n_psr]
    return [p[0] for p in selected], [p[1] for p in selected], len(disco_all)


def clone_psrs_with_distances(psrs, distances):
    clones = copy.deepcopy(psrs)
    for psr, dist in zip(clones, distances):
        sigma = float(psr.pdist[1])
        try:
            psr.pdist = (float(dist), sigma)
        except Exception:
            psr._pdist = (float(dist), sigma)
    return clones


def make_truth_distances(prior_mean, prior_sigma, rng):
    if TRUTH_DISTANCE_MODE == "prior_mean":
        return prior_mean.copy()
    if TRUTH_DISTANCE_MODE == "gaussian_prior_draw":
        return clipped_normal(prior_mean, prior_sigma * TRUTH_SIGMA_MULTIPLIER, rng, TRUTH_SIGMA_CLIP)
    raise ValueError(TRUTH_DISTANCE_MODE)


def cw_params_from_library(n_cw, log10_h=None, log10_mc=None):
    if n_cw > len(CW_LIBRARY):
        raise ValueError(f"N_CW={n_cw} exceeds CW_LIBRARY size={len(CW_LIBRARY)}")
    out = [dict(cw) for cw in CW_LIBRARY[:n_cw]]
    for cw in out:
        if log10_h is not None:
            cw["log10_h"] = log10_h
        if log10_mc is not None:
            cw["log10_mc"] = log10_mc
    return out


def cw_list_from_enterprise_params(enterprise_params, cw_block_names):
    keys = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi")
    return [{key: float(enterprise_params[f"{name}_{key}"]) for key in keys} for name in cw_block_names]


def distance_key(psr):
    return f"{psr.name}_cw_p_dist"


def set_distances(values, param_keys, disco_psrs, distances):
    out = np.array(values, dtype=float).copy()
    for psr, dist in zip(disco_psrs, distances):
        out[param_keys.index(distance_key(psr))] = dist
    return out


def min_mode_spacings(disco_psrs, cw_params_list):
    out = []
    for psr in disco_psrs:
        vals = [compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"], cw["log10_fgw"], psr.pos) for cw in cw_params_list]
        out.append(float(np.nanmin(vals)))
    return np.array(out)


def score_distances(distances):
    err_modes = (np.asarray(distances) - truth_dist) / mode_spacings
    return dict(
        n_recovered=int(np.sum(np.abs(err_modes) < 0.5)),
        frac_recovered=float(np.mean(np.abs(err_modes) < 0.5)),
        median_abs_modes=float(np.nanmedian(np.abs(err_modes))),
        max_abs_modes=float(np.nanmax(np.abs(err_modes))),
    )


def top_local_maxima(x, y, k):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(y)
    x, y = x[finite], y[finite]
    if len(x) == 0:
        return np.array([]), np.array([])
    if len(x) >= 3:
        idx = np.where(np.r_[False, (y[1:-1] >= y[:-2]) & (y[1:-1] >= y[2:]), False])[0]
    else:
        idx = np.array([], dtype=int)
    if len(idx) == 0:
        idx = np.array([int(np.nanargmax(y))])
    idx = idx[np.argsort(y[idx])[::-1]][:k]
    return x[idx], y[idx]
```

## Build Injection And Likelihood


```python
print("BUILD INJECTION AND LIKELIHOOD")
t0 = time.time()

ent_psrs, disco_psrs, total_loaded = select_pulsars(N_PSR)
prior_mean = np.array([p.pdist[0] for p in disco_psrs], dtype=float)
sigmas = np.array([p.pdist[1] for p in disco_psrs], dtype=float)
truth_dist = make_truth_distances(prior_mean, sigmas, rng)
ent_psrs_inj = clone_psrs_with_distances(ent_psrs, truth_dist)
disco_psrs_inj = clone_psrs_with_distances(disco_psrs, truth_dist)

if DATA_MODE == "pure":
    cw_params_list = cw_params_from_library(N_CW, LOG10_H, LOG10_MC)
    residual_map, _ = inject_noisefree_cw(disco_psrs_inj, cw_params_list)
    logl_fn, param_keys, base_values = build_pure_cw_likelihood(disco_psrs, residual_map, cw_params_list)
    enterprise_params = None
    cw_block_names = None
else:
    pta, cw_block_names, _ = build_enterprise_pta(
        ent_psrs_inj, N_CW, components=COMPONENTS,
        include_rn=INCLUDE_RN, rn_components=RN_COMPONENTS,
    )
    enterprise_params = generate_injection_params(
        pta, ent_psrs_inj, N_CW, cw_block_names,
        log10_h=LOG10_H,
        scenario=STOCHASTIC_SCENARIO,
        rng=rng,
        gwb_log10_A=GWB_LOG10_A,
        gwb_gamma=GWB_GAMMA,
        include_rn=INCLUDE_RN,
    )
    if LOG10_MC is not None:
        for name in cw_block_names:
            enterprise_params[f"{name}_log10_mc"] = LOG10_MC
    cw_params_list = cw_list_from_enterprise_params(enterprise_params, cw_block_names)
    sim_resids = simulate(pta, enterprise_params, sparse_cholesky=True)
    residual_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}
    logl_fn_disco, param_keys_disco, base_values_disco = build_disco_likelihood(
        disco_psrs, residual_map,
        num_cw=N_CW,
        enterprise_params=enterprise_params,
        cw_block_names=cw_block_names,
        components=COMPONENTS,
        include_gwb=INCLUDE_GWB,
        include_rn=INCLUDE_RN,
        rn_components=RN_COMPONENTS,
    )
    if INCLUDE_GWB:
        try:
            logl_fn_fast, param_keys_fast, base_values_fast = build_fast_scan_likelihood(
                disco_psrs, residual_map,
                num_cw=N_CW,
                enterprise_params=enterprise_params,
                cw_block_names=cw_block_names,
                components=COMPONENTS,
                include_rn=INCLUDE_RN,
                rn_components=RN_COMPONENTS,
            )
            truth_probe = set_distances(base_values_fast, param_keys_fast, disco_psrs, truth_dist)
            ll_fast = float(logl_fn_fast(jnp.asarray(truth_probe)))
            ll_disco = float(logl_fn_disco(jnp.asarray(truth_probe)))
            print(f"fast/disco truth lnL delta={ll_fast - ll_disco:.3e}")
            if param_keys_fast != param_keys_disco or abs(ll_fast - ll_disco) > 1e-6:
                print("WARNING: fast likelihood validation failed; using discovery likelihood.")
                logl_fn, param_keys, base_values = logl_fn_disco, param_keys_disco, base_values_disco
            else:
                logl_fn, param_keys, base_values = logl_fn_fast, param_keys_fast, base_values_fast
        except Exception as e:
            print(f"WARNING: fast likelihood failed ({e}); using discovery likelihood.")
            logl_fn, param_keys, base_values = logl_fn_disco, param_keys_disco, base_values_disco
    else:
        logl_fn, param_keys, base_values = logl_fn_disco, param_keys_disco, base_values_disco

mode_spacings = min_mode_spacings(disco_psrs, cw_params_list)
hardness = sigmas / mode_spacings
truth_values = set_distances(base_values, param_keys, disco_psrs, truth_dist)
prior_values = set_distances(base_values, param_keys, disco_psrs, prior_mean)
truth_lnL = float(logl_fn(jnp.asarray(truth_values)))
prior_mean_lnL = float(logl_fn(jnp.asarray(prior_values)))

print(f"loaded total pulsars: {total_loaded}")
print(f"selected N_PSR={N_PSR}, N_CW={N_CW}, DATA_MODE={DATA_MODE}")
print(f"truth lnL={truth_lnL:.3f}")
print(f"prior mean lnL={prior_mean_lnL:.3f}, delta={prior_mean_lnL - truth_lnL:.3f}")
print(f"mode spacing min/median/max={np.nanmin(mode_spacings):.3g}/{np.nanmedian(mode_spacings):.3g}/{np.nanmax(mode_spacings):.3g} kpc")
print(f"build seconds={time.time() - t0:.1f}")
```

    BUILD INJECTION AND LIKELIHOOD
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/B1855+09.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/B1855+09.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/B1937+21.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/B1937+21.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/B1953+29.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/B1953+29.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0023+0923.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0023+0923.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0030+0451.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0030+0451.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0125-2327.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0125-2327.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0340+4130.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0340+4130.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0406+3039.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0406+3039.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0437-4715.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0437-4715.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0509+0856.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0509+0856.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0557+1551.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0557+1551.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0605+3757.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0605+3757.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0610-2100.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0610-2100.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0613-0200.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0613-0200.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0614-3329.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0614-3329.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0636+5128.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0636+5128.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0636-3044.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0636-3044.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0645+5158.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0645+5158.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0709+0458.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0709+0458.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0711-6830.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0711-6830.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0740+6620.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0740+6620.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0751+1807.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0751+1807.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0900-3144.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0900-3144.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0931-1902.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0931-1902.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J0955-6150.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J0955-6150.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1012+5307.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1012+5307.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1012-4235.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1012-4235.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1017-7156.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1017-7156.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1022+1001.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1022+1001.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1024-0719.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1024-0719.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1036-8317.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1036-8317.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1045-4509.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1045-4509.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1101-6424.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1101-6424.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1103-5403.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1103-5403.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1125+7819.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1125+7819.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1125-5825.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1125-5825.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1125-6014.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1125-6014.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1216-6410.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1216-6410.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1312+0051.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1312+0051.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1327-0755.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1327-0755.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1421-4409.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1421-4409.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1431-5740.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1431-5740.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1435-6100.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1435-6100.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1446-4701.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1446-4701.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1453+1902.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1453+1902.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1455-3330.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1455-3330.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1525-5545.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1525-5545.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1543-5149.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1543-5149.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1545-4550.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1545-4550.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1547-5709.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1547-5709.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1600-3053.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1600-3053.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1603-7202.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1603-7202.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1614-2230.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1614-2230.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1629-6902.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1629-6902.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1630+3734.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1630+3734.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1640+2224.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1640+2224.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1643-1224.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1643-1224.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1652-4838.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1652-4838.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1653-2054.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1653-2054.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1658-5324.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1658-5324.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1705-1903.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1705-1903.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1708-3506.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1708-3506.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1713+0747.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1713+0747.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1719-1438.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1719-1438.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1721-2457.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1721-2457.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1730-2304.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1730-2304.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1732-5049.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1732-5049.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1737-0811.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1737-0811.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1738+0333.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1738+0333.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1741+1351.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1741+1351.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1744-1134.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1744-1134.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1745+1017.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1745+1017.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1747-4036.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1747-4036.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1751-2857.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1751-2857.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1756-2251.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1756-2251.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1757-5322.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1757-5322.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1801-1417.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1801-1417.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1802-2124.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1802-2124.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1804-2717.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1804-2717.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1811-2405.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1811-2405.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1824-2452A.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1824-2452A.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1825-0319.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1825-0319.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1832-0836.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1832-0836.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1843-1113.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1843-1113.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1853+1303.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1853+1303.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1902-5105.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1902-5105.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1903+0327.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1903+0327.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1903-7051.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1903-7051.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1909-3744.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1909-3744.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1910+1256.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1910+1256.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1911+1347.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1911+1347.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1918-0642.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1918-0642.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1923+2515.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1923+2515.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1933-6211.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1933-6211.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1944+0907.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1944+0907.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1946+3417.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1946+3417.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J1946-5403.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J1946-5403.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2010-1323.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2010-1323.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2017+0603.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2017+0603.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2033+1734.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2033+1734.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2039-3616.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2039-3616.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2043+1711.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2043+1711.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2124-3358.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2124-3358.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2129-5721.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2129-5721.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2145-0750.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2145-0750.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2150-0326.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2150-0326.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2214+3000.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2214+3000.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2222-0137.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2222-0137.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2229+2643.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2229+2643.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2234+0611.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2234+0611.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2234+0944.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2234+0944.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2241-5236.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2241-5236.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2302+4442.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2302+4442.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2317+1439.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2317+1439.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2322+2057.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2322+2057.feather.
    FeatherPulsar.read_feather: cannot find dmx in feather file /home/mattm/projects/HSYMT/data_products/J2322-2650.feather.
    FeatherPulsar.read_feather: cannot find _pdist in feather file /home/mattm/projects/HSYMT/data_products/J2322-2650.feather.


    Duplicate signal J0437-4715_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0437-4715_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0437-4715_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0437-4715_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0437-4715_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0437-4715_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0437-4715_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0437-4715_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0437-4715_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0030+0451_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0030+0451_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0030+0451_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0030+0451_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0030+0451_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0030+0451_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0030+0451_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0030+0451_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0030+0451_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1744-1134_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1744-1134_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1744-1134_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1744-1134_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1744-1134_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1744-1134_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1744-1134_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1744-1134_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1744-1134_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1909-3744_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1909-3744_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1909-3744_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1909-3744_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1909-3744_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1909-3744_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1909-3744_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1909-3744_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1909-3744_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1022+1001_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1022+1001_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1022+1001_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1022+1001_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1022+1001_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1022+1001_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1022+1001_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1022+1001_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1022+1001_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1713+0747_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1713+0747_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1713+0747_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1713+0747_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1713+0747_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1713+0747_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1713+0747_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1713+0747_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1713+0747_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0711-6830_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0711-6830_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0711-6830_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0711-6830_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0711-6830_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0711-6830_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0711-6830_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0711-6830_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0711-6830_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1012+5307_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1012+5307_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1012+5307_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1012+5307_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1012+5307_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1012+5307_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1012+5307_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1012+5307_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1012+5307_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1730-2304_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1730-2304_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1730-2304_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1730-2304_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1730-2304_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1730-2304_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1730-2304_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1730-2304_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1730-2304_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2145-0750_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J2145-0750_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2145-0750_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2145-0750_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J2145-0750_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2145-0750_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2145-0750_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J2145-0750_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2145-0750_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1630+3734_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1630+3734_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1630+3734_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1630+3734_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1630+3734_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1630+3734_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1630+3734_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1630+3734_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1630+3734_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1614-2230_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1614-2230_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1614-2230_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1614-2230_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1614-2230_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1614-2230_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1614-2230_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1614-2230_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1614-2230_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1737-0811_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1737-0811_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1737-0811_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1737-0811_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1737-0811_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1737-0811_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1737-0811_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1737-0811_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1737-0811_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1024-0719_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1024-0719_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1024-0719_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1024-0719_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1024-0719_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1024-0719_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1024-0719_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1024-0719_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1024-0719_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2322-2650_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J2322-2650_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2322-2650_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2322-2650_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J2322-2650_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2322-2650_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2322-2650_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J2322-2650_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2322-2650_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2222-0137_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J2222-0137_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2222-0137_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2222-0137_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J2222-0137_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2222-0137_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2222-0137_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J2222-0137_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2222-0137_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2124-3358_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J2124-3358_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2124-3358_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2124-3358_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J2124-3358_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2124-3358_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2124-3358_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J2124-3358_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2124-3358_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1643-1224_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1643-1224_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1643-1224_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1643-1224_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1643-1224_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1643-1224_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1643-1224_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1643-1224_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1643-1224_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0751+1807_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0751+1807_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0751+1807_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0751+1807_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0751+1807_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0751+1807_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0751+1807_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0751+1807_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0751+1807_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1045-4509_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1045-4509_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1045-4509_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1045-4509_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1045-4509_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1045-4509_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1045-4509_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1045-4509_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1045-4509_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1640+2224_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1640+2224_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1640+2224_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1640+2224_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1640+2224_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1640+2224_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1640+2224_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1640+2224_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1640+2224_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0613-0200_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0613-0200_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0613-0200_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0613-0200_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0613-0200_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0613-0200_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0613-0200_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0613-0200_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0613-0200_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1738+0333_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1738+0333_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1738+0333_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1738+0333_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1738+0333_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1738+0333_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1738+0333_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1738+0333_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1738+0333_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0023+0923_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0023+0923_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0023+0923_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0023+0923_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0023+0923_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0023+0923_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0023+0923_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0023+0923_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0023+0923_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1918-0642_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1918-0642_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1918-0642_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1918-0642_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1918-0642_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1918-0642_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1918-0642_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1918-0642_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1918-0642_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2043+1711_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J2043+1711_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2043+1711_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2043+1711_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J2043+1711_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2043+1711_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2043+1711_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J2043+1711_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2043+1711_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal B1855+09_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, B1855+09_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, B1855+09_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal B1855+09_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, B1855+09_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, B1855+09_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal B1855+09_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, B1855+09_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, B1855+09_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0614-3329_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0614-3329_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0614-3329_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0614-3329_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0614-3329_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0614-3329_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0614-3329_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0614-3329_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0614-3329_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1125+7819_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1125+7819_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1125+7819_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1125+7819_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1125+7819_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1125+7819_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1125+7819_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1125+7819_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1125+7819_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1933-6211_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1933-6211_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1933-6211_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1933-6211_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1933-6211_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1933-6211_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1933-6211_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1933-6211_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1933-6211_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0636-3044_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0636-3044_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0636-3044_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0636-3044_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0636-3044_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0636-3044_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0636-3044_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0636-3044_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0636-3044_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0605+3757_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0605+3757_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0605+3757_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0605+3757_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0605+3757_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0605+3757_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0605+3757_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0605+3757_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0605+3757_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0610-2100_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0610-2100_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0610-2100_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0610-2100_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0610-2100_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0610-2100_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0610-2100_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0610-2100_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0610-2100_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1756-2251_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1756-2251_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1756-2251_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1756-2251_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1756-2251_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1756-2251_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1756-2251_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1756-2251_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1756-2251_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1312+0051_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1312+0051_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1312+0051_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1312+0051_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1312+0051_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1312+0051_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1312+0051_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1312+0051_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1312+0051_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1853+1303_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1853+1303_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1853+1303_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1853+1303_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1853+1303_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1853+1303_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1853+1303_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1853+1303_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1853+1303_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2234+0611_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J2234+0611_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2234+0611_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2234+0611_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J2234+0611_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2234+0611_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J2234+0611_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J2234+0611_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J2234+0611_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0125-2327_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0125-2327_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0125-2327_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0125-2327_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0125-2327_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0125-2327_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0125-2327_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0125-2327_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0125-2327_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1658-5324_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J1658-5324_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1658-5324_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1658-5324_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J1658-5324_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1658-5324_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J1658-5324_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J1658-5324_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J1658-5324_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0900-3144_cw from objects <Enterprise Signal object cw[cw2_cos_gwtheta, cw2_gwphi, cw2_cos_inc, cw2_log10_mc, cw2_log10_fgw, cw2_log10_h, cw2_phase0, cw2_psi, J0900-3144_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0900-3144_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0900-3144_cw from objects <Enterprise Signal object cw[cw3_cos_gwtheta, cw3_gwphi, cw3_cos_inc, cw3_log10_mc, cw3_log10_fgw, cw3_log10_h, cw3_phase0, cw3_psi, J0900-3144_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0900-3144_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    
    Duplicate signal J0900-3144_cw from objects <Enterprise Signal object cw[cw4_cos_gwtheta, cw4_gwphi, cw4_cos_inc, cw4_log10_mc, cw4_log10_fgw, cw4_log10_h, cw4_phase0, cw4_psi, J0900-3144_cw_p_dist]> and <Enterprise Signal object cw[cw_cos_gwtheta, cw_gwphi, cw_cos_inc, cw_log10_mc, cw_log10_fgw, cw_log10_h, cw_phase0, cw_psi, J0900-3144_cw_p_dist]>.
    This functionality was added in v1.1.0 and may cause post v1.1.0 functionality to break.
    This may not cause other errors but it is recommended that you use a custom name for one of the duplicate signals.
    


    fast/disco truth lnL delta=0.000e+00
    loaded total pulsars: 116
    selected N_PSR=40, N_CW=4, DATA_MODE=stochastic
    truth lnL=158941.139
    prior mean lnL=-1679351.629, delta=-1838292.769
    mode spacing min/median/max=0.000219/0.000388/0.00116 kpc
    build seconds=47.2


## Shared Diagnostic Helpers


```python
def scan_bounds(i, n_points=None):
    lo = max(0.01, prior_mean[i] - MODE_PRIOR_SIGMA_WIDTH * sigmas[i])
    hi = prior_mean[i] + MODE_PRIOR_SIGMA_WIDTH * sigmas[i]
    if n_points is not None:
        return lo, hi, int(n_points)
    dL = mode_spacings[i]
    if np.isfinite(dL) and dL > 0:
        n = int(np.ceil((hi - lo) / dL * MODE_SAMPLES_PER_MODE)) + 1
    else:
        n = MAX_1D_SCAN_POINTS
    return lo, hi, int(np.clip(n, 25, MAX_1D_SCAN_POINTS))


def conditional_scan(i, context_distances, n_points=DIAGNOSTIC_SCAN_POINTS, required_points=None):
    lo, hi, n = scan_bounds(i, n_points)
    req = [truth_dist[i], prior_mean[i]]
    if required_points is not None:
        req += list(np.ravel(required_points))
    vals = set_distances(base_values, param_keys, disco_psrs, context_distances)
    x, y = scan_pulsar_distance(
        logl_fn, vals, param_keys, distance_key(disco_psrs[i]),
        lo, hi, n_points=n, chunk_size=None, n_components=COMPONENTS,
        required_points=req,
    )
    return np.asarray(x, dtype=float), np.asarray(y, dtype=float)


def scan_metrics(i, x, y):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(y)
    if not np.any(finite):
        return dict(contrast=np.nan, n_degenerate=np.nan, truth_is_global_max=False, truth_peak=np.nan, global_peak=np.nan)

    peak_x, peak_y = top_local_maxima(x, y, max(300, TOP_MODE_CANDIDATES))
    if len(peak_x) == 0:
        return dict(contrast=np.nan, n_degenerate=np.nan, truth_is_global_max=False, truth_peak=np.nan, global_peak=float(x[np.nanargmax(y)]))

    truth_peak_idx = int(np.argmin(np.abs(peak_x - truth_dist[i])))
    truth_peak = float(peak_x[truth_peak_idx])
    truth_peak_y = float(peak_y[truth_peak_idx])
    global_idx = int(np.nanargmax(y))
    global_peak = float(x[global_idx])
    truth_is_global_max = abs(global_peak - truth_peak) < 0.5 * mode_spacings[i] or truth_peak_y >= float(np.nanmax(y)) - 1e-6

    near = np.where(np.abs(x - truth_peak) <= 0.5 * mode_spacings[i])[0]
    if len(near) >= 3:
        left = near[0]
        right = near[-1]
        local_min = float(np.nanmin(y[left:right + 1]))
    else:
        local_min = float(np.nanmin(y[finite]))
    contrast = truth_peak_y - local_min
    n_degenerate = int(np.sum(peak_y >= truth_peak_y - 1.0))

    return dict(
        contrast=float(contrast),
        n_degenerate=n_degenerate,
        truth_is_global_max=bool(truth_is_global_max),
        truth_peak=truth_peak,
        global_peak=global_peak,
    )


def choose_peak_distance(i, context_distances, n_points=DIAGNOSTIC_SCAN_POINTS):
    x, y = conditional_scan(i, context_distances, n_points=n_points)
    px, py = top_local_maxima(x, y, max(300, TOP_MODE_CANDIDATES if hardness[i] <= 500 else 300))
    if len(px) == 0:
        return float(x[np.nanargmax(y)]), float(np.nanmax(y))
    return float(px[0]), float(py[0])


def make_candidate_for_pulsar(i):
    context_rng = np.random.default_rng(RNG_SEED + 1000 + i)
    contexts = [prior_mean.copy()]
    for _ in range(N_CONTEXT_RANDOM):
        contexts.append(clipped_normal(prior_mean, sigmas, context_rng, MODE_PRIOR_SIGMA_WIDTH))
    contexts.append(truth_dist.copy())

    peaks = []
    scores = []
    for ctx in contexts:
        x, y = conditional_scan(i, ctx, n_points=None)
        k = 300 if hardness[i] > 500 else TOP_MODE_CANDIDATES
        px, py = top_local_maxima(x, y, k)
        peaks.extend(px.tolist())
        scores.extend(py.tolist())
    if len(peaks) == 0:
        return np.array([prior_mean[i]]), np.array([0.0])
    df = pd.DataFrame({"distance": peaks, "score": scores})
    df["mode_bin"] = np.round(df["distance"] / mode_spacings[i]).astype(int)
    df = df.sort_values("score", ascending=False).drop_duplicates("mode_bin")
    k = 300 if hardness[i] > 500 else TOP_MODE_CANDIDATES
    df = df.head(k)
    return df["distance"].to_numpy(float), df["score"].to_numpy(float)


def build_empirical_candidates():
    print("Building empirical 1D candidates")
    t0 = time.time()
    all_modes = [None] * N_PSR
    all_scores = [None] * N_PSR
    if PARALLEL_CANDIDATE_SCANS:
        with cf.ThreadPoolExecutor(max_workers=MAX_CANDIDATE_WORKERS) as ex:
            futs = {ex.submit(make_candidate_for_pulsar, i): i for i in range(N_PSR)}
            for n, fut in enumerate(cf.as_completed(futs), 1):
                i = futs[fut]
                all_modes[i], all_scores[i] = fut.result()
                if n % 5 == 0 or n == N_PSR:
                    print(f"  candidates {n}/{N_PSR}")
    else:
        for i in range(N_PSR):
            all_modes[i], all_scores[i] = make_candidate_for_pulsar(i)
            if (i + 1) % 5 == 0 or i + 1 == N_PSR:
                print(f"  candidates {i + 1}/{N_PSR}")
    print(f"candidate seconds={time.time() - t0:.1f}")
    return all_modes, all_scores


def weighted_choice(modes, scores, rng):
    modes = np.asarray(modes, dtype=float)
    if len(modes) == 0:
        return np.nan
    scores = np.asarray(scores, dtype=float)
    if len(scores) != len(modes) or not np.all(np.isfinite(scores)):
        return float(rng.choice(modes))
    w = np.exp(scores - np.nanmax(scores))
    w = w / np.sum(w) if np.sum(w) > 0 else np.ones(len(modes)) / len(modes)
    return float(rng.choice(modes, p=w))


def build_starts(all_modes=None, all_scores=None, n_mode_starts=FAST_MODE_STARTS, n_prior_starts=FAST_PRIOR_STARTS):
    starts = [("prior_mean", prior_mean.copy()), ("truth", truth_dist.copy())]
    if all_modes is not None:
        d_best1d = np.array([m[0] if len(m) else prior_mean[i] for i, m in enumerate(all_modes)], dtype=float)
        starts.append(("best_1d", d_best1d))
        for k in range(n_mode_starts):
            d = np.array([
                weighted_choice(all_modes[i], all_scores[i], rng) if len(all_modes[i]) else prior_mean[i]
                for i in range(N_PSR)
            ], dtype=float)
            starts.append((f"mode_combo_{k}", d))
    for k in range(n_prior_starts):
        starts.append((f"prior_draw_{k}", clipped_normal(prior_mean, sigmas, rng, MODE_PRIOR_SIGMA_WIDTH)))
    return starts


def make_full_optimizer(maxiter):
    return make_distance_optimizer(
        logl_fn, param_keys, base_values, disco_psrs,
        prior_means=prior_mean, prior_sigmas=sigmas,
        n_sigma=MODE_PRIOR_SIGMA_WIDTH, objective="logpost",
        maxiter=maxiter, factr=OPT_FACTR,
    )[0]


def run_gradient_optimizer_early(all_modes=None, all_scores=None):
    print("Diagnostic multi-start optimizer")
    t0 = time.time()
    starts = build_starts(all_modes, all_scores)
    slow_opt = make_full_optimizer(INITIAL_GRADIENT_MAXITER)
    fast_opt = None
    results = []
    use_fast = False
    for k, (label, d0) in enumerate(starts):
        if k >= 10 and len(results) >= 3:
            top3 = sorted(results, key=lambda r: r["logpost"], reverse=True)[:3]
            if max(r["logpost"] for r in top3) - min(r["logpost"] for r in top3) < 0.1:
                use_fast = True
        if use_fast and fast_opt is None:
            fast_opt = make_full_optimizer(FAST_GRADIENT_MAXITER)
            print("  top-3 logpost converged; remaining starts use fast maxiter")
        opt = fast_opt if use_fast else slow_opt
        d_opt, lnL, logpost, info = opt(d0)
        rec = dict(label=label, d_opt=d_opt, lnL=float(lnL), logpost=float(logpost), **info, **score_distances(d_opt))
        results.append(rec)
        print(f"  {k + 1:02d}/{len(starts)} {label}: rec={rec['n_recovered']}/{N_PSR}, logpost={logpost:.2f}, nit={info.get('nit')}")
    results.sort(key=lambda r: r["logpost"], reverse=True)
    print(f"gradient seconds={time.time() - t0:.1f}")
    return results


def display_score(label, distances):
    s = score_distances(distances)
    print(f"{label}: {s['n_recovered']}/{N_PSR} recovered, median_abs_modes={s['median_abs_modes']:.3g}, max_abs_modes={s['max_abs_modes']:.3g}")
    return s
```

## Sanity Check


```python
print("SANITY CHECK")
t0 = time.time()
summary = {}

display_score("truth", truth_dist)
display_score("prior_mean", prior_mean)
print(f"truth lnL={truth_lnL:.3f}, prior lnL={prior_mean_lnL:.3f}, delta={prior_mean_lnL - truth_lnL:.3f}")

all_pulsar_modes = all_pulsar_mode_scores = None
best_gradient_dist = prior_mean.copy()
gradient_results = []

if RUN_SANITY_GRADIENT:
    all_pulsar_modes, all_pulsar_mode_scores = build_empirical_candidates()
    gradient_results = run_gradient_optimizer_early(all_pulsar_modes, all_pulsar_mode_scores)
    best_gradient_dist = gradient_results[0]["d_opt"]
    display(pd.DataFrame(gradient_results[:10])[["label", "n_recovered", "lnL", "logpost", "nit", "success", "median_abs_modes", "max_abs_modes"]])
else:
    print("Sanity gradient skipped; using prior_mean as best_gradient_dist.")

summary["sanity"] = score_distances(best_gradient_dist)
print(f"sanity seconds={time.time() - t0:.1f}")
```

    SANITY CHECK
    truth: 40/40 recovered, median_abs_modes=0, max_abs_modes=0
    prior_mean: 0/40 recovered, median_abs_modes=108, max_abs_modes=1.1e+03
    truth lnL=158941.139, prior lnL=-1679351.629, delta=-1838292.769
    Building empirical 1D candidates
      candidates 5/40
      candidates 10/40
      candidates 15/40
      candidates 20/40
      candidates 25/40
      candidates 30/40
      candidates 35/40
      candidates 40/40
    candidate seconds=111.9
    Diagnostic multi-start optimizer
      01/28 prior_mean: rec=0/40, logpost=-769468.44, nit=100
      02/28 truth: rec=40/40, logpost=158953.96, nit=33
      03/28 best_1d: rec=40/40, logpost=158953.96, nit=33
      04/28 mode_combo_0: rec=40/40, logpost=158953.96, nit=33
      05/28 mode_combo_1: rec=40/40, logpost=158953.96, nit=33
      06/28 mode_combo_2: rec=40/40, logpost=158953.96, nit=33
      07/28 mode_combo_3: rec=40/40, logpost=158953.96, nit=33
      08/28 mode_combo_4: rec=40/40, logpost=158953.96, nit=33
      09/28 mode_combo_5: rec=40/40, logpost=158953.96, nit=33
      10/28 mode_combo_6: rec=40/40, logpost=158953.96, nit=33
      top-3 logpost converged; remaining starts use fast maxiter
      11/28 mode_combo_7: rec=40/40, logpost=158953.96, nit=33
      12/28 mode_combo_8: rec=40/40, logpost=158953.96, nit=33
      13/28 mode_combo_9: rec=40/40, logpost=158953.96, nit=33
      14/28 mode_combo_10: rec=40/40, logpost=158953.96, nit=33
      15/28 mode_combo_11: rec=40/40, logpost=158953.96, nit=33
      16/28 mode_combo_12: rec=40/40, logpost=158953.96, nit=33
      17/28 mode_combo_13: rec=40/40, logpost=158953.96, nit=33
      18/28 mode_combo_14: rec=40/40, logpost=158953.96, nit=33
      19/28 mode_combo_15: rec=40/40, logpost=158953.96, nit=33
      20/28 mode_combo_16: rec=40/40, logpost=158953.96, nit=33
      21/28 mode_combo_17: rec=40/40, logpost=158953.96, nit=33
      22/28 mode_combo_18: rec=40/40, logpost=158953.96, nit=33
      23/28 mode_combo_19: rec=40/40, logpost=158953.96, nit=33
      24/28 prior_draw_0: rec=0/40, logpost=-738086.53, nit=100
      25/28 prior_draw_1: rec=0/40, logpost=-962890.32, nit=100
      26/28 prior_draw_2: rec=1/40, logpost=-709229.05, nit=100
      27/28 prior_draw_3: rec=0/40, logpost=-722298.27, nit=100
      28/28 prior_draw_4: rec=2/40, logpost=-537523.41, nit=100
    gradient seconds=123.8



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>label</th>
      <th>n_recovered</th>
      <th>lnL</th>
      <th>logpost</th>
      <th>nit</th>
      <th>success</th>
      <th>median_abs_modes</th>
      <th>max_abs_modes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>truth</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>1</th>
      <td>best_1d</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>2</th>
      <td>mode_combo_0</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>3</th>
      <td>mode_combo_1</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>4</th>
      <td>mode_combo_2</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>5</th>
      <td>mode_combo_3</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>6</th>
      <td>mode_combo_4</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>7</th>
      <td>mode_combo_5</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>8</th>
      <td>mode_combo_6</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
    <tr>
      <th>9</th>
      <td>mode_combo_7</td>
      <td>40</td>
      <td>158978.702974</td>
      <td>158953.963557</td>
      <td>33</td>
      <td>True</td>
      <td>0.001769</td>
      <td>0.024192</td>
    </tr>
  </tbody>
</table>
</div>


    sanity seconds=236.1


## Experiment 1 - Likelihood Surface Characterization


```python
print("EXPERIMENT 1 - LIKELIHOOD SURFACE CHARACTERIZATION")
t0 = time.time()
exp1_df = pd.DataFrame()

if RUN_EXP1:
    rows = []
    for i in REPRESENTATIVE_INDICES:
        truth_ctx = truth_dist.copy()
        prior_ctx = prior_mean.copy()
        x_t, y_t = conditional_scan(i, truth_ctx, n_points=DIAGNOSTIC_SCAN_POINTS)
        x_p, y_p = conditional_scan(i, prior_ctx, n_points=DIAGNOSTIC_SCAN_POINTS)
        mt = scan_metrics(i, x_t, y_t)
        mp = scan_metrics(i, x_p, y_p)
        rows.append(dict(
            idx=i,
            pulsar=disco_psrs[i].name,
            sigma_over_mode=float(hardness[i]),
            contrast_at_truth_context=mt["contrast"],
            contrast_at_prior_context=mp["contrast"],
            n_degenerate_at_truth=mt["n_degenerate"],
            n_degenerate_at_prior=mp["n_degenerate"],
            truth_is_global_max=mt["truth_is_global_max"],
            prior_context_truth_global=mp["truth_is_global_max"],
        ))
        print(f"  {disco_psrs[i].name}: truth contrast={mt['contrast']:.2f}, prior contrast={mp['contrast']:.2f}")
    exp1_df = pd.DataFrame(rows)
    display(exp1_df)
else:
    print("Skipped")

summary["exp1"] = exp1_df
print(f"exp1 seconds={time.time() - t0:.1f}")
```

    EXPERIMENT 1 - LIKELIHOOD SURFACE CHARACTERIZATION


      J0437-4715: truth contrast=50543.16, prior contrast=47785.89
      J1022+1001: truth contrast=1112.40, prior contrast=906.83
      J1630+3734: truth contrast=2508.38, prior contrast=2117.04
      J1640+2224: truth contrast=219622.85, prior contrast=210149.09
      J1312+0051: truth contrast=190154.20, prior contrast=215684.30



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>idx</th>
      <th>pulsar</th>
      <th>sigma_over_mode</th>
      <th>contrast_at_truth_context</th>
      <th>contrast_at_prior_context</th>
      <th>n_degenerate_at_truth</th>
      <th>n_degenerate_at_prior</th>
      <th>truth_is_global_max</th>
      <th>prior_context_truth_global</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>J0437-4715</td>
      <td>0.863615</td>
      <td>50543.163430</td>
      <td>47785.888941</td>
      <td>1</td>
      <td>1</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>4</td>
      <td>J1022+1001</td>
      <td>79.269657</td>
      <td>1112.396214</td>
      <td>906.831361</td>
      <td>1</td>
      <td>2</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>10</td>
      <td>J1630+3734</td>
      <td>98.000816</td>
      <td>2508.380040</td>
      <td>2117.038519</td>
      <td>1</td>
      <td>5</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>20</td>
      <td>J1640+2224</td>
      <td>376.555290</td>
      <td>219622.849515</td>
      <td>210149.088319</td>
      <td>1</td>
      <td>1</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>34</td>
      <td>J1312+0051</td>
      <td>777.087703</td>
      <td>190154.204693</td>
      <td>215684.296373</td>
      <td>1</td>
      <td>1</td>
      <td>True</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>


    exp1 seconds=2.3


## Experiment 2 - Iterative Conditional Mode Recovery


```python
print("EXPERIMENT 2 - ITERATIVE CONDITIONAL MODE RECOVERY")
t0 = time.time()
exp2_scores = {}

def iterative_conditional(order, label):
    current = prior_mean.copy()
    rows = []
    for step, i in enumerate(order, 1):
        peak, peak_lnL = choose_peak_distance(i, current, n_points=DIAGNOSTIC_SCAN_POINTS)
        current[i] = peak
        recovered = abs((current[i] - truth_dist[i]) / mode_spacings[i]) < 0.5
        rows.append(dict(step=step, idx=i, pulsar=disco_psrs[i].name, recovered=bool(recovered), err_modes=(current[i] - truth_dist[i]) / mode_spacings[i]))
        if step % 10 == 0 or step == len(order):
            s = score_distances(current)
            print(f"  {label} step {step}: {s['n_recovered']}/{N_PSR}")
    return current, pd.DataFrame(rows)

if RUN_EXP2:
    easy_order = list(np.argsort(hardness))
    hard_order = list(np.argsort(hardness)[::-1])
    d_easy, exp2_easy_steps = iterative_conditional(easy_order, "easy-first")
    d_hard, exp2_hard_steps = iterative_conditional(hard_order, "hard-first")
    exp2_scores["easy_first"] = display_score("iterative easy-first", d_easy)
    exp2_scores["hard_first"] = display_score("iterative hard-first", d_hard)
else:
    print("Skipped")

summary["exp2"] = exp2_scores
print(f"exp2 seconds={time.time() - t0:.1f}")
```

    EXPERIMENT 2 - ITERATIVE CONDITIONAL MODE RECOVERY
      easy-first step 10: 8/40
      easy-first step 20: 14/40
      easy-first step 30: 24/40
      easy-first step 40: 33/40
      hard-first step 10: 7/40
      hard-first step 20: 13/40
      hard-first step 30: 21/40
      hard-first step 40: 27/40
    iterative easy-first: 33/40 recovered, median_abs_modes=0, max_abs_modes=472
    iterative hard-first: 27/40 recovered, median_abs_modes=0, max_abs_modes=742
    exp2 seconds=18.1


## Experiment 3 - Truth-Seeded Progressive Corruption


```python
print("EXPERIMENT 3 - TRUTH-SEEDED PROGRESSIVE CORRUPTION")
t0 = time.time()
exp3_df = pd.DataFrame()

if RUN_EXP3:
    opt = make_full_optimizer(FAST_GRADIENT_MAXITER)
    rows = []
    hard_idxs = list(np.argsort(hardness)[::-1][:PROGRESSIVE_N_HARD])
    for step, i in enumerate(hard_idxs, 1):
        start = truth_dist.copy()
        start[i] = prior_mean[i]
        d_opt, lnL, logpost, info = opt(start)
        sc = score_distances(d_opt)
        err_modes = (d_opt - truth_dist) / mode_spacings
        cascaded = np.where(np.abs(err_modes) >= 0.5)[0].tolist()
        rows.append(dict(
            step=step, idx=i, pulsar=disco_psrs[i].name,
            sigma_over_mode=float(hardness[i]),
            recovered=sc["n_recovered"], lnL=float(lnL), logpost=float(logpost),
            nit=info.get("nit"), cascaded=cascaded,
        ))
        print(f"  corrupt {disco_psrs[i].name}: recovered={sc['n_recovered']}/{N_PSR}, cascaded={len(cascaded)}")
    exp3_df = pd.DataFrame(rows)
    display(exp3_df[["step", "pulsar", "sigma_over_mode", "recovered", "nit", "cascaded"]])
else:
    print("Skipped")

summary["exp3"] = exp3_df
print(f"exp3 seconds={time.time() - t0:.1f}")
```

    EXPERIMENT 3 - TRUTH-SEEDED PROGRESSIVE CORRUPTION
      corrupt J1312+0051: recovered=39/40, cascaded=1
      corrupt J1125+7819: recovered=39/40, cascaded=1
      corrupt J0900-3144: recovered=39/40, cascaded=1
      corrupt J1853+1303: recovered=39/40, cascaded=1
      corrupt J1658-5324: recovered=39/40, cascaded=1
      corrupt J1756-2251: recovered=39/40, cascaded=1
      corrupt J1640+2224: recovered=39/40, cascaded=1
      corrupt J0605+3757: recovered=39/40, cascaded=1
      corrupt J1738+0333: recovered=39/40, cascaded=1
      corrupt B1855+09: recovered=39/40, cascaded=1



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>step</th>
      <th>pulsar</th>
      <th>sigma_over_mode</th>
      <th>recovered</th>
      <th>nit</th>
      <th>cascaded</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>J1312+0051</td>
      <td>777.087703</td>
      <td>39</td>
      <td>85</td>
      <td>[34]</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>J1125+7819</td>
      <td>489.689239</td>
      <td>39</td>
      <td>56</td>
      <td>[28]</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>J0900-3144</td>
      <td>487.546564</td>
      <td>39</td>
      <td>38</td>
      <td>[39]</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>J1853+1303</td>
      <td>464.040500</td>
      <td>39</td>
      <td>66</td>
      <td>[35]</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>J1658-5324</td>
      <td>399.784943</td>
      <td>39</td>
      <td>100</td>
      <td>[38]</td>
    </tr>
    <tr>
      <th>5</th>
      <td>6</td>
      <td>J1756-2251</td>
      <td>380.165004</td>
      <td>39</td>
      <td>84</td>
      <td>[33]</td>
    </tr>
    <tr>
      <th>6</th>
      <td>7</td>
      <td>J1640+2224</td>
      <td>376.555290</td>
      <td>39</td>
      <td>100</td>
      <td>[20]</td>
    </tr>
    <tr>
      <th>7</th>
      <td>8</td>
      <td>J0605+3757</td>
      <td>357.287058</td>
      <td>39</td>
      <td>100</td>
      <td>[31]</td>
    </tr>
    <tr>
      <th>8</th>
      <td>9</td>
      <td>J1738+0333</td>
      <td>325.708119</td>
      <td>39</td>
      <td>48</td>
      <td>[22]</td>
    </tr>
    <tr>
      <th>9</th>
      <td>10</td>
      <td>B1855+09</td>
      <td>319.289801</td>
      <td>39</td>
      <td>66</td>
      <td>[26]</td>
    </tr>
  </tbody>
</table>
</div>


    exp3 seconds=51.6


## Experiment 4 - Block Coordinate Descent


```python
print("EXPERIMENT 4 - BLOCK COORDINATE DESCENT")
t0 = time.time()
exp4_df = pd.DataFrame()

# Cache one compiled objective per block. The current full distance vector is
# passed as a runtime argument, avoiding a fresh JAX compile for every update.
_block_optimizer_cache = {}


def get_block_optimizer(block_indices, maxiter=FAST_GRADIENT_MAXITER):
    block_indices = tuple(int(i) for i in block_indices)
    cache_key = (block_indices, int(maxiter))
    if cache_key in _block_optimizer_cache:
        return _block_optimizer_cache[cache_key]

    block_keys = [distance_key(disco_psrs[i]) for i in block_indices]
    block_param_indices_np = np.array([param_keys.index(k) for k in block_keys], dtype=int)
    block_param_indices = jnp.array(block_param_indices_np)
    mu = jnp.asarray(prior_mean[list(block_indices)], dtype=jnp.float64)
    sig = jnp.asarray(sigmas[list(block_indices)], dtype=jnp.float64)
    z_lo = jnp.maximum(jnp.full(len(block_indices), -MODE_PRIOR_SIGMA_WIDTH), (0.01 - mu) / sig)
    z_hi = jnp.full(len(block_indices), MODE_PRIOR_SIGMA_WIDTH)
    z_lo_np = np.asarray(z_lo, dtype=float)
    z_hi_np = np.asarray(z_hi, dtype=float)
    ftol = OPT_FACTR * np.finfo(float).eps

    @jax.jit
    def neg_obj(z, fixed_values):
        d = mu + sig * z
        x = fixed_values.at[block_param_indices].set(d)
        return -(logl_fn(x) - 0.5 * jnp.sum(z ** 2))

    grad = jax.jit(jax.grad(neg_obj, argnums=0))

    def opt(d_start):
        d_start = np.asarray(d_start, dtype=float)
        fixed_values = jnp.asarray(set_distances(base_values, param_keys, disco_psrs, d_start), dtype=jnp.float64)
        z0 = (d_start[list(block_indices)] - np.asarray(mu)) / np.asarray(sig)
        z0 = np.clip(z0, z_lo_np + 1e-10, z_hi_np - 1e-10)

        def fun(z):
            return float(neg_obj(jnp.asarray(z, dtype=jnp.float64), fixed_values))

        def jac(z):
            return np.asarray(grad(jnp.asarray(z, dtype=jnp.float64), fixed_values), dtype=float)

        res = scipy.optimize.minimize(
            fun, z0, jac=jac, method="L-BFGS-B",
            bounds=list(zip(z_lo_np, z_hi_np)),
            options={"maxiter": maxiter, "ftol": ftol, "gtol": 1e-10},
        )
        out = d_start.copy()
        out[list(block_indices)] = np.asarray(mu) + np.asarray(sig) * res.x
        return out, res

    _block_optimizer_cache[cache_key] = opt
    return opt


def random_blocks():
    idx = np.arange(N_PSR)
    rng.shuffle(idx)
    return [idx[i:i + BLOCK_SIZE].tolist() for i in range(0, N_PSR, BLOCK_SIZE)]


def sorted_blocks():
    idx = np.argsort(hardness)
    return [idx[i:i + BLOCK_SIZE].tolist() for i in range(0, N_PSR, BLOCK_SIZE)]


def sky_blocks():
    remaining = set(range(N_PSR))
    blocks = []
    pos = np.array([p.pos for p in disco_psrs])
    while remaining:
        seed = min(remaining)
        rem = np.array(sorted(remaining))
        dots = pos[rem] @ pos[seed]
        group = rem[np.argsort(dots)[::-1][:BLOCK_SIZE]].tolist()
        blocks.append(group)
        remaining -= set(group)
    return blocks


def run_block_cd(blocks, label):
    current = prior_mean.copy()
    for cycle in range(BLOCK_CYCLES):
        for b, block in enumerate(blocks, 1):
            opt = get_block_optimizer(block)
            current, res = opt(current)
            if b == 1 or b == len(blocks):
                sc = score_distances(current)
                print(f"    {label} cycle {cycle + 1} block {b}/{len(blocks)}: {sc['n_recovered']}/{N_PSR}, nit={res.nit}")
        sc = score_distances(current)
        print(f"  {label} cycle {cycle + 1}: {sc['n_recovered']}/{N_PSR}")
    return current


if RUN_EXP4:
    all_block_sets = {
        "random_8": random_blocks,
        "sky_8": sky_blocks,
        "sigma_sorted_8": sorted_blocks,
    }
    requested = globals().get("EXP4_STRATEGIES", list(all_block_sets))
    block_sets = {name: all_block_sets[name]() for name in requested}
    rows = []
    block_results = {}
    for label, blocks in block_sets.items():
        print(f"  strategy {label}: {len(blocks)} blocks, {BLOCK_CYCLES} cycle(s)")
        d = run_block_cd(blocks, label)
        block_results[label] = d
        rows.append(dict(strategy=label, **score_distances(d)))
    exp4_df = pd.DataFrame(rows).sort_values("n_recovered", ascending=False)
    display(exp4_df)
else:
    print("Skipped")

summary["exp4"] = exp4_df
print(f"exp4 seconds={time.time() - t0:.1f}")
```

    EXPERIMENT 4 - BLOCK COORDINATE DESCENT
      strategy sigma_sorted_8: 5 blocks, 1 cycle(s)
        sigma_sorted_8 cycle 1 block 1/5: 0/40, nit=50
        sigma_sorted_8 cycle 1 block 5/5: 0/40, nit=17
      sigma_sorted_8 cycle 1: 0/40



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>strategy</th>
      <th>n_recovered</th>
      <th>frac_recovered</th>
      <th>median_abs_modes</th>
      <th>max_abs_modes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>sigma_sorted_8</td>
      <td>0</td>
      <td>0.0</td>
      <td>210.710489</td>
      <td>1100.838699</td>
    </tr>
  </tbody>
</table>
</div>


    exp4 seconds=197.1


## Experiment 5 - Failure Diagnosis For Stuck Pulsars


```python
print("EXPERIMENT 5 - FAILURE DIAGNOSIS FOR STUCK PULSARS")
t0 = time.time()
exp5_df = pd.DataFrame()

if RUN_EXP5:
    err_best = (best_gradient_dist - truth_dist) / mode_spacings
    if FAILED_PULSAR_INDICES is None:
        failed = np.where(np.abs(err_best) >= 0.5)[0].tolist()
        failed_source = "best non-truth diagnostic solution"
    else:
        failed = [int(i) for i in FAILED_PULSAR_INDICES]
        failed_source = "FAILED_PULSAR_INDICES override"
    print(f"failed pulsars from {failed_source}: {len(failed)}")
    if len(failed) == 0:
        print("No failed pulsars to diagnose. Use FAILED_PULSAR_INDICES from notebook 05, or rerun sanity without truth-seeded dominance.")
    rows = []
    for n, i in enumerate(failed, 1):
        x_t, y_t = conditional_scan(i, truth_dist, n_points=DIAGNOSTIC_SCAN_POINTS)
        x_b, y_b = conditional_scan(i, best_gradient_dist, n_points=DIAGNOSTIC_SCAN_POINTS)
        mt = scan_metrics(i, x_t, y_t)
        mb = scan_metrics(i, x_b, y_b)
        joint_failure = mt["truth_is_global_max"]
        locked_in = not mb["truth_is_global_max"]
        rows.append(dict(
            idx=i,
            pulsar=disco_psrs[i].name,
            err_modes=float(err_best[i]),
            truth_context_truth_global=mt["truth_is_global_max"],
            best_context_truth_global=mb["truth_is_global_max"],
            contrast_truth_context=mt["contrast"],
            contrast_best_context=mb["contrast"],
            joint_optimization_failure=bool(joint_failure),
            locked_in=bool(locked_in),
        ))
        print(f"  {n}/{len(failed)} {disco_psrs[i].name}: joint={joint_failure}, locked={locked_in}")
    exp5_df = pd.DataFrame(rows)
    display(exp5_df)
else:
    print("Skipped")

summary["exp5"] = exp5_df
print(f"exp5 seconds={time.time() - t0:.1f}")
```

    EXPERIMENT 5 - FAILURE DIAGNOSIS FOR STUCK PULSARS
    failed pulsars from best non-truth diagnostic solution: 0
    No failed pulsars to diagnose. Use FAILED_PULSAR_INDICES from notebook 05, or rerun sanity without truth-seeded dominance.



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>


    exp5 seconds=0.0


## Investigation Summary


```python
print("INVESTIGATION SUMMARY")
print("=====================")

if isinstance(summary.get("exp1"), pd.DataFrame) and len(summary["exp1"]):
    avg_deg = float(summary["exp1"]["n_degenerate_at_truth"].mean())
    avg_contrast = float(summary["exp1"]["contrast_at_truth_context"].mean())
    if avg_contrast > 10:
        contrast_label = "sharp"
    elif avg_contrast > 2:
        contrast_label = "moderate"
    else:
        contrast_label = "flat"
    print(f"Exp 1 - Mode contrast: {contrast_label} at truth context, {avg_deg:.1f} degenerate modes avg")
else:
    print("Exp 1 - Mode contrast: skipped")

if summary.get("exp2"):
    e = summary["exp2"].get("easy_first", {}).get("n_recovered", 0)
    h = summary["exp2"].get("hard_first", {}).get("n_recovered", 0)
    print(f"Exp 2 - Iterative 1D: {e}/{N_PSR} recovered (easy-first), {h}/{N_PSR} (hard-first)")
else:
    print("Exp 2 - Iterative 1D: skipped")

if isinstance(summary.get("exp3"), pd.DataFrame) and len(summary["exp3"]):
    frag = summary["exp3"].sort_values("recovered").head(5)["pulsar"].tolist()
    print(f"Exp 3 - Fragile pulsars: {frag}")
else:
    print("Exp 3 - Fragile pulsars: skipped")

if isinstance(summary.get("exp4"), pd.DataFrame) and len(summary["exp4"]):
    best = summary["exp4"].sort_values("n_recovered", ascending=False).iloc[0]
    print(f"Exp 4 - Block CD: {int(best['n_recovered'])}/{N_PSR} recovered (best block strategy: {best['strategy']})")
else:
    print("Exp 4 - Block CD: skipped")

if isinstance(summary.get("exp5"), pd.DataFrame) and len(summary["exp5"]):
    n_joint = int(summary["exp5"]["joint_optimization_failure"].sum())
    n_locked = int(summary["exp5"]["locked_in"].sum())
    n_fail = len(summary["exp5"])
    print(f"Exp 5 - Failure diagnosis: {n_joint}/{n_fail} are joint-optimization failures, {n_locked}/{n_fail} are locked-in")
else:
    print("Exp 5 - Failure diagnosis: skipped/no failed pulsars")
```

    INVESTIGATION SUMMARY
    =====================
    Exp 1 - Mode contrast: sharp at truth context, 1.0 degenerate modes avg
    Exp 2 - Iterative 1D: 33/40 recovered (easy-first), 27/40 (hard-first)
    Exp 3 - Fragile pulsars: ['J1312+0051', 'J1125+7819', 'J0900-3144', 'J1853+1303', 'J1658-5324']
    Exp 4 - Block CD: 0/40 recovered (best block strategy: sigma_sorted_8)
    Exp 5 - Failure diagnosis: skipped/no failed pulsars

