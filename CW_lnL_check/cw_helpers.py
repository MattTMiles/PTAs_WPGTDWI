"""
Shared helper functions for CW likelihood distance-dependence checks.

Uses enterprise for data simulation and discovery for likelihood evaluation,
following the pipeline from HSYMT/data_likelihood_sandbox.ipynb.
"""

import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
os.environ["JAX_DISABLE_MMAP_CACHE"] = "1"

import numpy as np
import jax
jax.config.update('jax_enable_x64', True)
import jax.numpy as jnp
import functools

import discovery as ds
from discovery.deterministic import make_phase_connected_binary

import glob
import scipy.linalg as sl
import scipy.sparse as ss
from sksparse.cholmod import cholesky
from scipy.signal import find_peaks

from enterprise.signals import signal_base, selections, utils
import enterprise.signals.parameter as parameter
from enterprise.signals import white_signals, gp_signals
from enterprise.signals.gp_signals import get_timing_model_basis, BasisGP
from enterprise.signals.parameter import function
from enterprise_extensions.blocks import common_red_noise_block
from enterprise_extensions import deterministic as ent_det
from enterprise_extensions import load_feathers

FEATHER_DIR = "/home/mattm/projects/HSYMT/data_products/"
KPC_OVER_C = 1.0267e11  # kpc in light-seconds


# ============================================================
# Enterprise: simulation infrastructure
# ============================================================

@function
def _tm_prior(weights, toas, variance=1e40):
    return weights * variance * len(toas)


def _TimingModel(prior_variance=1e-14):
    basis = get_timing_model_basis(False, True)
    prior = _tm_prior(variance=prior_variance)
    BaseClass = BasisGP(prior, basis, coefficients=False, name="linear_timing_model")

    class TM(BaseClass):
        signal_type = "basis"
        signal_name = "linear timing model"
        signal_id = "linear_timing_model"

    return TM


def simulate(pta, params, sparse_cholesky=True):
    """Simulate residuals using enterprise PTA model."""
    delays, ndiags, fmats, phis = (
        pta.get_delay(params=params),
        pta.get_ndiag(params=params),
        pta.get_basis(params=params),
        pta.get_phi(params=params),
    )

    gpresiduals = []
    if pta._commonsignals:
        if sparse_cholesky:
            cf = cholesky(ss.csc_matrix(phis))
            gp = np.zeros(phis.shape[0])
            gp[cf.P()] = np.dot(cf.L().toarray(), np.random.randn(phis.shape[0]))
        else:
            gp = np.dot(sl.cholesky(phis, lower=True), np.random.randn(phis.shape[0]))

        i = 0
        for fmat in fmats:
            j = i + fmat.shape[1]
            gpresiduals.append(np.dot(fmat, gp[i:j]))
            i = j
    else:
        for fmat, phi in zip(fmats, phis):
            if phi is None:
                gpresiduals.append(0)
            elif phi.ndim == 1:
                gpresiduals.append(np.dot(fmat, np.sqrt(phi) * np.random.randn(phi.shape[0])))
            else:
                raise NotImplementedError

    whiteresiduals = []
    for delay, ndiag in zip(delays, ndiags):
        if ndiag is None:
            whiteresiduals.append(0)
        elif isinstance(ndiag, signal_base.ShermanMorrison):
            n = np.diag(ndiag._nvec)
            for j, s in zip(ndiag._jvec, ndiag._slices):
                n[s, s] += j
            whiteresiduals.append(
                delay + np.dot(sl.cholesky(n, lower=True), np.random.randn(n.shape[0]))
            )
        elif ndiag.ndim == 1:
            whiteresiduals.append(delay + np.sqrt(ndiag) * np.random.randn(ndiag.shape[0]))
        else:
            raise NotImplementedError

    return [np.array(g + w) for g, w in zip(gpresiduals, whiteresiduals)]


def build_enterprise_pta(
    psrs, num_cw, components=30, log10_equad=-8.0,
    include_rn=True, rn_components=30,
):
    """Build enterprise PTA model for simulation.

    Parameters
    ----------
    include_rn : bool
        If True, add per-pulsar achromatic power-law red noise.
    rn_components : int
        Number of Fourier components for the intrinsic red-noise basis.

    Returns (pta, cw_block_names, Tspan).
    """
    for psr in psrs:
        psr._pdist = psr.pdist
        psr.residuals = np.array(psr.toas) * 0.0

    tmin = [p.toas.min() for p in psrs]
    tmax = [p.toas.max() for p in psrs]
    Tspan = np.max(tmax) - np.min(tmin)

    selection = selections.Selection(selections.by_backend)
    efac = parameter.Constant(1)
    equad = parameter.Constant(log10_equad)

    log10_A_gw = parameter.Uniform(-18, -11)("gwb_log10_A")
    gamma_gw = parameter.Uniform(0, 7)("gwb_gamma")

    cw_block_names = ["cw" if i == 0 else f"cw{i+1}" for i in range(num_cw)]

    tm = _TimingModel()
    models = []
    for p in psrs:
        s = tm
        s += white_signals.MeasurementNoise(efac=efac, selection=selection)
        s += white_signals.TNEquadNoise(log10_tnequad=equad, selection=selection)

        if include_rn:
            log10_A_rn = parameter.Uniform(-20, -11)
            gamma_rn = parameter.Uniform(0, 7)
            pl_rn = utils.powerlaw(log10_A=log10_A_rn, gamma=gamma_rn)
            rn = gp_signals.FourierBasisGP(
                spectrum=pl_rn, components=rn_components, Tspan=Tspan, name="rednoise"
            )
            s += rn
        s += common_red_noise_block(
            psd="powerlaw", prior="log-uniform", components=components, orf="hd", name="gwb"
        )

        for name in cw_block_names:
            s += ent_det.cw_block_circ(
                amp_prior="log-uniform",
                dist_prior=None,
                skyloc=None,
                log10_fgw=None,
                psrTerm=True,
                phase_connected=True,
                discoclone=False,
                name=name,
                evolve=True,
            )
        models.append(s(p))

    pta = signal_base.PTA(models)
    pta.set_default_params({})
    return pta, cw_block_names, Tspan


# ============================================================
# Injection parameter generation
# ============================================================

def generate_injection_params(
    pta, psrs, num_cw, cw_block_names, log10_h,
    scenario="single", rng=None,
    log10_h_range=None, log10_fgw_range=None, log10_mc_range=None,
    gwb_log10_A=None, gwb_gamma=None, include_rn=True,
):
    """Generate a complete injection parameter dictionary.

    Parameters
    ----------
    scenario : str
        "single"         -- one CW (num_cw should be 1)
        "close_freq"     -- all CWs within +/-0.1 dex of f_gw = 10^-8 Hz
        "close_sky"      -- all CWs within ~15 deg on sky
        "well_separated" -- random sky + frequency
        "realistic"      -- per-source draws of log10_h, log10_fgw, log10_mc
                            from supplied ranges (tuples); random sky.
    log10_h_range, log10_fgw_range, log10_mc_range : tuple or None
        (lo, hi) ranges for per-source draws in scenario="realistic".
    gwb_log10_A, gwb_gamma : float or None
        Override GWB amplitude / spectral index (default log10_A=-17.5, gamma=4.333).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    params = {}

    # Per-pulsar intrinsic red noise -- realistic MSP range.
    # Keep discovery-compatible enterprise name "rednoise".
    if include_rn:
        params.update(
            {p: rng.uniform(-15, -13) for p in pta.param_names if "rednoise_log10_A" in p}
        )
        params.update(
            {p: rng.uniform(1.5, 5.0) for p in pta.param_names if "rednoise_gamma" in p}
        )

    # Pulsar distances from EM priors
    for psr in psrs:
        params.update(
            {p: psr.pdist[0] for p in pta.param_names if psr.name + "_cw_p_dist" in p}
        )

    # GWB
    params.update({
        "gwb_gamma": 4.333 if gwb_gamma is None else gwb_gamma,
        "gwb_log10_A": -17.5 if gwb_log10_A is None else gwb_log10_A,
    })

    # CW sources
    if scenario == "close_freq":
        base_freq = -8.0
        for name in cw_block_names:
            params.update({
                f"{name}_cos_inc": rng.uniform(-1.0, 1.0),
                f"{name}_log10_fgw": base_freq + rng.uniform(-0.1, 0.1),
                f"{name}_log10_h": log10_h,
                f"{name}_log10_mc": 9.0,
                f"{name}_gwphi": rng.uniform(0, 2 * np.pi),
                f"{name}_cos_gwtheta": rng.uniform(-1, 1),
                f"{name}_phase0": rng.uniform(0, 2 * np.pi),
                f"{name}_psi": rng.uniform(0, np.pi),
            })

    elif scenario == "close_sky":
        # ~6 degree diameter cluster (Fornax-like):
        # At cos_gwtheta=0.3 (theta~72.5 deg, sin~0.95),
        # 3 deg half-angle -> delta_cos ~ 0.05, delta_phi ~ 0.055 rad
        base_cos_gwtheta = 0.3
        base_gwphi = 2.5
        for name in cw_block_names:
            params.update({
                f"{name}_cos_inc": rng.uniform(-1.0, 1.0),
                f"{name}_log10_fgw": rng.uniform(-8.5, -7.5),
                f"{name}_log10_h": log10_h,
                f"{name}_log10_mc": 9.0,
                f"{name}_gwphi": base_gwphi + rng.uniform(-0.055, 0.055),
                f"{name}_cos_gwtheta": np.clip(
                    base_cos_gwtheta + rng.uniform(-0.05, 0.05), -1, 1
                ),
                f"{name}_phase0": rng.uniform(0, 2 * np.pi),
                f"{name}_psi": rng.uniform(0, np.pi),
            })

    elif scenario == "well_separated":
        for name in cw_block_names:
            params.update({
                f"{name}_cos_inc": rng.uniform(-1.0, 1.0),
                f"{name}_log10_fgw": rng.uniform(-8.5, -7.5),
                f"{name}_log10_h": log10_h,
                f"{name}_log10_mc": 9.0,
                f"{name}_gwphi": rng.uniform(0, 2 * np.pi),
                f"{name}_cos_gwtheta": rng.uniform(-1, 1),
                f"{name}_phase0": rng.uniform(0, 2 * np.pi),
                f"{name}_psi": rng.uniform(0, np.pi),
            })

    elif scenario == "realistic":
        h_lo, h_hi = log10_h_range if log10_h_range is not None else (-14.3, -13.8)
        f_lo, f_hi = log10_fgw_range if log10_fgw_range is not None else (-8.3, -7.7)
        mc_lo, mc_hi = log10_mc_range if log10_mc_range is not None else (8.5, 9.5)
        for name in cw_block_names:
            params.update({
                f"{name}_cos_inc": rng.uniform(-1.0, 1.0),
                f"{name}_log10_fgw": rng.uniform(f_lo, f_hi),
                f"{name}_log10_h": rng.uniform(h_lo, h_hi),
                f"{name}_log10_mc": rng.uniform(mc_lo, mc_hi),
                f"{name}_gwphi": rng.uniform(0, 2 * np.pi),
                f"{name}_cos_gwtheta": rng.uniform(-1, 1),
                f"{name}_phase0": rng.uniform(0, 2 * np.pi),
                f"{name}_psi": rng.uniform(0, np.pi),
            })

    else:  # "single" or default
        for name in cw_block_names:
            params.update({
                f"{name}_cos_inc": rng.uniform(-1.0, 1.0),
                f"{name}_log10_fgw": -8.0,
                f"{name}_log10_h": log10_h,
                f"{name}_log10_mc": 9.0,
                f"{name}_gwphi": rng.uniform(0, 2 * np.pi),
                f"{name}_cos_gwtheta": rng.uniform(-1, 1),
                f"{name}_phase0": rng.uniform(0, 2 * np.pi),
                f"{name}_psi": rng.uniform(0, np.pi),
            })

    return params


# ============================================================
# Discovery: likelihood construction
# ============================================================

class MultiSourceDelay:
    """CW delay model for N sources sharing a single pulsar distance."""

    global_params = (
        "cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
        "log10_fgw", "log10_h", "phase0", "psi",
    )

    def __init__(self, psr, n_sources, include_pterm):
        self.psr = psr
        self.n_sources = n_sources
        self.include_pterm = include_pterm
        self.single_source = make_phase_connected_binary(pulsarterm=include_pterm)
        self._suffixes = ["" if i == 0 else f"_{i+1}" for i in range(n_sources)]

        # vmap over CW params (axes 2-9), keep toas/pos/p_dist fixed
        self._vmapped_source = jax.vmap(
            self.single_source,
            in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0, None),
        )

        params = []
        for suffix in self._suffixes:
            for key in self.global_params:
                name = f"cw_{key}{suffix}"
                if name not in params:
                    params.append(name)
        if include_pterm:
            params.append(f"{psr.name}_cw_p_dist")
        self._pdist_key = f"{psr.name}_cw_p_dist"
        self.params = params

    def __call__(self, params):
        # Stack CW params into arrays of shape (n_sources,)
        cw_vals = {
            key: jnp.stack([params[f"cw_{key}{suffix}"] for suffix in self._suffixes])
            for key in self.global_params
        }
        p_dist = params[self._pdist_key] if self.include_pterm else 1.0

        all_delays = self._vmapped_source(
            self.psr.toas, self.psr.pos,
            cw_vals["cos_gwtheta"], cw_vals["gwphi"],
            cw_vals["cos_inc"], cw_vals["log10_mc"],
            cw_vals["log10_fgw"], cw_vals["log10_h"],
            cw_vals["phase0"], cw_vals["psi"],
            p_dist,
        )
        return jnp.sum(all_delays, axis=0)


def _enterprise_to_disco_params(enterprise_params, cw_block_names):
    """Map enterprise CW parameter names to discovery naming convention."""
    disco_suffixes = ["" if i == 0 else f"_{i+1}" for i in range(len(cw_block_names))]
    temp_dict = dict(enterprise_params)
    for name, suffix in zip(cw_block_names, disco_suffixes):
        for key, val in enterprise_params.items():
            if key.startswith(f"{name}_"):
                base = key.replace(name, "cw", 1)
                disco_key = f"{base}{suffix}"
                temp_dict[disco_key] = val
    return temp_dict


def build_disco_likelihood(
    disco_psrs, residual_map, num_cw, enterprise_params, cw_block_names,
    components=30, include_pterm=True, include_gwb=True, log10_equad=-8.0,
    include_rn=True, rn_components=30,
):
    """Build discovery ArrayLikelihood.

    Returns (logl_fn, param_keys, base_values).
    logl_fn is JIT-compiled and takes a 1-D jnp array.
    """
    noisedict = {}
    for psr in disco_psrs:
        noisedict[psr.name + "_KAT_MKBF_efac"] = 1.0
        noisedict[psr.name + "_KAT_MKBF_log10_ecorr"] = -8.0
        noisedict[psr.name + "_KAT_MKBF_log10_t2equad"] = log10_equad

    noise_terms = {
        psr.name: ds.makenoise_measurement(psr, noisedict=noisedict) for psr in disco_psrs
    }
    timing_terms = {
        psr.name: ds.makegp_timing(psr, variance=1e-14) for psr in disco_psrs
    }

    T = ds.getspan(disco_psrs)
    if include_rn:
        common_gp = ds.makecommongp_fourier(
            disco_psrs, ds.powerlaw, rn_components, T, name="rednoise"
        )
    else:
        common_gp = None
    if include_gwb:
        global_gp = ds.makeglobalgp_fourier(
            disco_psrs, ds.powerlaw, ds.hd_orf, components, T, name="gwb"
        )
    else:
        global_gp = None

    temp_dict = _enterprise_to_disco_params(enterprise_params, cw_block_names)

    pulsar_likes = [
        ds.PulsarLikelihood([
            np.array(residual_map[psr.name], copy=True),
            noise_terms[psr.name],
            timing_terms[psr.name],
            MultiSourceDelay(psr, num_cw, include_pterm=include_pterm),
        ])
        for psr in disco_psrs
    ]

    fml = ds.ArrayLikelihood(pulsar_likes, commongp=common_gp, globalgp=global_gp)
    logl = fml.logL

    order_map = {k: i for i, k in enumerate(logl.params)}
    sorted_items = sorted(
        temp_dict.items(), key=lambda kv: order_map.get(kv[0], float("inf"))
    )
    sorted_dict = {k: v for k, v in sorted_items if k in logl.params}
    param_keys = list(sorted_dict.keys())
    base_vals = jnp.array([sorted_dict[k] for k in param_keys], dtype=jnp.float64)

    def logl_wrapped(x_array):
        return logl({k: v for k, v in zip(param_keys, x_array)})

    logl_fn = jax.jit(logl_wrapped)

    return logl_fn, param_keys, base_vals


# ============================================================
# Fast-scan likelihood: caches outer Cholesky for distance scans
# ============================================================
#
# During a distance scan, GP hyperparameters are frozen at truth.
# The outer Cholesky of (Phi_global^-1 + block_diag(T^T Sigma^-1 T))
# — the (N_psr * 2 * N_comp)^2 Cholesky that dominates cost at large
# N_psr — is therefore constant across scan points. Factor it once at
# construction, then each eval reuses the cached factor.
#
# Limitations:
#  - All GP hyperparameters must be frozen at the values in
#    `enterprise_params`; callers cannot vary rednoise/GWB amplitudes.
#  - The inner per-pulsar Cholesky inside `vsm.kernelterms` is still
#    recomputed each eval, but its cost is O(N_psr * n_comp^3), small
#    compared to the cached (N_psr * n_comp)^3 outer factor.
#
# At N_psr=116, n_comp=14 this gives ~50-100x per-eval speedup.

def build_fast_scan_likelihood(
    disco_psrs, residual_map, num_cw, enterprise_params, cw_block_names,
    components=30, include_pterm=True, log10_equad=-8.0,
    include_rn=True, rn_components=30,
):
    """Fast discovery logL with cached outer Cholesky. See module docstring.

    Requires globalgp (GWB GP). GP hyperparameters are frozen at the values
    in `enterprise_params`.

    Returns (logl_fn, param_keys, base_values) just like build_disco_likelihood.
    """
    from discovery import matrix as dsm

    noisedict = {}
    for psr in disco_psrs:
        noisedict[psr.name + "_KAT_MKBF_efac"] = 1.0
        noisedict[psr.name + "_KAT_MKBF_log10_ecorr"] = -8.0
        noisedict[psr.name + "_KAT_MKBF_log10_t2equad"] = log10_equad

    noise_terms = {
        psr.name: ds.makenoise_measurement(psr, noisedict=noisedict) for psr in disco_psrs
    }
    timing_terms = {
        psr.name: ds.makegp_timing(psr, variance=1e-14) for psr in disco_psrs
    }

    T = ds.getspan(disco_psrs)
    if include_rn:
        common_gp = ds.makecommongp_fourier(
            disco_psrs, ds.powerlaw, rn_components, T, name="rednoise"
        )
    else:
        common_gp = None
    global_gp = ds.makeglobalgp_fourier(
        disco_psrs, ds.powerlaw, ds.hd_orf, components, T, name="gwb"
    )

    temp_dict = _enterprise_to_disco_params(enterprise_params, cw_block_names)

    pulsar_likes = [
        ds.PulsarLikelihood([
            np.array(residual_map[psr.name], copy=True),
            noise_terms[psr.name],
            timing_terms[psr.name],
            MultiSourceDelay(psr, num_cw, include_pterm=include_pterm),
        ])
        for psr in disco_psrs
    ]

    # Build ArrayLikelihood and force logL construction so fml.vsm, fml.ys exist.
    fml = ds.ArrayLikelihood(pulsar_likes, commongp=common_gp, globalgp=global_gp)
    _ = fml.logL  # sets up vsm, ys

    # Replicate the setup from ArrayLikelihood.logL (likelihood.py:640-658)
    P_var_inv = getattr(global_gp, "Phi_inv", None) or global_gp.Phi.make_inv()
    Fs = getattr(global_gp, "Fs", None)
    if Fs is None:
        F = getattr(global_gp, "F", None)
        Fs = F if isinstance(F, (list, tuple)) else [F] * len(pulsar_likes)
    kterms = fml.vsm.make_kernelterms(fml.ys, Fs)

    npsr = len(Fs)
    ngp = Fs[0].shape[1]

    # Build a dict of truth parameter values for the required params
    truth_dict = {k: temp_dict[k] for k in kterms.params if k in temp_dict}
    for k in P_var_inv.params:
        if k in temp_dict:
            truth_dict[k] = temp_dict[k]

    # Evaluate once at truth params to obtain terms[2] = c (fixed when GP hyperparams are fixed)
    # and Pinv, ldP at truth (both fixed).
    truth_full = dict(truth_dict)
    # Add any placeholder keys that kterms may need (CW params, distances).
    for k in kterms.params:
        if k not in truth_full and k in temp_dict:
            truth_full[k] = temp_dict[k]

    terms_truth = kterms(truth_full)
    c_truth = terms_truth[2]   # shape (npsr, ngp, ngp), GP-hyperparam-dependent only

    Pinv_truth, ldP_truth = P_var_inv(truth_full)

    # Build the outer Cholesky once at construction time.
    FtNmF_outer = dsm.jsp.linalg.block_diag(*c_truth)
    cf_cached = dsm.matrix_factor(Pinv_truth + FtNmF_outer)

    # Precompute the outer log-det so it's a compile-time constant.
    ldP_sum_cached = jnp.sum(ldP_truth)
    logdet_cf_cached = dsm.matrix_norm * jnp.sum(jnp.log(jnp.diag(cf_cached[0])))

    # Build the sorted param_keys matching the original logl.params, then map
    # them to the original full parameter dict.
    full_param_set = sorted(set(kterms.params) | set(P_var_inv.params))
    order_map = {k: i for i, k in enumerate(full_param_set)}
    sorted_items = sorted(
        temp_dict.items(), key=lambda kv: order_map.get(kv[0], float("inf"))
    )
    sorted_dict = {k: v for k, v in sorted_items if k in full_param_set}
    param_keys = list(sorted_dict.keys())
    base_vals = jnp.array([sorted_dict[k] for k in param_keys], dtype=jnp.float64)

    def fast_loglike(x_array):
        params = {k: v for k, v in zip(param_keys, x_array)}
        terms = kterms(params)
        p0 = jnp.sum(terms[0])
        FtNmy = terms[1].reshape(npsr * ngp)
        # Use cached outer Cholesky — skips the (N_psr*ngp)^3 factorization.
        logp = p0 + 0.5 * (
            FtNmy.T @ dsm.matrix_solve(cf_cached, FtNmy)
            - ldP_sum_cached - logdet_cf_cached
        )
        return logp

    logl_fn = jax.jit(fast_loglike)

    return logl_fn, param_keys, base_vals


# ============================================================
# Distance scanning and peak analysis
# ============================================================

def compute_mode_spacing(cos_gwtheta, gwphi, log10_fgw, psr_pos):
    """Distance mode spacing dL (kpc) for one CW-pulsar pair."""
    gwtheta = np.arccos(cos_gwtheta)
    f_gw = 10.0 ** log10_fgw
    omhat = np.array([
        -np.sin(gwtheta) * np.cos(gwphi),
        -np.sin(gwtheta) * np.sin(gwphi),
        -np.cos(gwtheta),
    ])
    omhat_dot_pos = np.dot(omhat, psr_pos)
    denom = abs(1.0 + omhat_dot_pos)
    if denom < 1e-4:
        return np.inf
    return 1.0 / (f_gw * KPC_OVER_C * denom)


_BATCHED_SCAN_CACHE = {}


def _grid_with_required_points(scan_min, scan_max, n_points, required_points=None):
    """Build sorted scan grid and force specific points into the grid."""
    scan_vals = np.linspace(scan_min, scan_max, n_points)
    if required_points is None:
        return scan_vals

    required = np.asarray(required_points, dtype=float)
    required = required[(required >= scan_min) & (required <= scan_max)]
    if len(required) == 0:
        return scan_vals
    return np.unique(np.concatenate([scan_vals, required]))


def scan_pulsar_distance(logl_fn, base_values, param_keys, psr_dist_key,
                         scan_min, scan_max, n_points=3000, chunk_size=None,
                         n_components=30, required_points=None):
    """Scan one pulsar's distance holding everything else fixed.

    Returns (scan_values, logl_values).
    """
    # Adaptive chunk size to fit GPU memory.
    # The GP marginalisation requires Cholesky of an (N_psr * 2 * N_components)^2 matrix
    # per evaluation. With float64, memory per eval ≈ mat_size^2 * 8 bytes.
    # Target ~2GB max for the batched Cholesky intermediates.
    if chunk_size is None:
        # With the cached-Cholesky fast path, per-eval memory is dominated by
        # the (n_psr*2*n_comp)^2 solve rhs, not a fresh Cholesky — so we can
        # use much bigger batches than the slow path allowed.
        n_dist_params = sum(1 for k in param_keys if k.endswith("_cw_p_dist"))
        n_psr = max(n_dist_params, 2)
        mat_dim = n_psr * 2 * n_components
        mem_per_eval = mat_dim * 8 * 4  # a few vectors of size mat_dim
        target_mem = 4e9
        chunk_size = max(32, min(1024, int(target_mem / mem_per_eval)))

    scan_idx = param_keys.index(psr_dist_key)
    scan_vals = _grid_with_required_points(
        scan_min, scan_max, n_points, required_points=required_points
    )

    # Cache the batched/jit'd function per logl_fn id so repeated scans across
    # pulsars don't re-trigger XLA compilation (tens of seconds per compile
    # at N_psr=116). jitted fns don't accept attributes so we key on id.
    cache_key = (id(logl_fn), chunk_size)
    batched_fn = _BATCHED_SCAN_CACHE.get(cache_key)
    if batched_fn is None:
        batched_fn = jax.jit(jax.vmap(logl_fn))
        _BATCHED_SCAN_CACHE[cache_key] = batched_fn
    template = jnp.repeat(base_values[None, :], chunk_size, axis=0)

    pieces = []
    n_scan = len(scan_vals)
    n_chunks = int(np.ceil(n_scan / chunk_size))
    for idx in range(n_chunks):
        start = idx * chunk_size
        end = min(start + chunk_size, n_scan)
        current = jnp.array(scan_vals[start:end])
        pad = chunk_size - current.shape[0]
        padded = current if pad == 0 else jnp.pad(current, (0, pad), constant_values=current[-1])
        x_block = template.at[:, scan_idx].set(padded)
        block_vals = batched_fn(x_block)
        pieces.append(np.array(block_vals[: end - start]))

    return scan_vals, np.concatenate(pieces)


_BATCHED_2D_SCAN_CACHE = {}


def scan_pulsar_pair_2d(logl_fn, base_values, param_keys,
                        dist_key_i, dist_key_j, grid_i, grid_j,
                        chunk_size=None, n_components=30):
    """Evaluate lnL on a 2D pulsar-distance grid.

    Returns (grid_i, grid_j, lnL_surface) with
    lnL_surface[i, j] = lnL(grid_i[i], grid_j[j]).
    """
    grid_i = np.asarray(grid_i, dtype=float)
    grid_j = np.asarray(grid_j, dtype=float)
    ni, nj = len(grid_i), len(grid_j)
    n_total = ni * nj

    if chunk_size is None:
        n_dist_params = sum(1 for k in param_keys if k.endswith("_cw_p_dist"))
        n_psr = max(n_dist_params, 2)
        mat_dim = n_psr * 2 * n_components
        mem_per_eval = mat_dim * 8 * 4
        target_mem = 4e9
        chunk_size = max(32, min(1024, int(target_mem / mem_per_eval)))

    idx_i = param_keys.index(dist_key_i)
    idx_j = param_keys.index(dist_key_j)
    Di_mesh, Dj_mesh = np.meshgrid(grid_i, grid_j, indexing="ij")
    Di_flat = Di_mesh.ravel()
    Dj_flat = Dj_mesh.ravel()

    cache_key = (id(logl_fn), chunk_size)
    batched_fn = _BATCHED_2D_SCAN_CACHE.get(cache_key)
    if batched_fn is None:
        batched_fn = jax.jit(jax.vmap(logl_fn))
        _BATCHED_2D_SCAN_CACHE[cache_key] = batched_fn

    template = jnp.repeat(base_values[None, :], chunk_size, axis=0)
    pieces = []
    n_chunks = int(np.ceil(n_total / chunk_size))
    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, n_total)
        di = jnp.array(Di_flat[start:end])
        dj = jnp.array(Dj_flat[start:end])
        pad = chunk_size - (end - start)
        if pad:
            di = jnp.pad(di, (0, pad), constant_values=di[-1])
            dj = jnp.pad(dj, (0, pad), constant_values=dj[-1])
        x_block = template.at[:, idx_i].set(di).at[:, idx_j].set(dj)
        vals = batched_fn(x_block)
        pieces.append(np.array(vals[: end - start]))

    return grid_i, grid_j, np.concatenate(pieces).reshape(ni, nj)


def scan_pair_window_2d(logl_fn, base_values, param_keys,
                        psr_i, psr_j, cw_params_list,
                        center_i=None, center_j=None,
                        half_width_modes=3.0, points_per_mode=8,
                        min_points=41, chunk_size=None,
                        n_components=30, include_truth=True):
    """Resolve a local 2D distance window around current pair estimates.

    This is intended for coordinate-ascent updates. It does not do a
    low-resolution full-prior pre-scan, because that aliases badly when the
    prior contains many distance modes.
    """
    center_i = psr_i.pdist[0] if center_i is None else center_i
    center_j = psr_j.pdist[0] if center_j is None else center_j

    dL_i = min(
        compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"], cw["log10_fgw"], psr_i.pos)
        for cw in cw_params_list
    )
    dL_j = min(
        compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"], cw["log10_fgw"], psr_j.pos)
        for cw in cw_params_list
    )

    n_i = int(max(min_points, np.ceil(2 * half_width_modes * points_per_mode + 1)))
    n_j = int(max(min_points, np.ceil(2 * half_width_modes * points_per_mode + 1)))

    req_i = [psr_i.pdist[0]] if include_truth else None
    req_j = [psr_j.pdist[0]] if include_truth else None
    grid_i = _grid_with_required_points(
        max(0.01, center_i - half_width_modes * dL_i),
        center_i + half_width_modes * dL_i,
        n_i,
        required_points=req_i,
    )
    grid_j = _grid_with_required_points(
        max(0.01, center_j - half_width_modes * dL_j),
        center_j + half_width_modes * dL_j,
        n_j,
        required_points=req_j,
    )

    dist_key_i = f"{psr_i.name}_cw_p_dist"
    dist_key_j = f"{psr_j.name}_cw_p_dist"
    grid_i, grid_j, lnL = scan_pulsar_pair_2d(
        logl_fn, base_values, param_keys,
        dist_key_i, dist_key_j, grid_i, grid_j,
        chunk_size=chunk_size, n_components=n_components,
    )
    mi, mj = np.unravel_index(np.argmax(lnL), lnL.shape)

    return {
        "grid_i": grid_i,
        "grid_j": grid_j,
        "lnL": lnL,
        "peak_i": float(grid_i[mi]),
        "peak_j": float(grid_j[mj]),
        "peak_lnL": float(lnL[mi, mj]),
        "dL_i": float(dL_i),
        "dL_j": float(dL_j),
    }


def analyze_peaks(scan_vals, logl_vals, true_dist, mode_spacing=None):
    """Find peaks in a distance scan and compare to the truth peak.

    Returns a dict with:
      - n_peaks, truth_logl, best_wrong_logl, delta_lnL_truth_vs_best_wrong
      - peak_distances, peak_heights  (arrays)
      - peak_delta_lnL: per-peak Δln L relative to truth peak
    """
    # Adaptive prominence: use a fraction of the total range
    prom = max(0.5, 0.01 * (np.max(logl_vals) - np.min(logl_vals)))
    peaks, props = find_peaks(logl_vals, prominence=prom, distance=3)

    if len(peaks) == 0:
        truth_ll = float(np.interp(true_dist, scan_vals, logl_vals))
        return {
            "n_peaks": 0, "truth_logl": truth_ll,
            "best_wrong_logl": np.nan,
            "delta_lnL_truth_vs_best_wrong": np.nan,
            "peak_distances": np.array([]),
            "peak_heights": np.array([]),
            "peak_delta_lnL": np.array([]),
        }

    peak_dists = scan_vals[peaks]
    peak_heights = logl_vals[peaks]

    # Identify truth peak
    truth_idx = np.argmin(np.abs(peak_dists - true_dist))
    truth_ll = peak_heights[truth_idx]

    # Non-truth peaks: exclude the peak closest to truth
    # (use mode_spacing if available to define "close")
    exclude_radius = 0.02 if mode_spacing is None else 0.3 * mode_spacing
    non_truth = np.abs(peak_dists - true_dist) > exclude_radius

    if np.any(non_truth):
        best_wrong_ll = np.max(peak_heights[non_truth])
    else:
        best_wrong_ll = np.nan

    delta = truth_ll - peak_heights
    delta[truth_idx] = 0.0  # truth peak itself

    return {
        "n_peaks": len(peaks),
        "truth_logl": truth_ll,
        "best_wrong_logl": best_wrong_ll,
        "delta_lnL_truth_vs_best_wrong": truth_ll - best_wrong_ll if not np.isnan(best_wrong_ll) else np.nan,
        "peak_distances": peak_dists,
        "peak_heights": peak_heights,
        "peak_delta_lnL": delta,
    }


# ============================================================
# Loading helpers
# ============================================================

def load_pulsars(n_pulsars=10):
    """Load enterprise and discovery pulsars (first n_pulsars)."""
    ent_psrs = load_feathers.load_feathers_from_folder(FEATHER_DIR)[:n_pulsars]
    disco_psrs = [
        ds.Pulsar.read_feather(f)
        for f in sorted(glob.glob(FEATHER_DIR + "*.feather"))
    ][:n_pulsars]
    for psr in disco_psrs:
        psr.toaerrs = np.full_like(psr.toas, 1e-6, dtype=np.float64)
    return ent_psrs, disco_psrs


# ============================================================
# Best wrong mode within EM distance prior
# ============================================================

def find_best_wrong_mode_in_prior(
    logl_fn, base_values, param_keys, psr_name,
    true_dist, dist_sigma, mode_spacing,
    n_sigma=3.0, n_points=None, n_components=30,
):
    """Scan one pulsar's distance within its EM prior and find the best non-truth peak.

    Parameters
    ----------
    true_dist : float
        EM distance mean (kpc).
    dist_sigma : float
        EM distance uncertainty (kpc).
    mode_spacing : float
        Smallest mode spacing across all CWs for this pulsar (kpc).
    n_sigma : float
        Number of sigma to scan either side of true_dist.
    n_points : int or None
        If None, adaptively choose ~4 points per mode within the prior.

    Returns
    -------
    dict with:
      best_wrong_dist  : distance of the best non-truth peak (kpc)
      best_wrong_logl  : log-likelihood at that peak
      truth_logl       : log-likelihood at the truth-nearest peak
      delta_lnL        : truth_logl - best_wrong_logl
      n_peaks_in_prior : number of peaks found within the prior range
      scan_d, scan_ll  : the full scan arrays (for plotting)
    """
    dist_key = f"{psr_name}_cw_p_dist"
    scan_min = max(1e-3, true_dist - n_sigma * dist_sigma)
    scan_max = true_dist + n_sigma * dist_sigma

    if n_points is None:
        prior_width = scan_max - scan_min
        n_modes = prior_width / mode_spacing if mode_spacing > 0 else 100
        n_points = int(np.clip(4 * n_modes, 200, 2000))

    scan_d, scan_ll = scan_pulsar_distance(
        logl_fn, base_values, param_keys, dist_key,
        scan_min, scan_max, n_points, n_components=n_components,
    )

    analysis = analyze_peaks(scan_d, scan_ll, true_dist, mode_spacing=mode_spacing)

    return {
        "best_wrong_dist": analysis["peak_distances"][
            np.argmax(analysis["peak_heights"][
                np.abs(analysis["peak_distances"] - true_dist) > 0.3 * mode_spacing
            ])
        ] if analysis["n_peaks"] > 1 and np.any(
            np.abs(analysis["peak_distances"] - true_dist) > 0.3 * mode_spacing
        ) else np.nan,
        "best_wrong_logl": analysis["best_wrong_logl"],
        "truth_logl": analysis["truth_logl"],
        "delta_lnL": analysis["delta_lnL_truth_vs_best_wrong"],
        "n_peaks_in_prior": analysis["n_peaks"],
        "scan_d": scan_d,
        "scan_ll": scan_ll,
    }


def compute_joint_best_wrong_in_prior(
    logl_fn, base_values, param_keys,
    ent_psrs, disco_psrs, cw_params_or_inj,
    cw_block_names=None, n_sigma=3.0, n_points=None, n_components=30,
    keep_scan_arrays=False,
):
    """For each pulsar, find its best wrong mode within the EM prior, then
    evaluate the joint likelihood with ALL pulsars at their respective best
    wrong modes.

    Parameters
    ----------
    cw_params_or_inj : list of dicts OR dict
        If list of dicts: noise-free CW param dicts (each has cos_gwtheta, gwphi, log10_fgw).
        If dict: enterprise injection params with cw_block_names prefixed keys.
    cw_block_names : list or None
        Required if cw_params_or_inj is an enterprise injection dict.

    Returns
    -------
    dict with:
      ll_truth       : joint log-likelihood at truth distances
      ll_best_wrong  : joint log-likelihood at per-pulsar best wrong modes
      delta_lnL      : ll_truth - ll_best_wrong
      per_pulsar     : dict mapping psr name -> per-pulsar scan results
    """
    ll_truth = float(logl_fn(base_values))

    # Extract CW source params for mode spacing calculation
    if isinstance(cw_params_or_inj, list):
        # Noise-free: list of CW param dicts
        cw_source_params = cw_params_or_inj
    else:
        # Enterprise injection dict: extract from cw_block_names
        cw_source_params = []
        for name in cw_block_names:
            cw_source_params.append({
                "cos_gwtheta": cw_params_or_inj[f"{name}_cos_gwtheta"],
                "gwphi": cw_params_or_inj[f"{name}_gwphi"],
                "log10_fgw": cw_params_or_inj[f"{name}_log10_fgw"],
            })

    wrong_vals = np.array(base_values)
    per_pulsar = {}

    for pe, pd in zip(ent_psrs, disco_psrs):
        dist_key = f"{pe.name}_cw_p_dist"
        if dist_key not in param_keys:
            continue

        true_dist = pe.pdist[0]
        dist_sigma = pe.pdist[1]

        # Smallest mode spacing across all CWs for this pulsar
        all_dL = [
            compute_mode_spacing(cp["cos_gwtheta"], cp["gwphi"], cp["log10_fgw"], pd.pos)
            for cp in cw_source_params
        ]
        min_dL = min(all_dL)

        result = find_best_wrong_mode_in_prior(
            logl_fn, base_values, param_keys, pe.name,
            true_dist, dist_sigma, min_dL,
            n_sigma=n_sigma, n_points=n_points, n_components=n_components,
        )

        # Distilled summary (keeps interpretability without megabytes of scans)
        best_wrong_offset = (result["best_wrong_dist"] - true_dist
                             if not np.isnan(result["best_wrong_dist"]) else np.nan)
        summary = {
            "truth_logl": result["truth_logl"],
            "best_wrong_logl": result["best_wrong_logl"],
            "delta_lnL": result["delta_lnL"],
            "best_wrong_dist": result["best_wrong_dist"],
            "best_wrong_offset_kpc": best_wrong_offset,
            "best_wrong_offset_modes": (best_wrong_offset / min_dL
                                        if not np.isnan(best_wrong_offset) else np.nan),
            "n_peaks_in_prior": result["n_peaks_in_prior"],
            "mode_spacing": min_dL,
            "true_dist": true_dist,
            "dist_sigma": dist_sigma,
            "prior_width_kpc": 2 * n_sigma * dist_sigma,
            "prior_width_modes": 2 * n_sigma * dist_sigma / min_dL if min_dL > 0 else np.nan,
        }
        if keep_scan_arrays:
            summary["scan_d"] = result["scan_d"]
            summary["scan_ll"] = result["scan_ll"]
        per_pulsar[pe.name] = summary

        # Set this pulsar to its best wrong mode in the joint vector
        if not np.isnan(result["best_wrong_dist"]):
            idx = param_keys.index(dist_key)
            wrong_vals[idx] = result["best_wrong_dist"]

    ll_best_wrong = float(logl_fn(jnp.array(wrong_vals)))

    return {
        "ll_truth": ll_truth,
        "ll_best_wrong": ll_best_wrong,
        "delta_lnL": ll_truth - ll_best_wrong,
        "per_pulsar": per_pulsar,
    }


# ============================================================
# Noise-free CW residual injection (bypasses enterprise)
# ============================================================

def inject_noisefree_cw(disco_psrs, cw_params_list):
    """Inject CW signals directly using discovery's phase_connected_binary.

    Parameters
    ----------
    disco_psrs : list of discovery Pulsar objects
    cw_params_list : list of dicts, each with keys:
        cos_gwtheta, gwphi, cos_inc, log10_mc, log10_fgw, log10_h, phase0, psi

    Returns
    -------
    residual_map : dict mapping pulsar name -> residual array (pure CW, no noise)
    true_dists : dict mapping pulsar name -> true distance used
    """
    cw_func = make_phase_connected_binary(pulsarterm=True)
    vmapped_cw = jax.vmap(
        cw_func,
        in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0, None),
    )

    # Stack CW params into arrays for vectorised evaluation
    keys = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
            "log10_fgw", "log10_h", "phase0", "psi")
    stacked = {k: jnp.array([cw[k] for cw in cw_params_list]) for k in keys}

    residual_map = {}
    true_dists = {}
    for psr in disco_psrs:
        # (n_sources, n_toas) -> sum over sources -> (n_toas,)
        all_delays = vmapped_cw(
            psr.toas, psr.pos,
            stacked["cos_gwtheta"], stacked["gwphi"],
            stacked["cos_inc"], stacked["log10_mc"],
            stacked["log10_fgw"], stacked["log10_h"],
            stacked["phase0"], stacked["psi"],
            psr.pdist[0],
        )
        residual_map[psr.name] = np.array(jnp.sum(all_delays, axis=0))
        true_dists[psr.name] = psr.pdist[0]

    return residual_map, true_dists


def build_noisefree_likelihood(disco_psrs, residual_map, num_cw, cw_params_list,
                               sigma_toa=1e-6, include_gwb=False,
                               include_red_noise=True, log10_equad=-8.0):
    """Build a simple white-noise-only likelihood for the noise-free case.

    Uses the same discovery ArrayLikelihood infrastructure but with:
    - Fixed white noise (efac=1, equad=-8)
    - A negligible-amplitude common red GP by default to keep discovery on its
      JAX-compatible ArrayLikelihood path
    - No GWB GP by default
    - CW delay model with pulsar term

    Returns (logl_fn, param_keys, base_values).
    """
    noisedict = {}
    for psr in disco_psrs:
        noisedict[psr.name + "_KAT_MKBF_efac"] = 1.0
        noisedict[psr.name + "_KAT_MKBF_log10_ecorr"] = -8.0
        noisedict[psr.name + "_KAT_MKBF_log10_t2equad"] = log10_equad

    noise_terms = {
        psr.name: ds.makenoise_measurement(psr, noisedict=noisedict) for psr in disco_psrs
    }
    timing_terms = {
        psr.name: ds.makegp_timing(psr, variance=1e-14) for psr in disco_psrs
    }

    T = ds.getspan(disco_psrs)
    if include_red_noise:
        common_gp = ds.makecommongp_fourier(disco_psrs, ds.powerlaw, 30, T, name="rednoise")
    else:
        common_gp = None
    if include_gwb:
        global_gp = ds.makeglobalgp_fourier(
            disco_psrs, ds.powerlaw, ds.hd_orf, 30, T, name="gwb"
        )
    else:
        global_gp = None

    pulsar_likes = [
        ds.PulsarLikelihood([
            np.array(residual_map[psr.name], copy=True),
            noise_terms[psr.name],
            timing_terms[psr.name],
            MultiSourceDelay(psr, num_cw, include_pterm=True),
        ])
        for psr in disco_psrs
    ]

    fml = ds.ArrayLikelihood(pulsar_likes, commongp=common_gp, globalgp=global_gp)
    logl = fml.logL

    # Build parameter vector: CW params + distances + noise GP params at benign values
    disco_suffixes = ["" if i == 0 else f"_{i+1}" for i in range(num_cw)]
    param_dict = {}

    for idx, (suffix, cw_p) in enumerate(zip(disco_suffixes, cw_params_list)):
        for key in MultiSourceDelay.global_params:
            param_dict[f"cw_{key}{suffix}"] = cw_p[key]

    for psr in disco_psrs:
        param_dict[f"{psr.name}_cw_p_dist"] = psr.pdist[0]

    # Red noise (per-pulsar in discovery) at negligible amplitudes
    if include_red_noise:
        for psr in disco_psrs:
            param_dict[f"{psr.name}_rednoise_log10_A"] = -20.0
            param_dict[f"{psr.name}_rednoise_gamma"] = 4.0
    # GWB at negligible amplitude (only if included in likelihood)
    if include_gwb:
        param_dict["gwb_log10_A"] = -20.0
        param_dict["gwb_gamma"] = 4.333

    order_map = {k: i for i, k in enumerate(logl.params)}
    sorted_items = sorted(
        param_dict.items(), key=lambda kv: order_map.get(kv[0], float("inf"))
    )
    sorted_dict = {k: v for k, v in sorted_items if k in logl.params}
    param_keys = list(sorted_dict.keys())
    base_vals = jnp.array([sorted_dict[k] for k in param_keys], dtype=jnp.float64)

    def logl_wrapped(x_array):
        return logl({k: v for k, v in zip(param_keys, x_array)})

    logl_fn = jax.jit(logl_wrapped)

    return logl_fn, param_keys, base_vals


def build_pure_cw_likelihood(disco_psrs, residual_map, cw_params_list,
                             log10_equad=-8.0):
    """Build JIT-safe white-noise + timing-model + CW likelihood, no GP terms.

    This avoids ``ArrayLikelihood(commongp=None)`` because that route can call
    SciPy solves inside a JIT trace.  Instead each pulsar's kernel-product
    function is built once with the callable CW residual model; discovery then
    uses its JAX-safe Woodbury path during evaluation.
    """
    num_cw = len(cw_params_list)

    noisedict = {}
    for psr in disco_psrs:
        noisedict[psr.name + "_KAT_MKBF_efac"] = 1.0
        noisedict[psr.name + "_KAT_MKBF_log10_ecorr"] = -8.0
        noisedict[psr.name + "_KAT_MKBF_log10_t2equad"] = log10_equad

    pulsar_likes = [
        ds.PulsarLikelihood([
            np.array(residual_map[psr.name], copy=True),
            ds.makenoise_measurement(psr, noisedict=noisedict),
            ds.makegp_timing(psr, variance=1e-14),
            MultiSourceDelay(psr, num_cw, include_pterm=True),
        ])
        for psr in disco_psrs
    ]

    kp_fns = [psl.N.make_kernelproduct(psl.y) for psl in pulsar_likes]
    all_params = sorted(set.union(*[set(kp.params) for kp in kp_fns]))

    disco_suffixes = ["" if i == 0 else f"_{i+1}" for i in range(num_cw)]
    param_dict = {}
    for suffix, cw_p in zip(disco_suffixes, cw_params_list):
        for key in MultiSourceDelay.global_params:
            param_dict[f"cw_{key}{suffix}"] = cw_p[key]
    for psr in disco_psrs:
        param_dict[f"{psr.name}_cw_p_dist"] = psr.pdist[0]

    order_map = {k: i for i, k in enumerate(all_params)}
    sorted_items = sorted(
        param_dict.items(), key=lambda kv: order_map.get(kv[0], float("inf"))
    )
    sorted_dict = {k: v for k, v in sorted_items if k in all_params}
    param_keys = list(sorted_dict.keys())
    base_vals = jnp.array([sorted_dict[k] for k in param_keys], dtype=jnp.float64)

    def logl_wrapped(x_array):
        params = {k: v for k, v in zip(param_keys, x_array)}
        return sum(kp(params) for kp in kp_fns)

    return jax.jit(logl_wrapped), param_keys, base_vals


# ============================================================
# Gradient distance optimization helpers
# ============================================================

def make_distance_optimizer(logl_fn, param_keys, base_values, disco_psrs,
                            prior_means=None, prior_sigmas=None,
                            n_sigma=3.0, objective="logpost", maxiter=200,
                            factr=None):
    """Build gradient-based optimizer for pulsar distances only.

    Optimizes scaled coordinates z_i = (d_i - prior_mean_i) / prior_sigma_i.
    Returns (optimize, fun_np, jac_np, neg_objective_z), where optimize takes
    and returns kpc.
    """
    import scipy.optimize

    if objective not in ("lnL", "logpost"):
        raise ValueError("objective must be 'lnL' or 'logpost'")

    n_psr = len(disco_psrs)
    if prior_means is None:
        prior_means = [psr.pdist[0] for psr in disco_psrs]
    if prior_sigmas is None:
        prior_sigmas = [psr.pdist[1] for psr in disco_psrs]

    dist_keys = [f"{psr.name}_cw_p_dist" for psr in disco_psrs]
    missing = [k for k in dist_keys if k not in param_keys]
    if missing:
        raise KeyError(f"Missing distance parameters in param_keys: {missing}")

    dist_indices_np = np.array([param_keys.index(k) for k in dist_keys], dtype=int)
    dist_indices = jnp.array(dist_indices_np)
    base_values = jnp.asarray(base_values, dtype=jnp.float64)
    mu = jnp.asarray(prior_means, dtype=jnp.float64)
    sig = jnp.asarray(prior_sigmas, dtype=jnp.float64)
    if mu.shape != (n_psr,) or sig.shape != (n_psr,):
        raise ValueError("prior_means and prior_sigmas must match disco_psrs length")
    if np.any(np.asarray(sig) <= 0):
        raise ValueError("prior_sigmas must be positive")

    z_lo = jnp.full(n_psr, -float(n_sigma), dtype=jnp.float64)
    z_hi = jnp.full(n_psr, float(n_sigma), dtype=jnp.float64)
    z_lo = jnp.maximum(z_lo, (0.01 - mu) / sig)
    z_lo_np = np.asarray(z_lo, dtype=np.float64)
    z_hi_np = np.asarray(z_hi, dtype=np.float64)
    bad_bounds = np.where(z_lo_np >= z_hi_np)[0]
    if len(bad_bounds):
        raise ValueError(f"Infeasible distance bounds for pulsar indices: {bad_bounds.tolist()}")

    def z_to_d(z):
        return mu + sig * z

    def d_to_z(d):
        return (d - mu) / sig

    @jax.jit
    def neg_objective_z(z):
        d = z_to_d(z)
        x = base_values.at[dist_indices].set(d)
        lnL = logl_fn(x)
        if objective == "logpost":
            return -(lnL - 0.5 * jnp.sum(z ** 2))
        return -lnL

    neg_grad_z = jax.jit(jax.grad(neg_objective_z))

    def fun_np(z_np):
        return float(neg_objective_z(jnp.asarray(z_np, dtype=jnp.float64)))

    def jac_np(z_np):
        return np.asarray(
            neg_grad_z(jnp.asarray(z_np, dtype=jnp.float64)), dtype=np.float64
        )

    z_test = d_to_z(mu)
    _ = fun_np(np.asarray(z_test, dtype=np.float64))
    _ = jac_np(np.asarray(z_test, dtype=np.float64))

    def optimize(d_start_kpc):
        d_start_kpc = np.asarray(d_start_kpc, dtype=np.float64)
        z0 = np.asarray(d_to_z(jnp.asarray(d_start_kpc)), dtype=np.float64)
        z0 = np.clip(z0, z_lo_np + 1e-10, z_hi_np - 1e-10)
        ftol = 1e-14 if factr is None else float(factr) * np.finfo(float).eps
        try:
            res = scipy.optimize.minimize(
                fun_np, z0, jac=jac_np, method="L-BFGS-B",
                bounds=list(zip(z_lo_np, z_hi_np)),
                options={"maxiter": maxiter, "ftol": ftol, "gtol": 1e-10},
            )
            d_opt = np.asarray(z_to_d(jnp.asarray(res.x)), dtype=np.float64)
            x_opt = np.array(base_values, dtype=np.float64, copy=True)
            x_opt[dist_indices_np] = d_opt
            lnL_opt = float(logl_fn(jnp.asarray(x_opt, dtype=jnp.float64)))
            z_opt = (d_opt - np.asarray(prior_means, dtype=np.float64)) / np.asarray(prior_sigmas, dtype=np.float64)
            logprior = -0.5 * np.sum(z_opt ** 2)
            jac = np.asarray(getattr(res, "jac", np.full_like(res.x, np.nan)), dtype=np.float64)
            return d_opt, lnL_opt, lnL_opt + logprior, dict(
                success=bool(res.success),
                message=str(res.message),
                nfev=getattr(res, "nfev", None),
                njev=getattr(res, "njev", None),
                nit=getattr(res, "nit", None),
                grad_norm=float(np.linalg.norm(jac)),
                z=np.asarray(res.x, dtype=np.float64),
            )
        except Exception as e:
            return d_start_kpc.copy(), np.nan, np.nan, dict(
                success=False, message=str(e), nfev=0, njev=0, nit=0,
                grad_norm=np.nan, z=z0,
            )

    optimize.neg_objective_z = neg_objective_z
    optimize.neg_grad_z = neg_grad_z
    optimize.z_to_d = z_to_d
    optimize.d_to_z = d_to_z
    optimize.z_bounds = (z_lo_np.copy(), z_hi_np.copy())
    optimize.dist_indices = dist_indices_np.copy()

    return optimize, fun_np, jac_np, neg_objective_z


def enumerate_pulsar_modes(psr, cw_params_list, prior_mean, prior_sigma,
                           n_sigma=3.0, empirical_reference=None,
                           snr_weights=None, score_threshold=None,
                           max_candidates=30):
    """Enumerate candidate distance modes for one pulsar.

    Candidates are pooled from each CW mode family and ranked by multi-CW
    consistency.  Higher score means closer agreement with all mode families.
    """
    empirical_reference = empirical_reference or {}
    n_cw = len(cw_params_list)
    weights = np.ones(n_cw, dtype=float) if snr_weights is None else np.asarray(snr_weights, dtype=float)
    if weights.shape != (n_cw,):
        raise ValueError("snr_weights must match cw_params_list length")
    if np.sum(weights) <= 0:
        weights = np.ones(n_cw, dtype=float)
    weights = weights / np.sum(weights)

    lo = max(0.01, float(prior_mean) - float(n_sigma) * float(prior_sigma))
    hi = float(prior_mean) + float(n_sigma) * float(prior_sigma)

    per_cw_spacings = []
    families = []
    pooled = [float(prior_mean)]
    for k, cw in enumerate(cw_params_list):
        dL = float(compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"], cw["log10_fgw"], psr.pos))
        per_cw_spacings.append(dL)
        if not np.isfinite(dL) or dL <= 0:
            families.append(np.array([float(prior_mean)]))
            continue

        d_ref = float(empirical_reference.get(k, prior_mean))
        n_min = int(np.floor((lo - d_ref) / dL)) - 1
        n_max = int(np.ceil((hi - d_ref) / dL)) + 1
        centers = d_ref + np.arange(n_min, n_max + 1, dtype=float) * dL
        centers = centers[(centers >= lo) & (centers <= hi)]
        if len(centers) == 0:
            centers = np.array([np.clip(d_ref, lo, hi)], dtype=float)
        families.append(centers)
        pooled.extend(centers.tolist())

    raw_candidates = np.unique(np.clip(np.asarray(pooled, dtype=float), lo, hi))
    scores = []
    for d in raw_candidates:
        score = 0.0
        for k, (centers, dL) in enumerate(zip(families, per_cw_spacings)):
            if not np.isfinite(dL) or dL <= 0:
                continue
            delta = np.min(np.abs(d - centers))
            score -= weights[k] * (delta / dL) ** 2
        scores.append(score)
    scores = np.asarray(scores, dtype=float)

    if score_threshold is not None:
        keep = scores >= float(score_threshold)
        keep |= np.isclose(raw_candidates, float(prior_mean), rtol=0, atol=1e-12)
        raw_candidates = raw_candidates[keep]
        scores = scores[keep]

    order = np.argsort(scores)[::-1]
    raw_candidates = raw_candidates[order][:max_candidates]
    scores = scores[order][:max_candidates]
    return raw_candidates, scores, per_cw_spacings
