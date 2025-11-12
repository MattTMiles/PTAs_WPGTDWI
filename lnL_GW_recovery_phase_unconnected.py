import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
os.environ["JAX_DISABLE_MMAP_CACHE"] = "1"
os.environ["XLA_FLAGS"] = "--xla_gpu_autotune_level=2"

import sys
import numpy as np
import discovery as ds
import numpy as np
import glob
import matplotlib.pyplot as plt
import jax
jax.config.update('jax_enable_x64', True)
jax.clear_caches()
import jax.numpy as jnp
import functools
from tqdm import tqdm

from enterprise.signals import signal_base

import numpy as np
import scipy.linalg as sl
import scipy.sparse as ss
from sksparse.cholmod import cholesky
from enterprise.signals import signal_base
from enterprise.signals.gp_signals import get_timing_model_basis, BasisGP
from enterprise.signals.parameter import function
from enterprise_extensions import load_feathers


import argparse
import time

    
from enterprise.signals import selections
import enterprise.signals.parameter as parameter
from enterprise.signals import white_signals
from enterprise_extensions.blocks import common_red_noise_block
from enterprise.signals import gp_signals    
from enterprise.signals import utils
from enterprise_extensions import deterministic
from discovery.deterministic import make_phase_connected_binary, make_phase_unconnected_binary
from discovery import const
from pathlib import Path


# Set data product directory. The simulation code requires pulsar timing models and data to initialise. 
# In this directory, we have stored these files for a set of pulsars (mimicking the IPTA DR2 dataset, but nonetheless synthetic).

feather_dir = "./data_products/"

CURRENT_LOGL = None
PARAM_KEYS = ()
gwamp_index = None
chunk_template = None

parser = argparse.ArgumentParser(description="GW amplitude recovery scan.")
parser.add_argument(
    "start_realisation",
    type=int,
    nargs="?",
    default=0,
    help="Index of the first realisation to run (default: 0).",
)
parser.add_argument(
    "--n-real",
    type=int,
    default=1,
    help="Number of consecutive realisations to run in a single process.",
)
parser.add_argument(
    "--save-prefix",
    default="logls_batch",
    help="Prefix for the aggregated output file (.npz).",
)
parser.add_argument(
    "--save-individual",
    action="store_true",
    help="Also write each realisation as its own .npy file.",
)
parser.add_argument(
    "--scan-param",
    #choices=["gwb_log10_A", "cw_log10_h", "cw_log10_h_2", "J0030+0451_cw_p_dist"],
    default="gwb_log10_A",
    help="Which amplitude parameter to scan over.",
)
parser.add_argument(
    "--pulsar-term",
    action="store_true",
    help="Include pulsar term when building the deterministic CW delay.",
)

parser.add_argument(
    "--seed",
    type=int,
    default=None,
    help="Seed for NumPy RNG",
)
parser.add_argument(
    "--output-root",
    default="lnLs_GWAmp",
    help="Directory where per-run folders and artefacts are written.",
)
parser.add_argument(
    "--n-cw-sources", 
    type=int, 
    default=1, 
    help="How many phase-connected CW signals to inject/model.",
    )
parser.add_argument(
    "--n-pulsars", 
    type=int, 
    default=10, 
    help="How many pulsars to include in the model.",
    )

args = parser.parse_args()

Npulsars = args.n_pulsars
num_cw = args.n_cw_sources

cw_block_names = ["cw" if i == 0 else f"cw{i+1}" for i in range(num_cw)]
disco_suffixes = ["" if i == 0 else f"_{i+1}" for i in range(num_cw)]

base_seed = args.seed if args.seed is not None else np.random.SeedSequence().entropy
rng = np.random.default_rng(base_seed)
if args.seed is not None:
    np.random.seed(args.seed)

scan_param = args.scan_param
_SCAN_RANGES = {
    "gwb_log10_A": (-17.0, -11.0),
    "cw_log10_h": (-17.0, -12.5),
    "cw_log10_h_2": (-17.0, -12.5),
}
_SCAN_RANGES.update({f"cw_log10_h{suffix}": (-17.0, -12.5) for suffix in disco_suffixes})


def simulate(pta, params, sparse_cholesky=True):
    """Simulate code with enterprise (instead of libstempo/PINT)"""
    delays, ndiags, fmats, phis = (pta.get_delay(params=params),
                                   pta.get_ndiag(params=params),
                                   pta.get_basis(params=params),
                                   pta.get_phi(params=params))

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

        assert len(gp) == i
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
            # this code is very slow...
            n = np.diag(ndiag._nvec)
            for j,s in zip(ndiag._jvec, ndiag._slices):
                n[s,s] += j
            whiteresiduals.append(delay + np.dot(sl.cholesky(n, lower=True), np.random.randn(n.shape[0])))
        elif ndiag.ndim == 1:
            whiteresiduals.append(delay + np.sqrt(ndiag) * np.random.randn(ndiag.shape[0]))
        else:
            raise NotImplementedError

    return [np.array(g + w) for g, w in zip(gpresiduals, whiteresiduals)]

def set_residuals(psr, y):
    psr._residuals[psr._isort] = y

@function
def tm_prior(weights, toas, variance=1e40):
    return weights * variance * len(toas)

def TimingModel(coefficients=False, name="linear_timing_model",
                use_svd=False, normed=True, prior_variance=1e40):
    """Class factory for marginalized linear timing model signals."""

    basis = get_timing_model_basis(use_svd, normed)
    prior = tm_prior(variance=prior_variance)

    BaseClass = BasisGP(prior, basis, coefficients=coefficients, name=name)

    class TimingModel(BaseClass):
        signal_type = "basis"
        signal_name = "linear timing model"
        signal_id = name + "_svd" if use_svd else name

    return TimingModel

# Load the pulsars data products as Enterprise Pulsar objects
psrs = load_feathers.load_feathers_from_folder(feather_dir)

_SCAN_RANGES.update({f"{psr.name}_cw_p_dist": (0.1, 5.0) for psr in psrs})
scan_min, scan_max = _SCAN_RANGES[scan_param]
gridsteps = 5000
scan_values = jnp.linspace(scan_min, scan_max, gridsteps)


# The _pdist attribute is used by enterprise to store the pulsar distance. Here, we back this up to allow us to vary it later.
for psr in psrs:
    psr._pdist = psr.pdist
    psr.residuals = np.array(psr.toas) * 0.0

tmin = [p.toas.min() for p in psrs]
tmax = [p.toas.max() for p in psrs]
Tspan = np.max(tmax) - np.min(tmin)

selection = selections.Selection(selections.by_backend)
efac = parameter.Constant(1)
equad = parameter.Constant(-8)
# ecorr = parameter.Constant(-8)

log10_A_red = parameter.Uniform(-18, -11)
gamma_red = parameter.Uniform(0, 7)
log10_A_gw = parameter.Uniform(-18,-11)('gwb_log10_A')
gamma_gw = parameter.Uniform(0,7)('gwb_gamma')

components = 30

tm = TimingModel(coefficients=False, name="linear_timing_model",
                use_svd=False, normed=True, prior_variance=1e-14)
models = []
for p in psrs:
    s = tm

    ef = white_signals.MeasurementNoise(efac=efac, selection=selection)
    s += ef

    eq = white_signals.TNEquadNoise(log10_tnequad=equad, selection=selection)
    s += eq

    # ec = gp_signals.EcorrBasisModel(log10_ecorr=ecorr,selection=selection)
    # s += ec

    pl = utils.powerlaw(log10_A=log10_A_red, gamma=gamma_red)
    rn = gp_signals.FourierBasisGP(spectrum=pl, components=components, Tspan=Tspan, name="rednoise")
    s += rn

    ## HD Correlated signal
    crn = common_red_noise_block(psd='powerlaw', prior='log-uniform',
                        components=components, orf='hd', name='gwb')


    s += crn

    for name in cw_block_names:
        s += deterministic.cw_block_circ(
            amp_prior="log-uniform",
            dist_prior=None,
            skyloc=None,
            log10_fgw=None,
            psrTerm=True,
            phase_connected=False,
            discoclone=False,
            name=name,
        )
    # s += cw 

    # s += cw2

    models.append(s(p))
    print(p.name)
pta = signal_base.PTA(models)

pta.set_default_params({})

if args.seed is not None:
    np.random.seed(args.seed)

def get_correct_p_phase(psr, gwphi):
    """Get correct pulsar phase offset for given pulsar and GW phi."""
    #const.c = 2.99792458e8  # speed of light [m/s]
    omega_gw = 2.0 * jnp.pi * 10**enterprise_params['cw_log10_fgw']
    dec = jnp.arcsin(enterprise_params['cw_sindec'])
    ra = enterprise_params['cw_ra']
    omega_hat = jnp.array([ -jnp.cos(dec) * jnp.cos(ra), 
                           -jnp.cos(dec) * jnp.sin(ra),
                           -jnp.sin(dec)
                         ])
    pos = jnp.array(psr.pos)
    p_dist_m = psr.pdist[0] * const.kpc  # Convert kpc to meters
    phi_psr = (p_dist_m / const.c) * omega_gw  * (1.0 + jnp.dot(omega_hat, pos))

    return phi_psr

enterprise_params ={}
enterprise_params.update({p: np.random.uniform(-18,-11) for p in pta.param_names if 'rednoise_log10_A' in p})
enterprise_params.update({p: np.random.uniform(1,7) for p in pta.param_names if 'rednoise_gamma' in p})

enterprise_params.update({'gwb_gamma': 4.333, 'gwb_log10_A': -14.5})
for name in cw_block_names:
    enterprise_params.update(
        {
            f"{name}_cos_inc": np.random.uniform(-1.0,1.0),
            f"{name}_log10_fgw": np.random.uniform(-9,-7.0),
            #f"{name}_log10_fgw": -8,
            f"{name}_log10_h": np.random.uniform(-18.0, -11.0),
            #f"{name}_log10_h": -12.5,
            f"{name}_log10_mc": np.random.uniform(6.0, 10.0),
            f"{name}_gwphi": np.random.uniform(0, 2 * np.pi),
            f"{name}_cos_gwtheta": np.random.uniform(-1, 1),
            f"{name}_phase0": np.random.uniform(0, 2 * np.pi),
            f"{name}_psi": np.random.uniform(0, np.pi),
        }
    )

for psr in psrs:
    enterprise_params.update({p: psr.pdist[0] for p in pta.param_names if psr.name+'_cw_p_dist' in p})
    enterprise_params.update({p: get_correct_p_phase(psr, enterprise_params['cw_gwphi']) for p in pta.param_names if psr.name+'_cw_p_phase' in p})

temp_dict = dict(enterprise_params)
reverse_temp_dict = {k: k for k in enterprise_params}

for name, suffix in zip(cw_block_names, disco_suffixes):
    for key, val in enterprise_params.items():
        if key.startswith(f"{name}_"):
            base = key.replace(name, "cw", 1)
            disco_key = f"{base}{suffix}"
            temp_dict[disco_key] = val
            reverse_temp_dict[disco_key] = key

start_real = args.start_realisation
n_real = args.n_real


disco_psrs = [ds.Pulsar.read_feather(f) for f in sorted(glob.glob(feather_dir + "*.feather"))][:Npulsars] 
for psr in disco_psrs:
    psr.toaerrs = np.full_like(psr.toas, 1e-6, dtype=np.float64)

# white-noise parameters fixed once
noisedict = {
    psr.name + "_KAT_MKBF_efac": 1.0
    for psr in disco_psrs
}
noisedict.update({
    psr.name + "_KAT_MKBF_log10_ecorr": -8.0
    for psr in disco_psrs
})
noisedict.update({
    psr.name + "_KAT_MKBF_log10_t2equad": -8.0
    for psr in disco_psrs
})

noise_terms = {
    psr.name: ds.makenoise_measurement(psr, noisedict=noisedict)
    for psr in disco_psrs
}
timing_terms = {
    psr.name: ds.makegp_timing(psr, variance=1e-14)
    for psr in disco_psrs
}


class MultiSourceDelay:
    global_params = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi")

    def __init__(self, psr, n_sources, include_pterm):
        self.psr = psr
        self.n_sources = n_sources
        self.include_pterm = include_pterm
        self.single_source = make_phase_unconnected_binary(pulsarterm=self.include_pterm)
        params = []

        # CW global parameters shared across pulsars
        for idx in range(n_sources):
            suffix = "" if idx == 0 else f"_{idx+1}"
            for key in self.global_params:
                name = f"cw_{key}{suffix}"
                if name not in params:
                    params.append(name)

        # Pulsar-specific distance offsets (one per source if requested)
        if self.include_pterm:
            for idx in range(n_sources):
                suffix = "" if idx == 0 else f"_{idx+1}"
                params.append(f"{self.psr.name}_cw_p_dist")
                params.append(f"{self.psr.name}_cw_p_phase")

        self.params = params

    def __call__(self, params):
        total = 0.0
        for idx in range(self.n_sources):
            suffix = "" if idx == 0 else f"_{idx+1}"
            kwargs = {
                "cos_gwtheta": params[f"cw_cos_gwtheta{suffix}"],
                "gwphi": params[f"cw_gwphi{suffix}"],
                "cos_inc": params[f"cw_cos_inc{suffix}"],
                "log10_mc": params[f"cw_log10_mc{suffix}"],
                "log10_fgw": params[f"cw_log10_fgw{suffix}"],
                "log10_h": params[f"cw_log10_h{suffix}"],
                "phase0": params[f"cw_phase0{suffix}"],
                "psi": params[f"cw_psi{suffix}"],
            }
            if self.include_pterm:
                kwargs["p_dist"] = params[f"{self.psr.name}_cw_p_dist"]
                kwargs["p_phase"] = params[f"{self.psr.name}_cw_p_phase"]
            total += self.single_source(self.psr.toas, self.psr.pos, **kwargs)
        return total
    
T = ds.getspan(disco_psrs)

common_gp = ds.makecommongp_fourier(disco_psrs, ds.powerlaw, 30, T, name="rednoise")
global_gp = ds.makeglobalgp_fourier(
    disco_psrs,
    ds.powerlaw,
    ds.hd_orf,
    30,
    T,
    name="gwb",
)

chunk_size = 128

def build_likelihood(residual_map):
    """Construct ArrayLikelihood and return logl callable with parameter metadata."""
    pulsar_likes = [
        ds.PulsarLikelihood(
            [
                np.array(residual_map[psr.name], copy=True),
                noise_terms[psr.name],
                timing_terms[psr.name],
                MultiSourceDelay(psr, num_cw, include_pterm=args.pulsar_term),
            ]
        )
        for psr in disco_psrs
    ]
    fml = ds.ArrayLikelihood(
        pulsar_likes,
        commongp=common_gp,
        globalgp=global_gp,
    )
    logl = fml.logL

    order_map = {k: i for i, k in enumerate(logl.params)}
    sorted_items = sorted(temp_dict.items(), key=lambda kv: order_map.get(kv[0], float("inf")))
    sorted_dict = {k: v for k, v in sorted_items if k in logl.params}
    param_keys = list(sorted_dict.keys())
    base_vals = jnp.array([sorted_dict[k] for k in param_keys], dtype=jnp.float64)

    return logl, param_keys, base_vals

def eval_grid(eval_fn, template, scan_idx):
    n_chunks = int(np.ceil(scan_values.shape[0] / chunk_size))
    pieces = []
    for idx in range(n_chunks):
        start = idx * chunk_size
        end = min(start + chunk_size, scan_values.shape[0])
        current = scan_values[start:end]
        pad = chunk_size - current.shape[0]
        padded = current if pad == 0 else jnp.pad(current, (0, pad), constant_values=current[-1])
        x_block = template.at[:, scan_idx].set(padded)
        block_vals = eval_fn(x_block)
        pieces.append(np.array(block_vals[: end - start]))
    return np.concatenate(pieces)

all_logls = np.empty((n_real, gridsteps), dtype=np.float64)
start_time = time.time()

param_keys = None
scan_param_index = None
pterm_tag = "psrTerm" if args.pulsar_term else "earthTerm"
combined_prefix = f"{args.save_prefix}_{scan_param}_{pterm_tag}"
out_dir = Path(args.output_root) / args.save_prefix / scan_param / pterm_tag
out_dir.mkdir(parents=True, exist_ok=True)

for offset in tqdm(range(n_real), desc="Realisations", unit="run"):
    if args.seed is not None:
        np.random.seed(args.seed + offset + 1)
    else:
        np.random.seed(None)

    sim_resids = simulate(pta, enterprise_params, sparse_cholesky=True)
    name_to_resid = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim_resids)}

    missing = [psr.name for psr in disco_psrs if psr.name not in name_to_resid]
    if missing:
        raise RuntimeError(f"Pulsars missing from simulate(): {missing}")

    logl, keys, base_values = build_likelihood(name_to_resid)

    if param_keys is None:
        param_keys = keys
        scan_param_index = keys.index(scan_param)
    elif keys != param_keys:
        raise RuntimeError("Parameter ordering changed between realisations.")

    def logl_wrapped(x_array):
        params = {k: v for k, v in zip(param_keys, x_array)}
        return logl(params)
    

    batched_logl = jax.jit(jax.vmap(logl_wrapped), donate_argnums=(0,))
    template = jnp.repeat(base_values[None, :], chunk_size, axis=0)

    logls = eval_grid(batched_logl, template, scan_param_index)
    all_logls[offset] = logls

    if args.save_individual:
        ireal = start_real + offset
        np.save(out_dir / f"{combined_prefix}_realisation_{ireal}.npy", logls)

elapsed = time.time() - start_time
out_path = out_dir / f"{combined_prefix}_{start_real}_{start_real + n_real - 1}.npz"
np.savez(
    out_path,
    scan_values=np.array(scan_values),
    logls=all_logls,
    scan_param=scan_param,
    pulsar_term=args.pulsar_term,
)
print(f"Saved {n_real} realisations to {out_path} in {elapsed:.1f} s")

plot_key = reverse_temp_dict.get(scan_param, scan_param)

fig, ax = plt.subplots(figsize=(8, 5))
for curve in all_logls:
    ax.plot(scan_values, curve, alpha=0.4, linewidth=0.8)

likelihood_deltas = []
low_idx = int(np.argmin(scan_values))
for curve in all_logls:
    peak_idx = int(np.argmax(curve))
    likelihood_deltas.append(curve[peak_idx] - curve[low_idx])

legend_handles = []
legend_labels = []

if plot_key in enterprise_params:
    inj_line = ax.axvline(
        enterprise_params[plot_key], color="k", linestyle="--", label="injected value", alpha=0.4,
    )
    legend_handles.append(inj_line)
    legend_labels.append("injected value")

if likelihood_deltas:
    if len(likelihood_deltas) > 1:
        delta_value = float(np.mean(likelihood_deltas))
        delta_label = f"delta_logL_avg(peak-low): {delta_value:.3f}"
    else:
        delta_label = f"delta_logL(peak-low): {likelihood_deltas[0]:.3f}"
    dummy_handle = ax.plot([], [], " ", label=delta_label)[0]
    legend_handles.append(dummy_handle)
    legend_labels.append(dummy_handle.get_label())

if legend_handles:
    ax.legend(legend_handles, legend_labels)
ax.set_xlabel(scan_param)
ax.set_ylabel("Log-Likelihood")
ax.set_title(f"LogL realisations vs {combined_prefix}_{start_real}_{start_real + n_real - 1} phase_connected")
ax.grid(True, linestyle="--", alpha=0.3)
fig.tight_layout()
fig.savefig(out_dir / f"{combined_prefix}_{start_real}_{start_real + n_real - 1}.png", dpi=150)
plt.close(fig)
