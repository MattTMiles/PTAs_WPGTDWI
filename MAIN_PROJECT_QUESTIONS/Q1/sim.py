"""Utilities to assemble Q1 ring-PTA simulations and discovery log-likelihoods.

The goal is to mirror the exploratory script in the project root while making
it reusable for MAP/optimisation runs inside MAIN_PROJECT_QUESTIONS/Q1.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import jax
import jax.numpy as jnp
import numpy as np

import discovery as ds
from enterprise.signals import selections, signal_base
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import utils
from enterprise.signals.gp_signals import get_timing_model_basis, BasisGP
from enterprise.signals.parameter import function
import enterprise.signals.parameter as parameter
from enterprise_extensions.blocks import common_red_noise_block
from enterprise_extensions import load_feathers
from enterprise_extensions import deterministic
from sksparse.cholmod import cholesky
import scipy.linalg as sl
import scipy.sparse as ss
from discovery.deterministic import make_phase_connected_binary
from discovery import deterministic as det

# Silence XLA preallocation to avoid large grabs on first JIT call
import os
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
os.environ.setdefault("JAX_DISABLE_MMAP_CACHE", "1")
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=2")
jax.config.update("jax_enable_x64", True)


@dataclass
class DistanceTier:
    """Simple container for distance-prior inflation/deflation."""

    name: str
    scale: float  # multiply sigma by this factor (scale < 1 tightens)


@dataclass
class Q1Config:
    """Configuration knobs for Q1 scenarios."""

    n_pulsars: int = 30
    n_cw_sources: int = 1
    ring_radius_deg: float = 20.0
    same_loc: bool = False
    same_freq: bool = False
    include_pterm: bool = True
    inj_evolve: bool = True
    rec_evolve: bool = True
    scenario: str = "A"  # A: WN, B: WN+RN, C: WN+RN+GWB
    distance_tier: DistanceTier = dataclasses.field(default_factory=lambda: DistanceTier("real", 1.0))
    seed: int | None = None
    feather_dir: Path = dataclasses.field(
        default_factory=lambda: Path(__file__).resolve().parents[2] / "data_products"
    )


def make_rng(seed: int | None) -> np.random.Generator:
    if seed is None:
        return np.random.default_rng(np.random.SeedSequence().entropy)
    return np.random.default_rng(seed)


def unit_vector(phi: float, theta: float) -> np.ndarray:
    st = np.sin(theta)
    return np.array([st * np.cos(phi), st * np.sin(phi), np.cos(theta)])


def angles_from_vector(vec: np.ndarray) -> Tuple[float, float]:
    vec = vec / np.linalg.norm(vec)
    theta = np.arccos(np.clip(vec[2], -1.0, 1.0))
    phi = np.mod(np.arctan2(vec[1], vec[0]), 2.0 * np.pi)
    return phi, theta


def build_ring_positions(cw_dirs: Sequence[np.ndarray], n_psr: int, radius: float) -> np.ndarray:
    """Place pulsars on a ring around the mean CW direction."""
    avg = np.mean(cw_dirs, axis=0)
    norm = np.linalg.norm(avg)
    if norm < 1e-8:
        avg = cw_dirs[0]
    else:
        avg = avg / norm

    # pick a reference perpendicular to avg
    ref = np.array([0.0, 0.0, 1.0]) if abs(avg[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    u = np.cross(avg, ref)
    u = u / np.linalg.norm(u)
    v = np.cross(avg, u)

    positions = np.zeros((n_psr, 3))
    for idx in range(n_psr):
        beta = 2.0 * np.pi * idx / n_psr
        tangent = np.cos(beta) * u + np.sin(beta) * v
        vec = np.cos(radius) * avg + np.sin(radius) * tangent
        positions[idx] = vec / np.linalg.norm(vec)
    return positions


def apply_ring_positions(psr_list: Sequence, ring_vecs: np.ndarray) -> None:
    """Mutate pulsar objects (enterprise or discovery) to new sky positions."""
    for psr, vec in zip(psr_list, ring_vecs):
        phi, theta = angles_from_vector(vec)
        psr.pos = vec.copy()
        psr.phi = phi
        psr.theta = theta


def draw_cw_population(cfg: Q1Config, rng: np.random.Generator) -> List[Dict[str, float]]:
    cw_draws: List[Dict[str, float]] = []
    for idx in range(cfg.n_cw_sources):
        suffix = "" if idx == 0 else f"_{idx+1}"
        cw_draws.append(
            {
                "suffix": suffix,
                "cos_inc": rng.uniform(-1.0, 1.0),
                "log10_h": -12.0,
                "log10_mc": 10.0,
                "phase0": rng.uniform(0.0, 2.0 * np.pi),
                "psi": rng.uniform(0.0, np.pi),
            }
        )
        if cfg.same_loc:
            cw_draws[-1]["gwphi"] = np.pi
            cw_draws[-1]["cos_gwtheta"] = 0.0
        else:
            cw_draws[-1]["gwphi"] = rng.uniform(0.0, 2.0 * np.pi)
            cw_draws[-1]["cos_gwtheta"] = rng.uniform(-1.0, 1.0)

        if cfg.same_freq:
            cw_draws[-1]["log10_fgw"] = -8.0
        else:
            cw_draws[-1]["log10_fgw"] = rng.uniform(-9.0, -7.5)
    return cw_draws


def compute_ring_radius(cw_draws: Sequence[Dict[str, float]], cfg: Q1Config) -> float:
    """Choose a ring radius that comfortably surrounds all CW directions."""
    ring_radius = np.deg2rad(cfg.ring_radius_deg)
    cw_dirs = [unit_vector(draw["gwphi"], np.arccos(draw["cos_gwtheta"])) for draw in cw_draws]
    avg = np.mean(cw_dirs, axis=0)
    avg /= np.linalg.norm(avg)
    distances = np.array([np.arccos(np.clip(np.dot(avg, d), -1.0, 1.0)) for d in cw_dirs])
    required = distances.max() + np.deg2rad(5.0)
    return max(ring_radius, required)


def load_enterprise_pulsars(cfg: Q1Config) -> List:
    psrs = load_feathers.load_feathers_from_folder(str(cfg.feather_dir))
    return psrs[: cfg.n_pulsars]


def load_discovery_pulsars(cfg: Q1Config):
    files = sorted(Path(cfg.feather_dir).glob("*.feather"))
    return [ds.Pulsar.read_feather(f) for f in files[: cfg.n_pulsars]]


def set_toaerrs(psrs: Sequence, sigma: float = 1e-6) -> None:
    """Force uniform TOA uncertainties to avoid outliers during injection."""
    for psr in psrs:
        psr.toaerrs = np.full_like(psr.toaerrs, sigma, dtype=np.float64)


def apply_distance_tier(psrs: Sequence, tier: DistanceTier) -> Dict[str, float]:
    """Update pulsar distance priors based on tier; return true distances."""
    true = {}
    for psr in psrs:
        true[psr.name + "_cw_p_dist"] = psr.pdist[0]
        # enterprise uses [mean, sigma]; rescale sigma
        if tier.scale == 0.0:
            psr.pdist = [psr.pdist[0], 1e-6]
        else:
            psr.pdist = [psr.pdist[0], psr.pdist[1] * tier.scale]
        psr._pdist = psr.pdist
        psr.residuals = np.array(psr.toas) * 0.0
    return true


@function
def tm_prior(weights, toas, variance=1e40):
    return weights * variance * len(toas)


def TimingModel(coefficients=False, name="linear_timing_model", use_svd=False, normed=True, prior_variance=1e40):
    """Create a timing model GP with the same defaults as the root script."""
    basis = get_timing_model_basis(use_svd, normed)
    prior = tm_prior(variance=prior_variance)
    BaseClass = BasisGP(prior, basis, coefficients=coefficients, name=name)

    class _TimingModel(BaseClass):
        signal_type = "basis"
        signal_name = "linear timing model"
        signal_id = name + "_svd" if use_svd else name

    return _TimingModel


def build_enterprise_pta(psrs: Sequence, cfg: Q1Config, cw_draws: Sequence[Dict[str, float]]):
    """Construct enterprise PTA model mirroring the root script."""
    selection = selections.Selection(selections.by_backend)
    efac = parameter.Constant(1)
    equad = parameter.Constant(-8)

    log10_A_red = parameter.Uniform(-18, -11)
    gamma_red = parameter.Uniform(0, 7)
    components = 30

    tmin = [p.toas.min() for p in psrs]
    tmax = [p.toas.max() for p in psrs]
    Tspan = np.max(tmax) - np.min(tmin)

    ring_radius = compute_ring_radius(cw_draws, cfg)
    cw_dirs = [unit_vector(draw["gwphi"], np.arccos(draw["cos_gwtheta"])) for draw in cw_draws]
    ring_vecs = build_ring_positions(cw_dirs, cfg.n_pulsars, ring_radius)
    apply_ring_positions(psrs, ring_vecs)

    tm = TimingModel(coefficients=False, name="linear_timing_model", use_svd=False, normed=True, prior_variance=1e-14)

    models = []
    for p in psrs:
        s = tm
        s += white_signals.MeasurementNoise(efac=efac, selection=selection)
        s += white_signals.TNEquadNoise(log10_tnequad=equad, selection=selection)

        if cfg.scenario in ("B", "C"):
            pl = utils.powerlaw(log10_A=log10_A_red, gamma=gamma_red)
            rn = gp_signals.FourierBasisGP(spectrum=pl, components=components, Tspan=Tspan, name="rednoise")
            s += rn

        if cfg.scenario == "C":
            crn = common_red_noise_block(psd="powerlaw", prior="log-uniform", components=components, orf="hd", name="gwb")
            s += crn

        for idx in range(cfg.n_cw_sources):
            name = "cw" if idx == 0 else f"cw{idx+1}"
            s += deterministic.cw_block_circ(
                amp_prior="log-uniform",
                dist_prior=None,
                skyloc=None,
                log10_fgw=None,
                psrTerm=True,
                phase_connected=True,
                discoclone=False,
                name=name,
                evolve=cfg.inj_evolve,
            )
        models.append(s(p))
    pta = signal_base.PTA(models)
    return pta, ring_vecs


def draw_noise_params(pta, cfg: Q1Config, rng: np.random.Generator, psrs: Sequence | None = None) -> Dict[str, float]:
    """Sample noise parameters conditioned on scenario."""
    params: Dict[str, float] = {}
    if cfg.scenario == "A":
        params.update({p: -18 for p in pta.param_names if "rednoise_log10_A" in p})
        params.update({p: 3 for p in pta.param_names if "rednoise_gamma" in p})
        params.update({"gwb_gamma": 4.333, "gwb_log10_A": -18.5})
    else:
        params.update({p: rng.uniform(-17.0, -13.0) for p in pta.param_names if "rednoise_log10_A" in p})
        params.update({p: rng.uniform(2.0, 6.0) for p in pta.param_names if "rednoise_gamma" in p})
        params.update({"gwb_gamma": 4.333, "gwb_log10_A": -14.5})
    psr_list = psrs if psrs is not None else getattr(pta, "pulsars", [])
    for psr in psr_list:
        name = getattr(psr, "name", psr)
        pdist0 = getattr(psr, "pdist", [1.0])[0] if hasattr(psr, "pdist") else 1.0
        params.update({p: pdist0 for p in pta.param_names if f"{name}_cw_p_dist" in p})
    return params


def inject_cw_params(params: Dict[str, float], cw_draws: Sequence[Dict[str, float]], cfg: Q1Config) -> Dict[str, float]:
    params = dict(params)
    for idx, draw in enumerate(cw_draws):
        name = "cw" if idx == 0 else f"cw{idx+1}"
        params.update(
            {
                f"{name}_cos_inc": draw["cos_inc"],
                f"{name}_log10_fgw": draw["log10_fgw"],
                f"{name}_log10_h": draw["log10_h"],
                f"{name}_log10_mc": draw["log10_mc"],
                f"{name}_gwphi": draw["gwphi"],
                f"{name}_cos_gwtheta": draw["cos_gwtheta"],
                f"{name}_phase0": draw["phase0"],
                f"{name}_psi": draw["psi"],
            }
        )
    return params


def simulate_enterprise_residuals(pta, params: Dict[str, float], sparse_cholesky: bool = True) -> List[np.ndarray]:
    """Reimplementation of simulate() from the root script with minimal tweaks."""
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
            whiteresiduals.append(delay)
        elif isinstance(ndiag, signal_base.ShermanMorrison):
            n = np.diag(ndiag._nvec)
            for j, s in zip(ndiag._jvec, ndiag._slices):
                n[s, s] += j
            whiteresiduals.append(delay + np.dot(sl.cholesky(n, lower=True), np.random.randn(n.shape[0])))
        elif ndiag.ndim == 1:
            whiteresiduals.append(delay + np.sqrt(ndiag) * np.random.randn(ndiag.shape[0]))
        else:
            raise NotImplementedError
    return [np.array(g + w) for g, w in zip(gpresiduals, whiteresiduals)]


def build_discovery_loglike(
    disco_psrs: Sequence,
    cw_draws: Sequence[Dict[str, float]],
    cfg: Q1Config,
    ring_vecs: np.ndarray | None = None,
    residual_map: Dict[str, np.ndarray] | None = None,
):
    """Return discovery logL callable + base param ordering."""
    for psr in disco_psrs:
        psr.toaerrs = np.full_like(psr.toas, 1e-6, dtype=np.float64)
        psr.pdist = 1
        if residual_map is not None and psr.name in residual_map:
            # Keep the discovery Pulsar object's residuals in sync with the simulated data
            psr.residuals = np.array(residual_map[psr.name], copy=True)

    if ring_vecs is None:
        radius = compute_ring_radius(cw_draws, cfg)
        ring_vecs = build_ring_positions(
            [unit_vector(d["gwphi"], np.arccos(d["cos_gwtheta"])) for d in cw_draws],
            cfg.n_pulsars,
            radius,
        )
    apply_ring_positions(disco_psrs, ring_vecs)

    noise_terms = {psr.name: ds.makenoise_measurement(psr, noisedict={
        psr.name + "_KAT_MKBF_efac": 1.0,
        psr.name + "_KAT_MKBF_log10_ecorr": -8.0,
        psr.name + "_KAT_MKBF_log10_t2equad": -8.0,
    }) for psr in disco_psrs}
    timing_terms = {psr.name: ds.makegp_timing(psr, variance=1e-14) for psr in disco_psrs}

    class MultiSourceDelay:
        global_params = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi")

        def __init__(self, psr, n_sources, include_pterm, include_evolve):
            self.psr = psr
            self.n_sources = n_sources
            self.include_pterm = include_pterm
            self.include_evolve = include_evolve
            self.single_source = make_phase_connected_binary(pulsarterm=self.include_pterm, evolve=self.include_evolve)
            params = []
            for idx in range(n_sources):
                suffix = "" if idx == 0 else f"_{idx+1}"
                for key in self.global_params:
                    name = f"cw_{key}{suffix}"
                    if name not in params:
                        params.append(name)
                if self.include_pterm:
                    params.append(f"{self.psr.name}_cw_p_dist")
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
                total += self.single_source(self.psr.toas, self.psr.pos, **kwargs)
            return total

    T = ds.getspan(disco_psrs)
    common_gp = None
    global_gp = None
    if cfg.scenario in ("B", "C"):
        common_gp = ds.makecommongp_fourier(disco_psrs, ds.powerlaw, 30, T, name="rednoise")
    if cfg.scenario == "C":
        global_gp = ds.makeglobalgp_fourier(disco_psrs, ds.powerlaw, ds.hd_orf, 30, T, name="gwb")

    delays = [
        MultiSourceDelay(psr, cfg.n_cw_sources, include_pterm=cfg.include_pterm, include_evolve=cfg.rec_evolve)
        for psr in disco_psrs
    ]
    pulsar_likes = [
        ds.PulsarLikelihood(
            [
                np.array(np.zeros_like(psr.toas) if residual_map is None else residual_map[psr.name], copy=True),
                noise_terms[psr.name],
                timing_terms[psr.name],
                delay,
            ]
        )
        for psr, delay in zip(disco_psrs, delays)
    ]

    # Force JAX-friendly logL by keeping y callable inside make_kernelproduct
    for pl in pulsar_likes:
        if callable(pl.y):
            kp = pl.N.make_kernelproduct(pl.y)

            def _logl(params, kp=kp):
                return kp(params)

            _logl.params = getattr(kp, "params", [])
            pl.logL = _logl

    fml = ds.ArrayLikelihood(pulsar_likes, commongp=common_gp, globalgp=global_gp)
    logl = fml.logL
    param_keys = list(getattr(logl, "params", ()))
    if not param_keys:
        gathered = []
        for delay in delays:
            gathered.extend(getattr(delay, "params", ()))
        if common_gp is not None and hasattr(common_gp, "params"):
            gathered.extend(list(common_gp.params))
        if global_gp is not None and hasattr(global_gp, "params"):
            gathered.extend(list(global_gp.params))
        param_keys = list(dict.fromkeys(gathered))  # preserve order, drop dups
    return logl, param_keys


def build_discovery_loglike_marginalized(
    disco_psrs: Sequence,
    cw_draws: Sequence[Dict[str, float]],
    cfg: Q1Config,
    dist_priors: Dict[str, Tuple[float, float]],
    ring_vecs: np.ndarray | None = None,
    residual_map: Dict[str, np.ndarray] | None = None,
    include_cw: bool = True,
):
    """Return discovery logL with analytic distance marginalisation (no distance params).

    The pulsar term is evaluated at the prior mean distance and attenuated by
    kappa = exp(-0.5 * (omega * (1-cos_mu) * kpc/c * sigma_d)^2).
    This removes the distance parameters from the model while retaining their
    effect through the coherence factor.
    """
    for psr in disco_psrs:
        psr.toaerrs = np.full_like(psr.toas, 1e-6, dtype=np.float64)
        if residual_map is not None and psr.name in residual_map:
            psr.residuals = np.array(residual_map[psr.name], copy=True)

    if ring_vecs is None:
        radius = compute_ring_radius(cw_draws, cfg)
        ring_vecs = build_ring_positions(
            [unit_vector(d["gwphi"], np.arccos(d["cos_gwtheta"])) for d in cw_draws],
            cfg.n_pulsars,
            radius,
        )
    apply_ring_positions(disco_psrs, ring_vecs)

    noise_terms = {
        psr.name: ds.makenoise_measurement(
            psr,
            noisedict={
                psr.name + "_KAT_MKBF_efac": 1.0,
                psr.name + "_KAT_MKBF_log10_ecorr": -8.0,
                psr.name + "_KAT_MKBF_log10_t2equad": -8.0,
            },
        )
        for psr in disco_psrs
    }
    timing_terms = {psr.name: ds.makegp_timing(psr, variance=1e-14) for psr in disco_psrs}

    class MarginalizedDelay:
        global_params = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc", "log10_fgw", "log10_h", "phase0", "psi")

        def __init__(self, psr, n_sources, include_pterm, include_evolve):
            self.psr = psr
            self.n_sources = n_sources
            self.include_pterm = include_pterm
            self.include_evolve = include_evolve
            self.earth_term = make_phase_connected_binary(pulsarterm=False, evolve=self.include_evolve)
            self.full_term = make_phase_connected_binary(pulsarterm=True, evolve=self.include_evolve)
            params = []
            for idx in range(n_sources):
                suffix = "" if idx == 0 else f"_{idx+1}"
                for key in self.global_params:
                    name = f"cw_{key}{suffix}"
                    if name not in params:
                        params.append(name)
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
                # Earth term only (pulsar term off)
                earth = self.earth_term(self.psr.toas, self.psr.pos, **kwargs)
                if not self.include_pterm:
                    total += earth
                    continue

                # Pulsar term at mean distance
                d0, sigma = dist_priors.get(self.psr.name, (None, None))
                if d0 is None:
                    raise ValueError(f"Missing distance prior for pulsar {self.psr.name}")
                full = self.full_term(self.psr.toas, self.psr.pos, p_dist=d0, **kwargs)
                pulsar_part = full - earth

                # Coherence factor kappa = exp(-0.5 * (w0 * (1-cos_mu) * kpc/c * sigma)^2)
                gwtheta = jnp.arccos(kwargs["cos_gwtheta"])
                fplus, fcross, cos_mu = det.fpcmu_fast(self.psr.pos, gwtheta, kwargs["gwphi"])
                w0 = jnp.pi * (10.0 ** kwargs["log10_fgw"])
                # Use full GW phase factor ~ 2*pi*f * (1 + cos_mu) * d/c
                sigma_phi = 2.0 * w0 * (1.0 - cos_mu) * (det.const.kpc / det.const.c) * sigma
                # Coherence factor kappa = exp(-0.5 * [2*w0 * (1-cos_mu) * (kpc/c) * sigma]^2)
                kappa = jnp.exp(-0.5 * sigma_phi ** 2) if sigma is not None else 1.0

                total += earth + kappa * pulsar_part
            return total

    T = ds.getspan(disco_psrs)
    common_gp = None
    global_gp = None
    if cfg.scenario in ("B", "C"):
        common_gp = ds.makecommongp_fourier(disco_psrs, ds.powerlaw, 30, T, name="rednoise")
    if cfg.scenario == "C":
        global_gp = ds.makeglobalgp_fourier(disco_psrs, ds.powerlaw, ds.hd_orf, 30, T, name="gwb")

    delays = []
    if include_cw:
        delays = [
            MarginalizedDelay(psr, cfg.n_cw_sources, include_pterm=cfg.include_pterm, include_evolve=cfg.rec_evolve)
            for psr in disco_psrs
        ]
    else:
        delays = [None for _ in disco_psrs]

    pulsar_likes = []
    for psr, delay in zip(disco_psrs, delays):
        signals = [
            np.array(np.zeros_like(psr.toas) if residual_map is None else residual_map[psr.name], copy=True),
            noise_terms[psr.name],
            timing_terms[psr.name],
        ]
        if include_cw:
            signals.append(delay)
        pulsar_likes.append(ds.PulsarLikelihood(signals))

    # Force JAX-friendly logL by keeping y callable inside make_kernelproduct
    for pl in pulsar_likes:
        if callable(pl.y):
            kp = pl.N.make_kernelproduct(pl.y)

            def _logl(params, kp=kp):
                return kp(params)

            _logl.params = getattr(kp, "params", [])
            pl.logL = _logl

    fml = ds.ArrayLikelihood(pulsar_likes, commongp=common_gp, globalgp=global_gp)
    logl = fml.logL
    param_keys = list(getattr(logl, "params", ()))
    if not param_keys:
        gathered = []
        for delay in delays:
            if delay is None:
                continue
            gathered.extend(getattr(delay, "params", ()))
        if common_gp is not None and hasattr(common_gp, "params"):
            gathered.extend(list(common_gp.params))
        if global_gp is not None and hasattr(global_gp, "params"):
            gathered.extend(list(global_gp.params))
        param_keys = list(dict.fromkeys(gathered))  # preserve order, drop dups
    return logl, param_keys
