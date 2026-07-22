"""Track B — B1 CORE: the noisy, targeted (counterpart-conditioned) machinery.

B1 tests the CONDITIONAL pipeline, not cold-start certification (Track B terminal verdict:
cold start is an INFORMATION-THEORETIC NO-GO, f=6.9e-7, break-even sky prior 0.188 deg/src).
The live scenario is TARGETED: the loud sources' SKY is supplied externally (an EM
counterpart clears the 22-arcsec registration tolerance by >=1e4). Everything else the
pipeline must find. Everything before B1 was zero-noise Asimov; B1 draws real noise.

WHAT THIS MODULE ADDS over the banked Track B machinery
-------------------------------------------------------
1. **A per-pulsar PULSAR-TERM MASK as a RUNTIME argument.** `MaskedDelay` evaluates both the
   full (Earth+pulsar) and Earth-only waveform and returns
       delay_p = d_earth + m_p * (d_full - d_earth),      m_p in {0,1}
   so ONE compiled graph serves (a) the Earth-only source fit, (b) the phase-up fit with the
   pulsar terms of a certified subset switched on, and (c) the full model. This is what the
   B1.3 loop-gain experiment needs and it is why B1 does not reuse `stagec_fisher`'s two
   separate (full / earth) amortised builds.

2. **DATA as a runtime argument all the way down.** `stagec_fisher.build_fisher_amortised`
   already does this for the array likelihood, but `trackB_split.Split` BAKES `data_obs` into
   its jitted per-pulsar evaluators (`self.data` is read at trace time). With >=20 noise
   realisations that would mean >=20 compiles of everything. `B1Split` threads (data, pmask)
   through as runtime args, so the ~500 s compile is paid ONCE for the whole ensemble.

3. **Real noise draws**: heterogeneous white (the exact `N` diagonal the likelihood uses) +
   per-pulsar red-noise GP (30 comp, frozen hyperparams) + the HD-correlated GWB GP (14 comp,
   frozen at NG15). All three are in the likelihood's noise covariance, so all three are
   drawn. The linear timing-model GP is NOT drawn (it is a marginalisation device, not noise).

ARM A vs ARM B (the distance-truth convention) — decided with Matt, 2026-07-08
------------------------------------------------------------------------------
`generate_injection_params` (cw_helpers.py:228) injects every true distance AT the EM prior
mean, so n_true == 0 for all 116 pulsars. Every Track B deliverable to date inherits that.
It (i) makes the Gaussian prior tail maximally favour the true fringe and (ii) leaves the
array pre-registered before any fitting, because the pipeline's own starting distances are L0.

  * ARM A  (continuity): L_true = L0. Directly comparable to the census ceilings
    (J0711 0.953 / J1713 0.989 / J1909 0.998). NOT the headline.
  * ARM B  (HEADLINE): n_true drawn per pulsar from the prior-weighted distribution over the
    SAMPLED fringe centres, and the within-fringe offset drawn uniform in (-1/2, +1/2) dL --
    so neither the integer nor the sub-fringe position is pinned to a convention. This is the
    first Track B computation in which truth is not at the prior mean.

The A-vs-B delta prices what registration-from-the-prior-mean was worth.

Env: jug-gpu, XLA autotune off (shared GPU), x64, jax persistent cache. Matt commits; never
commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, hashlib
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")

import discovery as ds
from discovery import matrix as dsm
from discovery.deterministic import make_phase_connected_binary

from cw_helpers import (load_pulsars, build_enterprise_pta, compute_mode_spacing,
                        _enterprise_to_disco_params)
from stagec_fisher import (make_geometry_injection, N_COMPONENTS, LOG10_EQUAD,
                           CW_PARAM_BASE, KPC_TO_PC, cw_param_keys)
import stagec_anchor_a2 as A2
import trackB_estimator as TE

CWT = "/home/mattm/projects/HSYMT/CW_transition"
CENSUS = ["J0711-6830", "J1713+0747", "J1909-3744"]
ANCHOR = "J0437-4715"                # the sole prior-pinned (K<=3) pulsar (R POSTMORTEM)
NEXT_TIER = ["J2222-0137", "J1744-1134", "J0030+0451", "J1600-3053", "J1017-7156"]

B_CERT = 512                         # certification fringe budget -- NEVER lower
PMASK = "__pmask"                    # runtime per-pulsar pulsar-term mask key
SMASK = "__smask"                    # runtime per-SOURCE feed weight (the source-side twin of
                                     # PMASK; FORGE-G soft source side). ABSENT -> weight 1 for
                                     # every source, so every pre-FORGE_G caller is bit-identical.
                                     # FORGE-G2: a JITTED RUNTIME argument (TURBINE) -- pattern
                                     # changes never recompile; only None <-> array retraces.

# Double-where guard (TURBINE, FORGE-G2). A CARRIED (w_s == 0) source's waveform inputs are
# replaced by this known-evaluable point BEFORE the waveform is evaluated, so neither the
# primal nor any gradient ever touches the carried source's actual params (which may sit in
# the non-evaluable coalescence corner or be junk -- SPARK-3 §2). Values: the benign corner
# of the generative prior box (lowest chirp mass x lowest frequency = longest coalescence
# time), zero angles.
SMASK_SAFE = dict(cos_gwtheta=0.0, gwphi=1.0, cos_inc=0.0, log10_mc=8.5,
                  log10_fgw=-8.0, log10_h=-15.0, phase0=0.0, psi=0.0)

# population config (the saved-seed draw used by every Track B deliverable)
POP = dict(ncw=16, seed=3000, population=(3, -13.25, -14.25))
N_LOUD = 3


# ============================================================
# 1. Masked CW delay: pulsar term switched on per pulsar at RUNTIME
# ============================================================
class MaskedDelay:
    """delay_p = d_earth + m_p * (d_full - d_earth), m_p read from params[PMASK][ipsr].

    m_p = 1 -> exactly MultiSourceDelay(include_pterm=True) (gated 0.00 vs stagec_fisher).
    m_p = 0 -> exactly trackB_estimator.EarthDelay          (gated 0.00).
    Cost: two waveform evaluations per pulsar (the pulsar term is not separately exposed by
    `make_phase_connected_binary`). This ~2x on the delay is what buys ONE compiled graph for
    the whole B1.3 phase-up ladder instead of one compile per certified set.
    """
    _G = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
          "log10_fgw", "log10_h", "phase0", "psi")

    def __init__(self, psr, n_sources, ipsr):
        self.psr = psr
        self.n_sources = n_sources
        self.ipsr = int(ipsr)
        self._full = jax.vmap(make_phase_connected_binary(pulsarterm=True),
                              in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0, None))
        self._earth = jax.vmap(make_phase_connected_binary(pulsarterm=False),
                               in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0))
        self._suf = ["" if i == 0 else f"_{i+1}" for i in range(n_sources)]
        params = []
        for s in self._suf:
            for k in self._G:
                nm = f"cw_{k}{s}"
                if nm not in params:
                    params.append(nm)
        self._pdist_key = f"{psr.name}_cw_p_dist"
        params.append(self._pdist_key)
        params.append(PMASK)
        self.params = params

    def __call__(self, params):
        cw = {k: jnp.stack([params[f"cw_{k}{s}"] for s in self._suf]) for k in self._G}
        L = params[self._pdist_key]
        m = params[PMASK][self.ipsr]
        # Per-SOURCE feed weight w_s (SMASK). w_s == 1 -> the source contributes its full delay;
        # w_s == 0 -> the source is ABSENT from the template (exactly, in fp64) -- the source-side
        # analogue of switching a pulsar term off. ABSENT key -> w_s == 1 for every source, so the
        # 1.0*x == x identity leaves every pre-FORGE_G caller bit-identical. FORGE-G2: w may be a
        # TRACED runtime argument -- nothing here may branch on its values.
        w = params.get(SMASK, None)
        if w is not None:
            gate1 = w > 0.0
            # DOUBLE-WHERE (the TURBINE guard). The outer where below already makes a carried
            # source's VALUE exact zero, but autodiff pushes a (zero) cotangent THROUGH the
            # waveform evaluated at the carried source's ACTUAL params; if those are
            # non-evaluable (junk reservoir, coalescence corner -- SPARK-3 §2), 0 * NaN == NaN
            # poisons the whole gradient. Sanitising the INPUTS means the waveform and its
            # grads are only ever evaluated at SMASK_SAFE for carried sources; a fed source's
            # inputs pass through where() unchanged, so the primal stays bit-exact.
            cw = {k: jnp.where(gate1, cw[k], SMASK_SAFE[k]) for k in self._G}
        args = (self.psr.toas, self.psr.pos, cw["cos_gwtheta"], cw["gwphi"], cw["cos_inc"],
                cw["log10_mc"], cw["log10_fgw"], cw["log10_h"], cw["phase0"], cw["psi"])
        full = self._full(*args, L)                 # (n_src, n_toa)
        earth = self._earth(*args)                   # (n_src, n_toa)
        if w is not None:
            # jnp.where SELECTS 0 for a not-fed source rather than MULTIPLYING by 0. Selection
            # makes w_s == 0 TRUE absence -- a carried source cannot NaN the template. A FED
            # source (w_s > 0) is still evaluated at its actual params and must be valid.
            gate = gate1[:, jnp.newaxis]
            wcol = w[:, jnp.newaxis]
            d_full = jnp.sum(jnp.where(gate, wcol * full, 0.0), axis=0)
            d_earth = jnp.sum(jnp.where(gate, wcol * earth, 0.0), axis=0)
        else:
            d_full = jnp.sum(full, axis=0)
            d_earth = jnp.sum(earth, axis=0)
        return d_earth + m * (d_full - d_earth)


# ============================================================
# 2. Amortised likelihood with (data, pmask) as runtime args
#    (mirrors stagec_fisher.build_fisher_amortised; stagec_fisher is left untouched)
# ============================================================
def build_b1_amortised(disco_psrs, num_cw, enterprise_params, cw_block_names,
                       components=N_COMPONENTS, log10_equad=LOG10_EQUAD, rn_components=30):
    noisedict = {}
    for psr in disco_psrs:
        noisedict[psr.name + "_KAT_MKBF_efac"] = 1.0
        noisedict[psr.name + "_KAT_MKBF_log10_ecorr"] = -8.0
        noisedict[psr.name + "_KAT_MKBF_log10_t2equad"] = log10_equad

    noise_terms = {p.name: ds.makenoise_measurement(p, noisedict=noisedict) for p in disco_psrs}
    timing_terms = {p.name: ds.makegp_timing(p, variance=1e-14) for p in disco_psrs}

    T = ds.getspan(disco_psrs)
    common_gp = ds.makecommongp_fourier(disco_psrs, ds.powerlaw, rn_components, T, name="rednoise")
    global_gp = ds.makeglobalgp_fourier(disco_psrs, ds.powerlaw, ds.hd_orf, components, T, name="gwb")

    temp_dict = _enterprise_to_disco_params(enterprise_params, cw_block_names)
    msds = [MaskedDelay(p, num_cw, i) for i, p in enumerate(disco_psrs)]
    data_keys = [f"{p.name}__data" for p in disco_psrs]

    ys = []
    for p, msd, dkey in zip(disco_psrs, msds, data_keys):
        def y(params, _msd=msd, _dkey=dkey):
            return params[_dkey] - _msd(params)
        y.params = sorted(set(msd.params) | {dkey})
        ys.append(y)

    # FAST PATH for mask == 1. MaskedDelay evaluates BOTH the full and the Earth-only waveform to
    # form earth + m*(full-earth); when the mask is identically one the residual is just
    # data - full(L), so the E-step's 116 x B inner loop pays a 2x tax for a mask it never varies.
    # These residuals are NOT in the likelihood graph -- only B1Split's per-fringe evaluator uses
    # them, and only after `enable_fast_full()`, which is gated bit-exactly against the masked path.
    from cw_helpers import MultiSourceDelay as _MSD
    msds_full = [_MSD(p, num_cw, include_pterm=True) for p in disco_psrs]
    ys_full = []
    for p, msd, dkey in zip(disco_psrs, msds_full, data_keys):
        def yf(params, _msd=msd, _dkey=dkey):
            return params[_dkey] - _msd(params)
        yf.params = sorted(set(msd.params) | {dkey})
        ys_full.append(yf)

    pulsar_likes = [
        ds.PulsarLikelihood([np.zeros(len(p.toas)), noise_terms[p.name],
                             timing_terms[p.name], msd])
        for p, msd in zip(disco_psrs, msds)
    ]
    fml = ds.ArrayLikelihood(pulsar_likes, commongp=common_gp, globalgp=global_gp)
    _ = fml.logL

    P_var_inv = getattr(global_gp, "Phi_inv", None) or global_gp.Phi.make_inv()
    Fs = getattr(global_gp, "Fs", None)
    if Fs is None:
        F = getattr(global_gp, "F", None)
        Fs = F if isinstance(F, (list, tuple)) else [F] * len(pulsar_likes)

    kterms = fml.vsm.make_kernelterms(ys, Fs)
    npsr = len(Fs)
    ngp = Fs[0].shape[1]

    truth_full = {k: temp_dict[k] for k in kterms.params if k in temp_dict}
    for k in P_var_inv.params:
        if k in temp_dict:
            truth_full[k] = temp_dict[k]

    zero_data = {dk: jnp.zeros(len(p.toas)) for dk, p in zip(data_keys, disco_psrs)}
    ones_mask = jnp.ones(npsr)
    terms_truth = kterms({**truth_full, **zero_data, PMASK: ones_mask})
    c_truth = terms_truth[2]                       # distance/data/mask INDEPENDENT
    Pinv_truth, ldP_truth = P_var_inv(truth_full)
    cf_cached = dsm.matrix_factor(Pinv_truth + dsm.jsp.linalg.block_diag(*c_truth))
    ldP_sum_cached = jnp.sum(ldP_truth)
    logdet_cf_cached = dsm.matrix_norm * jnp.sum(jnp.log(jnp.diag(cf_cached[0])))

    dist_keys = [f"{p.name}_cw_p_dist" for p in disco_psrs]
    theta_keys = dist_keys + cw_param_keys(num_cw)
    n_dist = len(dist_keys)
    theta_truth = jnp.array([temp_dict[k] for k in theta_keys], dtype=jnp.float64)
    frozen = {k: temp_dict[k] for k in truth_full if k not in theta_keys}

    def logL(theta_arr, data_tuple, pmask):
        params = dict(frozen)
        for k, v in zip(theta_keys, theta_arr):
            params[k] = v
        for dk, d in zip(data_keys, data_tuple):
            params[dk] = d
        params[PMASK] = pmask
        terms = kterms(params)
        p0 = jnp.sum(terms[0])
        FtNmy = terms[1].reshape(npsr * ngp)
        return p0 + 0.5 * (FtNmy.T @ dsm.matrix_solve(cf_cached, FtNmy)
                           - ldP_sum_cached - logdet_cf_cached)

    def inject_delay(theta_arr, pmask):
        """Per-pulsar CW delay at theta (the noise-free signal)."""
        params = dict(frozen)
        for k, v in zip(theta_keys, theta_arr):
            params[k] = v
        params[PMASK] = pmask
        return tuple(msd(params) for msd in msds)

    return dict(
        logL=jax.jit(logL),
        logL_batch_theta=jax.jit(jax.vmap(logL, in_axes=(0, None, None))),
        inject_delay=jax.jit(inject_delay),
        theta_keys=theta_keys, n_dist=n_dist, n_theta=len(theta_keys),
        theta_truth=np.asarray(theta_truth), npsr=npsr, ngp_gwb=ngp,
        disco_psrs=disco_psrs,
        internals=dict(fml=fml, vsm=fml.vsm, kterms=kterms, P_var_inv=P_var_inv,
                       cf_cached=cf_cached, ldP_sum_cached=ldP_sum_cached,
                       logdet_cf_cached=logdet_cf_cached, c_truth=c_truth,
                       Fs_gwb=Fs, ys=ys, ys_full=ys_full, msds=msds, msds_full=msds_full,
                       frozen=frozen, data_keys=data_keys,
                       npsr=npsr, ngp_gwb=ngp, global_gp=global_gp, common_gp=common_gp,
                       Pinv_gwb=Pinv_truth, truth_full=truth_full),
    )


def _white_diag(N, i):
    """Per-TOA white variance out of vsm.Ns[i] (a WoodburyKernel_novar around the white
    NoiseMatrix1D_novar). Fails loudly rather than silently drawing the wrong covariance."""
    inner = getattr(N, "N", None)
    vec = getattr(inner, "N", None)
    if not isinstance(vec, np.ndarray) or vec.ndim != 1:
        raise RuntimeError(f"pulsar {i}: cannot locate a 1-D white-noise diagonal inside "
                           f"{type(N).__name__} (inner {type(inner).__name__}). The B1 noise "
                           f"draw must match the likelihood's covariance exactly -- aborting.")
    return vec


# ============================================================
# 3. Low-rank split with (data, pmask) runtime  -- the E-step engine
# ============================================================
class B1Split:
    """`trackB_split.Split`, with data + pulsar-term mask threaded as RUNTIME args.

    Frozen pieces (per-pulsar rednoise Woodbury cf_p, TtNmF_p, the outer cf_cached, ldN)
    depend only on the noise/GP hyperparameters -> built once, shared by every realisation.
    `a_const` is data- AND mask-independent (it collects ldN + the rednoise logdet); this is
    gated, not assumed.
    """

    def __init__(self, amo, theta_truth, data_ref, mask_ref):
        it = amo["internals"]
        self.amo = amo
        self.theta_keys = amo["theta_keys"]; self.n_dist = amo["n_dist"]
        self.frozen = it["frozen"]; self.data_keys = it["data_keys"]
        self.ys = it["ys"]; self._ys_masked = it["ys"]; self._ys_full = it["ys_full"]
        self._fast_full = False
        self.cf_cached = it["cf_cached"]; self.ldP_cached = it["ldP_sum_cached"]
        self.logdet_cached = it["logdet_cf_cached"]
        self.npsr = it["npsr"]; self.ngp_gwb = it["ngp_gwb"]
        vsm = it["vsm"]
        Ns, Fs_rn, Ts = vsm.Ns, vsm.Fs, it["Fs_gwb"]

        NmFs, ldNs = zip(*[N.solve_2d(F) for N, F in zip(Ns, Fs_rn)])
        if dsm.regularize_FtNmF:
            FtNmFs = [dsm.zero_out_negative_eigenvalues(F.T @ NmF) for F, NmF in zip(Fs_rn, NmFs)]
        else:
            FtNmFs = [F.T @ NmF for F, NmF in zip(Fs_rn, NmFs)]
        TtNmFs = [T.T @ NmF for T, NmF in zip(Ts, NmFs)]
        self.FtNmF = jnp.asarray(np.array(FtNmFs))
        self.TtNmF = jnp.asarray(np.array(TtNmFs))
        self.ldN = float(sum(ldNs))
        self.solve1d = [N.make_solve_1d() for N in Ns]
        self.Ft = [jnp.asarray(np.asarray(F.T)) for F in Fs_rn]
        self.Tt = [jnp.asarray(np.asarray(T.T)) for T in Ts]
        # vsm.Ns[p] is a WoodburyKernel_novar = white NoiseMatrix1D_novar (.N.N, the per-TOA
        # variance) Woodbury-combined with the linear timing-model GP (.F, prior .P). B1 draws
        # the WHITE part only; the timing GP is a marginalisation device, not a noise process.
        self.N_diag = [_white_diag(N, i) for i, N in enumerate(Ns)]
        self.tm_basis = [np.asarray(getattr(N, "F", np.zeros((0, 0)))) for N in Ns]
        self.tm_prior = [np.asarray(getattr(getattr(N, "P", None), "N", np.zeros(0))) for N in Ns]
        self.Fs_rn = [np.asarray(F) for F in Fs_rn]
        self.Ts_gwb = [np.asarray(T) for T in Ts]

        Pinv, ldP = vsm.P_var.make_inv()(self.frozen)
        if np.asarray(Pinv).ndim != 2:
            raise RuntimeError(f"rednoise Pinv ndim={np.asarray(Pinv).ndim}, expected 2 "
                               f"(npsr, ncomp) diagonal -- the RN draw assumes a diagonal prior.")
        self.Phi_rn_diag = np.asarray(1.0 / Pinv)            # (npsr, ncomp_rn) rednoise prior var
        if Pinv.ndim == 2:
            i1, i2 = jnp.diag_indices(self.FtNmF.shape[1], ndim=2)
            self.cf = dsm.matrix_factor(self.FtNmF.at[:, i1, i2].add(Pinv))
        else:
            self.cf = dsm.matrix_factor(self.FtNmF + Pinv)

        self.smask = None                       # per-source feed weight (FORGE-G); None -> all-on
                                                # FORGE-G2: a runtime arg of the jitted evaluators
        self._ppab = jax.jit(self._per_pulsar_ab_impl)
        self.a_const = self.a_const_from(theta_truth, data_ref, mask_ref)
        self._Minv = None
        self._ab_fns = {}
        self.AI = None

    def a_const_from(self, theta, data_tuple, pmask):
        """a_const by difference: lnL - (sum a_dep + 0.5[b M^-1 b - ldP - logdet]).
        Invariant to (theta, data, pmask) -- it collects ldN + the rednoise logdet. GATED."""
        a_dep, b = self._ppab(jnp.asarray(theta), data_tuple, jnp.asarray(pmask), self.smask)
        rest = float(jnp.sum(a_dep)) + 0.5 * float(
            b.reshape(-1) @ dsm.matrix_solve(self.cf_cached, b.reshape(-1))
            - self.ldP_cached - self.logdet_cached)
        return float(self.amo["logL"](jnp.asarray(theta), data_tuple, jnp.asarray(pmask))) - rest

    def enable_fast_full(self, on=True):
        """Swap the per-fringe residual to the pterm-only (mask-free) callable. ONLY valid when
        every subsequent call uses pmask == 1; the caller must gate this (see B1Marg)."""
        if on == self._fast_full:
            return
        self._fast_full = on
        self.ys = self._ys_full if on else self._ys_masked
        self._ab_fns = {}                       # per-pulsar jitted evaluators are now stale
        self._ppab = jax.jit(self._per_pulsar_ab_impl)

    def set_smask(self, w):
        """Set the per-source feed weight (FORGE-G soft source side). FORGE-G2 (TURBINE): smask
        is a RUNTIME argument of the jitted evaluators, so changing the PATTERN is a plain
        argument change -- no cache clear, no recompile. Only switching None <-> array retraces
        (once per direction: the params pytree gains/loses the SMASK key), and both traces stay
        cached thereafter. w == None restores the incumbent all-on behaviour: the SMASK key is
        absent from params, so the trace is bit-identical to the pre-FORGE-G path.

        The pterm-only fast-full residual (ys_full = MultiSourceDelay) does NOT honour SMASK, so a
        non-None smask forces the masked residual path -- the caller (B1Marg) must not re-enable
        fast-full while a smask is live. Gated by the soft-source module's threshold-inf limit."""
        if w is not None:
            w = jnp.asarray(np.asarray(w, dtype=np.float64))
            if self._fast_full:
                self.enable_fast_full(False)     # MultiSourceDelay ignores SMASK -- force masked path
        self.smask = w

    # ---- core algebra ----
    def _params(self, theta_arr, data_tuple, pmask, smask=None):
        params = dict(self.frozen)
        for k, v in zip(self.theta_keys, theta_arr):
            params[k] = v
        for dk, d in zip(self.data_keys, data_tuple):
            params[dk] = d
        params[PMASK] = pmask
        if smask is not None:
            params[SMASK] = smask
        return params

    def _per_pulsar_ab_impl(self, theta_arr, data_tuple, pmask, smask=None):
        params = self._params(theta_arr, data_tuple, pmask, smask)
        ytNmy = []; FtNmy = []; TtNmy = []
        for p in range(self.npsr):
            yp = self.ys[p](params)
            Nmy, _ = self.solve1d[p](yp)
            ytNmy.append(yp @ Nmy)
            FtNmy.append(self.Ft[p] @ Nmy)
            TtNmy.append(self.Tt[p] @ Nmy)
        ytNmy = jnp.stack(ytNmy); FtNmy = jnp.stack(FtNmy); TtNmy = jnp.stack(TtNmy)
        sol = dsm.matrix_solve(self.cf, FtNmy)
        a_dep = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))
        b = TtNmy - jnp.sum(self.TtNmF * sol[:, jnp.newaxis, :], axis=2)
        return a_dep, b

    def lnL(self, theta_arr, data_tuple, pmask):
        a_dep, b = self._ppab(jnp.asarray(theta_arr), data_tuple, jnp.asarray(pmask), self.smask)
        bflat = b.reshape(-1)
        return float(jnp.sum(a_dep) + self.a_const
                     + 0.5 * (bflat @ dsm.matrix_solve(self.cf_cached, bflat)
                              - self.ldP_cached - self.logdet_cached))

    def _ensure_minv(self):
        if self._Minv is None:
            ng = self.ngp_gwb; M = self.npsr * ng
            self._Minv = np.asarray(dsm.matrix_solve(self.cf_cached, jnp.eye(M)))
            self._Minv_pp = np.stack([self._Minv[p * ng:(p + 1) * ng, p * ng:(p + 1) * ng]
                                      for p in range(self.npsr)])

    def _pulsar_ab_fn(self, p):
        """(theta_base, Lvals, data, pmask, smask) -> (a_dep (B,), b (B,ngp)) for pulsar p.
        Everything is a runtime arg (incl. smask, FORGE-G2), so this compiles ONCE and serves
        every realisation and every feed pattern."""
        if p in self._ab_fns:
            return self._ab_fns[p]
        Ft = self.Ft[p]; Tt = self.Tt[p]; s1d = self.solve1d[p]; ys_p = self.ys[p]
        cf_p = (self.cf[0][p], self.cf[1]); TtNmF_p = self.TtNmF[p]; AI_p = int(self.AI[p])

        @jax.jit
        def run(theta_base, Lvals, data_tuple, pmask, smask):
            def one(L):
                params = self._params(theta_base.at[AI_p].set(L), data_tuple, pmask, smask)
                yp = ys_p(params); Nmy, _ = s1d(yp)
                return yp @ Nmy, Ft @ Nmy, Tt @ Nmy
            ytNmy, FtNmy, TtNmy = jax.vmap(one)(Lvals)
            sol = dsm.matrix_solve(cf_p, FtNmy.T).T
            a_dep = -0.5 * (ytNmy - jnp.sum(FtNmy * sol, axis=1))
            b = TtNmy - sol @ TtNmF_p.T
            return a_dep, b
        self._ab_fns[p] = run
        return run

    def estep(self, theta_base, EV, AI, data_tuple, pmask):
        """LNL (npsr, B): pulsar p's distance swept over EV[p], all others at theta_base."""
        self.AI = np.asarray(AI); self._ensure_minv()
        tb = jnp.asarray(theta_base); pm = jnp.asarray(pmask)
        a_base, b_base = self._ppab(tb, data_tuple, pm, self.smask)
        a_base = np.asarray(a_base); b_base = np.asarray(b_base)
        bflat = b_base.reshape(-1)
        u = self._Minv @ bflat
        Q_base = float(bflat @ u)
        u_p = u.reshape(self.npsr, self.ngp_gwb)
        sum_a = float(a_base.sum())
        const = self.a_const - 0.5 * (self.ldP_cached + self.logdet_cached)
        npsr, B = EV.shape
        LNL = np.empty((npsr, B))
        for p in range(npsr):
            a_pf, b_pf = self._pulsar_ab_fn(p)(tb, jnp.asarray(EV[p]), data_tuple, pm, self.smask)
            a_pf = np.asarray(a_pf); b_pf = np.asarray(b_pf)
            db = b_pf - b_base[p]
            Q_pf = Q_base + 2.0 * (db @ u_p[p]) + np.einsum('bi,ij,bj->b', db, self._Minv_pp[p], db)
            LNL[p] = const + (sum_a - a_base[p] + a_pf) + 0.5 * Q_pf
        return LNL

    def lnref(self, theta_base, data_tuple, pmask):
        """Exact lnL at theta_base (assembled from a/b, i.e. the split's own value)."""
        return self.lnL(theta_base, data_tuple, pmask)


# ============================================================
# 4. Real noise draws (white + per-pulsar RN GP + HD-correlated GWB GP)
# ============================================================

# The GWB square root is BANKED, not recomputed. See load_or_build_L_gwb.
L_GWB_BANK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "b1_L_gwb.npz")


def lgwb_fingerprint(L):
    """Content hash of a GWB square root. Two runs agree iff they draw the same noise."""
    return hashlib.sha256(np.ascontiguousarray(L, dtype=np.float64).tobytes()).hexdigest()[:16]


def load_or_build_L_gwb(Phi_gwb, path=L_GWB_BANK, strict=False):
    """Return (L_gwb, provenance) with Phi_gwb = L L^T.

    THE THREAD-COUNT HAZARD (CHORUS 2026-07-13, §0.1; found by its g1 gate).
    Phi_gwb is the HD-correlated GWB prior covariance and its spectrum is NEAR-DEGENERATE.
    `np.linalg.eigh` fixes the eigenvector basis INSIDE a degenerate subspace arbitrarily,
    and LAPACK's choice depends on the BLAS THREAD COUNT. Two runs at the same seed but
    different `--cpus-per-task` therefore draw a ROTATED — different, though
    equal-in-distribution — GWB realisation. Every banked noisy realisation in this repo
    (FORGE, IGNITE, IGNITE-2, D4, SURFACE, CHORUS, ANCHOR) was drawn at the convention
    `--cpus-per-task=8`, and reproduces bit-identically ONLY there. Until this fix,
    *cpus-per-task was part of the seed* — an undeclared input to every result.

    THE FIX: bank the square root itself. `L_gwb` is a frozen artifact, loaded from disk;
    no BLAS call sits between the seed and the draw, so `draw(seed)` is bit-identical at
    ANY thread count, on any machine. Gated by `trackB_lgwb_gate.py`.

    SCOPE OF INFERENCE. The banked artifact must be the one the existing banks were drawn
    under, i.e. generated at `--cpus-per-task=8` ON ACCRE and validated against ANCHOR's g1
    replay (80 banked `ig_nullN_*` realisations, bit-identical). A file generated anywhere
    else reproduces itself forever but NOT the bank. This box (cronus, 24 cores) cannot
    produce the canonical basis; the fallback path below is therefore a WARNING, not a
    silent recompute.
    """
    if os.path.exists(path):
        z = np.load(path)
        L = np.asarray(z["L_gwb"], dtype=np.float64)
        fp = lgwb_fingerprint(L)
        banked_fp = str(z["fingerprint"]) if "fingerprint" in z.files else fp
        if fp != banked_fp:
            raise RuntimeError(
                f"L_gwb bank {path} is CORRUPT: fingerprint {fp} != recorded {banked_fp}")
        if L.shape != Phi_gwb.shape:
            raise RuntimeError(
                f"L_gwb bank {path} has shape {L.shape}, problem needs {Phi_gwb.shape} "
                "— the bank belongs to a different array/GP configuration.")
        prov = f"BANKED {os.path.basename(path)} fp={fp} cpus={z['cpus'] if 'cpus' in z.files else '?'}"
        return L, prov

    # ---- fallback: no bank on disk. Recompute, and say so loudly. ----
    nthreads = os.environ.get("OMP_NUM_THREADS") or os.environ.get("OPENBLAS_NUM_THREADS") or "unset"
    w, V = np.linalg.eigh(Phi_gwb)
    w = np.clip(w, 0.0, None)
    L = V * np.sqrt(w)
    fp = lgwb_fingerprint(L)
    msg = (f"[NoiseDrawer] WARNING: no banked L_gwb at {path}. Recomputed by np.linalg.eigh "
           f"at BLAS threads={nthreads} -> fingerprint {fp}. THIS DRAW IS THREAD-COUNT "
           f"DEPENDENT and is NOT guaranteed to reproduce any banked realisation. "
           f"Bank L_gwb (see load_or_build_L_gwb.__doc__) before quoting any number from it.")
    if strict:
        raise RuntimeError(msg)
    print(msg, file=sys.stderr)
    return L, f"RECOMPUTED-UNSAFE threads={nthreads} fp={fp}"


class NoiseDrawer:
    """Draws the three stochastic components the likelihood's covariance actually contains.

    white_p ~ N(0, diag(N_p))                        N_p = efac^2 (toaerr^2 + equad^2)
    rn_p    ~ F_rn,p @ N(0, diag(Phi_rn,p))          30 comp, frozen per-pulsar powerlaw
    gwb     ~ [T_p] @ N(0, Phi_gwb)                  14 comp, HD-correlated across pulsars

    The linear timing-model GP is deliberately NOT drawn: it is a marginalisation device
    (prior variance 1e-14 on the normalised design matrix), not a physical noise process.
    """

    def __init__(self, sp, amo, lgwb_path=L_GWB_BANK, strict=False):
        self.sp = sp
        it = amo["internals"]
        self.npsr = sp.npsr
        self.ngp = sp.ngp_gwb
        Pinv_gwb = np.asarray(it["Pinv_gwb"])          # (npsr*ngp, npsr*ngp)
        Phi_gwb = np.linalg.inv(0.5 * (Pinv_gwb + Pinv_gwb.T))
        Phi_gwb = 0.5 * (Phi_gwb + Phi_gwb.T)
        self.L_gwb, self.lgwb_prov = load_or_build_L_gwb(Phi_gwb, lgwb_path, strict=strict)
        self.Phi_gwb = Phi_gwb
        self.sig_white = [np.sqrt(nd) for nd in sp.N_diag]
        self.sig_rn = np.sqrt(np.clip(sp.Phi_rn_diag, 0.0, None))   # (npsr, ncomp_rn)

    def draw(self, seed, components=("white", "rn", "gwb")):
        rng = np.random.default_rng(seed)
        M = self.npsr * self.ngp
        out = []
        c_gwb = self.L_gwb @ rng.standard_normal(M) if "gwb" in components else np.zeros(M)
        for p in range(self.npsr):
            n = len(self.sig_white[p])
            r = np.zeros(n)
            if "white" in components:
                r = r + self.sig_white[p] * rng.standard_normal(n)
            if "rn" in components:
                r = r + self.sp.Fs_rn[p] @ (self.sig_rn[p] * rng.standard_normal(self.sig_rn.shape[1]))
            if "gwb" in components:
                r = r + self.sp.Ts_gwb[p] @ c_gwb[p * self.ngp:(p + 1) * self.ngp]
            out.append(jnp.asarray(r))
        return tuple(out)


# ============================================================
# 5. Problem construction (pop draw, targeted)
# ============================================================
def build_b1_problem(verbose=True):
    t0 = time.time()
    ent_psrs, disco_psrs = load_pulsars(116)
    names = [p.name for p in ent_psrs]
    priors = A2.load_real_priors(names)
    ncw = POP["ncw"]
    pta, cwb, _ = build_enterprise_pta(ent_psrs, ncw, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)
    inj = make_geometry_injection(pta, ent_psrs, ncw, cwb, seed=POP["seed"],
                                  population=POP["population"])
    amo = build_b1_amortised(disco_psrs, ncw, inj, cwb)
    if verbose:
        print(f"  build_b1_amortised: {time.time()-t0:.1f}s", flush=True)

    temp = _enterprise_to_disco_params(inj, cwb)
    theta_truth = np.array([temp[k] for k in amo["theta_keys"]], dtype=np.float64)
    npsr = len(ent_psrs)
    L0 = np.array([pe.pdist[0] for pe in ent_psrs])
    idx_of = {k: i for i, k in enumerate(amo["theta_keys"])}
    AI = np.array([idx_of[f"{pe.name}_cw_p_dist"] for pe in ent_psrs])

    cw_list = [{k: inj[f"{nm}_{k}"] for k in CW_PARAM_BASE} for nm in cwb]
    dL_truth = np.array([min(compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"],
                                                  cw["log10_fgw"], disco_psrs[a].pos)
                             for cw in cw_list) for a in range(npsr)])
    ones = jnp.ones(npsr); zeros = jnp.zeros(npsr)

    # compile warmers
    data_zero = amo["inject_delay"](jnp.asarray(theta_truth), ones)
    _ = amo["logL"](jnp.asarray(theta_truth), data_zero, ones)
    _ = amo["logL"](jnp.asarray(theta_truth), data_zero, zeros)

    sp = B1Split(amo, theta_truth, data_zero, np.ones(npsr))
    nd = NoiseDrawer(sp, amo)

    P = dict(ent_psrs=ent_psrs, disco_psrs=disco_psrs, names=names, priors=priors,
             pta=pta, cwb=cwb, inj=inj, ncw=ncw, npsr=npsr, amo=amo, sp=sp, nd=nd,
             theta_truth=theta_truth, L0=L0, AI=AI, dL_truth=dL_truth,
             n_dist=amo["n_dist"], n_theta=amo["n_theta"],
             mask_one=np.ones(npsr), mask_zero=np.zeros(npsr),
             census_idx=np.array([names.index(n) for n in CENSUS]),
             anchor_idx=names.index(ANCHOR))
    P["EV_truth"] = build_EV(P, dL_truth)
    if verbose:
        print(f"  build_b1_problem total {time.time()-t0:.1f}s "
              f"(npsr={npsr}, ncw={ncw}, n_theta={P['n_theta']})", flush=True)
    return P


def build_EV(P, dL):
    """Per-pulsar B_CERT-point evaluation grid over +-3 sigma_EM about L0 (A2.eval_grid)."""
    EV = np.empty((P["npsr"], B_CERT))
    pri = P["priors"]
    for a, pe in enumerate(P["ent_psrs"]):
        R = 3.0 * max(pri["feather"][a], pri["script"][a], pri["lit"][a])
        EV[a] = A2.eval_grid(pe.pdist[0], dL[a], R)
    return EV


# ============================================================
# 6. Fringe posterior (vectorised, prior-column = literature)
# ============================================================
class FringeTables:
    """Per-pulsar unique-fringe reduction + Gaussian prior tail, for a given (EV, dL)."""

    def __init__(self, P, EV, dL, prior_key="lit"):
        self.npsr = P["npsr"]
        self.L0 = P["L0"]; self.dL = dL; self.EV = EV
        sig = P["priors"][prior_key]
        self.sig = sig
        self.order = []; self.redidx = []; self.uk = []; self.lnprior = []; self.K_counted = []
        for a in range(self.npsr):
            k = np.round((EV[a] - self.L0[a]) / dL[a]).astype(int)
            o = np.argsort(k, kind="stable"); self.order.append(o)
            uk, i0 = np.unique(k[o], return_index=True)
            self.uk.append(uk); self.redidx.append(i0)
            offs = uk * dL[a]
            self.lnprior.append(-offs ** 2 / (2.0 * sig[a] ** 2))
            self.K_counted.append(2 * int(np.floor(3.0 * sig[a] / dL[a] + 1e-9)))
        self.K_counted = np.array(self.K_counted)

    def posterior(self, LNL):
        """LNL (npsr,B) -> dict of per-pulsar arrays.

        q_max     : posterior mass on the MAP fringe (the pipeline's TRUTH-BLIND confidence)
        map_k     : MAP fringe index
        lik_map_k : likelihood-only (no prior) best fringe index
        peaks     : list of per-fringe peak lnL
        """
        npsr = self.npsr
        out = dict(q_max=np.zeros(npsr), map_k=np.zeros(npsr, int),
                   lik_map_k=np.zeros(npsr, int), map_L=np.zeros(npsr),
                   entropy=np.zeros(npsr))
        self.q = [None] * npsr
        self.peak = [None] * npsr
        for a in range(npsr):
            pk = np.maximum.reduceat(LNL[a][self.order[a]], self.redidx[a])
            logw = (pk - pk.max()) + self.lnprior[a]
            w = np.exp(logw - logw.max()); w /= w.sum()
            j = int(np.argmax(w))
            out["q_max"][a] = w[j]
            out["map_k"][a] = self.uk[a][j]
            out["map_L"][a] = self.L0[a] + self.uk[a][j] * self.dL[a]
            out["lik_map_k"][a] = self.uk[a][int(np.argmax(pk))]
            ww = w[w > 0]
            out["entropy"][a] = float(-np.sum(ww * np.log(ww)))
            self.q[a] = w; self.peak[a] = pk
        return out

    def p_of_fringe(self, n_true):
        """q_p(n_true) per pulsar (needs posterior() first). NaN if n_true not sampled."""
        out = np.full(self.npsr, np.nan)
        rank_post = np.full(self.npsr, -1, int)
        rank_lik = np.full(self.npsr, -1, int)
        for a in range(self.npsr):
            hit = np.where(self.uk[a] == n_true[a])[0]
            if len(hit) == 0:
                continue
            j = hit[0]
            out[a] = self.q[a][j]
            rank_post[a] = int(np.sum(self.q[a] > self.q[a][j]))
            rank_lik[a] = int(np.sum(self.peak[a] > self.peak[a][j]))
        return out, rank_post, rank_lik


# ============================================================
# 7. Arm-B truth draw
# ============================================================
def draw_true_distances(P, seed):
    """n_true ~ prior-weighted over SAMPLED fringe centres; within-fringe offset ~ U(-1/2,1/2).

    Returns (L_true (npsr,), n_true (npsr,), u_true (npsr,)).
    Drawing n_true only over fringes the B_CERT grid actually samples guarantees the truth is
    inside the search window (otherwise 'wrong certification' would be trivially guaranteed for
    the wide-prior pulsars whose window holds K >> B_CERT fringes).
    """
    rng = np.random.default_rng(seed)
    EV = P["EV_truth"]; dL = P["dL_truth"]; L0 = P["L0"]; sig = P["priors"]["lit"]
    n_true = np.zeros(P["npsr"], int); u_true = np.zeros(P["npsr"])
    for a in range(P["npsr"]):
        k = np.unique(np.round((EV[a] - L0[a]) / dL[a]).astype(int))
        lw = -(k * dL[a]) ** 2 / (2.0 * sig[a] ** 2)
        w = np.exp(lw - lw.max()); w /= w.sum()
        n_true[a] = int(rng.choice(k, p=w))
        u_true[a] = rng.uniform(-0.5, 0.5)
    L_true = L0 + (n_true + u_true) * dL
    bad = L_true <= 1e-3
    if bad.any():                       # keep distances physical; re-pin those to the fringe
        L_true[bad] = L0[bad] + n_true[bad] * dL[bad]
        u_true[bad] = 0.0
    return L_true, n_true, u_true


def n_true_on_grid(L_true, L0, dL):
    """The fringe index of the grid (anchored at L0, spacing dL) that CONTAINS L_true."""
    return np.round((L_true - L0) / dL).astype(int)


# ============================================================
# 8. GATES
# ============================================================
def _finite(name, *arrs):
    for a in arrs:
        a = np.asarray(a)
        if not np.all(np.isfinite(a)):
            raise RuntimeError(f"NON-FINITE VALUES in {name}: "
                               f"{int(np.sum(~np.isfinite(a)))}/{a.size}")


def mode_gate():
    print("=== B1 CORE GATES ===", flush=True)
    ok = True
    P = build_b1_problem()
    amo = P["amo"]; sp = P["sp"]; tt = P["theta_truth"]
    npsr = P["npsr"]
    one = jnp.ones(npsr); zero = jnp.zeros(npsr)
    data0 = amo["inject_delay"](jnp.asarray(tt), one)

    # ---- G0: mask=1 == stagec_fisher's full likelihood ----
    print("\n-- G0: masked lnL(mask=1) == stagec_fisher.fisher_logL (full pterm) --", flush=True)
    from stagec_fisher import build_fisher_amortised
    t0 = time.time()
    amo_ref = build_fisher_amortised(P["disco_psrs"], P["ncw"], P["inj"], P["cwb"])
    print(f"   reference build {time.time()-t0:.0f}s", flush=True)
    assert list(amo_ref["theta_keys"]) == list(amo["theta_keys"]), "theta layout mismatch"
    data_ref = amo_ref["inject_data"](jnp.asarray(tt))
    dd = max(float(np.max(np.abs(np.asarray(a) - np.asarray(b)))) for a, b in zip(data0, data_ref))
    print(f"   inject_delay(mask=1) vs inject_data: max|diff| = {dd:.3e} s", flush=True)
    rng = np.random.default_rng(0)
    scale = TE.phi_scale(P); plo, phi_hi = TE.phi_bounds(P)
    thetas = [tt.copy()]
    for _ in range(3):
        th = tt.copy()
        th[P["n_dist"]:] = np.clip(th[P["n_dist"]:] + rng.normal(size=len(tt) - P["n_dist"]) * scale * 0.05,
                                   plo, phi_hi)
        thetas.append(th)
    w0 = 0.0
    for th in thetas:
        a = float(amo["logL"](jnp.asarray(th), data_ref, one))
        b = float(amo_ref["fisher_logL"](jnp.asarray(th), data_ref))
        w0 = max(w0, abs(a - b))
    g0 = w0 < 1e-8; ok &= g0
    print(f"   [G0] max|diff| = {w0:.3e} -> {'PASS' if g0 else 'FAIL'} (<1e-8)", flush=True)

    # ---- G1: mask=0 == EarthDelay likelihood ----
    print("\n-- G1: masked lnL(mask=0) == EarthDelay (pterm off) --", flush=True)
    t0 = time.time()
    amo_e = build_fisher_amortised(P["disco_psrs"], P["ncw"], P["inj"], P["cwb"],
                                   msd_factory=TE.EarthDelay)
    print(f"   earth reference build {time.time()-t0:.0f}s", flush=True)
    w1 = 0.0
    for th in thetas:
        a = float(amo["logL"](jnp.asarray(th), data_ref, zero))
        b = float(amo_e["fisher_logL"](jnp.asarray(th), data_ref))
        w1 = max(w1, abs(a - b))
    g1 = w1 < 1e-8; ok &= g1
    print(f"   [G1] max|diff| = {w1:.3e} -> {'PASS' if g1 else 'FAIL'} (<1e-8)", flush=True)

    # partial mask sanity: mask with only census on lies strictly between
    m_part = np.zeros(npsr); m_part[P["census_idx"]] = 1.0
    lp = float(amo["logL"](jnp.asarray(tt), data_ref, jnp.asarray(m_part)))
    l1 = float(amo["logL"](jnp.asarray(tt), data_ref, one))
    l0 = float(amo["logL"](jnp.asarray(tt), data_ref, zero))
    print(f"   partial mask (census-3 on): lnL0={l0:.2f} < lnL_part={lp:.2f} < lnL1={l1:.2f} "
          f"-> {'OK' if l0 < lp < l1 else 'CHECK'}", flush=True)

    # ---- G2: a_const is data- and mask-independent ----
    print("\n-- G2: split a_const invariant to (theta, data, mask) --", flush=True)
    nz = P["nd"].draw(seed=101)
    data_n = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(data0, nz))
    ac_alt = sp.a_const_from(thetas[1], data_n, m_part)
    g2 = abs(sp.a_const - ac_alt) < 1e-6 * max(1.0, abs(sp.a_const)); ok &= g2
    print(f"   a_const(truth, zero-noise, mask=1)     = {sp.a_const:.6f}", flush=True)
    print(f"   a_const(perturbed, noisy, mask=census) = {ac_alt:.6f}  "
          f"|diff|={abs(sp.a_const-ac_alt):.3e} -> {'PASS' if g2 else 'FAIL'}", flush=True)

    # ---- G3: split.lnL == amo.logL over (data, mask) combos ----
    print("\n-- G3: B1Split.lnL == masked logL (data + mask runtime) --", flush=True)
    w3 = 0.0
    for dat, lab in [(data_ref, "zero-noise"), (data_n, "noisy")]:
        for m, ml in [(one, "mask=1"), (zero, "mask=0"), (jnp.asarray(m_part), "mask=census")]:
            for th in thetas[:2]:
                a = float(amo["logL"](jnp.asarray(th), dat, m))
                b = sp.lnL(th, dat, m)
                w3 = max(w3, abs(a - b))
    g3 = w3 < 1e-8; ok &= g3
    print(f"   [G3] max|diff| over {{data}}x{{mask}}x{{theta}} = {w3:.3e} -> "
          f"{'PASS' if g3 else 'FAIL'} (<1e-8)", flush=True)

    # ---- G4: E-step == banked Split.estep (zero-noise, mask=1) ----
    print("\n-- G4: B1Split.estep == trackB_split.Split.estep (zero-noise, mask=1) --", flush=True)
    import trackB_split as TS
    sp_ref = TS.Split(amo_ref, data_ref, tt)
    EV = P["EV_truth"]; AI = P["AI"]
    L_ref = sp_ref.estep(tt, EV, AI)
    t0 = time.time(); L_b1 = sp.estep(tt, EV, AI, data_ref, one); dt = time.time() - t0
    _finite("estep", L_b1)
    w4 = float(np.max(np.abs(L_ref - L_b1)))
    nbad = int(np.sum(np.max(np.abs(L_ref - L_b1), axis=1) > 1e-8))
    g4 = w4 < 1e-8; ok &= g4
    print(f"   [G4] max|dlnL| = {w4:.3e} ({nbad}/{npsr} pulsars >1e-8) -> "
          f"{'PASS' if g4 else 'FAIL'}; warm estep {dt:.2f}s", flush=True)

    # ---- G5: E-step at truth reproduces the A2 census ceilings (zero noise) ----
    print("\n-- G5: zero-noise E-step at truth == A2 census ceilings --", flush=True)
    FT = FringeTables(P, EV, P["dL_truth"])
    post = FT.posterior(L_b1)
    pt, rp, rl = FT.p_of_fringe(np.zeros(npsr, int))
    g5 = True
    for nm in CENSUS:
        a = P["names"].index(nm)
        d = abs(pt[a] - TE.CEILING[nm]); g5 &= d < 0.02
        print(f"   {nm:14s} P(true)={pt[a]:.4f}  ceiling={TE.CEILING[nm]:.3f}  |diff|={d:.4f}", flush=True)
    ncert = int(np.sum(post["q_max"] > 0.9))
    ncert_true = int(np.sum((post["q_max"] > 0.9) & (post["map_k"] == 0)))
    print(f"   truth-blind certified (q_max>0.9): {ncert}   of which on the TRUE fringe: {ncert_true}", flush=True)
    ok &= g5
    print(f"   [G5] -> {'PASS' if g5 else 'FAIL'}", flush=True)

    # ---- G6: noise draw statistics ----
    print("\n-- G6: noise draw covariance sanity --", flush=True)
    nd = P["nd"]
    ndraw = 24
    a0 = P["names"].index(CENSUS[0]); a1 = P["names"].index(CENSUS[1])
    W = np.array([np.asarray(nd.draw(seed=1000 + i, components=("white",))[a0]) for i in range(ndraw)])
    emp = W.var(axis=0, ddof=1); exp_ = np.asarray(sp.N_diag[a0])
    r_white = float(np.median(emp / exp_))
    # GWB: HD cross-correlation between two pulsars over many draws
    ndg = 400
    G0 = np.zeros(ndg); G1 = np.zeros(ndg)
    ng = sp.ngp_gwb
    rngc = np.random.default_rng(7)
    C = nd.Phi_gwb
    hd_pred = C[a0 * ng, a1 * ng] / np.sqrt(C[a0 * ng, a0 * ng] * C[a1 * ng, a1 * ng])
    Lg = nd.L_gwb
    for i in range(ndg):
        c = Lg @ rngc.standard_normal(npsr * ng)
        G0[i] = c[a0 * ng]; G1[i] = c[a1 * ng]
    emp_corr = float(np.corrcoef(G0, G1)[0, 1])
    hd_true = float(ds.hd_orf(P["disco_psrs"][a0].pos, P["disco_psrs"][a1].pos))
    hd_auto0 = float(ds.hd_orf(P["disco_psrs"][a0].pos, P["disco_psrs"][a0].pos))
    print(f"   white: median(empirical var / N_diag) over {ndraw} draws = {r_white:.3f} "
          f"(expect ~1, {ndraw}-draw scatter ~{np.sqrt(2/(ndraw-1)):.2f})", flush=True)
    print(f"   gwb  : corr({CENSUS[0]},{CENSUS[1]}) empirical {emp_corr:+.4f} vs Phi-implied "
          f"{hd_pred:+.4f}; hd_orf(cross)={hd_true:+.4f} hd_orf(auto)={hd_auto0:+.4f} "
          f"-> normalised ORF {hd_true/hd_auto0:+.4f}", flush=True)
    g6 = (abs(r_white - 1.0) < 0.35) and (abs(emp_corr - hd_pred) < 0.15)
    ok &= g6
    print(f"   [G6] -> {'PASS' if g6 else 'FAIL'}", flush=True)

    # ---- G7: noisy data is finite; lnL finite; noise RMS vs signal RMS ----
    print("\n-- G7: noisy-data finiteness + signal/noise scale --", flush=True)
    nz_full = nd.draw(seed=2001)
    _finite("noise draw", *[np.asarray(x) for x in nz_full])
    sig_rms = np.median([float(np.std(np.asarray(d))) for d in data0])
    noi_rms = np.median([float(np.std(np.asarray(n))) for n in nz_full])
    dn = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(data0, nz_full))
    ln = float(amo["logL"](jnp.asarray(tt), dn, one))
    lz = float(amo["logL"](jnp.asarray(tt), data_ref, one))
    g7 = np.isfinite(ln)
    ok &= g7
    w_rms = np.median([float(np.sqrt(np.mean(v))) for v in sp.N_diag])
    r_rms = np.median([float(np.std(sp.Fs_rn[p] @ (nd.sig_rn[p] * np.random.default_rng(p).standard_normal(nd.sig_rn.shape[1]))))
                       for p in range(npsr)])
    tm_rms = np.median([float(np.sqrt(np.mean(np.sum(sp.tm_basis[p] ** 2 * sp.tm_prior[p][None, :], axis=1))))
                        for p in range(npsr) if sp.tm_prior[p].size])
    print(f"   median per-pulsar CW rms = {sig_rms*1e9:.1f} ns; total noise rms = {noi_rms*1e9:.1f} ns", flush=True)
    print(f"   components: white {w_rms*1e9:.1f} ns, rednoise {r_rms*1e9:.1f} ns "
          f"(gwb = remainder); timing-GP prior rms {tm_rms*1e9:.3g} ns (NOT drawn: "
          f"marginalisation device)", flush=True)
    print(f"   lnL(truth | zero-noise) = {lz:.2f}; lnL(truth | noisy) = {ln:.2f}  "
          f"-> {'PASS' if g7 else 'FAIL'}", flush=True)

    # ---- G8: arm-B truth draw is inside the search window ----
    print("\n-- G8: arm-B distance truth draw --", flush=True)
    Lt, nt, ut = draw_true_distances(P, seed=5000)
    _finite("L_true", Lt)
    inwin = np.array([nt[a] in FT.uk[a] for a in range(npsr)])
    off_pc = (Lt - P["L0"]) * KPC_TO_PC
    g8 = bool(inwin.all()) and bool((Lt > 0).all())
    ok &= g8
    print(f"   |n_true| median {int(np.median(np.abs(nt)))}, max {int(np.abs(nt).max())}; "
          f"sampled-in-window {int(inwin.sum())}/{npsr}", flush=True)
    print(f"   |L_true - L0| median {np.median(np.abs(off_pc)):.2f} pc, max {np.abs(off_pc).max():.1f} pc; "
          f"n_true==0 for {int(np.sum(nt==0))} pulsars", flush=True)
    print(f"   census n_true: " + " ".join(f"{nm.split('-')[0].split('+')[0]}={nt[P['names'].index(nm)]}"
                                           for nm in CENSUS)
          + f";  {ANCHOR}={nt[P['anchor_idx']]}", flush=True)
    print(f"   [G8] -> {'PASS' if g8 else 'FAIL'}", flush=True)

    print(f"\n=== B1 CORE GATES: {'ALL PASS' if ok else 'FAIL'} ===", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    sys.exit(mode_gate())
