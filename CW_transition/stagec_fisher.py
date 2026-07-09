"""Stage C — Deliverable D2: amortised zero-noise Fisher + conditional/marginal sigma_L.

Amortisation (D2 carry-forward note 2): the fast-scan residual is a callable
`y(params) = data - delay(params)`; only `data` is a baked constant. Here we make
`data` a RUNTIME argument, so ONE compiled Hessian graph serves every CW-geometry
draw in the ensemble — pay the ~350-500 s XLA compile once, then ~1-3 s warm per draw.

Zero-noise Fisher (D2 carry-forward note 1): evaluate on data = the injected CW at
truth (no noise draw), but with the real heterogeneous noise covariance (white +
frozen GWB GP + frozen RN GP) still weighting the likelihood. Then residual(truth)=0
exactly, so H = -Fisher, negative-semidefinite by construction. Verified via the
negative-curvature count (expect N_psr/N_psr).

sigma_L (D2 body):
  conditional_i = 1/sqrt(F_dd[i,i])                      (CW params known)
  marginal_i    = sqrt( pinv(F_marg)[i,i] ),  F_marg = F_dd - F_dc F_cc^{-1} F_cd
                                                          (CW params marginalised;
                                                           Schur complement)
Both converted to PARSEC and compared to each pulsar's EM prior (psr.pdist[1]).
F_cc^{-1} (and F_marg pinv) use the toy's eigh pseudo-inverse discipline (Stage A.1):
rcond in {1e-8, 1e-10, 1e-12}, require sigma_L stable across them.

Run in jug-gpu:
  python stagec_fisher.py gate        # amortised H == D1 dense H on 5-psr zero-noise
  python stagec_fisher.py d2 --ndraws 10
"""

import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse
import sys
import time

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
import discovery as ds
from discovery import matrix as dsm
from cw_helpers import (
    load_pulsars, build_enterprise_pta, generate_injection_params,
    MultiSourceDelay, _enterprise_to_disco_params, compute_mode_spacing,
)

# nb-05 realistic operating point (matches stagec_hessian.py)
SEED = 1234
N_COMPONENTS = 14
LOG10_EQUAD = -6.0
GWB_LOG10_A = -14.6
GWB_GAMMA = 13.0 / 3.0
LOG10_H_RANGE = (-14.0, -13.5)
LOG10_FGW_RANGE = (-8.0, -7.5)
LOG10_MC_RANGE = (8.5, 9.5)
KPC_TO_PC = 1000.0

CW_PARAM_BASE = [
    "cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
    "log10_fgw", "log10_h", "phase0", "psi",
]


def cw_param_keys(n_cw):
    suffixes = ["" if i == 0 else f"_{i+1}" for i in range(n_cw)]
    return [f"cw_{base}{suffix}" for suffix in suffixes for base in CW_PARAM_BASE]


def device_peak_gib():
    st = jax.local_devices()[0].memory_stats() or {}
    return st.get("peak_bytes_in_use", -1) / 2**30


# ============================================================
# Amortised zero-noise Fisher likelihood
# ============================================================

def build_fisher_amortised(disco_psrs, num_cw, enterprise_params, cw_block_names,
                           components=N_COMPONENTS, log10_equad=LOG10_EQUAD,
                           rn_components=30, include_pterm=True, msd_factory=None):
    """Replicate build_fast_scan_likelihood but with the residual DATA as a runtime arg.

    Returns dict with:
      fisher_hessian(theta_arr, data_tuple) -> (n_theta, n_theta) Hessian (=-Fisher
          at truth on zero-noise data), JIT-compiled, one shape for all draws.
      theta_keys : ordered [dist keys + cw keys]
      n_dist     : number of distance params
      inject_data(theta_arr) -> data_tuple : per-pulsar CW delay at theta (the
          zero-noise "data"); residual(theta)=0 when data=inject_data(theta).
    """
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

    # Per-pulsar CW delay models (same objects the fast-scan uses). A msd_factory
    # (p, num_cw) -> delay may be supplied to swap the delay model, e.g. Track B's
    # Earth-term-only (fully-decohered) delay for the coherence anneal; its dist keys
    # stay in theta_keys but are inert (the delay ignores p_dist).
    if msd_factory is None:
        msds = [MultiSourceDelay(p, num_cw, include_pterm=include_pterm) for p in disco_psrs]
    else:
        msds = [msd_factory(p, num_cw) for p in disco_psrs]

    # Residual callables with the DATA read from params (key f"{name}__data")
    data_keys = [f"{p.name}__data" for p in disco_psrs]
    ys = []
    for p, msd, dkey in zip(disco_psrs, msds, data_keys):
        def y(params, _msd=msd, _dkey=dkey):
            return params[_dkey] - _msd(params)
        y.params = sorted(set(msd.params) | {dkey})
        ys.append(y)

    pulsar_likes = [
        ds.PulsarLikelihood([
            np.zeros(len(p.toas)),  # placeholder residual; real data via ys callable below
            noise_terms[p.name],
            timing_terms[p.name],
            msd,
        ])
        for p, msd in zip(disco_psrs, msds)
    ]
    fml = ds.ArrayLikelihood(pulsar_likes, commongp=common_gp, globalgp=global_gp)
    _ = fml.logL  # sets up vsm

    P_var_inv = getattr(global_gp, "Phi_inv", None) or global_gp.Phi.make_inv()
    Fs = getattr(global_gp, "Fs", None)
    if Fs is None:
        F = getattr(global_gp, "F", None)
        Fs = F if isinstance(F, (list, tuple)) else [F] * len(pulsar_likes)

    # vary-path kernelterms with OUR data-carrying residual callables
    kterms = fml.vsm.make_kernelterms(ys, Fs)

    npsr = len(Fs)
    ngp = Fs[0].shape[1]

    # frozen (noise/GP) truth values needed by kterms + P_var_inv
    truth_full = {k: temp_dict[k] for k in kterms.params if k in temp_dict}
    for k in P_var_inv.params:
        if k in temp_dict:
            truth_full[k] = temp_dict[k]

    # cf_cached (frozen GWB), exactly as build_fast_scan_likelihood
    # c = terms[2] is data-independent; evaluate at truth with zero data placeholders.
    zero_data = {dk: jnp.zeros(len(p.toas)) for dk, p in zip(data_keys, disco_psrs)}
    terms_truth = kterms({**truth_full, **zero_data})
    c_truth = terms_truth[2]
    Pinv_truth, ldP_truth = P_var_inv(truth_full)
    cf_cached = dsm.matrix_factor(Pinv_truth + dsm.jsp.linalg.block_diag(*c_truth))
    ldP_sum_cached = jnp.sum(ldP_truth)
    logdet_cf_cached = dsm.matrix_norm * jnp.sum(jnp.log(jnp.diag(cf_cached[0])))

    # theta = [dist keys + cw keys]
    dist_keys = [f"{p.name}_cw_p_dist" for p in disco_psrs]
    cw_keys = cw_param_keys(num_cw)
    theta_keys = dist_keys + cw_keys
    n_dist = len(dist_keys)
    theta_truth = jnp.array([temp_dict[k] for k in theta_keys], dtype=jnp.float64)

    # frozen params that are NOT theta and NOT data (white + RN + GWB)
    frozen = {k: temp_dict[k] for k in truth_full if k not in theta_keys}

    def fisher_logL(theta_arr, data_tuple):
        params = dict(frozen)
        for k, v in zip(theta_keys, theta_arr):
            params[k] = v
        for dk, d in zip(data_keys, data_tuple):
            params[dk] = d
        terms = kterms(params)
        p0 = jnp.sum(terms[0])
        FtNmy = terms[1].reshape(npsr * ngp)
        return p0 + 0.5 * (
            FtNmy.T @ dsm.matrix_solve(cf_cached, FtNmy)
            - ldP_sum_cached - logdet_cf_cached
        )

    def inject_data(theta_arr):
        """Per-pulsar CW delay at theta = the zero-noise data (residual(theta)=0)."""
        params = dict(frozen)
        for k, v in zip(theta_keys, theta_arr):
            params[k] = v
        return tuple(msd(params) for msd in msds)

    fisher_hessian = jax.jit(jax.hessian(fisher_logL, argnums=0))
    fisher_logL_jit = jax.jit(fisher_logL)

    # HVP-based Hessian assembly (dense jax.hessian OOMs at large N_CW): one
    # unit-tangent column at a time, chunked, so memory ~ one gradient not the
    # full forward-over-reverse stack. Assembles the full (n_theta, n_theta) H.
    n_theta = len(theta_keys)

    def _hvp_single(theta_arr, data_tuple, v):
        g = lambda th: jax.grad(fisher_logL, argnums=0)(th, data_tuple)
        return jax.jvp(g, (theta_arr,), (v,))[1]

    hvp_batched = jax.jit(jax.vmap(_hvp_single, in_axes=(None, None, 0)))

    # batched logL over a stack of theta (for fringe scans): data fixed
    fisher_logL_batched = jax.jit(jax.vmap(fisher_logL, in_axes=(0, None)))

    return dict(
        fisher_hessian=fisher_hessian, fisher_logL=fisher_logL_jit,
        fisher_logL_batched=fisher_logL_batched,
        hvp_batched=hvp_batched, n_theta=n_theta,
        theta_keys=theta_keys, n_dist=n_dist, theta_truth=theta_truth,
        inject_data=jax.jit(inject_data), disco_psrs=disco_psrs,
        # frozen internals for the low-rank QuickCW split (trackB_split.py). Additive:
        # existing consumers ignore this key. All are distance-independent EXCEPT ys/msds
        # (the residual/delay callables, which carry the L_p dependence the split exploits).
        internals=dict(
            fml=fml, vsm=fml.vsm, kterms=kterms, P_var_inv=P_var_inv,
            cf_cached=cf_cached, ldP_sum_cached=ldP_sum_cached,
            logdet_cf_cached=logdet_cf_cached, c_truth=c_truth,
            Fs_gwb=Fs, ys=ys, msds=msds, frozen=frozen, data_keys=data_keys,
            npsr=npsr, ngp_gwb=ngp,
        ),
    )


# ============================================================
# Geometry-draw problem builder (fixed noise, varying CW geometry)
# ============================================================

def make_geometry_injection(pta, ent_psrs, num_cw, cw_block_names, seed,
                            h_range=None, population=None):
    """Realistic per-source CW geometry draw at FIXED noise (frozen GWB + fixed RN).

    h_range : (lo,hi) tuple overriding LOG10_H_RANGE. Pass (h,h) for EQUAL STRAIN.
    population : (k_loud, log10_h_loud, log10_h_faint) or None. If set, the first
        k_loud sources get log10_h_loud, the rest log10_h_faint (overrides h_range).
    """
    rng = np.random.default_rng(seed)
    np.random.seed(seed)
    inj = generate_injection_params(
        pta, ent_psrs, num_cw, cw_block_names, log10_h=None,
        scenario="realistic", rng=rng,
        log10_h_range=(h_range if h_range is not None else LOG10_H_RANGE),
        log10_fgw_range=LOG10_FGW_RANGE, log10_mc_range=LOG10_MC_RANGE,
        gwb_log10_A=GWB_LOG10_A, gwb_gamma=GWB_GAMMA,
    )
    if population is not None:
        k_loud, h_loud, h_faint = population
        for i, name in enumerate(cw_block_names):
            inj[f"{name}_log10_h"] = h_loud if i < k_loud else h_faint
    # Freeze red noise across draws (so noise covariance is identical): overwrite
    # the per-draw RN draws with a fixed reference (seed-0 draw). Keeps one compiled
    # shape AND a fixed noise floor; only CW geometry varies.
    ref_rng = np.random.default_rng(0)
    for k in list(inj.keys()):
        if "rednoise_log10_A" in k:
            inj[k] = ref_rng.uniform(-15, -13)
        elif "rednoise_gamma" in k:
            inj[k] = ref_rng.uniform(1.5, 5.0)
    return inj


# ============================================================
# sigma_L extraction (conditional + marginal, eigh pinv discipline)
# ============================================================

def eigh_pinv(A, rcond):
    """Symmetric pseudo-inverse via eigh, dropping eigenvalues < rcond*max|eig|."""
    w, V = np.linalg.eigh(0.5 * (A + A.T))
    wmax = np.max(np.abs(w))
    keep = np.abs(w) > rcond * wmax
    winv = np.where(keep, 1.0 / np.where(keep, w, 1.0), 0.0)
    return (V * winv) @ V.T


def sigma_L_from_hessian(H, n_dist, rcond=1e-10):
    """conditional + marginal per-pulsar sigma_L (distance param unit = kpc).

    Following the toy (Stage A/A.1) and the D2 spec, both are per-pulsar and
    CONDITIONAL ON THE OTHER DISTANCES; they differ only in the CW block:
      conditional_i = 1/sqrt(F_dd[i,i])                    (CW params known)
      marginal_i    = 1/sqrt(F_dd[i,i] - (F_dc F_cc^{-1} F_cd)_ii)
                                                           (CW params marginalised
                                                            = diagonal of the Schur
                                                            complement F_marg)
    Only the CW block F_cc (8*N_CW square) is inverted, via the toy's eigh
    pseudo-inverse (rcond*max|eig| cut). No 116x116 inverse -> rcond-stable.
    (The full-covariance variant sqrt(diag(inv(F_marg))) additionally
    marginalises over the OTHER distances; that inverse is near-singular for a
    single CW and rcond-unstable, so it is NOT the reported metric.)
    """
    F = -np.asarray(H)  # Fisher = -Hessian (NSD Hessian -> PSD Fisher)
    F = 0.5 * (F + F.T)
    F_dd = F[:n_dist, :n_dist]
    F_dc = F[:n_dist, n_dist:]
    F_cc = F[n_dist:, n_dist:]

    diag_dd = np.diag(F_dd)
    with np.errstate(divide="ignore", invalid="ignore"):
        cond = 1.0 / np.sqrt(diag_dd)  # kpc

    F_cc_inv = eigh_pinv(F_cc, rcond)
    # per-pulsar Schur diagonal: diag(F_dc F_cc^{-1} F_cd)
    schur_diag = np.einsum("ic,cd,id->i", F_dc, F_cc_inv, F_dc)
    marg_info = diag_dd - schur_diag           # PSD => >= 0
    marg_info = np.clip(marg_info, 0.0, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        marg = 1.0 / np.sqrt(marg_info)         # kpc; inf where fully degenerate
    return cond, marg, marg_info


# ============================================================
# CLI
# ============================================================

def _build_pta_and_psrs(n_psr, num_cw):
    ent_psrs, disco_psrs = load_pulsars(n_psr)
    pta, cw_block_names, _ = build_enterprise_pta(
        ent_psrs, num_cw, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
    return ent_psrs, disco_psrs, pta, cw_block_names


def run_gate():
    """Amortised zero-noise Fisher H == D1 dense Hessian on zero-noise data (5 psr, 1 CW)."""
    from stagec_hessian import build_stagec_problem, stagec_hessian
    print("gate: 5 psr, 1 CW, zero-noise data", flush=True)
    ent_psrs, disco_psrs, pta, cwb = _build_pta_and_psrs(5, 1)
    inj = make_geometry_injection(pta, ent_psrs, 1, cwb, seed=SEED)

    amo = build_fisher_amortised(disco_psrs, 1, inj, cwb)
    theta = amo["theta_truth"]
    data = amo["inject_data"](theta)
    # residual(truth)=0 check
    ll = float(amo["fisher_logL"](theta, data))
    H_amo = np.array(amo["fisher_hessian"](theta, data))
    n_dist = amo["n_dist"]
    negc = int((np.diag(H_amo)[:n_dist] < 0).sum())
    print(f"  amortised: logL(truth,zero-noise)={ll:.6e}  neg-curv dist diag={negc}/{n_dist}")

    # D1 dense Hessian on the SAME zero-noise data: build a fast-scan problem whose
    # residual_map is the pure CW injection at truth.
    from cw_helpers import inject_noisefree_cw, build_fast_scan_likelihood
    cw_list = [{k: inj[f"cw_{k}"] if f"cw_{k}" in inj else inj[f"{cwb[0]}_{k}"]
                for k in CW_PARAM_BASE}]
    resid_map, _ = inject_noisefree_cw(disco_psrs, cw_list)
    logl_fn, pk, bv = build_fast_scan_likelihood(disco_psrs, resid_map, 1, inj, cwb,
                                                 components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
    # reorder amortised H to fast-scan param order
    sel = [pk.index(k) for k in amo["theta_keys"]]
    import jax as _jax
    H_d1_full = np.array(_jax.jit(_jax.hessian(lambda x: logl_fn(x)))(bv))
    H_d1 = H_d1_full[np.ix_(sel, sel)]
    Hmax = np.max(np.abs(H_d1))
    rel = np.max(np.abs(H_amo - H_d1)) / Hmax
    print(f"  amortised vs D1 dense (zero-noise, reordered): max rel dev = {rel:.3e} "
          f"({'PASS' if rel < 1e-8 else 'FAIL'} vs 1e-8)")
    negc_d1 = int((np.diag(H_d1)[:n_dist] < 0).sum())
    print(f"  D1 dense neg-curv dist diag={negc_d1}/{n_dist}")
    return 0 if rel < 1e-8 else 1


def run_d2(ndraws):
    print(f"=== D2: 116 psr, 1 CW, {ndraws} geometry draws ===", flush=True)
    ent_psrs, disco_psrs, pta, cwb = _build_pta_and_psrs(116, 1)
    em_prior_pc = np.array([p.pdist[1] for p in ent_psrs]) * KPC_TO_PC
    psr_names = [p.name for p in ent_psrs]

    # Build the amortised Fisher ONCE (fixed noise). Reuse across draws.
    inj0 = make_geometry_injection(pta, ent_psrs, 1, cwb, seed=1000)
    t0 = time.time()
    amo = build_fisher_amortised(disco_psrs, 1, inj0, cwb)
    n_dist = amo["n_dist"]
    theta_keys = amo["theta_keys"]
    # compile once on the reference draw
    theta0 = amo["theta_truth"]; data0 = amo["inject_data"](theta0)
    _ = np.array(amo["fisher_hessian"](theta0, data0))
    print(f"build+compile once: {time.time()-t0:.1f}s, peak {device_peak_gib():.2f} GiB", flush=True)

    rconds = [1e-8, 1e-10, 1e-12]
    cond_draws, marg_draws = [], []
    negc_list = []
    marg_rcond = {rc: [] for rc in rconds}
    warm_times = []

    for d in range(ndraws):
        inj = make_geometry_injection(pta, ent_psrs, 1, cwb, seed=2000 + d)
        temp = _enterprise_to_disco_params(inj, cwb)
        theta = jnp.array([temp[k] for k in theta_keys], dtype=jnp.float64)
        data = amo["inject_data"](theta)
        t = time.time()
        H = np.array(amo["fisher_hessian"](theta, data))
        warm_times.append(time.time() - t)
        negc_list.append(int((np.diag(H)[:n_dist] < 0).sum()))
        cond, marg, _ = sigma_L_from_hessian(H, n_dist, rcond=1e-10)
        cond_draws.append(cond * KPC_TO_PC)   # pc
        marg_draws.append(marg * KPC_TO_PC)   # pc
        for rc in rconds:
            _, m_rc, _ = sigma_L_from_hessian(H, n_dist, rcond=rc)
            marg_rcond[rc].append(m_rc * KPC_TO_PC)
        ndegen = int(np.sum(~np.isfinite(marg)))
        print(f"  draw {d}: warm {warm_times[-1]:.2f}s neg-curv {negc_list[-1]}/{n_dist} "
              f"median marg sigma_L={np.median(marg[np.isfinite(marg)]*KPC_TO_PC):.4g} pc "
              f"degenerate(marg=inf)={ndegen}", flush=True)

    cond_draws = np.array(cond_draws)   # (ndraws, npsr) pc
    marg_draws = np.array(marg_draws)
    # per-pulsar medians + 16/84 across geometry draws
    # per-pulsar median/percentiles over draws, ignoring non-finite (degenerate) draws
    def nan_pct(a, q):
        a = np.where(np.isfinite(a), a, np.nan)
        return np.nanpercentile(a, q, axis=0)
    cond_med = np.median(cond_draws, axis=0)
    marg_med = nan_pct(marg_draws, 50)
    marg_lo = nan_pct(marg_draws, 16)
    marg_hi = nan_pct(marg_draws, 84)
    ratio = marg_draws / cond_draws
    n_degen_per_draw = np.sum(~np.isfinite(marg_draws), axis=1)

    # headline: pulsars with marginal sigma_L below EM prior (per-draw count + median)
    below = marg_draws < em_prior_pc[None, :]
    below_count_per_draw = below.sum(axis=1)

    # rcond sensitivity of per-pulsar median marginal sigma_L (finite pulsars only)
    def med_over_draws(arr):
        a = np.where(np.isfinite(arr), arr, np.nan)
        return np.nanmedian(a, axis=0)
    rcond_med = {rc: med_over_draws(np.array(marg_rcond[rc])) for rc in rconds}
    fin = np.isfinite(rcond_med[rconds[0]]) & np.isfinite(rcond_med[rconds[-1]]) \
        & (rcond_med[rconds[1]] > 0)
    rcond_spread = float(np.max(np.abs(rcond_med[rconds[0]][fin] - rcond_med[rconds[-1]][fin])
                                / rcond_med[rconds[1]][fin]))

    np.savez(
        "/home/mattm/projects/HSYMT/CW_transition/stagec_d2_results.npz",
        psr_names=np.array(psr_names), em_prior_pc=em_prior_pc,
        cond_draws=cond_draws, marg_draws=marg_draws,
        cond_med=cond_med, marg_med=marg_med, marg_lo=marg_lo, marg_hi=marg_hi,
        ratio=ratio, below_count_per_draw=below_count_per_draw,
        negc_list=np.array(negc_list), warm_times=np.array(warm_times),
    )

    # ---- report ----
    print("\n================ D2 SUMMARY ================")
    print(f"neg-curv count per draw (expect {n_dist}/{n_dist}): "
          f"min {min(negc_list)} max {max(negc_list)} "
          f"({'all NSD' if min(negc_list)==n_dist else 'SOME POSITIVE — see note'})")
    print(f"warm Hessian eval: median {np.median(warm_times):.2f}s "
          f"(compile paid once) — amortisation working")

    # 10-best table by median marginal sigma_L
    order = np.argsort(marg_med)[:10]
    print("\n10 best pulsars (by median marginal sigma_L across draws):")
    print(f"{'pulsar':14s} {'marg_pc':>10s} {'cond_pc':>10s} {'EM_pc':>10s} "
          f"{'marg/cond':>10s} {'marg<EM?':>8s}")
    for i in order:
        flag = "YES" if marg_med[i] < em_prior_pc[i] else "no"
        print(f"{psr_names[i]:14s} {marg_med[i]:10.4g} {cond_med[i]:10.4g} "
              f"{em_prior_pc[i]:10.4g} {np.median(ratio[:,i]):10.4g} {flag:>8s}")

    print(f"\ndegenerate (marginal sigma_L = inf, distance fully absorbed by CW): "
          f"median {int(np.median(n_degen_per_draw))} of {n_dist} per draw "
          f"(range {n_degen_per_draw.min()}-{n_degen_per_draw.max()})")
    print(f"HEADLINE: pulsars with MARGINAL sigma_L < EM prior: "
          f"median {int(np.median(below_count_per_draw))} of {n_dist} "
          f"(range {below_count_per_draw.min()}-{below_count_per_draw.max()} across draws)")
    r = ratio[np.isfinite(ratio)].flatten()
    print(f"marg/cond distribution (finite, all psr x draws): median {np.median(r):.4g}, "
          f"16/84 {np.percentile(r,16):.4g}/{np.percentile(r,84):.4g}, "
          f"min {r.min():.4g} max {r.max():.4g}")
    print(f"rcond sensitivity of median marginal sigma_L "
          f"({rconds[0]:.0e} vs {rconds[-1]:.0e}): max rel spread = {rcond_spread:.3e}")
    print("saved stagec_d2_results.npz")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["gate", "d2"])
    ap.add_argument("--ndraws", type=int, default=10)
    args = ap.parse_args()
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    if args.cmd == "gate":
        return run_gate()
    return run_d2(args.ndraws)


if __name__ == "__main__":
    sys.exit(main())
