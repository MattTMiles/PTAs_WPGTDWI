"""Track B — B0.1/B0.2: cold-start DAEM fringe-identification estimator.

Deterministic-Annealing EM (Ueda & Nakano 1998), annealing in pulsar-term COHERENCE
(Stage A.2). Reaches the Anchor-Census P(true) ceilings from a COLD start (sources drawn
from the prior; no truth beyond the EM priors). Spec: trackB_estimator_spec.md.

Design (see spec):
  * MODEL: N_CW sources (8 params each) + 116 distances; amortised zero-noise fast-scan
    logL (stagec_fisher.build_fisher_amortised). GP hyperparams frozen at truth.
  * COHERENCE ANNEAL: two amortised likelihoods with IDENTICAL theta layout ---
    amo_earth (include_pterm=False, fully-decohered limit) and amo_full (pterm on,
    coherent limit). lnL_beta = (1-beta)*lnL_earth + beta*lnL_full, beta=exp(-taubar/2 tc)
    in (0,1]. beta->0 unimodal Earth-term source recovery; beta->1 full A2 likelihood.
    The dist keys are inert in amo_earth (delay ignores p_dist), so the theta vector and
    the compiled shape are shared.
  * E-STEP (per pulsar, star-topology exact): scan B fringe centres on amo_full; the A2
    Bayesian softmax verbatim with the beta coherence weight:
        q_p(n) prop exp[ beta*(lnL_p(n)-lnL_p*) - (L_p(n)-L0_p)^2/(2 sigma_EM^2) ]
    P(true)_p = q_p(n_true). At beta=1 this is classify_full's P_true exactly.
  * M-STEP (soft/hard fringe-weighted): hard-assign each pulsar's MAP fringe centre
    (snapped to local conditional max), ascend lnL_beta on the 8*N_CW source params
    (Adam, param-scaled). Distances held at MAP; their earth-term gradient is ~0, the
    full-term gradient sharp -> at low beta theta moves by the Earth term (annealing).
  * MULTI-START: >=8 cold prior draws, staged coarse->kill->polish. Basin stats.

Run in jug-gpu:
  python trackB_estimator.py validate-estep      # E-step at truth reproduces A2 pop ceilings
  python trackB_estimator.py split-cost          # precompute-split E-step cost before/after
  python trackB_estimator.py asimov --config pop # B0.2 gate (pop / eq8 / eq16)
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import (load_pulsars, build_enterprise_pta, compute_mode_spacing,
                        _enterprise_to_disco_params, generate_injection_params)
from stagec_fisher import (build_fisher_amortised, make_geometry_injection,
                           N_COMPONENTS, LOG10_EQUAD, CW_PARAM_BASE, KPC_TO_PC,
                           GWB_LOG10_A, GWB_GAMMA, LOG10_FGW_RANGE, LOG10_MC_RANGE)
import stagec_anchor_a2 as A2
from discovery.deterministic import make_phase_connected_binary

CWT = "/home/mattm/projects/HSYMT/CW_transition"


class EarthDelay:
    """Earth-term-only multi-source CW delay (the fully-decohered coherence limit).

    Mirrors MultiSourceDelay but calls the pterm=False waveform (p_dist is keyword-only
    there, so it is NOT passed positionally). Keeps a nominal p_dist param key so the
    theta layout matches the full model; the key is inert (delay ignores it)."""
    _GLOBALS = ("cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
                "log10_fgw", "log10_h", "phase0", "psi")

    def __init__(self, psr, n_sources):
        self.psr = psr; self.n_sources = n_sources
        src = make_phase_connected_binary(pulsarterm=False)
        self._vmapped = jax.vmap(src, in_axes=(None, None, 0, 0, 0, 0, 0, 0, 0, 0))
        self._suffixes = ["" if i == 0 else f"_{i+1}" for i in range(n_sources)]
        params = []
        for suf in self._suffixes:
            for key in self._GLOBALS:
                nm = f"cw_{key}{suf}"
                if nm not in params:
                    params.append(nm)
        params.append(f"{psr.name}_cw_p_dist")   # inert, keeps theta layout aligned
        self.params = params

    def __call__(self, params):
        cw = {k: jnp.stack([params[f"cw_{k}{s}"] for s in self._suffixes]) for k in self._GLOBALS}
        d = self._vmapped(self.psr.toas, self.psr.pos, cw["cos_gwtheta"], cw["gwphi"],
                          cw["cos_inc"], cw["log10_mc"], cw["log10_fgw"], cw["log10_h"],
                          cw["phase0"], cw["psi"])
        return jnp.sum(d, axis=0)
EQUAL_H = -13.75
B = A2.B  # 512 evals per pulsar

# config -> (ncw, injection kwargs, seed)
CONFIGS = {
    "pop":  dict(ncw=16, seed=3000, population=(3, -13.25, -14.25)),
    "eq8":  dict(ncw=8,  seed=2000, h_range=(EQUAL_H, EQUAL_H)),
    "eq16": dict(ncw=16, seed=2000, h_range=(EQUAL_H, EQUAL_H)),
}


# ============================================================
# Problem construction: two amortised likelihoods sharing theta layout
# ============================================================
def build_problem(config, verbose=True, build_earth=True):
    ncw = CONFIGS[config]["ncw"]
    ent_psrs, disco_psrs = load_pulsars(116)
    names = [p.name for p in ent_psrs]
    priors = A2.load_real_priors(names)
    pta, cwb, _ = build_enterprise_pta(ent_psrs, ncw, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)

    # truth injection for THIS config
    kw = {k: v for k, v in CONFIGS[config].items() if k not in ("ncw",)}
    inj = make_geometry_injection(pta, ent_psrs, ncw, cwb, **kw)

    # amortised full likelihood (pterm on) built at the truth injection (frozen noise)
    t0 = time.time()
    amo_full = build_fisher_amortised(disco_psrs, ncw, inj, cwb)
    # amortised earth-only likelihood: SAME injection/noise, pterm off. Build a private
    # variant by monkey-swapping include_pterm on the MultiSourceDelay via a flag.
    # build_earth=False (split/probe work) skips this second compile (~30-60s saved); only
    # trackB_estimator's own make_mstep coherence-anneal reads amo_earth.
    amo_earth = (build_fisher_amortised(disco_psrs, ncw, inj, cwb, msd_factory=EarthDelay)
                 if build_earth else None)
    if verbose:
        print(f"  build amo_full+amo_earth: {time.time()-t0:.1f}s", flush=True)

    # truth theta + observed (zero-noise Asimov) data on the FULL model
    temp = _enterprise_to_disco_params(inj, cwb)
    theta_truth = np.array([temp[k] for k in amo_full["theta_keys"]], dtype=np.float64)
    data_obs = amo_full["inject_data"](jnp.asarray(theta_truth))

    # per-pulsar fringe grid dL (kpc), AI index map, EV eval grid, L0 prior means
    cw_list = [{k: inj[f"{nm}_{k}"] for k in CW_PARAM_BASE} for nm in cwb]
    npsr = len(ent_psrs)
    dL = np.array([min(compute_mode_spacing(cw["cos_gwtheta"], cw["gwphi"], cw["log10_fgw"],
                       disco_psrs[a].pos) for cw in cw_list) for a in range(npsr)])
    idx_of = {k: i for i, k in enumerate(amo_full["theta_keys"])}
    AI = np.array([idx_of[f"{pe.name}_cw_p_dist"] for pe in ent_psrs])
    L0 = np.array([pe.pdist[0] for pe in ent_psrs])
    EV = np.empty((npsr, B))
    for a, pe in enumerate(ent_psrs):
        R = 3.0 * max(priors["feather"][a], priors["script"][a], priors["lit"][a])
        EV[a] = A2.eval_grid(pe.pdist[0], dL[a], R)

    vmap2_full = A2.make_vmap2(amo_full)
    # warm compiles
    _ = amo_full["fisher_logL"](jnp.asarray(theta_truth), data_obs)
    if amo_earth is not None:
        _ = amo_earth["fisher_logL"](jnp.asarray(theta_truth), data_obs)
    _ = vmap2_full(jnp.zeros((1, B, amo_full["n_theta"])), A2.stack_data([data_obs]))

    n_dist = amo_full["n_dist"]
    n_cw_params = amo_full["n_theta"] - n_dist
    return dict(config=config, ncw=ncw, ent_psrs=ent_psrs, disco_psrs=disco_psrs, names=names,
                pta=pta, cwb=cwb, priors=priors, inj=inj,
                amo_full=amo_full, amo_earth=amo_earth, vmap2_full=vmap2_full,
                theta_truth=theta_truth, data_obs=data_obs,
                dL=dL, AI=AI, L0=L0, EV=EV, n_dist=n_dist, n_cw_params=n_cw_params, npsr=npsr)


# ============================================================
# F-statistic seed scan (TRUTH-BLIND): coarse freq x sky matched filter,
# Earth-term single-source likelihood profiled over amplitude/phase/psi/inc.
# Places the LOUD sources near truth before the DAEM handoff. No injection-truth
# info enters the grid, the thresholds, or the peak selection.
# ============================================================
FGW_BAND = (-8.0, -7.5)              # log10 Hz (10-32 nHz), the D2 operating band
SEED_MC = 9.0                        # fixed reference chirp mass for the scan
SEED_NSIDE = 4                       # healpix sky grid (192 px, ~14.7 deg; M-step polishes)


def freq_grid(T):
    """log10_fgw grid at 1/(3T) linear spacing across the band."""
    df = 1.0 / (3.0 * T)
    f = np.arange(10 ** FGW_BAND[0], 10 ** FGW_BAND[1] + df, df)
    return np.log10(f)


def sky_grid(nside):
    """healpix pixel centres -> (cos_gwtheta, gwphi)."""
    import healpy as hp
    npix = hp.nside2npix(nside)
    th, ph = hp.pix2ang(nside, np.arange(npix))   # colatitude, longitude
    return np.cos(th), ph


def build_earth_single(P):
    """N_CW=1 Earth-term amortised likelihood (EarthDelay), sharing P's frozen noise
    and fed P's observed data. Used only for the seed scan."""
    ent, disco, cwb16 = P["ent_psrs"], P["disco_psrs"], P["cwb"]
    pta1, cwb1, _ = build_enterprise_pta(ent, 1, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
    inj1 = make_geometry_injection(pta1, ent, 1, cwb1, seed=1000, h_range=(EQUAL_H, EQUAL_H))
    earth1 = build_fisher_amortised(disco, 1, inj1, cwb1, msd_factory=EarthDelay)
    return earth1


# free = [cos_inc, log10_h, phase0, psi]; bounds for the profile
_FREE_LO = np.array([-1.0, -16.0, 0.0, 0.0])
_FREE_HI = np.array([1.0, -12.0, 2 * np.pi, np.pi])


def make_fstat(P, earth1):
    """Return a jitted, vmapped profiled-statistic over grid points.
    Each point (cos_gwtheta, gwphi, log10_fgw) -> max over [cos_inc,log10_h,phase0,psi]
    of the Earth-term single-source lnL (matched filter, amplitude/phase profiled)."""
    fl = earth1["fisher_logL"]; data = P["data_obs"]
    L0 = jnp.asarray(P["L0"])
    flo = jnp.asarray(_FREE_LO); fhi = jnp.asarray(_FREE_HI)
    fscale = jnp.asarray([0.5, 1.0, 3.14, 1.57])

    def lnl_at(cosgt, gwphi, lfgw, free):
        cw = jnp.array([cosgt, gwphi, free[0], SEED_MC, lfgw, free[1], free[2], free[3]])
        return fl(jnp.concatenate([L0, cw]), data)

    def profile_point(cosgt, gwphi, lfgw):
        g = jax.grad(lambda fr: lnl_at(cosgt, gwphi, lfgw, fr))
        free = jnp.array([0.0, -14.0, 3.14, 1.57])
        m = jnp.zeros(4); v = jnp.zeros(4)
        b1, b2, eps = 0.9, 0.999, 1e-8
        def body(t, c):
            free, m, v = c
            gr = jnp.nan_to_num(g(free))
            m = b1 * m + (1 - b1) * gr; v = b2 * v + (1 - b2) * gr * gr
            mh = m / (1 - b1 ** t); vh = v / (1 - b2 ** t)
            free = jnp.clip(free + 0.05 * fscale * mh / (jnp.sqrt(vh) + eps), flo, fhi)
            return (free, m, v)
        free, _, _ = jax.lax.fori_loop(1, 41, body, (free, m, v))
        return lnl_at(cosgt, gwphi, lfgw, free), free

    return jax.jit(jax.vmap(profile_point))


def seed_scan(P, verbose=True):
    """Run the profiled matched-filter over the full (freq x sky) grid, return the
    statistic map + per-point best free params + the grid axes. TRUTH-BLIND."""
    earth1 = build_earth_single(P)
    fstat = make_fstat(P, earth1)
    T = 6.992e8  # array span (getspan); fixed here to avoid recompute
    lfgws = freq_grid(T); cosgt, gwphi = sky_grid(SEED_NSIDE)
    nf, ns = len(lfgws), len(cosgt)
    # baseline (no-signal) lnL for the detection statistic 2F = 2(lnL_max - lnL_0)
    ll0 = float(earth1["fisher_logL"](jnp.concatenate([jnp.asarray(P["L0"]),
                  jnp.array([0., 0., 0., SEED_MC, np.mean(lfgws), -18.0, 0., 0.])]), P["data_obs"]))
    # build the full grid (nf*ns points), chunk through the vmap
    CG, GP, LF = [], [], []
    for i in range(nf):
        for j in range(ns):
            CG.append(cosgt[j]); GP.append(gwphi[j]); LF.append(lfgws[i])
    CG = np.array(CG); GP = np.array(GP); LF = np.array(LF)
    npt = len(CG)
    stat = np.empty(npt); free = np.empty((npt, 4))
    t0 = time.time(); CH = 512
    for c0 in range(0, npt, CH):
        c1 = min(c0 + CH, npt)
        s, fr = fstat(jnp.asarray(CG[c0:c1]), jnp.asarray(GP[c0:c1]), jnp.asarray(LF[c0:c1]))
        stat[c0:c1] = np.asarray(s); free[c0:c1] = np.asarray(fr)
    twoF = 2.0 * (stat - ll0)
    if verbose:
        print(f"  seed_scan: grid {nf} freq x {ns} sky = {npt} pts, {time.time()-t0:.1f}s; "
              f"2F max {twoF.max():.1f} median {np.median(twoF):.1f}", flush=True)
    return dict(twoF=twoF.reshape(nf, ns), free=free.reshape(nf, ns, 4),
                lfgws=lfgws, cosgt=cosgt, gwphi=gwphi, ll0=ll0, npt=npt)


# ============================================================
# E-step: per-pulsar fringe posterior at current theta estimate
# ============================================================
def scan_current(P, theta_full):
    """lnL_p(n) for all pulsars at the current theta estimate (distances in theta_full
    are the 'other' distances held fixed; each pulsar's own is swept over EV)."""
    LNL = A2.scan_all(P["vmap2_full"], np.asarray(theta_full)[None],
                      A2.stack_data([P["data_obs"]]), P["EV"][None], P["AI"])[0]
    return LNL  # (npsr, B)


def fringe_posterior(P, LNL, ll_at_truthfringe, priors_col, beta):
    """Per-pulsar softmax over sampled fringe centres (A2 classify_full math + beta).

    Returns dict of arrays (npsr,): P_true, map_offset (kpc from L0), map_lnL,
    map_evalL, tv_ready q stored as (offsets,weights) list, K_counted, dlnL.
    ll_at_truthfringe: lnL at the TRUE fringe (for P_true). For a cold estimate we do
    not know the truth fringe; P_true is reported only at eval/gate time using the true
    distance's nearest fringe (see p_true_of).
    """
    npsr = P["npsr"]
    out = dict(P_true=np.zeros(npsr), map_evalL=np.zeros(npsr), map_lnL=np.zeros(npsr),
               map_k=np.zeros(npsr, int), K_counted=np.zeros(npsr, int), dlnL=np.zeros(npsr),
               map_qmax=np.zeros(npsr), qlist=[None] * npsr)
    for a, pe in enumerate(P["ent_psrs"]):
        evalL, lnL = P["EV"][a], LNL[a]
        L0 = pe.pdist[0]; dLa = P["dL"][a]; sig = priors_col[a]
        k_all = np.round((evalL - L0) / dLa).astype(int)
        kmax_win = int(np.floor(3.0 * sig / dLa + 1e-9))
        # unique fringes among sampled points, per-fringe peak lnL (like classify_full)
        order = np.argsort(k_all)
        ks, ll_s = k_all[order], lnL[order]
        uk, idx0 = np.unique(ks, return_index=True)
        lnL_u = np.maximum.reduceat(ll_s, idx0)
        offs_u = uk * dLa
        # coherence-annealed Bayesian weights over ALL sampled fringes
        logw = beta * (lnL_u - lnL_u.max()) - offs_u ** 2 / (2.0 * sig ** 2)
        w = np.exp(logw - logw.max()); w /= w.sum()
        jmap = int(np.argmax(w))
        out["map_evalL"][a] = L0 + uk[jmap] * dLa
        out["map_k"][a] = uk[jmap]
        out["map_lnL"][a] = lnL_u[jmap]
        out["map_qmax"][a] = float(w[jmap])
        out["qlist"][a] = (uk.copy(), offs_u.copy(), lnL_u.copy(), w.copy())
        out["K_counted"][a] = 2 * kmax_win
    return out


def p_true_at(P, LNL, priors_col, beta=1.0):
    """P(true fringe) per pulsar at the current LNL, using the TRUE distance's fringe as
    'true' (gate/eval only). Mirrors A2 classify_full's P_true with the beta weight."""
    npsr = P["npsr"]; Pt = np.zeros(npsr); dlnL = np.zeros(npsr); Kc = np.zeros(npsr, int)
    for a, pe in enumerate(P["ent_psrs"]):
        evalL, lnL = P["EV"][a], LNL[a]
        L0 = pe.pdist[0]; dLa = P["dL"][a]; sig = priors_col[a]
        kmax_win = int(np.floor(3.0 * sig / dLa + 1e-9)); Kc[a] = 2 * kmax_win
        k_all = np.round((evalL - L0) / dLa).astype(int)
        sel = (np.abs(k_all) >= 1) & (np.abs(k_all) <= kmax_win)
        if kmax_win >= 1 and sel.any():
            kk, lw = k_all[sel], lnL[sel]
            o = np.argsort(kk); kk, lw = kk[o], lw[o]
            uk, i0 = np.unique(kk, return_index=True)
            lnL_u = np.maximum.reduceat(lw, i0); offs = uk * dLa
            # weight of true fringe (k=0) relative to wrong fringes, at current LNL
            ll_true_here = float(lnL[np.abs(k_all).argmin()])
            logw = beta * (lnL_u - ll_true_here) - offs ** 2 / (2.0 * sig ** 2)
            Pt[a] = 1.0 / (1.0 + float(np.exp(logw).sum()))
            dlnL[a] = ll_true_here - float(lnL_u.max())
        else:
            Pt[a] = 1.0; dlnL[a] = np.inf
    return Pt, dlnL, Kc


# ============================================================
# M-step: ascend the coherence-blended lnL on the source params
# ============================================================
def make_mstep(P):
    """M-step = Adam ascent on the coherence-blended lnL, with the WHOLE Adam loop
    inside one JITted lax.fori_loop (60 device dispatches -> 1 compiled call). The
    per-cw-param scale multiplies the step so mixed-scale params (angles O(1) vs logs)
    move commensurately. nsteps is dynamic (traced) -> no recompile per rung."""
    fl_full = P["amo_full"]["fisher_logL"]
    fl_earth = P["amo_earth"]["fisher_logL"]
    data = P["data_obs"]
    step_scale = jnp.asarray(phi_scale(P))
    plo, phi_hi = phi_bounds(P)
    plo = jnp.asarray(plo); phi_hi = jnp.asarray(phi_hi)

    def obj(phi, base_L, beta):
        # M-step ascends the FULL likelihood at the current (confidence-gated) fringes.
        # Its maximum is truth (zero-noise Asimov) -> truth is a fixed point (no earth-
        # model bias). Annealing lives in the E-step temperature beta + the confidence
        # gate on distance reassignment, NOT in a decohered-model blend. beta unused here.
        theta = jnp.concatenate([base_L, phi])
        return fl_full(theta, data)

    grad_obj = jax.grad(obj, argnums=0)
    b1, b2, eps = 0.9, 0.999, 1e-8

    @jax.jit
    def adam_run(phi0, base_L, beta, nsteps, lr):
        def body(t, carry):
            phi, m, v = carry
            g = grad_obj(phi, base_L, beta)
            g = jnp.nan_to_num(g, nan=0.0, posinf=0.0, neginf=0.0)    # singular grad -> no move
            m = b1 * m + (1 - b1) * g
            v = b2 * v + (1 - b2) * g * g
            mh = m / (1 - b1 ** t); vh = v / (1 - b2 ** t)
            step = lr * step_scale * mh / (jnp.sqrt(vh) + eps)
            step = jnp.clip(step, -0.5 * step_scale, 0.5 * step_scale)  # trust region
            phi = jnp.clip(phi + step, plo, phi_hi)                     # ASCEND + keep valid
            return (phi, m, v)
        z = jnp.zeros_like(phi0)
        phi, _, _ = jax.lax.fori_loop(1, nsteps + 1, body, (jnp.clip(phi0, plo, phi_hi), z, z))
        return phi, obj(phi, base_L, beta)

    def mstep(phi, base_L, beta, nsteps, lr):
        phi_out, f = adam_run(jnp.asarray(phi), jnp.asarray(base_L),
                              float(beta), int(nsteps), float(lr))
        return np.asarray(phi_out), float(f)

    return mstep, adam_run


def snap_distances(P, phi, base_L):
    """Parabolic within-fringe refine of each distance on amo_full lnL (3-pt)."""
    fl = P["amo_full"]["fisher_logL"]; data = P["data_obs"]
    base_L = np.asarray(base_L).copy()
    for a in range(P["npsr"]):
        d = P["dL"][a] / 8.0
        Ls = np.array([base_L[a] - d, base_L[a], base_L[a] + d])
        ys = []
        for L in Ls:
            th = np.concatenate([base_L, phi]); th[a] = L
            ys.append(float(fl(jnp.asarray(th), data)))
        y0, y1, y2 = ys
        denom = (y0 - 2 * y1 + y2)
        if denom < 0:  # concave -> interior max
            shift = 0.5 * (y0 - y2) / denom * d
            base_L[a] = base_L[a] + np.clip(shift, -d, d)
    return base_L


# ============================================================
def _unit(cosgt, gwphi):
    s = np.sqrt(np.clip(1 - cosgt ** 2, 0, 1))
    return np.array([s * np.cos(gwphi), s * np.sin(gwphi), cosgt])


def pick_seeds(scan, twoF_thresh, min_sep_deg=25.0, min_df_bins=1.5):
    """Non-maximum suppression on the (freq x sky) 2F map: greedily take peaks above
    threshold, each separated from stronger accepted peaks by > min_sep in SKY (any frequency).
    STEP-2 REPAIR (2026-07-07): SKY-ONLY exclusion (was sky-AND-freq). The single-source F-stat
    map throws cross-frequency sidelobes of the loud sources (2F 285-350) that sit at ~same sky /
    different freq; the old sky-AND-freq rule kept them (freq differed), burying weaker true
    sources -- loud2's true peak sank to rank-13. Sky-only exclusion at 25deg kills those sidelobes
    and surfaces loud2 (11.88 deg, rank ~11); GATE 3/3 loud captured on the pop draw. min_df_bins
    retained only as a param for callers; not used in the sky-only rule."""
    twoF = scan["twoF"]; lfgws = scan["lfgws"]; cosgt = scan["cosgt"]; gwphi = scan["gwphi"]
    nf, ns = twoF.shape
    flat = [(twoF[i, j], i, j) for i in range(nf) for j in range(ns) if twoF[i, j] > twoF_thresh]
    flat.sort(reverse=True)
    df = np.median(np.diff(lfgws))
    seeds = []
    for val, i, j in flat:
        u = _unit(cosgt[j], gwphi[j])
        ok = True
        for s in seeds:
            ang = np.degrees(np.arccos(np.clip(np.dot(u, s["u"]), -1, 1)))
            if ang < min_sep_deg:                       # sky-only sidelobe suppression
                ok = False; break
        if ok:
            fr = scan["free"][i, j]
            seeds.append(dict(twoF=val, lfgw=lfgws[i], cosgt=cosgt[j], gwphi=gwphi[j],
                              cos_inc=fr[0], log10_h=fr[1], phase0=fr[2], psi=fr[3],
                              u=u, log10_mc=SEED_MC))
    return seeds


def seed_cw_vec(s):
    return np.array([s["cosgt"], s["gwphi"], s["cos_inc"], s["log10_mc"],
                     s["lfgw"], s["log10_h"], s["phase0"], s["psi"]])


def seeds_to_phi(P, seeds, filler_seed=9100):
    """Build an N_CW phi: seeds in the first slots, faint cold-prior filler in the rest."""
    ncw = P["ncw"]
    phi = draw_cold_phi(P, seed=filler_seed).copy()
    for k in range(min(len(seeds), ncw)):
        phi[8 * k:8 * k + 8] = seed_cw_vec(seeds[k])
    # faint filler: keep filler amplitudes low so they don't fight the seeds
    for k in range(len(seeds), ncw):
        phi[8 * k + 5] = -15.0
    return phi


def loud_truth(P):
    """The population's loud sources (first k_loud) truth (cos_gwtheta, gwphi, log10_fgw)."""
    inj = P["inj"]; cwb = P["cwb"]
    k = CONFIGS[P["config"]].get("population", (P["ncw"],))[0]
    return [dict(cosgt=inj[f"{cwb[i]}_cos_gwtheta"], gwphi=inj[f"{cwb[i]}_gwphi"],
                 lfgw=inj[f"{cwb[i]}_log10_fgw"], u=_unit(inj[f"{cwb[i]}_cos_gwtheta"],
                 inj[f"{cwb[i]}_gwphi"])) for i in range(k)]


# ============================================================
# DAEM driver (single start)
# ============================================================
def run_daem_start(P, phi0, priors_col, mstep, beta_ladder, iters_per_rung,
                   mstep_steps, lr, tol_phi=1e-3, verbose=False, on_iter=None):
    phi = np.asarray(phi0, float).copy()
    base_L = P["L0"].copy()                     # distances start at the EM prior mean (k=0)
    scale = phi_scale(P)
    q_thresh = 0.5                              # only reassign a fringe the posterior is sure of
    LNL = None; post = None; niter = 0
    for beta in beta_ladder:
        for it in range(iters_per_rung):
            theta = np.concatenate([base_L, phi])
            LNL = scan_current(P, theta)                       # E-step
            post = fringe_posterior(P, LNL, None, priors_col, beta)
            conf = post["map_qmax"] > q_thresh                 # confident pulsars only
            base_L = np.where(conf, post["map_evalL"], base_L) # else hold (no noisy MAP)
            phi_new, fval = mstep(phi, base_L, beta, mstep_steps, lr)   # M-step
            dphi = np.max(np.abs(phi_new - phi) / scale)
            phi = phi_new; niter += 1
            if verbose:
                print(f"    beta={beta:.3f} it={it} f={fval:.3f} dphi={dphi:.2e}", flush=True)
            if on_iter is not None:
                on_iter(beta, it, niter, phi, base_L, LNL, fval, dphi)
            if beta == beta_ladder[-1] and dphi < tol_phi:
                break
    # final E-step + snap
    theta = np.concatenate([base_L, phi])
    LNL = scan_current(P, theta)
    post = fringe_posterior(P, LNL, None, priors_col, 1.0)
    base_L = post["map_evalL"].copy()
    base_L = snap_distances(P, phi, base_L)
    theta = np.concatenate([base_L, phi])
    LNL = scan_current(P, theta)
    fval = float(P["amo_full"]["fisher_logL"](jnp.asarray(theta), P["data_obs"]))
    return dict(phi=phi, base_L=base_L, LNL=LNL, fval=fval, niter=niter)


def phi_scale(P):
    """Per-cw-param scale for the dphi convergence norm (angle O(1), log O(1) too here)."""
    ncw = P["ncw"]
    base = np.array([1.0, 3.14, 1.0, 0.5, 0.3, 0.3, 3.14, 1.57])  # per 8-param block
    return np.tile(base, ncw)


# per-8-param-block bounds: cos_gwtheta, gwphi, cos_inc, log10_mc, log10_fgw,
# log10_h, phase0, psi  (generous physical ranges; clipping prevents waveform nan)
_PHI_LO = np.array([-1.0, 0.0, -1.0, 7.0, -8.7, -16.0, 0.0, 0.0])
_PHI_HI = np.array([1.0, 2 * np.pi, 1.0, 10.0, -7.0, -12.0, 2 * np.pi, np.pi])


def phi_bounds(P):
    return np.tile(_PHI_LO, P["ncw"]), np.tile(_PHI_HI, P["ncw"])


def draw_cold_phi(P, seed):
    """Cold source params from the generative prior (no truth)."""
    rng = np.random.default_rng(seed); np.random.seed(seed)
    inj = generate_injection_params(P["pta"], P["ent_psrs"], P["ncw"], P["cwb"],
                                    log10_h=None, scenario="realistic", rng=rng,
                                    log10_h_range=(EQUAL_H, EQUAL_H),
                                    log10_fgw_range=LOG10_FGW_RANGE, log10_mc_range=LOG10_MC_RANGE,
                                    gwb_log10_A=GWB_LOG10_A, gwb_gamma=GWB_GAMMA)
    temp = _enterprise_to_disco_params(inj, P["cwb"])
    theta = np.array([temp[k] for k in P["amo_full"]["theta_keys"]], dtype=np.float64)
    return theta[P["n_dist"]:]  # cw slice only


# ============================================================
# Multi-start (staged) + eval
# ============================================================
def multistart(P, priors_col, n_start=8, coarse_rungs=None, full_ladder=None,
               iters_per_rung=3, mstep_steps=40, lr=8e-3, verbose=True, seed_phi=None):
    if full_ladder is None:
        full_ladder = [0.05, 0.1, 0.2, 0.4, 0.7, 1.0]
    if coarse_rungs is None:
        coarse_rungs = full_ladder[:3]
    mstep, _ = make_mstep(P)
    # start 0 = F-stat seed (if given); remaining = cold prior draws for basin diversity
    init_phis = ([seed_phi] if seed_phi is not None else []) + \
                [draw_cold_phi(P, seed=7000 + s) for s in range(n_start - (seed_phi is not None))]
    # STAGE 1: coarse EM for all starts
    starts = []
    t0 = time.time()
    for s, phi0 in enumerate(init_phis):
        r = run_daem_start(P, phi0, priors_col, mstep, coarse_rungs, iters_per_rung,
                           mstep_steps, lr, verbose=bool(os.environ.get("TB_VERBOSE")))
        starts.append(r)
        if verbose:
            print(f"  [coarse] start {s}: f={r['fval']:.2f}", flush=True)
    fvals = np.array([r["fval"] for r in starts])
    fbest = fvals.max()
    # STAGE 2: kill clearly-losing basins (gap > 50 nats below best), polish survivors
    keep = np.where(fvals > fbest - 50.0)[0]
    survivors = []
    for s in keep:
        phi0 = starts[s]["phi"]
        r = run_daem_start(P, phi0, priors_col, mstep, full_ladder, iters_per_rung,
                           mstep_steps, lr, verbose=bool(os.environ.get("TB_VERBOSE")))
        survivors.append(r)
        if verbose:
            print(f"  [polish] start {s}: f={r['fval']:.2f} niter={r['niter']}", flush=True)
    fpol = np.array([r["fval"] for r in survivors])
    best = survivors[int(fpol.argmax())]
    frac_best = float(np.mean(fpol > fpol.max() - 1.0))
    wall = time.time() - t0
    return dict(best=best, survivors=survivors, coarse_fvals=fvals,
                polish_fvals=fpol, n_start=n_start, n_survive=len(keep),
                frac_best=frac_best, wall=wall)


CEILING = {"J0711-6830": 0.953, "J1713+0747": 0.989, "J1909-3744": 0.998}


def eval_gate(P, LNL, priors_col):
    Pt, dlnL, Kc = p_true_at(P, LNL, priors_col, beta=1.0)
    names = P["names"]
    certified = [names[a] for a in range(P["npsr"]) if Pt[a] > 0.9]
    strict = [names[a] for a in range(P["npsr"]) if Pt[a] > 0.99]
    return dict(P_true=Pt, dlnL=dlnL, K_counted=Kc, certified=certified, strict=strict)


# ============================================================
# CLI modes
# ============================================================
def mode_validate_estep():
    """E-step at TRUTH theta must reproduce the A2 population ceilings (lit prior)."""
    print("=== validate-estep: pop draw, E-step at truth vs A2 ceilings ===", flush=True)
    P = build_problem("pop")
    LNL = scan_current(P, P["theta_truth"])
    ev = eval_gate(P, LNL, P["priors"]["lit"])
    Pt = ev["P_true"]; names = P["names"]
    print(f"\n  {'pulsar':14s} {'P_true(est)':>12s} {'A2 ceiling':>12s} {'|diff|':>8s}")
    ok = True
    for nm, ceil in CEILING.items():
        a = names.index(nm); diff = abs(Pt[a] - ceil)
        ok &= diff < 0.05
        print(f"  {nm:14s} {Pt[a]:12.3f} {ceil:12.3f} {diff:8.3f}")
    cert = ev["certified"]
    print(f"\n  certified (P>0.9): {cert}")
    print(f"  count P>0.9 = {len(cert)} (A2 pop bayes = 3)")
    verdict = "PASS" if (ok and set(cert) == set(CEILING)) else "CHECK"
    print(f"\n  [validate-estep] {verdict} (3 ceiling pulsars within 0.05 & certified set matches)")
    return 0


def mode_split_cost():
    """Precompute-split E-step cost: cached-outer-Cholesky (after) vs rebuild (before)."""
    print("=== split-cost: E-step per-fringe likelihood, before/after the QuickCW split ===", flush=True)
    P = build_problem("pop", verbose=True)
    from cw_helpers import build_disco_likelihood, inject_noisefree_cw
    # AFTER: amortised cached outer Cholesky, batched 116x512 fringe scan
    for _ in range(2):  # warm
        LNL = scan_current(P, P["theta_truth"])
    t = time.time(); LNL = scan_current(P, P["theta_truth"]); t_after = time.time() - t
    n_eval = P["npsr"] * B
    print(f"  AFTER  (cached outer Cholesky, amortised): {t_after:.3f}s for {n_eval} evals "
          f"= {t_after/n_eval*1e3:.3f} ms/eval", flush=True)
    # BEFORE: full discovery logL rebuilding the outer Cholesky each eval (no cache).
    inj = P["inj"]; cwb = P["cwb"]; disco = P["disco_psrs"]
    cw_list = [{k: inj[f"{nm}_{k}"] for k in CW_PARAM_BASE} for nm in cwb]
    resid_map, _ = inject_noisefree_cw(disco, cw_list)
    logl_fn, pk, bv = build_disco_likelihood(disco, resid_map, P["ncw"], inj, cwb,
                                             components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
    a_dist = pk.index(f"{P['ent_psrs'][0].name}_cw_p_dist")
    bv = np.asarray(bv)
    _ = float(logl_fn(jnp.asarray(bv)))  # warm compile
    n_before = 64
    t = time.time()
    for i in range(n_before):
        b2 = bv.copy(); b2[a_dist] = P["EV"][0][i % B]
        _ = float(logl_fn(jnp.asarray(b2)))
    t_before = (time.time() - t) / n_before
    print(f"  BEFORE (rebuild outer Cholesky per eval, build_disco_likelihood): "
          f"{t_before*1e3:.3f} ms/eval", flush=True)
    print(f"  SPLIT SPEEDUP = {t_before/(t_after/n_eval):.1f}x per eval "
          f"(cached outer Cholesky + batched scan)", flush=True)
    return 0


def report_seed_recovery(P, seeds, thresh):
    """Seeding gate 3a: how many loud/faint truth sources a seed captured, offsets."""
    loud = loud_truth(P); df = np.median(np.diff(freq_grid(6.992e8)))
    names_loud = [f"loud{i}" for i in range(len(loud))]
    print(f"  seeds above 2F>{thresh:.0f}: {len(seeds)} "
          f"(2F {[round(s['twoF'],1) for s in sorted(seeds,key=lambda x:-x['twoF'])[:6]]})", flush=True)
    nrec = 0
    print(f"  {'loud src':9s} {'ang off(deg)':>12s} {'df off(bins)':>12s} {'seed 2F':>9s} {'captured':>9s}")
    for i, lt in enumerate(loud):
        if not seeds:
            print(f"  {names_loud[i]:9s} {'--':>12s}"); continue
        angs = [np.degrees(np.arccos(np.clip(np.dot(lt["u"], s["u"]), -1, 1))) for s in seeds]
        dfs = [abs(lt["lfgw"] - s["lfgw"]) / df for s in seeds]
        cost = [a + 30 * d for a, d in zip(angs, dfs)]  # nearest in sky+freq
        j = int(np.argmin(cost))
        cap = angs[j] < 25.0 and dfs[j] < 2.0            # within M-step capture range
        nrec += cap
        print(f"  {names_loud[i]:9s} {angs[j]:12.1f} {dfs[j]:12.2f} {seeds[j]['twoF']:9.1f} {str(cap):>9s}")
    print(f"  LOUD recovered: {nrec}/{len(loud)}", flush=True)
    return nrec


def mode_seed_scan(config, thresh):
    print(f"=== seed-scan (F-stat seeding gate 3a): config={config} ===", flush=True)
    P = build_problem(config)
    scan = seed_scan(P)
    seeds = pick_seeds(scan, twoF_thresh=thresh)
    report_seed_recovery(P, seeds, thresh)
    return 0


def mode_asimov(config, n_start, save, seed=True):
    print(f"=== B0.2 Asimov gate: config={config} (cold-start DAEM"
          f"{' + F-stat seed' if seed else ''}) ===", flush=True)
    P = build_problem(config)
    prior_col = P["priors"]["lit"]
    seed_phi = None
    if seed:
        scan = seed_scan(P)
        seeds = pick_seeds(scan, twoF_thresh=25.0)
        nrec = report_seed_recovery(P, seeds, 25.0)
        seed_phi = seeds_to_phi(P, seeds)
    ms = multistart(P, prior_col, n_start=n_start, verbose=True, seed_phi=seed_phi)
    best = ms["best"]
    ev = eval_gate(P, best["LNL"], prior_col)
    Pt = ev["P_true"]; names = P["names"]

    # source recovery quality: converged phi vs truth
    phi_truth = P["theta_truth"][P["n_dist"]:]
    dphi = np.abs(best["phi"] - phi_truth) / phi_scale(P)
    print(f"\n  source recovery: max scaled |phi-truth| = {dphi.max():.3f}, "
          f"median {np.median(dphi):.3f}", flush=True)

    print(f"\n  {'pulsar':14s} {'P_true(cold)':>12s} {'A2 ceiling':>12s} {'|diff|':>8s}")
    ok = True
    for nm, ceil in CEILING.items():
        a = names.index(nm); diff = abs(Pt[a] - ceil)
        ok &= diff < 0.05
        print(f"  {nm:14s} {Pt[a]:12.3f} {ceil:12.3f} {diff:8.3f}")
    print(f"\n  certified (P>0.9): {ev['certified']}")
    print(f"  strict   (P>0.99): {ev['strict']}")
    extra = set(ev["certified"]) - set(CEILING)
    print(f"  extra beyond ceiling: {sorted(extra)}")
    print(f"\n  basin: {ms['n_survive']}/{ms['n_start']} survived coarse cut; "
          f"frac at best optimum = {ms['frac_best']:.2f}; wall {ms['wall']:.1f}s; "
          f"best niter {best['niter']}", flush=True)
    if config == "pop":
        verdict = "PASS" if (ok and not extra) else "FAIL"
        print(f"\n  [B0.2 {config}] GATE {verdict}", flush=True)
    if save:
        np.savez(f"{CWT}/trackB_asimov_{config}.npz",
                 names=np.array(names), P_true=Pt, dlnL=ev["dlnL"], K_counted=ev["K_counted"],
                 phi=best["phi"], phi_truth=phi_truth, base_L=best["base_L"],
                 coarse_fvals=ms["coarse_fvals"], polish_fvals=ms["polish_fvals"],
                 frac_best=ms["frac_best"], wall=ms["wall"], niter=best["niter"])
        print(f"  saved trackB_asimov_{config}.npz", flush=True)
    return 0


def mode_debug(config, seed, steps, iters, lr, ladder):
    """Single cold start with per-iter diagnostics vs truth (Asimov). Tune the M-step."""
    print(f"=== debug: config={config} seed={seed} steps={steps} iters={iters} lr={lr} "
          f"ladder={ladder} ===", flush=True)
    P = build_problem(config)
    prior_col = P["priors"]["lit"]
    mstep, _ = make_mstep(P)
    phi_truth = P["theta_truth"][P["n_dist"]:]
    scale = phi_scale(P)
    idx = {nm: P["names"].index(nm) for nm in CEILING}
    # loud sources = first 3 blocks (population sets first k_loud loud)
    n_loud = CONFIGS[config].get("population", (0,))[0] if config == "pop" else P["ncw"]
    loud_sl = slice(0, 8 * max(n_loud, 1))

    def on_iter(beta, it, niter, phi, base_L, LNL, fval, dphi):
        Pt, _, _ = p_true_at(P, LNL, prior_col, beta=1.0)
        eloud = np.max(np.abs(phi[loud_sl] - phi_truth[loud_sl]) / scale[loud_sl])
        ce = " ".join(f"{nm.split('-')[0].split('+')[0]}={Pt[idx[nm]]:.2f}" for nm in CEILING)
        ncert = int((Pt > 0.9).sum())
        print(f"  b={beta:.3f} it{it} f={fval:.1f} dphi={dphi:.2e} loudErr={eloud:.2f} "
              f"ncert={ncert} [{ce}]", flush=True)

    if seed < 0:   # init AT truth (machinery diagnostic: does it HOLD the ceilings?)
        phi0 = phi_truth.copy()
        print("  [init from TRUTH source params]", flush=True)
    else:
        phi0 = draw_cold_phi(P, seed=seed)
    t0 = time.time()
    r = run_daem_start(P, phi0, prior_col, mstep, ladder, iters,
                       steps, lr, verbose=False, on_iter=on_iter)
    ev = eval_gate(P, r["LNL"], prior_col)
    Pt = ev["P_true"]
    print(f"\n  FINAL ({time.time()-t0:.0f}s, niter={r['niter']}):", flush=True)
    for nm, ceil in CEILING.items():
        print(f"    {nm:14s} P_true={Pt[P['names'].index(nm)]:.3f} (ceiling {ceil})", flush=True)
    print(f"    certified P>0.9: {ev['certified']}", flush=True)
    print(f"    extra beyond ceiling: {sorted(set(ev['certified'])-set(CEILING))}", flush=True)
    return 0


def mode_perturb(config):
    """Certification TOLERANCE on source params: perturb truth loud sources by a range
    of scaled magnitudes, run the E-step only (no EM), report ceiling P_true vs delta.
    Measures whether the A2 ceilings (computed AT truth) survive finite source error."""
    print(f"=== perturb: certification tolerance on source params (config={config}) ===", flush=True)
    P = build_problem(config)
    prior_col = P["priors"]["lit"]
    phi_truth = P["theta_truth"][P["n_dist"]:].copy()
    scale = phi_scale(P)
    names = P["names"]
    n_loud = CONFIGS[config].get("population", (P["ncw"],))[0]
    loud_sl = slice(0, 8 * n_loud)
    # per-param subset: also test FREQUENCY-only vs ALL-param perturbation
    rng = np.random.default_rng(1234)
    print(f"  {'delta(scaled)':>13s} {'J0711':>7s} {'J1713':>7s} {'J1909':>7s}  {'ncert':>5s}  mode", flush=True)
    for mode, which in [("all-loud", None), ("freq-only", 4)]:
        for delta in [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
            phi = phi_truth.copy()
            if which is None:
                pert = rng.normal(size=8 * n_loud) * scale[loud_sl] * delta
                phi[loud_sl] += pert
            else:  # perturb only log10_fgw of each loud source
                for k in range(n_loud):
                    phi[8 * k + which] += rng.normal() * scale[8 * k + which] * delta
            theta = np.concatenate([P["L0"], phi])
            LNL = scan_current(P, theta)
            Pt, _, _ = p_true_at(P, LNL, prior_col, beta=1.0)
            vals = [Pt[names.index(nm)] for nm in CEILING]
            print(f"  {delta:13.1e} {vals[0]:7.3f} {vals[1]:7.3f} {vals[2]:7.3f}  "
                  f"{int((Pt>0.9).sum()):5d}  {mode}", flush=True)
    return 0


P1_CHUNK = 4      # fixed batched-draw size (padded) -> ONE scan_all compile shape


def p1_eval_trials(P, trial_thetas, verbose=False):
    """For each trial full-theta, run the E-step (distances swept), return per trial:
    joint lnL profiled over distances (each pulsar at its best-fringe centre), the
    registration count (#pulsars q_max>0.5), and the TRUE-fringe registration count
    (#pulsars q_max>0.5 AND best fringe == truth k=0). Chunked at a FIXED padded
    draw-size (D=P1_CHUNK) so scan_all compiles one shape and reuses it."""
    fl_b = P["amo_full"]["fisher_logL_batched"]
    ndist = P["n_dist"]; nth = P["amo_full"]["n_theta"]; nT = len(trial_thetas)
    data_stack = A2.stack_data([P["data_obs"]] * P1_CHUNK)
    EVc = np.repeat(P["EV"][None], P1_CHUNK, 0)
    jlnl = np.full(nT, -np.inf); nreg = np.zeros(nT, int); nreg_true = np.zeros(nT, int)
    t0 = time.time()
    for c0 in range(0, nT, P1_CHUNK):
        idx = list(range(c0, min(c0 + P1_CHUNK, nT)))
        th = np.zeros((P1_CHUNK, nth)); th[:len(idx)] = np.stack([trial_thetas[i] for i in idx])
        LNL = A2.scan_all(P["vmap2_full"], th, data_stack, EVc, P["AI"])   # (P1_CHUNK,116,B)
        best_theta = th.copy()
        for j, i in enumerate(idx):
            post = fringe_posterior(P, LNL[j], None, P["priors"]["lit"], 1.0)
            nreg[i] = int(np.sum(post["map_qmax"] > 0.5))
            nreg_true[i] = int(np.sum((post["map_qmax"] > 0.5) & (post["map_k"] == 0)))
            best_theta[j, :ndist] = post["map_evalL"]
        jl = np.asarray(fl_b(jnp.asarray(best_theta), P["data_obs"]))
        for j, i in enumerate(idx):
            jlnl[i] = jl[j]
        if verbose and (c0 // P1_CHUNK) % 5 == 0:
            print(f"    ...{c0+len(idx)}/{nT} trials ({time.time()-t0:.0f}s)", flush=True)
    return jlnl, nreg, nreg_true


def mode_p1(config):
    """P1 registration map: scan loud0's sky around truth (others at truth), distances
    profiled; is there a unique JOINT-lnL / true-registration needle at truth, and how
    wide vs the 1-rad pulsar-term-phase prediction?"""
    print(f"=== P1 registration map (config={config}) ===", flush=True)
    P = build_problem(config)
    tt = P["theta_truth"].copy(); ndist = P["n_dist"]
    cg0, gp0 = tt[ndist + 0], tt[ndist + 1]                      # loud0 truth sky
    lt = loud_truth(P)
    print(f"  loud0 truth: cos_gwtheta={cg0:.4f} gwphi={gp0:.4f}", flush=True)

    def run_grid(half_deg, n, tag):
        # square grid in (cos_gwtheta, gwphi) spanning +/- half_deg about truth
        dr = np.radians(half_deg)
        cgs = np.clip(cg0 + np.linspace(-dr, dr, n), -1, 1)
        gps = gp0 + np.linspace(-dr, dr, n)
        trials, coords = [], []
        for cg in cgs:
            for gp in gps:
                th = tt.copy(); th[ndist + 0] = cg; th[ndist + 1] = gp
                trials.append(th); coords.append((cg, gp))
        jlnl, nreg, nreg_true = p1_eval_trials(P, trials, verbose=True)
        jlnl = jlnl.reshape(n, n); nreg_true = nreg_true.reshape(n, n)
        i0 = np.unravel_index(np.argmax(jlnl), jlnl.shape)
        best_cg, best_gp = cgs[i0[0]], gps[i0[1]]
        ang_off = np.degrees(np.arccos(np.clip(np.dot(
            _unit(best_cg, best_gp), _unit(cg0, gp0)), -1, 1)))
        print(f"\n  [{tag}] grid {n}x{n} +/-{half_deg} deg: joint lnL max={jlnl.max():.1f} "
              f"at ang {ang_off:.4f} deg from truth; max true-reg={nreg_true.max()}/116 "
              f"(at truth-cell {nreg_true[n//2,n//2]})", flush=True)
        # 1-D profile through truth row for width
        row = jlnl[:, n // 2]; row = row - row.max()
        half = np.where(row > -0.5)[0]
        if len(half) > 1:
            w = np.degrees(2 * dr) * (half[-1] - half[0]) / (n - 1)
            print(f"  [{tag}] joint-lnL FWHM(0.5nat) along cos_gwtheta ~ {w:.4f} deg", flush=True)
        return best_cg, best_gp, jlnl, nreg_true

    print("  (first grid pays a one-time scan_all compile for the D=4 batched shape)", flush=True)
    run_grid(20.0, 13, "coarse")
    run_grid(2.0, 13, "zoom1")
    run_grid(0.05, 13, "zoom2")   # resolve the ~0.006 deg needle
    return 0


def mode_p0(config):
    """P0 miscalibration anatomy: at perturbed sources, is each pulsar's fringe posterior
    CONFIDENT-WRONG (concentrated on a shifted fringe) or diffuse? Report q_true = q(n_true)
    vs q_max = max_n q(n) and the MAP fringe offset, for the ceiling pulsars + the array."""
    print(f"=== P0 miscalibration anatomy (config={config}) ===", flush=True)
    P = build_problem(config)
    prior_col = P["priors"]["lit"]
    phi_truth = P["theta_truth"][P["n_dist"]:].copy()
    scale = phi_scale(P); names = P["names"]
    n_loud = CONFIGS[config].get("population", (P["ncw"],))[0]
    loud_sl = slice(0, 8 * n_loud)
    rng = np.random.default_rng(1234)
    for delta in [0.0, 1e-3, 3e-3, 1e-2]:
        phi = phi_truth.copy()
        phi[loud_sl] += rng.normal(size=8 * n_loud) * scale[loud_sl] * delta
        LNL = scan_current(P, np.concatenate([P["L0"], phi]))
        post = fringe_posterior(P, LNL, None, prior_col, 1.0)
        Pt, _, _ = p_true_at(P, LNL, prior_col, beta=1.0)   # q_true per pulsar
        qmax = post["map_qmax"]; kmap = post["map_k"]
        # array-level: confident (qmax>0.9) AND wrong (map fringe != true k=0)
        conf_wrong = int(np.sum((qmax > 0.9) & (kmap != 0)))
        conf_right = int(np.sum((qmax > 0.9) & (kmap == 0)))
        print(f"\n  delta={delta:.0e}: array qmax>0.9 -> {conf_right} on TRUE fringe, "
              f"{conf_wrong} CONFIDENT-WRONG (shifted); diffuse(qmax<0.5)="
              f"{int(np.sum(qmax<0.5))}", flush=True)
        print(f"    {'pulsar':12s} {'q_true':>7s} {'q_max':>7s} {'map_k':>6s} {'verdict':>16s}")
        for nm in CEILING:
            a = names.index(nm)
            v = "TRUE-confident" if (qmax[a] > 0.9 and kmap[a] == 0) else \
                ("CONFIDENT-WRONG" if qmax[a] > 0.9 else "diffuse")
            print(f"    {nm:12s} {Pt[a]:7.3f} {qmax[a]:7.3f} {kmap[a]:6d} {v:>16s}", flush=True)
    return 0


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["validate-estep", "split-cost", "asimov", "debug", "seed-scan", "perturb", "p0", "p1"])
    ap.add_argument("--config", default="pop", choices=list(CONFIGS))
    ap.add_argument("--n_start", type=int, default=8)
    ap.add_argument("--seed", type=int, default=7000)
    ap.add_argument("--steps", type=int, default=40)
    ap.add_argument("--iters", type=int, default=3)
    ap.add_argument("--lr", type=float, default=1e-2)
    ap.add_argument("--thresh", type=float, default=25.0)
    ap.add_argument("--no-seed", action="store_true")
    ap.add_argument("--ladder", default="0.05,0.1,0.2,0.4,0.7,1.0")
    ap.add_argument("--no-save", action="store_true")
    a = ap.parse_args()
    ladder = [float(x) for x in a.ladder.split(",")]
    if a.mode == "validate-estep":
        sys.exit(mode_validate_estep())
    elif a.mode == "split-cost":
        sys.exit(mode_split_cost())
    elif a.mode == "debug":
        sys.exit(mode_debug(a.config, a.seed, a.steps, a.iters, a.lr, ladder))
    elif a.mode == "seed-scan":
        sys.exit(mode_seed_scan(a.config, a.thresh))
    elif a.mode == "perturb":
        sys.exit(mode_perturb(a.config))
    elif a.mode == "p0":
        sys.exit(mode_p0(a.config))
    elif a.mode == "p1":
        sys.exit(mode_p1(a.config))
    else:
        sys.exit(mode_asimov(a.config, a.n_start, save=not a.no_save, seed=not a.no_seed))
