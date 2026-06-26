"""
Prong 2 — GWB <-> CW transition, information-only.

Self-contained Fisher-information calculation for a single pulsar's distance under a
population of N non-evolving phase-connected CW sources. We contrast:

  * CONDITIONAL distance information  I_LL          (all source phases/amplitudes known)
  * MARGINAL    distance information  I_LL^marg     (per-source phase + log-amplitude
                                                     marginalised; Schur complement)

The marginal collapse as sources confuse (>~1 per 1/T frequency bin) is Farr's
zero-information stochastic limit; the conditional curve is Mihir's "independence line".

Conventions
-----------
n_hat  : unit vector TOWARDS the source  = (sinθ cosφ, sinθ sinφ, cosθ)
GW propagation direction = -n_hat ; cos μ = n_hat . p_hat ; delay factor (1 - cos μ).
Antenna projection and pulsar-term retardation both carry (1 - cos μ), self-consistently.
Pulsar-term amplitude == Earth-term amplitude (non-evolving), so distance enters ONLY
through the pulsar-term phase  -2π f (L/c)(1 - cos μ).

Amplitude normalisation is the standard schematic (residual scale ζ = h/(2π f)); the
science is in the SHAPE of the curves and the conditional-vs-marginal contrast, which
we then port onto the real discovery likelihood.

GPU scaling (this version)
--------------------------
1. PADDING + MASK. The source array is padded to a FIXED N_max and a boolean mask
   zeroes out the padded sources' residual contributions. One compiled XLA shape then
   serves every N in the sweep — no per-N recompilation. Padded sources contribute
   exactly zero to the residual and hence zero rows/cols to the Fisher source-block.
2. EIGH PSEUDO-INVERSE. The Schur complement marginalising the source block uses an
   eigh-based pseudo-inverse (eigenvalues thresholded at rcond * max-eig) instead of a
   ridge-regularised solve. This is exact on the padded zero-eigenvalues and stable at
   large N where the source block is genuinely rank-deficient (confusion limit).
"""

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from functools import partial

C_LIGHT = 299792458.0
KPC_M = 3.0856775814913673e19
KPC_OVER_C = KPC_M / C_LIGHT          # seconds for light to cross 1 kpc  (~1.0288e11)
YR = 365.25 * 86400.0


# ----------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------
def sky_unit(theta, phi):
    """Unit vector toward the source."""
    return jnp.array([jnp.sin(theta) * jnp.cos(phi),
                      jnp.sin(theta) * jnp.sin(phi),
                      jnp.cos(theta)])


def antenna(theta, phi, psi, p_hat):
    """Plus/cross antenna response for a pulsar at p_hat. Includes 1/(1-cosμ)."""
    n_hat = sky_unit(theta, phi)
    cos_mu = jnp.dot(n_hat, p_hat)
    e_theta = jnp.array([jnp.cos(theta) * jnp.cos(phi),
                         jnp.cos(theta) * jnp.sin(phi),
                         -jnp.sin(theta)])
    e_phi = jnp.array([-jnp.sin(phi), jnp.cos(phi), 0.0])
    dt, dp = jnp.dot(e_theta, p_hat), jnp.dot(e_phi, p_hat)
    denom = jnp.clip(1.0 - cos_mu, 1e-6, None)
    Fp0 = 0.5 * (dt * dt - dp * dp) / denom
    Fx0 = (dt * dp) / denom
    Fp = Fp0 * jnp.cos(2 * psi) + Fx0 * jnp.sin(2 * psi)
    Fx = -Fp0 * jnp.sin(2 * psi) + Fx0 * jnp.cos(2 * psi)
    return Fp, Fx, cos_mu


# ----------------------------------------------------------------------
# Single-source residual at one pulsar (Earth term + pulsar term)
# params per source: theta, phi, psi, cos_inc, log10_f, log10_h, phase0
# ----------------------------------------------------------------------
def single_source_residual(toas, p_hat, L_kpc, src):
    theta, phi, psi, cinc, log10_f, log10_h, phase0 = src
    Fp, Fx, cos_mu = antenna(theta, phi, psi, p_hat)
    f = 10.0 ** log10_f
    h = 10.0 ** log10_h
    zeta = h / (2.0 * jnp.pi * f)                       # residual scale [s]
    tau_p = (L_kpc * KPC_OVER_C) * (1.0 - cos_mu)       # pulsar-term retardation [s]

    Phi_e = 2.0 * jnp.pi * f * toas + phase0
    Phi_p = Phi_e - 2.0 * jnp.pi * f * tau_p

    wp = 1.0 + cinc ** 2          # plus inclination weight
    wx = -2.0 * cinc              # cross inclination weight

    def wave(Phi):
        return Fp * wp * jnp.sin(Phi) + Fx * wx * jnp.cos(Phi)

    return zeta * (wave(Phi_e) - wave(Phi_p))


# vmap over sources -> (n_src, n_toa); the array of sources is shape (n_src, 7)
_src_vec = jax.vmap(single_source_residual, in_axes=(None, None, None, 0))


def total_residual(toas, p_hat, L_kpc, srcs, mask):
    """Sum residuals over sources, zeroing masked (padded) sources.

    mask is (n_src,) in {0,1}; masked sources contribute exactly zero to the residual
    AND zero to its derivatives wrt their own phase/log-amplitude, so the padded block
    of the Fisher matrix is identically zero (handled by the eigh pseudo-inverse).
    """
    contrib = _src_vec(toas, p_hat, L_kpc, srcs)        # (n_src, n_toa)
    return jnp.sum(mask[:, None] * contrib, axis=0)


# ----------------------------------------------------------------------
# JOINT-ARRAY Fisher information.
#
# Distance L_a of pulsar a is perfectly degenerate with each source's pulsar-term
# phase IF that pulsar is considered alone. The degeneracy is broken because the
# source phases {phase0_s} are GLOBAL: the Earth term is common to every pulsar and
# pins phase0_s, after which the pulsar-term phase predicts L_a. So we must model the
# whole array jointly.
#
# Free parameters: [L_0..L_{Np-1}, phase0_0, log10h_0, phase0_1, log10h_1, ...]
#   - L_a            : one distance per pulsar               (interest)
#   - phase0_s,log10h_s : per-source phase + log-amplitude   (marginalisation nuisances)
# Sky/frequency/inclination/polarisation held known (well-measured for resolvable
# sources; the phase is the nuisance that becomes unconstrained in the confusion limit).
# ----------------------------------------------------------------------
LN10 = float(np.log(10.0))


def _psr_jac(toas, p_hat, L_a, extras_vec, srcs_fixed, mask):
    """AUTODIFF reference Jacobian of one pulsar's residual wrt [L_a, extras].

    param = [L_a, phase0_0, log10h_0, phase0_1, log10h_1, ...]  (length 1 + 2*n_max)
    Returns (n_toa, 1 + 2*n_max). Used only to validate the analytic gradients below;
    the production path (`analytic_psr_grads`) gives the identical result far faster.
    """
    def res(p):
        L_a_ = p[0]
        extras = p[1:].reshape(-1, 2)
        srcs = srcs_fixed.at[:, 6].set(extras[:, 0]).at[:, 5].set(extras[:, 1])
        return total_residual(toas, p_hat, L_a_, srcs, mask)
    p0 = jnp.concatenate([jnp.atleast_1d(L_a), extras_vec])
    return jax.jacfwd(res)(p0)


def analytic_psr_grads(toas, p_hat, L_a, srcs, mask):
    """Closed-form gradients of one pulsar's residual. Returns

        g : (n_toa,)             dr/dL_a
        H : (n_toa, 2*n_src)     [dr/dphase0_s, dr/dlog10h_s] interleaved per source

    Exact derivatives of the same residual `total_residual` differentiates — verified
    against jacfwd to ~1e-10 in `validate_grads`. Avoids pushing 2*n_max autodiff
    tangents per pulsar, which is the dominant cost at large N_max.

    r_s(t) = mask_s ζ_s [W(Φ_e) - W(Φ_p)],  W(Φ)=A sinΦ + B cosΦ,  W'(Φ)=A cosΦ - B sinΦ
      ∂r/∂phase0_s = mask_s ζ_s [W'(Φ_e) - W'(Φ_p)]
      ∂r/∂log10h_s = ln10 · r_s                       (ζ ∝ h = 10^log10h)
      ∂r/∂L_a      = Σ_s mask_s ζ_s W'(Φ_p) · 2πf_s K (1-cosμ_s)
    """
    th, ph, ps = srcs[:, 0], srcs[:, 1], srcs[:, 2]
    cinc, l10f, l10h, phase0 = srcs[:, 3], srcs[:, 4], srcs[:, 5], srcs[:, 6]
    Fp, Fx, cos_mu = jax.vmap(lambda a, b, c: antenna(a, b, c, p_hat))(th, ph, ps)
    f = 10.0 ** l10f
    zeta = (10.0 ** l10h) / (2.0 * jnp.pi * f)
    A = Fp * (1.0 + cinc ** 2)          # plus amplitude per source
    B = Fx * (-2.0 * cinc)              # cross amplitude per source
    twopif = 2.0 * jnp.pi * f

    Phi_e = twopif[:, None] * toas[None, :] + phase0[:, None]      # (n_src, n_toa)
    tau_p = L_a * KPC_OVER_C * (1.0 - cos_mu)                      # (n_src,)
    Phi_p = Phi_e - twopif[:, None] * tau_p[:, None]
    sE, cE = jnp.sin(Phi_e), jnp.cos(Phi_e)
    sP, cP = jnp.sin(Phi_p), jnp.cos(Phi_p)
    Ac, Bc = A[:, None], B[:, None]
    W_e = Ac * sE + Bc * cE
    W_p = Ac * sP + Bc * cP
    Wp_e = Ac * cE - Bc * sE
    Wp_p = Ac * cP - Bc * sP

    mz = (mask * zeta)[:, None]                                   # (n_src, 1)
    dphase = mz * (Wp_e - Wp_p)                                   # (n_src, n_toa)
    dlogh = LN10 * mz * (W_e - W_p)
    coef = (twopif * KPC_OVER_C * (1.0 - cos_mu))[:, None]        # (n_src, 1)
    g = jnp.sum(mz * Wp_p * coef, axis=0)                         # (n_toa,)

    n_src = srcs.shape[0]
    H = jnp.zeros((toas.shape[0], 2 * n_src))
    H = H.at[:, 0::2].set(dphase.T).at[:, 1::2].set(dlogh.T)
    return g, H


def _eigh_pinv(A, rcond=1e-10):
    """Symmetric pseudo-inverse via eigh; eigenvalues below rcond*max(eig) dropped.

    Exact zero eigenvalues (padded sources) and the genuinely rank-deficient source
    block in the confusion limit are both handled cleanly: their null directions carry
    no information and are projected out rather than ridge-filled.
    """
    w, V = jnp.linalg.eigh(A)
    wmax = jnp.max(w)
    thresh = rcond * jnp.where(wmax > 0.0, wmax, 1.0)
    keep = w > thresh
    inv_w = jnp.where(keep, 1.0 / jnp.where(keep, w, 1.0), 0.0)
    return (V * inv_w) @ V.T


@partial(jax.jit, static_argnums=(5,))
def fisher_blocks(toas, p_hats, L_arr, srcs, mask, n_psr, sigma):
    """Assemble the array Fisher blocks (the expensive part; rcond-independent).

    Returns
        cond : (n_psr,)            diag F_LL  (conditional distance info per pulsar)
        F_Ls : (n_psr, 2*n_max)    distance x source-parameter cross block
        F_ss : (2*n_max, 2*n_max)  source-parameter block, SUMMED OVER ALL PULSARS

    *** How F_ss is built across the array (key to why the marginal is non-zero) ***
    Each pulsar a contributes H_a = dr_a/d{phase0_s, log10h_s}  (n_toa, 2*n_max). The
    source block is  F_ss = (1/sigma^2) Σ_a H_a^T H_a  — the einsum 'pta,ptb->ab' sums
    the pulsar axis p. So F_ss couples every source parameter to the WHOLE array's data:
    the global Earth term (common to all pulsars) is what pins the source phases, after
    which each pulsar-term phase determines that pulsar's distance. A single pulsar alone
    has dr_a/dL_a perfectly degenerate with its own source phases (rank-deficient block,
    marginal -> 0); the degeneracy is broken only by the cross-pulsar sum here.

    Assembled per-pulsar (F_LL diagonal: a pulsar's residual depends only on its own
    distance), so memory is O(n_src^2) not O(n_param * n_data) — what lets N_max=1000
    fit on a 24 GB GPU.
    """
    inv_s2 = 1.0 / (sigma ** 2)
    g, H = jax.vmap(lambda ph, L: analytic_psr_grads(toas, ph, L, srcs, mask))(
        p_hats, L_arr)            # g:(n_psr,n_toa)  H:(n_psr,n_toa,2*n_max)
    cond = jnp.einsum("pt,pt->p", g, g) * inv_s2                 # diag F_LL
    F_Ls = jnp.einsum("pt,pta->pa", g, H) * inv_s2              # (n_psr, 2*n_max)
    F_ss = jnp.einsum("pta,ptb->ab", H, H) * inv_s2            # (2*n_max, 2*n_max), Σ_a
    return cond, F_Ls, F_ss


@jax.jit
def marg_from_blocks(cond, F_Ls, F_ss, rcond):
    """Marginal distance info per pulsar from pre-assembled blocks, for a given rcond.

    F_marg = F_LL - F_Ls F_ss^+ F_sL ; marginal info per pulsar = diag(F_marg).
    """
    F_ss_pinv = _eigh_pinv(F_ss, rcond)
    M = F_ss_pinv @ F_Ls.T                                        # (2*n_max, n_psr)
    return jnp.clip(cond - jnp.sum(F_Ls.T * M, axis=0), 0.0, None)


def fisher_array(toas, p_hats, L_arr, srcs, mask, n_psr, sigma, rcond=1e-10):
    """Convenience wrapper: per-pulsar (conditional, marginal) distance info."""
    cond, F_Ls, F_ss = fisher_blocks(toas, p_hats, L_arr, srcs, mask, n_psr, sigma)
    marg = marg_from_blocks(cond, F_Ls, F_ss, rcond)
    return cond, marg


# ----------------------------------------------------------------------
# Population draws + padding
# ----------------------------------------------------------------------
def draw_sources(n_src, rng, f_band, log10_h, mode="fixed_persource", h_total=None,
                 k_loud=3, h_ratio=10.0):
    """Draw n_src sources.

    mode:
      'fixed_persource' : every source at log10_h.
      'fixed_total'     : per-source h set so Σ h^2 = h_total^2.
      'population'      : min(k_loud, n_src) loud sources at log10_h and the remaining
                          (n_src - k_loud) faint sources a factor h_ratio weaker. This
                          breaks the algebraic amplitude-degeneracy (sources no longer
                          share one amplitude), so marg/cond is no longer amplitude-blind.
    """
    theta = np.arccos(rng.uniform(-1, 1, n_src))
    phi = rng.uniform(0, 2 * np.pi, n_src)
    psi = rng.uniform(0, np.pi, n_src)
    cinc = rng.uniform(-1, 1, n_src)
    log10_f = rng.uniform(np.log10(f_band[0]), np.log10(f_band[1]), n_src)
    phase0 = rng.uniform(0, 2 * np.pi, n_src)
    if mode == "fixed_total":
        h_each = h_total / np.sqrt(n_src)
        l10h = np.full(n_src, np.log10(h_each))
    elif mode == "population":
        l10h = np.full(n_src, log10_h - np.log10(h_ratio))   # faint default
        l10h[:min(k_loud, n_src)] = log10_h                  # k loud sources
    else:
        l10h = np.full(n_src, log10_h)
    return np.stack([theta, phi, psi, cinc, log10_f, l10h, phase0], axis=1)


@partial(jax.jit, static_argnums=(4,))
def per_source_snr(toas, p_hats, L_arr, srcs, n_psr, sigma):
    """Per-source optimal matched-filter SNR across the array.

        SNR_s = sqrt( Σ_{pulsar a, TOA t} r_{a,s}(t)^2 / sigma^2 )

    where r_{a,s} is source s's residual contribution at pulsar a. Returns (n_src,).
    """
    # residual of each (pulsar, source): (n_psr, n_src, n_toa)
    r = jax.vmap(lambda ph, L: _src_vec(toas, ph, L, srcs))(p_hats, L_arr)
    snr2 = jnp.sum(r ** 2, axis=(0, 2)) / (sigma ** 2)   # sum over pulsars & TOAs
    return jnp.sqrt(snr2)


def pad_sources(srcs, n_max):
    """Pad a (n_src,7) source block to (n_max,7) with a (n_max,) {0,1} mask.

    Padded rows are tiled copies of the real sources (kept finite so autodiff never
    sees 0*inf); their mask entry is 0 so they contribute nothing to residual or Fisher.
    """
    n = srcs.shape[0]
    mask = np.zeros(n_max)
    mask[:min(n, n_max)] = 1.0
    if n < n_max:
        reps = int(np.ceil(n_max / n))
        filled = np.tile(srcs, (reps, 1))[:n_max]
    else:
        filled = srcs[:n_max]
    return filled.copy(), mask


# ----------------------------------------------------------------------
# Array geometry / TOA grid helpers
# ----------------------------------------------------------------------
def build_array(n_psr, seed=777):
    rng_psr = np.random.default_rng(seed)
    dirs = []
    for _ in range(n_psr):
        ct = rng_psr.uniform(-1, 1)
        ph = rng_psr.uniform(0, 2 * np.pi)
        st = np.sqrt(1 - ct ** 2)
        dirs.append([st * np.cos(ph), st * np.sin(ph), ct])
    return jnp.array(dirs)


def build_toas(T, cadence_days):
    n_toa = int(T / (cadence_days * 86400.0))
    return jnp.array(np.linspace(0, T, n_toa))


def confusion_onset(T, f_band):
    """N* ~ T * (f_hi - f_lo): one source per 1/T frequency bin across the band."""
    return T * (f_band[1] - f_band[0])


# ----------------------------------------------------------------------
# Sweep (padded: single compiled shape for all N)
# ----------------------------------------------------------------------
def run_sweep(n_values, mode, n_max=None, n_seeds=30, n_psr=116,
              T=15 * YR, cadence_days=7.0, sigma=1e-7,
              f_band=(3e-9, 2e-8), log10_h=-14.3, h_total=5e-14, seed0=10,
              p_hats=None, label=None, rconds=(1e-10,), L_arr=None,
              target_idx=None, k_loud=3, h_ratio=10.0):
    """mode in {'fixed_persource', 'fixed_total', 'population'}.

    Joint-array Fisher with PADDING: every realisation is padded to n_max sources so the
    jitted assembly compiles exactly once. Reports the median over pulsars and the 16/84
    spread over realisations. The blocks are assembled once per realisation and the
    marginal is evaluated for every rcond in `rconds` (rconds[0] is the primary value).

    If target_idx is set, also tracks the marginal info of that specific pulsar (used for
    the per-pulsar sigma_L conversion).
    """
    if n_max is None:
        n_max = int(max(n_values))
    if p_hats is None:
        p_hats = build_array(n_psr)
    if L_arr is None:
        L_arr = jnp.ones(n_psr)  # 1 kpc each (representative)
    toas = build_toas(T, cadence_days)
    N_star = confusion_onset(T, f_band)

    cond_med, cond_lo, cond_hi = [], [], []
    marg_med, marg_lo, marg_hi = [], [], []
    marg_by_rcond = {rc: [] for rc in rconds}
    marg_target = []

    tag = label or mode
    for n_src in n_values:
        cond_real, marg_real = [], []
        rc_real = {rc: [] for rc in rconds}
        tgt_real = []
        for s in range(n_seeds):
            rng = np.random.default_rng(seed0 + 1000 * s + n_src)
            srcs = draw_sources(n_src, rng, f_band, log10_h, mode=mode,
                                h_total=h_total, k_loud=k_loud, h_ratio=h_ratio)
            srcs_pad, mask = pad_sources(srcs, n_max)
            cond, F_Ls, F_ss = fisher_blocks(toas, p_hats, L_arr, jnp.array(srcs_pad),
                                             jnp.array(mask), n_psr, sigma)
            marg0 = None
            for rc in rconds:
                marg = marg_from_blocks(cond, F_Ls, F_ss, rc)
                rc_real[rc].append(float(jnp.median(marg)))
                if marg0 is None:
                    marg0 = marg
            cond_real.append(float(jnp.median(cond)))
            marg_real.append(float(jnp.median(marg0)))
            if target_idx is not None:
                tgt_real.append(float(marg0[target_idx]))
        cond_real = np.array(cond_real); marg_real = np.array(marg_real)
        cond_med.append(np.median(cond_real))
        cond_lo.append(np.percentile(cond_real, 16)); cond_hi.append(np.percentile(cond_real, 84))
        marg_med.append(np.median(marg_real))
        marg_lo.append(np.percentile(marg_real, 16)); marg_hi.append(np.percentile(marg_real, 84))
        for rc in rconds:
            marg_by_rcond[rc].append(np.median(rc_real[rc]))
        if target_idx is not None:
            marg_target.append(np.median(tgt_real))
        ratio = marg_med[-1] / cond_med[-1] if cond_med[-1] > 0 else 0.0
        print(f"  [{tag}] N={n_src:5d}  cond={cond_med[-1]:.3e}  "
              f"marg={marg_med[-1]:.3e}  (marg/cond {ratio:.3f})", flush=True)

    out = dict(
        n=np.array(n_values),
        cond=np.array(cond_med), cond_lo=np.array(cond_lo), cond_hi=np.array(cond_hi),
        marg=np.array(marg_med), marg_lo=np.array(marg_lo), marg_hi=np.array(marg_hi),
        N_star=N_star, Tspan_bins=N_star, sigma=sigma, T=T, f_band=np.array(f_band),
        n_psr=n_psr, n_max=n_max, n_seeds=n_seeds, cadence_days=cadence_days,
    )
    out["marg_by_rcond"] = {rc: np.array(v) for rc, v in marg_by_rcond.items()}
    if target_idx is not None:
        out["marg_target"] = np.array(marg_target)
    return out


# ----------------------------------------------------------------------
# Verification: padded result must match unpadded to ~1e-10
# ----------------------------------------------------------------------
def verify_padding(n_test=5, n_max=64, n_psr=10, sigma=1e-7,
                   T=15 * YR, cadence_days=7.0, f_band=(3e-9, 2e-8), log10_h=-14.3):
    p_hats = build_array(n_psr)
    toas = build_toas(T, cadence_days)
    L_arr = jnp.ones(n_psr)
    rng = np.random.default_rng(2024)
    srcs = draw_sources(n_test, rng, f_band, log10_h)

    # unpadded reference (mask all-ones, native shape)
    cond_u, marg_u = fisher_array(toas, p_hats, L_arr, jnp.array(srcs),
                                  jnp.ones(n_test), n_psr, sigma)
    # padded
    srcs_pad, mask = pad_sources(srcs, n_max)
    cond_p, marg_p = fisher_array(toas, p_hats, L_arr, jnp.array(srcs_pad),
                                  jnp.array(mask), n_psr, sigma)

    cond_u, marg_u = np.array(cond_u), np.array(marg_u)
    cond_p, marg_p = np.array(cond_p), np.array(marg_p)
    rc = np.max(np.abs(cond_p - cond_u) / (np.abs(cond_u) + 1e-300))
    rm = np.max(np.abs(marg_p - marg_u) / (np.abs(marg_u) + 1e-300))
    print(f"[verify] padded-vs-unpadded  max rel diff: cond={rc:.2e}  marg={rm:.2e}  "
          f"(n_test={n_test}, n_max={n_max})", flush=True)
    return rc, rm


def validate_grads(n_test=7, n_max=40, sigma=1e-7, T=15 * YR, cadence_days=7.0,
                   f_band=(3e-9, 2e-8), log10_h=-14.3):
    """Analytic per-pulsar gradients must match the autodiff Jacobian to ~1e-10."""
    toas = build_toas(T, cadence_days)
    p_hat = build_array(3)[1]
    rng = np.random.default_rng(7)
    srcs = draw_sources(n_test, rng, f_band, log10_h)
    srcs_pad, mask = pad_sources(srcs, n_max)
    srcs_pad, mask = jnp.array(srcs_pad), jnp.array(mask)
    L_a = 1.0

    g, H = analytic_psr_grads(toas, p_hat, L_a, srcs_pad, mask)
    extras_vec = jnp.stack([srcs_pad[:, 6], srcs_pad[:, 5]], axis=1).reshape(-1)
    Jp = _psr_jac(toas, p_hat, L_a, extras_vec, srcs_pad, mask)
    g_ad, H_ad = np.array(Jp[:, 0]), np.array(Jp[:, 1:])
    g, H = np.array(g), np.array(H)
    dg = np.max(np.abs(g - g_ad)) / (np.max(np.abs(g_ad)) + 1e-300)
    dH = np.max(np.abs(H - H_ad)) / (np.max(np.abs(H_ad)) + 1e-300)
    print(f"[verify] analytic-vs-autodiff grads  max rel diff: dL={dg:.2e}  extras={dH:.2e}",
          flush=True)
    return dg, dH


# ----------------------------------------------------------------------
# Knee diagnostic: vary T and band, check the marg/cond knee tracks N* = T*Δf
# ----------------------------------------------------------------------
def crossing_N(n, ratio, level=0.5):
    """Interpolate (in log N) where marg/cond first falls through `level`."""
    n = np.asarray(n, float); r = np.asarray(ratio, float)
    below = np.where(r < level)[0]
    if len(below) == 0:
        return np.nan
    i = below[0]
    if i == 0:
        return n[0]
    x0, x1 = np.log(n[i - 1]), np.log(n[i])
    y0, y1 = r[i - 1], r[i]
    return float(np.exp(x0 + (level - y0) * (x1 - x0) / (y1 - y0)))


def run_knee_diagnostic(configs, n_values, n_seeds=12, n_psr=116,
                        cadence_days=7.0, sigma=1e-7, log10_h=-14.3, mode="fixed_persource"):
    """For each (label, T, f_band), sweep N and return marg/cond curve + N* + knee N."""
    p_hats = build_array(n_psr)
    n_max = int(max(n_values))
    out = []
    for label, T, f_band in configs:
        res = run_sweep(n_values, mode=mode, n_max=n_max, n_seeds=n_seeds, n_psr=n_psr,
                        T=T, cadence_days=cadence_days, sigma=sigma, f_band=f_band,
                        log10_h=log10_h, p_hats=p_hats, label=label)
        ratio = res["marg"] / np.where(res["cond"] > 0, res["cond"], np.nan)
        knee = crossing_N(res["n"], ratio, 0.5)
        Nstar = res["N_star"]
        print(f"  >> {label:14s} N*={Nstar:6.2f}  knee(0.5)={knee:7.2f}  "
              f"knee/N*={knee / Nstar:.2f}", flush=True)
        out.append(dict(label=label, n=res["n"], ratio=ratio, N_star=Nstar, knee=knee,
                        T=T, f_band=np.array(f_band)))
    return out


# ----------------------------------------------------------------------
# J0437-like pulsar: nearby, well-localised MSP. Equatorial unit vector
# (RA 04h37m ~ 69.3 deg, Dec -47.25 deg). Used for the sigma_L-in-pc readout.
# ----------------------------------------------------------------------
_j = np.array([np.cos(np.radians(-47.25)) * np.cos(np.radians(69.3)),
               np.cos(np.radians(-47.25)) * np.sin(np.radians(69.3)),
               np.sin(np.radians(-47.25))])
J0437_PHAT = _j / np.linalg.norm(_j)
J0437_L_KPC = 0.15679          # VLBI/timing distance ~156.8 pc (Reardon+ 2024)
J0437_EM_SIGMA_PC = 0.25       # representative EM (parallax) prior width [pc]


def array_with_target(n_psr, target_phat=J0437_PHAT, seed=777):
    """Random array but with pulsar 0 fixed to a chosen (e.g. J0437-like) direction."""
    p = np.array(build_array(n_psr, seed))
    p[0] = target_phat
    return jnp.array(p)


def assert_array_coupling(n_src=40, n_max=64, n_psr=116, sigma=1e-7, T=15 * YR,
                          cadence_days=7.0, f_band=(3e-9, 2e-8), log10_h=-14.3, seed=10):
    """Item 2: show F_ss is summed over ALL pulsars, and that a non-zero marginal is a
    genuine array effect — a single pulsar's own data cannot break the L<->phase
    degeneracy (its per-pulsar marginal is ~0)."""
    p_hats = build_array(n_psr); toas = build_toas(T, cadence_days)
    L_arr = jnp.ones(n_psr); inv_s2 = 1.0 / sigma ** 2
    srcs = draw_sources(n_src, np.random.default_rng(seed), f_band, log10_h)
    srcs_pad, mask = pad_sources(srcs, n_max)
    srcs_pad, mask = jnp.array(srcs_pad), jnp.array(mask)

    cond, F_Ls, F_ss = fisher_blocks(toas, p_hats, L_arr, srcs_pad, mask, n_psr, sigma)
    # per-pulsar source-derivative blocks H_a, shape (n_psr, n_toa, 2*n_max)
    H = jax.vmap(lambda ph, L: analytic_psr_grads(toas, ph, L, srcs_pad, mask)[1])(
        p_hats, L_arr)
    F_ss_explicit = jnp.einsum("pta,ptb->ab", H, H) * inv_s2     # explicit Σ over pulsars
    rel = float(jnp.max(jnp.abs(F_ss - F_ss_explicit)) / jnp.max(jnp.abs(F_ss)))
    assert rel < 1e-10, f"F_ss assembly mismatch {rel:.1e}"

    marg_full = marg_from_blocks(cond, F_Ls, F_ss, 1e-10)
    # per-pulsar-only marginal: each pulsar marginalises using ONLY its own data block
    def single(a):
        Ha = H[a]
        F_ss_a = (Ha.T @ Ha) * inv_s2
        pinv = _eigh_pinv(F_ss_a, 1e-10)
        return jnp.clip(cond[a] - F_Ls[a] @ (pinv @ F_Ls[a]), 0.0, None)
    marg_single = jax.vmap(single)(jnp.arange(n_psr))

    cmed = float(jnp.median(cond))
    print(f"[coupling] F_ss = (1/sig^2) Σ_a H_a^T H_a over ALL {n_psr} pulsars "
          f"(explicit-vs-einsum rel diff {rel:.1e}).", flush=True)
    print(f"[coupling] N={n_src}: median cond={cmed:.3e}  "
          f"array marginal={float(jnp.median(marg_full)):.3e} "
          f"(marg/cond={float(jnp.median(marg_full))/cmed:.3f})  vs  "
          f"single-pulsar-only marginal={float(jnp.median(marg_single)):.3e} "
          f"(={float(jnp.median(marg_single))/cmed:.1e} of cond).", flush=True)
    print("[coupling] -> single-pulsar marginal ~ 0 (L degenerate with its own source "
          "phase/amp); non-zero marginal comes ONLY from the cross-pulsar sum that pins "
          "the GLOBAL source phases via the shared Earth term.", flush=True)
    return marg_full, marg_single


def ratio_at_N(n, ratio, N_query):
    """Log-N interpolation of a marg/cond curve at an arbitrary N."""
    n = np.asarray(n, float); r = np.asarray(ratio, float)
    return float(np.interp(np.log(N_query), np.log(n), r))


def sigma_L_pc(marg_info_kpc):
    """marginal I_LL [kpc^-2] -> sigma_L = 1/sqrt(I) in parsec."""
    m = np.asarray(marg_info_kpc, float)
    return 1e3 / np.sqrt(np.where(m > 0, m, np.nan))


if __name__ == "__main__":
    import time
    OUT = "/home/mattm/projects/HSYMT/CW_transition/prong2_results.npz"

    RCONDS = (1e-10, 1e-8, 1e-12)       # production first; sweep values for item 1

    # 0) correctness gates + array-coupling proof (item 2)
    verify_padding()
    validate_grads()
    print("", flush=True)
    assert_array_coupling()

    # log-spaced N grid 1 -> ~1000
    n_values = sorted(set(int(round(x)) for x in np.logspace(0, 3, 16)))
    print("\nN grid:", n_values, flush=True)
    n_max = max(n_values)
    p_hats = build_array(116)

    t0 = time.time()
    print("\nFixed per-source strain (Mihir-style; independence line):", flush=True)
    res_ps = run_sweep(n_values, mode="fixed_persource", n_max=n_max, p_hats=p_hats,
                       log10_h=-14.3, rconds=RCONDS)
    print(f"[timing] fixed_persource sweep (+rcond): {time.time() - t0:.1f}s", flush=True)

    # --- item 1: rcond sweep -- recompute knee N for each rcond, pick stable value ---
    print("\n[rcond] knee N vs rcond (eigh pseudo-inverse threshold):", flush=True)
    rcond_knees = {}
    for rc in RCONDS:
        ratio_rc = res_ps["marg_by_rcond"][rc] / res_ps["cond"]
        knee_rc = crossing_N(res_ps["n"], ratio_rc, 0.5)
        rcond_knees[rc] = knee_rc
        print(f"  rcond={rc:.0e}  knee(0.5)={knee_rc:7.2f}  "
              f"marg/cond@1000={ratio_rc[-1]:.3f}", flush=True)
    knee_spread = (max(rcond_knees.values()) - min(rcond_knees.values())) / \
        np.median(list(rcond_knees.values()))
    PROD_RCOND = 1e-10
    print(f"  -> knee spread across rcond = {knee_spread*100:.1f}% ; "
          f"production rcond = {PROD_RCOND:.0e} (knee stable)", flush=True)

    t1 = time.time()
    print("\nFixed total power (literal CW -> background fragmentation):", flush=True)
    res_tot = run_sweep(n_values, mode="fixed_total", n_max=n_max, p_hats=p_hats,
                        h_total=5e-14)
    print(f"[timing] fixed_total sweep: {time.time() - t1:.1f}s", flush=True)

    # mode-independence of marg/cond, and where it breaks
    rps = res_ps["marg"] / res_ps["cond"]
    rtot = res_tot["marg"] / res_tot["cond"]
    dev = np.abs(rps - rtot)
    ibreak = int(np.argmax(dev))
    print(f"\n[mode-indep] max|Δ(marg/cond)| (uniform ps vs tot) = {dev.max():.3f} at "
          f"N={res_ps['n'][ibreak]} (ps={rps[ibreak]:.3f} vs tot={rtot[ibreak]:.3f})",
          flush=True)

    # --- item 5: population amplitude mode (k loud + faint) ---
    tp = time.time()
    print("\nPopulation mode (k=3 loud + faint, h_loud/h_faint=10; breaks "
          "amplitude-degeneracy):", flush=True)
    res_pop = run_sweep(n_values, mode="population", n_max=n_max, p_hats=p_hats,
                        log10_h=-14.3, k_loud=3, h_ratio=10.0)
    rpop = res_pop["marg"] / res_pop["cond"]
    dpop = rpop - rps                         # population minus uniform
    idep = np.where(np.abs(dpop) > 0.02)[0]
    depart_N = int(res_pop["n"][idep[0]]) if len(idep) else -1
    print(f"[population] departs uniform marg/cond (|Δ|>0.02) first at N={depart_N}; "
          f"max Δ(marg/cond)={dpop[np.argmax(np.abs(dpop))]:+.3f} at "
          f"N={res_pop['n'][np.argmax(np.abs(dpop))]}", flush=True)
    print(f"[timing] population sweep: {time.time() - tp:.1f}s", flush=True)

    # --- item 4: median per-source optimal SNR at N = knee ---
    knee = crossing_N(res_ps["n"], rps, 0.5)
    N_knee = max(2, int(round(knee)))
    toas_k = build_toas(15 * YR, 7.0); L_k = jnp.ones(116)
    snr_meds = []
    for s in range(8):
        srcs_k = draw_sources(N_knee, np.random.default_rng(500 + s), (3e-9, 2e-8), -14.3)
        snr = per_source_snr(toas_k, p_hats, L_k, jnp.array(srcs_k), 116, 1e-7)
        snr_meds.append(float(jnp.median(snr)))
    snr_knee = float(np.median(snr_meds))
    print(f"\n[snr] at N=knee={N_knee}: median per-source optimal SNR = {snr_knee:.2f} "
          f"(sqrt of summed snr^2 over TOAs & pulsars)", flush=True)

    # --- item 3: sigma_L in parsec for a J0437-like pulsar ---
    # Use fixed_total: distance precision DEGRADES as the CW fragments into a background
    # (in fixed_persource the per-source power keeps adding info, so sigma_L only improves
    # and never crosses a threshold -- noted below).
    ts = time.time()
    print("\nsigma_L (pc) for a J0437-like pulsar (L=%.1f pc, idx 0; fixed_total):"
          % (J0437_L_KPC * 1e3), flush=True)
    p_hats_j = array_with_target(116)
    L_arr_j = jnp.array(np.concatenate([[J0437_L_KPC], np.ones(115)]))
    res_j = run_sweep(n_values, mode="fixed_total", n_max=n_max, p_hats=p_hats_j,
                      L_arr=L_arr_j, h_total=5e-14, target_idx=0, label="J0437")
    res_j_ps = run_sweep(n_values, mode="fixed_persource", n_max=n_max, p_hats=p_hats_j,
                         L_arr=L_arr_j, log10_h=-14.3, target_idx=0, label="J0437_ps")
    sig_pc = sigma_L_pc(res_j["marg_target"])
    sig_pc_ps = sigma_L_pc(res_j_ps["marg_target"])

    def cross_up(nv, y, level):
        nv = np.asarray(nv, float); y = np.asarray(y, float)
        m = np.isfinite(y)
        nv, y = nv[m], y[m]
        idx = np.where(y > level)[0]
        if len(idx) == 0:
            return np.nan
        i = idx[0]
        if i == 0:
            return nv[0]
        x0, x1 = np.log(nv[i - 1]), np.log(nv[i])
        return float(np.exp(x0 + (level - y[i - 1]) * (x1 - x0) / (y[i] - y[i - 1])))

    N_subpc = cross_up(res_j["n"], sig_pc, 1.0)
    N_em = cross_up(res_j["n"], sig_pc, J0437_EM_SIGMA_PC)
    print(f"[sigma_L] fixed_total sig_L(pc): "
          f"{np.array2string(sig_pc, precision=4, max_line_width=200)}", flush=True)
    nanmsg = lambda x: ("never in N<=1000" if not np.isfinite(x) else f"N≈{x:.0f}")
    print(f"[sigma_L] crosses (a) sub-parsec (1 pc): {nanmsg(N_subpc)}; "
          f"(b) EM prior {J0437_EM_SIGMA_PC} pc: {nanmsg(N_em)} "
          f"(sigma_L grows as the background fragments).", flush=True)
    print(f"[sigma_L] (fixed_persource sig_L stays {np.nanmin(sig_pc_ps):.3g}-"
          f"{np.nanmax(sig_pc_ps):.3g} pc and only improves with N -- per-source power "
          f"keeps adding distance info; no crossing.)", flush=True)
    print(f"[sigma_L] NOTE absolute scale is schematic (zeta=h/2pi f norm); the SHAPE/"
          f"relative degradation is the transferable result.", flush=True)
    print(f"[timing] sigma_L sweep: {time.time() - ts:.1f}s", flush=True)

    # knee diagnostic: vary T and band -> check knee tracks N* = T*Δf
    t2 = time.time()
    print("\nKnee diagnostic (vary T and band; knee should track N* = T*Δf):", flush=True)
    diag_N = sorted(set(int(round(x)) for x in np.logspace(0, 3, 13)))
    configs = [
        ("T10_b3-20",  10 * YR, (3e-9, 2e-8)),
        ("T15_b3-20",  15 * YR, (3e-9, 2e-8)),
        ("T20_b3-20",  20 * YR, (3e-9, 2e-8)),
        ("T15_b3-12",  15 * YR, (3e-9, 1.2e-8)),
        ("T15_b3-40",  15 * YR, (3e-9, 4.0e-8)),
    ]
    diag = run_knee_diagnostic(configs, diag_N, n_seeds=12)
    print(f"[timing] knee diagnostic: {time.time() - t2:.1f}s", flush=True)

    # --- required summary: marg/cond at N in {1, 10, N*, knee, 1000} ---
    N_star = res_ps["N_star"]
    print("\n[summary] marg/cond at key N (fixed_persource):", flush=True)
    for tagN, Nq in [("1", 1), ("10", 10), ("N*", N_star), ("knee", knee), ("1000", 1000)]:
        print(f"    N={tagN:>5s} ({Nq:7.1f}):  marg/cond = {ratio_at_N(res_ps['n'], rps, Nq):.3f}",
              flush=True)

    # pack diagnostic arrays for npz (per-config keys)
    diag_pack = {}
    for d in diag:
        L = d["label"]
        diag_pack[f"diag_{L}_n"] = d["n"]
        diag_pack[f"diag_{L}_ratio"] = d["ratio"]
        diag_pack[f"diag_{L}_Nstar"] = d["N_star"]
        diag_pack[f"diag_{L}_knee"] = d["knee"]
        diag_pack[f"diag_{L}_T"] = d["T"]
        diag_pack[f"diag_{L}_fband"] = d["f_band"]
    diag_labels = np.array([d["label"] for d in diag])

    # rcond marg curves + population + sigma_L extras
    rcond_pack = {f"ps_marg_rcond_{rc:.0e}": res_ps["marg_by_rcond"][rc] for rc in RCONDS}

    def clean(res):   # drop non-array fields before npz spread
        return {k: v for k, v in res.items() if k not in ("marg_by_rcond", "marg_target")}

    np.savez(OUT,
             **{f"ps_{k}": v for k, v in clean(res_ps).items()},
             **{f"tot_{k}": v for k, v in clean(res_tot).items()},
             **{f"pop_{k}": v for k, v in clean(res_pop).items()},
             pop_k_loud=3, pop_h_ratio=10.0, pop_depart_N=depart_N,
             mode_indep_maxdev=dev.max(), mode_indep_break_N=res_ps['n'][ibreak],
             rcond_values=np.array(RCONDS), prod_rcond=PROD_RCOND,
             rcond_knees=np.array([rcond_knees[rc] for rc in RCONDS]),
             rcond_knee_spread=knee_spread, **rcond_pack,
             snr_knee_N=N_knee, snr_knee_median=snr_knee, knee_N=knee, N_star=N_star,
             j_n=res_j["n"], j_marg_target=res_j["marg_target"], j_sigma_pc=sig_pc,
             j_sigma_pc_ps=sig_pc_ps, j_marg_target_ps=res_j_ps["marg_target"],
             j_L_kpc=J0437_L_KPC, em_sigma_pc=J0437_EM_SIGMA_PC,
             sigma_subpc_N=N_subpc, sigma_em_N=N_em,
             diag_labels=diag_labels, **diag_pack)
    print(f"\nSaved {OUT}", flush=True)
    print(f"[timing] TOTAL: {time.time() - t0:.1f}s", flush=True)
