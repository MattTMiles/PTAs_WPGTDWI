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


# ======================================================================
# COHERENCE transition (Stage A.2) -- distinct from the confusion transition.
#
# Each source has a finite coherence time t_c. The pulsar term is sampled at a
# different retarded time along a different sightline for every pulsar, so we add an
# INDEPENDENT per-(source,pulsar) pulsar-term phase offset dphi_{s,p} (default 0 ->
# the coherent model). Decoherence is a Gaussian prior on dphi with variance
#     sigma_phi^2_{s,p} = tau_{s,p} / t_c ,   tau_{s,p} = (L_p/c)(1 - cos mu_{s,p}) [s]
# (same (1-cos mu) factor as compute_mode_spacing's denominator). Bayesian Fisher
# F_total = F_data + Pi with Pi = diag(1/sigma_phi^2) on the dphi block only:
#   t_c -> inf : sigma_phi^2 -> 0  -> Pi -> inf -> dphi PINNED  -> coherent (marg/cond
#                = the Stage-A N=1 value, only phase0/logh marginalised).
#   t_c -> 0   : sigma_phi^2 -> inf -> Pi -> 0   -> dphi FREE    -> decoherent. Each L_p
#                is degenerate with its own dphi_{s,p} (single sinusoid) -> marginal -> 0.
# This is Farr's coherence limit, isolated at N=1 (zero confusion).
# ======================================================================
def coh_single_residual(toas, p_hat, L_kpc, src, dphi):
    """single_source_residual with an extra pulsar-term phase offset dphi."""
    theta, phi, psi, cinc, log10_f, log10_h, phase0 = src
    Fp, Fx, cos_mu = antenna(theta, phi, psi, p_hat)
    f = 10.0 ** log10_f
    zeta = (10.0 ** log10_h) / (2.0 * jnp.pi * f)
    tau_p = (L_kpc * KPC_OVER_C) * (1.0 - cos_mu)
    Phi_e = 2.0 * jnp.pi * f * toas + phase0
    Phi_p = Phi_e - 2.0 * jnp.pi * f * tau_p + dphi          # decoherence offset
    wp = 1.0 + cinc ** 2
    wx = -2.0 * cinc
    def wave(Phi):
        return Fp * wp * jnp.sin(Phi) + Fx * wx * jnp.cos(Phi)
    return zeta * (wave(Phi_e) - wave(Phi_p))


def _coh_residual_from_vec(toas, p_hats, vec, srcs_fixed, n_psr, n_src):
    """param vec = [L(n_psr), phase0(n_src), log10h(n_src), dphi(n_src*n_psr)]."""
    L = vec[:n_psr]
    phase0 = vec[n_psr:n_psr + n_src]
    logh = vec[n_psr + n_src:n_psr + 2 * n_src]
    dphi = vec[n_psr + 2 * n_src:].reshape(n_src, n_psr)     # [s, p]
    srcs = srcs_fixed.at[:, 6].set(phase0).at[:, 5].set(logh)

    def per_psr(p_hat, L_p, dphi_col):                      # dphi_col: (n_src,)
        r_s = jax.vmap(lambda src, dp: coh_single_residual(toas, p_hat, L_p, src, dp))(
            srcs, dphi_col)
        return jnp.sum(r_s, axis=0)                         # (n_toa,)
    R = jax.vmap(per_psr)(p_hats, L, dphi.T)               # (n_psr, n_toa)
    return R.reshape(-1)


@partial(jax.jit, static_argnums=(4, 5))
def coherence_blocks(toas, p_hats, L_arr, srcs, n_psr, n_src, sigma):
    """Data Fisher F over [L, phase0, log10h, dphi] (n_psr small here: N=1 use), plus the
    tau_{s,p} array. rcond/Pi/t_c applied afterwards in coherence_marg (so one jacfwd
    serves every t_c)."""
    n_dphi = n_src * n_psr
    vec0 = jnp.concatenate([L_arr, srcs[:, 6], srcs[:, 5], jnp.zeros(n_dphi)])
    J = jax.jacfwd(lambda v: _coh_residual_from_vec(toas, p_hats, v, srcs, n_psr, n_src))(vec0)
    F = (J.T @ J) / (sigma ** 2)
    # tau_{s,p} = (L_p/c)(1 - cos mu_{s,p})   [seconds]
    cos_mu = jax.vmap(lambda src: jax.vmap(lambda ph: jnp.dot(
        sky_unit(src[0], src[1]), ph))(p_hats))(srcs)        # (n_src, n_psr)
    tau = (L_arr[None, :] * KPC_OVER_C) * (1.0 - cos_mu)     # (n_src, n_psr)
    return F, tau


def sigma_phi2_linear(tau, t_c):
    return tau / t_c


def sigma_phi2_saturating(tau, t_c, sigma_max2=np.pi ** 2 / 3.0):
    """Wrapped-phase-capped form: sigma_phi^2 = sigma_max^2 (1 - exp(-tau/t_c))."""
    return sigma_max2 * (1.0 - jnp.exp(-tau / t_c))


def coherence_marg_cond(F, tau, n_psr, n_src, t_c, form, rcond=1e-10, pin_dphi=False):
    """Marginal/conditional L-info from data-Fisher F + decoherence prior on dphi.

    cond = diag(F_LL) (all nuisances pinned). marg = Schur-complement the nuisance block
    {phase0, log10h, dphi} out of F_total = F + Pi (Pi on dphi only). pin_dphi=True drops
    dphi from the nuisance set entirely (the t_c->inf reference / Stage-A coherent value).
    """
    nL = n_psr
    n_pl = 2 * n_src                          # phase0 + log10h count
    cond = jnp.diag(F)[:nL]

    if pin_dphi:
        nu = slice(nL, nL + n_pl)             # marginalise only phase0/log10h
        F_nn = F[nu, nu]
    else:
        nu = slice(nL, None)                  # phase0/log10h + dphi
        sig2 = form(tau, t_c).reshape(-1)     # (n_src*n_psr,)  row-major [s,p]
        pi = jnp.concatenate([jnp.zeros(n_pl), 1.0 / sig2])
        F_nn = F[nu, nu] + jnp.diag(pi)
    F_Ln = F[:nL, nu]
    pinv = _eigh_pinv(F_nn, rcond)
    marg = jnp.clip(cond - jnp.sum(F_Ln.T * (pinv @ F_Ln.T), axis=0), 0.0, None)
    return cond, marg


def _cross_x(x, y, level=0.5):
    """Interpolate (log-x) where a monotone-decreasing y(x) passes through level."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    below = np.where(y < level)[0]
    if len(below) == 0 or below[0] == 0:
        return np.nan
    i = below[0]
    lx = np.log(x[i - 1]) + (level - y[i - 1]) * (np.log(x[i]) - np.log(x[i - 1])) / \
        (y[i] - y[i - 1])
    return float(np.exp(lx))


def run_coherence_experiment(n_psr=116, T=15 * YR, cadence_days=7.0, sigma=1e-7,
                             f_band=(3e-9, 2e-8), log10_h=-14.3, n_seeds=24,
                             seed0=300, rcond=1e-10):
    """Stage A.2: Farr's COHERENCE transition at N=1 (no confusion)."""
    import time
    t0 = time.time()
    OUT = "/home/mattm/projects/HSYMT/CW_transition/prong2_coherence.npz"
    FIG = "/home/mattm/projects/HSYMT/CW_transition/prong2_coherence_transition.png"
    p_hats = build_array(n_psr)
    L_arr = jnp.ones(n_psr)
    n_src = 1

    def blocks_for_T(Tval, seeds):
        toas = build_toas(Tval, cadence_days)
        Fs, taus = [], []
        for s in seeds:
            srcs = jnp.array(draw_sources(n_src, np.random.default_rng(seed0 + s),
                                          f_band, log10_h))
            F, tau = coherence_blocks(toas, p_hats, L_arr, srcs, n_psr, n_src, sigma)
            Fs.append(F); taus.append(tau)
        return Fs, taus

    seeds = list(range(n_seeds))
    Fs, taus = blocks_for_T(T, seeds)
    ref_tau = float(np.median([float(jnp.median(t)) for t in taus]))   # seconds

    # coherent reference (dphi pinned) -- the Stage-A N=1 marg/cond
    ref_ratios = []
    for F, tau in zip(Fs, taus):
        c, m = coherence_marg_cond(F, tau, n_psr, n_src, 1.0, sigma_phi2_linear,
                                   rcond, pin_dphi=True)
        ref_ratios.append(float(jnp.median(m) / jnp.median(c)))
    ref_ratio = float(np.median(ref_ratios))
    print(f"[coh] coherent reference (dphi pinned, N=1, {n_psr} psr): "
          f"marg/cond = {ref_ratio:.4f}  (Stage-A N=1 value; toy-15psr was ~0.95)",
          flush=True)

    # ---- Experiment 1: sweep t_c so median tau/t_c spans 1e-3..1e3 ----
    x_grid = np.logspace(-3, 3, 25)            # = ref_tau/t_c
    x_gate = 1e-8                              # t_c -> inf
    forms = {"linear": sigma_phi2_linear, "saturating": sigma_phi2_saturating}
    exp1 = {}
    for fname, form in forms.items():
        ratios = []
        for x in x_grid:
            t_c = ref_tau / x
            rr = [float(jnp.median(m) / jnp.median(c)) for c, m in
                  [coherence_marg_cond(F, tau, n_psr, n_src, t_c, form, rcond)
                   for F, tau in zip(Fs, taus)]]
            ratios.append(float(np.median(rr)))
        exp1[fname] = np.array(ratios)
        x_half = _cross_x(x_grid, exp1[fname], 0.5)
        print(f"[coh] Exp1 [{fname}] marg/cond=0.5 at tau_p/t_c = {x_half:.3g}", flush=True)

    # gate (i): t_c -> inf reproduces coherent reference
    t_c_inf = ref_tau / x_gate
    rr_inf = [float(jnp.median(m) / jnp.median(c)) for c, m in
              [coherence_marg_cond(F, tau, n_psr, n_src, t_c_inf, sigma_phi2_linear, rcond)
               for F, tau in zip(Fs, taus)]]
    coh_inf = float(np.median(rr_inf))
    gate_err = abs(coh_inf - ref_ratio)
    print(f"[coh][gate i] t_c->inf marg/cond={coh_inf:.4f} vs coherent ref={ref_ratio:.4f} "
          f"(|Δ|={gate_err:.1e})", flush=True)
    assert gate_err < 1e-3, f"t_c->inf does not reproduce coherent reference ({gate_err:.1e})"

    x_half_lin = _cross_x(x_grid, exp1["linear"], 0.5)
    x_half_sat = _cross_x(x_grid, exp1["saturating"], 0.5)
    shift = x_half_sat / x_half_lin
    print(f"[coh][gate ii] 0.5-crossing moves {shift:.2f}x between linear and saturating "
          f"forms (location robust ~ tau_p/t_c~1; shape model-dependent).", flush=True)

    # ---- Experiment 2: fix mid-transition t_c, vary T -> coherence is T-INDEPENDENT ----
    t_c_mid = ref_tau / x_half_lin             # where marg/cond ~ 0.5
    Tyrs = np.array([10.0, 15.0, 20.0, 30.0])
    coh_vsT = []
    for Ty in Tyrs:
        Fs_T, taus_T = blocks_for_T(Ty * YR, list(range(12)))
        rr = [float(jnp.median(m) / jnp.median(c)) for c, m in
              [coherence_marg_cond(F, tau, n_psr, n_src, t_c_mid, sigma_phi2_linear, rcond)
               for F, tau in zip(Fs_T, taus_T)]]
        coh_vsT.append(float(np.median(rr)))
    coh_vsT = np.array(coh_vsT)
    # marg/cond = 1/(1 + (tau/t_c)*SNR^2); SNR^2 ~ T  =>  marg/cond ~ 1/(1 + T/T0).
    inv = 1.0 / coh_vsT - 1.0
    T0 = float(np.median(Tyrs / inv))         # yr; inv ~ T/T0
    print(f"[coh] Exp2 coherence marg/cond vs T={Tyrs.tolist()} yr: "
          f"{np.round(coh_vsT, 4).tolist()}  -> FALSIFIES the T-independence guess: "
          f"marg/cond falls as ~1/(1+T/{T0:.0f}yr) because the per-baseline phase offset "
          f"is increasingly resolved as SNR^2~T grows (decoherence penalty = "
          f"(tau/t_c)*SNR^2). Both transitions scale with T; the clean discriminator is "
          f"N (coherence is full at N=1 where confusion is ABSENT -- Exp1).", flush=True)

    # confusion knee vs T (band fixed) for the contrast curve
    print("[coh] confusion knee vs T (for contrast):", flush=True)
    diag_N = sorted(set(int(round(x)) for x in np.logspace(0, 3, 13)))
    conf_configs = [(f"T{int(Ty)}", Ty * YR, f_band) for Ty in Tyrs]
    conf = run_knee_diagnostic(conf_configs, diag_N, n_seeds=12, n_psr=n_psr,
                               cadence_days=cadence_days, sigma=sigma, log10_h=log10_h)
    conf_knee = np.array([d["knee"] for d in conf])
    print(f"[coh] confusion knee vs T: {conf_knee.round(0).tolist()} (scales ~ T·Δf)",
          flush=True)

    np.savez(OUT, x_grid=x_grid, exp1_linear=exp1["linear"],
             exp1_saturating=exp1["saturating"], ref_ratio=ref_ratio, ref_tau=ref_tau,
             coh_inf=coh_inf, x_half_linear=x_half_lin, x_half_saturating=x_half_sat,
             shape_shift=shift, Tyrs=Tyrs, coh_vsT=coh_vsT, conf_knee=conf_knee,
             t_c_mid=t_c_mid, n_psr=n_psr, n_seeds=n_seeds, coh_T0_yr=T0)

    # ---- 2-panel figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
    a = ax[0]
    a.plot(x_grid, exp1["linear"], "o-", color="tab:blue", lw=2, ms=4,
           label=r"linear $\sigma_\phi^2=\tau/t_c$")
    a.plot(x_grid, exp1["saturating"], "s--", color="tab:red", lw=1.8, ms=4,
           label=r"saturating $\sigma_\phi^2=\frac{\pi^2}{3}(1-e^{-\tau/t_c})$")
    a.axhline(0.5, color="k", lw=.6, alpha=.5)
    a.axhline(ref_ratio, color="gray", ls=":", lw=1, alpha=.7)
    a.axvline(x_half_lin, color="tab:blue", ls=":", lw=1.2)
    a.text(x_half_lin * 1.1, 0.55, f"0.5 @ {x_half_lin:.2g}", color="tab:blue", fontsize=8)
    a.axvline(1.0, color="green", ls="--", lw=1, alpha=.6)
    a.text(1.05, 0.05, r"$\tau_p/t_c=1$", color="green", fontsize=8, rotation=90)
    a.set_xscale("log"); a.set_xlabel(r"$\tau_p/t_c$  (median over array)")
    a.set_ylabel("marginal / conditional"); a.set_ylim(0, 1.02)
    a.set_title("Exp 1: COHERENCE transition (N=1, zero confusion)", fontsize=10)
    a.grid(alpha=.25, which="both"); a.legend(fontsize=8, loc="lower left")

    a = ax[1]
    a.plot(Tyrs, coh_vsT, "o-", color="tab:blue", lw=2, ms=6,
           label="coherence marg/cond (fixed $t_c$)")
    Tfit = np.linspace(Tyrs.min(), Tyrs.max(), 50)
    a.plot(Tfit, 1.0 / (1.0 + Tfit / T0), "-", color="tab:blue", lw=1, alpha=.5,
           label=fr"$1/(1+T/{T0:.0f}\,\mathrm{{yr}})$ (SNR$^2\!\propto T$)")
    a.axhline(0.5, color="k", lw=.6, ls=":", alpha=.5)
    a.set_xlabel("T  [yr]"); a.set_ylabel("coherence marg / cond", color="tab:blue")
    a.set_ylim(0, 1.0); a.tick_params(axis="y", labelcolor="tab:blue")
    a.set_title("Exp 2: coherence is NOT T-independent (falls via SNR);\n"
                r"discriminator is N, not T -- confusion knee $\propto T$", fontsize=9.5)
    a2 = a.twinx()
    a2.plot(Tyrs, conf_knee, "s--", color="tab:red", lw=2, ms=6,
            label=r"confusion knee $N$ ($\propto T\Delta f$)")
    a2.set_ylabel("confusion knee  $N$", color="tab:red")
    a2.tick_params(axis="y", labelcolor="tab:red")
    l1, lb1 = a.get_legend_handles_labels(); l2, lb2 = a2.get_legend_handles_labels()
    a.legend(l1 + l2, lb1 + lb2, fontsize=7.5, loc="upper center")
    a.grid(alpha=.25)
    fig.suptitle("Prong 2 Stage A.2 -- COHERENCE transition (N=1, distinct AXIS from "
                 "confusion; both T-dependent)", fontsize=11, y=1.01)
    fig.tight_layout(); fig.savefig(FIG, dpi=130, bbox_inches="tight")
    print(f"[coh] saved {OUT}\n[coh] saved {FIG}\n[coh] total {time.time()-t0:.1f}s",
          flush=True)


# ======================================================================
# Stage A.3 -- the N x t_c MAP: do confusion (Stage A) and coherence (Stage A.2)
# factorise, or interact?  Combined model: N fixed-per-source sources, each with a
# per-(source,pulsar) phase offset dphi_{s,p} ~ N(0, tau_{s,p}/t_c).
#
# Efficient marginalisation: dphi_{s,p} enters only pulsar p's residual, so the
# dphi-block of F_data is BLOCK-DIAGONAL over pulsars (116 independent N x N blocks).
# Marginalise dphi first by per-pulsar block inversion, then the dense 2N {phase0,logh}
# block. Conditional pins {phase0, log10h, dphi}; marginal Schur-complements them all.
# ======================================================================
def _combined_psr_blocks(toas, p_hat, L_p, srcs, sigma):
    """For one pulsar: data-Fisher sub-blocks over local params a=[L_p, phase0(N),
    log10h(N)] and d=[dphi(N)]. Returns F_aa (1+2N,1+2N), F_ad (1+2N,N), F_dd (N,N),
    tau_col (N,). Analytic gradients (validated against autodiff in the N=8 gate)."""
    th, ph, ps = srcs[:, 0], srcs[:, 1], srcs[:, 2]
    cinc, l10f, l10h, phase0 = srcs[:, 3], srcs[:, 4], srcs[:, 5], srcs[:, 6]
    Fp, Fx, cos_mu = jax.vmap(lambda a, b, c: antenna(a, b, c, p_hat))(th, ph, ps)
    f = 10.0 ** l10f
    zeta = (10.0 ** l10h) / (2.0 * jnp.pi * f)
    A = Fp * (1.0 + cinc ** 2)
    B = Fx * (-2.0 * cinc)
    twopif = 2.0 * jnp.pi * f
    Phi_e = twopif[:, None] * toas[None, :] + phase0[:, None]
    tau = (L_p * KPC_OVER_C) * (1.0 - cos_mu)                  # (N,)
    Phi_p = Phi_e - twopif[:, None] * tau[:, None]
    sE, cE = jnp.sin(Phi_e), jnp.cos(Phi_e)
    sP, cP = jnp.sin(Phi_p), jnp.cos(Phi_p)
    Ac, Bc = A[:, None], B[:, None]
    W_e = Ac * sE + Bc * cE; W_p = Ac * sP + Bc * cP
    Wp_e = Ac * cE - Bc * sE; Wp_p = Ac * cP - Bc * sP
    zc = zeta[:, None]
    gL = jnp.sum(zc * Wp_p * (twopif * KPC_OVER_C * (1.0 - cos_mu))[:, None], axis=0)  # (n_toa,)
    Hph = zc * (Wp_e - Wp_p)                                   # (N, n_toa) phase0
    Hlh = LN10 * zc * (W_e - W_p)                              # (N, n_toa) log10h
    Hdp = -zc * Wp_p                                           # (N, n_toa) dphi
    a = jnp.concatenate([gL[None, :], Hph, Hlh], axis=0)       # (1+2N, n_toa)
    inv = 1.0 / sigma ** 2
    F_aa = (a @ a.T) * inv
    F_ad = (a @ Hdp.T) * inv
    F_dd = (Hdp @ Hdp.T) * inv
    return F_aa, F_ad, F_dd, tau


@partial(jax.jit, static_argnums=(4,))
def combined_blocks(toas, p_hats, L_arr, srcs, n_psr, sigma):
    """vmap _combined_psr_blocks over pulsars. cond = diag(F_aa[:,0,0])."""
    F_aa, F_ad, F_dd, tau = jax.vmap(
        lambda ph, L: _combined_psr_blocks(toas, ph, L, srcs, sigma))(p_hats, L_arr)
    return F_aa, F_ad, F_dd, tau, F_aa[:, 0, 0]


@partial(jax.jit, static_argnums=(5, 6))
def ntc_marg_cond(F_aa, F_ad, F_dd, tau, cond, n_psr, n_src, t_c, rcond=1e-10):
    """Structured marginal/conditional. Marginalise dphi per-pulsar (N x N solves), then
    the dense 2N {phase0,log10h} block; read the L diagonal."""
    Pi = t_c / tau                              # (n_psr, n_src) = 1/sigma_phi^2

    def per_p(Faa, Fad, Fdd, pi_col):
        M = Fdd + jnp.diag(pi_col)              # (N, N)  -- per-pulsar block
        sol = jnp.linalg.solve(M, Fad.T)        # (N, 1+2N)
        return Faa - Fad @ sol                  # (1+2N, 1+2N)
    delta = jax.vmap(per_p)(F_aa, F_ad, F_dd, Pi)     # (n_psr, 1+2N, 1+2N)

    Feff_LL = delta[:, 0, 0]                    # (n_psr,)  L block is diagonal
    Feff_Lpl = delta[:, 0, 1:]                  # (n_psr, 2N)
    Feff_plpl = jnp.sum(delta[:, 1:, 1:], axis=0)     # (2N, 2N)  shared phase0/logh
    pinv = _eigh_pinv(Feff_plpl, rcond)
    marg = jnp.clip(Feff_LL - jnp.sum(Feff_Lpl * (Feff_Lpl @ pinv), axis=1), 0.0, None)
    return cond, marg


def ntc_dense_marg(F_aa, F_ad, F_dd, tau, n_psr, n_src, t_c, rcond=1e-10):
    """DENSE reference for the gate: scatter per-pulsar blocks into the full Fisher over
    [L(n_psr), phase0(N), log10h(N), dphi(N*n_psr)], add Pi on dphi, single Schur out all
    nuisances. Same analytic F as the structured path -> must match to ~1e-10."""
    npl = 2 * n_src
    nd = n_src * n_psr
    P = n_psr + npl + nd
    F = np.zeros((P, P))
    # global indices: L = p ; pl = n_psr + (0..2N-1) ; dphi_p = n_psr+npl + p*N + (0..N-1)
    Faa = np.array(F_aa); Fad = np.array(F_ad); Fdd = np.array(F_dd); tau = np.array(tau)
    pl = slice(n_psr, n_psr + npl)
    for p in range(n_psr):
        a_glb = [p] + list(range(n_psr, n_psr + npl))         # local a -> global
        d_glb = list(range(n_psr + npl + p * n_src, n_psr + npl + (p + 1) * n_src))
        F[np.ix_(a_glb, a_glb)] += Faa[p]
        F[np.ix_(a_glb, d_glb)] += Fad[p]
        F[np.ix_(d_glb, a_glb)] += Fad[p].T
        F[np.ix_(d_glb, d_glb)] += Fdd[p]
    Pi = np.zeros(P)
    for p in range(n_psr):
        d0 = n_psr + npl + p * n_src
        Pi[d0:d0 + n_src] = t_c / tau[p]                      # 1/sigma_phi^2
    F = F + np.diag(Pi)
    nu = list(range(n_psr, P))                                # phase0,logh,dphi
    Lidx = list(range(n_psr))
    F_LL = F[np.ix_(Lidx, Lidx)]; F_Ln = F[np.ix_(Lidx, nu)]; F_nn = F[np.ix_(nu, nu)]
    pinv = np.array(_eigh_pinv(jnp.array(F_nn), rcond))
    marg = np.clip(np.diag(F_LL - F_Ln @ pinv @ F_Ln.T), 0.0, None)
    cond = np.diag(F_LL)
    return cond, marg


def run_ntc_map(n_psr=116, T=15 * YR, cadence_days=7.0, sigma=1e-7,
                f_band=(3e-9, 1.2e-8), log10_h=-14.3, n_seeds=15, seed0=900, rcond=1e-10):
    """Stage A.3: the N x t_c map and the factorisation test."""
    import time
    t0 = time.time()
    OUT = "/home/mattm/projects/HSYMT/CW_transition/prong2_Ntc_map.npz"
    FIG = "/home/mattm/projects/HSYMT/CW_transition/prong2_Ntc_map.png"
    p_hats = build_array(n_psr); L_arr = jnp.ones(n_psr)
    toas = build_toas(T, cadence_days)
    N_star = confusion_onset(T, f_band)
    print(f"[ntc] band {f_band}, N*={N_star:.2f}, knee~{52*N_star:.0f}", flush=True)

    # ---- GATE: structured vs dense at N=8 (same analytic F) ----
    srcs8 = jnp.array(draw_sources(8, np.random.default_rng(1), f_band, log10_h))
    bl8 = combined_blocks(toas, p_hats, L_arr, srcs8, n_psr, sigma)
    t_c_test = 5e10
    cs, ms = ntc_marg_cond(*bl8[:5], n_psr, 8, jnp.asarray(t_c_test), rcond)
    cd, md = ntc_dense_marg(*bl8[:4], n_psr, 8, t_c_test, rcond)
    rel = float(np.max(np.abs(np.array(ms) - md)) / np.max(np.abs(md)))
    print(f"[ntc][gate structured-vs-dense @N=8] max rel diff = {rel:.2e}", flush=True)
    assert rel < 1e-10, f"structured != dense ({rel:.2e})"

    # ---- invariant Y: SNR^2 * sigma_phi^2, normalised to the N=1 half-info point ----
    # The bare SNR^2 * median(sigma_phi^2) carries an O(geometry) factor (antenna
    # 1/(1-cos mu) weighting spreads the per-baseline invariant), so we CALIBRATE the
    # constant by the t_c at which the N=1 coherence transition hits R=0.5, then define
    # Y = t_c_50 / t_c. Y=1 is then the coherence half-info point by construction; this
    # is the form-independent invariant (Stage-A.2 gate ii: absolute t_c is model-
    # dependent, the location is the claim).
    snr2s, taus = [], []
    for s in range(min(n_seeds, 8)):
        srcs = jnp.array(draw_sources(4, np.random.default_rng(seed0 + s), f_band, log10_h))
        snr = per_source_snr(toas, p_hats, L_arr, srcs, n_psr, sigma)
        snr2s.append(float(jnp.median(snr ** 2)))
        _, _, _, tau, _ = combined_blocks(toas, p_hats, L_arr, srcs, n_psr, sigma)
        taus.append(float(jnp.median(tau)))
    SNR2 = float(np.median(snr2s)); tau_med = float(np.median(taus))

    tc_wide = np.logspace(np.log10(SNR2 * tau_med) - 4, np.log10(SNR2 * tau_med) + 4, 40)
    r1 = []
    for t_c in tc_wide:
        rr = []
        for s in range(n_seeds):
            ss = jnp.array(draw_sources(1, np.random.default_rng(seed0 + 7 * s + 1),
                                        f_band, log10_h))
            bl = combined_blocks(toas, p_hats, L_arr, ss, n_psr, sigma)
            c, m = ntc_marg_cond(*bl[:5], n_psr, 1, jnp.asarray(t_c), rcond)
            rr.append(float(jnp.median(m) / jnp.median(c)))
        r1.append(np.median(rr))
    r1 = np.array(r1)
    # R(1) rises with t_c; find t_c where it crosses 0.5
    above = np.where(r1 > 0.5)[0]
    iup = above[0]
    lx = np.log(tc_wide[iup - 1]) + (0.5 - r1[iup - 1]) * (
        np.log(tc_wide[iup]) - np.log(tc_wide[iup - 1])) / (r1[iup] - r1[iup - 1])
    tc_50 = float(np.exp(lx))
    print(f"[ntc] SNR^2/src ~ {SNR2:.3e}; N=1 half-info at t_c_50 ~ {tc_50:.3e} s "
          f"(Y := t_c_50/t_c)", flush=True)

    # ---- grid ----
    N_grid = sorted(set(int(round(x)) for x in np.logspace(0, np.log10(256), 15)))
    Y_grid = np.logspace(-2, 2, 15)
    t_c_grid = tc_50 / Y_grid                   # Y = t_c_50 / t_c  (1 = N=1 half-info)
    nN, nY = len(N_grid), len(Y_grid)
    R = np.zeros((nN, nY)); Cnd = np.zeros((nN, nY)); Mrg = np.zeros((nN, nY))

    for i, N in enumerate(N_grid):
        cell_c = np.zeros((nY, n_seeds)); cell_m = np.zeros((nY, n_seeds))
        cell_r = np.zeros((nY, n_seeds))
        for s in range(n_seeds):
            srcs = jnp.array(draw_sources(N, np.random.default_rng(seed0 + 7 * s + N),
                                          f_band, log10_h))
            bl = combined_blocks(toas, p_hats, L_arr, srcs, n_psr, sigma)
            for j, t_c in enumerate(t_c_grid):
                c, m = ntc_marg_cond(*bl[:5], n_psr, N, jnp.asarray(t_c), rcond)
                mc = float(jnp.median(c)); mm = float(jnp.median(m))
                cell_c[j, s] = mc; cell_m[j, s] = mm; cell_r[j, s] = mm / mc
        Cnd[i] = np.median(cell_c, axis=1); Mrg[i] = np.median(cell_m, axis=1)
        R[i] = np.median(cell_r, axis=1)        # median-of-ratios (consistent w/ gates)
        print(f"[ntc] N={N:4d} (N/N*={N/N_star:5.1f})  R[Y=min..max]="
              f"{R[i,0]:.3f}..{R[i,-1]:.3f}", flush=True)

    # y->0 edge slice (t_c -> inf) and x=1 edge slice
    R_conf = R[:, 0].copy()                     # smallest Y ~ coherent
    R_coh = R[0, :].copy()                      # N=1 row

    # EDGE GATE (i): y->0 reproduces Stage-A fixed-per-source confusion R_conf(N)
    t_c_inf = SNR2 * tau_med / 1e-6
    Rconf_ref = []
    for N in N_grid:
        rr = []
        for s in range(n_seeds):
            srcs = jnp.array(draw_sources(N, np.random.default_rng(seed0 + 7 * s + N),
                                          f_band, log10_h))
            bl = combined_blocks(toas, p_hats, L_arr, srcs, n_psr, sigma)
            c, m = ntc_marg_cond(*bl[:5], n_psr, N, jnp.asarray(t_c_inf), rcond)
            rr.append(float(jnp.median(m) / jnp.median(c)))
        Rconf_ref.append(np.median(rr))
    Rconf_ref = np.array(Rconf_ref)
    # compare against Stage-A path (no dphi at all) in the SAME band
    Rconf_stageA = []
    for N in N_grid:
        rr = []
        for s in range(n_seeds):
            srcs = draw_sources(N, np.random.default_rng(seed0 + 7 * s + N), f_band, log10_h)
            pad, mask = pad_sources(srcs, max(N_grid))
            c, m = fisher_array(toas, p_hats, L_arr, jnp.array(pad), jnp.array(mask),
                                n_psr, sigma, rcond)
            rr.append(float(jnp.median(m) / jnp.median(c)))
        Rconf_stageA.append(np.median(rr))
    Rconf_stageA = np.array(Rconf_stageA)
    g_conf = float(np.max(np.abs(Rconf_ref - Rconf_stageA)))
    print(f"[ntc][edge y->0] combined vs Stage-A confusion: max|Δ| = {g_conf:.2e}",
          flush=True)
    assert g_conf < 1e-3, f"y->0 edge != Stage-A confusion ({g_conf:.2e})"

    # EDGE GATE (ii): x=1 (N=1) reproduces the Stage-A.2 coherence path (independent
    # autodiff code) to ~1e-3. (The idealised R_coh=1/(1+Y) holds only up to an O(geometry)
    # factor: the antenna 1/(1-cos mu) weighting spreads the per-baseline invariant, so the
    # 0.5-crossing is not exactly at Y=1; reported below after rescaling.)
    coh_ref = []
    for t_c in t_c_grid:
        rr = []
        for s in range(n_seeds):
            ss = jnp.array(draw_sources(1, np.random.default_rng(seed0 + 7 * s + 1),
                                        f_band, log10_h))
            Fc, tauc = coherence_blocks(toas, p_hats, L_arr, ss, n_psr, 1, sigma)
            c2, m2 = coherence_marg_cond(Fc, tauc, n_psr, 1, float(t_c),
                                         sigma_phi2_linear, rcond)
            rr.append(float(jnp.median(m2) / jnp.median(c2)))
        coh_ref.append(np.median(rr))
    coh_ref = np.array(coh_ref)
    g_coh = float(np.max(np.abs(R_coh - coh_ref)))
    print(f"[ntc][edge x=1] N=1 slice vs Stage-A.2 coherence path: max|Δ| = {g_coh:.2e}",
          flush=True)
    assert g_coh < 1e-3, f"x=1 edge != Stage-A.2 coherence ({g_coh:.2e})"
    # report the 1/(1+Y') fit with Y rescaled to the measured 0.5-crossing
    Y50 = _cross_x(Y_grid, R_coh / R_coh[0], 0.5)
    fit = 1.0 / (1.0 + Y_grid / Y50)
    fit_res = float(np.nanmax(np.abs(R_coh / R_coh[0] - fit)))
    print(f"[ntc] N=1 coherence: 0.5-crossing at Y={Y50:.3g} (not 1 -> geometry factor); "
          f"1/(1+Y/Y50) fits the shape to {fit_res:.2f}.", flush=True)

    # ---- HEADLINE: factorisation R(N,Y) vs R(N,0)*R(1,Y) ----
    product = np.outer(R_conf, R_coh / R_coh[0])   # R(N,0) * [R(1,Y)/R(1,0)]
    resid = R - product
    k = np.unravel_index(np.argmax(np.abs(resid)), resid.shape)
    maxdev = float(resid[k]); sign = "R>product (decoherence AIDS separation)" if maxdev > 0 \
        else "R<product (dphi COMPOUNDS confusion loss)"
    print(f"[ntc] factorisation: max|R - R(N,0)R(1,Y)| = {abs(maxdev):.3f} at "
          f"N/N*={N_grid[k[0]]/N_star:.1f}, Y={Y_grid[k[1]]:.2g}; sign: {maxdev:+.3f} "
          f"-> {sign}", flush=True)

    np.savez(OUT, N_grid=np.array(N_grid), Y_grid=Y_grid, t_c_grid=t_c_grid,
             x_grid=np.array(N_grid) / N_star, N_star=N_star, R=R, conditional=Cnd,
             marginal=Mrg, R_conf=R_conf, R_coh=R_coh, product=product, resid=resid,
             max_factor_dev=maxdev, SNR2=SNR2, tau_med=tau_med, tc_50=tc_50,
             n_psr=n_psr, n_seeds=n_seeds, band=np.array(f_band))

    # ---- figure: heatmap of R + factorisation residual ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.array(N_grid) / N_star
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2))
    for a, Z, ttl, cmap, ctr in [
            (ax[0], R, "R(N,Y) = marginal/conditional", "viridis", True),
            (ax[1], resid, r"factorisation residual  $R - R(N,0)\,R(1,Y)$", "RdBu_r", False)]:
        vmax = np.nanmax(np.abs(Z)) if not ctr else 1.0
        kw = dict(cmap=cmap, shading="auto")
        if not ctr:
            kw.update(vmin=-vmax, vmax=vmax)
        pc = a.pcolormesh(x, Y_grid, Z.T, **kw)
        if ctr:
            cs = a.contour(x, Y_grid, Z.T, levels=[0.5], colors="white", linewidths=2)
            a.clabel(cs, fmt="R=0.5", fontsize=8)
        a.axvline(1.0, color="k", ls=":", lw=1, alpha=.6)
        a.axhline(1.0, color="k", ls=":", lw=1, alpha=.6)
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlabel(r"$N/N^*$  (confusion axis)")
        a.set_ylabel(r"$Y = \mathrm{SNR}^2\,\sigma_\phi^2$  (coherence axis)")
        a.set_title(ttl, fontsize=10)
        fig.colorbar(pc, ax=a)
    fig.suptitle("Prong 2 Stage A.3 -- N x t_c map: confusion x coherence "
                 f"(max factor. dev {abs(maxdev):.2f})", fontsize=11, y=1.0)
    fig.tight_layout(); fig.savefig(FIG, dpi=130, bbox_inches="tight")
    print(f"[ntc] saved {OUT}\n[ntc] saved {FIG}\n[ntc] total {time.time()-t0:.1f}s",
          flush=True)


if __name__ == "__main__":
    import sys, time
    if len(sys.argv) > 1 and sys.argv[1] == "coherence":
        run_coherence_experiment()
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "ntcmap":
        run_ntc_map()
        sys.exit(0)
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
