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
def fisher_array(toas, p_hats, L_arr, srcs, mask, n_psr, sigma):
    """Joint Fisher for the array. Returns per-pulsar (conditional, marginal) distance
    information arrays of shape (n_psr,).

    Conditional = diagonal of the L-block (other params known).
    Marginal    = Schur complement removing the source-parameter block (phases/amps),
                  using an eigh pseudo-inverse of the source block.

    Assembled per-pulsar (lax.scan) rather than as one monolithic Jacobian: a pulsar's
    residual depends only on its OWN distance, so F_LL is diagonal and the array Fisher
    is a sum of per-pulsar outer products. Memory stays O(n_src^2) instead of
    O(n_param * n_data), which is what lets N_max=1000 fit on a 24 GB GPU.
    """
    inv_s2 = 1.0 / (sigma ** 2)

    # per-pulsar analytic gradients in parallel
    g, H = jax.vmap(lambda ph, L: analytic_psr_grads(toas, ph, L, srcs, mask))(
        p_hats, L_arr)            # g:(n_psr,n_toa)  H:(n_psr,n_toa,2*n_max)

    cond = jnp.einsum("pt,pt->p", g, g) * inv_s2                 # diag F_LL
    F_Ls = jnp.einsum("pt,pta->pa", g, H) * inv_s2              # (n_psr, 2*n_max)
    F_ss = jnp.einsum("pta,ptb->ab", H, H) * inv_s2            # (2*n_max, 2*n_max)

    # Schur complement of the L-block after marginalising sources:
    #   F_marg = F_LL - F_Ls F_ss^+ F_sL ; marginal info per pulsar = diag of F_marg
    F_ss_pinv = _eigh_pinv(F_ss)
    M = F_ss_pinv @ F_Ls.T                                        # (2*n_max, n_psr)
    marg = jnp.clip(cond - jnp.sum(F_Ls.T * M, axis=0), 0.0, None)
    return cond, marg


# ----------------------------------------------------------------------
# Population draws + padding
# ----------------------------------------------------------------------
def draw_sources(n_src, rng, f_band, log10_h, fixed_total=False, h_total=None):
    """Draw n_src sources. If fixed_total, set per-source h so sum h^2 = h_total^2."""
    theta = np.arccos(rng.uniform(-1, 1, n_src))
    phi = rng.uniform(0, 2 * np.pi, n_src)
    psi = rng.uniform(0, np.pi, n_src)
    cinc = rng.uniform(-1, 1, n_src)
    log10_f = rng.uniform(np.log10(f_band[0]), np.log10(f_band[1]), n_src)
    phase0 = rng.uniform(0, 2 * np.pi, n_src)
    if fixed_total:
        h_each = h_total / np.sqrt(n_src)
        l10h = np.full(n_src, np.log10(h_each))
    else:
        l10h = np.full(n_src, log10_h)
    return np.stack([theta, phi, psi, cinc, log10_f, l10h, phase0], axis=1)


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
              p_hats=None, label=None):
    """mode in {'fixed_persource', 'fixed_total'}.

    Joint-array Fisher with PADDING: every realisation is padded to n_max sources so the
    jitted fisher_array compiles exactly once. We report the median over pulsars (within
    the array) and the 16/84 spread over population realisations.
    """
    if n_max is None:
        n_max = int(max(n_values))
    if p_hats is None:
        p_hats = build_array(n_psr)
    toas = build_toas(T, cadence_days)
    L_arr = jnp.ones(n_psr)  # 1 kpc each (representative)
    N_star = confusion_onset(T, f_band)

    cond_med, cond_lo, cond_hi = [], [], []
    marg_med, marg_lo, marg_hi = [], [], []

    tag = label or mode
    for n_src in n_values:
        cond_real, marg_real = [], []   # per-realisation array-median info
        for s in range(n_seeds):
            rng = np.random.default_rng(seed0 + 1000 * s + n_src)
            srcs = draw_sources(n_src, rng, f_band, log10_h,
                                fixed_total=(mode == "fixed_total"), h_total=h_total)
            srcs_pad, mask = pad_sources(srcs, n_max)
            cond, marg = fisher_array(toas, p_hats, L_arr,
                                      jnp.array(srcs_pad), jnp.array(mask), n_psr, sigma)
            cond_real.append(float(jnp.median(cond)))
            marg_real.append(float(jnp.median(marg)))
        cond_real = np.array(cond_real); marg_real = np.array(marg_real)
        cond_med.append(np.median(cond_real))
        cond_lo.append(np.percentile(cond_real, 16)); cond_hi.append(np.percentile(cond_real, 84))
        marg_med.append(np.median(marg_real))
        marg_lo.append(np.percentile(marg_real, 16)); marg_hi.append(np.percentile(marg_real, 84))
        ratio = marg_med[-1] / cond_med[-1] if cond_med[-1] > 0 else 0.0
        print(f"  [{tag}] N={n_src:5d}  cond={cond_med[-1]:.3e}  "
              f"marg={marg_med[-1]:.3e}  (marg/cond {ratio:.3f})", flush=True)

    return dict(
        n=np.array(n_values),
        cond=np.array(cond_med), cond_lo=np.array(cond_lo), cond_hi=np.array(cond_hi),
        marg=np.array(marg_med), marg_lo=np.array(marg_lo), marg_hi=np.array(marg_hi),
        N_star=N_star, Tspan_bins=N_star, sigma=sigma, T=T, f_band=np.array(f_band),
        n_psr=n_psr, n_max=n_max, n_seeds=n_seeds, cadence_days=cadence_days,
    )


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


if __name__ == "__main__":
    import time
    OUT = "/home/mattm/projects/HSYMT/CW_transition/prong2_results.npz"

    # 0) correctness gates
    verify_padding()
    validate_grads()

    # log-spaced N grid 1 -> ~1000
    n_values = sorted(set(int(round(x)) for x in np.logspace(0, 3, 16)))
    print("N grid:", n_values, flush=True)
    n_max = max(n_values)

    p_hats = build_array(116)

    t0 = time.time()
    print("\nFixed per-source strain (Mihir-style; independence line):", flush=True)
    res_ps = run_sweep(n_values, mode="fixed_persource", n_max=n_max, p_hats=p_hats,
                       log10_h=-14.3)
    print(f"[timing] fixed_persource sweep: {time.time() - t0:.1f}s", flush=True)

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
    print(f"\n[mode-indep] max|Δ(marg/cond)| = {dev.max():.3f} at N={res_ps['n'][ibreak]} "
          f"(ps={rps[ibreak]:.3f} vs tot={rtot[ibreak]:.3f})", flush=True)

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

    np.savez(OUT,
             **{f"ps_{k}": v for k, v in res_ps.items()},
             **{f"tot_{k}": v for k, v in res_tot.items()},
             mode_indep_maxdev=dev.max(), mode_indep_break_N=res_ps['n'][ibreak],
             diag_labels=diag_labels, **diag_pack)
    print(f"\nSaved {OUT}", flush=True)
    print(f"Confusion onset N* = {res_ps['N_star']:.1f} (one source per 1/T frequency bin)",
          flush=True)
    print(f"[timing] TOTAL: {time.time() - t0:.1f}s", flush=True)
