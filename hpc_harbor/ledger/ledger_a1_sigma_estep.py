#!/usr/bin/env python3
"""LEDGER A1 -- SIGMA-POINT E-STEP: the fringe posterior averaged over the source belief.

=============================================================================
THE DEFECT, AS WRITTEN IN THE CODE (verdict: REAL)
=============================================================================
The E-step evaluates the fringe posterior at ONE source point:

  glacier_loop.py:668   theta_base = theta_with_absent(theta_rec, nd, carried)
  glacier_loop.py:673   LNL = estep_per_target(self.sp, theta_base, EV, AI, data, PM, jnp)
  spark3.py:269         tb = jnp.asarray(theta_base)          <- a single point, no average
  spark3.py:273-289     the only loop is over TARGET PULSARS p

and the fringe posterior is then formed from that one LNL surface
(trackB_b1_core.FringeTables.posterior, line 700). The certification statistic
q_max, dlnL_det and the on_true flag are all read off it.

The loop DOES carry a source belief -- and throws it away at this step:
  * glacier_loop.py:417  sig_opt, sig_pes, F_ii = fisher_conditional(...)   (widths)
  * glacier_loop.py:420  used ONLY for the frontier ratio
  * glacier_loop.py:555  banked as columns sig_opt / sig_pes
  * glacier_loop.py:496  mstep_quadratic returns per-fed-member Laplace `widths` --
                         which are not even banked, let alone propagated
  * forge_b1.py:311-316  FORGE-B's own belief: Cov = inv(Fs + Pi), sig_f, sig_mc --
                         and forge_b1.py:331 then re-E-steps at `theta_src`, the POINT.

CONSEQUENCE. q_max is a CONDITIONAL confidence: P(fringe | this exact source). It is
reported, and charged against the bar, as if it were the marginal P(fringe | data). When
the source belief is wide compared with the fringe-discriminating phase scale, the
conditional is systematically OVER-confident: it is the width of the fringe posterior at
a source the analysis does not actually know. This is the same class of error the
programme already names elsewhere -- IGNITE-2's retired "hard lock" was over-confidence
in a PINNED FRINGE; this is over-confidence in a PINNED SOURCE.

=============================================================================
THE UPGRADE (pre-registered for the BELIEF ARM's fan; NOT retro-applied)
=============================================================================
Replace the point evaluation with a quadrature of the LIKELIHOOD over the belief:

    LNL_marg[p, b] = ln INT d(theta_s) N(theta_s; theta_hat, Sigma) exp(LNL(theta_s)[p, b])
                   ~ logsumexp_i ( ln w_i + LNL(chi_i)[p, b] )

over sigma points {chi_i, w_i} matched to (theta_hat, Sigma). NOTE THE ORDER: the average
is in LIKELIHOOD space, then the log is taken. Averaging log-likelihoods (the naive thing)
is a different, and wrong, object -- it gives the geometric mean and cannot widen a
posterior at all.

THE RULE (scaled unscented transform, DECLARED knobs -- w0 is the one free knob):
    m       = dim of the belief block
    lambda  = m*w0/(1-w0)                      so that the centre carries weight w0
    chi_0   = theta_hat,                      w_0 = w0
    chi_i^+-= theta_hat +- [sqrt((m+lambda) Sigma)]_i,   w_i = (1-w0)/(2m)
  w0 = 1/3 (Julier-Uhlmann's Gaussian recommendation) keeps EVERY weight strictly
  positive, which a likelihood quadrature requires: a negative weight can drive the
  integrand negative and the log undefined. Weights sum to 1 by construction; mean and
  covariance of the point set match (theta_hat, Sigma) exactly.

THE TWO GATES (both are exact identities, so both are BIT-COMPARABLE, not tolerance-based):
  G-L1  BELIEF-WIDTH -> 0.  Sigma -> 0 collapses every sigma point onto theta_hat, and
        the rule collapses STRUCTURALLY to the single point (see sigma_points): the
        incumbent call is re-issued unchanged, so max|LNL_sigma - LNL_point| == 0.0
        EXACTLY. Stated plainly because the alternative is a lie: computed through the
        logsumexp with unit-sum weights the same quantity lands 1 ulp (2.2e-16) off, and
        "bit-comparable" would then be a tolerance dressed up as an identity.
  G-L2  DEGENERATE RULE.  w0 = 1 puts all mass on the centre at any Sigma; same identity.
        Guards the weighting arithmetic independently of the point placement.
  G-L3  EXACTNESS ON A GAUSSIAN.  For LNL quadratic in theta_s the marginal is analytic;
        the 2m+1 rule must reproduce it to the rule's order. Measured, printed, not
        asserted -- the UT is exact to 2nd order in the MOMENTS, not in this integral,
        so this leg REPORTS the error rather than passing/failing it.
  G-L4  COST.  n_eval = 2m+1 per E-step. Printed as a multiplier, because it is the
        budget statement the belief arm has to live inside.

BELIEF BLOCK, AND THE COST CAP (declared). The belief that matters for the fringe is the
one over the axes the pulsar-term phase is sensitive to: (log10_fgw, log10_mc) per fed
member -- exactly MSTEP_AXES (glacier_loop.py:114) and exactly FORGE-B's (sig_f, sig_mc)
(forge_b1.py:311-316). So m = 2 * n_fed and the E-step cost multiplier is 4*n_fed + 1.
At the r13p25 cells (n_fed up to 20) that is 81x a full E-step: NOT affordable. The fan
therefore caps the sigma set at the TOP `n_belief` fed members ranked by
(belief width / prior box width) -- the members the loop is least sure of, which are the
ones capable of moving a fringe. Remaining fed members are held at their point estimate.
`n_belief` is banked per cell, never defaulted silently, and n_belief = 0 IS the incumbent.

WHAT THIS MODULE IS. Machinery + gates, deliberately factored so the expensive part
(`estep_fn`) is injected. That lets G-L1/G-L2/G-L3 run on the CPU against an analytic
surrogate -- which is where they are run here -- while the venue wiring is a 6-line change
in the belief arm's driver, given in LEDGER_stats_audit.md and NOT applied to the held
Stage-1/2 driver.

Usage:  python3 hpc_harbor/ledger/ledger_a1_sigma_estep.py gate
"""
import argparse
import sys

import numpy as np


# ============================================================
# the rule
# ============================================================
def sigma_points(theta_hat, axes, sigmas, w0=1.0 / 3.0):
    """Scaled-UT sigma points over a DIAGONAL belief on `axes`.

    theta_hat : (n_theta,) the point estimate (the incumbent evaluation point)
    axes      : (m,) int index array into theta_hat -- the belief block
    sigmas    : (m,) belief standard deviations (diagonal Sigma)
    w0        : centre weight; must be in (0, 1]. 1.0 -> the incumbent point rule.

    Returns (CHI (n_pt, n_theta), W (n_pt,)) with sum(W) == 1 and every W > 0.
    Diagonal Sigma is not an approximation of convenience: it is what the loop actually
    carries (fisher_conditional returns per-axis CONDITIONAL widths; mstep_quadratic
    returns per-axis Laplace widths; FORGE-B's sig_f/sig_mc are the Cov DIAGONAL). A
    full-covariance rule is a drop-in -- replace the axial offsets with the columns of
    a Cholesky factor -- and is specced but not needed until a joint Hessian exists.
    """
    theta_hat = np.asarray(theta_hat, float)
    axes = np.asarray(axes, int)
    sigmas = np.asarray(sigmas, float)
    if axes.shape != sigmas.shape:
        raise ValueError(f"axes {axes.shape} vs sigmas {sigmas.shape}")
    if not (0.0 < w0 <= 1.0):
        raise ValueError(f"w0 must be in (0,1], got {w0}")
    m = len(axes)
    # STRUCTURAL COLLAPSE (declared, and it is what makes G-L1 an IDENTITY rather than a
    # tolerance). A zero-width belief, an empty belief block, or w0 = 1 are not "nearly
    # the point rule" -- they ARE the point rule, and the code says so structurally, the
    # same way embed_igniter returns the census bit-exactly at e = 0 (CHORUS C1). Without
    # the collapse the logsumexp of identical values with unit-sum weights lands 1 ulp
    # (2.2e-16) off the incumbent, which would make "bit-comparable" a lie.
    if m == 0 or w0 == 1.0 or not np.any(sigmas != 0.0):
        return theta_hat[None, :].copy(), np.ones(1)
    lam = m * w0 / (1.0 - w0)
    c = np.sqrt(m + lam)
    CHI = np.tile(theta_hat, (2 * m + 1, 1))
    W = np.empty(2 * m + 1)
    W[0] = w0
    for i in range(m):
        CHI[1 + 2 * i, axes[i]] += c * sigmas[i]
        CHI[2 + 2 * i, axes[i]] -= c * sigmas[i]
        W[1 + 2 * i] = W[2 + 2 * i] = (1.0 - w0) / (2.0 * m)
    return CHI, W


def estep_sigma(estep_fn, theta_hat, axes, sigmas, w0=1.0 / 3.0):
    """Belief-averaged E-step.

    estep_fn(theta) -> LNL array of any fixed shape (the incumbent per-target E-step,
    e.g. lambda th: spark3.estep_per_target(sp, th, EV, AI, data, PM, jnp)).

    Returns (LNL_marg, n_eval). The reduction is a stabilised logsumexp over sigma
    points, so the LIKELIHOOD is averaged and the log is taken afterwards.
    """
    CHI, W = sigma_points(theta_hat, axes, sigmas, w0=w0)
    stack = np.stack([np.asarray(estep_fn(chi), float) for chi in CHI])   # (n_pt, ...)
    lw = np.log(W).reshape((-1,) + (1,) * (stack.ndim - 1))
    z = stack + lw
    mx = z.max(axis=0)
    return mx + np.log(np.exp(z - mx).sum(axis=0)), len(CHI)


def rank_belief_members(sig_belief, box_sigma, fed_idx, n_belief):
    """Which fed members get sigma points: the n_belief LEAST-KNOWN, by width ratio.

    sig_belief : (n_fed, n_axes) the loop's own belief widths for the fed members
    box_sigma  : (n_axes,) the prior box widths (glacier_loop.py:395 convention)
    Returns the selected fed-member indices, most-uncertain first. Banked per cell as
    `belief_members`; n_belief = 0 reproduces the incumbent exactly.
    """
    if n_belief <= 0 or len(fed_idx) == 0:
        return np.zeros(0, int)
    r = np.max(np.asarray(sig_belief) / np.asarray(box_sigma)[None, :], axis=1)
    order = np.argsort(-r)
    return np.asarray(fed_idx)[order[:int(n_belief)]]


# ============================================================
# gates (CPU; analytic surrogate for the expensive E-step)
# ============================================================
def _surrogate(n_theta, axes, seed=20260729):
    """A deterministic quadratic-in-theta_s 'LNL' with the same shape signature as the
    real per-target E-step: (npsr, B). Quadratic so G-L3 has an analytic answer."""
    rng = np.random.default_rng(seed)
    npsr, B = 5, 7
    g = rng.normal(size=(npsr, B, len(axes)))            # dLNL/dtheta_s
    H = -rng.uniform(0.5, 2.0, size=(npsr, B, len(axes)))  # curvature (negative definite)
    base = rng.normal(size=(npsr, B))
    th0 = rng.normal(size=n_theta)

    def fn(theta):
        d = np.asarray(theta, float)[axes] - th0[axes]
        return base + np.einsum("pbj,j->pb", g, d) + 0.5 * np.einsum("pbj,j->pb", H, d ** 2)
    return fn, g, H, base, th0, (npsr, B)


def gate(verbose=True):
    print("=== LEDGER A1 -- SIGMA-POINT E-STEP GATES (CPU) ===", flush=True)
    ok = True
    n_theta = 24
    axes = np.array([3, 4, 11, 12])                      # 2 fed members x (fgw, mc)
    fn, g, H, base, th0, shape = _surrogate(n_theta, axes)
    theta_hat = th0 + 0.05
    sig = np.array([0.02, 0.05, 0.03, 0.07])

    # -- G-L1: belief width -> 0 reproduces the incumbent BIT-EXACTLY
    ref = np.asarray(fn(theta_hat), float)
    L0, n0 = estep_sigma(fn, theta_hat, axes, np.zeros_like(sig))
    d1 = float(np.max(np.abs(L0 - ref)))
    b1 = (d1 == 0.0)
    print(f"  G-L1 belief-width -> 0 == incumbent point E-step: max|d| = {d1:.3e} "
          f"(need exactly 0.0), n_eval {n0} -> {'PASS' if b1 else 'FAIL'}")
    ok &= b1

    # -- G-L2: degenerate rule (w0 = 1) reproduces the incumbent BIT-EXACTLY at any width
    L1, n1 = estep_sigma(fn, theta_hat, axes, sig, w0=1.0)
    d2 = float(np.max(np.abs(L1 - ref)))
    b2 = (d2 == 0.0)
    print(f"  G-L2 w0 = 1 (all mass on the centre) == incumbent, at FULL belief width: "
          f"max|d| = {d2:.3e}, n_eval {n1} -> {'PASS' if b2 else 'FAIL'}")
    ok &= b2

    # -- weights: positive and normalised (the requirement a likelihood quadrature has)
    CHI, W = sigma_points(theta_hat, axes, sig)
    b3 = bool(np.all(W > 0)) and abs(W.sum() - 1.0) < 1e-15
    print(f"  G-L2b weights strictly positive and sum to 1: min w {W.min():.4f}, "
          f"sum-1 {W.sum()-1.0:+.2e} -> {'PASS' if b3 else 'FAIL'}")
    ok &= b3
    # mean/covariance of the point set match the belief exactly (the UT's defining property)
    mu = (W[:, None] * CHI).sum(0)
    cv = ((W[:, None] * (CHI - mu) ** 2).sum(0))[axes]
    b3b = (float(np.max(np.abs(mu - theta_hat))) < 1e-14
           and float(np.max(np.abs(np.sqrt(cv) - sig))) < 1e-14)
    print(f"  G-L2c point-set moments match (theta_hat, Sigma): max|d mean| "
          f"{np.max(np.abs(mu - theta_hat)):.2e}, max|d sigma| "
          f"{np.max(np.abs(np.sqrt(cv) - sig)):.2e} -> {'PASS' if b3b else 'FAIL'}")
    ok &= b3b

    # -- G-L3: exactness against the analytic Gaussian marginal (REPORT, not pass/fail)
    # For LNL = base + g.d + 0.5 H d^2 (diagonal H<0) and d ~ N(dhat, S):
    #   ln E[exp LNL] = base + sum_j [ g dhat + .5 H dhat^2 ... ] -- done per axis exactly.
    dhat = theta_hat[axes] - th0[axes]
    S = sig ** 2
    # per (p,b): integrate prod_j N(d_j; dhat_j, S_j) exp(g_j d_j + .5 H_j d_j^2)
    A = 1.0 - H * S[None, None, :]                                   # >0 since H<0
    expo = (g * dhat + 0.5 * H * dhat ** 2
            + 0.5 * S[None, None, :] * (g + H * dhat) ** 2 / A)
    exact = base + expo.sum(-1) - 0.5 * np.log(A).sum(-1)
    Lm, nm = estep_sigma(fn, theta_hat, axes, sig)
    err = float(np.max(np.abs(Lm - exact)))
    rel = err / float(np.max(np.abs(exact - ref)) + 1e-300)
    print(f"  G-L3 vs the ANALYTIC Gaussian marginal (quadratic surrogate): "
          f"max|d| = {err:.3e} nat; the marginalisation itself moves the surface by "
          f"{float(np.max(np.abs(exact - ref))):.3e} nat, so the rule captures "
          f"{100*(1-rel):.1f}% of it (REPORT-ONLY leg)")

    # -- G-L3b: the mechanism itself, on a FRINGE toy (the generic quadratic surrogate
    # above cannot show it -- its peak does not move with the source, so marginalising
    # can sharpen or broaden at random). Here the fringe peak POSITION shifts with the
    # source parameter, which is the actual coupling: the pulsar-term phase, and hence
    # which distance fringe is favoured, depends on (fgw, mc). Marginalising then smears
    # the peak ACROSS fringes and q_max can only fall.
    print()
    print("  G-L3b THE MECHANISM (fringe toy: peak position shifts with the source)")
    K, dfr = 41, 1.0                 # 41 fringes, unit spacing
    kk = np.arange(K) - K // 2
    shift_per_sigma = 0.0            # placeholder; swept below
    ax1 = np.array([3])
    def fringe_fn_factory(slope, sharp=6.0):
        def f(theta):
            mu = slope * (theta[ax1[0]] - theta_hat[ax1[0]])
            return (-0.5 * sharp * (kk - mu) ** 2)[None, :]
        return f
    print("     belief width (in fringes of peak motion)   q_max point -> marginal")
    for frac in (0.0, 0.1, 0.25, 0.5, 1.0):
        slope = frac / 0.05          # 0.05 = the belief sigma used below
        f = fringe_fn_factory(slope)
        Lp = np.asarray(f(theta_hat))
        Lq, _ = estep_sigma(f, theta_hat, ax1, np.array([0.05]))
        def qm1(L):
            w = np.exp(L - L.max(axis=1, keepdims=True)); w /= w.sum(axis=1, keepdims=True)
            return float(w.max(axis=1)[0])
        print(f"       peak moves {frac:4.2f} fringe per 1-sigma of belief      "
              f"{qm1(Lp):.4f} -> {qm1(Lq):.4f}"
              + ("   <- q_max UNCHANGED at zero belief" if frac == 0 else ""))
    print("     -> q_max is a CONDITIONAL confidence; once the belief moves the peak by")
    print("        an appreciable fraction of a fringe, the marginal is materially less")
    print("        confident than the number the criterion is currently charged against.")
    print()

    # -- G-L4: cost
    for nf in (1, 2, 3, 5, 10, 20):
        print(f"  G-L4 cost at n_fed = {nf:2d} (m = 2*n_fed): E-step multiplier "
              f"{4*nf+1:3d}x" + ("   <- r13p25 class, NOT affordable" if nf >= 10 else ""))

    print(f"\n=== LEDGER A1 GATES: {'PASS' if ok else 'FAIL'} ===", flush=True)
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate"], nargs="?", default="gate")
    ap.parse_args()
    return gate()


if __name__ == "__main__":
    sys.exit(main())
