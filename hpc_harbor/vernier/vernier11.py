#!/usr/bin/env python3
"""
VERNIER-11 -- population priors. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier7/vernier8 for the measured census ladder.

Model. For one source with everything known except amplitude A, the matched-filter
output is x ~ N(rho_true, 1) and the likelihood ratio at trial amplitude is
    L(rho) = exp(rho x - rho^2/2)
so the evidence for "a source is present" is
    Z = integral p(rho) exp(rho x - rho^2/2) d rho ,   lnB = ln Z.
Everything about the prior enters through p(rho). That is the whole experiment.

    FLAT prior      : uniform in log10 A over W dex (the toys' default)
    POPULATION prior: the census luminosity function, dN/dh ~ h^-gamma
                      => p(ln h) ~ h^-(gamma-1)

The safety question is whether the resulting change in lnB is a DIFFERENTIAL
effect (signal gains more than noise) or a COMMON-MODE Occam shift (both gain
equally, in which case recalibrating on the null erases it).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import vernier7 as S
import vernier8 as E

REPORTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "reports")
REPORTS = os.path.abspath(REPORTS)

W_DEX = 5.0                 # flat prior width in log10 amplitude
GAMMA = 4.0                 # dN/dh ~ h^-gamma (Euclidean)
RHO_REF = 20.0              # rank-1 array SNR, as everywhere in the series
NG = 900


def rho_grid(rho_ref=RHO_REF, w=W_DEX):
    lo = np.log10(rho_ref) - w
    return 10.0 ** np.linspace(lo, np.log10(rho_ref) + 0.3, NG)


def prior_flat(g):
    p = np.ones_like(g) / g                      # uniform in ln rho
    return p / np.trapezoid(p, g)


def prior_pop(g, gamma=GAMMA, shift_dex=0.0):
    """p(ln h) ~ h^-(gamma-1); shift_dex slides the LF's amplitude scale."""
    gs = g / 10.0 ** shift_dex
    p = gs ** -(gamma - 1.0) / g                 # convert to density in rho
    p = np.where(np.isfinite(p), p, 0.0)
    return p / np.trapezoid(p, g)


def prior_informed(g, rho_exp, width_dex=0.30, shift_dex=0.0, slope=None):
    """
    The prior a population analysis ACTUALLY supplies: where the loudest source
    should be. Log-normal in amplitude, centred on the order-statistic
    expectation with the census's measured spread (VERNIER-3: top-3 spread
    median 0.258 dex). shift_dex misspecifies the centre; slope adds a tilt.
    """
    lg = np.log10(g)
    c = np.log10(rho_exp) + shift_dex
    p = np.exp(-0.5 * ((lg - c) / width_dex) ** 2) / g
    if slope is not None:
        p = p * g ** slope
    return p / np.trapezoid(p, g)


def lnB(x, g, p):
    """ln evidence for a source, marginalised over amplitude."""
    lw = g * x - 0.5 * g ** 2
    m = lw.max()
    return m + np.log(np.trapezoid(p * np.exp(lw - m), g))


def run(n_trial=20000, seed=0, gamma=GAMMA, shift=0.0, ranks=25):
    rng = np.random.default_rng(seed)
    g = rho_grid()
    pf, pp = prior_flat(g), prior_pop(g, gamma, shift)
    rho_true = RHO_REF * 10.0 ** S.LAD[:ranks]
    # NULL: no source
    xn = rng.normal(0.0, 1.0, n_trial)
    bn_f = np.array([lnB(x, g, pf) for x in xn])
    bn_p = np.array([lnB(x, g, pp) for x in xn])
    # SIGNAL: one draw per rank, many realisations
    bs_f, bs_p, rk = [], [], []
    per = max(n_trial // ranks, 200)
    for i, r in enumerate(rho_true):
        xs = rng.normal(r, 1.0, per)
        bs_f += [lnB(x, g, pf) for x in xs]
        bs_p += [lnB(x, g, pp) for x in xs]
        rk += [i + 1] * per
    return (np.array(bn_f), np.array(bn_p), np.array(bs_f), np.array(bs_p),
            np.array(rk), rho_true)
