#!/usr/bin/env python3
"""
VERNIER-15 -- deconfusion. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier8.

This is an IDENTIFIABILITY question, not a sensitivity one. Put M sources in the
SAME frequency bin at comparable brightness and ask how many can be separated.

COUNTING PREDICTION
  Per bin each pulsar supplies one complex Fourier coefficient = 2 real numbers.
  Earth terms alone -> 2 N_psr measurements.
  Because f_P != f_E, each pulsar's PULSAR term lands in its OWN separate bin, so
  with known distances it is an extra independent complex datum -> 4 N_psr.
  Each source costs ~7 parameters (sky x2, A, phase, psi, cos iota, f).
      M_max(Earth only) ~ 2 N / 7        M_max(both) ~ 4 N / 7      ratio ~ 2
Separability is measured from the FISHER matrix of the full M-source model: a
source is resolved when its Fisher-predicted sky and amplitude errors stay inside
tolerance and the matrix stays invertible.
"""
import os
import numpy as np
import vernier_toy as V
import vernier8 as E

REPORTS = V.REPORTS
# Within ONE bin the Earth term has NO frequency dependence -- that is what "same
# bin" means -- so f is not a free parameter for the Earth-only model and including
# it makes the Fisher singular by construction. With pulsar terms the phase carries
# 2 pi f tau with a millennia-long lever arm, so f BECOMES measurable: that is
# itself one of the extra degrees of freedom the pulsar term supplies.
NPAR_E = 6            # th, ph, A, psi, cos iota, phase
NPAR_P = 7            # ... + f


def basis(omega):
    """Polarisation basis orthogonal to the propagation direction."""
    z = np.array([0.0, 0.0, 1.0])
    u = np.cross(z, omega)
    n = np.linalg.norm(u)
    u = u / n if n > 1e-12 else np.array([1.0, 0.0, 0.0])
    v = np.cross(omega, u)
    return u, v


def antenna(p, omega):
    """F+, Fx for pulsar direction p and source direction omega."""
    u, v = basis(omega)
    d = 1.0 + float(p @ omega)
    if abs(d) < 1e-6:
        d = 1e-6
    pu, pv = float(p @ u), float(p @ v)
    return (pu ** 2 - pv ** 2) / (2 * d), (pu * pv) / d


def response(p, th, ph, A, psi, ci, phase, f, tau=None):
    """Complex response of one pulsar to one source (Earth term, or pulsar term)."""
    om = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    Fp, Fc = antenna(p, om)
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    Fp, Fc = Fp * c2 + Fc * s2, -Fp * s2 + Fc * c2
    ap, ac = (1 + ci ** 2) / 2.0, ci
    ph_tot = phase if tau is None else phase - 2 * np.pi * f * tau
    return A * (Fp * ap + 1j * Fc * ac) * np.exp(1j * ph_tot)


def data_vector(P, srcs, taus=None, f=1e-8):
    """
    Stacked complex data: Earth bin for every pulsar, plus pulsar-term bins.

    taus must be shaped (n_psr, n_src): the light-travel delay depends on the
    angle between THAT pulsar and THAT source. An earlier version averaged tau
    over sources, which made every source share one delay per pulsar -- the
    pulsar-term block then reduced to the Earth block times a common phase and
    carried no independent information at all (rank 5 where 14 was needed).
    """
    out = []
    for a, p in enumerate(P):
        out.append(sum(response(p, *s, tau=None) for s in srcs))
    if taus is not None:
        for a, p in enumerate(P):
            out.append(sum(response(p, *s, tau=taus[a][i])
                           for i, s in enumerate(srcs)))
    return np.array(out)


def fisher(P, srcs, taus=None, f=1e-8, eps=1e-5):
    """Fisher matrix of the full M-source model (real+imag treated as data)."""
    M = len(srcs)
    npar = NPAR_P if taus is not None else NPAR_E
    base = data_vector(P, srcs, taus, f)
    J = np.zeros((2 * base.size, M * npar))
    for i, s in enumerate(srcs):
        for k in range(npar):
            sp = [list(x) for x in srcs]
            step = eps * (abs(sp[i][k]) + 1.0)
            sp[i][k] += step
            d = (data_vector(P, sp, taus, f) - base) / step
            J[:, i * npar + k] = np.concatenate([d.real, d.imag])
    return J.T @ J


def draw(rng, n_psr, M, L_lo=300.0, L_hi=5000.0, f=1e-8):
    P = V.iso_dirs(n_psr, rng)
    L = 10 ** rng.uniform(np.log10(L_lo), np.log10(L_hi), n_psr)
    srcs = []
    for _ in range(M):
        srcs.append([np.arccos(rng.uniform(-1, 1)), rng.uniform(0, 2 * np.pi),
                     1.0, rng.uniform(0, np.pi), rng.uniform(-1, 1),
                     rng.uniform(0, 2 * np.pi), f])
    om = np.array([[np.sin(s[0]) * np.cos(s[1]), np.sin(s[0]) * np.sin(s[1]),
                    np.cos(s[0])] for s in srcs])
    taus = np.array([[V.tau_of(L[a], max(1 - float(P[a] @ o), 5e-2)) for o in om]
                     for a in range(n_psr)])          # (n_psr, n_src) -- per PAIR
    return P, srcs, taus, f


def rank_test(rng, n_psr, M, use_pterm, tol=1e-10):
    """
    Identifiability by the EFFECTIVE RANK of the Fisher matrix.

    pinv is the wrong tool here: it truncates degenerate directions and returns a
    SMALL pseudo-error for them, so a singular Fisher reads as a precise
    measurement. An earlier version of this file did that and produced a
    non-monotonic Earth-only curve (0 resolved at M=10, 30 at M=30). The honest
    test is whether the Fisher is full rank: M sources are separable only if
    rank(F) = M * npar, i.e. no exactly degenerate parameter combination.
    """
    P, srcs, taus, f = draw(rng, n_psr, M)
    tt = taus if use_pterm else None
    npar = NPAR_P if use_pterm else NPAR_E
    F = fisher(P, srcs, tt, f)
    # Parameters have wildly different units: d/df in the pulsar term carries a
    # factor tau ~ 1e11 s, so one eigenvalue swamps the spectrum and a relative
    # rank tolerance discards everything else. Normalise to correlation form
    # first -- the rank is then scale-invariant, which is what identifiability
    # actually means.
    d = np.sqrt(np.clip(np.diag(F), 1e-300, None))
    F = F / np.outer(d, d)
    ev = np.linalg.eigvalsh(F)
    ev = ev[ev > 0]
    rank = int(np.sum(ev > tol * ev.max())) if ev.size else 0
    need = M * npar
    ndata = 2 * n_psr * (2 if use_pterm else 1)
    return rank, need, ndata
