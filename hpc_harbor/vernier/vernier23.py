#!/usr/bin/env python3
"""
VERNIER-23 -- drop the certification requirement. EXPLORATORY, NOT pre-registered.
Clean room: imports only vernier_toy / vernier15.

HYPOTHESIS. V15's separability gain (8 -> 20) used KNOWN distances. But the lever
arm is the light-travel TIME, not knowledge of it. Add the N_psr distances as FREE
parameters with their existing priors and MARGINALISE them (Schur complement of the
Fisher). If the source-parameter block stays full rank, separation is an
array-level measurement and the per-pulsar certification bar never enters.

COUNTING PREDICTION
  data = 4 N_psr (Earth block + pulsar-term block)
  (a) Earth only      : M*6 <= 2 N        -> M_max = 2N/6
  (b) distances known : M*7 <= 4 N        -> M_max = 4N/7
  (c) marginalised    : M*7 + N <= 4 N    -> M_max = 3N/7
so (c) should sit BETWEEN (a) and (b), not collapse to (a).
"""
import os
import numpy as np

import vernier_toy as V
import vernier15 as D

REPORTS = V.REPORTS


def fisher_with_distances(P, L, srcs, f, eps=1e-5):
    """
    Jacobian over [source params ..., pulsar distances ...].
    Distances enter only through the pulsar-term block.
    """
    n_psr = P.shape[0]
    M = len(srcs)
    npar = D.NPAR_P
    om = np.array([[np.sin(s[0]) * np.cos(s[1]), np.sin(s[0]) * np.sin(s[1]),
                    np.cos(s[0])] for s in srcs])
    def taus_of(Lv):
        return np.array([[V.tau_of(Lv[a], max(1 - float(P[a] @ o), 5e-2))
                          for o in om] for a in range(n_psr)])
    base = D.data_vector(P, srcs, taus_of(L), f)
    ncol = M * npar + n_psr
    J = np.zeros((2 * base.size, ncol))
    for i in range(M):
        for k in range(npar):
            sp = [list(x) for x in srcs]
            st = eps * (abs(sp[i][k]) + 1.0)
            sp[i][k] += st
            d = (D.data_vector(P, sp, taus_of(L), f) - base) / st
            J[:, i * npar + k] = np.concatenate([d.real, d.imag])
    for a in range(n_psr):
        Lv = L.copy(); st = eps * (abs(Lv[a]) + 1.0); Lv[a] += st
        d = (D.data_vector(P, srcs, taus_of(Lv), f) - base) / st
        J[:, M * npar + a] = np.concatenate([d.real, d.imag])
    return J.T @ J, M * npar, n_psr


def rank_marginalised(rng, n_psr, M, sigma_L_pc, tol=1e-10, f=1e-8):
    """Rank of the SOURCE block after marginalising the pulsar distances."""
    P, srcs, taus, _ = D.draw(rng, n_psr, M, f=f)
    L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), n_psr)
    F, ns, nd = fisher_with_distances(P, L, srcs, f)
    # distance priors enter the distance block
    Fdd = F[ns:, ns:] + np.eye(nd) / sigma_L_pc ** 2
    Fss = F[:ns, :ns]; Fsd = F[:ns, ns:]
    try:
        Fm = Fss - Fsd @ np.linalg.solve(Fdd, Fsd.T)
    except np.linalg.LinAlgError:
        return 0, ns
    d = np.sqrt(np.clip(np.diag(Fm), 1e-300, None))
    Fm = Fm / np.outer(d, d)
    ev = np.linalg.eigvalsh(Fm); ev = ev[ev > 0]
    rank = int(np.sum(ev > tol * ev.max())) if ev.size else 0
    return rank, ns
