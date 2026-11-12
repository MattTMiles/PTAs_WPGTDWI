#!/usr/bin/env python3
"""
VERNIER-30 -- integer ambiguity resolution (LAMBDA/ILS) on the fringe comb.

EXPLORATORY, NOT pre-registered.  cronus, CPU, numpy only.
Clean room: imports only vernier_toy (for the physical constants / fringe spacing
and the geometry draw).  No campaign code, no jax.

THE PROBLEM, STATED AS GNSS
---------------------------
One pulsar, M sources at distinct frequencies.  Write the distance offset from the
prior mean as x (pc).  Source k gives that pulsar a fringe spacing (carrier
wavelength)

    lam_k = c / (f_k (1 - cos mu_k))                        [vernier_toy.dL_of]

and an AMBIGUOUS pulsar-term phase.  In cycles:

    y_k = x / lam_k - n_k + eps_k ,   n_k in Z ,  eps_k ~ N(0, sig_k^2)   (cycles)
    0   = x        + nu           ,               nu    ~ N(0, sigL^2)    (pc)

y_k is the measured fractional pulsar-term phase (wrapped), n_k the fringe index,
sig_k = 1/(2 pi rho_k) with rho_k the per-(pulsar,source) SNR, sigL the distance
prior width.  This is EXACTLY single-baseline multi-frequency carrier-phase
ambiguity resolution: lam = carrier wavelength, the distance prior = the code /
pseudorange solution, the M sources = a multi-constellation GNSS receiver.

The system is EXACTLY DETERMINED (M phase obs + 1 prior vs 1 real + M integers),
so the float solution is closed form:

    xhat = 0 ,  nhat = -y ,
    Q_xx = sigL^2 ,  Q_xn[k] = sigL^2 / lam_k ,
    Q_nn[j,k] = sigL^2/(lam_j lam_k) + delta_jk sig_k^2      (rank-1 + diagonal)

The rank-1 part is the whole story: Q_nn is catastrophically ill-conditioned
(cond ~ (sigL/lam/sig_eps)^2 ~ 1e10), and the integer combinations that are
FIXABLE are exactly those nearly orthogonal to the 1/lam direction -- the GNSS
WIDE LANES.  That is the vernier, in Teunissen's coordinates.  Decorrelation is
therefore not cosmetic here; it is the method.

THREE ENGINES, RNG-MATCHED on identical realisations
----------------------------------------------------
 (a) SNAP    -- the node sampler's step: seed at the float (x=0), snap every
                fringe index to its MAP value, closed-form conditional x, iterate;
                then MODE-HOP proposals (shift one source's integer by +-1..H,
                re-snap, accept uphill).  This is snap_distances / MAP-fringe snap.
 (b) GRID    -- brute comb grid: dense sub-fringe scan of x over +-3 sigL, integers
                snapped and x re-solved in closed form at every node.  Sub-fringe
                sampling is DELIBERATE -- a fringe-CENTRE grid pins the sub-fringe
                offset u=0 and fabricates fringe evidence (the banked artefact).
                This engine is exact-by-exhaustion and is the ceiling.
 (c) ILS     -- LAMBDA: LDL' factorisation of Q_nn, integer Gauss + permutation
                decorrelation (Z-transform), depth-first search-and-shrink for the
                two best integer vectors, conditional back-substitution for x.
                Ported from the classic Teunissen/RTKLIB formulation.

All three minimise the SAME objective
    s(n) = (n - nhat)' Q_nn^-1 (n - nhat)  ==  -2 lnL_profiled(n)
so their solutions are directly comparable and dlnL = (s_2 - s_1)/2 exactly.

THE RATIO TEST vs OUR GATE
--------------------------
GNSS accepts a fix when R = s_2/s_1 > mu (ratio test).  The campaign accepts when
the fringe contrast dlnL(best vs runner-up) = (s_2 - s_1)/2 exceeds the trials bar
ln K, K = number of fringes in the prior window.  One is a RATIO of the two
residuals, the other a DIFFERENCE.  Experiment 4 asks whether they order the same
realisations the same way, and which is safer at matched acceptance.

Outputs: reports/vernier30.npz, reports/vernier30.log (tee'd by the caller).
Staged, never committed here.
"""
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vernier_toy as V                                     # noqa: E402

REPORTS = V.REPORTS

# ------------------------------------------------------------------ knobs
BAND_NHZ   = V.BAND_NHZ            # (3, 20) nHz campaign band
LOG10MC    = (8.5, 9.5)
SIGL_PC    = 100.0                 # default distance prior width (pc); campaign
                                   # tight-lane pulsars sit at 89-155 pc
RHO_DEF    = 4.0                   # default per-(pulsar,source) SNR
M_DEF      = 4
N_REAL     = 400
GRID_PER_FRINGE = 8                # sub-fringe grid density (NOT fringe centres)
GRID_CAP   = 4_000_000
HOPS       = 6                     # mode-hop radius per source for engine (a)
LOOPMAX    = 200_000               # ILS search iteration guard


# =================================================================== LAMBDA
def ldl(Q):
    """Q = L' diag(D) L with L unit lower triangular (Teunissen's LD form)."""
    n = Q.shape[0]
    A = np.array(Q, dtype=float, copy=True)
    L = np.zeros((n, n))
    D = np.zeros(n)
    for i in range(n - 1, -1, -1):
        D[i] = A[i, i]
        if D[i] <= 0.0:
            raise np.linalg.LinAlgError("Q not positive definite in LDL")
        L[i, :i + 1] = A[i, :i + 1] / np.sqrt(A[i, i])
        for j in range(i):
            A[j, :j + 1] -= L[i, :j + 1] * L[i, j]
        L[i, :i + 1] /= L[i, i]
    return L, D


def _gauss(L, Z, i, j):
    mu = int(np.rint(L[i, j]))
    if mu != 0:
        L[i:, j] -= mu * L[i:, i]
        Z[:, j] -= mu * Z[:, i]


def _perm(L, D, j, delta, Z):
    n = len(D)
    eta = D[j] / delta
    lam = D[j + 1] * L[j + 1, j] / delta
    D[j] = eta * D[j + 1]
    D[j + 1] = delta
    if j > 0:
        a0 = L[j, :j].copy()
        a1 = L[j + 1, :j].copy()
        L[j, :j] = -L[j + 1, j] * a0 + a1
        L[j + 1, :j] = eta * a0 + lam * a1
    L[j + 1, j] = lam
    if j + 2 < n:
        L[j + 2:, [j, j + 1]] = L[j + 2:, [j + 1, j]]
    Z[:, [j, j + 1]] = Z[:, [j + 1, j]]


def reduction(L, D):
    """Integer Gauss + permutation decorrelation.  Returns (Z, L, D) in place."""
    n = len(D)
    Z = np.eye(n)
    j = k = n - 2
    while j >= 0:
        if j <= k:
            for i in range(j + 1, n):
                _gauss(L, Z, i, j)
        delta = D[j] + L[j + 1, j] ** 2 * D[j + 1]
        if delta + 1e-12 < D[j + 1]:
            _perm(L, D, j, delta, Z)
            k = j
            j = n - 2
        else:
            j -= 1
    return Z, L, D


def ils_search(zs, L, D, m=2, loopmax=LOOPMAX):
    """Depth-first search-and-shrink over the decorrelated lattice.

    Returns (zn (n,m), s (m,), n_visited).  Candidates sorted by s ascending.
    """
    n = len(zs)
    S = np.zeros((n, n))
    dist = np.zeros(n)
    zb = np.zeros(n)
    z = np.zeros(n)
    step = np.zeros(n)
    zn = np.zeros((n, m))
    s = np.zeros(m)

    def sgn(v):
        return 1.0 if v >= 0.0 else -1.0

    k = n - 1
    dist[k] = 0.0
    zb[k] = zs[k]
    z[k] = np.rint(zb[k])
    y = zb[k] - z[k]
    step[k] = sgn(y)

    maxdist = 1e99
    nn = 0
    imax = 0
    visited = 0
    for c in range(loopmax):
        visited += 1
        newdist = dist[k] + y * y / D[k]
        if newdist < maxdist:
            if k != 0:
                k -= 1
                dist[k] = newdist
                S[k, :k + 1] = S[k + 1, :k + 1] + (z[k + 1] - zb[k + 1]) * L[k + 1, :k + 1]
                zb[k] = zs[k] + S[k, k]
                z[k] = np.rint(zb[k])
                y = zb[k] - z[k]
                step[k] = sgn(y)
            else:
                if nn < m:
                    if nn == 0 or newdist > s[imax]:
                        imax = nn
                    zn[:, nn] = z
                    s[nn] = newdist
                    nn += 1
                else:
                    if newdist < s[imax]:
                        zn[:, imax] = z
                        s[imax] = newdist
                        imax = int(np.argmax(s))
                    maxdist = s[imax]
                z[0] += step[0]
                y = zb[0] - z[0]
                step[0] = -step[0] - sgn(step[0])
        else:
            if k == n - 1:
                break
            k += 1
            z[k] += step[k]
            y = zb[k] - z[k]
            step[k] = -step[k] - sgn(step[k])

    order = np.argsort(s[:nn])
    return zn[:, :nn][:, order], s[:nn][order], visited


def lambda_fix(nhat, Qnn, m=2):
    """Full LAMBDA: decorrelate, ILS, back-transform.  Returns (nfix (M,m), s, nvis)."""
    L, D = ldl(Qnn)
    Z, L, D = reduction(L, D)
    zs = Z.T @ nhat
    zn, s, nvis = ils_search(zs, L, D, m=m)
    # a = Z^-T z   (Z is unimodular integer)
    Zinv_T = np.linalg.inv(Z.T)
    nfix = np.rint(Zinv_T @ zn)
    return nfix, s, nvis


# ================================================== the toy: one realisation
def draw_realisation(rng, M, rho, sigL, band=BAND_NHZ, null=False):
    """One pulsar, M sources.  Returns a dict with lam, y, sig, truth."""
    p_hat = V.iso_dirs(1, rng)[0]
    f0 = 10.0 ** rng.uniform(np.log10(band[0] * 1e-9), np.log10(band[1] * 1e-9), M)
    s_hat = V.iso_dirs(M, rng)
    omc = np.maximum(1.0 - s_hat @ p_hat, 1e-2)
    lam = V.dL_of(f0, omc)                                   # pc, per source
    rho_k = np.atleast_1d(rho)
    if rho_k.size == 1:
        rho_k = np.full(M, float(rho))
    sig = 1.0 / (2.0 * np.pi * rho_k)                        # cycles

    x_true = rng.normal(0.0, sigL)                           # true offset from prior mean
    if null:
        # NO coherent comb: fractional phases are uniform, nothing to fix.
        y = rng.uniform(-0.5, 0.5, M)
        n_true = np.full(M, np.nan)
    else:
        phase = x_true / lam + rng.normal(0.0, sig)
        y = phase - np.rint(phase)                           # wrapped fractional cycles
        n_true = np.rint(x_true / lam - y)                   # the integer consistent with y
    return dict(lam=lam, y=y, sig=sig, sigL=sigL, x_true=x_true, n_true=n_true,
                M=M, rho=rho_k)


def float_solution(R):
    """Closed-form float solution + covariance (exactly determined system)."""
    lam, sig, sigL, y = R["lam"], R["sig"], R["sigL"], R["y"]
    xhat = 0.0
    nhat = -y
    Qnn = sigL ** 2 * np.outer(1.0 / lam, 1.0 / lam) + np.diag(sig ** 2)
    Qxn = sigL ** 2 / lam
    return xhat, nhat, Qnn, Qxn, sigL ** 2


def obj(R, n):
    """s(n) = min_x [ sum_k (x/lam_k - y_k - n_k)^2/sig_k^2 + x^2/sigL^2 ].

    Vectorised over leading axes of n.  Returns (s, x_hat_conditional).
    """
    lam, sig, sigL, y = R["lam"], R["sig"], R["sigL"], R["y"]
    w = 1.0 / sig ** 2
    b = (y + n) @ (w / lam)
    a = np.sum(w / lam ** 2) + 1.0 / sigL ** 2
    x = b / a
    r = x[..., None] / lam - y - n
    s = np.sum(r ** 2 * w, axis=-1) + x ** 2 / sigL ** 2
    return s, x


# ============================================================ engine (a) SNAP
def engine_snap(R, hops=HOPS, max_iter=40):
    """Node-sampler step: MAP-fringe snap to convergence, then mode-hop."""
    lam, y = R["lam"], R["y"]
    x = 0.0                                                  # the float landing
    n = np.rint(x / lam - y)
    for _ in range(max_iter):
        s, x = obj(R, n)
        n_new = np.rint(x / lam - y)
        if np.array_equal(n_new, n):
            break
        n = n_new
    s_best, x_best = obj(R, n)
    n_best = n.copy()
    s_run = np.inf                                           # runner-up for the gate
    for k in range(R["M"]):
        for d in range(-hops, hops + 1):
            if d == 0:
                continue
            n_try = n_best.copy()
            n_try[k] += d
            # re-snap the others
            for _ in range(max_iter):
                s_t, x_t = obj(R, n_try)
                n_re = np.rint(x_t / lam - y)
                n_re[k] = n_try[k]
                if np.array_equal(n_re, n_try):
                    break
                n_try = n_re
            s_t, x_t = obj(R, n_try)
            if s_t < s_best:
                s_run, = (s_best,)
                s_best, x_best, n_best = s_t, x_t, n_try.copy()
            elif s_t < s_run and not np.array_equal(n_try, n_best):
                s_run = s_t
    return n_best, x_best, s_best, s_run


def engine_snap_restart(R, n_restart=64, hops=2, max_iter=40):
    """Fair-budget SNAP: n_restart seeds drawn from the PRIOR (what a sampler
    chain explores over time), each snapped to its fixed point, best kept.
    The bare single-seed engine_snap is the one-shot node-sampler step; this is
    the same move given a search budget comparable to the other engines."""
    lam, y, sigL = R["lam"], R["y"], R["sigL"]
    rng = np.random.default_rng(abs(int(1e6 * y[0])) % (2 ** 31))
    seeds = np.concatenate([[0.0], rng.normal(0.0, sigL, n_restart - 1)])
    n_all = np.rint(seeds[:, None] / lam - y)
    for _ in range(max_iter):
        s, xs = obj(R, n_all)
        n_new = np.rint(xs[:, None] / lam - y)
        if np.array_equal(n_new, n_all):
            break
        n_all = n_new
    s, xs = obj(R, n_all)
    i1 = int(np.argmin(s))
    n1 = n_all[i1]
    diff = np.any(n_all != n1, axis=1)
    s2 = float(np.min(s[diff])) if diff.any() else np.inf
    return n1, float(xs[i1]), float(s[i1]), s2, n_restart


# ============================================================ engine (b) GRID
def engine_grid(R, per_fringe=GRID_PER_FRINGE, cap=GRID_CAP):
    """Brute comb grid, sub-fringe sampled, integers snapped + x re-solved."""
    lam, sigL, y = R["lam"], R["sigL"], R["y"]
    half = 3.0 * sigL
    step = lam.min() / per_fringe
    npts = int(min(2.0 * half / step, cap)) + 1
    xg = np.linspace(-half, half, npts)
    ng = np.rint(xg[:, None] / lam - y)
    s, xs = obj(R, ng)
    i1 = int(np.argmin(s))
    n1 = ng[i1]
    # runner-up = best node whose integer vector differs from the winner
    diff = np.any(ng != n1, axis=1)
    s2 = float(np.min(s[diff])) if diff.any() else np.inf
    return n1, float(xs[i1]), float(s[i1]), s2, npts


# ============================================================= engine (c) ILS
def engine_ils(R, m=2):
    xhat, nhat, Qnn, Qxn, Qxx = float_solution(R)
    nfix, s_lam, nvis = lambda_fix(nhat, Qnn, m=m)
    # rescore on the true objective (identical up to the exactly-determined float,
    # but recomputing keeps the three engines on one scale by construction)
    s = np.array([obj(R, nfix[:, j])[0] for j in range(nfix.shape[1])])
    order = np.argsort(s)
    nfix, s = nfix[:, order], s[order]
    s1, xh = obj(R, nfix[:, 0])
    s2 = float(s[1]) if s.size > 1 else np.inf
    return nfix[:, 0], float(xh), float(s1), s2, nvis


# ================================================================ experiments
ENGINES = ("snap", "snapR", "grid", "ils")


def run_cell(rng_seed, M, rho, sigL, n_real, engines=ENGINES,
             null=False, grid_per_fringe=GRID_PER_FRINGE):
    """RNG-MATCHED: one realisation stream, every engine sees identical data."""
    rng = np.random.default_rng(rng_seed)
    out = {e: dict(succ=[], s1=[], s2=[], t=0.0, extra=[]) for e in engines}
    K_tr = []
    for _ in range(n_real):
        R = draw_realisation(rng, M, rho, sigL, null=null)
        K_tr.append(6.0 * sigL / R["lam"].min())
        for e in engines:
            t0 = time.perf_counter()
            if e == "snap":
                n, x, s1, s2 = engine_snap(R)
                ex = 0
            elif e == "snapR":
                n, x, s1, s2, ex = engine_snap_restart(R)
            elif e == "grid":
                n, x, s1, s2, ex = engine_grid(R, per_fringe=grid_per_fringe)
            else:
                n, x, s1, s2, ex = engine_ils(R)
            out[e]["t"] += time.perf_counter() - t0
            ok = (not null) and np.array_equal(n, R["n_true"])
            out[e]["succ"].append(bool(ok))
            out[e]["s1"].append(float(s1))
            out[e]["s2"].append(float(s2))
            out[e]["extra"].append(ex)
    for e in engines:
        for kk in ("succ", "s1", "s2", "extra"):
            out[e][kk] = np.asarray(out[e][kk])
        out[e]["t_per"] = out[e]["t"] / n_real
    return out, np.asarray(K_tr)


def fmt_cell(tag, out, K_tr, engines=ENGINES):
    line = f"  {tag:>22s} |"
    for e in engines:
        line += f" {e}:{100*out[e]['succ'].mean():5.1f}%/{1e3*out[e]['t_per']:7.3f}ms |"
    return line


# --------------------------------------------------------------- gate compare
def gate_stats(out, K_tr, engine="ils"):
    """Ratio test (s2/s1) vs contrast gate ((s2-s1)/2 vs ln K)."""
    s1 = out[engine]["s1"]
    s2 = out[engine]["s2"]
    succ = out[engine]["succ"]
    ratio = np.where(s1 > 0, s2 / np.maximum(s1, 1e-300), np.inf)
    contrast = 0.5 * (s2 - s1)
    bar = np.log(np.maximum(K_tr, 2.0))
    return ratio, contrast, bar, succ


def roc_at_matched(ratio, contrast, bar, succ, accept_rates):
    """Wrong-fix rate of each criterion at MATCHED acceptance rate.

    Gate statistic is the MARGIN (contrast - ln K), i.e. the gate's own quantity;
    sliding its threshold off zero is what makes the acceptance rates matchable.
    """
    margin = contrast - bar
    rows = []
    for ar in accept_rates:
        if ar <= 0 or ar >= 1:
            continue
        q = 1.0 - ar
        thr_r = np.quantile(ratio[np.isfinite(ratio)], q) if np.isfinite(ratio).any() else np.inf
        thr_c = np.quantile(margin, q)
        acc_r = ratio >= thr_r
        acc_c = margin >= thr_c
        wf_r = (acc_r & ~succ).sum() / max(acc_r.sum(), 1)
        wf_c = (acc_c & ~succ).sum() / max(acc_c.sum(), 1)
        agree = (acc_r == acc_c).mean()
        rows.append((ar, acc_r.mean(), wf_r, acc_c.mean(), wf_c, agree,
                     float(thr_r), float(thr_c)))
    return np.array(rows)


# ======================================================================== main
def main():
    t_all = time.perf_counter()
    out = {}
    print("=" * 78)
    print("VERNIER-30  --  integer ambiguity resolution (LAMBDA/ILS) on the comb")
    print("=" * 78)

    # ---------------------------------------------------------- SANITY: LAMBDA
    print("\n--- S0. LAMBDA machinery sanity ------------------------------------")
    rs = np.random.default_rng(7)
    worst = 0.0
    condQ = []
    for _ in range(200):
        R = draw_realisation(rs, 5, 4.0, SIGL_PC)
        _, nhat, Qnn, _, _ = float_solution(R)
        condQ.append(np.linalg.cond(Qnn))
        L, D = ldl(Qnn)
        rec = L.T @ np.diag(D) @ L
        worst = max(worst, np.max(np.abs(rec - Qnn)) / np.max(np.abs(Qnn)))
        Lc, Dc = ldl(Qnn)
        Z, Lz, Dz = reduction(Lc, Dc)
        detZ = abs(np.linalg.det(Z))
        assert abs(detZ - 1.0) < 1e-6, f"Z not unimodular: det={detZ}"
        Qz = Z.T @ Qnn @ Z
        recz = Lz.T @ np.diag(Dz) @ Lz
        worst = max(worst, np.max(np.abs(recz - Qz)) / np.max(np.abs(Qz)))
    print(f"  LDL' + Z-transform round-trip max rel err : {worst:.2e}   (200 draws, M=5)")
    print(f"  cond(Q_nn) median                         : {np.median(condQ):.3e}")
    print(f"  -> Q_nn is rank-1-dominant; decorrelation is the method, not cosmetics.")
    out["sanity_ldl_err"] = worst
    out["sanity_condQ"] = np.asarray(condQ)

    # ILS vs the exhaustive engine.  The reference is the DENSE sub-fringe grid
    # over the whole +-3 sigL prior window (engine b) at 64 pts/fringe, which
    # enumerates every integer vector the prior admits.  A +-R per-component box
    # around nhat is NOT a valid reference: Q_nn is rank-1 dominant, so the
    # admissible integer vectors run along a lattice line of length ~sigL/lam
    # (hundreds of units), and any small box misses the global optimum.
    rs = np.random.default_rng(11)
    same, ils_le, ils_lt, grid_lt = 0, 0, 0, 0
    NS = 300
    for _ in range(NS):
        R = draw_realisation(rs, 4, 4.0, 10.0)
        n_i, _, s_i, _, _ = engine_ils(R)
        n_g, _, s_g, _, _ = engine_grid(R, per_fringe=64)
        if np.array_equal(n_i, n_g):
            same += 1
        if s_i <= s_g + 1e-9:
            ils_le += 1
        if s_i < s_g - 1e-9:
            ils_lt += 1
        if s_g < s_i - 1e-9:
            grid_lt += 1
    print(f"  ILS vs dense-grid exhaustive (M=4, sigL=10pc, 64 pts/fringe):")
    print(f"    identical integer vector                : {same}/{NS}")
    print(f"    ILS objective <= grid objective         : {ils_le}/{NS}")
    print(f"    ILS strictly BETTER than grid           : {ils_lt}/{NS}  "
          f"(grid missed the optimum)")
    print(f"    grid strictly BETTER than ILS           : {grid_lt}/{NS}  "
          f"(must be 0 -- ILS is exact)")
    out["sanity_ils_grid"] = np.array([same, ils_le, ils_lt, grid_lt, NS])

    # ------------------------------------------------------------- EXP 1: K
    print("\n--- E1. success / wall-clock vs K = number of sources ---------------")
    print(f"    (rho={RHO_DEF}, sigL={SIGL_PC} pc, {N_REAL} realisations, RNG-matched)")
    Ks = list(range(2, 9))
    e1 = {}
    for M in Ks:
        o, K_tr = run_cell(1000 + M, M, RHO_DEF, SIGL_PC, N_REAL)
        e1[M] = (o, K_tr)
        print(fmt_cell(f"K={M} (lnK={np.log(K_tr.mean()):.1f})", o, K_tr))
    for e in ENGINES:
        out[f"e1_{e}_succ"] = np.array([e1[M][0][e]["succ"].mean() for M in Ks])
        out[f"e1_{e}_t"] = np.array([e1[M][0][e]["t_per"] for M in Ks])
    out["e1_K"] = np.array(Ks)

    # ------------------------------------------------------------ EXP 2: SNR
    print("\n--- E2. success vs per-(pulsar,source) SNR rho ----------------------")
    print(f"    (K={M_DEF}, sigL={SIGL_PC} pc, {N_REAL} realisations)")
    rhos = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
    e2 = {}
    for r in rhos:
        o, K_tr = run_cell(2000 + int(10 * r), M_DEF, r, SIGL_PC, N_REAL)
        e2[r] = (o, K_tr)
        print(fmt_cell(f"rho={r:g}", o, K_tr))
    for e in ENGINES:
        out[f"e2_{e}_succ"] = np.array([e2[r][0][e]["succ"].mean() for r in rhos])
        out[f"e2_{e}_t"] = np.array([e2[r][0][e]["t_per"] for r in rhos])
    out["e2_rho"] = np.array(rhos)

    # ---------------------------------------------------- EXP 3: prior width
    print("\n--- E3. success / wall-clock vs prior width sigL --------------------")
    print(f"    (K={M_DEF}, rho={RHO_DEF}, {N_REAL} realisations)")
    sigs = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0]
    e3 = {}
    for sg in sigs:
        nr = N_REAL if sg <= 300.0 else (100 if sg <= 1000.0 else 40)
        o, K_tr = run_cell(3000 + int(sg), M_DEF, RHO_DEF, sg, nr)
        e3[sg] = (o, K_tr)
        print(fmt_cell(f"sigL={sg:g}pc (lnK={np.log(K_tr.mean()):.1f})", o, K_tr))
    for e in ENGINES:
        out[f"e3_{e}_succ"] = np.array([e3[sg][0][e]["succ"].mean() for sg in sigs])
        out[f"e3_{e}_t"] = np.array([e3[sg][0][e]["t_per"] for sg in sigs])
    out["e3_sigL"] = np.array(sigs)
    out["e3_gridpts"] = np.array([e3[sg][0]["grid"]["extra"].mean() for sg in sigs])
    out["e3_ilsvisit"] = np.array([e3[sg][0]["ils"]["extra"].mean() for sg in sigs])

    # ------------------------------------------- EXP 4: ratio test vs our gate
    print("\n--- E4. ratio test vs contrast-vs-bar certification gate ------------")
    pool_r, pool_c, pool_b, pool_s, pool_s1 = [], [], [], [], []
    for M in (3, 4, 5, 6):
        for r in (2.0, 4.0, 8.0):
            for sg in (10.0, 100.0, 1000.0):        # vary sigL so ln K VARIES:
                                                    # the gate's per-realisation
                                                    # bar must do real work here
                o, K_tr = run_cell(4000 + 100 * M + 10 * int(r) + int(sg), M, r, sg,
                                   120, engines=("ils",))
                ratio, contrast, bar, succ = gate_stats(o, K_tr, "ils")
                pool_r.append(ratio); pool_c.append(contrast)
                pool_b.append(bar); pool_s.append(succ)
                pool_s1.append(o["ils"]["s1"])
    ratio = np.concatenate(pool_r); contrast = np.concatenate(pool_c)
    bar = np.concatenate(pool_b); succ = np.concatenate(pool_s)
    s1p = np.concatenate(pool_s1)
    fin = np.isfinite(ratio) & np.isfinite(contrast)
    ratio, contrast, bar, succ, s1p = (ratio[fin], contrast[fin], bar[fin],
                                       succ[fin], s1p[fin])
    print(f"    pooled n = {succ.size}, overall ILS success {100*succ.mean():.1f}%")
    print(f"    trials bar ln K spans {bar.min():.2f} - {bar.max():.2f} "
          f"(sd {bar.std():.2f}) -- the bar is not a constant in this pool")

    # native operating points
    acc_gate = contrast >= bar
    print(f"    NATIVE gate  (contrast >= ln K)  : accept {100*acc_gate.mean():5.1f}%  "
          f"wrong-fix {100*((acc_gate & ~succ).sum()/max(acc_gate.sum(),1)):5.2f}%")
    for mu in (1.5, 2.0, 3.0):
        acc_r = ratio >= mu
        print(f"    NATIVE ratio (R >= {mu:.1f})           : accept {100*acc_r.mean():5.1f}%  "
              f"wrong-fix {100*((acc_r & ~succ).sum()/max(acc_r.sum(),1)):5.2f}%  "
              f"agree-with-gate {100*(acc_r == acc_gate).mean():5.1f}%")

    rows = roc_at_matched(ratio, contrast, bar, succ,
                          [0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90])
    print("\n    MATCHED-ACCEPTANCE comparison (which criterion is safer):")
    print("      target  ratio-acc  ratio-WF  gate-acc  gate-WF   agree")
    for a, ar_, wr, ac_, wc, ag, tr, tc in rows:
        print(f"      {100*a:5.0f}%   {100*ar_:7.1f}%  {100*wr:7.2f}%  "
              f"{100*ac_:7.1f}%  {100*wc:7.2f}%  {100*ag:6.1f}%")
    out["e4_ratio"] = ratio; out["e4_contrast"] = contrast
    out["e4_bar"] = bar; out["e4_succ"] = succ; out["e4_matched"] = rows

    # how much can each criterion accept while staying under 1% / 5% wrong-fix?
    print("\n    Max acceptance rate holding the wrong-fix rate under a target:")
    margin = contrast - bar
    for wf_tgt in (0.01, 0.05):
        best = {}
        for name, stat in (("ratio", ratio), ("gate-margin", margin)):
            ok = 0.0
            for ar in np.linspace(0.01, 0.99, 99):
                thr = np.quantile(stat[np.isfinite(stat)], 1.0 - ar)
                acc = stat >= thr
                wf = (acc & ~succ).sum() / max(acc.sum(), 1)
                if wf <= wf_tgt:
                    ok = ar
            best[name] = ok
        rel = (f"{best['gate-margin']/best['ratio']:.1f}x" if best["ratio"] > 0
               else "ratio never reaches this target")
        print(f"      WF <= {100*wf_tgt:4.1f}%  ->  ratio {100*best['ratio']:5.1f}%   "
              f"gate {100*best['gate-margin']:5.1f}%   ({rel})")
        out[f"e4_maxacc_wf{int(100*wf_tgt)}"] = np.array(
            [best["ratio"], best["gate-margin"]])

    # WHY the ratio test is the less safe statistic here: s_1 sits in its
    # denominator, so a realisation whose best residual is accidentally tiny
    # gets a huge R with no real evidence behind it.
    acc_r3 = ratio >= 3.0
    bad = acc_r3 & ~succ
    if bad.any() and (acc_r3 & succ).any():
        print(f"\n    mechanism: among R>=3 accepts, median s_1 = "
              f"{np.median(s1p[acc_r3 & succ]):.3f} (correct fixes) vs "
              f"{np.median(s1p[bad]):.3f} (WRONG fixes); "
              f"median contrast {np.median(contrast[acc_r3 & succ]):.2f} vs "
              f"{np.median(contrast[bad]):.2f} nat")
    out["e4_s1"] = s1p

    # rank correlation of the two statistics
    from scipy.stats import spearmanr, kendalltau
    rho_s = spearmanr(ratio, contrast).statistic
    tau_k = kendalltau(ratio, contrast).statistic
    print(f"\n    Spearman(ratio, contrast) = {rho_s:.4f}   Kendall tau = {tau_k:.4f}")
    out["e4_spearman"] = rho_s; out["e4_kendall"] = tau_k

    # ------------------------------------------------- EXP 5: null / FPR
    print("\n--- E5. NULL control (no coherent comb): false-acceptance -----------")
    o, K_tr = run_cell(5150, M_DEF, RHO_DEF, SIGL_PC, 400, engines=("ils",), null=True)
    rN, cN, bN, _ = gate_stats(o, K_tr, "ils")
    finN = np.isfinite(rN) & np.isfinite(cN)
    rN, cN, bN = rN[finN], cN[finN], bN[finN]
    print(f"    gate  (contrast >= ln K) fires on pure noise : "
          f"{100*(cN >= bN).mean():5.2f}%")
    for mu in (1.5, 2.0, 3.0):
        print(f"    ratio (R >= {mu:.1f})           fires on pure noise : "
              f"{100*(rN >= mu).mean():5.2f}%")
    out["e5_null_ratio"] = rN; out["e5_null_contrast"] = cN; out["e5_null_bar"] = bN

    # --------------------------------------------------- EXP 6: PAR (bonus)
    print("\n--- E6. partial ambiguity resolution: are wide lanes fixable? -------")
    rs = np.random.default_rng(6060)
    fixable = []
    for _ in range(300):
        R = draw_realisation(rs, 5, RHO_DEF, SIGL_PC)
        _, nhat, Qnn, _, _ = float_solution(R)
        L, D = ldl(Qnn)
        Z, Lz, Dz = reduction(L, D)
        # bootstrapped success rate per decorrelated component (Teunissen)
        from scipy.special import erf
        P = erf(1.0 / (2.0 * np.sqrt(2.0 * Dz)))
        fixable.append(P)
    Pb = np.array(fixable)
    print("    bootstrapped per-component fix probability, decorrelated (M=5):")
    print("      component :   " + "  ".join(f"z{j+1}" for j in range(5)))
    print("      median P  :   " + "  ".join(f"{v:.3f}" for v in np.median(Pb, 0)))
    print(f"    full-vector bootstrap success (product) median: "
          f"{np.median(np.prod(Pb, 1)):.4f}")
    print(f"    >=1 component with P>0.99: {100*(Pb > 0.99).any(1).mean():.1f}% of draws; "
          f"mean count {np.mean((Pb > 0.99).sum(1)):.2f}/5")
    o6, _ = run_cell(6061, 5, RHO_DEF, SIGL_PC, 300, engines=("ils",))
    print(f"    Teunissen bootstrap PREDICTS full-fix success "
          f"{np.median(np.prod(Pb, 1)):.3f}; MEASURED ILS success on a matched "
          f"cell = {o6['ils']['succ'].mean():.3f}")
    out["e6_Pboot"] = Pb
    out["e6_measured"] = o6["ils"]["succ"].mean()

    out["runtime_s"] = time.perf_counter() - t_all
    np.savez_compressed(os.path.join(REPORTS, "vernier30.npz"), **out)
    print(f"\ntotal runtime {out['runtime_s']:.1f} s -> reports/vernier30.npz")


if __name__ == "__main__":
    main()
