"""SIEVE-C T4 -- CLOSURE OBSERVABLES.  Zero-noise / Asimov only.

Toy: 6 pulsars x 3 sources (top-6 census pulsars, real positions/distances).

PART 1 -- the algebra, verified numerically to machine precision.
  Observable: the per-(pulsar, source) PULSAR-TERM phase relative to the Earth
  term.  The measured phase carries two nuisance classes the brief names:

      phi_pk = g_pk(sky, L_p, f_k)  +  Phi_k  +  c_p   [ + 2 pi f_k tau_p ]

    g_pk  = -2 pi f_k L_p (1 - shat_k . phat_p) / c    (the SIGNAL geometry)
    Phi_k = per-SOURCE phase offset   (unknown initial orbital phase)
    c_p   = per-PULSAR ACHROMATIC phase offset (a "clock-like" term that is the
            same phase for every source)
    tau_p = per-PULSAR TIME offset (a genuine clock error) -> phase 2 pi f_k tau_p,
            which is CHROMATIC (grows with f_k).  Kept separate on purpose.

  Closure classes:
    C1 source-pair  (fixed p) : phi_pk - phi_pk'                 kills c_p
    C2 pulsar-pair  (fixed k) : phi_pk - phi_p'k                 kills Phi_k
    C3 double       (closure) : phi_pk - phi_pk' - phi_p'k + phi_p'k'
                                                                 kills BOTH
    C4 chromatic double       : the C3 combination formed on phi/f rather than
                                phi -- kills Phi_k/f_k and tau_p, NOT c_p.

PART 2 -- THE MONEY TEST.  Inject the WRONG-SKY manufacture mode (STOP-#2 class:
  a bright template at a SCRAMBLED sky) and ask whether a closure-consistency
  statistic REJECTS it where the raw statistic accepts.  Both statistics are
  given the SAME freedom (all 6 pulsar distances free, all nuisances free) so
  the comparison is apples-to-apples: the only difference is that the closure
  statistic is nuisance-free BY CONSTRUCTION, whereas the raw statistic has to
  spend its own freedom absorbing Phi_k and c_p.
"""

import os
import sys
import json
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sieve_toy as S

OUT = os.path.join(S.ROOT, "reports", "sieve_t4_closure.npz")
NPSR, NSRC = 6, 3
RNG_ALGEBRA = 20260729          # deterministic nuisance draws for the ALGEBRA gate only
N_SCRAM = 16                    # deterministic scrambled-sky ensemble


def make_sources():
    """3 sources, band-spread in f, distinct skies.  Deterministic."""
    return [
        dict(costheta=+0.35, phi=1.90, log10_h=-13.2, log10_f=-8.00,
             cosi=+0.40, psi=0.70, phi0=1.10),
        dict(costheta=-0.55, phi=4.10, log10_h=-13.4, log10_f=-7.82,
             cosi=-0.30, psi=2.05, phi0=3.60),
        dict(costheta=+0.80, phi=0.55, log10_h=-13.5, log10_f=-7.66,
             cosi=+0.75, psi=1.30, phi0=5.20),
    ]


def geom_phase(spine, src, L_m=None):
    """g_pk = -2 pi f_k L_p (1 - shat_k . phat_p)/c   (radians, unwrapped)."""
    shat = S.sky_unit(src["costheta"], src["phi"])
    f = 10.0 ** src["log10_f"]
    L = spine["L_m"] if L_m is None else L_m
    return -2.0 * np.pi * f * L * (1.0 - spine["phat"] @ shat) / S.C_LIGHT


def phase_matrix(spine, srcs, L_m=None, Phi=None, c=None, tau=None):
    """(npsr, nsrc) observable phase matrix with the nuisances applied."""
    n_p, n_s = spine["npsr"], len(srcs)
    P = np.zeros((n_p, n_s))
    for k, s in enumerate(srcs):
        P[:, k] = geom_phase(spine, s, L_m=L_m)
        if Phi is not None:
            P[:, k] += Phi[k]
        if c is not None:
            P[:, k] += c
        if tau is not None:
            P[:, k] += 2.0 * np.pi * (10.0 ** s["log10_f"]) * tau
    return P


# ---- closure operators -------------------------------------------------
def C1(P):
    """source-pair differences, (npsr, nsrc-1): phi_pk - phi_p,k+1"""
    return P[:, :-1] - P[:, 1:]


def C2(P):
    """pulsar-pair differences, (npsr-1, nsrc): phi_pk - phi_p+1,k"""
    return P[:-1, :] - P[1:, :]


def C3(P):
    """double differences (the closure), (npsr-1, nsrc-1)"""
    return P[:-1, :-1] - P[:-1, 1:] - P[1:, :-1] + P[1:, 1:]


def C4(P, fk):
    """the chromatic closure: C3 applied to phi/f (kills tau_p and Phi_k/f_k)"""
    return C3(P / fk[None, :])


def C3b(P):
    """batched C3 -- operates on the LAST two axes, so it accepts (..., np, ns)."""
    return (P[..., :-1, :-1] - P[..., :-1, 1:]
            - P[..., 1:, :-1] + P[..., 1:, 1:])


def main():
    spine = S.build_spine(npsr=NPSR, T_yr=18.0, cadence_per_yr=20.0)
    srcs = make_sources()
    fk = np.array([10.0 ** s["log10_f"] for s in srcs])
    rng = np.random.default_rng(RNG_ALGEBRA)

    print(f"[T4] pulsars: {spine['names']}")
    print(f"[T4] f_k = {fk}")

    # ==================================================================
    # PART 1 -- ALGEBRA GATE
    # ==================================================================
    print("\n[T4-P1] which nuisance does each closure class cancel?")
    P0 = phase_matrix(spine, srcs)
    Phi = rng.uniform(0, 2 * np.pi, NSRC)
    cp = rng.uniform(0, 2 * np.pi, NPSR)
    tau = rng.uniform(-1.0, 1.0, NPSR) * 1e7          # s, a clock-like TIME offset

    cases = {
        "Phi_k only": dict(Phi=Phi),
        "c_p only": dict(c=cp),
        "Phi_k + c_p": dict(Phi=Phi, c=cp),
        "tau_p only": dict(tau=tau),
        "Phi_k + c_p + tau_p": dict(Phi=Phi, c=cp, tau=tau),
    }
    algebra = {}
    hdr = f"  {'nuisance injected':22s} {'C1 src-pair':>13s} {'C2 psr-pair':>13s} " \
          f"{'C3 closure':>13s} {'C4 chromatic':>14s}"
    print(hdr)
    for nm, kw in cases.items():
        P = phase_matrix(spine, srcs, **kw)
        d = [float(np.abs(C1(P) - C1(P0)).max()),
             float(np.abs(C2(P) - C2(P0)).max()),
             float(np.abs(C3(P) - C3(P0)).max()),
             float(np.abs(C4(P, fk) - C4(P0, fk)).max())]
        algebra[nm] = d
        print(f"  {nm:22s} " + " ".join(f"{v:13.3e}" for v in d))

    # ==================================================================
    # PART 2 -- THE MONEY TEST: the WRONG-SKY manufacture mode
    # ==================================================================
    # The observables are phases, so they are known only MOD 2 pi -- that is what
    # a fringe measurement IS, and it is the whole reason a wrong sky can be
    # absorbed by re-hopping pulsar distances.  Everything below therefore runs
    # in the COMPLEX (wrapped) domain, z_pk = exp(i phi_pk), with the pulsar
    # distances free over their FULL canonical +-3 sigma prior range (which spans
    # 1 to ~1270 fringes for these pulsars -- so fringe hopping is ALLOWED, as it
    # is in the real problem).
    print("\n[T4-P2] money test: does closure REJECT a bright template at a "
          "scrambled sky?")
    rho = np.zeros((NPSR, NSRC))
    for k, s in enumerate(srcs):
        r = S.residuals(spine, s, front="plane")
        rho[:, k] = np.sqrt(np.sum((r / spine["sigma"][:, None]) ** 2, axis=1))
    W = np.clip(rho, 1e-6, None) ** 2          # weight ~ 1/sigma_phi^2 = rho^2
    print(f"  per-(p,k) SNR {rho.min():.1f}-{rho.max():.1f}")
    nfr = 3.0 * spine["Lerr_m"] / (S.C_LIGHT / fk.max())
    print(f"  +-3sigma distance prior spans {nfr.min():.1f}-{nfr.max():.1f} "
          f"fringes per pulsar -- fringe hopping ALLOWED")

    def coef(srcs_try):
        """d(phase)/d(L_p) per (p,k):  a_pk = -2 pi f_k (1 - cos mu_pk)/c."""
        A = np.zeros((NPSR, NSRC))
        for k, s in enumerate(srcs_try):
            shat = S.sky_unit(s["costheta"], s["phi"])
            A[:, k] = -2.0 * np.pi * (10.0 ** s["log10_f"]) * \
                (1.0 - spine["phat"] @ shat) / S.C_LIGHT
        return A

    def wrap(x):
        return (x + np.pi) % (2.0 * np.pi) - np.pi

    def stat_batch(Pm, P_obs, statistic):
        """Statistic for a BATCH of model phase matrices Pm (nb, NPSR, NSRC).

        raw_blind : wrapped residual, only the per-SOURCE phase Phi_k profiled
                    out (the pipeline's actual freedom -- it has no c_p term)
        raw_nuis  : both Phi_k and c_p profiled out (the fair comparison)
        closure   : the C3 double difference -- immune to BOTH by construction
        """
        R = wrap(P_obs[None, :, :] - Pm)
        if statistic == "closure":
            Zo = C3(P_obs)[None, :, :]
            Zm = C3b(Pm)
            s2 = (1.0 / W[:-1, :-1] + 1.0 / W[:-1, 1:] +
                  1.0 / W[1:, :-1] + 1.0 / W[1:, 1:])
            return np.sum(2.0 * (1.0 - np.cos(wrap(Zo - Zm))) / s2, axis=(1, 2))
        Z = np.exp(1j * R)
        do_c = (statistic == "raw_nuis")
        a = np.zeros((Z.shape[0], NPSR))
        b = np.zeros((Z.shape[0], NSRC))
        for _ in range(40):                       # weighted two-way circular centring
            b = np.angle(np.sum(W[None] * Z * np.exp(-1j * a[:, :, None]), axis=1))
            if do_c:
                a = np.angle(np.sum(W[None] * Z * np.exp(-1j * b[:, None, :]), axis=2))
        E = wrap(R - a[:, :, None] - b[:, None, :])
        return np.sum(W[None] * 2.0 * (1.0 - np.cos(E)), axis=(1, 2))

    def fit_at_sky(srcs_try, P_obs, statistic, n_sweep=6, n_grid=3001):
        """Minimise the statistic over the 6 pulsar distances by coordinate
        descent with a DENSE, VECTORISED 1-D scan per pulsar over its full
        +-3 sigma prior (so every fringe in the prior is visited, not just the
        local one), then a golden refine in the winning cell."""
        A = coef(srcs_try)
        g0 = A * spine["L_m"][:, None]
        lo, hi = -3.0 * spine["Lerr_m"], 3.0 * spine["Lerr_m"]
        dL = np.zeros(NPSR)
        best = np.inf
        for _ in range(n_sweep):
            for p in range(NPSR):
                grid = np.linspace(lo[p], hi[p], n_grid)
                D = np.tile(dL, (n_grid, 1)); D[:, p] = grid
                Pm = g0[None] + A[None] * D[:, :, None]
                vals = stat_batch(Pm, P_obs, statistic)
                j = int(np.argmin(vals))
                a0 = grid[max(j - 1, 0)]
                b0 = grid[min(j + 1, n_grid - 1)]
                sub = np.linspace(a0, b0, 401)              # refine in the cell
                D = np.tile(dL, (401, 1)); D[:, p] = sub
                Pm = g0[None] + A[None] * D[:, :, None]
                v2 = stat_batch(Pm, P_obs, statistic)
                dL[p] = float(sub[int(np.argmin(v2))])
            v = float(stat_batch((g0 + A * dL[:, None])[None], P_obs,
                                 statistic)[0])
            if v > best - 1e-12:
                best = min(best, v)
                break
            best = v
        return best, dL

    rs = np.random.default_rng(918273)
    scram_skies = [(float(rs.uniform(-1, 1)), float(rs.uniform(0, 2 * np.pi)))
                   for _ in range(N_SCRAM)]

    results = {}
    for arm, c_inj in (("clean", None),
                       ("clock", rs.uniform(0, 2 * np.pi, NPSR))):
        # arm "clock": an UNMODELLED per-pulsar achromatic clock term is present
        # in the data.  This is the systematic the raw (c_p-blind) statistic
        # cannot represent -- and the one C3 cancels by construction.
        P_obs = phase_matrix(spine, srcs, c=c_inj)
        stats = ["raw_blind", "raw_nuis", "closure"]
        tru = {st: fit_at_sky(srcs, P_obs, st)[0] for st in stats}
        print(f"\n  ARM '{arm}' (per-pulsar clock injected: "
              f"{'YES' if c_inj is not None else 'no'})")
        print(f"    TRUE sky : " +
              "  ".join(f"{st}={tru[st]:.4e}" for st in stats))
        ens = {st: [] for st in stats}
        for j, (ct, ph) in enumerate(scram_skies):
            st_srcs = [dict(s) for s in srcs]
            st_srcs[0]["costheta"] = ct
            st_srcs[0]["phi"] = ph
            for st in stats:
                ens[st].append(fit_at_sky(st_srcs, P_obs, st)[0])
            print(f"     scram {j:2d}: " +
                  "  ".join(f"{st}={ens[st][-1]:.4e}" for st in stats))
        for st in stats:
            ens[st] = np.array(ens[st])
        results[arm] = dict(true=tru, ens={k: v.tolist() for k, v in ens.items()})
        print(f"    --- ARM '{arm}' SUMMARY ---")
        for st in stats:
            n_beat = int(np.sum(ens[st] <= tru[st] + 1e-9))
            print(f"      {st:10s}  truth={tru[st]:.4e}  scram min={ens[st].min():.4e}"
                  f"  median={np.median(ens[st]):.4e}"
                  f"  MANUFACTURE (scram beats truth): {n_beat}/{N_SCRAM}")

    np.savez(OUT,
             algebra=json.dumps(algebra),
             names=np.array(spine["names"]),
             fk=fk, rho=rho, nfringe_3sig=nfr,
             money=json.dumps(results),
             scram_skies=np.array(scram_skies),
             srcs=json.dumps(srcs), zero_noise=True)
    print(f"\n[T4] banked -> {OUT}")


if __name__ == "__main__":
    main()
