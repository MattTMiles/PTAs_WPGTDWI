"""SIEVE-C T1 -- WAVEFRONT CURVATURE.  Zero-noise / Asimov only (cronus is
RECOMPUTED-UNSAFE for noisy realisations; nothing here draws noise).

Inject a CURVED (spherical) wavefront at D = {16, 50, 200} Mpc at loud strain;
recover with FLAT (plane-wave) templates.  Readouts:
  (a) recovered source sky / distance bias vs D
  (b) dlnL cost of flatness vs D
  (c) at 16 Mpc: profile lnL over D with curvature IN the template
      -- the GW-parallax check.

Each readout is run in THREE pulsar-distance regimes, because the answer
depends entirely on which one the pipeline is in:
    pinned  : L_p held at truth (no distance freedom at all)
    prior   : L_p free under the canonical Gaussian prior (best_distances.txt)
    free    : L_p free, unconstrained (local fit -- no fringe hopping)
"""

import os
import sys
import json
import numpy as np
from scipy.optimize import minimize, minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sieve_toy as S

OUT = os.path.join(S.ROOT, "reports", "sieve_t1_curvature.npz")

# ---- the injected source (loud; band-centre; a real-sky-class geometry) ----
SRC0 = dict(costheta=0.35, phi=1.9, log10_h=-13.2, log10_f=-8.0,
            cosi=0.4, psi=0.7, phi0=1.1)
D_LIST_MPC = [16.0, 50.0, 200.0]
SKEYS = ["costheta", "phi", "log10_h", "log10_f", "cosi", "psi", "phi0"]


def _unpack(x, spine, src0, nd):
    s = dict(src0)
    for i, k in enumerate(SKEYS):
        s[k] = x[i]
    s["costheta"] = np.clip(s["costheta"], -1.0, 1.0)
    s["cosi"] = np.clip(s["cosi"], -1.0, 1.0)
    if nd:
        # distances carried as offsets in units of the TRUE-sky fringe spacing
        s["L_m"] = spine["L_m"] + x[len(SKEYS):] * spine["_dLfid"]
    return s


def _psr_nll(spine, data, s, front):
    """Per-pulsar negative lnL contributions -- the pulsar-term residual of
    pulsar p depends ONLY on L_p, so the distance block is EXACTLY separable
    at fixed source parameters.  That separability is what makes the fits below
    reliable (a 23-D simplex is not)."""
    r = S.residuals(spine, s, front=front)
    d = (data - r) / spine["sigma"][:, None]
    return 0.5 * np.sum(d * d, axis=1)


def fit_flat(spine, data, src0, regime, front="plane", D_m=None,
             n_sweep=12, x0=None):
    """Maximise lnL (+ Gaussian distance prior) by ALTERNATING:
        (i)  7-D source-parameter simplex at fixed distances
        (ii) EXACT per-pulsar 1-D distance minimisation (separable)
    `regime`: pinned | prior | free.
    """
    nd = spine["npsr"] if regime in ("prior", "free") else 0
    sig_frac = spine["Lerr_m"] / spine["_dLfid"]
    xs = np.array([src0[k] for k in SKEYS], float) if x0 is None else np.asarray(x0[:7], float)
    u = np.zeros(spine["npsr"]) if x0 is None else np.asarray(x0[7:], float).copy()

    def build(xs_, u_):
        s = dict(src0)
        for i, k in enumerate(SKEYS):
            s[k] = xs_[i]
        s["costheta"] = float(np.clip(s["costheta"], -1.0, 1.0))
        s["cosi"] = float(np.clip(s["cosi"], -1.0, 1.0))
        if D_m is not None:
            s["D_m"] = D_m
        if nd:
            s["L_m"] = spine["L_m"] + u_ * spine["_dLfid"]
        return s

    def total(xs_, u_):
        s = build(xs_, u_)
        v = -S.lnL(spine, data, s, front=front)
        if regime == "prior":
            v += 0.5 * float(np.sum((u_ / sig_frac) ** 2))
        return v

    prev = np.inf
    for _ in range(n_sweep):
        r = minimize(lambda z: total(z, u), xs, method="Nelder-Mead",
                     options=dict(maxiter=3000, maxfev=6000,
                                  xatol=1e-11, fatol=1e-11))
        xs = r.x
        if nd:
            # (ii) exact separable 1-D sweep, one pulsar at a time, on a local
            # bracket of +-1 fringe (no fringe hopping -- this is a LOCAL fit).
            for p in range(spine["npsr"]):
                grid = u[p] + np.linspace(-1.0, 1.0, 81)
                best, bu = np.inf, u[p]
                for g in grid:
                    uu = u.copy(); uu[p] = g
                    v = _psr_nll(spine, data, build(xs, uu), front)[p]
                    if regime == "prior":
                        v += 0.5 * (g / sig_frac[p]) ** 2
                    if v < best:
                        best, bu = v, g
                lo, hi = bu - 0.05, bu + 0.05

                def f1(g):
                    uu = u.copy(); uu[p] = g
                    v = _psr_nll(spine, data, build(xs, uu), front)[p]
                    if regime == "prior":
                        v += 0.5 * (g / sig_frac[p]) ** 2
                    return v
                rr = minimize_scalar(f1, bracket=None, bounds=(lo, hi),
                                     method="bounded",
                                     options=dict(xatol=1e-13))
                u[p] = rr.x
        cur = total(xs, u)
        if abs(prev - cur) < 1e-13 * max(1.0, abs(cur)):
            break
        prev = cur
    return total(xs, u), build(xs, u), nd, np.concatenate([xs, u])


def main():
    spine = S.build_spine(npsr=16, T_yr=18.0, cadence_per_yr=20.0)
    spine["_dLfid"] = S.fringe_spacing(spine, SRC0)
    lam_pc = S.C_LIGHT / 10.0 ** SRC0["log10_f"] / S.PC
    print(f"[T1] GW wavelength lambda = {lam_pc:.4f} pc; fringe spacing "
          f"{(spine['_dLfid']/S.PC).min():.3f}-{(spine['_dLfid']/S.PC).max():.3f} pc")

    rows, prof = [], {}

    # ---------------------------------------------------------------
    # G0 -- the identity gate: is the curved front EXACTLY a per-pulsar
    #       effective distance offset dL_eff = L^2 (1+cos mu) / (2 D) ?
    # ---------------------------------------------------------------
    print("\n[T1-G0] curvature <-> distance identity gate")
    gate = []
    for Dm in D_LIST_MPC:
        D = Dm * S.MPC
        src_c = dict(SRC0, D_m=D)
        d_sph = S.residuals(spine, src_c, front="sphere")
        d_frs = S.residuals(spine, src_c, front="fresnel")
        dLe = S.curvature_dL_eff(spine, SRC0, D)
        src_shift = dict(SRC0, L_m=spine["L_m"] + dLe)
        d_shift = S.residuals(spine, src_shift, front="plane")
        s2 = S.snr2(spine, d_sph)
        c_sph_shift = -S.lnL(spine, d_sph, src_shift, front="plane")
        c_sph_frs = -S.lnL(spine, d_sph, src_c, front="fresnel")
        cyc = np.abs(dLe / spine["_dLfid"])
        print(f"  D={Dm:6.1f} Mpc  SNR={np.sqrt(s2):8.2f}  "
              f"max|dL_eff|/fringe={cyc.max():.4f}  "
              f"dlnL(sphere vs dL_eff-shifted flat)={c_sph_shift:.3e}  "
              f"dlnL(sphere vs fresnel)={c_sph_frs:.3e}")
        gate.append(dict(D_Mpc=Dm, snr=float(np.sqrt(s2)),
                         max_cycles=float(cyc.max()),
                         dlnL_sphere_vs_shift=float(c_sph_shift),
                         dlnL_sphere_vs_fresnel=float(c_sph_frs),
                         dL_eff_pc=(dLe / S.PC).tolist(),
                         cycles=cyc.tolist()))
        _ = d_frs, d_shift

    # ---------------------------------------------------------------
    # (a) + (b) -- recover CURVED data with FLAT templates
    # ---------------------------------------------------------------
    print("\n[T1-a/b] flat-template recovery of curved data")
    for Dm in D_LIST_MPC:
        D = Dm * S.MPC
        src_c = dict(SRC0, D_m=D)
        data = S.residuals(spine, src_c, front="sphere")
        snr = np.sqrt(S.snr2(spine, data))
        warm = None
        for regime in ("pinned", "prior", "free"):
            tot, sbest, nd, xfull = fit_flat(spine, data, SRC0, regime, x0=warm)
            if regime == "prior":
                warm = xfull                          # warm-start `free` from `prior`
            # always quote the PURE mismatch cost, prior term stripped
            cost = -S.lnL(spine, data, sbest, front="plane")
            shat_t = S.sky_unit(SRC0["costheta"], SRC0["phi"])
            shat_r = S.sky_unit(sbest["costheta"], sbest["phi"])
            dsky = np.degrees(np.arccos(np.clip(shat_t @ shat_r, -1, 1)))
            dh = sbest["log10_h"] - SRC0["log10_h"]
            df = sbest["log10_f"] - SRC0["log10_f"]
            dL = ((np.asarray(sbest.get("L_m", spine["L_m"])) - spine["L_m"]) / S.PC
                  if nd else np.zeros(spine["npsr"]))
            row = dict(D_Mpc=Dm, regime=regime, snr=float(snr),
                       dlnL_cost=float(cost), dsky_deg=float(dsky),
                       dlog10h=float(dh), dlog10f=float(df),
                       max_dL_pc=float(np.abs(dL).max()),
                       rms_dL_pc=float(np.sqrt(np.mean(dL ** 2))))
            rows.append(row)
            print(f"  D={Dm:6.1f} {regime:7s} dlnL_cost={cost:12.4e} "
                  f"dsky={dsky:8.4f} deg  dlog10h={dh:+.2e} "
                  f"dlog10f={df:+.2e}  max|dL|={np.abs(dL).max():.4f} pc")

    # ---------------------------------------------------------------
    # (c) -- GW parallax: profile lnL over D with curvature IN the template
    # ---------------------------------------------------------------
    print("\n[T1-c] GW-parallax profile at D_true = 16 Mpc "
          "(curvature IN the template)")
    D_true = 16.0 * S.MPC
    data = S.residuals(spine, dict(SRC0, D_m=D_true), front="sphere")
    D_grid = np.array([4., 8., 12., 16., 24., 40., 80., 200., 1000., 1e5]) * S.MPC
    for regime in ("pinned", "prior", "free"):
        vals, warm = [], None
        for Dg in D_grid:                     # warm-chained along the grid
            tot, sbest, nd, xfull = fit_flat(spine, data, SRC0, regime,
                                             front="sphere", D_m=Dg, x0=warm)
            warm = xfull
            vals.append(-tot)                 # PROFILE of the (penalised) objective
        vals = np.array(vals)
        prof[regime] = vals
        best = D_grid[int(np.argmax(vals))] / S.MPC
        print(f"  {regime:7s} peak at D={best:9.1f} Mpc; "
              f"lnL(16 Mpc)-lnL(1e5 Mpc) = {vals[3] - vals[-1]:+.4e} nat")
        print("           " + "  ".join(f"{v:+.2e}" for v in vals - vals.max()))

    np.savez(OUT,
             gate=json.dumps(gate), rows=json.dumps(rows),
             D_grid_Mpc=D_grid / S.MPC,
             prof_pinned=prof["pinned"], prof_prior=prof["prior"],
             prof_free=prof["free"],
             names=np.array(spine["names"]), L_kpc=spine["L_kpc"],
             Lerr_kpc=spine["Lerr_kpc"], sigma_s=spine["sigma"],
             fringe_pc=spine["_dLfid"] / S.PC,
             src0=json.dumps(SRC0), zero_noise=True)
    print(f"\n[T1] banked -> {OUT}")


if __name__ == "__main__":
    main()
