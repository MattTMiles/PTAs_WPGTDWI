#!/usr/bin/env python3
"""LEDGER A2 -- THE RED-NOISE MONITOR: is the loop leaking source power into per-pulsar RN?

=============================================================================
THE AUDIT (verdict: REAL -- the per-pulsar RN parameters are FROZEN, never refit)
=============================================================================
Traced end to end:

  glacier_pop.build_glacier_problem (docstring, glacier_pop.py:530-535):
      "RN frozen at the seed-0 reference draw exactly as make_geometry_injection does
       (fixed noise covariance across skies)"
  trackB_b1_core.build_b1_amortised:190-192  common_gp = makecommongp_fourier(..., "rednoise")
  trackB_b1_core.build_b1_amortised:252      theta_keys = dist_keys + cw_param_keys(num_cw)
  trackB_b1_core.build_b1_amortised:256      frozen = {k: temp_dict[k] for k in truth_full
                                                       if k not in theta_keys}
      -> every `*_rednoise_log10_A` and `*_rednoise_gamma` lands in `frozen`. They are
         NOT in theta, so no stage of the loop -- Fisher (a), frontier (b), drain (c),
         E-step (d), M-step (e), scoreboard (f) -- can ever move them.
  trackB_b1_core.B1Split.__init__:357        Pinv, ldP = vsm.P_var.make_inv()(self.frozen)
  trackB_b1_core.B1Split.__init__:358        self.Phi_rn_diag = 1.0 / Pinv
      -> the RN prior variance, and the per-pulsar Cholesky `self.cf` built from it, are
         cut ONCE at construction and reused by every iteration of every cell.

AND THE DATA IS DRAWN FROM THE SAME OBJECT:
  trackB_b1_core.NoiseDrawer:590             self.sig_rn = sqrt(sp.Phi_rn_diag)
  trackB_b1_core.NoiseDrawer.draw:609        r += Fs_rn[p] @ (sig_rn[p] * N(0,1))

So the model's RN and the data's RN share hyperparameters EXACTLY. There is no
mis-specification, and there is no refit. That is a strong idealisation and it should be
named: **the campaign never pays for not knowing the noise.** ANCHOR SS4 already measured
that deliberate RN mis-specification does not move the certification floor (the tail is
template-dominated), so this is not a floor defect. The exposure it leaves open is
different and is what this monitor watches:

  THE LEAK. Feeding a source removes power from the residual; the M-step then WANDERS the
  fed member off truth by 0.01-0.44 dex (S4.15, measured). Mismatched source power has to
  go somewhere. With the GWB amplitude fitted (BackgroundFit) and per-pulsar RN frozen,
  the drain absorbs the common-mode part -- but the PER-PULSAR part has nowhere to go
  except the residual, where it inflates the effective per-pulsar noise the fringe
  likelihood is being scored against. Frozen RN cannot show that as a fitted parameter
  moving; it shows up only as posterior RN power drifting away from its prior. Which is
  exactly what this monitor measures.

REFIT-IN-LOOP IS A DESIGN CHANGE AND IS **SPEC ONLY** (per the charter and the brief):
see LEDGER_stats_audit.md SS A2-SPEC. It would move every banked floor.

=============================================================================
THE MONITOR (additive; it cannot change any banked verdict)
=============================================================================
The per-pulsar RN GP coefficients have an EXACT Gaussian conditional posterior at the
current template state -- and every object it needs is already held by B1Split, so the
monitor costs one 1-D noise solve + one matvec + one triangular solve per pulsar per
iteration. No new factorisation, no extra likelihood call, no RNG draw (so it cannot
perturb a realisation, and a resumed checkpoint is unaffected).

  y_p       = data_p - delay_p(theta_rec)              (B1Split.ys[p], the split's own)
  b_p       = F_p^T N_p^-1 y_p                          (B1Split.Ft[p], solve1d[p])
  Lambda_p  = F_p^T N_p^-1 F_p + Phi_p^-1               (B1Split.cf -- ALREADY FACTORED)
  E[c_p]    = Lambda_p^-1 b_p ,  Cov[c_p] = Lambda_p^-1

  S_p  =  (1/n_k) SUM_k  ( E[c_pk]^2 + Cov[c_p]_kk ) / Phi_pk         <- THE STATISTIC

S_p is the posterior RN power in units of its own prior: the whitened per-pulsar red-noise
occupancy. Under a correct template it shrinks toward the prior; when source mismatch
leaks into pulsar p's red band it RISES. Its exact sampling spread under the posterior is
analytic for a Gaussian (no Monte Carlo):

  Var(S_p) = (1/n_k^2) SUM_k ( 2 Cov_kk^2 + 4 E[c_pk]^2 Cov_kk ) / Phi_pk^2

THE BAND, AND THE FLAG (pre-registered here, forward):
  reference   S_p(0)         -- the iteration-0 value in the SAME cell (self-referenced,
                               so a venue's own RN occupancy is never compared across
                               venues, which would be meaningless)
  band        sqrt(Var(S_p)(0) + Var(S_p)(i))
  FLAG        |S_p(i) - S_p(0)| > 3 * band                (per pulsar, per iteration)
  CELL FLAG   any pulsar flagged, OR the array-median |drift| exceeding 3 * its own band
              (the common-mode leg -- a leak that touches every pulsar a little is
              exactly what a per-pulsar 3-sigma test would miss)
The flag is REPORT-ONLY at introduction. It does not STOP a cell and it does not enter
any certification. Promotion to a STOP condition is a pre-registered bar change and is
therefore Matt's call (charter #1), not this agent's.

DECLARED LIMITATION. S_p is CONDITIONAL on the frozen GWB block and the current template;
it is a drift detector, not an inference of log10_A_rn. The profiled per-pulsar RN
amplitude (the same band-selective diagonal-rescale identity BackgroundFit uses for the
GWB, applied to Phi_rn instead) is the natural upgrade and is specced in the report; it
costs one extra ncomp x ncomp Cholesky per pulsar per grid point and is affordable, but
it is a NEW fitted quantity and so is not introduced silently mid-campaign.

Usage:
    python3 hpc_harbor/ledger/ledger_a2_rn_monitor.py gate      # CPU, synthetic, no venue
    # in a driver:  from ledger_a2_rn_monitor import rn_occupancy, rn_drift_flags
"""
import argparse
import sys

import numpy as np


# ============================================================
# the statistic
# ============================================================
def rn_occupancy(sp, theta, data_tuple, pmask, jnp=None, dsm=None):
    """Per-pulsar whitened RN occupancy S_p and its analytic variance, at this state.

    `sp` is a trackB_b1_core.B1Split. Reads only objects it already holds:
    sp.ys, sp.solve1d, sp.Ft, sp.cf, sp.Phi_rn_diag, sp.smask. Consumes no RNG.

    Returns (S (npsr,), VarS (npsr,), Emc2 (npsr, ncomp)) -- the last is the per-mode
    posterior mean-square in prior units, banked so a leak can be localised in FREQUENCY
    (a source leak is band-limited; a genuine RN change is not).
    """
    if jnp is None:
        import jax.numpy as jnp                                    # noqa: PLC0415
    if dsm is None:
        from discovery import matrix as dsm                        # noqa: PLC0415
    params = sp._params(jnp.asarray(theta), data_tuple, jnp.asarray(pmask), sp.smask)
    npsr = sp.npsr
    Phi = np.asarray(sp.Phi_rn_diag, float)                        # (npsr, ncomp)
    nk = Phi.shape[1]
    b = []
    for p in range(npsr):
        yp = sp.ys[p](params)
        Nmy, _ = sp.solve1d[p](yp)
        b.append(np.asarray(sp.Ft[p] @ Nmy, float))
    b = np.stack(b)                                                # (npsr, ncomp)
    mean = np.asarray(dsm.matrix_solve(sp.cf, jnp.asarray(b)), float)   # Lambda^-1 b
    # Cov diagonal: solve Lambda X = I once per pulsar (ncomp is 30-95; this is cheap)
    eye = jnp.broadcast_to(jnp.eye(nk), (npsr, nk, nk))
    cov = np.asarray(dsm.matrix_solve(sp.cf, eye), float)
    cd = np.einsum("pkk->pk", cov) if cov.ndim == 3 else np.stack(
        [np.diag(cov[p]) for p in range(npsr)])
    Emc2 = (mean ** 2 + cd) / Phi
    S = Emc2.mean(axis=1)
    VarS = ((2.0 * cd ** 2 + 4.0 * mean ** 2 * cd) / Phi ** 2).sum(axis=1) / nk ** 2
    return S, VarS, Emc2


def rn_drift_flags(S, VarS, S0, VarS0, kappa=3.0):
    """The pre-registered band + flag. Returns (drift, band, flag_psr, flag_cell)."""
    S = np.asarray(S, float); S0 = np.asarray(S0, float)
    drift = S - S0
    band = np.sqrt(np.asarray(VarS, float) + np.asarray(VarS0, float))
    flag_psr = np.abs(drift) > kappa * band
    med = float(np.median(drift))
    med_band = float(np.median(band)) / np.sqrt(max(len(drift), 1))
    flag_cell = bool(flag_psr.any() or abs(med) > kappa * med_band)
    return drift, band, flag_psr, flag_cell


# ============================================================
# gate (CPU, synthetic -- no venue, no GPU, no pulsar load)
# ============================================================
class _FakeSplit:
    """Minimal stand-in exercising the exact algebra path rn_occupancy uses."""

    def __init__(self, npsr=6, nk=10, ntoa=40, seed=20260729):
        rng = np.random.default_rng(seed)
        self.npsr, self.nk, self.ntoa = npsr, nk, ntoa
        self.Phi_rn_diag = rng.uniform(0.5, 2.0, size=(npsr, nk))
        self.F = rng.normal(size=(npsr, ntoa, nk))
        self.Ninv = rng.uniform(0.5, 1.5, size=(npsr, ntoa))
        self.smask = None
        self.cf = None            # the real B1Split carries the factored Lambda here
        self._y = np.zeros((npsr, ntoa))
        self.Ft = [self.F[p].T for p in range(npsr)]
        self.solve1d = [(lambda y, _p=p: (self.Ninv[_p] * y, 0.0)) for p in range(npsr)]
        self.ys = [(lambda params, _p=p: self._y[_p]) for p in range(npsr)]
        self.Lam = np.stack([self.F[p].T @ (self.Ninv[p][:, None] * self.F[p])
                             + np.diag(1.0 / self.Phi_rn_diag[p]) for p in range(npsr)])

    def _params(self, *a, **k):
        return {}


class _FakeDsm:
    def __init__(self, split):
        self.s = split

    def matrix_solve(self, cf, rhs):
        r = np.asarray(rhs)
        if r.ndim == 2:
            return np.stack([np.linalg.solve(self.s.Lam[p], r[p])
                             for p in range(self.s.npsr)])
        return np.stack([np.linalg.solve(self.s.Lam[p], r[p]) for p in range(self.s.npsr)])


class _FakeJnp:
    asarray = staticmethod(np.asarray)
    eye = staticmethod(np.eye)

    @staticmethod
    def broadcast_to(a, shp):
        return np.broadcast_to(a, shp)


def gate():
    print("=== LEDGER A2 -- RN MONITOR GATES (CPU, synthetic) ===", flush=True)
    ok = True
    sp = _FakeSplit()
    dsm = _FakeDsm(sp)
    jnp = _FakeJnp()

    # G-R1: zero residual -> S is the pure shrinkage floor, and it is FINITE and < 1
    sp._y[:] = 0.0
    S0, V0, _ = rn_occupancy(sp, np.zeros(3), (), np.ones(sp.npsr), jnp=jnp, dsm=dsm)
    b1 = bool(np.all(np.isfinite(S0)) and np.all(S0 > 0) and np.all(S0 < 1.0))
    print(f"  G-R1 zero-residual occupancy finite and below prior (shrinkage): "
          f"S in [{S0.min():.4f}, {S0.max():.4f}] -> {'PASS' if b1 else 'FAIL'}")
    ok &= b1

    # G-R2: injecting RN-band power into ONE pulsar raises ONLY that pulsar's S, and the
    # flag fires on it and on no other. This is the leak the monitor exists to catch.
    rng = np.random.default_rng(7)
    c_inj = np.sqrt(sp.Phi_rn_diag[2]) * rng.normal(size=sp.nk) * 4.0
    sp._y[2] = sp.F[2] @ c_inj
    S1, V1, _ = rn_occupancy(sp, np.zeros(3), (), np.ones(sp.npsr), jnp=jnp, dsm=dsm)
    drift, band, fp, fc = rn_drift_flags(S1, V1, S0, V0)
    b2 = bool(fp[2] and not fp[np.arange(sp.npsr) != 2].any() and fc)
    print(f"  G-R2 single-pulsar RN-band leak: drift {drift[2]:+.4f} vs band "
          f"{band[2]:.4f} ({drift[2]/band[2]:+.1f} sigma); flagged pulsars "
          f"{list(np.where(fp)[0])} (need [2]) -> {'PASS' if b2 else 'FAIL'}")
    ok &= b2

    # G-R3: NO leak -> no flag (the false-positive leg; the monitor must be quiet when
    # nothing has happened, else it is noise in the fan's log)
    sp._y[:] = 0.0
    S2, V2, _ = rn_occupancy(sp, np.zeros(3), (), np.ones(sp.npsr), jnp=jnp, dsm=dsm)
    _, _, fp2, fc2 = rn_drift_flags(S2, V2, S0, V0)
    b3 = (not fp2.any()) and (not fc2)
    print(f"  G-R3 unchanged state raises no flag: flagged {int(fp2.sum())}, cell flag "
          f"{fc2} -> {'PASS' if b3 else 'FAIL'}")
    ok &= b3

    # G-R4: COMMON-MODE leak (small, every pulsar) -- the leg a per-pulsar 3-sigma test
    # misses by construction and the median leg is there to catch
    for p in range(sp.npsr):
        sp._y[p] = sp.F[p] @ (np.sqrt(sp.Phi_rn_diag[p])
                              * np.random.default_rng(100 + p).normal(size=sp.nk) * 0.9)
    S3, V3, _ = rn_occupancy(sp, np.zeros(3), (), np.ones(sp.npsr), jnp=jnp, dsm=dsm)
    d3, bnd3, fp3, fc3 = rn_drift_flags(S3, V3, S0, V0)
    print(f"  G-R4 common-mode leak: per-pulsar flags {int(fp3.sum())}/{sp.npsr}, "
          f"median drift {np.median(d3):+.4f}, CELL flag {fc3} "
          f"-> {'PASS' if fc3 else 'FAIL'} (the median leg is what fires here)")
    ok &= bool(fc3)

    print(f"\n=== LEDGER A2 GATES: {'PASS' if ok else 'FAIL'} ===", flush=True)
    return 0 if ok else 1


# ============================================================
# venue replay -- the monitor RUN on banked fan cells, post hoc
# ============================================================
def replay(rung, arm, sky, t_label, wscale, n_src, venue_T=None):
    """Compute the monitor on a banked cell's OWN checkpoints, without touching the
    driver. theta_rec and fed_mask are banked per iteration (glacier_loop.py:566, 557), so
    the whole iteration history can be re-scored after the fact. Reads banks; writes only
    a new ledger_a2_* bank. No GPU is required -- rn_occupancy never calls
    estep_per_target, so the XLA-CPU vm.max_map_count hazard (CPU-lane map-count note)
    does not apply and this runs on the free `batch` partition."""
    import os as _os, sys as _sys, glob as _glob
    _os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    HSYMT = "/home/mattm/projects/HSYMT"
    for q in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
              "hpc_harbor/atlas", "hpc_harbor/spark", "hpc_harbor/glacier"):
        _sys.path.insert(0, f"{HSYMT}/{q}")
    import glacier_loop as GL
    import glacier_pop as GP
    from glacier_pop import bank_npz, lane_tag
    import jax.numpy as jnp
    from discovery import matrix as dsm

    T = venue_T or t_label
    wtag = "" if wscale == 1.0 else "_w" + f"{wscale:g}".replace(".", "p")
    stem = (f"gl1_cell_{arm}_s{sky}_T{T}_{GL.TIER}" if rung is None else
            f"gl2_{rung}{wtag}_cell_{arm}_s{sky}_T{T}_{GL.TIER}")
    cks = sorted(_glob.glob(f"{GP.OUT}/{stem}_i*__*.npz"))
    if not cks:
        raise SystemExit(f"no banked checkpoints for {stem}")
    print(f"[LEDGER-A2] replay {stem}: {len(cks)} iterations, lane {lane_tag()}", flush=True)

    jax, jnpm, C, B1, TE, IG, F, FL = GL._stack()
    if rung is None:
        pop = GP.draw_population(GP.SEED_POP_BASE + sky, n_src=n_src,
                                 band_log10f=GL.BAND_CAMPAIGN)
    else:
        pop, _ = GP.draw_population_conditional(sky, rung, n_src=n_src,
                                                band_log10f=GL.BAND_CAMPAIGN)
    slots, *_ = GL.embed_igniter(pop, GL.E_IGNITER[arm], GL.VENUE_SPAN_S[T])
    pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
    G = GP.build_glacier_problem(T, pop_s, verbose=True)
    G["slots"] = slots
    amo = G["amo"]; nd = amo["n_dist"]
    ones = jnpm.ones(amo["npsr"])
    sb = GL.CertScoreboard(G, amo, jnpm, C, prior_key=GL.TIER)
    noise_seed = GL.SEED_NOISE_BASE + 1000 * sky
    L_true, _ = sb.draw_truth(IG, noise_seed + 10_000_000, tier=GL.TIER)
    th_t = np.asarray(amo["theta_truth"], float).copy(); th_t[sb.AI] = L_true
    clean = amo["inject_delay"](jnpm.asarray(th_t), ones)
    ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=GL.LGWB_BANKS[T], strict=True)
    noise = ndraw.draw(noise_seed, components=("white", "rn"), white_scale=wscale)
    data = tuple(jnpm.asarray(np.asarray(a) + np.asarray(b)) for a, b in zip(clean, noise))

    S_all, V_all, E_all, its = [], [], [], []
    for ck in cks:
        z = np.load(ck, allow_pickle=True)
        it = int(z["iter"])
        carried = np.where(~np.asarray(z["fed_mask"], bool))[0]
        th = GL.theta_with_absent(np.asarray(z["theta_rec"], float), nd, carried)
        S, V, E = rn_occupancy(sb.sp, th, data, ones, jnp=jnpm, dsm=dsm)
        S_all.append(S); V_all.append(V); E_all.append(E); its.append(it)
        print(f"  iter {it}: median S {np.median(S):.5f}  max S {S.max():.5f}", flush=True)
    S_all = np.stack(S_all); V_all = np.stack(V_all)
    o = np.argsort(its); S_all, V_all, its = S_all[o], V_all[o], np.array(its)[o]
    drifts, flags, cellflags = [], [], []
    for i in range(len(its)):
        d, b, fp, fc = rn_drift_flags(S_all[i], V_all[i], S_all[0], V_all[0])
        drifts.append(d); flags.append(fp); cellflags.append(fc)
        print(f"  iter {its[i]}: median drift {np.median(d):+.5f}, flagged "
              f"{int(fp.sum())}/{len(fp)} pulsars, CELL FLAG {fc}", flush=True)
    # NB: bank_npz's own first parameter is named `stem`, so the cell's stem must NOT be
    # passed as `stem=` -- it collides. Banked as `cell_stem`.
    bank_npz(f"ledger_a2_rn_{stem}", cell_stem=stem, iters=its, S=S_all, VarS=V_all,
             Emc2_last=E_all[-1], drift=np.stack(drifts), flag_psr=np.stack(flags),
             flag_cell=np.array(cellflags), kappa=3.0,
             note=("LEDGER A2 RN monitor, REPLAYED from banked checkpoints. Additive and "
                   "report-only: no banked column is rewritten, no verdict moves."))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "replay"], nargs="?", default="gate")
    ap.add_argument("--rung", default=None)
    ap.add_argument("--arm", default="none")
    ap.add_argument("--sky", type=int, default=1)
    ap.add_argument("--t", type=int, default=30)
    ap.add_argument("--wscale", type=float, default=1.0)
    ap.add_argument("--n", type=int, default=256)
    a = ap.parse_args()
    if a.mode == "replay":
        return replay(a.rung, a.arm, a.sky, a.t, a.wscale, a.n)
    return gate()


if __name__ == "__main__":
    sys.exit(main())
