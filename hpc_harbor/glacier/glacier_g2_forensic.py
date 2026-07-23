"""GLACIER g2a-ii FORENSIC -- why does the band-matched profile peg at the grid floor?

Jobs 12645909 and 12646058 both read log10 Ahat = -15.6 [EDGE HIT] against A_eff = -14.54,
and the banked profile is PURE OCCAM: lnl(-15.6) - lnl(-14.6) = +3.8 nats, monotone, with
no visible quadratic (signal) gain anywhere on the 2-dex grid. Two hypotheses produce
exactly that curve, and they demand opposite remedies:

  H-BUG   : the __smask key is dropped between BackgroundFit._data_terms and MaskedDelay,
            the template silently contains every source at truth, and the zero-noise
            residual is identically zero (quad == 0 exactly). Remedy: fix the plumbing.
  H-FAINT : the residual is real (the full drawn population) but its P-weighted in-band
            power sits below the array's in-band noise level at this (T, n), so the ML
            amplitude is genuinely at/below the grid floor -- the gate as designed asks
            for more information than the venue holds. Remedy: redesign the gate venue,
            not the machinery.

This module MEASURES which. Zero-noise Asimov throughout -- no noise realisation is drawn,
no L_gwb basis is touched, so the run is gate-class, not campaign-class.

  F1  template/residual identities: template(smask=0) must be exactly 0 per pulsar;
      the residual under smask=0 must equal the injected data (per-pulsar RMS printed).
  F2  FtNmy under smask = zeros / ones / ABSENT key: ones vs absent must agree bit-exactly
      (the 1.0*x == x identity); ones must be ~0 (zero-noise: data == template at truth);
      zeros must carry the full signal projection. ||FtNmy(zeros)|| ~ 0 -> H-BUG.
  F3  profile decomposition on a grid extended 2 dex below the smoke grid: quad(A) and
      occam(A) printed SEPARATELY. H-BUG reads quad == 0 everywhere; H-FAINT reads
      quad > 0 with d(quad)/dlnA < d(occam)/dlnA at the floor.
  F4  POSITIVE CONTROL: a synthetic in-band GP realisation drawn at the NG15 reference
      amplitude from the EXACT banded prior (correlated draw from inv(Pinv0) restricted
      to in-band rows, fixed seed 4242), pushed through the identical fit path. The
      recovered (ahat, sigma) measures the venue's information content directly:
      recovers ~A0 with small sigma -> machinery sound, and sigma is the precision any
      amplitude gate can honestly demand at this (T, n).

Run:  glacier_g2_forensic.py --t 15 --n 32     (the smoke venue, ~15 min)
      glacier_g2_forensic.py --t 30 --n 256    (the main-gate venue)
Banks g2forensic_T{t}_n{n}__{lane}.npz via glacier_pop.bank_npz (lane-tagged, refusal on
cross-lane collision). No gate verdict is cut here -- this is measurement, and the g2a-ii
remedy is designed from its numbers, pre-registered, then re-gated.
"""
import argparse, os, sys, time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import glacier_pop as GP


def data_terms_nokey(bf, theta_arr, data_tuple, pmask):
    """BackgroundFit._data_terms MINUS the __smask key: the ABSENT-key limb of F2.
    Deliberately reproduces the pre-FORGE-G identity (every source in the template)."""
    jnp = bf.jnp
    params = dict(bf.frozen)
    for k, v in zip(bf.theta_keys, theta_arr):
        params[k] = v
    for dk, d in zip(bf.data_keys, data_tuple):
        params[dk] = d
    params["__pmask"] = pmask
    terms = bf.kterms(params)
    return float(jnp.sum(terms[0])), np.asarray(terms[1]).reshape(-1)


def quad_occam(bf, a_log10, FtNmy):
    """One profile point, decomposed: (quad, occam) with occam = ldP + logdet(Pinv+C).
    lnL = p0 + 0.5*(quad - occam); only differences across the grid are meaningful."""
    jnp, dsm = bf.jnp, bf.dsm
    Pinv, ldP_sum = bf._pinv_band(a_log10)
    if not hasattr(bf, "_C"):
        bf._C = np.asarray(dsm.jsp.linalg.block_diag(*bf.c_truth))
    cf = dsm.matrix_factor(jnp.asarray(Pinv + bf._C))
    ld_cf = dsm.matrix_norm * float(jnp.sum(jnp.log(jnp.diag(cf[0]))))
    quad = float(FtNmy @ np.asarray(dsm.matrix_solve(cf, jnp.asarray(FtNmy))))
    return quad, ldP_sum + ld_cf


def profile_from_terms(bf, FtNmy, grid):
    """The profile shape from a fixed FtNmy (p0 is grid-constant and dropped): returns
    lnl-up-to-constant, quads, occams, ahat, sigma, edge flag."""
    quads = np.empty(len(grid)); occs = np.empty(len(grid))
    for i, a in enumerate(grid):
        quads[i], occs[i] = quad_occam(bf, a, FtNmy)
    lnl = 0.5 * (quads - occs)
    k = int(np.argmax(lnl))
    if 0 < k < len(grid) - 1:
        a, b, c = lnl[k - 1], lnl[k], lnl[k + 1]
        d = grid[1] - grid[0]
        off = 0.5 * (a - c) / (a - 2 * b + c)
        ahat = grid[k] + off * d
        curv = (a - 2 * b + c) / d ** 2
        sig = 1.0 / np.sqrt(max(-curv, 1e-300))
        edge = False
    else:
        ahat, sig, edge = grid[k], np.inf, True
    return lnl, quads, occs, float(ahat), float(sig), bool(edge)


def control_draw(bf, seed=4242):
    """A correlated in-band GP realisation at the reference amplitude a0: invert the full
    banded prior Pinv0 once, restrict to in-band rows, Cholesky, draw. Out-of-band
    coefficients are ZERO (the control lives exactly where the fit scales), declared."""
    Phi0 = np.linalg.inv(bf.Pinv0)
    sub = Phi0[np.ix_(bf.mask_full, bf.mask_full)]
    sub = 0.5 * (sub + sub.T)
    jit = 1e-12 * float(np.mean(np.diag(sub)))
    for _ in range(3):
        try:
            L = np.linalg.cholesky(sub + jit * np.eye(sub.shape[0]))
            break
        except np.linalg.LinAlgError:
            jit *= 100.0
    else:
        raise GP.CampaignStop("control prior submatrix is not PSD at jitter 1e-8*diag -- "
                              "the banded prior itself is suspect; STOP.")
    rng = np.random.default_rng(seed)
    c = np.zeros(bf.npsr * bf.ngp)
    c[bf.mask_full] = L @ rng.standard_normal(int(bf.mask_full.sum()))
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t", type=int, default=15)
    ap.add_argument("--n", type=int, default=32)
    ap.add_argument("--seed", type=int, default=GP.SEED_POP_BASE)
    ap.add_argument("--flo", type=float, default=GP.LOG10_FGW_RANGE[0],
                    help="log10 f_lo of the GENERATIVE band (remedy-a study: -8.7)")
    ap.add_argument("--fhi", type=float, default=GP.LOG10_FGW_RANGE[1])
    args = ap.parse_args()
    GP.check_affinity()
    T_label, n_src = args.t, args.n
    band = (args.flo, args.fhi)
    band_default = tuple(band) == tuple(GP.LOG10_FGW_RANGE)

    print(f"\n=== GLACIER g2a-ii FORENSIC (T={T_label}, ncw={n_src}, band "
          f"[{args.flo}, {args.fhi}]{'' if band_default else ' <- REMEDY-A STUDY'}, "
          f"zero-noise, lane {GP.lane_tag()}) ===", flush=True)
    pop = GP.draw_population(args.seed, n_src=n_src, band_log10f=band)
    a_eff, _ = GP.a_eff_projection(pop["src"], band_log10f=band)
    G = GP.build_glacier_problem(T_label, pop)
    amo = G["amo"]
    import jax.numpy as jnp
    theta = np.asarray(amo["theta_truth"], float)
    ones = jnp.ones(amo["npsr"])
    t0 = time.time()
    data = amo["inject_delay"](jnp.asarray(theta), ones)
    print(f"  Asimov injection: {time.time()-t0:.1f}s")
    bf = GP.BackgroundFit(amo, band_log10f=band)
    it = amo["internals"]

    # ---- F1: template and residual identities under smask = 0 ----
    params = dict(bf.frozen)
    for k, v in zip(bf.theta_keys, theta):
        params[k] = v
    for dk, d in zip(bf.data_keys, data):
        params[dk] = d
    params["__pmask"] = ones
    params["__smask"] = jnp.zeros(n_src)
    tmax = max(float(jnp.max(jnp.abs(m(params)))) for m in it["msds"])
    rms = np.array([float(jnp.sqrt(jnp.mean(d ** 2))) for d in data])
    print(f"\n[F1] max|template(smask=0)| over all pulsars = {tmax:.3e}  "
          f"(H-BUG at the MaskedDelay level iff > 0)")
    print(f"[F1] injected-population delay RMS: median {np.median(rms)*1e9:.1f} ns, "
          f"max {rms.max()*1e9:.1f} ns, min {rms.min()*1e9:.1f} ns  "
          f"(the residual the fit should see)")

    # ---- F2: FtNmy under zeros / ones / absent ----
    t0 = time.time()
    p0_z, F_z = bf._data_terms(theta, data, ones, np.zeros(n_src))
    p0_o, F_o = bf._data_terms(theta, data, ones, np.ones(n_src))
    p0_n, F_n = data_terms_nokey(bf, theta, data, ones)
    nz, no, nn = [float(np.linalg.norm(v)) for v in (F_z, F_o, F_n)]
    print(f"\n[F2] ||FtNmy||: smask=zeros {nz:.6e} | smask=ones {no:.6e} | "
          f"ABSENT key {nn:.6e}   ({time.time()-t0:.1f}s)")
    print(f"[F2] ones vs absent max|diff| = {float(np.max(np.abs(F_o-F_n))):.3e} "
          f"(1.0*x==x identity; must be 0)")
    smask_honoured = nz > 1e3 * max(no, 1e-30)
    print(f"[F2] __smask reaches the template through kterms: "
          f"{'YES -- H-BUG REFUTED at this layer' if smask_honoured else 'NO -- H-BUG CONFIRMED'}")

    # ---- F3: decomposed profile, grid extended 2 dex below the smoke grid ----
    grid = np.round(np.arange(-17.6, -13.09, 0.1), 10)
    t0 = time.time()
    lnl, quads, occs, ahat, sig, edge = profile_from_terms(bf, F_z, grid)
    print(f"\n[F3] population profile over [{grid[0]}, {grid[-1]}] "
          f"({len(grid)} pts, {time.time()-t0:.1f}s):")
    i0 = int(np.argmin(np.abs(grid - (-14.6))))
    for i in range(0, len(grid), 4):
        print(f"     A={grid[i]:+.2f}  dquad={quads[i]-quads[0]:+12.4f}  "
              f"doccam={occs[i]-occs[0]:+12.4f}  dlnl={lnl[i]-lnl[0]:+12.4f}")
    print(f"[F3] quad rise floor->A0(-14.6): {quads[i0]-quads[0]:+.4f} nats*2 ; "
          f"occam rise: {occs[i0]-occs[0]:+.4f}")
    print(f"[F3] ahat = {ahat:.4f} +- {sig:.4f}{'  [EDGE]' if edge else ''}  "
          f"(A_eff drawn = {a_eff:.4f})")

    # ---- F4: positive control at the reference amplitude ----
    t0 = time.time()
    c = control_draw(bf)
    Fs = it["Fs_gwb"]
    data_ctrl = tuple(jnp.asarray(np.asarray(Fs[p]) @ c[p * bf.ngp:(p + 1) * bf.ngp])
                      for p in range(bf.npsr))
    rms_c = np.array([float(jnp.sqrt(jnp.mean(d ** 2))) for d in data_ctrl])
    _, F_c = bf._data_terms(theta, data_ctrl, ones, np.zeros(n_src))
    lnl_c, quads_c, occs_c, ahat_c, sig_c, edge_c = profile_from_terms(bf, F_c, grid)
    print(f"\n[F4] control draw (in-band GP at A0=-14.6, seed 4242): residual RMS median "
          f"{np.median(rms_c)*1e9:.1f} ns   ({time.time()-t0:.1f}s)")
    print(f"[F4] control ahat = {ahat_c:.4f} +- {sig_c:.4f}{'  [EDGE]' if edge_c else ''}")
    print(f"[F4] quad rise floor->A0: {quads_c[i0]-quads_c[0]:+.4f} ; "
          f"occam rise: {occs_c[i0]-occs_c[0]:+.4f}")
    print(f"[F4] THE VENUE'S INFORMATION: any amplitude gate at (T={T_label}, n={n_src}) "
          f"can honestly demand ~{sig_c:.2f} dex, no tighter.")

    btag = "" if band_default else f"_flo{str(args.flo).replace('-', 'm').replace('.', 'p')}"
    path = GP.bank_npz(f"g2forensic_T{T_label}_n{n_src}{btag}",
                       grid=grid, lnl_pop=lnl, quad_pop=quads, occam_pop=occs,
                       ahat_pop=ahat, sig_pop=sig, edge_pop=edge,
                       lnl_ctrl=lnl_c, quad_ctrl=quads_c, occam_ctrl=occs_c,
                       ahat_ctrl=ahat_c, sig_ctrl=sig_c, edge_ctrl=edge_c,
                       a_eff_log10=a_eff, a0_ref=bf.a0, T_label=T_label, n_src=n_src,
                       band_log10f=np.array(band), seed_pop=args.seed, seed_ctrl=4242,
                       ftnmy_norm_zeros=nz, ftnmy_norm_ones=no, ftnmy_norm_nokey=nn,
                       template_smask0_max=tmax, data_rms=rms, ctrl_rms=rms_c,
                       smask_honoured=smask_honoured,
                       orientation="lnl_* are 0.5*(quad-occam) up to a grid-constant; "
                                   "sig_ctrl is the venue's honest amplitude precision")
    print(f"\n  banked -> {path}")
    verdict = ("H-BUG" if not smask_honoured else
               ("H-FAINT" if (edge or abs(ahat - a_eff) > 0.15) else "FIT-OK"))
    print(f"\n=== FORENSIC READING (T={T_label}, n={n_src}): {verdict} ===", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
