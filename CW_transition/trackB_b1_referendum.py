"""Track B — B1 STEP 2: the TARGETED (f, mc) posterior referendum, at three conditioning tiers.

Deliverable R asked: with NO external information, does the fringe-marginalised posterior
concentrate at the needle? Answer f = Z_needle/(Z_needle+Z_plateau) = 6.9e-7, integrating the
6 SKY registration dims. R held frequency and chirp mass at truth and argued their Laplace
factors cancel between needle and plateau.

The pilot showed that argument is wrong for log10_mc (needle sigma ~1e-3 scaled vs plateau
~1.7), which makes R's blind verdict FIRMER, not weaker: including mc shrinks Z_needle.
But it also means the FORWARD claim -- "an EM counterpart supplies the 0.188 deg sky prior,
therefore certification becomes possible" (design-theorem lever iii) -- does NOT follow.
Removing the sky leaves (log10_fgw, log10_mc) as registration dims with their OWN volume
contest. This module runs it.

The object is R's, verbatim (`trackB_b1.B1Marg`): the count-once star-topology fringe-marginal
    lnL_marg(theta) = lnL_ref + sum_p LSE_n[ (LNL_p(n) - lnL_ref) + lnprior_p(n) ]
with the SKY of every loud source FIXED at truth (the counterpart baseline) and the ACTIVE
dims = (log10_fgw, log10_mc) x 3 loud = 6. Extrinsics are fixed at truth: the pilot's M2 shows
census P(true) is flat in them out to delta = 1e-2 scaled, i.e. they are not registration dims
and their Laplace factors DO cancel (the justification R asserted, here measured).

Three pre-registered conditioning tiers (sky given at truth throughout):

  Tier A  sky only.
          f box  = Earth-term posterior sigma (pilot M3).
          mc box = the model's generative prior, U(8.5, 9.5) -> half-width 1.0 scaled.

  Tier B  sky + EM period.  An electromagnetic counterpart with a photometric periodicity
          gives P_orb, hence f_gw = 2/P_orb. Adopted fractional precision sigma_P/P = 1e-3.
          JUSTIFICATION + why this is generous: the SMBHB-candidate periodicities in the
          literature are measured far more coarsely than this. OJ 287's ~12 yr cycle is the
          best-timed case (outburst epochs to ~days, ~1e-4 fractional) but it is a single,
          disputed object; the CRTS / PTF / Pan-STARRS candidate populations (Graham+ 2015,
          Charisi+ 2016, Liu+ 2019) have periods from only ~1.5-2 observed cycles, so
          sigma_P/P ~ 1e-2 at best, and most candidates do not survive red-noise reanalysis
          (Vaughan+ 2016). Adopting 1e-3 is therefore OPTIMISTIC by ~10x for a realistic
          counterpart and pessimistic only vs the single best-timed object. Reported as an
          adopted value with a sensitivity row over {1e-2, 1e-3, 1e-4}, not as a measurement.
          -> sigma(log10_fgw) = (1e-3)/ln10 dex.
          mc box = the model prior (a period gives NO chirp mass).

  Tier C  sky + period + luminosity distance.  The host redshift gives D_L; with the CW model's
          own amplitude relation
              log10 h = (5/3) log10 Mc + (2/3) log10 f - log10 D_L + const
          the array's measured h (pilot M3 posterior) plus D_L pins Mc:
              sigma(log10_mc) = (3/5) sqrt( sigma_logh^2 + ((2/3) sigma_logf)^2 + sigma_logDL^2 )
          with sigma_logDL = 0.005 dex (host z to ~1%). This is the tightest mc a counterpart
          can physically deliver -- it does NOT come from the pulsar terms.

Per tier: ln Z_needle by truth-placed local quadrature (LABELLED diagnostic truth use, as R1),
ln Z_box by tempered SMC with the needle sub-box excised (as R2), then
    f = Z_needle / (Z_needle + Z_box)
and the break-even (f, mc) box. PRE-REGISTERED (Matt, 2026-07-09): the loosest tier with
f >= 0.95 defines B1's legitimate conditioning.

Env: jug-gpu, autotune off, x64. Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_b1_core as C
import trackB_b1 as B1
import trackB_estimator as TE

CWT = "/home/mattm/projects/HSYMT/CW_transition"
I_MC, I_FGW = 3, 4
SIG_LOGDL = 0.005            # host redshift to ~1% -> log10 D_L, dex
SIG_P_OVER_P = 1e-3          # Tier B adopted EM-period fractional precision (see docstring)


# ============================================================
# tier boxes (SCALED units, half-widths, per loud source)
# ============================================================
def tier_boxes(P, tier, sig_p_over_p=SIG_P_OVER_P):
    """Return (half_f (n_loud,), half_mc (n_loud,)) in SCALED units, plus a provenance string."""
    scale = TE.phi_scale(P)
    m3 = np.load(f"{CWT}/b1_pilot_m3.npz")
    sig = m3["sigma_scaled"]
    sig_f_earth = np.array([sig[8 * k + I_FGW] for k in range(C.N_LOUD)])
    sig_h_earth = np.array([sig[8 * k + 5] for k in range(C.N_LOUD)]) * scale[5]     # dex
    sig_f_earth_dex = sig_f_earth * scale[I_FGW]

    # the generative prior on log10_mc is U(8.5, 9.5): half-width 0.5 dex
    half_mc_prior = np.full(C.N_LOUD, 0.5 / scale[I_MC])

    if tier == "A":
        half_f = 3.0 * sig_f_earth
        half_mc = half_mc_prior
        prov = ("f box = 3 x Earth-term posterior sigma (pilot M3); "
                "mc box = generative prior U(8.5,9.5)")
    elif tier == "B":
        sig_f_dex = sig_p_over_p / np.log(10.0)            # d log10 f = (df/f)/ln10
        half_f = np.full(C.N_LOUD, 3.0 * sig_f_dex / scale[I_FGW])
        half_f = np.minimum(half_f, 3.0 * sig_f_earth)     # never wider than the array's own
        half_mc = half_mc_prior
        prov = (f"f box = 3 x sigma from EM period (sigma_P/P = {sig_p_over_p:.0e} -> "
                f"{sig_f_dex:.2e} dex), capped at the Earth-term posterior; "
                f"mc box = generative prior")
    elif tier == "C":
        sig_f_dex = sig_p_over_p / np.log(10.0)
        half_f = np.full(C.N_LOUD, 3.0 * sig_f_dex / scale[I_FGW])
        half_f = np.minimum(half_f, 3.0 * sig_f_earth)
        sig_mc_dex = (3.0 / 5.0) * np.sqrt(sig_h_earth ** 2
                                           + ((2.0 / 3.0) * sig_f_dex) ** 2
                                           + SIG_LOGDL ** 2)
        half_mc = 3.0 * sig_mc_dex / scale[I_MC]
        half_mc = np.minimum(half_mc, half_mc_prior)
        prov = (f"f box as Tier B; mc box = 3 x sigma from h + f + D_L "
                f"(sigma_logh {np.round(sig_h_earth,4)} dex, sigma_logDL {SIG_LOGDL}) "
                f"-> sigma_mc {np.round(sig_mc_dex,4)} dex")
    else:
        raise ValueError(tier)
    return half_f, half_mc, prov


def active_index(k, j):
    """index into the 6-vector of active (f, mc) dims for loud source k."""
    return 2 * k + (0 if j == I_FGW else 1)


class TargetedMarg:
    """lnL_marg as a function of the 6 ACTIVE dims (log10_fgw, log10_mc) x 3 loud, in SCALED
    offsets from truth. Sky + extrinsics + all faint sources FIXED at truth."""

    def __init__(self, P, data, mask):
        self.P = P
        self.E = B1.B1Marg(P, data, mask)
        self.nd = P["n_dist"]
        self.src_truth = P["theta_truth"][self.nd:].copy()
        self.scale = TE.phi_scale(P)
        self.slots = np.array([8 * k + j for k in range(C.N_LOUD) for j in (I_FGW, I_MC)])
        self.sc = np.array([self.scale[j] for k in range(C.N_LOUD) for j in (I_FGW, I_MC)])
        self.D = len(self.slots)

    def src_from_x(self, X):
        X = np.atleast_2d(np.asarray(X, float))
        src = np.tile(self.src_truth, (X.shape[0], 1))
        src[:, self.slots] = self.src_truth[self.slots] + X * self.sc
        return src

    def __call__(self, X):
        return self.E.lnmarg(self.src_from_x(X))


def _lse(x):
    m = np.max(x)
    return m + np.log(np.sum(np.exp(x - m)))


# ============================================================
# R1-style: Z_needle by local quadrature (TRUTH-PLACED, labelled)
# ============================================================
def _fd_hessian(G, g0, hstep, D):
    pts = [np.zeros(D)]
    iph, imh = {}, {}
    for i in range(D):
        iph[i] = len(pts); x = np.zeros(D); x[i] = hstep[i]; pts.append(x.copy())
        imh[i] = len(pts); x = np.zeros(D); x[i] = -hstep[i]; pts.append(x.copy())
    ipp = {}
    for i in range(D):
        for j in range(i + 1, D):
            ipp[(i, j)] = len(pts); x = np.zeros(D)
            x[i] = hstep[i]; x[j] = hstep[j]; pts.append(x.copy())
    gv = G(np.array(pts)); g0b = gv[0]
    H = np.zeros((D, D))
    for i in range(D):
        H[i, i] = (gv[iph[i]] + gv[imh[i]] - 2 * g0b) / hstep[i] ** 2
    for i in range(D):
        for j in range(i + 1, D):
            H[i, j] = H[j, i] = (gv[ipp[(i, j)]] - gv[iph[i]] - gv[iph[j]] + g0b) / (hstep[i] * hstep[j])
    return -0.5 * (H + H.T)


def _coord_steps(G, g0, boxhalf, D, target=-0.4, lo=-1.5, hi=-0.1):
    """Per-coordinate FD step with dg in [lo, hi], closest to `target`. A step much larger than
    the needle's eigen-sigma samples the SHOULDER and fakes negative curvature -- so `target` is
    deliberately shallower than R1's -1 nat."""
    HLAD = np.geomspace(1e-8, 1.0, 22)
    pts = []
    for i in range(D):
        for h in HLAD:
            x = np.zeros(D); x[i] = min(h, 0.5 * boxhalf[i]); pts.append(x)
    dg = (G(np.array(pts)) - g0).reshape(D, len(HLAD))
    hstep = np.zeros(D)
    for i in range(D):
        hl = np.minimum(HLAD, 0.5 * boxhalf[i]); d = dg[i]
        good = np.where((d <= hi) & (d >= lo))[0]
        hstep[i] = hl[good[np.argmin(np.abs(d[good] - target))]] if len(good) else (
            hl[-1] if d[-1] > hi else hl[0])
    return hstep, dg


def z_needle(G, boxhalf, g0, tag):
    """Z_needle by truth-placed local quadrature. TRUTH USE IS DIAGNOSTIC AND LABELLED.

    Hardened after the Tier-C first pass produced 1/6 NEGATIVE curvature eigenvalue and failed the
    doubling gate (|dlnJ| = 0.566). Both symptoms had one cause: the coordinate FD step was chosen
    to give dg ~ -1 nat, i.e. 5-8x OUTSIDE the needle's eigen-sigma, so the second difference
    measured the shoulder. Hardening:
      * TWO-PASS Hessian: coordinate pass -> eigen-sigma -> re-choose steps at ~0.3 x the smallest
        eigen-sigma -> re-do. Richardson check H(h) vs H(h/2).
      * NEGATIVE eigenvalues are NEVER clipped. Clipping turns a non-curving direction into a very
        WIDE Gaussian and INFLATES Z_needle -- biased toward the needle winning, i.e. toward the
        answer that lets B1 proceed.
      * explicit 1-D profiles along every eigenvector, printed, so a negative direction is SEEN.
      * quadrature ladder 49 -> 97 -> 193 until the doubling gate passes.
    """
    D = G.D
    print(f"\n  -- Z_needle [{tag}] (truth-placed local quadrature; diagnostic truth use) --",
          flush=True)
    hstep, _ = _coord_steps(G, g0, boxhalf, D)
    print(f"     pass-1 coord steps (scaled): {np.array2string(hstep, precision=2)}", flush=True)
    H1 = _fd_hessian(G, g0, hstep, D)
    e1, _ = np.linalg.eigh(H1)
    sig_min = 1.0 / np.sqrt(max(e1.max(), 1e-300))
    # pass 2: every coordinate step at ~0.3 x the SHARPEST eigen-sigma (never above the pass-1 step)
    hstep2 = np.minimum(hstep, 0.3 * sig_min * np.ones(D))
    hstep2 = np.maximum(hstep2, 1e-9)
    print(f"     sharpest eigen-sigma {sig_min:.3e} -> pass-2 steps "
          f"{np.array2string(hstep2, precision=2)}", flush=True)
    Hpos = _fd_hessian(G, g0, hstep2, D)
    Hhalf = _fd_hessian(G, g0, 0.5 * hstep2, D)
    rich = float(np.max(np.abs(Hpos - Hhalf)) / max(np.max(np.abs(Hpos)), 1e-300))
    print(f"     Richardson |H(h)-H(h/2)|/max|H| = {rich:.3e} "
          f"-> {'STABLE' if rich < 0.1 else 'FD-NOISY'}", flush=True)

    evals, evecs = np.linalg.eigh(Hpos)
    npos = int((evals > 0).sum())
    print(f"     Hessian: {npos}/{D} positive-curvature eigs; "
          f"range [{evals.min():.3e}, {evals.max():.3e}]", flush=True)

    # ---- explicit 1-D profiles along every eigenvector (SEE the curvature) ----
    print(f"     1-D profiles of lnL_marg along each eigenvector (dg vs multiples of sigma_eff):",
          flush=True)
    sig_eff = np.array([1.0 / np.sqrt(e) if e > 0 else np.nan for e in evals])
    mult = np.array([-6, -3, -1, 1, 3, 6], float)
    prof = np.full((D, len(mult)), np.nan)
    for k in range(D):
        base = sig_eff[k] if evals[k] > 0 else 0.3 * np.min(boxhalf / (np.abs(evecs[:, k]) + 1e-30))
        s = mult * base
        cap = np.min(boxhalf / (np.abs(evecs[:, k]) + 1e-30))
        s = np.clip(s, -cap, cap)
        X = s[:, None] * evecs[:, k][None, :]
        prof[k] = G(X) - g0
        tagk = "POS" if evals[k] > 0 else "NEG/FLAT"
        print(f"       eig{k} lam={evals[k]:+.3e} [{tagk}] base={base:.2e}  dg = "
              + " ".join(f"{v:+8.2f}" for v in prof[k]), flush=True)

    # ---- log-scan: separate a genuine MICRO-DIP at truth from a macroscopic saddle ----
    # lnL_marg = lnL_ref + sum_p m_p with m_p >= 0. lnL_ref peaks at truth (zero-noise Asimov), but
    # every m_p GROWS as its pulsar de-registers, so a micro-offset can buy more fringe entropy than
    # it costs in lnL_ref. A negative curvature at 1e-6 is therefore physically possible while the
    # surface still falls by hundreds of nat at 1e-2. Distinguish the two by scanning decades.
    decades = np.geomspace(1e-7, 1.0, 15)
    print(f"     log-scan of dg along each eigenvector (|s| in scaled units, mean of +/-):",
          flush=True)
    print(f"       {'|s|':>9s} " + " ".join(f"{'eig'+str(k):>9s}" for k in range(D)), flush=True)
    logscan = np.full((len(decades), D), np.nan)
    for i, frac in enumerate(decades):
        row = []
        for k in range(D):
            cap = np.min(boxhalf / (np.abs(evecs[:, k]) + 1e-30))
            s = frac * cap
            v = G(np.array([+s * evecs[:, k], -s * evecs[:, k]])) - g0
            logscan[i, k] = 0.5 * (v[0] + v[1])
            row.append(logscan[i, k])
        print(f"       {decades[i]:9.1e} " + " ".join(f"{v:9.2f}" for v in row), flush=True)
    micro = logscan[0]                       # dg at |s| = 1e-7 x cap
    if npos < D:
        rising = [k for k in range(D) if evals[k] <= 0 and micro[k] > 0]
        print(f"     *** {D-npos} NON-POSITIVE curvature direction(s) at truth. Of these, "
              f"{len(rising)} actually RISE at |s|=1e-7 x cap (dg > 0). ***", flush=True)
        print(f"     Interpretation: negative curvature + dg<0 at every scanned scale = a MICRO-DIP "
              f"at truth inside a globally sharp peak (fringe-entropy gain beating lnL_ref loss at "
              f"sub-fringe offsets), NOT a saddle. Z_needle is then still a well-defined local "
              f"integral, but it is NOT a Laplace integral -- quadrature only.", flush=True)

    # ---- per-eigenvector quadrature: BRACKET the peak, then integrate inside it ----
    # The first Tier-A pass set the integration range from the eigenvalue (6 sigma, or the whole
    # box cap for the non-positive direction) and laid a uniform grid across it. For a peak of
    # width ~1e-5 inside a box of half-width ~4e-2 the trapezoid then returns nothing but the grid
    # spacing: sum lnJ fell by EXACTLY ln 2 per doubling (-53.0955 -> -53.7886 -> -54.4817), the
    # signature of an unresolved peak. Range must come from the FUNCTION, not from the curvature.
    DROP = 30.0                       # bracket where lnL_marg has fallen DROP nat below the peak

    def _edge(k, sign, cap):
        """smallest |s| along +/- evec_k where dg <= -DROP; bisect on a log-bracketed scan."""
        s = min(1e-7, 0.5 * cap)
        for _ in range(40):                      # expand
            if s >= cap:
                return cap, False                # never drops inside the box -> flag
            d = G((sign * s * evecs[:, k])[None])[0] - g0
            if d <= -DROP:
                break
            s *= 1.7
        lo, hi = s / 1.7, min(s, cap)
        for _ in range(18):                      # bisect
            mid = 0.5 * (lo + hi)
            d = G((sign * mid * evecs[:, k])[None])[0] - g0
            if d <= -DROP:
                hi = mid
            else:
                lo = mid
        return hi, True

    edges = np.zeros((D, 2)); closed = np.zeros(D, bool)
    for k in range(D):
        cap = np.min(boxhalf / (np.abs(evecs[:, k]) + 1e-30))
        em, okm = _edge(k, -1.0, cap)
        ep, okp = _edge(k, +1.0, cap)
        edges[k] = (em, ep); closed[k] = okm and okp
        print(f"       eig{k}: bracket [-{em:.3e}, +{ep:.3e}] (cap {cap:.2e}) "
              f"{'closed' if closed[k] else 'OPEN -> peak not contained in box'}", flush=True)

    def lnJ(nq):
        X = []; sv = []
        for k in range(D):
            s = np.concatenate([np.linspace(-edges[k, 0], 0, nq // 2, endpoint=False),
                                np.linspace(0, edges[k, 1], nq - nq // 2)])
            sv.append(s); X.append(s[:, None] * evecs[:, k][None, :])
        p = (G(np.concatenate(X)) - g0)
        out = []; o = 0
        for k in range(D):
            n = len(sv[k]); y = np.exp(p[o:o + n]); o += n
            out.append(np.log(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(sv[k]))))
        return np.array(out)

    Jprev = lnJ(129); dJ = np.inf
    for nq in (257, 513):
        Jcur = lnJ(nq)
        dJ = abs(Jcur.sum() - Jprev.sum())
        print(f"     quadrature nq={nq}: sum lnJ = {Jcur.sum():.4f}; doubling |dlnJ| = {dJ:.3f} "
              f"-> {'STABLE' if dJ < 0.2 else 'REFINE'}", flush=True)
        Jprev = Jcur
        if dJ < 0.2:
            break
    if dJ >= 0.2:
        print(f"     *** quadrature did NOT converge (|dlnJ| = {dJ:.3f}); ln Z_needle is a BOUND, "
              f"not a measurement ***", flush=True)
    if not closed.all():
        print(f"     *** {int((~closed).sum())} eigendirection(s) never fall {DROP:.0f} nat inside "
              f"the box: the needle is NOT contained in this tier's box ***", flush=True)
    lnZ_q = g0 + Jprev.sum()
    if npos == D:
        lnZ_lap = g0 + 0.5 * (D * np.log(2 * np.pi) - np.sum(np.log(evals)))
        print(f"     Laplace ln Z_needle = {lnZ_lap:.4f}; quadrature = {lnZ_q:.4f} "
              f"(diff {lnZ_q-lnZ_lap:+.3f} nat)", flush=True)
    else:
        lnZ_lap = np.nan
        print(f"     Laplace ln Z_needle = N/A (non-positive curvature); "
              f"quadrature = {lnZ_q:.4f}", flush=True)
    sig_needle = np.where(evals > 0, sig_eff, np.max(edges, axis=1) / 6.0)
    return (lnZ_q, lnZ_lap, evals, evecs, sig_needle, dJ, rich, npos, prof,
            edges, closed, logscan)


# ============================================================
# R2-style: Z_box by tempered SMC, needle excised
# ============================================================
def _resample(w, rng):
    N = len(w); pos = (rng.random() + np.arange(N)) / N
    return np.searchsorted(np.cumsum(w), pos)


def z_box(G, boxhalf, N, seed, needle_sig, n_mcmc=4, tag="", target_acc=0.25, max_mcmc=14):
    """Tempered SMC for ln Z over the tier box.

    STRENGTHENED (pre-registered, Matt 2026-07-09). The first Tier-A pass used a FIXED number of
    RW-MCMC sweeps per stage; acceptance collapsed to 0.13-0.14 at high beta and the two seeds
    disagreed by ~2.4 nat (R2's benchmark on the same machinery is 0.05 nat). Fix: keep sweeping
    at each stage until the sweep acceptance clears `target_acc`, up to `max_mcmc` sweeps, adapting
    the RW scale every sweep. Under-mixed particles at beta~1 bias logZ and inflate the seed
    scatter, which is exactly what the gate in mode_run tests.
    """
    rng = np.random.default_rng(seed)
    D = len(boxhalf)
    X = rng.uniform(-boxhalf, boxhalf, size=(N, D))
    L = G(X)
    beta = 0.0; logZ = 0.0; s = 0.5; stage = 0
    ess_hist = []; acc_hist = []
    while beta < 1.0 - 1e-12:
        def ess_at(db):
            lw = db * L; lw = lw - lw.max()
            w = np.exp(lw)
            return (w.sum() ** 2) / np.sum(w ** 2)
        hi = 1.0 - beta
        if ess_at(hi) >= 0.5 * N:
            db = hi
        else:
            lo = 0.0
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if ess_at(mid) >= 0.5 * N: lo = mid
                else: hi = mid
            db = lo
        lw = db * L
        logZ += _lse(lw) - np.log(N)
        w = np.exp(lw - _lse(lw)); ess = 1.0 / np.sum(w ** 2)
        idx = _resample(w, rng); X = X[idx]; L = L[idx]
        beta += db; stage += 1; ess_hist.append(ess)
        Cv = np.atleast_2d(np.cov(X.T))
        Cv = Cv + 1e-9 * np.eye(D) * max(np.trace(Cv) / D, 1e-30)
        try:
            Lc = np.linalg.cholesky(Cv)
        except np.linalg.LinAlgError:
            Lc = np.diag(np.sqrt(np.maximum(np.diag(Cv), 1e-30)))
        # beta-adaptive mixing: sweep until acceptance clears target_acc (cap max_mcmc), adapting
        # the RW scale after EVERY sweep so a collapsed acceptance can actually recover.
        accs = []
        nsw = 0
        while nsw < max_mcmc:
            Xp = X + s * (rng.normal(size=(N, D)) @ Lc.T)
            inbox = np.all(np.abs(Xp) <= boxhalf, axis=1)
            Lp = np.where(inbox, G(Xp), -np.inf)
            a = beta * (Lp - L)
            take = inbox & (np.log(rng.random(N)) < a)
            X[take] = Xp[take]; L[take] = Lp[take]
            a_i = float(take.mean()); accs.append(a_i); nsw += 1
            s = float(np.clip(s * np.exp(a_i - 0.234), 1e-6, 1.0))
            if nsw >= n_mcmc and np.mean(accs[-n_mcmc:]) >= target_acc:
                break
        acc = float(np.mean(accs[-n_mcmc:]))
        acc_hist.append((beta, acc, nsw))
        print(f"     [smc {tag} s{seed}] sweeps={nsw} "
              f"stage {stage}: beta={beta:.4f} ESS={ess:.0f} "
              f"acc={acc:.2f} Lmax={L.max():.1f} logZdens={logZ:.3f}", flush=True)
        np.savez(f"{CWT}/b1_ref_{tag}_s{seed}_ckpt.npz", beta=beta, logZ=logZ, X=X, L=L,
                 stage=stage, ess=np.array(ess_hist))
    excised = int(np.sum(np.all(np.abs(X) < 5.0 * needle_sig, axis=1)))
    lnV = float(np.sum(np.log(2 * boxhalf)))
    ah = np.array(acc_hist)                      # (stage, 3) = beta, acc, sweeps
    hi = ah[ah[:, 0] > 0.5]
    acc_hi = float(hi[:, 1].min()) if len(hi) else float(ah[:, 1].min())
    return dict(logZ_density=logZ, lnVbox=lnV, logZ=logZ + lnV, Lmax=float(L.max()),
                X=X, L=L, excised=excised, ess=np.array(ess_hist), stages=stage,
                acc_hist=ah, acc_hi=acc_hi, sweeps=int(ah[:, 2].sum()))


# ============================================================
def mode_run(args):
    print(f"=== B1 STEP 2: targeted (f, mc) referendum, tier {args.tier} ===", flush=True)
    P = C.build_b1_problem()
    tt = P["theta_truth"]
    mask1 = P["mask_one"]
    data = P["amo"]["inject_delay"](jnp.asarray(tt), jnp.asarray(mask1))
    if args.noise:
        nz = P["nd"].draw(seed=args.noise_seed)
        data = tuple(jnp.asarray(np.asarray(d) + np.asarray(n)) for d, n in zip(data, nz))
        print(f"  data = CW(truth) + noise draw seed {args.noise_seed}", flush=True)
    else:
        print(f"  data = zero-noise Asimov", flush=True)

    G = TargetedMarg(P, data, mask1)
    half_f, half_mc, prov = tier_boxes(P, args.tier, args.sig_p)
    boxhalf = np.empty(G.D)
    for k in range(C.N_LOUD):
        boxhalf[active_index(k, I_FGW)] = half_f[k]
        boxhalf[active_index(k, I_MC)] = half_mc[k]
    print(f"  TIER {args.tier}: {prov}", flush=True)
    print(f"  box half-widths (scaled): f {np.round(half_f,5)}  mc {np.round(half_mc,4)}", flush=True)
    print(f"  ln V_box (6 active dims) = {np.sum(np.log(2*boxhalf)):.4f}", flush=True)

    # gate the fast-full residual BEFORE any science uses it (bit-exactness, not a benchmark)
    rng_g = np.random.default_rng(11)
    Xg = rng_g.uniform(-1, 1, size=(3, G.D)) * 0.3
    G.E.gate_fast_full(G.src_from_x(np.vstack([np.zeros(G.D), Xg])))

    t0 = time.time(); g0 = float(G(np.zeros((1, G.D)))[0]); t_first = time.time() - t0
    lnl0 = float(P["amo"]["logL"](jnp.asarray(tt), data, jnp.asarray(mask1)))
    print(f"  lnL(truth) = {lnl0:.4f}; lnL_marg(truth) = {g0:.4f} "
          f"(fringe entropy {g0-lnl0:+.2f} nat); first eval {t_first:.0f}s", flush=True)
    t0 = time.time(); _ = G(np.zeros((G.E.T_CHUNK, G.D))); tw = time.time() - t0
    print(f"  warm: {tw/G.E.T_CHUNK*1e3:.0f} ms per lnL_marg", flush=True)

    (lnZn_q, lnZn_l, evals, evecs, sig_needle, dJ, rich, npos, prof,
     edges, closed, logscan) = z_needle(G, boxhalf, g0, f"tier{args.tier}")
    if args.needle_only:
        np.savez(f"{CWT}/b1_needlediag_tier{args.tier}.npz", g0=g0, lnZn_quad=lnZn_q,
                 lnZn_lap=lnZn_l, evals=evals, evecs=evecs, sig_needle=sig_needle,
                 doubling=dJ, richardson=rich, npos=npos, prof=prof, boxhalf=boxhalf,
                 edges=edges, closed=closed, logscan=logscan)
        print(f"\n  [needle-only] saved b1_needlediag_tier{args.tier}.npz", flush=True)
        return 0

    print(f"\n  -- Z_box [tier {args.tier}] (tempered SMC, needle excised) --", flush=True)
    runs = []
    sd = 0
    # SMC STRENGTHENING GATE (pre-registered, Matt 2026-07-09): seeds must agree to <= 0.3 nat AND
    # high-beta acceptance >= 0.25. If either fails, ADD SEEDS rather than quote. R2's benchmark on
    # this machinery is 0.05 nat, so a 2-nat scatter is a mixing failure, not irreducible noise.
    while sd < args.max_seeds:
        t0 = time.time()
        r = z_box(G, boxhalf, args.N, sd, sig_needle, n_mcmc=args.mcmc, tag=f"t{args.tier}",
                  target_acc=args.target_acc, max_mcmc=args.max_mcmc)
        r["wall"] = time.time() - t0
        runs.append(r); sd += 1
        print(f"     [smc s{sd-1}] ln Z_box = {r['logZ']:.3f} (density {r['logZ_density']:.3f} "
              f"+ lnV {r['lnVbox']:.3f}); stages {r['stages']} sweeps {r['sweeps']} "
              f"acc_hi {r['acc_hi']:.3f} Lmax {r['Lmax']:.1f} excised {r['excised']} "
              f"({r['wall']:.0f}s)", flush=True)
        if sd < args.seeds:
            continue
        lz_ = np.array([q["logZ"] for q in runs])
        spread = float(lz_.max() - lz_.min())
        acc_hi = float(min(q["acc_hi"] for q in runs))
        ok = (spread <= args.gate_spread) and (acc_hi >= args.target_acc)
        print(f"     [GATE smc] seed spread {spread:.3f} nat (<= {args.gate_spread}) ; "
              f"min high-beta acc {acc_hi:.3f} (>= {args.target_acc}) -> "
              f"{'PASS' if ok else 'FAIL -> adding a seed'}", flush=True)
        if ok:
            break
    lz = np.array([r["logZ"] for r in runs])
    nseed = len(lz)
    lnZp = lz.mean()
    lnZp_e = lz.std(ddof=1) / np.sqrt(nseed) if nseed > 1 else np.inf   # s.e.m. on the mean
    spread = float(lz.max() - lz.min())
    acc_hi = float(min(r["acc_hi"] for r in runs))
    gate_ok = (spread <= args.gate_spread) and (acc_hi >= args.target_acc)
    print(f"     ln Z_box = {lnZp:.3f} +- {lnZp_e:.3f} (s.e.m., {nseed} seeds); "
          f"spread {spread:.3f} nat; min high-beta acc {acc_hi:.3f}", flush=True)

    d = lnZn_q - lnZp
    f = 1.0 / (1.0 + np.exp(-d))
    sig_d = lnZp_e                                   # Z_needle is deterministic quadrature
    sig_f = f * (1.0 - f) * sig_d                    # delta method through the logistic
    print(f"\n  === REFERENDUM, TIER {args.tier} ===", flush=True)
    print(f"  ln Z_needle = {lnZn_q:.3f} (quad) / {lnZn_l} (Laplace)", flush=True)
    print(f"  ln Z_box    = {lnZp:.3f} +- {lnZp_e:.3f}", flush=True)
    print(f"  ln Z_needle - ln Z_box = {d:+.3f} +- {sig_d:.3f} nat", flush=True)
    print(f"  >>> f = {f:.6g} +- {sig_f:.4g}   (2-sigma band "
          f"[{max(f-2*sig_f,0):.4g}, {min(f+2*sig_f,1):.4g}])", flush=True)
    lam = np.exp(d / G.D)
    print(f"  BREAK-EVEN: shrink every active dim by {lam:.4g}x "
          f"-> f box +-{np.round(half_f*lam,6)}, mc box +-{np.round(half_mc*lam,6)} scaled", flush=True)
    print(f"              (mc break-even in dex: {np.round(half_mc*lam*0.5,6)})", flush=True)

    # PASS/FAIL WITH ERROR BARS (pre-registered, amendment 2)
    if not gate_ok or not np.isfinite(sig_f):
        verdict = ("NO VERDICT -- SMC gate FAILED (seed spread / high-beta acceptance). "
                   "Strengthen further; do NOT branch.")
    elif f - 2 * sig_f >= 0.95:
        verdict = "PASS (f - 2sigma >= 0.95): this tier is legitimate conditioning for B1"
    elif f + 2 * sig_f < 0.95:
        verdict = "FAIL (f + 2sigma < 0.95): the targeted posterior does NOT concentrate here"
    else:
        verdict = ("STRADDLE (f +- 2sigma spans 0.95): NO BRANCH NAMED. Strengthen "
                   "(more seeds / sweeps) and re-run.")
    print(f"  PRE-REGISTERED VERDICT: {verdict}", flush=True)
    if not gate_ok:
        print(f"  (quadrature doubling |dlnJ| = {dJ:.3f}; Richardson {rich:.2e}; "
              f"{npos}/{G.D} positive-curvature eigs)", flush=True)
    np.savez(f"{CWT}/b1_referendum_tier{args.tier}{'_noisy' if args.noise else ''}.npz",
             tier=args.tier, g0=g0, lnl0=lnl0, lnZn_quad=lnZn_q, lnZn_lap=lnZn_l,
             lnZbox=lnZp, lnZbox_err=lnZp_e, f=f, sigma_f=sig_f, nseed=nseed,
             seed_spread=spread, acc_hi=acc_hi, gate_ok=gate_ok, verdict=verdict,
             lnZbox_all=lz, boxhalf=boxhalf, evals=evals,
             sig_needle=sig_needle, doubling=dJ, richardson=rich, npos=npos, prof=prof,
             edges=edges, closed=closed, logscan=logscan, breakeven_lambda=lam,
             half_f=half_f, half_mc=half_mc, prov=prov, N=args.N, seeds=args.seeds)
    print(f"  saved b1_referendum_tier{args.tier}.npz", flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="A", choices=["A", "B", "C"])
    ap.add_argument("--noise", action="store_true")
    ap.add_argument("--noise_seed", type=int, default=90001)
    ap.add_argument("--N", type=int, default=256)
    ap.add_argument("--seeds", type=int, default=2, help="minimum seeds before the gate is tested")
    ap.add_argument("--max_seeds", type=int, default=5, help="add seeds until the gate passes")
    ap.add_argument("--mcmc", type=int, default=4, help="minimum MCMC sweeps per SMC stage")
    ap.add_argument("--max_mcmc", type=int, default=14, help="cap on sweeps per stage")
    ap.add_argument("--target_acc", type=float, default=0.25, help="high-beta acceptance floor")
    ap.add_argument("--gate_spread", type=float, default=0.3, help="max seed spread, nat")
    ap.add_argument("--sig_p", type=float, default=SIG_P_OVER_P)
    ap.add_argument("--needle_only", action="store_true",
                    help="run only the Z_needle diagnostic (Hessian + eigen-profiles), skip SMC")
    a = ap.parse_args()
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    sys.exit(mode_run(a))
