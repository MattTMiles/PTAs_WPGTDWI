#!/usr/bin/env python
"""SPARK -- the cascade's launch criterion. AVALANCHE's proposed successor, executed.

WHY THIS EXISTS
---------------
AVALANCHE (reports/AVALANCHE_cascade.md) stopped the multi-source growing-list campaign
pre-flight. One of its three grounds was that the cascade's [E]->[D] feedback edge does not
exist in code: the incumbent seeder

    TE.seed_scan -> TE.build_earth_single -> build_fisher_amortised(msd_factory=EarthDelay)

is the FULLY-DECOHERED limit BY CONSTRUCTION (TE.EarthDelay docstring: "Earth-term-only ...
the fully-decohered coherence limit"; make_phase_connected_binary(pulsarterm=False)). Certified
pulsar terms cannot enter it, so [D] returns the same source list at every iteration and the
loop is inert for a CODE reason, not a physics one.

SPARK builds the missing path and prices the cascade.

S1 -- BUILD the pulsar-term-coherent DETECTION statistic (in the DETECTOR, not the fitter).
   The mask machinery already exists but was never wired into a detector:
   trackB_b1_core.MaskedDelay -- delay_p = d_earth + m_p*(d_full - d_earth), m_p read from
   params[PMASK][ipsr] at RUNTIME (gated: m_p=1 == MultiSourceDelay(pterm=True) 0.00;
   m_p=0 == TE.EarthDelay 0.00). SPARK wires it into an ncw=1 detector so the certified set
   enters through (pmask, Ldist) as RUNTIME args -> ONE compiled graph for every certified set.

   SOFT / q-WEIGHTED per spec 3, NEVER hard-locked: m_p = q_p in [0,1] (the pulsar's own
   certification confidence), not a 0/1 lock. Uncertified pulsars sit at m_p=0 -- decohered,
   NOT pinned to a wrong MAP fringe (the hard-lock failure that IGNITE-2 retired).

   GATE g0a: SPARK's re-wired EarthDelay F-stat (data as a runtime arg) == TE.make_fstat
             verbatim. Proves the re-wiring did not move the incumbent.
   GATE g0b: THE USER'S GATE -- the coherent detector at ZERO certified pulsars (pmask=0)
             reproduces the EarthDelay F-stat. Bar: EMBER 2.2(b) adopted -- discrete exact,
             continuous < 1e-6 (NOT bit-identity: different builders, different graphs).
   GATE g0c: the data-driven loud-source list (F-stat + sky-exclusion NMS, TRUTH-BLIND in
             selection) is reproduced from the coherent detector's own pmask=0 map: 3/3 loud.

S2 -- THE ARITHMETIC. Gain vs cost, per igniter, on the faint reservoir.

   GAIN (measured): for each of the 13 faint sources (the reservoir, log10_h=-14.25), the
   detection statistic 2F at its TRUE (sky, freq) -- ORACLE-ANCHORED -- profiled over
   (cos_inc, log10_h, phase0, psi), at certified states:
       s0    : pmask = 0            (the incumbent seeder; the cascade's iteration 0)
       sC(g) : pmask = q_p on the igniter g's banked certified set, 0 elsewhere; the certified
               pulsars' distances at L_true (a correct certification KNOWS L_p -- Arm B:
               L_true = L0 + (n_true+u_true)*dL, so certification is what pins it)
       sMAX  : pmask = 1 on ALL 116 at L_true -- THE CEILING. Every pulsar in the array
               perfectly certified. Unreachable by construction; it bounds what coherence
               could EVER buy.

   COST (measured): the null-calibrated 2F floor at matched state + the lnK growth per added
   source (RECUT adopt(): zf<=0.20 -> Gumbel, else emp q95 + bootstrap; zero-fraction is a
   REQUIRED column).

   VERDICT: CASCADE-ALIVE if the gain clears the cost anywhere; else CASCADE-DEAD + shortfall.

THE A-FORTIORI STRUCTURE (why a DEAD verdict here is unconditional). Every convention is
chosen to FAVOUR the cascade:
   (i)   ORACLE-ANCHORED -- 2F is read at the faint source's TRUE sky/freq. A real detector
         must search; searching can only do worse (B0.2: cold-start lands ~0.5 scaled vs the
         ~1e-4 needed; the needle is ~1e-5 deg wide).
   (ii)  sMAX certifies ALL 116 pulsars at their TRUE distances -- the coherence ceiling.
   (iii) The floor is read WITHOUT a trials factor (no search -> no max-over-grid penalty).
         The searched floor is also measured and is strictly higher.
   (iv)  Measured on the BASE census (3 loud + 13 faint, no eccentric harmonics). The
         igniters' harmonics ADD confusion power to the data, which only makes faint
         detection harder; the igniter enters here through its CERTIFIED SET, which is the
         lever the cascade actually pulls (channels -> certs -> coherence).
If the gain fails to clear the cost under ALL FOUR, no real cascade can succeed.

SCOPE OF INFERENCE. 116-pulsar MOCK array (AXIS, 1440 MHz single-frequency, 10.15(a)); no
real TOA is touched, the residuals ARE the injected CW+CURN. Real pulsar positions, real TOA
uncertainties, real published distance priors, drawn white + per-pulsar red + GWB noise.
A verdict here is a statement about an approximately-real IPTA's geometry and noise budget.

CONVENTIONS. cpus-per-task=8 is part of the seed (NoiseDrawer eigh basis rotates with BLAS
threads). T=30 has NO canonical L_gwb bank -> force_recompute_lgwb (IGNITE-2's path).
Lean-npz: raw statistics banked per state, never verdicts.

Run:  spark.py g0      (S1 build + gate chain)
      spark.py s2      (the ledger; refuses to run unless g0 passed and is banked)
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, time, argparse, glob
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir",
                  os.environ.get("HARBOR_JAXCACHE_IN", "/home/mattm/.cache/jax_stagec_cache"))
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

HSYMT = "/home/mattm/projects/HSYMT"
for _p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite"):
    sys.path.insert(0, f"{HSYMT}/{_p}")

import discovery as ds
import trackB_b1_core as C
import trackB_estimator as TE
import ignite as IG
from stagec_fisher import (build_fisher_amortised, make_geometry_injection,
                           N_COMPONENTS, LOG10_EQUAD)
from cw_helpers import build_enterprise_pta, compute_mode_spacing

OUT = f"{HSYMT}/SPARK_results"
REPORTS = f"{HSYMT}/reports"
I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)

# ---- the cell SPARK runs at: the CHORUS igniter cell, verbatim from the banks ----
H_CELL = -13.25          # loud strain
T_LABEL = 30             # CHORUS T_label
TIER = "lit"
GEO_SRC = None           # geo_id = -1, the fiducial POP sky (banked igniters are geo_id=-1)
N_LOUD = 3               # C.N_LOUD; sources 0..2 are loud, 3..15 are THE RESERVOIR
N_FAINT = 13
H_ABSENT = -18.0         # "source not in the model" (TE.seed_scan's own ll0 convention)

# igniter variants: the CHORUS switch-on set (CHORUS 2.1 -- 1 ecc member -> e=0.5;
# 2+ -> e=0.3). NOT the refuted "one e=0.3 member" headline.
IGNITERS = ["m1e07", "m1e05", "m2e03"]

ALPHA = 0.05
ZF_GATE = 0.20           # RECUT recut_surface.py:56
C_SD = 2.80              # RECUT gate G7
Z_ALPHA = 2.9701952521018403
BOOT = 4000
BOOT_SEED = 20260716

GATE_CONT = 1e-6         # EMBER 2.2(b) adopted bar for continuous columns


def force_recompute_lgwb(C):
    """kindle_loop.py:47-71, inlined verbatim (importing kindle_loop would drag the whole
    loop's dependency set in for a 6-line wrapper).

    There is NO canonical b1_L_gwb bank for the T=30 geometry: the canonical bank is
    (3248,3248) for the T=15 ANCHOR array, build_ignite_problem(30) needs (5336,5336), and
    load_or_build_L_gwb hits its shape-mismatch raise. IGNITE-2/KINDLE's T=30 runs used the
    RECOMPUTE path (np.linalg.eigh of Phi_gwb at cpus=8), which 10.14(g)/10.16(e) state is
    the convention every banked noisy number was drawn under and is bit-identical AT cpus=8.
    Point the loader at a non-existent path so it takes that branch. cpus-per-task MUST be 8
    (the noise-draw thread-count hazard); the sbatch sets it and the wrapper prints the
    provenance so the log proves the path taken.
    """
    _orig = C.load_or_build_L_gwb

    def _recompute(Phi_gwb, path, strict=False):
        L, prov = _orig(Phi_gwb, "/nonexistent_spark_force_recompute.npz", strict=False)
        print(f"[SPARK] L_gwb via {prov}", flush=True)
        return L, prov

    C.load_or_build_L_gwb = _recompute


# ============================================================
# S1 -- the detectors
# ============================================================
def build_detectors(P, verbose=True):
    """The ncw=1 twins, built on IDENTICAL pulsars / frozen noise / GP component counts so
    the g0 gate is apples-to-apples.

      earth1 : msd_factory=TE.EarthDelay   -- the INCUMBENT seeder's likelihood.
      coh1   : C.build_b1_amortised        -- MaskedDelay; pmask is a RUNTIME arg.

    NOTE (a real inconsistency in the incumbent, reported not fixed): TE.build_earth_single
    hard-codes the DEFAULT components (N_COMPONENTS=14, rn=30) regardless of T. At T=30 the
    problem's own GP counts are span-scaled (ext_diag), so TE.build_earth_single's noise model
    does NOT match a T=30 problem's data. The seeder was only ever run at T=15, where this is
    a no-op. SPARK builds BOTH detectors at the PROBLEM's scaled counts, so the two detectors
    agree with each other AND with the data they are fed.
    """
    t0 = time.time()
    ent, disco = P["ent_psrs"], P["disco_psrs"]
    ed = P.get("ext_diag") or {}
    gwb_comp = int(ed.get("gwb_comp", N_COMPONENTS))
    rn_comp = int(ed.get("rn_comp", 30))
    if verbose:
        print(f"[SPARK] detectors at PROBLEM component counts: gwb={gwb_comp} rn={rn_comp}",
              flush=True)
    pta1, cwb1, _ = build_enterprise_pta(ent, 1, components=N_COMPONENTS,
                                         log10_equad=LOG10_EQUAD)
    inj1 = make_geometry_injection(pta1, ent, 1, cwb1, seed=1000,
                                   h_range=(TE.EQUAL_H, TE.EQUAL_H))
    earth1 = build_fisher_amortised(disco, 1, inj1, cwb1, components=gwb_comp,
                                    rn_components=rn_comp, msd_factory=TE.EarthDelay)
    if verbose:
        print(f"[SPARK]   earth1 (EarthDelay, incumbent) {time.time()-t0:.1f}s", flush=True)
    t1 = time.time()
    coh1 = C.build_b1_amortised(disco, 1, inj1, cwb1, components=gwb_comp,
                                rn_components=rn_comp)
    if verbose:
        print(f"[SPARK]   coh1   (MaskedDelay, SPARK)    {time.time()-t1:.1f}s", flush=True)
    return earth1, coh1


def _adam_profile(lnl_fn):
    """TE.make_fstat's profile, VERBATIM (same init, same 40 fori_loop steps, same clip,
    same fscale). Shared by both detectors so the gate isolates the MSD PATH and nothing else.
    """
    flo = jnp.asarray(TE._FREE_LO); fhi = jnp.asarray(TE._FREE_HI)
    fscale = jnp.asarray([0.5, 1.0, 3.14, 1.57])

    def profile(*carry):
        g = jax.grad(lambda fr: lnl_fn(*carry, fr))
        free = jnp.array([0.0, -14.0, 3.14, 1.57])
        m = jnp.zeros(4); v = jnp.zeros(4)
        b1, b2, eps = 0.9, 0.999, 1e-8

        def body(t, c):
            free, m, v = c
            gr = jnp.nan_to_num(g(free))
            m = b1 * m + (1 - b1) * gr
            v = b2 * v + (1 - b2) * gr * gr
            mh = m / (1 - b1 ** t); vh = v / (1 - b2 ** t)
            free = jnp.clip(free + 0.05 * fscale * mh / (jnp.sqrt(vh) + eps), flo, fhi)
            return (free, m, v)

        free, _, _ = jax.lax.fori_loop(1, 41, body, (free, m, v))
        return lnl_fn(*carry, free), free
    return profile


def make_fstat_earth(earth1):
    """TE.make_fstat with (Ld, data) as RUNTIME args (TE bakes P['data_obs'] / P['L0'] at
    trace time; SPARK draws its own realisations, so they must be runtime)."""
    fl = earth1["fisher_logL"]

    def lnl_at(Ld, data, cosgt, gwphi, lfgw, free):
        cw = jnp.array([cosgt, gwphi, free[0], TE.SEED_MC, lfgw, free[1], free[2], free[3]])
        return fl(jnp.concatenate([Ld, cw]), data)

    prof = _adam_profile(lnl_at)
    return jax.jit(jax.vmap(prof, in_axes=(None, None, 0, 0, 0)))


def make_fstat_coh(coh1):
    """SPARK's COHERENT detector. Identical to make_fstat_earth in every respect EXCEPT the
    msd path: MaskedDelay(pmask) instead of EarthDelay. (Ld, pmask, data) are runtime ->
    one compiled graph for every certified set."""
    fl = coh1["logL"]

    def lnl_at(Ld, pmask, data, cosgt, gwphi, lfgw, free):
        cw = jnp.array([cosgt, gwphi, free[0], TE.SEED_MC, lfgw, free[1], free[2], free[3]])
        return fl(jnp.concatenate([Ld, cw]), data, pmask)

    prof = _adam_profile(lnl_at)
    return jax.jit(jax.vmap(prof, in_axes=(None, None, None, 0, 0, 0)))


def ll0_earth(earth1, Ld, data, lfgw_ref):
    """No-signal reference (TE.seed_scan's convention: log10_h = -18)."""
    cw = jnp.array([0., 0., 0., TE.SEED_MC, lfgw_ref, H_ABSENT, 0., 0.])
    return float(earth1["fisher_logL"](jnp.concatenate([jnp.asarray(Ld), cw]), data))


def ll0_coh(coh1, Ld, pmask, data, lfgw_ref):
    cw = jnp.array([0., 0., 0., TE.SEED_MC, lfgw_ref, H_ABSENT, 0., 0.])
    return float(coh1["logL"](jnp.concatenate([jnp.asarray(Ld), cw]),
                              data, jnp.asarray(pmask)))


def scan_grid(fstat, args_prefix, CG, GP, LF, chunk=256):
    """Chunked evaluation of a profiled statistic over a flat grid."""
    npt = len(CG)
    stat = np.empty(npt); free = np.empty((npt, 4))
    for c0 in range(0, npt, chunk):
        c1 = min(c0 + chunk, npt)
        s, fr = fstat(*args_prefix, jnp.asarray(CG[c0:c1]), jnp.asarray(GP[c0:c1]),
                      jnp.asarray(LF[c0:c1]))
        stat[c0:c1] = np.asarray(s); free[c0:c1] = np.asarray(fr)
    return stat, free


def flat_grid(P):
    """The seeder's own grid, at the PROBLEM's span (TE.seed_scan hard-codes T=6.992e8, the
    T=15 span, 'to avoid recompute'; at T=30 that is the wrong df. SPARK uses the real span)."""
    T = float(ds.getspan(P["disco_psrs"]))
    lfgws = TE.freq_grid(T)
    cosgt, gwphi = TE.sky_grid(TE.SEED_NSIDE)
    CG, GP, LF = [], [], []
    for i in range(len(lfgws)):
        for j in range(len(cosgt)):
            CG.append(cosgt[j]); GP.append(gwphi[j]); LF.append(lfgws[i])
    return (np.array(CG), np.array(GP), np.array(LF), lfgws, cosgt, gwphi, T)


# ============================================================
# realisation draw (Arm B, matching IGNITE/CHORUS exactly)
# ============================================================
def draw_realisation(P, noise_seed, dist_seed, faint_present=True):
    """Arm-B realisation on the base census. faint_present=False -> the 13 faint sources are
    REMOVED from the injected data (log10_h -> -18) but the 3 loud remain: SPARK's null for
    the detection floor (the reservoir is absent; everything else is matched)."""
    G = IG.cell_geometry(P, GEO_SRC, H_CELL, TIER)
    dL, EV = G["dL"], G["EV"]
    nd, AI = P["n_dist"], P["AI"]
    L_true, n_true = IG.draw_true_distances_tier(P, dL, EV, seed=dist_seed, tier=TIER)
    theta_true = G["theta"].copy()
    theta_true[AI] = L_true
    if not faint_present:
        for i in range(N_LOUD, P["ncw"]):
            theta_true[nd + 8 * i + I_H] = H_ABSENT
    one = jnp.ones(P["npsr"])
    clean = P["amo"]["inject_delay"](jnp.asarray(theta_true), one)
    noise = P["nd"].draw(noise_seed)
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n)) for c, n in zip(clean, noise))
    return dict(data=data, theta_true=theta_true, L_true=L_true, n_true=n_true,
                dL=dL, EV=EV, FT=G["FT"], theta_cell=G["theta"])


def faint_targets(P, theta_true):
    """The reservoir's true (sky, freq, h) -- the ORACLE anchors."""
    nd = P["n_dist"]
    out = []
    for i in range(N_LOUD, P["ncw"]):
        b = nd + 8 * i
        out.append(dict(idx=i, cosgt=float(theta_true[b + I_COSGT]),
                        gwphi=float(theta_true[b + I_GWPHI]),
                        lfgw=float(theta_true[b + I_FGW]),
                        log10_h=float(theta_true[b + I_H]),
                        log10_mc=float(theta_true[b + I_MC])))
    return out


# ============================================================
# certified sets from the banked igniters (re-derived from RAW, per lean-npz discipline)
# ============================================================
def chorus_floor(tag):
    z = np.load(f"{REPORTS}/recut_chorus.npz", allow_pickle=True)
    tags = [str(t) for t in np.asarray(z["tags"]).ravel()]
    cols = [str(c) for c in np.asarray(z["cols"]).ravel()]
    t = np.asarray(z["table"])
    i = tags.index(tag)
    return {c: float(t[i, j]) for j, c in enumerate(cols)}


def certified_sets(igniter, tier=TIER):
    """Re-derive the certified set from RAW columns (dlnL_det, lnK, qmax) at the ADOPTED
    criterion-v2.2 floor -- 'do not trust the banked bool' (ember_anatomy.py:26).
    Returns one entry per banked realisation."""
    tag = f"{igniter}_{tier}"
    fl = chorus_floor(tag)
    floor = fl["floor"]
    fs = sorted(glob.glob(f"{HSYMT}/CHORUS_results/ch_sig_{igniter}_{tier}_*.npz"))
    rows = []
    for f in fs:
        z = np.load(f, allow_pickle=True)
        dlnL = np.asarray(z["dlnL_det"], float)
        lnK = np.asarray(z["lnK"], float)
        qmax = np.asarray(z["qmax"], float)
        on_true = np.asarray(z["on_true"], bool)
        bar = np.maximum(lnK, floor)
        det = dlnL > bar
        cert = det & (qmax > 0.9)
        rows.append(dict(path=f, cert=cert, qmax=qmax, on_true=on_true,
                         n_cert=int(cert.sum()),
                         n_cert_true=int((cert & on_true).sum()),
                         n_cert_wrong=int((cert & ~on_true).sum()),
                         noise_seed=int(z["noise_seed"]), dist_seed=int(z["dist_seed"]),
                         n_active=int(z["n_active"])))
    return dict(tag=tag, floor=floor, rows=rows,
                corr_banked=fl["corr"], wrong_banked=fl["wrong"])


# ============================================================
# floors (RECUT adopt(), verbatim semantics)
# ============================================================
def gumbel_floor(x):
    """RECUT's Gumbel fit. Guarded: a degenerate sample (all-equal, e.g. a 100% zero point
    mass) makes gumbel_r.fit's bracket search overflow. Return NaN rather than crash -- the
    zero-fraction gate sends any such cell to emp_q95 regardless, and a NaN in the banked
    'gumbel' column is the honest record that the fit was not defined here."""
    from scipy.stats import gumbel_r
    x = np.asarray(x, float)
    n = len(x)
    if n < 2 or not np.isfinite(x).all() or np.ptp(x) <= 0:
        return np.nan, np.nan, np.nan, np.nan, n
    try:
        mu, beta = gumbel_r.fit(x)
    except (OverflowError, RuntimeError, FloatingPointError):
        return np.nan, np.nan, np.nan, np.nan, n
    return mu + beta * Z_ALPHA, mu, beta, C_SD * beta / np.sqrt(n), n


def boot_sd(x, q=1 - ALPHA):
    rng = np.random.default_rng(BOOT_SEED)
    x = np.asarray(x, float); n = len(x)
    return float(np.std([np.quantile(rng.choice(x, n, replace=True), q) for _ in range(BOOT)]))


def adopt(x):
    """RECUT recut_surface.py:101 adopt(), semantics verbatim. Zero-fraction is a REQUIRED
    column, not a caveat (ANCHOR 4)."""
    x = np.asarray(x, float)
    zf = float(np.mean(x == 0.0))
    gu, mu, beta, gu_sd, n = gumbel_floor(x)
    emp = float(np.quantile(x, 1 - ALPHA))
    esd = boot_sd(x)
    if zf <= ZF_GATE:
        return dict(floor=float(gu), err=float(gu_sd), estimator="gumbel", zf=zf,
                    gumbel=float(gu), gumbel_sd=float(gu_sd), emp=emp, emp_sd=esd, n=n)
    return dict(floor=emp, err=esd, estimator="emp_q95", zf=zf,
                gumbel=float(gu), gumbel_sd=float(gu_sd), emp=emp, emp_sd=esd, n=n)


# ============================================================
# lnK growth per added source (pure numpy -- no GPU, no likelihood)
# ============================================================
def lnK_vs_nsrc(P, theta_cell):
    """dL = MIN over the source list of the mode spacing (forge_b1.apply_geometry:91-94), so
    adding a source can only SHRINK dL -> more fringes in the prior window -> K_counted grows.
    Measured exactly by restricting the min to the first n sources. CPU, seconds."""
    nd = P["n_dist"]
    src = theta_cell[nd:]
    cw = [dict(cos_gwtheta=src[8 * i + I_COSGT], gwphi=src[8 * i + I_GWPHI],
               log10_fgw=src[8 * i + I_FGW]) for i in range(P["ncw"])]
    rows = []
    for n in range(1, P["ncw"] + 1):
        dL = np.array([min(compute_mode_spacing(c["cos_gwtheta"], c["gwphi"], c["log10_fgw"],
                                                P["disco_psrs"][a].pos) for c in cw[:n])
                       for a in range(P["npsr"])])
        EV = C.build_EV(P, dL)
        FT = C.FringeTables(P, EV, dL, prior_key=TIER)
        K = np.asarray(FT.K_counted, float)
        rows.append(dict(n_src=n, K_sum=float(K.sum()), lnK_mean=float(np.log(np.maximum(K, 1)).mean()),
                         lnK_med=float(np.median(np.log(np.maximum(K, 1)))),
                         dL_med=float(np.median(dL))))
    return rows


# ============================================================
# MODE g0 -- S1 build + the gate chain
# ============================================================
def mode_g0():
    print("=" * 78, flush=True)
    print("SPARK S1 -- the pulsar-term-coherent DETECTION statistic: BUILD + GATE", flush=True)
    print("=" * 78, flush=True)
    force_recompute_lgwb(C)
    t0 = time.time()
    P = IG.build_ignite_problem(T_LABEL, verbose=True)
    print(f"[SPARK] problem T={T_LABEL} built {time.time()-t0:.1f}s "
          f"(npsr={P['npsr']} ncw={P['ncw']} n_theta={P['n_theta']})", flush=True)

    earth1, coh1 = build_detectors(P)
    fe = make_fstat_earth(earth1)
    fc = make_fstat_coh(coh1)

    R = draw_realisation(P, noise_seed=7560000, dist_seed=17560000, faint_present=True)
    data = R["data"]
    P["data_obs"] = data            # TE.make_fstat reads this at trace time
    # TE.loud_truth / TE.make_fstat read P["config"], a Track B build_problem key that
    # build_ignite_problem does not set. "pop" is the correct tag -- ASSERTED, not assumed:
    # TE.CONFIGS["pop"] must be C.POP, the population every Track B deliverable draws from.
    assert TE.CONFIGS["pop"] == C.POP, (
        f"STOP: TE.CONFIGS['pop']={TE.CONFIGS['pop']} != C.POP={C.POP}; the seeder's "
        f"population tag has drifted from the loop's and loud_truth would name the wrong "
        f"sources.")
    P["config"] = "pop"

    CG, GP, LF, lfgws, cosgt, gwphi, T = flat_grid(P)
    print(f"[SPARK] grid: {len(lfgws)} freq x {len(cosgt)} sky = {len(CG)} pts "
          f"(span {T/3.15576e7:.2f} yr)", flush=True)

    L0 = np.asarray(P["L0"], float)
    zero = np.zeros(P["npsr"])

    # ---- g0a: SPARK's re-wired EarthDelay F-stat == TE.make_fstat verbatim ----
    print("\n-- g0a: SPARK make_fstat_earth == TE.make_fstat (the incumbent, verbatim) --",
          flush=True)
    te_f = TE.make_fstat(P, earth1)
    sub = np.arange(0, len(CG), max(1, len(CG) // 512))[:512]
    t1 = time.time()
    s_te, _ = scan_grid(lambda *a: te_f(*a[2:]), (None, None), CG[sub], GP[sub], LF[sub])
    s_me, _ = scan_grid(fe, (jnp.asarray(L0), data), CG[sub], GP[sub], LF[sub])
    d_a = float(np.max(np.abs(s_te - s_me)))
    print(f"   {len(sub)} pts, max|SPARK - TE| = {d_a:.3e}  ({time.time()-t1:.1f}s)", flush=True)
    g0a = d_a < GATE_CONT
    print(f"   g0a {'PASS' if g0a else '*** FAIL'} (bar {GATE_CONT:.0e})", flush=True)

    # ---- g0b: THE USER'S GATE -- coherent at ZERO certified == EarthDelay ----
    print("\n-- g0b: coherent detector at pmask=0 == EarthDelay F-stat (THE gate) --",
          flush=True)
    ck = f"{OUT}/spark_g0_grids.npz"
    if os.path.exists(ck):
        # skip-on-exist resume (proven on the production SLURM path, EMBER 2.4)
        z = np.load(ck, allow_pickle=True)
        stat_e = z["stat_earth"]; free_e = z["free_earth"]
        stat_c0 = z["stat_coh_m0"]; free_c0 = z["free_coh_m0"]
        assert len(stat_e) == len(CG), (
            f"STOP: checkpoint grid {len(stat_e)} != this grid {len(CG)}; refusing to "
            f"score a stale checkpoint.")
        print(f"   [resume] both grids from {ck} ({len(stat_e)} pts)", flush=True)
    else:
        t1 = time.time()
        stat_e, free_e = scan_grid(fe, (jnp.asarray(L0), data), CG, GP, LF)
        print(f"   earth full grid {time.time()-t1:.1f}s", flush=True)
        t1 = time.time()
        stat_c0, free_c0 = scan_grid(fc, (jnp.asarray(L0), jnp.asarray(zero), data),
                                     CG, GP, LF)
        print(f"   coherent(pmask=0) full grid {time.time()-t1:.1f}s", flush=True)
    d_b = float(np.max(np.abs(stat_e - stat_c0)))
    d_bf = float(np.max(np.abs(free_e - free_c0)))
    print(f"   max|lnL_coh(m=0) - lnL_earth| = {d_b:.3e}", flush=True)
    print(f"   max|free_coh(m=0)  - free_earth| = {d_bf:.3e}", flush=True)
    g0b = d_b < GATE_CONT
    print(f"   g0b {'PASS' if g0b else '*** FAIL'} (bar {GATE_CONT:.0e}, EMBER 2.2b)",
          flush=True)

    # ---- CHECKPOINT the grids BEFORE the fragile readback step. The two full-grid scans are
    # ~30 min of GPU; a late exception must not discard them (EMBER 2.5: never lose banked
    # work to a step that ran after it).
    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark_g0_grids.npz", stat_earth=stat_e, stat_coh_m0=stat_c0,
             free_earth=free_e, free_coh_m0=free_c0, CG=CG, GP=GP, LF=LF,
             lfgws=lfgws, cosgt=cosgt, gwphi=gwphi, span_s=T,
             g0a_maxdev=d_a, g0b_maxdev=d_b, g0b_free_maxdev=d_bf)
    print(f"   [checkpoint] grids -> {OUT}/spark_g0_grids.npz", flush=True)

    # ---- g0c: the data-driven loud list from the coherent detector's own m=0 map ----
    print("\n-- g0c: loud-source list from the coherent detector (TRUTH-BLIND selection) --",
          flush=True)
    l0e = ll0_earth(earth1, L0, data, float(np.mean(lfgws)))
    l0c = ll0_coh(coh1, L0, zero, data, float(np.mean(lfgws)))
    print(f"   ll0: earth {l0e:.6f}  coh(m=0) {l0c:.6f}  |d| {abs(l0e-l0c):.3e}", flush=True)
    nf, ns = len(lfgws), len(cosgt)
    twoF = (2.0 * (stat_c0 - l0c)).reshape(nf, ns)
    scan = dict(twoF=twoF, free=free_c0.reshape(nf, ns, 4), lfgws=lfgws,
                cosgt=cosgt, gwphi=gwphi, ll0=l0c)
    import trackB_step2_seeder as S2S
    seeds = S2S.pick_seeds_sky(scan, twoF_thresh=25.0, excl_sky_deg=25.0)
    print(f"   2F: max {twoF.max():.1f}  median {np.median(twoF):.1f}  "
          f"seeds {len(seeds)}", flush=True)
    ncap = S2S.report_capture(P, seeds, "coh m=0", cap_deg=25.0)
    g0c = (ncap == 3)
    print(f"   g0c {'PASS' if g0c else '*** FAIL'} ({ncap}/3 loud captured)", flush=True)

    ok = bool(g0a and g0b and g0c)
    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark_g0.npz",
             g0a_maxdev=d_a, g0b_maxdev=d_b, g0b_free_maxdev=d_bf,
             g0c_ncap=ncap, gate_bar=GATE_CONT, passed=ok,
             stat_earth=stat_e, stat_coh_m0=stat_c0, twoF_coh_m0=twoF,
             ll0_earth=l0e, ll0_coh=l0c, lfgws=lfgws, cosgt=cosgt, gwphi=gwphi,
             npt=len(CG), span_s=T, T_label=T_LABEL, h=H_CELL, tier=TIER,
             seed_twoF=np.array([s["twoF"] for s in seeds]) if seeds else np.zeros(0),
             note=("g0a: SPARK's re-wired EarthDelay F-stat == TE.make_fstat verbatim. "
                   "g0b: coherent detector at pmask=0 == EarthDelay F-stat (EMBER 2.2b bar: "
                   "continuous <1e-6; NOT bit-identity -- different builders/graphs). "
                   "g0c: 3/3 loud captured from the coherent map, selection TRUTH-BLIND."))
    print(f"\n[SPARK] g0 {'PASSED' if ok else '*** FAILED'} -> {OUT}/spark_g0.npz", flush=True)
    return 0 if ok else 1


# ============================================================
# MODE s2 -- the arithmetic
# ============================================================
def mode_s2(n_null=100):
    gp = f"{OUT}/spark_g0.npz"
    if not os.path.exists(gp):
        print("*** STOP: g0 has not run. The detector is ungated; s2 refuses.", flush=True)
        return 2
    if not bool(np.load(gp, allow_pickle=True)["passed"]):
        print("*** STOP: g0 FAILED. s2 refuses to score an ungated detector.", flush=True)
        return 2

    print("=" * 78, flush=True)
    print("SPARK S2 -- the gain-vs-cost ledger", flush=True)
    print("=" * 78, flush=True)
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(T_LABEL, verbose=True)
    earth1, coh1 = build_detectors(P)
    fc = make_fstat_coh(coh1)

    NOISE0, DIST0 = 7560000, 17560000
    R = draw_realisation(P, NOISE0, DIST0, faint_present=True)
    tg = faint_targets(P, R["theta_true"])
    CGf = np.array([t["cosgt"] for t in tg])
    GPf = np.array([t["gwphi"] for t in tg])
    LFf = np.array([t["lfgw"] for t in tg])
    L0 = np.asarray(P["L0"], float)
    Lt = np.asarray(R["L_true"], float)
    zero = np.zeros(P["npsr"]); one = np.ones(P["npsr"])
    lref = float(np.mean(TE.freq_grid(float(ds.getspan(P["disco_psrs"])))))

    print(f"[SPARK] reservoir: {len(tg)} faint sources, log10_h "
          f"{min(t['log10_h'] for t in tg):.2f}..{max(t['log10_h'] for t in tg):.2f}", flush=True)

    def twoF_at(Ld, pmask, data):
        s, _ = scan_grid(fc, (jnp.asarray(Ld), jnp.asarray(pmask), data), CGf, GPf, LFf,
                         chunk=16)
        return 2.0 * (s - ll0_coh(coh1, Ld, pmask, data, lref))

    # ---------- GAIN ----------
    states = {}
    states["s0"] = dict(Ld=L0, pmask=zero, label="pmask=0 (incumbent seeder, iteration 0)")
    cert_info = {}
    for g in IGNITERS:
        cs = certified_sets(g)
        cert_info[g] = cs
        # the MODAL banked realisation's certified set (re-derived from raw)
        ncs = np.array([r["n_cert"] for r in cs["rows"]])
        j = int(np.argmax(ncs))          # the MOST-certified realisation: generous
        r = cs["rows"][j]
        pm = np.where(r["cert"], np.asarray(r["qmax"], float), 0.0)
        Ld = np.where(r["cert"], Lt, L0)
        states[f"sC_{g}"] = dict(Ld=Ld, pmask=pm, label=(
            f"{g}: q-weighted on the banked certified set, N_cert={r['n_cert']} "
            f"(max over {len(ncs)} reals; mean {ncs.mean():.2f}, banked corr "
            f"{cs['corr_banked']:.2f}); certified distances at L_true"))
        print(f"[SPARK] {g}: N_cert per real min/mean/max = "
              f"{ncs.min()}/{ncs.mean():.2f}/{ncs.max()}  banked corr={cs['corr_banked']:.2f} "
              f"wrong={cs['wrong_banked']:.3f}  floor={cs['floor']:.3f}", flush=True)
    states["sMAX"] = dict(Ld=Lt, pmask=one, label=(
        "THE CEILING: all 116 pulsars coherent at L_true (perfect array-wide certification; "
        "unreachable by construction)"))

    rows = {}
    for k, st in states.items():
        t0 = time.time()
        rows[k] = twoF_at(st["Ld"], st["pmask"], R["data"])
        print(f"[SPARK] 2F[{k:12s}] median {np.median(rows[k]):7.3f}  max {rows[k].max():7.3f}"
              f"   ({time.time()-t0:.1f}s)  {st['label']}", flush=True)

    # ---------- COST: the null-calibrated detection floor ----------
    print(f"\n[SPARK] nulls: {n_null} draws, faint reservoir ABSENT from the data "
          f"(3 loud + noise retained) -> 2F at the 13 reservoir positions", flush=True)
    null0, nullmax = [], []
    for i in range(n_null):
        Rn = draw_realisation(P, NOISE0 + 1000 + i, DIST0 + 1000 + i, faint_present=False)
        null0.append(twoF_at(L0, zero, Rn["data"]))
        nullmax.append(twoF_at(np.asarray(Rn["L_true"], float), one, Rn["data"]))
        if (i + 1) % 20 == 0:
            print(f"   null {i+1}/{n_null}", flush=True)
    null0 = np.concatenate(null0); nullmax = np.concatenate(nullmax)
    f0 = adopt(np.maximum(null0, 0.0))
    fmax = adopt(np.maximum(nullmax, 0.0))
    print(f"[SPARK] floor(s0)   {f0['floor']:.3f} +- {f0['err']:.3f} [{f0['estimator']}] "
          f"zf={f0['zf']:.3f} n={f0['n']}", flush=True)
    print(f"[SPARK] floor(sMAX) {fmax['floor']:.3f} +- {fmax['err']:.3f} "
          f"[{fmax['estimator']}] zf={fmax['zf']:.3f} n={fmax['n']}", flush=True)

    # ---------- COST: lnK growth per added source ----------
    print(f"\n[SPARK] lnK growth per added source (CPU, exact)", flush=True)
    lk = lnK_vs_nsrc(P, R["theta_cell"])
    for r in (lk[0], lk[2], lk[3], lk[-1]):
        print(f"   n_src={r['n_src']:2d}  K_sum={r['K_sum']:12.0f}  "
              f"lnK_med={r['lnK_med']:6.3f}  dL_med={r['dL_med']:.3e}", flush=True)

    # ---------- THE LEDGER ----------
    print("\n" + "=" * 78, flush=True)
    print("THE LEDGER -- gain vs cost, per igniter", flush=True)
    print("=" * 78, flush=True)
    bar0 = f0["floor"] + f0["err"]
    barM = fmax["floor"] + fmax["err"]
    print(f"{'state':14s} {'2F med':>9s} {'2F max':>9s} {'gain med':>9s} "
          f"{'bar':>8s} {'clears?':>8s} {'shortfall x':>12s}", flush=True)
    ledger = {}
    for k in states:
        v = rows[k]
        gain = float(np.median(v - rows["s0"]))
        bar = barM if k == "sMAX" else bar0
        clears = int((v > bar).sum())
        short = float(bar / max(np.max(v), 1e-9))
        ledger[k] = dict(med=float(np.median(v)), mx=float(np.max(v)), gain=gain,
                         bar=float(bar), clears=clears, shortfall=short)
        print(f"{k:14s} {np.median(v):9.3f} {v.max():9.3f} {gain:9.3f} "
              f"{bar:8.3f} {clears:8d} {short:12.2f}", flush=True)

    alive = any(ledger[k]["clears"] > 0 for k in ledger)
    verdict = "CASCADE-ALIVE" if alive else "CASCADE-DEAD"
    print(f"\nVERDICT: {verdict}", flush=True)
    if not alive:
        print(f"  Even sMAX (all 116 coherent at L_true, oracle-anchored, no trials factor) "
              f"clears 0/{len(tg)}.", flush=True)
        print(f"  Shortfall at the ceiling: {ledger['sMAX']['shortfall']:.2f}x in 2F "
              f"(need 2F>{barM:.2f}, best faint source reaches {ledger['sMAX']['mx']:.3f}).",
              flush=True)

    os.makedirs(OUT, exist_ok=True)
    bank = dict(verdict=verdict, alive=alive,
                faint_idx=np.array([t["idx"] for t in tg]),
                faint_cosgt=CGf, faint_gwphi=GPf, faint_lfgw=LFf,
                faint_h=np.array([t["log10_h"] for t in tg]),
                floor_s0=f0["floor"], floor_s0_err=f0["err"], floor_s0_zf=f0["zf"],
                floor_s0_est=f0["estimator"], floor_s0_n=f0["n"],
                floor_s0_gumbel=f0["gumbel"], floor_s0_emp=f0["emp"],
                floor_max=fmax["floor"], floor_max_err=fmax["err"], floor_max_zf=fmax["zf"],
                floor_max_est=fmax["estimator"], floor_max_n=fmax["n"],
                null_twoF_s0=null0, null_twoF_max=nullmax,
                lnK_n_src=np.array([r["n_src"] for r in lk]),
                lnK_K_sum=np.array([r["K_sum"] for r in lk]),
                lnK_med=np.array([r["lnK_med"] for r in lk]),
                lnK_dL_med=np.array([r["dL_med"] for r in lk]),
                h=H_CELL, tier=TIER, T_label=T_LABEL, n_null=n_null,
                noise_seed=NOISE0, dist_seed=DIST0,
                note=("RAW 2F per (state, faint source) banked -- NOT verdicts. "
                      "Oracle-anchored: 2F read at each faint source's TRUE sky/freq. "
                      "Floor: no trials factor (no search). sMAX = all 116 coherent at "
                      "L_true. All four conventions FAVOUR the cascade -> a DEAD verdict "
                      "is a fortiori."))
    for k in states:
        bank[f"twoF_{k}"] = rows[k]
        bank[f"pmask_{k}"] = np.asarray(states[k]["pmask"])
        bank[f"label_{k}"] = states[k]["label"]
    for g in IGNITERS:
        cs = cert_info[g]
        bank[f"ncert_{g}"] = np.array([r["n_cert"] for r in cs["rows"]])
        bank[f"ncert_true_{g}"] = np.array([r["n_cert_true"] for r in cs["rows"]])
        bank[f"ncert_wrong_{g}"] = np.array([r["n_cert_wrong"] for r in cs["rows"]])
        bank[f"floor_{g}"] = cs["floor"]
        bank[f"corr_banked_{g}"] = cs["corr_banked"]
    np.savez(f"{OUT}/spark_s2.npz", **bank)
    print(f"[SPARK] banked -> {OUT}/spark_s2.npz", flush=True)
    return 0


# ============================================================
# MODE s2b -- THE ORACLE LEDGER (the corrected arithmetic)
# ============================================================
# WHY s2 (the F-stat ledger, banked as spark_s2_fstat_trail.npz) DOES NOT ANSWER THE QUESTION.
# It ran and is kept as a trail, but two defects make its verdict uninterpretable, BOTH found
# by reading its own output rather than by any gate:
#
#  (1) LOUD-SIDELOBE LEAKAGE. The ncw=1 detector models ONE source against data holding 3 loud
#      + 13 faint. At a faint position the statistic is dominated by the LOUD sources' sidelobes
#      (this is exactly the pathology TE's sky-exclusion NMS exists to defeat --
#      trackB_estimator.py:412-415, "loud2's true peak sank to rank-13"). Measured: the null
#      floor came out at 145.3 (a chi2_4-like statistic should sit near 9.5) because the null
#      RETAINS the loud sources. The single faint source that "cleared" is a faint source lying
#      near a loud one; the detector was seeing the loud sidelobe, not the faint source. s2's
#      "CASCADE-ALIVE" was an artifact of that leakage and is RETRACTED here, not reported.
#      In the cascade the detected sources are FITTED and subtracted before the next [D]; the
#      statistic must model them.
#
#  (2) mc MIS-SPECIFICATION. TE.make_fstat hard-codes SEED_MC = 9.0 (a fixed reference chirp
#      mass, TE:163). For the EARTH term that is nearly harmless -- over T the frequency barely
#      evolves. For the PULSAR term it is fatal: the pterm looks back ~L/c (kyr) into the
#      binary's past, where the frequency was materially different, so its phase depends
#      strongly on mc. Measured consequence in s2: coherence at SEED_MC LOWERED the statistic
#      (sMAX median 15.9 vs s0's 21.1) -- the coherent model actively mismatches the data.
#      That is a real and reportable fact about a coherent detector that does not know mc, but
#      it is ANTI-cascade, and SPARK's a-fortiori structure requires every convention to FAVOUR
#      the cascade. A ceiling must be measured at the TRUE mc.
#
# s2b fixes both, and in doing so makes the statistic STRICTLY MORE GENEROUS than any real
# detector could be:
#   * LOUD SUBTRACTED, and the other 12 faint subtracted too: the model holds all 16 sources at
#     TRUTH and toggles ONLY source k. No confusion of any kind survives.
#   * ORACLE at truth on EVERY axis, mc included -- no profiling, no search, no SEED_MC.
#   * dlnL_k(pmask) = logL(all 16 at truth) - logL(same, source k removed), the Neyman-Pearson
#     statistic when everything is known. Nothing can beat it.
#   * The floor is null-calibrated at MATCHED state and carries the zero-fraction column.
# If THIS does not clear the floor, no detector can, and the cascade is arithmetically closed.
def mode_s2b(n_null=100):
    gp = f"{OUT}/spark_g0.npz"
    if not os.path.exists(gp) or not bool(np.load(gp, allow_pickle=True)["passed"]):
        print("*** STOP: g0 absent or failed; s2b refuses to score an ungated detector.",
              flush=True)
        return 2
    print("=" * 78, flush=True)
    print("SPARK S2b -- THE ORACLE LEDGER (loud subtracted, true mc, all-coherent ceiling)",
          flush=True)
    print("=" * 78, flush=True)
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(T_LABEL, verbose=True)
    lb = P["amo"]["logL_batch_theta"]          # jax.vmap(logL, in_axes=(0, None, None))
    nd, AI = P["n_dist"], P["AI"]
    faint = list(range(N_LOUD, P["ncw"]))

    NOISE0, DIST0 = 7560000, 17560000
    L0 = np.asarray(P["L0"], float)
    zero = np.zeros(P["npsr"]); one = np.ones(P["npsr"])

    def dlnL_all(R, pmask, Ld):
        """dlnL_k = logL(all 16 at truth) - logL(same, source k REMOVED), for all 13 faint k.
        ONE batched call: theta stack is (1 + 13, nth)."""
        th_on = R["theta_true"].copy(); th_on[AI] = np.asarray(Ld, float)
        stack = [th_on]
        for k in faint:
            t = th_on.copy(); t[nd + 8 * k + I_H] = H_ABSENT
            stack.append(t)
        lls = np.asarray(lb(jnp.asarray(np.stack(stack)), R["data"], jnp.asarray(pmask)))
        return lls[0] - lls[1:]

    # ---------- states ----------
    R = draw_realisation(P, NOISE0, DIST0, faint_present=True)
    Lt = np.asarray(R["L_true"], float)
    states = {"s0": dict(Ld=L0, pmask=zero,
                         label="pmask=0 (fully decohered; the incumbent seeder's limit)")}
    cert_info = {}
    for g in IGNITERS:
        cs = certified_sets(g); cert_info[g] = cs
        ncs = np.array([r["n_cert"] for r in cs["rows"]])
        r = cs["rows"][int(np.argmax(ncs))]          # the MOST-certified realisation: generous
        states[f"sC_{g}"] = dict(
            Ld=np.where(r["cert"], Lt, L0),
            pmask=np.where(r["cert"], np.asarray(r["qmax"], float), 0.0),
            label=(f"{g}: q-weighted on the banked certified set, N_cert={r['n_cert']} "
                   f"(MAX over {len(ncs)} reals; mean {ncs.mean():.2f} = banked corr "
                   f"{cs['corr_banked']:.2f} + wrong {cs['wrong_banked']:.3f})"))
        # READBACK GATE: re-derived n_cert must equal banked corr + wrong, exactly.
        dev = abs(ncs.mean() - (cs["corr_banked"] + cs["wrong_banked"]))
        print(f"[SPARK] {g}: N_cert min/mean/max={ncs.min()}/{ncs.mean():.3f}/{ncs.max()}  "
              f"banked corr+wrong={cs['corr_banked']+cs['wrong_banked']:.3f}  "
              f"|dev|={dev:.3e}  {'OK' if dev < 1e-9 else '*** READBACK MISMATCH'}", flush=True)
    states["sMAX"] = dict(Ld=Lt, pmask=one, label=(
        "THE CEILING: all 116 pulsars coherent at L_true (perfect array-wide certification)"))

    rows = {k: dlnL_all(R, st["pmask"], st["Ld"]) for k, st in states.items()}
    for k in states:
        v = rows[k]
        print(f"[SPARK] dlnL[{k:12s}] median {np.median(v):9.4f}  max {v.max():9.4f}  "
              f"{states[k]['label']}", flush=True)

    # ---------- floors: matched nulls (source k ABSENT from the DATA) ----------
    print(f"\n[SPARK] nulls: {n_null} draws, the 13 faint ABSENT from the DATA; the model "
          f"still toggles each faint slot -> the statistic's own null", flush=True)
    nulls = {k: [] for k in states}
    for i in range(n_null):
        Rn = draw_realisation(P, NOISE0 + 1000 + i, DIST0 + 1000 + i, faint_present=False)
        Ltn = np.asarray(Rn["L_true"], float)
        for k, st in states.items():
            Ld = Ltn if k == "sMAX" else (np.where(st["pmask"] > 0, Ltn, L0)
                                          if k != "s0" else L0)
            nulls[k].append(dlnL_all(Rn, st["pmask"], Ld))
        if (i + 1) % 25 == 0:
            print(f"   null {i+1}/{n_null}", flush=True)
    floors = {}
    for k in states:
        x = np.concatenate(nulls[k])
        floors[k] = adopt(np.maximum(x, 0.0))
        f = floors[k]
        print(f"[SPARK] floor[{k:12s}] {f['floor']:8.4f} +- {f['err']:.4f} "
              f"[{f['estimator']}] zf={f['zf']:.3f} n={f['n']}", flush=True)

    # ---------- lnK growth ----------
    print(f"\n[SPARK] lnK growth per added source (CPU, exact)", flush=True)
    lk = lnK_vs_nsrc(P, R["theta_cell"])
    for r in (lk[0], lk[2], lk[3], lk[-1]):
        print(f"   n_src={r['n_src']:2d}  K_sum={r['K_sum']:12.0f}  lnK_med={r['lnK_med']:6.3f}"
              f"  dL_med={r['dL_med']:.3e}", flush=True)
    dlnK = lk[-1]["lnK_med"] - lk[2]["lnK_med"]
    print(f"   lnK cost of growing the list 3 -> 16 sources: +{dlnK:.3f} nat (median pulsar)",
          flush=True)

    # ---------- THE LEDGER ----------
    print("\n" + "=" * 78, flush=True)
    print("THE ORACLE LEDGER -- gain vs cost", flush=True)
    print("=" * 78, flush=True)
    print(f"{'state':14s} {'dlnL med':>10s} {'dlnL max':>10s} {'gain med':>10s} "
          f"{'bar':>9s} {'clears':>7s} {'shortfall x':>12s}", flush=True)
    ledger = {}
    for k in states:
        v = rows[k]; f = floors[k]
        bar = f["floor"] + f["err"]
        gain = float(np.median(v - rows["s0"]))
        clears = int((v > bar).sum())
        short = float(bar / max(np.max(v), 1e-12))
        ledger[k] = dict(med=float(np.median(v)), mx=float(np.max(v)), gain=gain,
                         bar=float(bar), clears=clears, shortfall=short)
        print(f"{k:14s} {np.median(v):10.4f} {v.max():10.4f} {gain:10.4f} "
              f"{bar:9.4f} {clears:7d} {short:12.2f}", flush=True)

    alive = any(ledger[k]["clears"] > 0 for k in ledger)
    verdict = "CASCADE-ALIVE" if alive else "CASCADE-DEAD"
    print(f"\nVERDICT: {verdict}", flush=True)
    if not alive:
        m = ledger["sMAX"]
        print(f"  The CEILING (all 116 coherent at L_true, oracle at truth on every axis, "
              f"every other source subtracted) clears 0/{len(faint)}.", flush=True)
        print(f"  Shortfall at the ceiling: {m['shortfall']:.2f}x  (bar {m['bar']:.4f}, "
              f"best faint source {m['mx']:.4f}).", flush=True)
        print(f"  Coherence gain at the ceiling: {m['gain']:+.4f} nat (median) -- against a "
              f"lnK cost of +{dlnK:.3f} nat for the list growth it would have to pay for.",
              flush=True)

    bank = dict(verdict=verdict, alive=alive, faint_idx=np.array(faint),
                lnK_n_src=np.array([r["n_src"] for r in lk]),
                lnK_K_sum=np.array([r["K_sum"] for r in lk]),
                lnK_med=np.array([r["lnK_med"] for r in lk]),
                lnK_dL_med=np.array([r["dL_med"] for r in lk]),
                lnK_cost_3_to_16=dlnK,
                h=H_CELL, tier=TIER, T_label=T_LABEL, n_null=n_null,
                noise_seed=NOISE0, dist_seed=DIST0,
                note=("ORACLE ledger. RAW dlnL per (state, faint source) banked -- never "
                      "verdicts. dlnL_k = logL(all 16 at truth) - logL(same, k removed): the "
                      "Neyman-Pearson statistic with EVERY other source subtracted, at truth "
                      "on every axis INCLUDING mc, no search, no profiling. sMAX = all 116 "
                      "coherent at L_true. Floors null-calibrated at matched state (faint "
                      "ABSENT from the data), RECUT adopt(), zero-fraction banked. Every "
                      "convention favours the cascade -> a DEAD verdict is a fortiori. "
                      "Supersedes spark_s2_fstat_trail.npz, whose ncw=1 F-stat was dominated "
                      "by loud-sidelobe leakage (null floor 145.3) and mis-specified at "
                      "SEED_MC=9.0; that run's CASCADE-ALIVE is RETRACTED."))
    for k in states:
        bank[f"dlnL_{k}"] = rows[k]
        bank[f"null_dlnL_{k}"] = np.concatenate(nulls[k])
        bank[f"pmask_{k}"] = np.asarray(states[k]["pmask"])
        bank[f"label_{k}"] = states[k]["label"]
        for kk, vv in floors[k].items():
            bank[f"floor_{k}_{kk}"] = vv
    for g in IGNITERS:
        cs = cert_info[g]
        bank[f"ncert_{g}"] = np.array([r["n_cert"] for r in cs["rows"]])
        bank[f"ncert_true_{g}"] = np.array([r["n_cert_true"] for r in cs["rows"]])
        bank[f"ncert_wrong_{g}"] = np.array([r["n_cert_wrong"] for r in cs["rows"]])
        bank[f"floor_chorus_{g}"] = cs["floor"]
        bank[f"corr_banked_{g}"] = cs["corr_banked"]
    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark_s2b.npz", **bank)
    print(f"[SPARK] banked -> {OUT}/spark_s2b.npz", flush=True)
    return 0


# ============================================================
# MODE s2c -- THE ORACLE F-STAT LEDGER (the statistic that is actually a detector)
# ============================================================
# WHY s2b's STATISTIC WAS WRONG, found by its own null (banked: it crashed the Gumbel fit).
# s2b used dlnL_k = logL(k in the model AT TRUTH) - logL(k removed). Under the NULL (k absent
# from the data) inserting k at its TRUE amplitude can only LOWER the likelihood, so the null
# collapsed to an all-<=0 sample: a 100% zero point mass, no floor, and a scipy OverflowError.
# That is the correct behaviour of an incorrect object. dlnL-at-truth is not a DETECTION
# statistic -- it already knows the true amplitude, which under the null does not exist. A
# detector must PROFILE the unknown source amplitude/phase; only then does the null have the
# chi2-like spread a floor can be cut from.
#
# s2c is s2b with that repaired, and keeps every generous convention:
#   2F_k(pmask) = 2*[ max_{cos_inc, log10_h, phase0, psi} logL(theta_k(free), data, pmask)
#                     - logL(theta with k REMOVED, data, pmask) ]
#   * PROFILED over the 4 unknown extrinsics -> a real detector, with a real null.
#   * ORACLE-ANCHORED on sky, freq, AND log10_mc (mc at TRUTH, not TE's SEED_MC=9.0 -- the
#     mis-specification that made s2's coherent states LOSE signal).
#   * ALL OTHER 15 SOURCES AT TRUTH IN THE MODEL -> the 3 loud AND the other 12 faint are
#     subtracted. No loud-sidelobe leakage (the defect that produced s2's 145.3 floor).
#   * sMAX = all 116 coherent at L_true: the coherence ceiling.
#   * Floor null-calibrated at MATCHED state, RECUT adopt(), zero-fraction banked.
# The same Adam profile as TE.make_fstat (via _adam_profile) -- the incumbent's own optimiser.
def make_profile16(P):
    """Profile source k's 4 extrinsics on the FULL 16-source masked likelihood. k enters as a
    TRACED base offset so the 13 reservoir sources vmap in one call."""
    fl = P["amo"]["logL"]

    def lnl_at(th_base, kb, data, pmask, free):
        th = th_base.at[kb + I_COSINC].set(free[0])
        th = th.at[kb + I_H].set(free[1])
        th = th.at[kb + I_PH0].set(free[2])
        th = th.at[kb + I_PSI].set(free[3])
        return fl(th, data, pmask)

    prof = _adam_profile(lnl_at)
    return jax.jit(jax.vmap(prof, in_axes=(None, 0, None, None)))


def mode_s2c(n_null=100):
    gp = f"{OUT}/spark_g0.npz"
    if not os.path.exists(gp) or not bool(np.load(gp, allow_pickle=True)["passed"]):
        print("*** STOP: g0 absent or failed; s2c refuses to score an ungated detector.",
              flush=True)
        return 2
    print("=" * 78, flush=True)
    print("SPARK S2c -- THE ORACLE F-STAT LEDGER", flush=True)
    print("=" * 78, flush=True)
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(T_LABEL, verbose=True)
    prof = make_profile16(P)
    lb = P["amo"]["logL_batch_theta"]
    nd, AI = P["n_dist"], P["AI"]
    faint = list(range(N_LOUD, P["ncw"]))
    KB = jnp.asarray([nd + 8 * k for k in faint])

    NOISE0, DIST0 = 7560000, 17560000
    L0 = np.asarray(P["L0"], float)
    zero = np.zeros(P["npsr"]); one = np.ones(P["npsr"])

    def twoF_all(R, pmask, Ld):
        """2F_k for all 13 reservoir sources: profiled max minus the k-removed reference."""
        th_on = R["theta_true"].copy(); th_on[AI] = np.asarray(Ld, float)
        off = []
        for k in faint:
            t = th_on.copy(); t[nd + 8 * k + I_H] = H_ABSENT
            off.append(t)
        ll_off = np.asarray(lb(jnp.asarray(np.stack(off)), R["data"], jnp.asarray(pmask)))
        mx, _ = prof(jnp.asarray(th_on), KB, R["data"], jnp.asarray(pmask))
        return 2.0 * (np.asarray(mx) - ll_off)

    R = draw_realisation(P, NOISE0, DIST0, faint_present=True)
    Lt = np.asarray(R["L_true"], float)
    states = {"s0": dict(Ld=L0, pmask=zero,
                         label="pmask=0 (fully decohered; the incumbent seeder's limit)")}
    cert_info = {}
    for g in IGNITERS:
        cs = certified_sets(g); cert_info[g] = cs
        ncs = np.array([r["n_cert"] for r in cs["rows"]])
        r = cs["rows"][int(np.argmax(ncs))]
        states[f"sC_{g}"] = dict(
            Ld=np.where(r["cert"], Lt, L0),
            pmask=np.where(r["cert"], np.asarray(r["qmax"], float), 0.0),
            label=(f"{g}: q-weighted certified set, N_cert={r['n_cert']} (MAX over "
                   f"{len(ncs)} reals; mean {ncs.mean():.3f})"))
        dev = abs(ncs.mean() - (cs["corr_banked"] + cs["wrong_banked"]))
        print(f"[SPARK] {g}: N_cert min/mean/max={ncs.min()}/{ncs.mean():.3f}/{ncs.max()}  "
              f"banked corr+wrong={cs['corr_banked']+cs['wrong_banked']:.3f}  |dev|={dev:.3e}  "
              f"{'OK' if dev < 1e-9 else '*** READBACK MISMATCH'}", flush=True)
    states["sMAX"] = dict(Ld=Lt, pmask=one,
                          label="THE CEILING: all 116 coherent at L_true")
    cert_masks = {k: (states[k]["pmask"] > 0) for k in states}

    rows = {k: twoF_all(R, st["pmask"], st["Ld"]) for k, st in states.items()}
    for k in states:
        v = rows[k]
        print(f"[SPARK] 2F[{k:12s}] median {np.median(v):9.4f}  max {v.max():9.4f}  "
              f"{states[k]['label']}", flush=True)

    print(f"\n[SPARK] nulls: {n_null} draws, the 13 faint ABSENT from the DATA "
          f"(3 loud retained IN THE MODEL and subtracted) -> the statistic's own null",
          flush=True)
    nulls = {k: [] for k in states}
    for i in range(n_null):
        Rn = draw_realisation(P, NOISE0 + 1000 + i, DIST0 + 1000 + i, faint_present=False)
        Ltn = np.asarray(Rn["L_true"], float)
        for k, st in states.items():
            Ld = L0 if k == "s0" else np.where(cert_masks[k], Ltn, L0)
            nulls[k].append(twoF_all(Rn, st["pmask"], Ld))
        if (i + 1) % 25 == 0:
            print(f"   null {i+1}/{n_null}", flush=True)
    floors = {}
    for k in states:
        x = np.concatenate(nulls[k])
        floors[k] = adopt(np.maximum(x, 0.0))
        f = floors[k]
        print(f"[SPARK] floor[{k:12s}] {f['floor']:8.4f} +- {f['err']:.4f} "
              f"[{f['estimator']}] zf={f['zf']:.3f} n={f['n']}", flush=True)

    print(f"\n[SPARK] lnK growth per added source (CPU, exact)", flush=True)
    lk = lnK_vs_nsrc(P, R["theta_cell"])
    for r in (lk[0], lk[2], lk[3], lk[-1]):
        print(f"   n_src={r['n_src']:2d}  K_sum={r['K_sum']:12.0f}  lnK_med={r['lnK_med']:6.3f}"
              f"  dL_med={r['dL_med']:.3e}", flush=True)
    dlnK = lk[-1]["lnK_med"] - lk[2]["lnK_med"]
    print(f"   lnK cost, list 3 -> 16 sources: +{dlnK:.3f} nat (median pulsar)", flush=True)

    print("\n" + "=" * 78, flush=True)
    print("THE ORACLE F-STAT LEDGER -- gain vs cost", flush=True)
    print("=" * 78, flush=True)
    print(f"{'state':14s} {'2F med':>9s} {'2F max':>9s} {'gain med':>9s} {'bar':>9s} "
          f"{'clears':>7s} {'shortfall x':>12s}", flush=True)
    ledger = {}
    for k in states:
        v = rows[k]; f = floors[k]
        bar = f["floor"] + f["err"]
        gain = float(np.median(v - rows["s0"]))
        clears = int((v > bar).sum())
        short = float(bar / max(np.max(v), 1e-12))
        ledger[k] = dict(med=float(np.median(v)), mx=float(np.max(v)), gain=gain,
                         bar=float(bar), clears=clears, shortfall=short)
        print(f"{k:14s} {np.median(v):9.4f} {v.max():9.4f} {gain:9.4f} {bar:9.4f} "
              f"{clears:7d} {short:12.2f}", flush=True)

    alive = any(ledger[k]["clears"] > 0 for k in ledger)
    verdict = "CASCADE-ALIVE" if alive else "CASCADE-DEAD"
    print(f"\nVERDICT: {verdict}", flush=True)
    m = ledger["sMAX"]
    if not alive:
        print(f"  The CEILING -- all 116 coherent at L_true, oracle on sky/freq/mc, every "
              f"other source subtracted, no search -- clears 0/{len(faint)}.", flush=True)
        print(f"  Shortfall at the ceiling: {m['shortfall']:.2f}x (bar {m['bar']:.4f}, best "
              f"reservoir source {m['mx']:.4f}).", flush=True)
    print(f"  Coherence gain at the CEILING: {m['gain']:+.4f} (median 2F), "
          f"i.e. {m['med']/max(ledger['s0']['med'],1e-9):.2f}x the decohered statistic.",
          flush=True)
    for g in IGNITERS:
        print(f"  Gain at the ACHIEVABLE certified set ({g}): "
              f"{ledger['sC_'+g]['gain']:+.4f} median 2F", flush=True)

    bank = dict(verdict=verdict, alive=alive, faint_idx=np.array(faint),
                lnK_n_src=np.array([r["n_src"] for r in lk]),
                lnK_K_sum=np.array([r["K_sum"] for r in lk]),
                lnK_med=np.array([r["lnK_med"] for r in lk]),
                lnK_dL_med=np.array([r["dL_med"] for r in lk]),
                lnK_cost_3_to_16=dlnK, h=H_CELL, tier=TIER, T_label=T_LABEL,
                n_null=n_null, noise_seed=NOISE0, dist_seed=DIST0,
                note=("ORACLE F-STAT ledger -- SPARK's headline. RAW 2F per (state, faint "
                      "source) banked, never verdicts. 2F_k = 2*[max over (cos_inc,log10_h,"
                      "phase0,psi) logL - logL(k removed)], sky/freq/mc ORACLE-ANCHORED at "
                      "truth, all other 15 sources at truth IN THE MODEL (loud subtracted), "
                      "same Adam profile as TE.make_fstat. sMAX = all 116 coherent at L_true. "
                      "Floors null-calibrated at matched state (faint ABSENT from data), "
                      "RECUT adopt(), zero-fraction banked. Every convention favours the "
                      "cascade -> a DEAD verdict is a fortiori. SUPERSEDES spark_s2b (its "
                      "dlnL-at-truth was not a detector: null had a 100% zero point mass) "
                      "and spark_s2_fstat_trail (loud-sidelobe leakage + SEED_MC "
                      "mis-specification; its CASCADE-ALIVE is RETRACTED)."))
    for k in states:
        bank[f"twoF_{k}"] = rows[k]
        bank[f"null_twoF_{k}"] = np.concatenate(nulls[k])
        bank[f"pmask_{k}"] = np.asarray(states[k]["pmask"])
        bank[f"label_{k}"] = states[k]["label"]
        for kk, vv in floors[k].items():
            bank[f"floor_{k}_{kk}"] = vv
    for g in IGNITERS:
        cs = cert_info[g]
        bank[f"ncert_{g}"] = np.array([r["n_cert"] for r in cs["rows"]])
        bank[f"ncert_true_{g}"] = np.array([r["n_cert_true"] for r in cs["rows"]])
        bank[f"ncert_wrong_{g}"] = np.array([r["n_cert_wrong"] for r in cs["rows"]])
        bank[f"floor_chorus_{g}"] = cs["floor"]
        bank[f"corr_banked_{g}"] = cs["corr_banked"]
    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark_s2c.npz", **bank)
    print(f"[SPARK] banked -> {OUT}/spark_s2c.npz", flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["g0", "s2", "s2b", "s2c"])
    ap.add_argument("--n-null", type=int, default=100)
    a = ap.parse_args()
    FN = {"g0": lambda: mode_g0(), "s2": lambda: mode_s2(a.n_null),
          "s2b": lambda: mode_s2b(a.n_null), "s2c": lambda: mode_s2c(a.n_null)}
    sys.exit(FN[a.mode]())
