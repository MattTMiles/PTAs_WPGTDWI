#!/usr/bin/env python
"""SPARK-3 -- ARROW 2, THE HONEST MIDDLE: deconfusion -> certification, at a venue where it CAN
move, with the instrument fixed.

WHAT SPARK-2 LEFT (reports/SPARK2_second_arrow.md). Its endpoints STAND and are the frame:
the faint reservoir UNMODELLED erases certification (offender 0.000) while TRUTH-modelled gives
4.435; the trials cost of modelling the reservoir is +0.578 nat and SATURATES. What it did NOT
deliver is the middle -- the soft-modelled ladder at reachable rungs -- and its 5 fix items are
this file's work plan, in its order.

FIX 0 -- THE VENUE (mode `venue`; SPARK-2 5(5)). SPARK-2 measured N_cert = 0 in EVERY state at
  its cell (h = -13.25, T = 30, lit, 3+13). NOTHING CAN CROSS A BAR NOBODY REACHES: with the
  count identically zero, arrow 2's ledger can only return zeros, and no rung difference is
  measurable. Re-cut from the banked SURFACE raw statistics (recut_surface.py verbatim), that
  cell reads corr = 0.567/real, verdict `below` -- SPARK-2 ran arrow 2 at a BELOW-ONSET cell.
  The venue moves to STORY App A's RESERVED above-onset cells (SURFACE Pair B, the loop cells
  KINDLE reserved -- this is their purpose).

FIX 1 -- THE PER-TARGET E-STEP (mode `gate1`; SPARK-2 1, the structural fix). `B1Split.estep`
  sweeps EVERY pulsar's distance under ONE GLOBAL pmask, so a pulsar at m_p = 0 has its OWN
  pulsar term off -> its distance is inert -> its fringe row is FLAT -> it cannot certify. The
  certified-coherent E-step arrow 2 needs -- TARGET's own term always live, the certified set
  coherent at q, uncertified off -- is one call PER TARGET. Built here (~116x, batched).

FIX 2 -- THE DOUBLE-PERTURBATION PATHOLOGY (mode `gate2`; SPARK-2 4b). Rungs 2 and 5 died
  all-NaN even after phi_bounds clipping. Diagnosed here rather than clipped further. EMBER's
  `mc_scan` measured "non-evaluable 0% across the entire mc box" for the LOUD sources; SPARK-2
  is the counter-example at FAINT. Faint soft draws are treated as UNTESTED MACHINERY and gated
  BEFORE the grid.

FIX 3 -- THE FISHER, TWO-SIDED (SPARK-2 4a). The monolithic joint is MEASURED unaffordable (two
  builds failed to return on an A100-80GB in 1 h). Pre-registered replacement, both banked as
  columns on every unit:
      OPTIMISTIC = the CONDITIONAL width  (sigma_cond = 1/sqrt(F_ii) <= sigma_marg ALWAYS:
                   it OVERSTATES how well the reservoir is known -> favours arrow 2)
      PESSIMISTIC = the PRIOR width       (no tightening at all: the loop learns nothing)
  The honest answer lives BETWEEN. If both bounds give the SAME verdict, the verdict stands
  without the joint. If they straddle -> STRADDLED, and the budgeted chunked-JVP cost to split
  them is specced with its wallclock cap stated up front.

CONVENTIONS. cpus-per-task=8 (the noise-draw thread hazard: NoiseDrawer's eigh basis rotates
with BLAS threads). Offender + floor: RECUT verbatim (recut_surface.py:75-78/:101 adopt());
zero-fraction is a REQUIRED column beside every floor. Lean-npz: raw per-pulsar dlnL/lnK/qmax
and raw null offender vectors banked per unit, NEVER verdicts. In-process rung looping (one
build, all rungs -- SPARK-2 4c's cache lesson). Machinery REUSED not reimplemented: ignite.py's
problem build, surface.py's structure axis, spark.py's detector + floors.

SCOPE. MOCK spine (AXIS 1440 MHz, 10.15(a)) -- no real TOA is touched; the residuals ARE the
injected CW+CURN. A verdict here is a statement about an approximately-real IPTA's geometry and
noise budget.

Modes:
  venue          FIX 0: does N_cert > 0 EXIST at the reserved cells? (CPU, banks only)
  gate1          FIX 1: the per-target E-step + its gate chain
  gate2          FIX 2: the faint soft-draw anatomy (non-evaluable fraction per rung)
  unit --real R  one realisation, ALL rungs in-process: (a) + (b) + the truth control
  ledger         (c): collate -> the ledger + the pre-registered verdict
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")
import sys, time, argparse, glob
import numpy as np

HSYMT = os.environ.get("HSYMT_ROOT", "/home/mattm/projects/HSYMT")
for _p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
           "hpc_harbor/surface", "hpc_harbor/spark"):
    sys.path.insert(0, f"{HSYMT}/{_p}")

# ============================================================
# THE VENUE -- STORY App A's RESERVED above-onset cells (SURFACE Pair B)
# ============================================================
# KINDLE reserved these as "the first genuinely above-onset, low-impurity seed set the
# programme has ever had, and both untouched by the floor fix". FIX 0 confirms N_cert > 0
# EXISTS here before a single GPU-second is spent on the grid.
VENUES = {
    "A": dict(h=-12.75, T=40, tier="vlbi", k=5, story_corr=4.07),
    "B": dict(h=-13.00, T=40, tier="vlbi", k=5, story_corr=3.57),
}
SPARK2_CELL = dict(h=-13.25, T=30, tier="lit", k=3)   # the control: where N_cert was 0

OUT = f"{HSYMT}/SPARK3_results"
QBAR = 0.9
RUNGS = [0, 2, 5, 8]
N_SOFT = 8
SOFT_SEED = 20260717
H_FAINT = -14.25
H_ABSENT = -18.0

# the loop-actionable axes per faint source; sky is frozen (ignite.py:102 LOOP_PARS -- the loop
# cannot act on it), so tightening is measured on what the loop can actually use.
I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)
FAINT_PARS = (I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI)
# the prior box per axis -- also the PESSIMISTIC (no-tightening) Fisher bound of FIX 3
BOX = {I_COSINC: 1.0, I_MC: 0.87, I_FGW: 0.5, I_H: 1.0, I_PH0: 3.14, I_PSI: 1.57}


def hkey(h):
    return f"{abs(h):.2f}".replace(".", "")


def cell_tag(h, T, tier, k):
    return f"h{hkey(h)}_T{T}_{tier}_k{k}"


# ============================================================
# FIX 0 -- THE VENUE. CPU only, banks only, no new realisations.
# ============================================================
def mode_venue():
    """Truth-control at the reserved cells: does N_cert > 0 EXIST under the standard pipeline?

    SPARK-2 5(5): "N_cert = 0 in BOTH states at both rungs means the current ladder may be
    ENTIRELY BELOW the certification bar at this cell -- check a louder cell or the truth
    control FIRST, before spending on a grid that can only return zeros." This is that check,
    and it is done from the banked SURFACE raw statistics: no GPU, no new draws.
    """
    sys.path.insert(0, f"{HSYMT}/CW_transition")
    from recut_surface import (offender, gumbel_floor, emp_quantile, boot_sd, adopt, QBAR as RQ)
    SR = f"{HSYMT}/SURFACE_results"

    def one(h, T, tier, k, lab, story=None):
        tag = cell_tag(h, T, tier, k)
        o = []
        for f in sorted(glob.glob(f"{SR}/sf_nullN_{tag}_g*_n*.npz")):
            d = np.load(f, allow_pickle=True)
            o.append(offender(np.nan_to_num(d["dlnL_det"], posinf=1e30), d["lnK"], d["qmax"]))
        o = np.array(o, float)
        gu, mu, beta, sd, n = gumbel_floor(o)
        zf = float((o == 0.0).mean())
        fl, err, est = adopt(gu, sd, o, zf)
        D, K, Q, OT, ncert = [], [], [], [], []
        for f in sorted(glob.glob(f"{SR}/sf_sig_{tag}_g*_n*.npz")):
            d = np.load(f, allow_pickle=True)
            dl = np.nan_to_num(d["dlnL_det"], posinf=1e30)
            lnK, q, ot = d["lnK"], d["qmax"], d["on_true"]
            D.append(dl); K.append(lnK); Q.append(q); OT.append(ot)
            cert = (dl > np.maximum(lnK, fl)) & (q > RQ)
            ncert.append((int(cert.sum()), int((cert & ot).sum()), int((cert & ~ot).sum())))
        ncert = np.array(ncert)
        nc = ncert[:, 0]
        print(f"{lab}", flush=True)
        print(f"   floor {fl:9.3f} +- {err:.3f} [{est}]  zf {zf:.2f}   nullN N = {len(o)}")
        print(f"   realisations {len(nc)}   N_cert/real: mean {nc.mean():.3f}  "
              f"median {np.median(nc):.1f}  max {nc.max()}")
        print(f"   >>> realisations with N_cert > 0 : {int((nc>0).sum())}/{len(nc)}  "
              f"({(nc>0).mean()*100:.0f}%)   LIVE = {(nc>0).any()}")
        print(f"   correct/real {ncert[:,1].mean():.3f}  wrong/real {ncert[:,2].mean():.3f}"
              + (f"   [STORY App A: {story}/real]" if story else ""))
        print(f"   N_cert vector: {list(nc)}\n", flush=True)
        return dict(h=h, T=T, tier=tier, k=k, floor=fl, err=err, est=est, zf=zf, n_null=len(o),
                    ncert=nc, corr=ncert[:, 1], wrong=ncert[:, 2], dlnL=np.array(D),
                    lnK=np.array(K), qmax=np.array(Q), on_true=np.array(OT))

    print("=" * 80)
    print("FIX 0 -- THE VENUE. Truth-control re-cut from the SURFACE banks (CPU, no new draws).")
    print("Nothing can cross a bar nobody reaches: SPARK-2 measured N_cert = 0 in every state.")
    print("=" * 80)
    R = {}
    for key, v in VENUES.items():
        R[key] = one(v["h"], v["T"], v["tier"], v["k"],
                     f"RESERVED {key}: ({v['h']}, {v['T']}, {v['tier']}, {v['k']}+{16-v['k']})",
                     story=v["story_corr"])
    R["C"] = one(SPARK2_CELL["h"], SPARK2_CELL["T"], SPARK2_CELL["tier"], SPARK2_CELL["k"],
                 "CONTROL:    (-13.25, 30, lit, 3+13) = SPARK-2's S2b venue")

    print("=" * 80)
    for key in ("A", "B"):
        live = (R[key]["ncert"] > 0).any()
        print(f"  VENUE {key}: N_cert > 0 EXISTS = {live}  "
              f"({int((R[key]['ncert']>0).sum())}/{len(R[key]['ncert'])} realisations)")
    print(f"  CONTROL (SPARK-2's cell): N_cert > 0 in "
          f"{int((R['C']['ncert']>0).sum())}/{len(R['C']['ncert'])} realisations")
    print("=" * 80, flush=True)

    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark3_venue.npz",
             **{f"{n}_{kk}": R[n][kk] for n in R for kk in
                ("h", "T", "tier", "k", "floor", "err", "est", "zf", "n_null", "ncert",
                 "corr", "wrong", "dlnL", "lnK", "qmax", "on_true")},
             note=("FIX 0, the venue check (SPARK-2 5 item 5). Truth-control re-cut of STORY "
                   "App A's RESERVED above-onset cells (SURFACE Pair B) and of SPARK-2's S2b "
                   "cell, from the banked SURFACE raw statistics under the standard pipeline "
                   "(recut_surface.py offender/adopt/score semantics verbatim). RAW per-pulsar "
                   "dlnL/lnK/qmax/on_true banked per realisation, never verdicts. The question "
                   "is EXISTENCE, not size: SPARK-2 measured N_cert = 0 in every state at its "
                   "cell, so its ladder could only return zeros. CPU only; no new realisations "
                   "-- a re-cut from the SURFACE banks suffices for the check."))
    print(f"[SPARK3] banked -> {OUT}/spark3_venue.npz")
    return 0


# ============================================================
# the venue build -- machinery REUSED, not reimplemented
# ============================================================
def _stack():
    import jax
    jax.config.update("jax_enable_x64", True)
    jax.config.update("jax_compilation_cache_dir",
                      os.environ.get("HARBOR_JAXCACHE_IN",
                                     "/home/mattm/.cache/jax_stagec_cache"))
    jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
    import jax.numpy as jnp
    import trackB_b1_core as C
    import trackB_estimator as TE
    import ignite as IG
    import forge_b1 as F
    import surface as SF
    from spark import force_recompute_lgwb
    return jax, jnp, C, TE, IG, F, SF, force_recompute_lgwb


def build_venue(key, geo_id, verbose=True):
    """The problem + geometry at the reserved cell, for sky `geo_id`.

    THE L_gwb PATH, AND WHY IT IS force_recompute (not a choice -- a reproduction
    requirement). The canonical `b1_L_gwb.npz` bank is (3248, 3248) = the T = 15 ANCHOR
    array; the venue is T = 40 and needs a different shape, so `load_or_build_L_gwb` would
    hit its shape-mismatch RAISE. SURFACE PREDATES the banking fix (CHORUS 2026-07-13 §0.1):
    with no bank on disk it took the RECOMPUTE branch (np.linalg.eigh of Phi_gwb) at the
    convention `--cpus-per-task=8`. To reproduce SURFACE's banked realisations -- whose
    certified sets this campaign cohere on -- SPARK-3 must take the SAME branch at the SAME
    thread count. Hence force_recompute_lgwb, and hence GATE g3a (mode `replay`), which
    PROVES the reproduction rather than assuming it. cpus-per-task=8 is part of the seed (the
    NoiseDrawer thread-count hazard); the sbatch sets it and the log prints the fingerprint
    so the provenance is on the record.
    """
    jax, jnp, C, TE, IG, F, SF, force_recompute_lgwb = _stack()
    v = VENUES[key]
    force_recompute_lgwb(C)
    P = IG.build_ignite_problem(v["T"], verbose=verbose)
    geo_src = None if geo_id < 0 else F.load_geo_skies([geo_id])[0]
    G = SF.cell_geometry_s(P, C, F, geo_src, v["h"], v["tier"], v["k"])
    return P, G, v


# ============================================================
# FIX 1 -- THE PER-TARGET E-STEP
# ============================================================
def estep_per_target(sp, theta_base, EV, AI, data_tuple, PM, jnp):
    """The certified-coherent E-step arrow 2 needs: for TARGET pulsar p, p's OWN pulsar term
    is live, the certified set's terms are coherent at their q, and the uncertified are off.

    WHY THIS EXISTS (SPARK-2 1, the structural finding). `B1Split.estep` sweeps EVERY
    pulsar's distance in ONE call under a GLOBAL pmask, so a pulsar at m_p = 0 has its own
    term switched off -> its distance is inert -> its fringe row is FLAT -> it cannot
    certify. At rung 0 the certification statistic does not exist at all (0/116 live rows).
    So the rung's mask cannot be applied globally; each target needs its OWN mask.

    PM: (npsr, npsr) -- PM[p] is the pmask in force when the target is p.

    THE COST, AND WHY IT IS NOT 116x. SPARK-2 priced this as "~116x ... a real build, not a
    parameter ... the single largest thing standing between this addendum and its verdict",
    reading it as 116 calls to `estep()` (each sweeping all 116 pulsars, 115/116 of the work
    discarded). It need not be built that way. `estep`'s own cost splits into
        (i)  ONE base evaluation `_ppab` at the pmask   -- 116 pulsar-evals + the full
             2668-dim GWB solve; and
        (ii) 116 per-pulsar fringe SWEEPS `_pulsar_ab_fn(p)` -- B = 512 evals EACH, on a
             per-pulsar GP block, and this is the bulk of the work.
    Only (i) depends on the target's mask; (ii) for target p is IDENTICAL whatever the other
    pulsars' mask is, because `MaskedDelay` reads `m = params[PMASK][self.ipsr]` -- pulsar
    p's own entry ONLY (trackB_b1_core.py:127). So the per-target form re-does (i) 116 times
    and (ii) exactly ONCE per target, as the standard E-step already does. Both `_ppab` and
    `_pulsar_ab_fn(p)` take pmask as a RUNTIME arg, so nothing recompiles.
    Measured cost is reported by mode `gate1` rather than asserted here.
    """
    sp.AI = np.asarray(AI)
    sp._ensure_minv()
    tb = jnp.asarray(theta_base)
    npsr, B = EV.shape
    LNL = np.empty((npsr, B))
    const = sp.a_const - 0.5 * (sp.ldP_cached + sp.logdet_cached)
    for p in range(npsr):
        pm = jnp.asarray(PM[p])
        a_base, b_base = sp._ppab(tb, data_tuple, pm)
        a_base = np.asarray(a_base); b_base = np.asarray(b_base)
        bflat = b_base.reshape(-1)
        u = sp._Minv @ bflat
        Q_base = float(bflat @ u)
        u_p = u.reshape(sp.npsr, sp.ngp_gwb)
        sum_a = float(a_base.sum())
        a_pf, b_pf = sp._pulsar_ab_fn(p)(tb, jnp.asarray(EV[p]), data_tuple, pm)
        a_pf = np.asarray(a_pf); b_pf = np.asarray(b_pf)
        db = b_pf - b_base[p]
        Q_pf = (Q_base + 2.0 * (db @ u_p[p])
                + np.einsum('bi,ij,bj->b', db, sp._Minv_pp[p], db))
        LNL[p] = const + (sum_a - a_base[p] + a_pf) + 0.5 * Q_pf
    return LNL


def rung_masks(npsr, cert_idx, qvals, k, rng=None):
    """PM for rung k. PM[p][p] = 1 ALWAYS (the target's own term is live -- else its fringe
    row is flat and it cannot certify, SPARK-2 1); PM[p][j] = q_j for the k certified
    pulsars j != p; 0 elsewhere (uncertified = DECOHERED, not pinned to a wrong MAP fringe --
    the hard-lock failure IGNITE-2 retired).

    The rung's k pulsars are the k certified with the LARGEST dlnL -- the array's own
    data-driven ranking, drawn from THIS realisation's own certified set (not imported from
    another cell's igniter bank)."""
    PM = np.zeros((npsr, npsr))
    take = np.asarray(cert_idx[:k], int) if k > 0 else np.zeros(0, int)
    for p in range(npsr):
        if len(take):
            PM[p, take] = np.asarray(qvals)[take]
        PM[p, p] = 1.0
    # pm_rung: the rung's certified-set mask with NO target override. This is the state of the
    # ARRAY at this rung, and it is what (a)'s Fisher must be evaluated at -- the faint
    # sources' constraint is a property of which pulsars are coherent, NOT of which pulsar is
    # currently being scored. Using a PM row here would force that row's target live even when
    # uncertified and silently add a 117th coherent term to the Fisher's state.
    pm_rung = np.zeros(npsr)
    if len(take):
        pm_rung[take] = np.asarray(qvals)[take]
    return PM, take, pm_rung


def venue_realisation(P, G, C, IG, jnp, tier, noise_seed, dist_seed, no_cw=False,
                      faint_present=True, k_loud=5):
    """One Arm-B realisation at the venue, drawn EXACTLY as surface.py:191-215 draws it
    (same geometry, same seeds, same order) so it reproduces the banked sf_sig realisation.
    faint_present=False removes the reservoir from the DATA (log10_h -> H_ABSENT) but keeps
    the loud: the matched null for the deconfusion floor."""
    dL, EV = G["dL"], G["EV"]
    nd, AI = P["n_dist"], P["AI"]
    L_true, n_true = IG.draw_true_distances_tier(P, dL, EV, seed=dist_seed, tier=tier)
    theta_true = G["theta"].copy()
    theta_true[AI] = L_true
    if not faint_present:
        for i in range(k_loud, P["ncw"]):
            theta_true[nd + 8 * i + I_H] = H_ABSENT
    one = jnp.ones(P["npsr"])
    if no_cw:
        clean = tuple(jnp.zeros(len(p.toas)) for p in P["disco_psrs"])
    else:
        clean = P["amo"]["inject_delay"](jnp.asarray(theta_true), one)
    noise = P["nd"].draw(noise_seed)
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n)) for c, n in zip(clean, noise))
    return dict(data=data, theta_true=theta_true, L_true=L_true, n_true=n_true)


# ============================================================
# FIX 2 -- THE EVALUABILITY BOUNDARY (the double-perturbation pathology, diagnosed)
# ============================================================
TSUN = 4.925490947e-6          # discovery const.Tsun (s)
TREF_MJD = 51544.5             # deterministic.py:499 -- tref = J2000
YR = 365.25 * 86400.0

# the POPULATION's generative prior (stagec_fisher.py:57-58) -- what actually drew the truth
POP_MC_RANGE = (8.5, 9.5)
POP_FGW_RANGE = (-8.0, -7.5)
# TE.phi_bounds (trackB_estimator.py:512-513) -- the ESTIMATOR's SEARCH box. SPARK-2 clipped
# the soft draw to THIS, believing it to be "the prior support". It is not: it is 3.0 dex in
# mc against the population's 1.0, and 1.7 dex in fgw against the population's 0.5.
PHI_MC_RANGE = (7.0, 10.0)
PHI_FGW_RANGE = (-8.7, -7.0)


def t_coal(log10_mc, log10_fgw):
    """Time from tref (J2000) to coalescence, seconds.

    THE NaN, EXACTLY (discovery/deterministic.py:509-510, the model this programme injects
    and recovers with):

        term  = 1.0 - (256.0/5.0) * mc**(5/3) * w0**(8/3) * t
        omega = w0 * jnp.power(term, -3.0/8.0)

    `jnp.power(negative, -3/8)` is NaN. term < 0 <=> t > t_coal = 5/(256 mc^(5/3) w0^(8/3)),
    i.e. the binary has ALREADY COALESCED at epoch t and the chirping waveform does not exist.

    WHICH TERM CAN FIRE. For the EARTH term t = toas - tref > 0 (the TOAs postdate J2000), so
    it CAN cross. For the PULSAR term t = tp = toas_rel - (L/c)(1 - cos_mu), and L/c ~ kpc/c
    ~ 3.3 kyr, so tp is hugely NEGATIVE -> term = 1 + |..| > 1 -> always evaluable. The
    pathology is an EARTH-TERM pathology, and it therefore fires at EVERY pmask -- coherent or
    not. This is why it wiped rungs 2 and 5 regardless of the certified set.
    """
    mc = (10.0 ** np.asarray(log10_mc, float)) * TSUN
    w0 = np.pi * (10.0 ** np.asarray(log10_fgw, float))
    return 5.0 / (256.0 * mc ** (5.0 / 3.0) * w0 ** (8.0 / 3.0))


def t_max_of(P):
    """max over the array of (toa - tref): the largest epoch the Earth term must evaluate at.
    Computed from the ACTUAL TOAs, never assumed from T_label (the extension moves it)."""
    tref = 86400.0 * TREF_MJD
    return float(max(float(np.max(np.asarray(p.toas))) for p in P["disco_psrs"]) - tref)


def evaluable(log10_mc, log10_fgw, t_max):
    """The exact, ANALYTIC evaluability predicate. No GPU, no waveform call."""
    return t_coal(log10_mc, log10_fgw) > t_max


def score_from_LNL(P, FT, LNL, C, L_true, dL):
    """ignite.py:349-358 semantics, verbatim -- the banked certification convention, applied
    to whichever LNL it is handed."""
    C._finite("estep", LNL)
    post = FT.posterior(LNL)
    n_true_grid = C.n_true_on_grid(L_true, P["L0"], dL)
    qmax = post["q_max"]; mapk = post["map_k"]
    lnK = np.log(np.maximum(FT.K_counted.astype(float), 1.0))
    dlnL_det = np.array([
        (np.sort(FT.peak[a])[-1] - np.sort(FT.peak[a])[-2])
        if np.asarray(FT.peak[a]).size > 1 else np.inf for a in range(P["npsr"])])
    return dict(dlnL=dlnL_det, lnK=lnK, qmax=qmax, mapk=mapk,
                on_true=(mapk == n_true_grid), n_true_grid=n_true_grid)


# ============================================================
# FIX 3 -- THE FISHER, TWO-SIDED. The joint is measured unaffordable; bound it instead.
# ============================================================
def faint_fisher_bounds(P, R, pmask, Ld, jnp, faint_lo):
    """The two Fisher bounds arrow 2's tightening is quoted BETWEEN.

    THE JOINT IS NOT AFFORDABLE, AND THAT IS MEASURED, NOT ASSUMED (SPARK-2 4a). The marginal
    width needs the full 78x78 Hessian of lnL over the faint block. TWO builds failed to
    return on an A100-80GB inside a 1 h walltime: a monolithic `jax.hessian` (job 12583525,
    >27 min on the Hessian, cancelled) and a chunked batched-JVP-of-grad Hessian at CH=8 (job
    12583899, >17 min, cancelled). The likelihood carries a 2668-dim correlated-GP solve per
    evaluation; 78 JVPs of its gradient is not a few-GPU-hr object. SPARK-3 does NOT retry it.

    THE PRE-REGISTERED REPLACEMENT -- bound the answer on BOTH sides and let the verdict fall
    where it falls:

      OPTIMISTIC  sigma_cond = 1/sqrt(F_ii), the CONDITIONAL width (a 1-D curvature scan per
        axis; one batched logL call per chunk, seconds). It IGNORES covariance with the other
        faint parameters, and sigma_cond <= sigma_marg ALWAYS. So it OVERSTATES how well the
        reservoir is known -> OVERSTATES how well it can be modelled -> FAVOURS arrow 2. An
        EDGE-ZERO measured here is a fortiori EDGE-ZERO with the true joint widths.

      PESSIMISTIC  sigma_prior = the prior box width. The loop learns NOTHING about the
        reservoir; the soft model is a draw from the prior. This CANNOT overstate the
        tightening, because there is none.

    The honest answer lives BETWEEN. If both give the SAME verdict, the verdict stands
    WITHOUT the joint -- the joint could only land between two bounds that already agree. If
    they STRADDLE, the bounds do not decide, and the report says so and prices the chunked JVP
    rather than quoting a number it cannot defend.

    Both are banked as COLUMNS on every unit (SPARK-2's labelling convention), never merged.
    """
    nd = P["n_dist"]
    faint = list(range(faint_lo, P["ncw"]))
    idx = np.array([nd + 8 * k + j for k in faint for j in FAINT_PARS])
    th = R["theta_true"].copy(); th[P["AI"]] = np.asarray(Ld, float)
    lb = P["amo"]["logL_batch_theta"]
    data = R["data"]; pm = jnp.asarray(pmask)
    box = np.array([BOX[j] for _ in faint for j in FAINT_PARS])
    step = 1e-3 * box
    n = len(idx)
    stack = []
    for i in range(n):
        for d in (-1, 0, 1):
            t = th.copy(); t[idx[i]] += d * step[i]
            stack.append(t)
    stack = np.stack(stack)
    lls = []
    CH = 48
    for c0 in range(0, len(stack), CH):
        lls.append(np.asarray(lb(jnp.asarray(stack[c0:c0 + CH]), data, pm)))
    L = np.concatenate(lls).reshape(n, 3)
    d2 = (L[:, 0] - 2.0 * L[:, 1] + L[:, 2]) / (step ** 2)
    F_ii = -d2
    cond_sig = np.where(F_ii > 0, 1.0 / np.sqrt(np.maximum(F_ii, 1e-300)), np.inf)
    # OPTIMISTIC: the conditional width, capped at the prior box (an unconstrained axis is
    # known no better than its prior -- that is the correct statement of "we know nothing").
    sig_opt = np.where(np.isfinite(cond_sig), np.minimum(cond_sig, box), box)
    # PESSIMISTIC: the prior box itself. No tightening at any rung.
    sig_pes = box.copy()
    return dict(idx=idx, F_ii=F_ii, cond_sig=cond_sig, sig_opt=sig_opt, sig_pes=sig_pes,
                box=box, n_pos=int((F_ii > 0).sum()), n_par=n, faint=np.array(faint),
                cond_num=float(np.max(np.abs(F_ii)) / max(np.min(np.abs(F_ii)), 1e-300)),
                bound_labels=np.array(["OPTIMISTIC=conditional 1/sqrt(F_ii) (<= marginal; "
                                       "favours arrow 2)",
                                       "PESSIMISTIC=prior width (no tightening at all)"]))


def soft_faint_theta(P, R, sig, idx, rng, bounds, t_max, faint_lo, max_tries=64):
    """A SOFT draw of the faint reservoir: truth + N(0, sigma), REJECTED AND REDRAWN where the
    binary has already coalesced. Never a hard commit.

    WHY REJECT-AND-REDRAW RATHER THAN A WIDER CLIP (FIX 2; the correction to SPARK-2 4b).
    SPARK-2 clipped the draw to `TE.phi_bounds` and called it "truncated at its own prior
    support, exactly as a real posterior is". phi_bounds is NOT the prior support -- it is the
    ESTIMATOR's SEARCH box: mc [7, 10] against the population's [8.5, 9.5], and fgw
    [-8.7, -7.0] against the population's [-8.0, -7.5]. That extra half-dex in EACH is exactly
    where the binary coalesces inside the span, so the clip SPARK-2 added as the fix is what
    ADMITS the non-evaluable region. One coalesced source is enough to kill everything: the
    delay SUMS all 16 sources (MaskedDelay.__call__), so a single NaN source NaNs every
    pulsar's residual -> the E-step returns non-finite on ALL npsr*B entries. That is exactly
    SPARK-2's signature (59392/59392 = 116 x 512), and it fired on rungs 2 and 5 and not 0 and
    8 for no deeper reason than which draws the RNG happened to put across the boundary.

    IT IS A DOUBLE PERTURBATION -- BUT OF (mc x fgw), NOT KINDLE's (scrambled x fix_mc). The
    boundary t_coal(mc, fgw) > t_max is a JOINT curve. EMBER's `mc_scan` (ember_map.py:133)
    scans mc with fgw HELD AT TRUTH and therefore never sees the corner; that -- not "loud vs
    faint" -- is why it measured "non-evaluable 0% across the entire mc box". SPARK-2 read
    that 0% as a loud-source result that "does not transfer to the faint". The strain plays no
    part in the mechanism at all: `t_coal` is a function of (mc, fgw) only.

    So the fix generalises EMBER's precedent from a per-axis scan to the JOINT analytic
    predicate: draw, test t_coal > t_max exactly (no waveform call), redraw the offending
    SOURCES only. The rejection fraction is BANKED per rung -- it is a genuine statement about
    how soft "soft" can be, not a nuisance to suppress.
    """
    th = R["theta_true"].copy()
    lo, hi = bounds
    s = np.asarray(sig, float)
    n_par_src = len(FAINT_PARS)
    n_rej = 0
    drawn = np.clip(th[idx] + rng.normal(0.0, s), lo, hi)
    for si in range(len(idx) // n_par_src):
        sl = slice(si * n_par_src, (si + 1) * n_par_src)
        j_mc = list(FAINT_PARS).index(I_MC)
        j_fg = list(FAINT_PARS).index(I_FGW)
        for _ in range(max_tries):
            v = drawn[sl]
            if evaluable(v[j_mc], v[j_fg], t_max):
                break
            n_rej += 1
            drawn[sl] = np.clip(th[idx][sl] + rng.normal(0.0, s[sl]), lo[sl], hi[sl])
        else:
            # the draw cannot be made evaluable at this width: fall back to TRUTH for this
            # source and BANK it. Never silently accept a NaN-producing model.
            drawn[sl] = th[idx][sl]
    th[idx] = drawn
    return th, n_rej


def _sig_files(key, geo_id):
    v = VENUES[key]
    return sorted(glob.glob(f"{HSYMT}/SURFACE_results/sf_sig_"
                            f"{cell_tag(v['h'], v['T'], v['tier'], v['k'])}_g{geo_id}_n*.npz"))


def _lane_tag():
    """The artifact key must carry the LANE. Learned the hard way: an H200 replay silently
    OVERWROTE the A100 replay's PASS because the filename keyed only on (venue, geo), leaving
    the report citing 0.000e+00 while the npz on disk said FAILED. An artifact whose name does
    not distinguish the conditions it was produced under is not an artifact, it is a race."""
    import socket
    h = socket.gethostname().split(".")[0]
    try:
        import subprocess
        g = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                           capture_output=True, text=True, timeout=20).stdout.split("\n")[0]
        g = g.strip().replace(" ", "")
    except Exception:
        g = "unknownGPU"
    return f"{h}_{g}"


def mode_replay(key, geo_id=3):
    """GATE g3a -- THE REPRODUCTION GATE. Re-draw a BANKED SURFACE realisation at the venue
    from its own seeds and reproduce its banked raw columns.

    WHY THIS GATE IS NOT OPTIONAL. This campaign coheres on the certified sets of SURFACE's
    banked realisations at this cell. If SPARK-3's data draw is not SURFACE's draw, the
    certified set belongs to a DIFFERENT realisation than the data being scored, and every
    rung is then coherence on a counterpart that is not there -- EMBER's manufacturing regime
    entered by a bookkeeping error rather than a physics one. The L_gwb recompute branch is
    thread-count dependent (the NoiseDrawer hazard), so the reproduction must be PROVEN at
    cpus-per-task=8, not assumed. It also transitively gates the whole build: problem,
    structure axis, tier, distance draw, noise draw, E-step and scorer.
    Bar: the discrete columns EXACT; the continuous columns at EMBER 2.2(b)'s adopted bar
    (< 1e-6). Anything worse -> STOP.
    """
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    fs = _sig_files(key, geo_id)
    if not fs:
        print(f"*** STOP: no banked sf_sig for venue {key} g{geo_id}", flush=True)
        return 2
    z = np.load(fs[0], allow_pickle=True)
    ns, ds = int(z["noise_seed"]), int(z["dist_seed"])
    v = VENUES[key]
    print(f"[SPARK3] REPLAY venue {key} = ({v['h']}, {v['T']}, {v['tier']}, {v['k']}+"
          f"{16-v['k']})  g{geo_id}  noise {ns} dist {ds}", flush=True)
    print(f"[SPARK3] target bank: {os.path.basename(fs[0])}", flush=True)

    t0 = time.time()
    P, G, v = build_venue(key, geo_id)
    print(f"[SPARK3] build {time.time()-t0:.0f}s", flush=True)
    R = venue_realisation(P, G, C, IG, jnp, v["tier"], ns, ds, k_loud=v["k"])
    one = np.ones(P["npsr"])
    th = R["theta_true"].copy(); th[P["AI"]] = P["L0"]      # E-step scans about the PRIOR MEAN
    t0 = time.time()
    LNL = P["sp"].estep(th, G["EV"], P["AI"], R["data"], one)
    print(f"[SPARK3] estep (first, incl. compile) {time.time()-t0:.1f}s", flush=True)
    s = score_from_LNL(P, G["FT"], LNL, C, R["L_true"], G["dL"])

    d_dlnL = float(np.max(np.abs(np.nan_to_num(s["dlnL"], posinf=1e30)
                                 - np.nan_to_num(z["dlnL_det"], posinf=1e30))))
    d_lnK = float(np.max(np.abs(s["lnK"] - z["lnK"])))
    d_q = float(np.max(np.abs(s["qmax"] - z["qmax"])))
    n_mapk = int((s["mapk"] != z["mapk"]).sum())
    n_ot = int((s["on_true"] != z["on_true"]).sum())
    print(f"\n  GATE g3a -- replay vs bank")
    print(f"    max |dlnL_det  - banked| = {d_dlnL:.3e}")
    print(f"    max |lnK       - banked| = {d_lnK:.3e}")
    print(f"    max |qmax      - banked| = {d_q:.3e}")
    print(f"    mapk    disagreements    = {n_mapk} / {P['npsr']}   (discrete: must be 0)")
    print(f"    on_true disagreements    = {n_ot} / {P['npsr']}   (discrete: must be 0)")
    ok = (n_mapk == 0 and n_ot == 0 and d_q < 1e-6 and d_lnK < 1e-6 and d_dlnL < 1e-6)
    print(f"    ==> g3a {'PASSED' if ok else '*** FAILED'}", flush=True)
    os.makedirs(OUT, exist_ok=True)
    lane = _lane_tag()
    np.savez(f"{OUT}/spark3_replay_{key}_g{geo_id}_{lane}.npz",
             venue=key, geo_id=geo_id, noise_seed=ns, dist_seed=ds, bank=fs[0], passed=ok,
             lane=lane, lgwb_fingerprint=str(getattr(C, "_SPARK3_LGWB_FP", "unrecorded")),
             d_dlnL=d_dlnL, d_lnK=d_lnK, d_qmax=d_q, n_mapk_disagree=n_mapk,
             n_ontrue_disagree=n_ot,
             dlnL=s["dlnL"], lnK=s["lnK"], qmax=s["qmax"], mapk=s["mapk"],
             dlnL_bank=z["dlnL_det"], lnK_bank=z["lnK"], qmax_bank=z["qmax"],
             mapk_bank=z["mapk"],
             note=("GATE g3a, the REPRODUCTION gate. SPARK-3 coheres on the certified sets of "
                   "SURFACE's banked realisations at this cell, so its data draw MUST be "
                   "SURFACE's draw -- otherwise the certified set belongs to a different "
                   "realisation than the data scored. The venue is T=40; the canonical "
                   "b1_L_gwb bank is the T=15 shape, so this takes the RECOMPUTE branch, "
                   "which is thread-count dependent (the NoiseDrawer hazard) -- hence the "
                   "reproduction is PROVEN at cpus-per-task=8, not assumed. Raw columns "
                   "banked beside the bank's own, never a verdict."))
    print(f"[SPARK3] banked -> {OUT}/spark3_replay_{key}_g{geo_id}_{lane}.npz", flush=True)
    return 0 if ok else 1


def mode_gate1(key, geo_id=3):
    """FIX 1's GATE CHAIN + the measured cost of the per-target E-step.

      g1a  THE CODE-PATH GATE: the per-target machinery driven at PM[p] == ONE for every
           target must reproduce `sp.estep(..., one)` BIT-IDENTICALLY. Same arithmetic, just
           116 redundant base evaluations. This proves the refactor did not move the
           incumbent (SPARK's g0a in spirit). Bar: 0.000e+00.
      g1b  THE PHYSICS GATE (the brief's rung-0 gate): the per-target form at RUNG 0 -- each
           target's own term live, NOTHING else coherent -- against the banked standard
           E-step's q columns at pmask = ONE.

    WHAT g1b IS ACTUALLY TESTING, derived from the machinery BEFORE it is run. Under the
    per-target rung-0 mask, target p has m_p = 1 exactly as it does under pmask = ONE; the
    other pulsars are switched off. `MaskedDelay` reads only `params[PMASK][self.ipsr]`
    (:127), so pulsar p's own residual -- hence `a_pf` and `db` -- is IDENTICAL in the two.
    The ONLY route by which the other pulsars' mask can reach target p's fringe row is the
    HD-CORRELATED GWB coupling: `u = Minv @ bflat` mixes every pulsar's `b_base`, and the
    cross-term `2 * db @ u_p[p]` is L-dependent. The constant offsets (`sum_a`, `const`)
    cancel in the normalised posterior. So g1b measures ONE thing: how much the GWB
    correlation couples other pulsars' pulsar-term state into this pulsar's fringe posterior.
    A pass says the rung-0 per-target E-step IS the banked convention to within that
    coupling; a fail is not a bug but a MEASUREMENT -- it would mean the two conventions are
    different objects and the ladder must state which one it is anchored to.
    """
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    fs = _sig_files(key, geo_id)
    z = np.load(fs[0], allow_pickle=True)
    ns, ds = int(z["noise_seed"]), int(z["dist_seed"])
    P, G, v = build_venue(key, geo_id)
    R = venue_realisation(P, G, C, IG, jnp, v["tier"], ns, ds, k_loud=v["k"])
    npsr = P["npsr"]
    one = np.ones(npsr)
    th = R["theta_true"].copy(); th[P["AI"]] = P["L0"]

    # ---- the standard E-step, warm ----
    L_std = P["sp"].estep(th, G["EV"], P["AI"], R["data"], one)   # compile
    t0 = time.time()
    L_std = P["sp"].estep(th, G["EV"], P["AI"], R["data"], one)
    t_std = time.time() - t0
    print(f"[SPARK3] standard estep (warm): {t_std:.2f}s", flush=True)

    # ---- g1a: per-target machinery at PM == ONE ----
    PM_one = np.ones((npsr, npsr))
    t0 = time.time()
    L_pt_one = estep_per_target(P["sp"], th, G["EV"], P["AI"], R["data"], PM_one, jnp)
    t_pt = time.time() - t0
    d_g1a = float(np.max(np.abs(L_pt_one - L_std)))
    print(f"[SPARK3] per-target estep (warm): {t_pt:.2f}s   = {t_pt/max(t_std,1e-9):.1f}x the "
          f"standard  [SPARK-2 priced this at ~116x]", flush=True)
    print(f"\n  GATE g1a (code path, PM==ONE): max |LNL_per_target - LNL_standard| = "
          f"{d_g1a:.3e}   ==> {'PASSED' if d_g1a == 0.0 else '*** FAILED'}", flush=True)

    # ---- g1b: per-target rung 0 vs the banked standard ----
    PM_r0 = np.eye(npsr)
    L_pt_r0 = estep_per_target(P["sp"], th, G["EV"], P["AI"], R["data"], PM_r0, jnp)
    s_r0 = score_from_LNL(P, G["FT"], L_pt_r0, C, R["L_true"], G["dL"])
    s_std = score_from_LNL(P, G["FT"], L_std, C, R["L_true"], G["dL"])
    live_r0 = int((np.array([np.ptp(np.asarray(G["FT"].peak[a])) for a in range(npsr)])
                   > 1e-9).sum())
    d_q = float(np.max(np.abs(s_r0["qmax"] - np.asarray(z["qmax"]))))
    d_d = float(np.max(np.abs(np.nan_to_num(s_r0["dlnL"], posinf=1e30)
                              - np.nan_to_num(z["dlnL_det"], posinf=1e30))))
    n_mapk = int((s_r0["mapk"] != z["mapk"]).sum())
    print(f"\n  GATE g1b (rung 0 per-target vs the BANKED standard q columns)")
    print(f"    LIVE fringe rows at rung 0 : {live_r0}/{npsr}   "
          f"[SPARK-2's global-pmask rung 0: 0/116 -- the confound this fixes]")
    print(f"    max |qmax - banked|        = {d_q:.3e}")
    print(f"    max |dlnL - banked|        = {d_d:.3e}")
    print(f"    mapk disagreements         = {n_mapk}/{npsr}")
    print(f"    (the ONLY channel between them is the HD-correlated GWB coupling "
          f"2*db@u_p[p]; see docstring)", flush=True)

    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark3_gate1_{key}_g{geo_id}.npz",
             venue=key, geo_id=geo_id, noise_seed=ns, dist_seed=ds,
             t_standard=t_std, t_per_target=t_pt, cost_ratio=t_pt / max(t_std, 1e-9),
             g1a_maxdev=d_g1a, g1a_passed=(d_g1a == 0.0),
             g1b_d_qmax=d_q, g1b_d_dlnL=d_d, g1b_n_mapk_disagree=n_mapk,
             n_live_rung0=live_r0,
             qmax_r0=s_r0["qmax"], dlnL_r0=s_r0["dlnL"], mapk_r0=s_r0["mapk"],
             qmax_std=s_std["qmax"], dlnL_std=s_std["dlnL"], mapk_std=s_std["mapk"],
             qmax_bank=z["qmax"], dlnL_bank=z["dlnL_det"], mapk_bank=z["mapk"],
             note=("FIX 1: the per-target E-step + gate chain. g1a = the code-path gate "
                   "(PM==ONE must reproduce sp.estep bit-identically). g1b = the rung-0 gate "
                   "(target's own term live, nothing else coherent) vs the BANKED standard "
                   "convention at pmask=ONE. The two differ ONLY through the HD-correlated "
                   "GWB coupling 2*db@u_p[p]: MaskedDelay reads params[PMASK][ipsr] only, so "
                   "a_pf and db are identical. n_live_rung0 is the direct counter to SPARK-2 "
                   "1's global-pmask confound (0/116 live rows at rung 0). Cost ratio is "
                   "MEASURED, not asserted. Raw LNL-derived columns banked, never verdicts."))
    print(f"[SPARK3] banked -> {OUT}/spark3_gate1_{key}_g{geo_id}.npz", flush=True)
    return 0


def mode_venue_self(key, geo_id, n_null, seed_block=0):
    """DERIVE the venue's floor + certified sets FROM OUR OWN DATA, on whatever host we are on.

    WHY THIS EXISTS, AND WHY IT IS BETTER SCIENCE THAN THE REPLAY IT REPLACES. SPARK-3's ONLY
    dependence on SURFACE's banks was `certified_of`: it imported each realisation's certified
    set from `sf_sig_*`, which forced our data to BE SURFACE's data bit-for-bit, which forced
    gate g3a, which forced the canonical `L_gwb` basis, which forced the dgx lane. That chain
    existed to save ONE standard E-step per realisation. It is not a scientific requirement.

    BIT-IDENTITY WAS NEVER THE REQUIREMENT -- CONSISTENCY WAS. The GWB square root enters as
    `L_gwb @ z`; a rotated eigenvector basis (a different host, or a different BLAS thread
    count) yields a DIFFERENT but EQUAL-IN-DISTRIBUTION realisation. Nothing physical changes.
    And arrow 2's statistic is a WITHIN-UNIT RELATIVE one -- does an uncertified pulsar's
    margin cross its bar between rung 0 and rung 8, on the SAME data, with the floor re-cut
    from that state's OWN nulls. The basis cancels out of that comparison entirely. What must
    hold is that the data, the certified set, and the floor all come from ONE basis -- which is
    exactly what this mode guarantees by deriving all three here.

    So: one standard E-step at pmask = ONE (the banked certification convention, ignite.py:349)
    per g3 realisation, and `n_null` matched nullN offenders for the floor via RECUT `adopt()`.
    ~2 min of GPU. It makes the campaign HOST-INDEPENDENT and runs on any full-fp64 card.

    THE GATE IS DISTRIBUTIONAL, NOT BIT-WISE -- and that is the right instrument for the claim.
    A rotated basis cannot reproduce SURFACE's numbers exactly and SHOULD NOT: it is a
    different draw. What it must reproduce is the CELL: a floor near SURFACE's 122.461 (A) /
    44.397 (B) and a live count (N_cert > 0) at g3. Agreement there says the host computes the
    same physics; disagreement says it does not, and would refuse the lane on a real ground
    rather than on a fingerprint mismatch that only ever meant "different noise draw".
    """
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    v = VENUES[key]
    P, G, _ = build_venue(key, geo_id)
    FT, EV, dL = G["FT"], G["EV"], G["dL"]
    L0 = np.asarray(P["L0"], float)
    one = np.ones(P["npsr"])
    lane = _lane_tag()
    print(f"[SPARK3] venue-self {key} g{geo_id} on lane {lane}: deriving floor + certified "
          f"sets from OUR OWN data (no SURFACE bank, no bit-identity requirement)", flush=True)

    # ---- matched nullN floor, this host's basis, SURFACE's convention (standard E-step) ----
    offs = []
    for i in range(n_null):
        # PER-VENUE, PER-BLOCK null seeds. The first cut used one seed set for BOTH venues
        # -- and nullN has no CW, so both floors were cut from the SAME noise draws: two
        # CORRELATED estimates, not two independent ones. That is one measurement wearing two
        # hats, and it cannot be used as evidence about the host. Fixed; `seed_block` buys a
        # genuinely independent replicate of the same floor.
        ns = 40_000_000 + 900_000 + 100_000 * (0 if key == "A" else 1) \
             + 10_000 * seed_block + i
        Ltn, _ = IG.draw_true_distances_tier(P, dL, EV, seed=ns + 10_000_000, tier=v["tier"])
        noise = P["nd"].draw(ns)
        dn = tuple(jnp.asarray(np.asarray(b)) for b in noise)      # nullN: no CW in the data
        th = G["theta"].copy(); th[P["AI"]] = L0
        s = score_from_LNL(P, FT, P["sp"].estep(th, EV, P["AI"], dn, one), C, Ltn, dL)
        offs.append(offender(s["dlnL"], s["lnK"], s["qmax"])[0])
        if (i + 1) % 25 == 0:
            print(f"   null {i+1}/{n_null}", flush=True)
    from spark import adopt
    FL = adopt(np.array(offs))
    print(f"[SPARK3] floor {FL['floor']:.3f} +- {FL['err']:.3f} [{FL['estimator']}] "
          f"zf={FL['zf']:.2f}   [SURFACE banked: {'122.461' if key=='A' else '44.397'}]",
          flush=True)

    # ---- the certified set per g3 realisation, at OUR data ----
    fs = _sig_files(key, geo_id)
    rows = []
    for r, f in enumerate(fs):
        z = np.load(f, allow_pickle=True)
        ns, ds = int(z["noise_seed"]), int(z["dist_seed"])       # SURFACE's SEEDS, our basis
        R = venue_realisation(P, G, C, IG, jnp, v["tier"], ns, ds, k_loud=v["k"])
        th = R["theta_true"].copy(); th[P["AI"]] = L0
        s = score_from_LNL(P, FT, P["sp"].estep(th, EV, P["AI"], R["data"], one), C,
                           R["L_true"], dL)
        cert = (s["dlnL"] > np.maximum(s["lnK"], FL["floor"])) & (s["qmax"] > QBAR)
        rows.append(dict(real=r, noise_seed=ns, dist_seed=ds, n_cert=int(cert.sum()),
                         n_corr=int((cert & s["on_true"]).sum()),
                         dlnL=s["dlnL"], lnK=s["lnK"], qmax=s["qmax"], cert=cert,
                         on_true=s["on_true"]))
        print(f"   real {r}: N_cert {int(cert.sum())}  (banked SURFACE: "
              f"{int(((np.nan_to_num(z['dlnL_det'],posinf=1e30) > np.maximum(z['lnK'], 122.461 if key=='A' else 44.397)) & (z['qmax']>QBAR)).sum())})",
              flush=True)
    nc = np.array([r["n_cert"] for r in rows])
    print(f"\n  DISTRIBUTIONAL GATE (not bit-wise -- a rotated basis is a different draw and "
          f"must not match exactly)")
    print(f"    floor      {FL['floor']:8.3f}   vs SURFACE {'122.461' if key=='A' else '44.397'}")
    print(f"    N_cert/g3  mean {nc.mean():.2f}  min {nc.min()}  max {nc.max()}  "
          f"live {int((nc>0).sum())}/{len(nc)}   vs SURFACE g3 "
          f"{'10-14' if key=='A' else '8-15'}")
    print(f"    rung-8 reachable: {int((nc>=8).sum())}/{len(nc)}", flush=True)
    os.makedirs(OUT, exist_ok=True)
    path = (f"{OUT}/spark3_venueself_{key}_g{geo_id}_{lane}"
            + (f"_b{seed_block}" if seed_block else "") + ".npz")
    np.savez(path, venue=key, geo_id=geo_id, lane=lane, floor=FL["floor"], err=FL["err"],
             zf=FL["zf"], estimator=FL["estimator"], null_off=np.array(offs),
             n_null=n_null, seed_block=seed_block,
             real=np.array([r["real"] for r in rows]),
             noise_seed=np.array([r["noise_seed"] for r in rows]),
             dist_seed=np.array([r["dist_seed"] for r in rows]),
             n_cert=nc, n_corr=np.array([r["n_corr"] for r in rows]),
             dlnL=np.array([r["dlnL"] for r in rows]),
             lnK=np.array([r["lnK"] for r in rows]),
             qmax=np.array([r["qmax"] for r in rows]),
             cert=np.array([r["cert"] for r in rows]),
             on_true=np.array([r["on_true"] for r in rows]),
             note=("SPARK-3 SELF-DERIVED venue: floor + certified sets computed from OUR OWN "
                   "data on THIS host's GWB basis, at SURFACE's seeds and SURFACE's "
                   "certification convention (standard E-step, pmask=ONE, ignite.py:349). "
                   "This CUTS the campaign's only dependence on SURFACE's banks and with it "
                   "the bit-identity requirement: a rotated L_gwb basis is a DIFFERENT but "
                   "EQUAL-IN-DISTRIBUTION draw, and arrow 2's statistic is a WITHIN-UNIT "
                   "relative comparison whose floor is re-cut from its own nulls -- the basis "
                   "cancels. What must hold is that data, certified set and floor share ONE "
                   "basis, which this guarantees. Gate is DISTRIBUTIONAL (floor near SURFACE's, "
                   "count live), never bit-wise. Raw dlnL/lnK/qmax + raw null offenders banked."))
    print(f"[SPARK3] banked -> {path}", flush=True)
    return 0


def certified_of(key, geo_id, real_i, self_lane=None):
    """THIS realisation's OWN certified set, re-derived from RAW columns at the venue's
    adopted floor -- 'do not trust the banked bool' (ember_anatomy.py:26). The rung's coherent
    pulsars are drawn from here, NOT imported from another cell's igniter bank.

    `self_lane`: prefer a SELF-DERIVED venue bank (mode `venue_self`) produced on this host's
    own GWB basis. That is what makes the campaign host-independent -- see mode_venue_self."""
    if self_lane:
        zs = sorted(glob.glob(f"{OUT}/spark3_venueself_{key}_g{geo_id}_*.npz"))
        if zs:
            z = np.load(zs[0], allow_pickle=True)
            i = int(np.where(np.asarray(z["real"]) == real_i)[0][0])
            dlnL = np.asarray(z["dlnL"][i], float); qmax = np.asarray(z["qmax"][i], float)
            cert = np.asarray(z["cert"][i], bool)
            order = np.where(cert)[0][np.argsort(-dlnL[np.where(cert)[0]])]
            return dict(cert=cert, qmax=qmax, dlnL=dlnL, lnK=np.asarray(z["lnK"][i], float),
                        floor=float(z["floor"]), on_true=np.asarray(z["on_true"][i], bool),
                        order=order, n_cert=int(cert.sum()),
                        noise_seed=int(z["noise_seed"][i]), dist_seed=int(z["dist_seed"][i]),
                        path=zs[0])
        print(f"*** no self-derived venue bank for {key} g{geo_id}; run `venue_self` first",
              flush=True)
        raise SystemExit(2)
    v = VENUES[key]
    z = np.load(f"{OUT}/spark3_venue.npz", allow_pickle=True)
    floor = float(z[f"{key}_floor"])
    fs = _sig_files(key, geo_id)
    d = np.load(fs[real_i], allow_pickle=True)
    dlnL = np.nan_to_num(np.asarray(d["dlnL_det"], float), posinf=1e30)
    lnK = np.asarray(d["lnK"], float); qmax = np.asarray(d["qmax"], float)
    cert = (dlnL > np.maximum(lnK, floor)) & (qmax > QBAR)
    order = np.where(cert)[0][np.argsort(-dlnL[np.where(cert)[0]])]
    return dict(cert=cert, qmax=qmax, dlnL=dlnL, lnK=lnK, floor=floor,
                on_true=np.asarray(d["on_true"], bool), order=order,
                n_cert=int(cert.sum()), noise_seed=int(d["noise_seed"]),
                dist_seed=int(d["dist_seed"]), path=fs[real_i])


def mode_gate2(key, geo_id=3, real_i=0):
    """FIX 2 -- THE ANATOMY OF THE DOUBLE PERTURBATION, and the gate on the repaired draw.

    SPARK-2 4b called this "IGNITE-2/KINDLE's known double-perturbation pathology (scrambled +
    fix_mc -> degenerate covariance -> all-non-finite E-step) arriving by a new route", named
    `log10_mc` x `fgw` coalescence as a "remaining suspect", and prescribed reject-and-redraw.
    The suspect is the whole story, the KINDLE attribution is wrong, and none of it needs a
    GPU to establish -- the boundary is ANALYTIC (see `t_coal`). This mode measures it, then
    proves the analytic predicate against the actual E-step rather than trusting it.

      g2a  THE PREDICATE IS EXACT: an E-step at a draw the predicate calls NON-evaluable must
           return non-finite on ALL npsr*B entries; at a draw it calls evaluable, on none.
      g2b  SPARK-2's RECIPE REPRODUCED: under its exact clip (phi_bounds) the per-rung
           non-evaluable fraction is > 0 -- i.e. its rung 2/5 deaths are reproduced as a
           property of the DRAW BOX, not of the rung.
      g2c  THE REPAIRED DRAW IS CLEAN: under the joint analytic predicate + reject-and-redraw,
           every rung's N_SOFT draws are evaluable and the E-step is finite. The rejection
           fraction is BANKED per rung.
    """
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    P, G, v = build_venue(key, geo_id)
    t_max = t_max_of(P)
    print(f"\n[SPARK3] t_max (max TOA - J2000) = {t_max/YR:.2f} yr   "
          f"(T_label {v['T']}; computed from the ACTUAL TOAs, not T_label)", flush=True)

    # ---- the analytic map ----
    print("\n  THE EVALUABILITY BOUNDARY  t_coal = 5/(256 mc^5/3 w0^8/3)  "
          "[deterministic.py:509-510]")
    print(f"  {'box':<40} {'log10_mc':>9} {'log10_fgw':>10} {'t_coal(yr)':>11} {'':>14}")
    corners = [("POPULATION prior worst corner", POP_MC_RANGE[1], POP_FGW_RANGE[1]),
               ("POPULATION prior best corner", POP_MC_RANGE[0], POP_FGW_RANGE[0]),
               ("TE.phi_bounds worst corner", PHI_MC_RANGE[1], PHI_FGW_RANGE[1]),
               ("TE.phi_bounds mc-max @ pop fgw-max", PHI_MC_RANGE[1], POP_FGW_RANGE[1]),
               ("TE.phi_bounds fgw-max @ pop mc-max", POP_MC_RANGE[1], PHI_FGW_RANGE[1])]
    amap = []
    for lab, a, b in corners:
        tc = float(t_coal(a, b)) / YR
        ok = tc > t_max / YR
        amap.append((lab, a, b, tc, ok))
        print(f"  {lab:<40} {a:9.2f} {b:10.2f} {tc:11.2f} "
              f"{'EVALUABLE' if ok else '*** NON-EVALUABLE':>14}")
    print(f"\n  ==> the POPULATION's generative box is evaluable EVERYWHERE "
          f"(worst corner {amap[0][3]:.1f} yr vs t_max {t_max/YR:.1f} yr = "
          f"{amap[0][3]/(t_max/YR):.1f}x margin).")
    print(f"  ==> TE.phi_bounds -- what SPARK-2 clipped to -- is NOT the prior support and "
          f"admits coalescence.", flush=True)

    cs = certified_of(key, geo_id, real_i)
    R = venue_realisation(P, G, C, IG, jnp, v["tier"], cs["noise_seed"], cs["dist_seed"],
                          k_loud=v["k"])
    npsr = P["npsr"]
    lo_b, hi_b = TE.phi_bounds(P)
    nd = P["n_dist"]

    # ---- a CLEAN (solo) re-measure of the E-step cost. gate1's absolutes were taken with a
    # co-tenant job on the node, which is an undeclared input to a wallclock number; the RATIO
    # survived that (same process, back-to-back) but the absolutes did not. Cheap to redo here.
    _one = np.ones(npsr)
    _th = R["theta_true"].copy(); _th[P["AI"]] = P["L0"]
    P["sp"].estep(_th, G["EV"], P["AI"], R["data"], _one)          # compile
    _t0 = time.time(); P["sp"].estep(_th, G["EV"], P["AI"], R["data"], _one)
    t_std_solo = time.time() - _t0
    _t0 = time.time()
    estep_per_target(P["sp"], _th, G["EV"], P["AI"], R["data"], np.eye(npsr), jnp)
    t_pt_solo = time.time() - _t0
    print(f"[SPARK3] SOLO estep timing: standard {t_std_solo:.2f}s  per-target {t_pt_solo:.2f}s"
          f"  = {t_pt_solo/max(t_std_solo,1e-9):.2f}x   (gate1 measured 0.91/3.73 = 4.08x "
          f"WITH a co-tenant on the node)", flush=True)

    # ---- g2a: the predicate vs the actual E-step ----
    print("\n  GATE g2a -- the ANALYTIC predicate vs the ACTUAL E-step", flush=True)
    one = np.ones(npsr)
    res_a = []
    for lab, lmc, lfgw in [("evaluable   (pop worst corner)", POP_MC_RANGE[1], POP_FGW_RANGE[1]),
                           ("NON-evaluable (phi worst corner)", PHI_MC_RANGE[1], PHI_FGW_RANGE[1])]:
        th = R["theta_true"].copy(); th[P["AI"]] = P["L0"]
        th[nd + 8 * v["k"] + I_MC] = lmc            # perturb ONE faint source only
        th[nd + 8 * v["k"] + I_FGW] = lfgw
        pred = bool(evaluable(lmc, lfgw, t_max))
        LNL = P["sp"].estep(th, G["EV"], P["AI"], R["data"], one)
        nbad = int((~np.isfinite(LNL)).sum())
        tot = int(LNL.size)
        agree = (pred and nbad == 0) or ((not pred) and nbad == tot)
        res_a.append((lab, pred, nbad, tot, agree))
        print(f"    {lab:<34} predicted {'EVALUABLE' if pred else 'NON-EVALUABLE':<14} "
              f"-> non-finite {nbad}/{tot}   {'AGREE' if agree else '*** DISAGREE'}",
              flush=True)
    g2a = all(r[4] for r in res_a)
    print(f"    ==> g2a {'PASSED -- one coalesced source NaNs the WHOLE array' if g2a else '*** FAILED'}",
          flush=True)

    # ---- g2b/g2c: per-rung non-evaluable fraction, SPARK-2's box vs the repair ----
    print("\n  GATE g2b/g2c -- per-rung non-evaluable fraction of the SOFT draw", flush=True)
    print(f"    {'rung':>4} {'coh':>4} | {'SPARK-2 box (phi_bounds)':>26} | "
          f"{'REPAIRED (pop box + joint gate)':>32}")
    print(f"    {'':>4} {'':>4} | {'src non-eval':>13} {'draws dead':>12} | "
          f"{'rejections':>12} {'draws dead':>12} {'E-step':>7}")
    rows = []
    for k in RUNGS:
        if cs["n_cert"] < k:
            print(f"    {k:>4} {'--':>4} |  rung unreachable (N_cert = {cs['n_cert']})")
            continue
        PM, take, pm_rung = rung_masks(npsr, cs["order"], cs["qmax"], k)
        fi = faint_fisher_bounds(P, R, pm_rung, P["L0"], jnp, v["k"])
        idx = fi["idx"]
        # --- SPARK-2's recipe: clip to phi_bounds, no evaluability test ---
        rng = np.random.default_rng(SOFT_SEED + 1000 * real_i + k)
        BND2 = (lo_b[idx - nd], hi_b[idx - nd])
        nbad_src = 0; dead = 0
        for _ in range(N_SOFT):
            d = np.clip(R["theta_true"][idx] + rng.normal(0.0, fi["sig_opt"]), *BND2)
            bad = 0
            for si in range(len(idx) // len(FAINT_PARS)):
                sl = slice(si * len(FAINT_PARS), (si + 1) * len(FAINT_PARS))
                mcv = d[sl][list(FAINT_PARS).index(I_MC)]
                fgv = d[sl][list(FAINT_PARS).index(I_FGW)]
                bad += int(not evaluable(mcv, fgv, t_max))
            nbad_src += bad
            dead += int(bad > 0)
        # --- the repair: population box + the JOINT analytic gate + reject-and-redraw ---
        lo_p = lo_b[idx - nd].copy(); hi_p = hi_b[idx - nd].copy()
        for si in range(len(idx) // len(FAINT_PARS)):
            b0 = si * len(FAINT_PARS)
            lo_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[0]
            hi_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[1]
            lo_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[0]
            hi_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[1]
        rng2 = np.random.default_rng(SOFT_SEED + 1000 * real_i + k)
        nrej = 0; dead2 = 0; nbad_e = 0
        for _ in range(N_SOFT):
            th_mo, nr = soft_faint_theta(P, R, fi["sig_opt"], idx, rng2, (lo_p, hi_p),
                                         t_max, v["k"])
            nrej += nr
            th_mo[P["AI"]] = P["L0"]
            LNL = P["sp"].estep(th_mo, G["EV"], P["AI"], R["data"], one)
            nb = int((~np.isfinite(LNL)).sum())
            nbad_e += nb
            dead2 += int(nb > 0)
        print(f"    {k:>4} {len(take):>4} | {nbad_src:>13} {dead:>5}/{N_SOFT:<6} | "
              f"{nrej:>12} {dead2:>5}/{N_SOFT:<6} {'FINITE' if nbad_e==0 else '*** NaN':>7}",
              flush=True)
        rows.append(dict(rung=k, coh=len(take), s2_bad_src=nbad_src, s2_dead=dead,
                         rep_rej=nrej, rep_dead=dead2, rep_nonfinite=nbad_e))
    g2c = all(r["rep_dead"] == 0 and r["rep_nonfinite"] == 0 for r in rows)
    g2b = any(r["s2_dead"] > 0 for r in rows)
    print(f"\n    ==> g2b {'PASSED -- SPARK-2 box reproduces the deaths' if g2b else 'no deaths reproduced under SPARK-2 box'}")
    print(f"    ==> g2c {'PASSED -- the repaired draw is finite at EVERY rung' if g2c else '*** FAILED'}",
          flush=True)

    os.makedirs(OUT, exist_ok=True)
    np.savez(f"{OUT}/spark3_gate2_{key}_g{geo_id}_r{real_i}.npz",
             venue=key, geo_id=geo_id, real_i=real_i, t_max_s=t_max, t_max_yr=t_max / YR,
             g2a_passed=g2a, g2b_reproduced=g2b, g2c_passed=g2c,
             t_standard_solo=t_std_solo, t_per_target_solo=t_pt_solo,
             cost_ratio_solo=t_pt_solo / max(t_std_solo, 1e-9),
             amap_label=np.array([a[0] for a in amap]),
             amap_mc=np.array([a[1] for a in amap]), amap_fgw=np.array([a[2] for a in amap]),
             amap_tcoal_yr=np.array([a[3] for a in amap]),
             amap_evaluable=np.array([a[4] for a in amap]),
             rung=np.array([r["rung"] for r in rows]),
             s2_bad_src=np.array([r["s2_bad_src"] for r in rows]),
             s2_dead=np.array([r["s2_dead"] for r in rows]),
             rep_rej=np.array([r["rep_rej"] for r in rows]),
             rep_dead=np.array([r["rep_dead"] for r in rows]),
             rep_nonfinite=np.array([r["rep_nonfinite"] for r in rows]),
             n_soft=N_SOFT, pop_mc=np.array(POP_MC_RANGE), pop_fgw=np.array(POP_FGW_RANGE),
             phi_mc=np.array(PHI_MC_RANGE), phi_fgw=np.array(PHI_FGW_RANGE),
             note=("FIX 2. The all-NaN E-step is an EARTH-TERM COALESCENCE pathology, exact "
                   "and analytic: term = 1 - (256/5) mc^(5/3) w0^(8/3) t < 0 for t > t_coal, "
                   "and jnp.power(term,-3/8) = NaN (deterministic.py:509-510). The delay SUMS "
                   "all 16 sources, so ONE coalesced source NaNs every pulsar -> non-finite "
                   "on ALL npsr*B entries = SPARK-2's 59392/59392 signature. ROOT CAUSE: "
                   "SPARK-2 clipped the soft draw to TE.phi_bounds (mc[7,10], fgw[-8.7,-7.0]) "
                   "-- the ESTIMATOR's SEARCH box -- believing it the prior support. The "
                   "POPULATION's generative prior is mc[8.5,9.5], fgw[-8.0,-7.5] and is "
                   "evaluable everywhere at this t_max. It IS a double perturbation, but of "
                   "(mc x fgw) JOINTLY -- not KINDLE's (scrambled x fix_mc) degenerate "
                   "covariance. EMBER's mc_scan holds fgw at truth so its 1-D scan cannot see "
                   "the joint corner: that, not 'loud vs faint', is why it measured 0% "
                   "non-evaluable. t_coal is a function of (mc,fgw) ONLY -- the strain plays "
                   "no part. Raw per-rung rejection counts banked, never verdicts."))
    print(f"[SPARK3] banked -> {OUT}/spark3_gate2_{key}_g{geo_id}_r{real_i}.npz", flush=True)
    return 0 if (g2a and g2c) else 1


def offender(dlnL, lnK, qmax):
    """RECUT recut_surface.py:75-78 verbatim -- max over pulsars, else 0.0 (the zero point
    mass). Non-finite dlnL (a pulsar with a single fringe peak, K=1) is EXCLUDED and counted:
    an inf offender would silently poison the extreme-value fit."""
    fin = np.isfinite(dlnL)
    m = (dlnL > lnK) & (qmax > QBAR) & fin
    return (float(dlnL[m].max()) if m.any() else 0.0), int((~fin).sum())


def _rung_path(key, geo_id, real_i, k):
    return f"{OUT}/s3r_{key}_g{geo_id}_r{real_i}_k{k}.npz"


def mode_unit(key, geo_id, real_i, n_null, only_rung=None, self_lane=False):
    """PER-RUNG CHECKPOINTING, and why it is not an optimisation.

    A unit is 4 rungs x ~2.3 h under the co-tenancy this grid actually runs at (the per-target
    E-step is 3.62 s solo but 20.0 s with 4 jobs sharing the node -- MEASURED, see the report's
    co-tenancy retraction). Skip-on-exist at UNIT granularity therefore means a unit that dies
    at rung 8 loses rungs 0/2/5 too -- hours of finished work thrown away for nothing. Each
    rung is banked the moment it completes and skipped if already on disk; the unit-level npz
    is assembled from the rung files. Resume is then free at the granularity the work actually
    has.
    """
    """ONE realisation, ALL rungs IN-PROCESS (SPARK-2 5(2): one build, all rungs -- the cache
    lesson; SPARK-2 paid a cold build on each of 12 one-shot units).

    THE LADDER, per rung (SPARK-2 5(3): without the truth rung the ledger has no scale):
        UNMODELLED   the reservoir is absent from the model (log10_h -> H_ABSENT)
        SOFT @ OPT   soft-modelled at the OPTIMISTIC (conditional) widths     ] FIX 3's
        SOFT @ PES   soft-modelled at the PESSIMISTIC (prior) widths          ] two bounds
        TRUTH        the reservoir modelled at truth -- the ceiling, and the control
    SPARK-2's endpoints (UNMODELLED 0.000 vs TRUTH 4.435 at ITS cell) are the frame; the soft
    states are the honest middle this campaign exists to measure.

    Certification uses the PER-TARGET E-step (FIX 1) at every rung: the target's own term is
    always live, so its fringe row is never flat and the rung-0 count is not structurally 0.
    """
    os.makedirs(OUT, exist_ok=True)
    path = f"{OUT}/s3_{key}_g{geo_id}_r{real_i}.npz"
    if only_rung is None and os.path.exists(path):
        print(f"[SPARK3] skip-on-exist: {path}", flush=True)
        return 0
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    t0 = time.time()
    P, G, v = build_venue(key, geo_id)
    print(f"[SPARK3] build {time.time()-t0:.0f}s", flush=True)
    t_max = t_max_of(P)
    npsr, nd = P["npsr"], P["n_dist"]
    FT, EV, dL = G["FT"], G["EV"], G["dL"]
    L0 = np.asarray(P["L0"], float)
    cs = certified_of(key, geo_id, real_i, self_lane)
    print(f"[SPARK3] unit venue {key} g{geo_id} real {real_i}: noise {cs['noise_seed']} "
          f"dist {cs['dist_seed']}  N_cert(own) {cs['n_cert']}  floor {cs['floor']:.3f}",
          flush=True)
    R = venue_realisation(P, G, C, IG, jnp, v["tier"], cs["noise_seed"], cs["dist_seed"],
                          k_loud=v["k"])
    Lt = np.asarray(R["L_true"], float)
    lo_b, hi_b = TE.phi_bounds(P)

    def score_state(th, PM, data, L_true):
        LNL = estep_per_target(P["sp"], th, EV, P["AI"], data, PM, jnp)
        return score_from_LNL(P, FT, LNL, C, L_true, dL)

    out_rows = []
    for k in (RUNGS if only_rung is None else [only_rung]):
        if cs["n_cert"] < k:
            print(f"[SPARK3]   rung {k}: UNREACHABLE (this realisation's N_cert = "
                  f"{cs['n_cert']}) -- skipped and BANKED as such", flush=True)
            out_rows.append(dict(rung=k, reachable=False))
            continue
        rp = _rung_path(key, geo_id, real_i, k)
        if os.path.exists(rp):
            print(f"[SPARK3]   rung {k}: skip-on-exist ({os.path.basename(rp)})", flush=True)
            z = np.load(rp, allow_pickle=True)
            out_rows.append({kk: z[kk] for kk in z.files} | dict(rung=k, reachable=True))
            continue
        tk = time.time()
        PM, take, pm_rung = rung_masks(npsr, cs["order"], cs["qmax"], k)
        # ---- (a) the two Fisher bounds at this rung, at the ARRAY's state (pm_rung), not at
        # any target's row: the faint constraint does not depend on who is being scored.
        fi = faint_fisher_bounds(P, R, pm_rung, L0, jnp, v["k"])
        idx = fi["idx"]
        # the soft draw's support: the POPULATION's generative box on (mc,fgw) -- the true
        # prior support -- and phi_bounds elsewhere; the JOINT evaluability gate on top.
        lo_p = lo_b[idx - nd].copy(); hi_p = hi_b[idx - nd].copy()
        for si in range(len(idx) // len(FAINT_PARS)):
            b0 = si * len(FAINT_PARS)
            lo_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[0]
            hi_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[1]
            lo_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[0]
            hi_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[1]

        # ---- the four reservoir states ----
        th_un = R["theta_true"].copy(); th_un[P["AI"]] = L0
        for i in range(v["k"], P["ncw"]):
            th_un[nd + 8 * i + I_H] = H_ABSENT
        th_tr = R["theta_true"].copy(); th_tr[P["AI"]] = L0

        s_un = score_state(th_un, PM, R["data"], Lt)
        s_tr = score_state(th_tr, PM, R["data"], Lt)
        soft = {}
        nrej = {}
        for bnd, sig in (("opt", fi["sig_opt"]), ("pes", fi["sig_pes"])):
            rng = np.random.default_rng(SOFT_SEED + 1000 * real_i + 10 * k
                                        + (0 if bnd == "opt" else 5))
            acc, rj = [], 0
            for _ in range(N_SOFT):
                th_mo, nr = soft_faint_theta(P, R, sig, idx, rng, (lo_p, hi_p), t_max, v["k"])
                rj += nr
                th_mo[P["AI"]] = L0
                acc.append(score_state(th_mo, PM, R["data"], Lt))
            # SOFT marginalisation: average the per-pulsar statistic over the faint posterior
            soft[bnd] = dict(dlnL=np.mean([x["dlnL"] for x in acc], axis=0),
                             lnK=acc[0]["lnK"],
                             qmax=np.mean([x["qmax"] for x in acc], axis=0),
                             on_true=acc[0]["on_true"])
            nrej[bnd] = rj

        states = dict(un=s_un, opt=soft["opt"], pes=soft["pes"], tr=s_tr)
        offs = {n: offender(s["dlnL"], s["lnK"], s["qmax"])[0] for n, s in states.items()}
        print(f"[SPARK3]   rung {k} (coh {len(take)}): offender  un {offs['un']:.3f}  "
              f"soft@opt {offs['opt']:.3f}  soft@pes {offs['pes']:.3f}  truth {offs['tr']:.3f}"
              f"   [rej opt {nrej['opt']} pes {nrej['pes']}]", flush=True)

        # ---- nulls -> the certification floor per state ----
        NUL = {n: [] for n in states}
        for i in range(n_null):
            sd = cs["noise_seed"] + 500_000 + i
            Ltn, _ = IG.draw_true_distances_tier(P, dL, EV, seed=cs["dist_seed"] + 500_000 + i,
                                                 tier=v["tier"])
            noise = P["nd"].draw(sd)
            dn = tuple(jnp.asarray(np.asarray(b)) for b in noise)   # nullN: NO CW in the data
            NUL["un"].append(offender(*[score_state(th_un, PM, dn, Ltn)[q]
                                        for q in ("dlnL", "lnK", "qmax")])[0])
            NUL["tr"].append(offender(*[score_state(th_tr, PM, dn, Ltn)[q]
                                        for q in ("dlnL", "lnK", "qmax")])[0])
            for bnd, sig in (("opt", fi["sig_opt"]), ("pes", fi["sig_pes"])):
                rngn = np.random.default_rng(SOFT_SEED + 77_000 + 1000 * real_i + 10 * k + i
                                             + (0 if bnd == "opt" else 5))
                thm, _ = soft_faint_theta(P, R, sig, idx, rngn, (lo_p, hi_p), t_max, v["k"])
                thm[P["AI"]] = L0
                s = score_state(thm, PM, dn, Ltn)
                NUL[bnd].append(offender(s["dlnL"], s["lnK"], s["qmax"])[0])
            if (i + 1) % 20 == 0:
                print(f"      null {i+1}/{n_null}", flush=True)

        from spark import adopt
        FL = {n: adopt(np.array(NUL[n])) for n in states}
        row = dict(rung=k, reachable=True, coh=len(take), take=take)
        for n in states:
            bar = np.maximum(states[n]["lnK"], FL[n]["floor"])
            cert = (states[n]["dlnL"] > bar) & (states[n]["qmax"] > QBAR)
            row[f"ncert_{n}"] = int(cert.sum())
            row[f"margin_{n}"] = states[n]["dlnL"] - bar
            row[f"floor_{n}"] = FL[n]["floor"]; row[f"zf_{n}"] = FL[n]["zf"]
            row[f"est_{n}"] = FL[n]["estimator"]; row[f"ferr_{n}"] = FL[n]["err"]
            row[f"off_{n}"] = offs[n]
            row[f"dlnL_{n}"] = states[n]["dlnL"]; row[f"qmax_{n}"] = states[n]["qmax"]
            row[f"null_{n}"] = np.array(NUL[n])
            print(f"[SPARK3]     {n:>3}: floor {FL[n]['floor']:8.3f} +- {FL[n]['err']:.3f} "
                  f"[{FL[n]['estimator']}] zf {FL[n]['zf']:.2f}  N_cert {row[f'ncert_{n}']}",
                  flush=True)
        row["lnK"] = states["un"]["lnK"]
        row["sig_opt"] = fi["sig_opt"]; row["sig_pes"] = fi["sig_pes"]
        row["cond_sig"] = fi["cond_sig"]; row["F_ii"] = fi["F_ii"]
        row["n_pos"] = fi["n_pos"]; row["cond_num"] = fi["cond_num"]
        row["nrej_opt"] = nrej["opt"]; row["nrej_pes"] = nrej["pes"]
        np.savez(rp, **{kk: vv for kk, vv in row.items() if kk != "reachable"},
                 note=("SPARK-3 per-RUNG checkpoint. Banked the moment the rung completes so a "
                       "unit that dies keeps its finished rungs (a rung is ~2.3 h under the "
                       "grid's measured co-tenancy). Raw per-pulsar dlnL/lnK/qmax + raw null "
                       "offenders per state; both Fisher bounds as columns; never verdicts."))
        out_rows.append(row)
        print(f"[SPARK3]   rung {k} done in {time.time()-tk:.0f}s -> {os.path.basename(rp)}",
              flush=True)

    blob = dict(venue=key, geo_id=geo_id, real_i=real_i, h=v["h"], T_label=v["T"],
                tier=v["tier"], k_loud=v["k"], noise_seed=cs["noise_seed"],
                dist_seed=cs["dist_seed"], n_cert_own=cs["n_cert"], venue_floor=cs["floor"],
                cert_order=cs["order"], cert_qmax=cs["qmax"], cert_on_true=cs["on_true"],
                n_null=n_null, n_soft=N_SOFT, t_max_yr=t_max / YR,
                rungs=np.array([r["rung"] for r in out_rows]),
                reachable=np.array([r["reachable"] for r in out_rows]),
                note=("SPARK-3 arrow-2 unit at a RESERVED above-onset venue. RAW per-pulsar "
                      "dlnL/lnK/qmax + raw null offender vectors banked per (rung, state), "
                      "never verdicts. States: un = reservoir UNMODELLED; opt/pes = SOFT "
                      "Monte-Carlo marginalisation over N(truth, sigma) at FIX 3's OPTIMISTIC "
                      "(conditional) and PESSIMISTIC (prior-width) bounds -- BOTH banked as "
                      "columns, never merged; tr = TRUTH-modelled control (SPARK-2 5(3): "
                      "without it the ladder has no scale). E-step is PER-TARGET (FIX 1): the "
                      "target's own term is always live, so no fringe row is structurally "
                      "flat. Soft draw is bounded by the POPULATION's generative support on "
                      "(mc,fgw) + the JOINT analytic evaluability gate; rejection counts "
                      "banked (FIX 2). Nulls: nullN (no CW in the data), model unchanged."))
    for r in out_rows:
        if not r["reachable"]:
            continue
        for kk, vv in r.items():
            if kk in ("rung", "reachable"):
                continue
            blob[f"k{r['rung']}_{kk}"] = vv
    if only_rung is not None:
        print(f"[SPARK3] rung-split job: rung {only_rung} banked; the unit npz is assembled by "
              f"`ledger`, which reads the rung checkpoints.", flush=True)
        return 0
    np.savez(path, **blob)
    print(f"[SPARK3] banked -> {path}", flush=True)
    return 0


def marg_width_from_hess(H, box):
    """The DEFENSIBLE marginal width for a possibly-INDEFINITE faint-block Hessian.

    The chunked-JVP reveals the faint Fisher F = -H is indefinite at truth (weak sources ->
    truth is a saddle, noise-dominated on ~half the directions), so sigma_marg = sqrt(diag(inv
    F)) is not well-defined -- inverting an indefinite matrix gives negative variances. The
    posterior precision F + P (P = prior precision = 1/box^2) is the right object, but it too
    is non-PD where the likelihood RUNS AWAY (generalized eig lam <= -1): there the only bound
    is the hard prior box. So: whiten by the prior, keep the data-CONSTRAINING directions
    (lam > 0), cap the rest at the prior. This is the constrained-subspace marginal and is
    well-defined for any H. Reduces to the ordinary marginal when F is PD, and to the prior on
    an unconstrained/runaway block.
    """
    from scipy.linalg import eigh as _eigh
    F = -np.asarray(H, float)
    Pi2 = np.diag(1.0 / box)                      # prior^{-1/2} whitening
    M = Pi2 @ F @ Pi2                              # whitened data precision (eig = generalized lam)
    w, V = _eigh(0.5 * (M + M.T))
    post_prec = 1.0 + np.where(w > 0, w, 0.0)     # keep data directions; runaway/flat -> prior (=1)
    Cov = Pi2 @ (V @ np.diag(1.0 / post_prec) @ V.T) @ Pi2
    return np.sqrt(np.clip(np.diag(Cov), 0.0, None))


def mode_remarg(key, geo_id, real_i, rung, n_null, self_lane=True):
    """RE-SCORE the crossing test at the MARGINAL width (the split of the STRADDLE).

    The JVP (mode `jvp`) banked the faint-block Hessian; `marg_width_from_hess` turns it into
    the defensible marginal width sigma_marg (median ~3x the conditional at the venue -- BETWEEN
    the optimistic and pessimistic bounds, so the width alone does not decide the crossings).
    This re-scores the SOFT state at sigma_marg -- one more soft state beside opt/pes -- and
    banks margin_marg, so the reader can re-cut crossings (uncertified@rung0 -> certified@rung8)
    at the true marginal. Floor re-cut at matched state, as everywhere.
    """
    os.makedirs(OUT, exist_ok=True)
    path = f"{OUT}/s3marg_{key}_g{geo_id}_r{real_i}_k{rung}.npz"
    if os.path.exists(path):
        print(f"[SPARK3] skip-on-exist: {path}", flush=True)
        return 0
    jvpf = f"{OUT}/s3jvp_{key}_g{geo_id}_r{real_i}_k{rung}.npz"
    if not os.path.exists(jvpf):
        print(f"*** no JVP Hessian for r{real_i} k{rung} -- run mode jvp first", flush=True)
        return 2
    zj = np.load(jvpf, allow_pickle=True)
    if bool(zj["capped"]):
        print(f"*** JVP r{real_i} k{rung} was CAPPED (straddle-unresolved) -- no marginal width",
              flush=True)
        return 2
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    P, G, v = build_venue(key, geo_id)
    EV, dL = G["EV"], G["dL"]; FT = G["FT"]
    L0 = np.asarray(P["L0"], float); nd = P["n_dist"]; npsr = P["npsr"]
    cs = certified_of(key, geo_id, real_i, self_lane)
    R = venue_realisation(P, G, C, IG, jnp, v["tier"], cs["noise_seed"], cs["dist_seed"],
                          k_loud=v["k"])
    Lt = np.asarray(R["L_true"], float)
    PM, take, pm_rung = rung_masks(npsr, cs["order"], cs["qmax"], rung)
    fi = faint_fisher_bounds(P, R, pm_rung, L0, jnp, v["k"])
    idx = fi["idx"]
    box = fi["box"]
    sig_marg = marg_width_from_hess(np.asarray(zj["hess"], float), box)
    print(f"[SPARK3] remarg r{real_i} k{rung}: median sigma cond {np.median(fi['cond_sig']):.4g}"
          f" -> MARG {np.median(sig_marg):.4g} -> prior {np.median(box):.4g}", flush=True)
    lo_b, hi_b = TE.phi_bounds(P)
    t_max = t_max_of(P)
    lo_p = lo_b[idx - nd].copy(); hi_p = hi_b[idx - nd].copy()
    for si in range(len(idx) // len(FAINT_PARS)):
        b0 = si * len(FAINT_PARS)
        lo_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[0]
        hi_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[1]
        lo_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[0]
        hi_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[1]

    def score_state(th, data, L_true):
        return score_from_LNL(P, FT, estep_per_target(P["sp"], th, EV, P["AI"], data, PM, jnp),
                              C, L_true, dL)

    rng = np.random.default_rng(SOFT_SEED + 55_000 + 1000 * real_i + 10 * rung)
    acc, nrej = [], 0
    for _ in range(N_SOFT):
        th_mo, nr = soft_faint_theta(P, R, sig_marg, idx, rng, (lo_p, hi_p), t_max, v["k"])
        nrej += nr; th_mo[P["AI"]] = L0
        acc.append(score_state(th_mo, R["data"], Lt))
    s_marg = dict(dlnL=np.mean([x["dlnL"] for x in acc], axis=0), lnK=acc[0]["lnK"],
                  qmax=np.mean([x["qmax"] for x in acc], axis=0))
    NUL = []
    for i in range(n_null):
        Ltn, _ = IG.draw_true_distances_tier(P, dL, EV, seed=cs["dist_seed"] + 800_000 + i,
                                             tier=v["tier"])
        noise = P["nd"].draw(cs["noise_seed"] + 800_000 + i)
        dn = tuple(jnp.asarray(np.asarray(b)) for b in noise)
        rngn = np.random.default_rng(SOFT_SEED + 66_000 + 1000 * real_i + 10 * rung + i)
        thm, _ = soft_faint_theta(P, R, sig_marg, idx, rngn, (lo_p, hi_p), t_max, v["k"])
        thm[P["AI"]] = L0
        s = score_state(thm, dn, Ltn)
        NUL.append(offender(s["dlnL"], s["lnK"], s["qmax"])[0])
        if (i + 1) % 25 == 0:
            print(f"   null {i+1}/{n_null}", flush=True)
    from spark import adopt
    FL = adopt(np.array(NUL))
    bar = np.maximum(s_marg["lnK"], FL["floor"])
    cert = (s_marg["dlnL"] > bar) & (s_marg["qmax"] > QBAR)
    margin_marg = s_marg["dlnL"] - bar
    print(f"[SPARK3] remarg r{real_i} k{rung}: floor {FL['floor']:.3f} [{FL['estimator']}] "
          f"zf {FL['zf']:.2f}  N_cert(marg) {int(cert.sum())}", flush=True)
    np.savez(path, venue=key, real_i=real_i, rung=rung, sig_marg=sig_marg,
             median_sig_marg=float(np.median(sig_marg)),
             median_sig_cond=float(np.median(fi["cond_sig"])),
             dlnL_marg=s_marg["dlnL"], qmax_marg=s_marg["qmax"], lnK=s_marg["lnK"],
             floor_marg=FL["floor"], zf_marg=FL["zf"], est_marg=FL["estimator"],
             margin_marg=margin_marg, ncert_marg=int(cert.sum()), null_marg=np.array(NUL),
             nrej=nrej, n_null=n_null, n_soft=N_SOFT,
             note=("SPARK-3 STRADDLE SPLIT: soft state re-scored at the MARGINAL width "
                   "sigma_marg = marg_width_from_hess(JVP Hessian) -- the defensible "
                   "constrained-subspace marginal for the indefinite faint Fisher. margin_marg "
                   "is the certification margin under the TRUE marginal; the reader re-cuts "
                   "crossings (margin_marg<0 @rung0 -> >0 @rung8) from it. Floor at matched "
                   "state. Raw columns banked, never a verdict."))
    print(f"[SPARK3] banked -> {path}", flush=True)
    return 0


def mode_jvp(key, reals, rungs=(0, 8), ch=8, cap_s=2700, self_lane=True):
    """ARM (c): THE BUDGETED CHUNKED-JVP — the true MARGINAL width, to split the STRADDLE.

    STRADDLED (§5) is the disagreement between the OPTIMISTIC (conditional, 1/sqrt(F_ii)) and
    PESSIMISTIC (prior) Fisher bounds. The verdict between them needs the MARGINAL width --
    sigma_marg = sqrt(diag(inv(H))) where H is the full faint-block Hessian of lnL -- which
    SPARK-2 measured is not a few-GPU-hr object (monolithic jax.hessian >27 min; chunked
    batched-JVP-of-grad at CH=8 >17 min, both CANCELLED unfinished on an A100-80GB). This mode
    does exactly the cheaper of those two, budgeted:

      f(x) = lnL(theta with the faint block set to x, data, pm_rung)    x in R^{78}
      g    = grad f                                                     (one GP-solve per eval)
      H    = jac g   computed COLUMN-CHUNKED: H @ E_chunk = jvp(g, x, E_chunk), CH columns/JVP
      sigma_marg[i] = sqrt( pinv(H)[i,i] )

    THE 45-MIN HARD CAP (cap_s), enforced PER HESSIAN inside the process. The chunk loop times
    every JVP and BAILS the moment cumulative Hessian time crosses the cap -- including the
    first chunk, which pays compilation and may alone exceed it (SPARK-2's CH=8 was still
    running at 17 min). A capped Hessian banks `capped=True` and sigma_marg=NaN: the unit is
    reported straddle-UNRESOLVED, never silently dropped. "Unaffordable even chunked" over
    enough units is an AVALANCHE-shaped closure and an ACCEPTABLE outcome (Matt's decision).

    GROUPED per process (SPARK-2's cache lesson): one build, then all (real x rung) Hessians in
    `reals x rungs` looped in-process -- no per-Hessian cold compile. The GRADIENT graph
    compiles once and is reused across chunks and across Hessians (same shapes).
    """
    import jax
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    v = VENUES[key]
    geo_id = 3
    t_build = time.time()
    P, G, _ = build_venue(key, geo_id)
    EV, dL = G["EV"], G["dL"]
    L0 = np.asarray(P["L0"], float)
    nd = P["n_dist"]
    lb = P["amo"]["logL"]                       # scalar lnL(theta, data, pmask)
    print(f"[SPARK3] JVP build {time.time()-t_build:.0f}s; venue {key} reals {list(reals)} "
          f"rungs {list(rungs)} CH={ch} cap={cap_s}s", flush=True)
    os.makedirs(OUT, exist_ok=True)

    for real_i in reals:
        cs = certified_of(key, geo_id, real_i, self_lane)
        R = venue_realisation(P, G, C, IG, jnp, v["tier"], cs["noise_seed"], cs["dist_seed"],
                              k_loud=v["k"])
        data = R["data"]
        for rung in rungs:
            path = f"{OUT}/s3jvp_{key}_g{geo_id}_r{real_i}_k{rung}.npz"
            if os.path.exists(path):
                print(f"[SPARK3]   skip-on-exist {os.path.basename(path)}", flush=True)
                continue
            if cs["n_cert"] < rung:
                print(f"[SPARK3]   r{real_i} k{rung}: unreachable (N_cert {cs['n_cert']})",
                      flush=True)
                continue
            PM, take, pm_rung = rung_masks(P["npsr"], cs["order"], cs["qmax"], rung)
            fi = faint_fisher_bounds(P, R, pm_rung, L0, jnp, v["k"])
            idx = fi["idx"]; n = len(idx)
            th0 = R["theta_true"].copy(); th0[P["AI"]] = L0
            th0j = jnp.asarray(th0); pm = jnp.asarray(pm_rung); idxj = jnp.asarray(idx)
            x0 = jnp.asarray(th0[idx])

            def f(x):                            # lnL as a function of the faint block only
                th = th0j.at[idxj].set(x)
                return lb(th, data, pm)
            g = jax.jit(jax.grad(f))
            g(x0)                                # warm the gradient graph once (shared)

            # ---- Hessian, column-chunked, with the per-Hessian wall cap ----
            E = np.eye(n)
            cols = []; capped = False; th = time.time()
            for c0 in range(0, n, ch):
                Ec = jnp.asarray(E[c0:c0 + ch])          # (cch, n) basis directions
                # H @ e for each e in the chunk = jvp of the gradient in that direction;
                # vmap over the chunk's directions -> cch rows of H at once.
                def hvp(e):
                    return jax.jvp(g, (x0,), (e,))[1]
                Hc = np.asarray(jax.vmap(hvp)(Ec))       # (cch, n) rows of H
                cols.append(Hc)
                if time.time() - th > cap_s:
                    capped = True
                    print(f"[SPARK3]   r{real_i} k{rung}: CAP hit at col {c0+len(Hc)}/{n} "
                          f"({time.time()-th:.0f}s) -> straddle-UNRESOLVED", flush=True)
                    break
            t_hess = time.time() - th
            if capped:
                np.savez(path, venue=key, real_i=real_i, rung=rung, capped=True,
                         idx=idx, sig_cond=fi["cond_sig"], sig_opt=fi["sig_opt"],
                         sig_pes=fi["sig_pes"], sig_marg=np.full(n, np.nan),
                         t_hess=t_hess, ncols_done=int(sum(len(c) for c in cols)),
                         note=("ARM (c) chunked-JVP marginal Hessian, CAPPED at the 45-min "
                               "per-Hessian wall limit -> unit is straddle-UNRESOLVED. "
                               "sigma_marg=NaN. 'Unaffordable even chunked' is a pre-registered "
                               "acceptable outcome (AVALANCHE-shaped closure)."))
                print(f"[SPARK3]   banked (capped) -> {os.path.basename(path)}", flush=True)
                continue
            H = np.concatenate(cols, axis=0)             # (n, n)
            H = 0.5 * (H + H.T)                           # symmetrise (AD round-off)
            # Fisher = -H (H is the Hessian of lnL; curvature is negative-definite at a peak).
            Fmat = -H
            # marginal variance = diag(inv(Fisher)); guard non-PD with a pinv + floor
            try:
                Cov = np.linalg.inv(Fmat)
            except np.linalg.LinAlgError:
                Cov = np.linalg.pinv(Fmat)
            var = np.diag(Cov)
            sig_marg = np.where(var > 0, np.sqrt(var), np.inf)
            box = fi["box"]
            sig_marg = np.where(np.isfinite(sig_marg), np.minimum(sig_marg, box), box)
            npd = int((np.linalg.eigvalsh(Fmat) <= 0).sum())
            print(f"[SPARK3]   r{real_i} k{rung}: Hessian {t_hess:.0f}s  "
                  f"median sigma cond {np.median(fi['cond_sig']):.4g} -> marg "
                  f"{np.median(sig_marg):.4g} (pes {np.median(fi['sig_pes']):.4g})  "
                  f"non-PD eig {npd}/{n}", flush=True)
            np.savez(path, venue=key, real_i=real_i, rung=rung, capped=False, idx=idx,
                     n_coherent=len(take), noise_seed=cs["noise_seed"],
                     dist_seed=cs["dist_seed"], F_ii=fi["F_ii"],
                     sig_cond=fi["cond_sig"], sig_opt=fi["sig_opt"], sig_pes=fi["sig_pes"],
                     sig_marg=sig_marg, hess=H, n_nonpd=npd, t_hess=t_hess,
                     note=("ARM (c) chunked-JVP MARGINAL width. sigma_marg = sqrt(diag(inv(-H))) "
                           "for the 78x78 faint-block Hessian H of lnL, computed column-chunked "
                           "at CH=%d via batched JVP-of-grad, one build per process (SPARK-2 "
                           "cache lesson). This is the width BETWEEN the campaign's optimistic "
                           "(conditional, sig_cond) and pessimistic (prior, sig_pes) bounds -- "
                           "the object that splits the STRADDLE. Banked beside both bounds as "
                           "columns. cond <= marg <= prior expected." % ch))
            print(f"[SPARK3]   banked -> {os.path.basename(path)}", flush=True)
    return 0


def mode_scram(key, geo_id, real_i, rung, n_null, mode="L", self_lane=False):
    """THE SCRAMBLED-COUNTERPART ARM — the manufacturing control. Pre-registered by the brief:
    *"scrambled-counterpart sanity at one rung (the soft reservoir must not manufacture — STOP +
    anatomy if it does)"*.

    THE QUESTION, AND WHY NOTHING ABOUT SAFETY MAY BE SAID WITHOUT IT. The soft state models the
    faint reservoir at N(truth, sigma) and AVERAGES the per-pulsar statistic over N_SOFT draws.
    Averaging a statistic over a family of models is not free: it can RAISE the statistic simply
    by giving the fit more ways to accommodate noise. If that inflation is large enough to push
    pulsars over the bar, arrow 2's ladder would be measuring the soft model's flexibility, not
    deconfusion — and every rung difference in §4 would be an artefact.

    THE CONTROL, following SURFACE's own null recipe (`surface.py:176` `scramble_s`, reused not
    reimplemented). The DATA keeps the true sky; the RECOVERY model's loud counterpart is
    sky-scrambled (`mode="L"`, `n_scr = k_loud`) or all 16 are (`mode="A"`). **No signal sits at
    the assumed counterpart, so EVERY certification here is FALSE by construction** — this is
    the same arm EMBER used to find the manufacturing boundary (S8.5.3: motion under a wrong
    counterpart, sensitivity 1.00, p = 0.002).

    THE PRE-REGISTERED READOUT, stated before it runs:
      * `N_false(soft@opt)` and `N_false(soft@pes)` vs `N_false(unmodelled)` AT THE SAME RUNG.
      * **MANUFACTURES** iff soft-modelling RAISES the false-cert count above unmodelled. That
        fires the brief's STOP, and the anatomy — which pulsars, at what margin, under which
        bound — is banked for the report rather than summarised.
      * **CLEAN** iff it does not. Note the asymmetry: this arm can only ever license the
        statement *"the soft reservoir does not manufacture in this cell at this rung"*. It is
        one rung of one realisation at one venue, and EMBER's S8.5.2 warning applies verbatim —
        a null here is not evidence of safety elsewhere, and it must never be quoted as one.

    Floors are re-cut at MATCHED state (the scrambled model's own nulls), never inherited from
    the signal arm: a bar computed under a different model is not this arm's bar.
    """
    os.makedirs(OUT, exist_ok=True)
    path = f"{OUT}/s3scr_{key}_g{geo_id}_r{real_i}_k{rung}_{mode}.npz"
    if os.path.exists(path):
        print(f"[SPARK3] skip-on-exist: {path}", flush=True)
        return 0
    jax, jnp, C, TE, IG, F, SF, _ = _stack()
    P, G, v = build_venue(key, geo_id)
    t_max = t_max_of(P)
    npsr, nd = P["npsr"], P["n_dist"]
    FT, EV, dL = G["FT"], G["EV"], G["dL"]
    L0 = np.asarray(P["L0"], float)
    cs = certified_of(key, geo_id, real_i, self_lane)
    if cs["n_cert"] < rung:
        print(f"*** rung {rung} unreachable (N_cert = {cs['n_cert']})", flush=True)
        return 2
    R = venue_realisation(P, G, C, IG, jnp, v["tier"], cs["noise_seed"], cs["dist_seed"],
                          k_loud=v["k"])
    Lt = np.asarray(R["L_true"], float)
    lo_b, hi_b = TE.phi_bounds(P)

    # ---- the WRONG counterpart: SURFACE's recipe, verbatim (surface.py:219-222) ----
    n_scr = v["k"] if mode == "L" else P["ncw"]
    theta_rec = SF.scramble_s(G["theta"], nd, P["ncw"], cs["noise_seed"] + 13, n_scr)
    print(f"[SPARK3] SCRAMBLED arm: venue {key} g{geo_id} real {real_i} rung {rung} "
          f"mode null{mode} (n_scr={n_scr}); every cert here is FALSE by construction",
          flush=True)

    PM, take, pm_rung = rung_masks(npsr, cs["order"], cs["qmax"], rung)
    fi = faint_fisher_bounds(P, R, pm_rung, L0, jnp, v["k"])
    idx = fi["idx"]
    lo_p = lo_b[idx - nd].copy(); hi_p = hi_b[idx - nd].copy()
    for si in range(len(idx) // len(FAINT_PARS)):
        b0 = si * len(FAINT_PARS)
        lo_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[0]
        hi_p[b0 + list(FAINT_PARS).index(I_MC)] = POP_MC_RANGE[1]
        lo_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[0]
        hi_p[b0 + list(FAINT_PARS).index(I_FGW)] = POP_FGW_RANGE[1]

    def score_state(th, data, L_true):
        return score_from_LNL(P, FT, estep_per_target(P["sp"], th, EV, P["AI"], data, PM, jnp),
                              C, L_true, dL)

    # the four reservoir states, all on the SCRAMBLED recovery model
    th_un = theta_rec.copy(); th_un[P["AI"]] = L0
    for i in range(v["k"], P["ncw"]):
        th_un[nd + 8 * i + I_H] = H_ABSENT
    th_tr = theta_rec.copy(); th_tr[P["AI"]] = L0
    s_un = score_state(th_un, R["data"], Lt)
    s_tr = score_state(th_tr, R["data"], Lt)
    soft, nrej = {}, {}
    for bnd, sig in (("opt", fi["sig_opt"]), ("pes", fi["sig_pes"])):
        rng = np.random.default_rng(SOFT_SEED + 31_000 + 1000 * real_i + 10 * rung
                                    + (0 if bnd == "opt" else 5))
        acc, rj = [], 0
        for _ in range(N_SOFT):
            th_mo, nr = soft_faint_theta(P, R, sig, idx, rng, (lo_p, hi_p), t_max, v["k"])
            # the faint draw is soft; the COUNTERPART stays wrong -- that is the whole point
            th_mo[:nd] = theta_rec[:nd]
            for j in range(v["k"]):
                th_mo[nd + 8 * j: nd + 8 * (j + 1)] = theta_rec[nd + 8 * j: nd + 8 * (j + 1)]
            rj += nr
            th_mo[P["AI"]] = L0
            acc.append(score_state(th_mo, R["data"], Lt))
        soft[bnd] = dict(dlnL=np.mean([x["dlnL"] for x in acc], axis=0), lnK=acc[0]["lnK"],
                         qmax=np.mean([x["qmax"] for x in acc], axis=0),
                         on_true=acc[0]["on_true"])
        nrej[bnd] = rj

    states = dict(un=s_un, opt=soft["opt"], pes=soft["pes"], tr=s_tr)
    # ---- matched nulls -> this arm's OWN floor per state ----
    NUL = {n: [] for n in states}
    for i in range(n_null):
        Ltn, _ = IG.draw_true_distances_tier(P, dL, EV, seed=cs["dist_seed"] + 700_000 + i,
                                             tier=v["tier"])
        noise = P["nd"].draw(cs["noise_seed"] + 700_000 + i)
        dn = tuple(jnp.asarray(np.asarray(b)) for b in noise)
        for n, th in (("un", th_un), ("tr", th_tr)):
            s = score_state(th, dn, Ltn)
            NUL[n].append(offender(s["dlnL"], s["lnK"], s["qmax"])[0])
        for bnd, sig in (("opt", fi["sig_opt"]), ("pes", fi["sig_pes"])):
            rngn = np.random.default_rng(SOFT_SEED + 91_000 + 1000 * real_i + 10 * rung + i
                                         + (0 if bnd == "opt" else 5))
            thm, _ = soft_faint_theta(P, R, sig, idx, rngn, (lo_p, hi_p), t_max, v["k"])
            thm[:nd] = theta_rec[:nd]
            for j in range(v["k"]):
                thm[nd + 8 * j: nd + 8 * (j + 1)] = theta_rec[nd + 8 * j: nd + 8 * (j + 1)]
            thm[P["AI"]] = L0
            s = score_state(thm, dn, Ltn)
            NUL[bnd].append(offender(s["dlnL"], s["lnK"], s["qmax"])[0])
        if (i + 1) % 20 == 0:
            print(f"   null {i+1}/{n_null}", flush=True)

    from spark import adopt
    FL = {n: adopt(np.array(NUL[n])) for n in states}
    blob = dict(venue=key, geo_id=geo_id, real_i=real_i, rung=rung, mode=mode, n_scr=n_scr,
                coh=len(take), noise_seed=cs["noise_seed"], dist_seed=cs["dist_seed"],
                n_null=n_null, n_soft=N_SOFT, sig_opt=fi["sig_opt"], sig_pes=fi["sig_pes"],
                nrej_opt=nrej["opt"], nrej_pes=nrej["pes"])
    print(f"\n  {'state':>5} {'floor':>9} {'zf':>5} {'est':>8} {'N_FALSE':>8}")
    nf = {}
    for n in states:
        bar = np.maximum(states[n]["lnK"], FL[n]["floor"])
        cert = (states[n]["dlnL"] > bar) & (states[n]["qmax"] > QBAR)
        nf[n] = int(cert.sum())          # scrambled => EVERY cert is false by construction
        blob[f"nfalse_{n}"] = nf[n]
        blob[f"margin_{n}"] = states[n]["dlnL"] - bar
        blob[f"dlnL_{n}"] = states[n]["dlnL"]; blob[f"qmax_{n}"] = states[n]["qmax"]
        blob[f"null_{n}"] = np.array(NUL[n]); blob[f"floor_{n}"] = FL[n]["floor"]
        blob[f"zf_{n}"] = FL[n]["zf"]; blob[f"est_{n}"] = FL[n]["estimator"]
        print(f"  {n:>5} {FL[n]['floor']:9.3f} {FL[n]['zf']:5.2f} {FL[n]['estimator']:>8} "
              f"{nf[n]:8d}")
    blob["lnK"] = states["un"]["lnK"]
    manufactures = (nf["opt"] > nf["un"]) or (nf["pes"] > nf["un"])
    blob["manufactures"] = manufactures
    print(f"\n  false certs: unmodelled {nf['un']}  soft@opt {nf['opt']}  soft@pes {nf['pes']}  "
          f"truth {nf['tr']}")
    print(f"  ==> {'*** MANUFACTURES -- the brief STOP fires; see the banked anatomy' if manufactures else 'CLEAN at this rung (soft modelling does not raise the false-cert count)'}",
          flush=True)
    blob["note"] = ("SPARK-3 SCRAMBLED-COUNTERPART arm -- the manufacturing control the brief "
                    "pre-registered. Recovery counterpart sky-scrambled by SURFACE's own recipe "
                    "(surface.py scramble_s), so EVERY certification is FALSE by construction. "
                    "MANUFACTURES iff soft-modelling the reservoir RAISES the false-cert count "
                    "above unmodelled at the same rung -> the brief's STOP. Floors re-cut at "
                    "MATCHED state (this arm's own nulls), never inherited from the signal arm. "
                    "Raw per-pulsar dlnL/qmax/margin + raw null offenders per state banked; "
                    "both Fisher bounds as columns; never verdicts. SCOPE: one rung, one "
                    "realisation, one venue -- EMBER S8.5.2 applies verbatim, a null here is "
                    "NOT evidence of safety elsewhere.")
    np.savez(path, **blob)
    print(f"[SPARK3] banked -> {path}", flush=True)
    return 0


def mode_ledger():
    """(c) THE LEDGER + the PRE-REGISTERED verdict.

    THE QUESTION, exactly: does any UNCERTIFIED pulsar's margin-to-bar cross zero between
    rung 0 and rung 8 under the SOFT-MODELLED middle? Not "does the floor fall" (it can fall
    for counting reasons) and not "does the count rise at truth" (truth is not reachable) --
    a CROSSING, by a pulsar that was uncertified at rung 0, under a reservoir model the loop
    could actually build.

    THE VERDICT, pre-registered before the grid ran:
      EDGE-POSITIVE  crossings under BOTH Fisher bounds -> the conclusion is bound-independent
                     -> GLACIER is pre-registered.
      EDGE-ZERO      no crossing under the OPTIMISTIC bound. Since sigma_cond <= sigma_marg,
                     the optimistic bound OVERSTATES the tightening: no crossing there means
                     no crossing with any honest width -> arrow 2 is CLOSED at current
                     capability, and the shortfall is quoted in nats.
      STRADDLED      the bounds disagree -> the two bounds do not decide it, and the report
                     prices a budgeted chunked JVP rather than quoting a number it cannot
                     defend.
    """
    # READ THE PER-RUNG CHECKPOINTS, not the unit npz: a unit that died mid-grid has no unit
    # npz but its finished rungs are on disk and are perfectly good data. The crossing test
    # then uses whatever (real, 0) and (real, 8) pairs actually exist.
    fs = sorted(glob.glob(f"{OUT}/s3r_*_g*_r*_k*.npz"))
    if not fs:
        print("*** STOP: no rungs banked.", flush=True)
        return 2
    R = {}
    for f in fs:
        b = os.path.basename(f)[4:-4].split("_")       # <venue>_g<geo>_r<real>_k<rung>
        key, geo, real, k = b[0], int(b[1][1:]), int(b[2][1:]), int(b[3][1:])
        R[(key, geo, real, k)] = f
    print(f"=== SPARK-3 LEDGER ({len(R)} rungs banked over "
          f"{len({(a,b,c) for a,b,c,_ in R})} units) ===", flush=True)
    rows = []
    for (key, geo, real, k), f in sorted(R.items()):
        z = np.load(f, allow_pickle=True)
        rows.append(dict(venue=key, geo=geo, real=real, rung=k, coh=int(z["coh"]),
                         **{f"{a}_{n}": (int(z[f"{a}_{n}"]) if a == "ncert"
                                         else float(z[f"{a}_{n}"]))
                            for n in ("un", "opt", "pes", "tr")
                            for a in ("ncert", "floor", "off", "zf")},
                         sig_opt=float(np.median(z["sig_opt"])),
                         sig_pes=float(np.median(z["sig_pes"])),
                         nrej_opt=int(z["nrej_opt"]), nrej_pes=int(z["nrej_pes"])))
    print(f"\n{'venue':>5} {'real':>4} {'rung':>4} {'coh':>4} | {'N_un':>4} {'N_opt':>5} "
          f"{'N_pes':>5} {'N_tr':>4} | {'fl_un':>8} {'fl_opt':>8} {'fl_tr':>8} | "
          f"{'med sig_opt':>11} {'rej':>4}")
    for r in rows:
        print(f"{r['venue']:>5} {r['real']:>4} {r['rung']:>4} {r['coh']:>4} | "
              f"{r['ncert_un']:>4} {r['ncert_opt']:>5} {r['ncert_pes']:>5} {r['ncert_tr']:>4} | "
              f"{r['floor_un']:>8.3f} {r['floor_opt']:>8.3f} {r['floor_tr']:>8.3f} | "
              f"{r['sig_opt']:>11.4g} {r['nrej_opt']:>4}")

    # ---- (c) THE CROSSING LEDGER ----
    print(f"\n=== (c) THE CROSSING LEDGER: uncertified at rung 0 -> certified at rung 8? ===")
    cross = {"opt": 0, "pes": 0}
    n_pairs = 0
    units = sorted({(a, b, c) for a, b, c, _ in R})
    for (key, geo, real) in units:
        if (key, geo, real, 0) not in R or (key, geo, real, 8) not in R:
            continue
        n_pairs += 1
        z0 = np.load(R[(key, geo, real, 0)], allow_pickle=True)
        z8 = np.load(R[(key, geo, real, 8)], allow_pickle=True)
        for bnd in ("opt", "pes"):
            m0 = np.asarray(z0[f"margin_{bnd}"], float)
            m8 = np.asarray(z8[f"margin_{bnd}"], float)
            c = int(((m0 < 0) & (m8 > 0) & np.isfinite(m0) & np.isfinite(m8)).sum())
            cross[bnd] += c
    print(f"  units with BOTH rung 0 and rung 8 banked: {n_pairs}")
    print(f"  crossings (margin < 0 at rung 0, > 0 at rung 8):")
    print(f"     under the OPTIMISTIC (conditional) bound : {cross['opt']}")
    print(f"     under the PESSIMISTIC (prior)      bound : {cross['pes']}")
    if cross["opt"] > 0 and cross["pes"] > 0:
        verdict = "EDGE-POSITIVE"
    elif cross["opt"] == 0:
        verdict = "EDGE-ZERO"
    else:
        verdict = "STRADDLED"
    print(f"\n  VERDICT: {verdict}", flush=True)

    # the shortfall, in nats: how far the nearest uncertified pulsar still is at rung 8
    short = []
    for (key, geo, real, k), f in sorted(R.items()):
        if k != 8:
            continue
        z = np.load(f, allow_pickle=True)
        m8 = np.asarray(z["margin_opt"], float)
        m8 = m8[np.isfinite(m8) & (m8 < 0)]
        if len(m8):
            short.append(float(np.sort(m8)[-1]))
    if short:
        print(f"  SHORTFALL (best uncertified margin at rung 8, OPTIMISTIC bound): "
              f"median {np.median(short):.3f} nat, best {max(short):.3f} nat "
              f"over {len(short)} units", flush=True)
    np.savez(f"{OUT}/spark3_ledger.npz", verdict=verdict,
             cross_opt=cross["opt"], cross_pes=cross["pes"], n_pairs=n_pairs,
             shortfall_opt=np.array(short),
             venue=np.array([r["venue"] for r in rows]),
             real=np.array([r["real"] for r in rows]),
             rung=np.array([r["rung"] for r in rows]),
             coh=np.array([r["coh"] for r in rows]),
             **{f"{a}_{n}": np.array([r[f"{a}_{n}"] for r in rows])
                for n in ("un", "opt", "pes", "tr") for a in ("ncert", "floor", "off", "zf")},
             note=("SPARK-3 (c). EDGE-POSITIVE iff an uncertified pulsar's margin CROSSES the "
                   "bar between rung 0 and rung 8 under BOTH Fisher bounds; EDGE-ZERO iff no "
                   "crossing under the OPTIMISTIC (conditional) bound, which overstates the "
                   "tightening and is therefore a fortiori; STRADDLED iff they disagree. "
                   "Collated from the per-unit RAW banks."))
    print(f"[SPARK3] banked -> {OUT}/spark3_ledger.npz", flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["venue", "venue_self", "replay", "gate1", "gate2",
                                     "unit", "scram", "jvp", "remarg", "ledger"])
    ap.add_argument("--venue", default="A", choices=["A", "B"])
    ap.add_argument("--real", type=int, default=0)
    ap.add_argument("--geo", type=int, default=3)
    ap.add_argument("--rung", type=int, default=None,
                    help="run ONE rung of the unit (the rung-split grid); default: all 4")
    ap.add_argument("--n-null", type=int, default=100)
    ap.add_argument("--seed-block", type=int, default=0,
                    help="independent replicate of the venue floor (fresh null seeds)")
    ap.add_argument("--self", dest="self_lane", action="store_true",
                    help="use the SELF-DERIVED venue bank (host-independent; no SURFACE "
                         "bit-identity requirement) -- see mode_venue_self")
    ap.add_argument("--reals", default="", help="jvp: comma-list of realisations for this process")
    ap.add_argument("--scr-mode", default="L", choices=["L", "A"],
                    help="L = scramble the k_loud counterpart (nullL); A = all 16 (nullA)")
    a = ap.parse_args()
    if a.mode == "venue":
        sys.exit(mode_venue())
    if a.mode == "venue_self":
        sys.exit(mode_venue_self(a.venue, a.geo, a.n_null, a.seed_block))
    if a.mode == "replay":
        sys.exit(mode_replay(a.venue, a.geo))
    if a.mode == "gate1":
        sys.exit(mode_gate1(a.venue, a.geo))
    if a.mode == "gate2":
        sys.exit(mode_gate2(a.venue, a.geo, a.real))
    if a.mode == "unit":
        sys.exit(mode_unit(a.venue, a.geo, a.real, a.n_null, a.rung, a.self_lane))
    if a.mode == "scram":
        if a.rung is None:
            print("*** scram needs --rung", flush=True); sys.exit(2)
        sys.exit(mode_scram(a.venue, a.geo, a.real, a.rung, a.n_null, a.scr_mode,
                            a.self_lane))
    if a.mode == "jvp":
        reals = [int(x) for x in a.reals.split(",") if x != ""]
        sys.exit(mode_jvp(a.venue, reals))
    if a.mode == "remarg":
        if a.rung is None:
            print("*** remarg needs --rung", flush=True); sys.exit(2)
        sys.exit(mode_remarg(a.venue, a.geo, a.real, a.rung, a.n_null))
    if a.mode == "ledger":
        sys.exit(mode_ledger())
    print(f"*** mode {a.mode} not yet built", flush=True)
    sys.exit(2)
