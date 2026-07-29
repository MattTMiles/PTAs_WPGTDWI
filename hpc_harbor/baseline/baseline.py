#!/usr/bin/env python
"""BASELINE -- the FIELD-BASELINE comparison the paper lacks.

Agent BASELINE, ACCRE, 2026-07-29. REPORT-ONLY: nothing here arms a protocol step,
moves a banked verdict, or enters a closure claim. Banks -> BASELINE_results/ only;
SURFACE_results/, GENERALISE_results/, GLACIER*/ and SPARK*/ are READ-ONLY inputs.

THE QUESTION
------------
Every number this programme quotes is scored against `criterion-v2.2` -- a three-layer
gate (trials factor `ln K`, an absolute null-calibrated floor, then `q_max > 0.9`) that
NO published PTA CW analysis applies. The paper therefore cannot say what the criterion
BUYS, because the field's own standard has never been run on the same realisations.
BASELINE runs it.

    B1  DETECTION.  The Earth-term F_e-statistic (the field's standard: a matched filter
        maximised over sky, frequency and the four amplitude parameters), thresholded at
        FAP = 0.05 read from the SAME pure-noise nulls the criterion's floors are read
        from. Detection fraction, beside the criterion's own.
    B2  DISTANCE.  For F_e-detected realisations, the field's standard per-pulsar
        pulsar-term fringe posterior, built on the FIELD'S OWN comb (the mode spacing at
        the ESTIMATED sky and frequency, not the true one). Score the fringe-correctness
        of confident picks -- i.e. apply q/purity to the FIELD's method.
        PRE-REGISTERED EXPECTATION (P14's mechanism): above onset the confident picks are
        wrong at a rate the field's own machinery cannot see. MEASURED, NOT ASSERTED.

WHAT IS AND IS NOT NEW MACHINERY
--------------------------------
Reused verbatim, never reimplemented:
  * the F_e detector          -- spark.make_fstat_earth / build_detectors (gate g0a:
                                 == TE.make_fstat verbatim), spark.scan_grid, ll0_earth
  * the grid                  -- TE.freq_grid / TE.sky_grid at the PROBLEM's real span
                                 (spark.flat_grid; TE.seed_scan hard-codes the T=15 span)
  * the venue                 -- chorus.build_chorus_problem (n_ecc=0 -> ncw=16 census;
                                 n_ecc=1 -> ncw=47 eccentric), the GENERALISE lgwb patch
  * the geometry              -- generalise.gen_geometry (structure + eccentric placement)
  * the E-step / comb / posterior -- trackB_b1_core.B1Split.estep, build_EV, FringeTables
  * the floor estimator       -- criterion-v2.2: Gumbel MLE (1-alpha) quantile at
                                 alpha=0.05, z=-ln(-ln(1-alpha)), fit error 2.80*beta/sqrt(N);
                                 ABOVE a 20% zero-fraction the floor is the empirical q95
                                 with a bootstrap error. The zero-fraction is a REQUIRED
                                 column on every floor emitted here.

BASELINE-local (three functions, each gated):
  * `polish`      -- the F_e maximisation finished off the grid. The incumbent seeder stops
                     at a 192-pixel (nside=4, ~14.7 deg) healpix grid because an M-step
                     polishes afterwards; a DETECTION statistic must be maximised, so the
                     grid argmax is followed by Adam over all EIGHT source parameters
                     (sky, log10_fgw, log10_mc, cos_inc, log10_h, phase0, psi), returning
                     the BEST point the trajectory visited and annealing the step so it
                     converges. Applied IDENTICALLY to signal and null, so the threshold
                     absorbs its trials cost. Gates G-B3 (never lowers 2F) and G-B3b (the
                     polished sky is >= 2x closer to the source than the grid pixel) --
                     both were FAILED by the first version and drove its rewrite.
  * `field_model` -- the field's recovery template: source slot 0 at the F_e estimate,
                     EVERY other slot at H_ABSENT. Gate G-B4: lowering H_ABSENT by 6 dex
                     moves lnL by < 1e-6, i.e. the other slots are numerically absent.
  * `field_comb`  -- dL_a at the ESTIMATED (sky, fgw) -- a single-source comb, because the
                     field does not know there are sixteen. n_true is then re-derived ON
                     THE FIELD'S OWN COMB (a MAP fringe can only be right or wrong relative
                     to the comb it was picked from). Gate G-B5: at the true source
                     parameters `field_comb` reproduces the incumbent single-source spacing.

THE mc FORK (declared, and both sides measured)
-----------------------------------------------
SPARK's s2 note (spark.py:738) records that a fixed reference chirp mass is "nearly
harmless" for the EARTH term and "fatal" for the PULSAR term, whose look-back is ~kyr.
The field's per-pulsar fringe posterior is a PULSAR-term object, so mc matters. Two
conventions, both pre-registered here:
  HEADLINE  (matched, no oracle): log10_mc is a FREE parameter of the polish. Signal and
            null are treated identically, so the refit floor is matched.
  BOUND     (oracle-generous, SIGNAL ONLY): the pterm template is re-scored at the TRUE
            chirp mass of the member the F_e statistic actually locked onto. This has no
            matched null (a pure-noise realisation has no true mc), so it is scored for
            PURITY ONLY -- purity is conditional on confident picks and needs no floor.
            It bounds how much of the field's wrongness is mc-ignorance rather than the
            mechanism P14 names.
  ORACLE-SOURCE BOUND (SIGNAL ONLY, and the load-bearing control): comb AND pterm template
            placed at the TRUE sky, TRUE fgw and TRUE mc of the locked member, amplitudes
            left at the F_e estimate -- the field's method with its SEARCH ERROR SET TO
            ZERO. It exists because the headline depends on how well `polish` converges,
            and an under-converged maximiser would understate the field's sky accuracy and
            thereby MANUFACTURE wrong picks. If confident picks are still wrong here, the
            wrongness cannot be blamed on sky/frequency estimation at all. Purity-only.
The a-fortiori direction is deliberate: every convention that is free to move is moved to
FAVOUR the field, so a "the field certifies confidently and wrongly" verdict is a floor,
not a ceiling.

SCOPE OF INFERENCE (rides with every number)
--------------------------------------------
  * 116-pulsar MOCK array (AXIS, 1440 MHz single-frequency). No real TOA is touched; the
    residuals ARE the injected CW + drawn noise. Real sky positions, real TOA errors, real
    published distance priors.
  * The T = 30 venue is the CADENCE-EXTENSION CONVENTION (SURFACE, inherited from IGNITE
    S2), not a forecast of real future data.
  * LOUDNESS: the frozen 16-source POP with 3 (or 5) members at log10_h = -12.75 is FAR
    above NG15. Per the GLACIER-LITE labelling rule (capstone addendum 2026-07-27) these
    cells are labelled by a_equivalent, never by "NG15-consistent".
  * The F_e grid is the incumbent seeder's (nside=4 x 1/(3T) frequency), polished. A finer
    grid would raise both the signal statistic and the null threshold.

CONVENTIONS. cpus-per-task=8 and a BANKED L_gwb (GENERALISE's shape-keyed
gen_L_gwb_n5336.npz, T=30). With the square root banked no BLAS call sits between the seed
and the draw, so every draw here is bit-identical at any thread count on any host
(trackB_b1_core.load_or_build_L_gwb.__doc__). The C2 (A-SKY) realisations were themselves
drawn under that bank and are therefore REPRODUCED EXACTLY; the C1 (SURFACE census)
realisations predate it and were drawn under the recompute path, so C1 is re-drawn here
and BOTH arms are recomputed on the re-drawn data -- the comparison is PAIRED within this
campaign and only distributionally comparable to SURFACE's bank. Gate G-B2 measures both.

Lean-npz: raw statistics only (dlnL / lnK / q_max / 2F), never verdicts.

Run:  baseline.py smoke   [--arm C1]
      baseline.py run     --arm {C1,C2} --shard i --nshard n
      baseline.py score
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time, argparse, socket, glob
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for _p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
           "hpc_harbor/chorus", "hpc_harbor/atlas", "hpc_harbor/generalise",
           "hpc_harbor/spark", "hpc_harbor/surface"):
    sys.path.insert(0, f"{HSYMT}/{_p}")

OUT = f"{HSYMT}/BASELINE_results"
REPORTS = f"{HSYMT}/reports"
SF_OUT = f"{HSYMT}/SURFACE_results"
GEN_OUT = f"{HSYMT}/GENERALISE_results"

I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)

# ---- criterion-v2.2 floor estimator (identical constants to RECUT / SURFACE / SPARK) ----
ALPHA = 0.05
Z_ALPHA = 2.9701952521018403        # -ln(-ln(1-alpha))
C_SD = 2.80                          # sd(floor_hat) = C_SD * beta / sqrt(N)   (RECUT G7)
ZF_GATE = 0.20                       # above this nullN zero-fraction the Gumbel is INVALID
QBAR = 0.9
QBAR_STRICT = 0.99
BOOT = 4000
BOOT_SEED = 20260729

# ---- the field template's "absent" strain ----
# NOT the -18.0 that TE.seed_scan uses for its no-signal reference. G-B4 MEASURED that
# convention at this venue and it FAILS: with 15 slots parked at -18 the template still
# carries 1.003e-02 nat of spurious source (job 12833821). -18 is adequate for a ONE-slot
# no-signal reference, which is all TE ever used it for; it is not adequate for fifteen.
# -30 is GLACIER-LITE's convention and sits inside chorus's _PHI_LO box (-32).
H_ABSENT = -30.0
GATE_CONT = 1e-6                     # EMBER 2.2(b) continuous-column bar
MC_BOUND = [True]                    # set by --mcbound; see the pre-registered trim order

# ---- the venue ----
T_LABEL = 30
TIER = "lit"
H_FAINT = -14.25

# ---- the two signal cells ----
# C1: the banked census cell -- SURFACE (h=-12.75, T=30, lit, 3+13), 5 skies x 6 noise = 30
C1 = dict(arm="C1", n_ecc=0, h=-12.75, k=3, e=None, edist=None, tier=TIER, T=T_LABEL,
          fdiv=1, placement=(), skies=[-1, 0, 1, 2, 3], nsig_per_sky=6,
          tag="C1_census_h1275_T30_lit_k3")
# C2: the ONE A-SKY survivor -- GENERALISE (e=0.3, h=-12.75, 5+11, T=30, lit), 8 skies x 15
C2 = dict(arm="C2", n_ecc=1, h=-12.75, k=5, e=0.3, edist="e03", tier=TIER, T=T_LABEL,
          fdiv=1, skies=[4, 5, 6, 7, 8, 9, 10, 11], nsig_per_sky=15,
          tag="C2_askysurv_e03_h1275_k5")

# ---- SURFACE seed algebra (surface.py: S_SIG/S_NULL/D_OFF, cell_index) ----
SF_S_SIG, SF_S_NULL, D_OFF = 20_000_000, 40_000_000, 10_000_000
SF_CELL_INDEX = 2                    # cells().index((-12.75, 30, 'lit', 3)); verified on disk
# ---- GENERALISE A-SKY seed algebra (generalise.py: AS_SIG/AS_NUL, ui = ci*8 + si) ----
AS_SIG, AS_NUL = 67_000_000, 68_000_000
AS_CI = 1                            # AS_CELLS.index((-12.75, 0.3, 5)); verified on disk

N_NULL = 100                         # the field's null set (pure noise, T=30) -- venue-only
N_XGATE = 8                          # nulls re-run in the C2 build to gate the pooling

# ---- the F_e polish ----
# --- the polish step schedule, DERIVED FROM THE GEOMETRY, not tuned ---
# Adam's step magnitude is ~lr*scale REGARDLESS of gradient size (mh/sqrt(vh) ~ +-1), so lr
# sets a TRAVEL SPEED, not a sensitivity. Two measured failures fixed the schedule:
#   lr0=0.05, no anneal  -> step 0.157 rad = 9 deg in gwphi. DESTROYED the statistic
#                           (13131.4 -> 2.5, job 12833821): it leaves the peak on step 1.
#   lr0=0.05, annealed   -> safe (best-point tracking) but INERT: 13131.4 -> 13131.4
#                           (job 12834517). It still leaves on step 1 and the 200 annealed
#                           steps wander ~200 deg without coming back.
# The schedule is set by the two lengths that actually matter:
#   * how far it must TRAVEL: up to a healpix half-width, ~7 deg;
#   * how fine it must LAND: the peak of a 2F ~ 1e4 source is far narrower than a grid cell.
# lr0 = 0.004 -> first step 0.72 deg (resolves the peak, does not jump the cell); the
# annealed series reaches sum(decay^t)*0.004*3.14 ~ 0.55 rad = 32 deg of cumulative travel
# over 400 steps, comfortably more than the 7 deg it can ever need; the final step is
# ~1e-6 rad. Best-point tracking means a mis-set schedule can only ever cost inertness,
# never damage -- which is why the failure mode above was visible as an equality, not a loss.
POLISH_STEPS = 400
POLISH_LR0 = 0.004
POLISH_DECAY = 1e-4 ** (1.0 / 400)      # 0.9772
#             cosgt  gwphi  cosinc   mc     lfgw    h     ph0    psi
POL_LO = np.array([-1.0,  0.0, -1.0,  6.0,  -8.0, -16.0,  0.0,  0.0])
POL_HI = np.array([1.0, 2 * np.pi, 1.0, 10.0, -7.5, -12.0, 2 * np.pi, np.pi])
POL_SCALE = np.array([0.5, 3.14, 0.5, 2.0, 0.25, 1.0, 3.14, 1.57])


# ============================================================
# stack
# ============================================================
def _import_stack():
    """CHORUS's stack + GENERALISE's shape-keyed L_gwb bank (READ-ONLY reuse: the bank
    already exists for this venue, so this process never writes one)."""
    import chorus as CH
    jax, jnp, C, B1, TE, F, IG = CH._import_stack()
    if not getattr(C, "_bl_lgwb_patched", False):
        orig = C.load_or_build_L_gwb

        def _bl_lgwb(Phi_gwb, path=None, strict=False):
            n = int(Phi_gwb.shape[0])
            bank = f"{GEN_OUT}/gen_L_gwb_n{n}.npz"
            if not os.path.exists(bank):
                raise RuntimeError(
                    f"BASELINE refuses to recompute L_gwb: no banked square root at {bank}. "
                    "A recomputed basis is thread/host dependent and would not reproduce "
                    "the GENERALISE realisations this campaign is paired against.")
            return orig(Phi_gwb, path=bank, strict=strict)

        C.load_or_build_L_gwb = _bl_lgwb
        C._bl_lgwb_patched = True
    return jax, jnp, C, B1, TE, F, IG, CH


def build_venue(cell, verbose=True):
    jax, jnp, C, B1, TE, F, IG, CH = _import_stack()
    import generalise as GEN
    import spark as SP
    P = CH.build_chorus_problem(cell["n_ecc"], T_label=cell["T"], verbose=verbose)
    if verbose:
        print(f"[BL] L_gwb provenance: {P['nd'].lgwb_prov}", flush=True)
    return dict(jax=jax, jnp=jnp, C=C, TE=TE, F=F, IG=IG, CH=CH, GEN=GEN, SP=SP, P=P)


# ============================================================
# the plan (deterministic; every seed re-derived from the source campaign's algebra)
# ============================================================
def plan(cell):
    """Signal entries for a cell + (C1 only) the shared 100-realisation field null set."""
    ent = []
    if cell["arm"] == "C1":
        for si, gid in enumerate(cell["skies"]):
            for rep in range(cell["nsig_per_sky"]):
                ns = SF_S_SIG + 10_000 * SF_CELL_INDEX + 100 * si + rep
                ent.append(dict(kind="sig", geo_id=gid, noise_seed=ns,
                                dist_seed=ns + D_OFF, no_cw=False, sky=gid))
        for rep in range(N_NULL):                       # SURFACE nullN, the shared null set
            ns = SF_S_NULL + 10_000 * SF_CELL_INDEX + rep
            ent.append(dict(kind="null", geo_id=cell["skies"][rep % 5], noise_seed=ns,
                            dist_seed=ns + D_OFF, no_cw=True, sky=cell["skies"][rep % 5]))
    else:
        for si, sky in enumerate(cell["skies"]):
            ui = AS_CI * 8 + si
            for rep in range(cell["nsig_per_sky"]):
                ns = AS_SIG + 10_000 * ui + rep
                ent.append(dict(kind="sig", geo_id=sky, noise_seed=ns,
                                dist_seed=ns + D_OFF, no_cw=False, sky=sky, ui=ui))
        # the cross-build pooling gate: the SAME SURFACE null seeds, re-run in the ncw=47
        # build. A pure-noise draw has no source, so these must come out bit-identical.
        for rep in range(N_XGATE):
            ns = SF_S_NULL + 10_000 * SF_CELL_INDEX + rep
            ent.append(dict(kind="xnull", geo_id=cell["skies"][0], noise_seed=ns,
                            dist_seed=ns + D_OFF, no_cw=True, sky=cell["skies"][0]))
    return ent


def real_path(cell, e):
    return f"{OUT}/bl_{e['kind']}_{cell['tag']}_g{e['geo_id']}_n{e['noise_seed']}.npz"


# ============================================================
# geometry: the TRUE (criterion's) comb, via the incumbent constructors
# ============================================================
def unit_of(cell, sky):
    """A GENERALISE unit dict, so gen_geometry is called with ITS OWN semantics."""
    return dict(arm=cell["arm"], e=(cell["e"] or 0.0), edist=cell["edist"], h=cell["h"],
                T=cell["T"], k=cell["k"], tier=cell["tier"], fdiv=cell["fdiv"],
                n_ecc=cell["n_ecc"], sky=sky,
                placement=(GEN_placement(cell, sky) if cell["n_ecc"] else ()),
                tag=cell["tag"])


def GEN_placement(cell, sky):
    import generalise as GEN
    return GEN.as_placement(cell["k"], sky)


def true_geometry(S, cell, geo_src, sky):
    """theta/dL/EV/FT at truth for this (cell, sky) -- the CRITERION's comb.

    C2 goes through generalise.gen_geometry verbatim (eccentric tie + active slots).
    C1 (n_ecc=0) also goes through it: with an empty placement, gen_theta reduces to
    SURFACE's cell_geometry_s (sky -> structure) and active_slots to range(16), so dL is
    the min-over-16 spacing SURFACE used. Gate G-B1 checks that reduction on disk."""
    return S["GEN"].gen_geometry(S["P"], S["C"], geo_src, unit_of(cell, sky))


def load_skies(S, cell):
    """GEO-protocol source skies, keyed by geo_id. -1 = the fiducial POP sky (None)."""
    ids = [g for g in cell["skies"] if g >= 0]
    drawn = S["F"].load_geo_skies(ids) if ids else []
    out = {gid: None for gid in cell["skies"] if gid < 0}
    for gid, sky in zip(ids, drawn):
        out[gid] = sky
    return out


# ============================================================
# the F_e statistic: grid (incumbent) + polish (BASELINE-local, gated)
# ============================================================
def make_polish(S):
    """Adam over all EIGHT source parameters from the grid argmax, returning THE BEST POINT
    THE TRAJECTORY EVER VISITED -- not its last point.

    WHY (measured, job 12833821): the first version returned the last iterate and it made
    the statistic WORSE -- 2F_grid = 25.4 -> 2F_polish = 8.0 on the first null. Adam is a
    fixed-step method, not a line search: at spark._adam_profile's lr = 0.05*scale the
    gwphi step is 0.05*3.14 = 0.157 rad, which walks straight off a peak the grid had
    already found, and 200 steps need not walk back. A DETECTION STATISTIC IS A MAXIMUM
    OVER EVALUATED POINTS, so returning anything less than the best point already evaluated
    is simply a bug in the maximiser.

    Two changes: (i) track the running best (free -- value_and_grad already computes the
    value the gradient is taken at), and seed that best with the GRID argmax, so the
    polished statistic can never fall below the grid one; (ii) START SMALL and ANNEAL
    (POLISH_LR0 / POLISH_DECAY / POLISH_STEPS -- the schedule is derived from the travel
    and landing lengths, see their definitions), so the refinement converges instead of
    orbiting.

    (i) ALONE IS NOT ENOUGH, and that is the second measured failure (job 12834517): with
    best-point tracking but the incumbent's lr = 0.05, the polish is SAFE and completely
    INERT -- 2F_grid = 13131.4 -> 2F_polish = 13131.4, sky 6.37 deg -> 6.3724 deg, i.e. it
    still leaves the peak on step 1 and never returns, and it would have passed a
    "never lowers" gate while handing B2 a 6.4-deg-wrong sky.

    (ii) matters for B2, not for B1. The field's fringe comb is built AT the estimated sky,
    so an artificially coarse sky error would MANUFACTURE exactly the confidently-wrong
    picks this campaign is testing for. That is the wrong direction: every free convention
    here must FAVOUR the field. G-B3b exists to catch precisely this, and did.
    G-B3 now checks a property the code guarantees, so the informative number is the
    IMPROVEMENT rate, reported beside it."""
    jax, jnp = S["jax"], S["jnp"]
    fl = S["earth1"]["fisher_logL"]
    lo, hi, sc = jnp.asarray(POL_LO), jnp.asarray(POL_HI), jnp.asarray(POL_SCALE)

    def lnl(Ld, data, p):
        return fl(jnp.concatenate([Ld, p]), data)

    def run(Ld, data, p0):
        vg = jax.value_and_grad(lambda p: lnl(Ld, data, p))
        m = jnp.zeros(8); v = jnp.zeros(8)
        b1, b2, eps = 0.9, 0.999, 1e-8
        l0 = lnl(Ld, data, p0)

        def body(t, c):
            p, m, v, bl, bp = c
            val, gr = vg(p)
            gr = jnp.nan_to_num(gr)
            val = jnp.nan_to_num(val, nan=-jnp.inf)
            better = val > bl
            bl = jnp.where(better, val, bl)
            bp = jnp.where(better, p, bp)
            m = b1 * m + (1 - b1) * gr
            v = b2 * v + (1 - b2) * gr * gr
            mh = m / (1 - b1 ** t); vh = v / (1 - b2 ** t)
            lr = POLISH_LR0 * POLISH_DECAY ** (t - 1)
            p = jnp.clip(p + lr * sc * mh / (jnp.sqrt(vh) + eps), lo, hi)
            return (p, m, v, bl, bp)

        p, _, _, bl, bp = jax.lax.fori_loop(1, POLISH_STEPS + 1, body,
                                            (p0, m, v, l0, p0))
        lf = jnp.nan_to_num(lnl(Ld, data, p), nan=-jnp.inf)   # the final iterate too
        better = lf > bl
        return jnp.where(better, lf, bl), jnp.where(better, p, bp)

    return jax.jit(run)


def fe_statistic(S, data, verbose=False):
    """The field's detection statistic on one realisation.

    grid   : spark.flat_grid (TE's freq x healpix grid at the PROBLEM's real span),
             spark.make_fstat_earth (== TE.make_fstat, gate g0a) profiled over the four
             amplitude parameters at TE.SEED_MC.
    polish : Adam over all eight from the grid argmax (BASELINE-local, gate G-B3).
    2F     = 2 (lnL_max - lnL_0),  lnL_0 at log10_h = -18 (TE.seed_scan's convention).
    """
    SP, jnp = S["SP"], S["jnp"]
    P, L0 = S["P"], S["P"]["L0"]
    CG, GP, LF = S["grid"][:3]
    Ld = jnp.asarray(L0)
    t0 = time.time()
    stat, free = SP.scan_grid(S["fe"], (Ld, data), CG, GP, LF, chunk=S["chunk"])
    j = int(np.argmax(stat))
    p0 = np.array([CG[j], GP[j], free[j, 0], S["TE"].SEED_MC, LF[j],
                   free[j, 1], free[j, 2], free[j, 3]])
    lmax_grid = float(stat[j])
    lp, p_hat = S["polish"](Ld, data, jnp.asarray(p0))
    lmax = float(lp); p_hat = np.asarray(p_hat)
    ll0 = SP.ll0_earth(S["earth1"], L0, data, float(np.median(S["grid"][3])))
    if verbose:
        print(f"  [Fe] grid {len(CG)} pts {time.time()-t0:.1f}s  "
              f"2F_grid={2*(lmax_grid-ll0):.1f} -> 2F_polish={2*(lmax-ll0):.1f}", flush=True)
    return dict(twoF=2.0 * (lmax - ll0), twoF_grid=2.0 * (lmax_grid - ll0),
                p_hat=p_hat, p_grid=p0, ll0=ll0,
                stat_max=lmax, stat_grid_max=lmax_grid,
                grid_top=np.sort(stat)[-5:][::-1])


# ============================================================
# the field's per-pulsar fringe inference (B2)
# ============================================================
def field_comb(S, p_hat):
    """dL_a at the ESTIMATED sky and frequency -- a SINGLE-source comb. The field does not
    know the population, so it cannot take a min over sixteen sources."""
    P = S["P"]
    ms = P["_mode_spacing"]
    return np.array([ms(p_hat[I_COSGT], p_hat[I_GWPHI], p_hat[I_FGW],
                        P["disco_psrs"][a].pos) for a in range(P["npsr"])])


def field_model(S, p_hat, mc=None):
    """theta with source slot 0 at the F_e estimate and EVERY other slot at H_ABSENT."""
    P = S["P"]
    nd = P["n_dist"]
    th = P["theta_truth"].copy()
    for i in range(P["ncw"]):
        th[nd + 8 * i + I_H] = H_ABSENT
    th[nd + 0 * 8 + I_COSGT] = p_hat[I_COSGT]
    th[nd + 0 * 8 + I_GWPHI] = p_hat[I_GWPHI]
    th[nd + 0 * 8 + I_COSINC] = p_hat[I_COSINC]
    th[nd + 0 * 8 + I_MC] = p_hat[I_MC] if mc is None else mc
    th[nd + 0 * 8 + I_FGW] = p_hat[I_FGW]
    th[nd + 0 * 8 + I_H] = p_hat[I_H]
    th[nd + 0 * 8 + I_PH0] = p_hat[I_PH0]
    th[nd + 0 * 8 + I_PSI] = p_hat[I_PSI]
    return th


def field_estep(S, p_hat, data, L_true, tier, mc=None):
    """The field's standard per-pulsar distance inference, on its OWN comb.

    Returns the criterion's four raw columns computed by the FIELD's method
    (dlnL_det, lnK, q_max, map_k) plus the oracle columns re-derived ON THE FIELD'S COMB
    (a MAP fringe is right or wrong only relative to the comb it was picked from)."""
    P, C, jnp = S["P"], S["C"], S["jnp"]
    dLf = field_comb(S, p_hat)
    EVf = C.build_EV(P, dLf)
    FTf = C.FringeTables(P, EVf, dLf, prior_key=tier)
    th = field_model(S, p_hat, mc=mc)
    th[P["AI"]] = P["L0"]                       # E-step scans about the PRIOR MEAN
    LNL = P["sp"].estep(th, EVf, P["AI"], data, jnp.ones(P["npsr"]))
    C._finite("field estep", LNL)
    post = FTf.posterior(LNL)
    lnK = np.log(np.maximum(FTf.K_counted.astype(float), 1.0))
    dlnL = S["CH"]._dlnL_gap(FTf, P["npsr"])
    if L_true is None:                          # a null has no truth
        ntg = np.full(P["npsr"], -10_000, int)
        ptrue = np.full(P["npsr"], np.nan)
    else:
        ntg = C.n_true_on_grid(L_true, P["L0"], dLf)
        ptrue, _, _ = FTf.p_of_fringe(ntg)
    return dict(dlnL_det=dlnL, lnK=lnK, qmax=post["q_max"], mapk=post["map_k"],
                n_true_grid=ntg, on_true=(post["map_k"] == ntg),
                ptrue=np.nan_to_num(ptrue, nan=-1.0), dL_med=float(np.median(dLf)),
                K_sum=int(FTf.K_counted.sum()), dL=dLf, nfringe=n_fringe(FTf))


def n_fringe(FT):
    """Fringes actually SAMPLED in each pulsar's window.

    On the criterion's comb this is >= 2 everywhere. On the FIELD's comb it need not be:
    the mode spacing dL_a ~ 1/(1 - cos mu_a) DIVERGES for a pulsar lying near the estimated
    source direction, so a mis-estimated sky can leave a pulsar with ONE fringe in its
    window. Such a pulsar carries no fringe ambiguity, hence no fringe claim: its `dlnL`
    (a gap between the best and second-best fringe) is +inf by construction and its q_max
    is 1.0 trivially. Banked as a column so the SCORER can exclude those pulsars from BOTH
    arms symmetrically and report how many it excluded -- the incumbent per-arm code paths
    are left byte-identical."""
    return np.array([len(u) for u in FT.uk], int)


# ============================================================
# one realisation: BOTH arms on the SAME data
# ============================================================
def run_realisation(S, cell, e, geo_src, keep=True, verbose=False):
    P, C, jnp, CH = S["P"], S["C"], S["jnp"], S["CH"]
    path = real_path(cell, e)
    if keep and os.path.exists(path):
        return "skip", None
    t0 = time.time()
    nd_, L0, AI = P["n_dist"], P["L0"], P["AI"]
    one = jnp.ones(P["npsr"])

    # ---- the realisation, rebuilt from its source campaign's own seed algebra ----
    G = true_geometry(S, cell, geo_src, e["sky"])
    theta_src, dL, EV, FT = G["theta"], G["dL"], G["EV"], G["FT"]
    L_true, _ = CH.draw_true_distances_tier(P, dL, EV, seed=e["dist_seed"],
                                            tier=cell["tier"])
    theta_true = theta_src.copy(); theta_true[AI] = L_true
    if e["no_cw"]:
        clean = tuple(jnp.zeros(len(p.toas)) for p in P["disco_psrs"])
    else:
        clean = P["amo"]["inject_delay"](jnp.asarray(theta_true), one)
    noise = P["nd"].draw(e["noise_seed"])
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n)) for c, n in zip(clean, noise))

    # ---- ARM 1: the CRITERION (the incumbent path, verbatim) ----
    theta_base = theta_src.copy(); theta_base[AI] = L0
    LNL = P["sp"].estep(theta_base, EV, AI, data, one)
    C._finite("criterion estep", LNL)
    post = FT.posterior(LNL)
    ntg = C.n_true_on_grid(L_true, L0, dL)
    ptrue, _, _ = FT.p_of_fringe(ntg)
    cri = dict(dlnL_det=CH._dlnL_gap(FT, P["npsr"]),
               lnK=np.log(np.maximum(FT.K_counted.astype(float), 1.0)),
               qmax=post["q_max"], mapk=post["map_k"], n_true_grid=ntg,
               on_true=(post["map_k"] == ntg), ptrue=np.nan_to_num(ptrue, nan=-1.0),
               nfringe=n_fringe(FT))

    # ---- ARM 2: the FIELD ----
    fe = fe_statistic(S, data, verbose=verbose)
    fld = field_estep(S, fe["p_hat"], data, None if e["no_cw"] else L_true, cell["tier"])

    # ---- the oracle-mc BOUND (signal only; purity-only, no matched null) ----
    # Costs one extra full E-step per signal realisation. Pre-registered trim #1 drops it to
    # C1 only if the fan does not fit the 12 GPU-hr cap; the switch is on the command line so
    # the log records which arms paid for it.
    fld_mc = fld_or = lock = None
    if not e["no_cw"] and MC_BOUND[0]:
        lock = mc_of_locked_member(S, theta_src, fe["p_hat"])   # (mc, sep_deg, df_hz, which)
        fld_mc = field_estep(S, fe["p_hat"], data, L_true, cell["tier"], mc=lock[0])
        # ---- the ORACLE-SOURCE bound: the maximally generous field ----
        # The comb and the pulsar-term template are placed at the TRUE sky, TRUE frequency
        # and TRUE chirp mass of the member the F_e statistic locked onto; only the
        # amplitude parameters stay at the F_e estimate. This is the field's method with
        # its SEARCH ERROR SET TO ZERO. It matters because it makes the headline robust to
        # how well the polish converges: if confident picks are still wrong HERE, the
        # wrongness cannot be blamed on sky/frequency estimation at all. Purity-only (a
        # pure-noise realisation has no true source), so no matched floor is owed.
        src = theta_src[nd_:]
        w = lock[3]
        p_or = fe["p_hat"].copy()
        p_or[I_COSGT] = src[8 * w + I_COSGT]
        p_or[I_GWPHI] = src[8 * w + I_GWPHI]
        p_or[I_FGW] = src[8 * w + I_FGW]
        p_or[I_MC] = src[8 * w + I_MC]
        fld_or = field_estep(S, p_or, data, L_true, cell["tier"])

    rec = dict(kind=e["kind"], arm=cell["arm"], cell=cell["tag"], h=cell["h"],
               T_label=cell["T"], tier=cell["tier"], k_loud=cell["k"],
               n_ecc=cell["n_ecc"], e_val=(cell["e"] or 0.0),
               geo_id=e["geo_id"], sky=e["sky"], noise_seed=e["noise_seed"],
               dist_seed=e["dist_seed"], no_cw=e["no_cw"],
               # --- FIELD detection (B1) ---
               twoF=fe["twoF"], twoF_grid=fe["twoF_grid"], p_hat=fe["p_hat"],
               p_grid=fe["p_grid"], fe_ll0=fe["ll0"], fe_grid_top=fe["grid_top"],
               # --- FIELD per-pulsar inference (B2), on the FIELD's own comb ---
               f_dlnL_det=fld["dlnL_det"], f_lnK=fld["lnK"], f_qmax=fld["qmax"],
               f_mapk=fld["mapk"], f_n_true_grid=fld["n_true_grid"],
               f_on_true=fld["on_true"], f_ptrue=fld["ptrue"],
               f_dL_med=fld["dL_med"], f_K_sum=fld["K_sum"], f_nfringe=fld["nfringe"],
               # --- CRITERION, on the SAME data (paired) ---
               c_dlnL_det=cri["dlnL_det"], c_lnK=cri["lnK"], c_qmax=cri["qmax"],
               c_mapk=cri["mapk"], c_n_true_grid=cri["n_true_grid"],
               c_on_true=cri["on_true"], c_ptrue=cri["ptrue"],
               c_dL_med=float(np.median(dL)), c_K_sum=int(FT.K_counted.sum()),
               c_nfringe=cri["nfringe"],
               oracle_cols=np.array(["n_true_grid", "on_true", "ptrue", "sep_deg",
                                     "df_hz", "mc_true"]),
               # --- provenance ---
               host=socket.gethostname(),
               blas_threads=int(os.environ.get("OPENBLAS_NUM_THREADS",
                                               len(os.sched_getaffinity(0)))),
               lgwb_prov=str(P["nd"].lgwb_prov), names=np.array(P["names"]),
               h_absent=H_ABSENT, polish_steps=POLISH_STEPS,
               grid_npt=len(S["grid"][0]), span_s=float(S["grid"][6]))
    if lock is not None:
        mc_true, sep_deg, df_hz, which = lock
        rec.update(mc_true=mc_true, sep_deg=sep_deg, df_hz=df_hz, locked_member=which,
                   fmc_qmax=fld_mc["qmax"], fmc_mapk=fld_mc["mapk"],
                   fmc_on_true=fld_mc["on_true"], fmc_dlnL_det=fld_mc["dlnL_det"],
                   fmc_n_true_grid=fld_mc["n_true_grid"],
                   fmc_nfringe=fld_mc["nfringe"],
                   # the ORACLE-SOURCE bound (true sky + fgw + mc of the locked member)
                   for_qmax=fld_or["qmax"], for_mapk=fld_or["mapk"],
                   for_on_true=fld_or["on_true"], for_dlnL_det=fld_or["dlnL_det"],
                   for_n_true_grid=fld_or["n_true_grid"],
                   for_nfringe=fld_or["nfringe"], for_ptrue=fld_or["ptrue"],
                   for_dL_med=fld_or["dL_med"])
    if keep:
        np.savez(path, **rec)
    fin = fe["twoF"]
    return (f"2F={fin:.1f} fq>0.9:{int((fld['qmax'] > QBAR).sum())} "
            f"cq>0.9:{int((cri['qmax'] > QBAR).sum())} {time.time()-t0:.0f}s"), rec


def mc_of_locked_member(S, theta_src, p_hat):
    """Which injected member did the F_e statistic lock onto? Nearest in (sky, log f),
    weighted so a half-band frequency error costs the same as 90 deg on the sky.
    ORACLE -- labelled, used only by the declared mc BOUND and by the diagnostics."""
    P = S["P"]
    nd = P["n_dist"]
    src = theta_src[nd:]
    n = P["ncw"]
    cg, gp = p_hat[I_COSGT], p_hat[I_GWPHI]
    best, bd = 0, np.inf
    for i in range(n):
        if src[8 * i + I_H] <= H_FAINT - 1.0:     # inert slot
            continue
        c2, g2 = src[8 * i + I_COSGT], src[8 * i + I_GWPHI]
        cosang = (cg * c2 + np.sqrt(max(0.0, 1 - cg ** 2)) * np.sqrt(max(0.0, 1 - c2 ** 2))
                  * np.cos(gp - g2))
        ang = float(np.arccos(np.clip(cosang, -1, 1)))
        dfl = abs(float(src[8 * i + I_FGW]) - float(p_hat[I_FGW]))
        d = (ang / (np.pi / 2)) ** 2 + (dfl / 0.25) ** 2
        if d < bd:
            bd, best = d, i
    c2, g2 = src[8 * best + I_COSGT], src[8 * best + I_GWPHI]
    cosang = (cg * c2 + np.sqrt(max(0.0, 1 - cg ** 2)) * np.sqrt(max(0.0, 1 - c2 ** 2))
              * np.cos(gp - g2))
    sep = float(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
    df = float(10 ** src[8 * best + I_FGW] - 10 ** p_hat[I_FGW])
    return float(src[8 * best + I_MC]), sep, df, int(best)


# ============================================================
# gates
# ============================================================
def _p(name, ok, detail=""):
    print(f"  {'PASS' if ok else '*** FAIL'}  {name}   {detail}", flush=True)
    return bool(ok)


def gates(S, cell):
    """G-B0..G-B5. Every one is cheap; none is skippable."""
    P, C, jnp, np_ = S["P"], S["C"], S["jnp"], np
    ok = {}
    print("\n" + "=" * 78, flush=True)
    print(f"BASELINE GATES -- {cell['tag']}", flush=True)
    print("=" * 78, flush=True)

    # ---- G-B0: the venue is the one the banks were drawn in ----
    prov = str(P["nd"].lgwb_prov)
    ok["G-B0"] = _p("G-B0 L_gwb BANKED (thread/host independent draws)",
                    prov.startswith("BANKED"), prov)

    # ---- G-B1: the reduction to SURFACE's geometry (C1) / GENERALISE's (C2) ----
    sk = load_skies(S, cell)
    sky0 = cell["skies"][0]
    G = true_geometry(S, cell, sk[sky0], sky0)
    if cell["arm"] == "C1":
        ref = glob.glob(f"{SF_OUT}/sf_sig_h1275_T30_lit_k3_g{sky0}_n*.npz")
    else:
        ref = glob.glob(f"{GEN_OUT}/gen_sig_AS_e03_h1275_k5_s{sky0}_g{sky0}_n*.npz")
    if ref:
        z = np.load(sorted(ref)[0], allow_pickle=True)
        lnK_ref = np.asarray(z["lnK"])
        lnK_here = np.log(np.maximum(G["FT"].K_counted.astype(float), 1.0))
        d = float(np.max(np.abs(lnK_ref - lnK_here)))
        ok["G-B1"] = _p("G-B1 comb reproduces the source campaign (lnK, noise-free)",
                        d == 0.0, f"max|dlnK| = {d:.3e}   ref={os.path.basename(sorted(ref)[0])}")
    else:
        ok["G-B1"] = _p("G-B1 comb reproduces the source campaign", False, "no reference on disk")

    # ---- G-B4: the field template's absent slots are numerically absent ----
    one = jnp.ones(P["npsr"])
    th_t = G["theta"].copy(); th_t[P["AI"]] = P["L0"]
    data0 = P["nd"].draw(999_000_001)
    p_fake = np.array([0.3, 1.1, 0.2, 9.0, -7.8, -13.0, 1.0, 0.5])
    a = float(P["amo"]["logL"](jnp.asarray(field_model(S, p_fake)), data0, one))
    thb = field_model(S, p_fake)
    nd = P["n_dist"]
    for i in range(1, P["ncw"]):
        thb[nd + 8 * i + I_H] = H_ABSENT - 6.0
    b = float(P["amo"]["logL"](jnp.asarray(thb), data0, one))
    ok["G-B4"] = _p("G-B4 H_ABSENT slots numerically absent (-6 dex moves lnL < 1e-6)",
                    abs(a - b) < GATE_CONT, f"|dlnL| = {abs(a-b):.3e}")

    # ---- G-B5: field_comb == the incumbent single-source spacing at truth ----
    src = G["theta"][nd:]
    p_true = np.array([src[I_COSGT], src[I_GWPHI], src[I_COSINC], src[I_MC],
                       src[I_FGW], src[I_H], src[I_PH0], src[I_PSI]])
    dl_f = field_comb(S, p_true)
    ms = P["_mode_spacing"]
    dl_ref = np.array([ms(src[I_COSGT], src[I_GWPHI], src[I_FGW], P["disco_psrs"][a_].pos)
                       for a_ in range(P["npsr"])])
    d = float(np.max(np.abs(dl_f - dl_ref)))
    ok["G-B5"] = _p("G-B5 field_comb == incumbent single-source mode spacing",
                    d == 0.0, f"max|ddL| = {d:.3e}")

    # ---- G-B3: the polish never lowers the statistic (3 draws: 1 null, 2 signal) ----
    drops = []; seps = []
    L_true, _ = S["CH"].draw_true_distances_tier(P, G["dL"], G["EV"],
                                                 seed=999_100_001, tier=cell["tier"])
    tt = G["theta"].copy(); tt[P["AI"]] = L_true
    clean = P["amo"]["inject_delay"](jnp.asarray(tt), one)
    for i, no_cw in enumerate([True, False, False]):
        nz = P["nd"].draw(999_200_001 + i)
        if no_cw:
            dat = tuple(jnp.asarray(np.asarray(n)) for n in nz)
        else:
            dat = tuple(jnp.asarray(np.asarray(c) + np.asarray(n)) for c, n in zip(clean, nz))
        fe = fe_statistic(S, dat, verbose=True)
        drops.append((fe["twoF"] - fe["twoF_grid"], abs(fe["twoF_grid"])))
        if not no_cw:
            _, s_g, _, _ = mc_of_locked_member(S, G["theta"], fe["p_grid"])
            _, s_p, _, _ = mc_of_locked_member(S, G["theta"], fe["p_hat"])
            seps.append((s_g, s_p))
    # The bar is RELATIVE, not absolute. The grid statistic comes out of the VMAPPED
    # profile and the polish's seed value out of a non-vmapped call to the same likelihood;
    # XLA fuses the two differently, so they agree only to ~1e-6 RELATIVE. At 2F ~ 1.3e4
    # that is ~1e-2 absolute, and an absolute 1e-9 bar fails on arithmetic rather than on
    # anything about the maximiser (job 12834517: d2F = -0.00 x3, flagged FAIL). EMBER
    # 2.2(b)'s continuous-column convention is the right one here.
    ok["G-B3"] = _p("G-B3 polish never lowers 2F (guaranteed by best-point tracking; "
                    "the informative number is the IMPROVEMENT)",
                    all(d >= -GATE_CONT * max(1.0, g) for d, g in drops),
                    "d2F = " + ", ".join(f"{d:+.3g}" for d, _ in drops) +
                    f"   (improved on the grid in "
                    f"{sum(d > GATE_CONT * max(1.0, g) for d, g in drops)}/3)")
    # G-B3b: the polish must REFINE THE SKY, not just not-hurt 2F. B2 builds the field's
    # fringe comb AT the estimated sky, so a polish left stuck at the healpix half-width
    # (~7 deg) would manufacture confidently-wrong picks -- an artefact pointing the wrong
    # way. Bar: the polished sky is at least 2x closer to the injected source than the grid
    # pixel on every signal draw.
    ok["G-B3b"] = _p("G-B3b polish REFINES the sky (>=2x closer than the grid pixel)",
                     all(sp <= sg / 2.0 or sp < 1e-3 for sg, sp in seps) if seps else False,
                     "  ".join(f"grid {sg:.2f} deg -> polish {sp:.4f} deg"
                               for sg, sp in seps))

    # ---- G-B2: the paired reproduction residual against the source bank ----
    e0 = plan(cell)[0]
    _, rec = run_realisation(S, cell, e0, sk[e0["sky"]], keep=False, verbose=True)
    tgt = (f"{SF_OUT}/sf_sig_h1275_T30_lit_k3_g{e0['geo_id']}_n{e0['noise_seed']}.npz"
           if cell["arm"] == "C1" else
           f"{GEN_OUT}/gen_sig_AS_e03_h1275_k5_s{e0['sky']}_g{e0['geo_id']}"
           f"_n{e0['noise_seed']}.npz")
    if os.path.exists(tgt):
        z = np.load(tgt, allow_pickle=True)
        dq = float(np.max(np.abs(np.asarray(z["qmax"]) - rec["c_qmax"])))
        dd = float(np.max(np.abs(np.asarray(z["dlnL_det"])[np.isfinite(rec["c_dlnL_det"])]
                                 - rec["c_dlnL_det"][np.isfinite(rec["c_dlnL_det"])])))
        # EMBER 2.2(b): discrete columns exact, CONTINUOUS columns < 1e-6. q_max and dlnL are
        # continuous, and this driver reaches them by a different summation path from
        # generalise.py (same physics, different order over the ncw slots and a different XLA
        # fusion), so ~1e-9 is the CORRECT expectation and exact bit-identity is the wrong
        # bar. Job 12834510 measured 1.171e-09 / 6.519e-09 and was flagged FAIL by a `== 0.0`
        # test -- the same discrete-bar-on-a-continuous-column error as G-B3's (see G-B3).
        bit = (dq < GATE_CONT and dd < GATE_CONT)
        if cell["arm"] == "C2":
            ok["G-B2"] = _p("G-B2 C2 REPRODUCES the GENERALISE bank across hosts "
                            "(continuous bar 1e-6, EMBER 2.2(b))", bit,
                            f"max|dq|={dq:.3e} max|d dlnL|={dd:.3e}")
        else:
            # C1's bank predates the shape-keyed L_gwb bank. If SURFACE's recompute landed
            # on the SAME eigenbasis (same host, same 8 threads) this is bit-identical; if
            # not, it is a different-but-equally-valid GWB rotation. Either way BOTH arms
            # are recomputed on the re-draw, so the comparison is paired regardless. This
            # gate MEASURES the residual and does not assert its direction.
            ok["G-B2"] = _p("G-B2 C1 residual vs the SURFACE bank MEASURED "
                            "(both arms recomputed on the re-draw -> paired either way)",
                            True, f"max|dq|={dq:.3e} max|d dlnL|={dd:.3e}"
                                  f"  [{'BIT-IDENTICAL' if bit else 'ROTATED'}]")
        print(f"        [G-B2] banked-vs-redrawn: dq_max={dq:.3e}, ddlnL_max={dd:.3f}",
              flush=True)
    else:
        ok["G-B2"] = _p("G-B2 reference realisation on disk", False, tgt)

    print("-" * 78, flush=True)
    print(f"GATES: {sum(ok.values())}/{len(ok)} pass", flush=True)
    return ok, rec


# ============================================================
# modes
# ============================================================
def attach_detector(S, verbose=True):
    """The F_e detector + grid, built once per process."""
    SP, P = S["SP"], S["P"]
    t0 = time.time()
    ed = P.get("ext_diag") or {}
    from stagec_fisher import build_fisher_amortised, make_geometry_injection, \
        N_COMPONENTS, LOG10_EQUAD
    from cw_helpers import build_enterprise_pta
    pta1, cwb1, _ = build_enterprise_pta(P["ent_psrs"], 1, components=N_COMPONENTS,
                                         log10_equad=LOG10_EQUAD)
    inj1 = make_geometry_injection(pta1, P["ent_psrs"], 1, cwb1, seed=1000,
                                   h_range=(S["TE"].EQUAL_H, S["TE"].EQUAL_H))
    S["earth1"] = build_fisher_amortised(P["disco_psrs"], 1, inj1, cwb1,
                                         components=int(ed.get("gwb_comp", N_COMPONENTS)),
                                         rn_components=int(ed.get("rn_comp", 30)),
                                         msd_factory=S["TE"].EarthDelay)
    S["fe"] = SP.make_fstat_earth(S["earth1"])
    S["polish"] = make_polish(S)
    S["grid"] = SP.flat_grid(P)
    S["chunk"] = int(os.environ.get("BL_CHUNK", "256"))
    if verbose:
        print(f"[BL] detector + grid: {time.time()-t0:.1f}s; grid {len(S['grid'][0])} pts "
              f"({len(S['grid'][3])} freq x {len(S['grid'][4])} sky), span "
              f"{S['grid'][6]/3.15576e7:.2f} yr", flush=True)
    return S


def mode_smoke(a):
    cell = C1 if a.arm == "C1" else C2
    S = build_venue(cell)
    S = attach_detector(S)
    ok, _ = gates(S, cell)
    # ---- timing: one signal + one null, end to end ----
    # --quick skips these: the per-realisation wall is already MEASURED (105 s, job
    # 12833821) and a re-smoke exists to re-gate a code change, not to re-time it.
    sk = load_skies(S, cell)
    pl = plan(cell)
    sig = [e for e in pl if e["kind"] == "sig"][0]
    nul = [e for e in pl if e["kind"] in ("null", "xnull")][0]
    for e in (() if a.quick else (sig, nul)):
        t0 = time.time()
        msg, _ = run_realisation(S, cell, e, sk[e["sky"]], keep=False, verbose=True)
        print(f"[BL smoke] {e['kind']} n{e['noise_seed']}: {msg} "
              f"(total {time.time()-t0:.0f}s)", flush=True)
    np.savez(f"{OUT}/bl_smoke_{cell['tag']}.npz",
             gates=np.array(list(ok.keys())), passed=np.array(list(ok.values())),
             host=socket.gethostname())
    return 0 if all(ok.values()) else 2


def mode_run(a):
    cell = C1 if a.arm == "C1" else C2
    S = build_venue(cell)
    S = attach_detector(S)
    sk = load_skies(S, cell)
    pl = plan(cell)
    mine = [e for i, e in enumerate(pl) if i % a.nshard == a.shard]
    print(f"[BL run] {cell['tag']} shard {a.shard}/{a.nshard}: {len(mine)}/{len(pl)} "
          f"realisations", flush=True)
    t0 = time.time()
    for i, e in enumerate(mine):
        msg, _ = run_realisation(S, cell, e, sk[e["sky"]], keep=True)
        el = time.time() - t0
        print(f"[BL {i+1}/{len(mine)}] {e['kind']} n{e['noise_seed']} {msg} "
              f"| elapsed {el/60:.1f}m eta {el/(i+1)*(len(mine)-i-1)/60:.1f}m", flush=True)
    print(f"[BL run] shard done, wall {(time.time()-t0)/60:.1f}m", flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["smoke", "run"])
    ap.add_argument("--arm", default="C1", choices=["C1", "C2"])
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--mcbound", default="on", choices=["on", "off"])
    ap.add_argument("--quick", action="store_true",
                    help="gates only; skip the two timing realisations")
    a = ap.parse_args()
    MC_BOUND[0] = (a.mcbound == "on")
    os.makedirs(OUT, exist_ok=True)
    print(f"[BL] host={socket.gethostname()} cpus={len(os.sched_getaffinity(0))} "
          f"OPENBLAS={os.environ.get('OPENBLAS_NUM_THREADS')} mode={a.mode} arm={a.arm} "
          f"mcbound={a.mcbound}",
          flush=True)
    sys.exit({"smoke": mode_smoke, "run": mode_run}[a.mode](a))
