#!/usr/bin/env python3
"""GLACIER-LITE -- the certification loop in the idealized corner (REPORT-ONLY).

THE QUESTION. GLACIER measured ignition FAILURE at realistic noise; ORACLE-IGNITION
measured that CONDITIONING alone cannot buy it either (276 cells, 0 certs, at realistic
noise). Both held the noise fixed. This campaign moves that one lever: it runs the FULL
CLOSED LOOP in the idealized regime the published single-shot literature uses -- white
noise only at 20-50 ns, a few sources conditioned on truth -- where ignition is expected
to be easy. Several groups have published SINGLE-SHOT idealized results
(arXiv:2503.23017, arXiv:2512.10729, arXiv:2512.10795); nobody has run the CLOSED LOOP.
The lap-by-lap compounding curve is the deliverable.

Nothing here arms any protocol step, feeds the loop protocol, or enters closure claims.

WHAT IS KEPT REAL (explicitly NOT idealized)
--------------------------------------------
  * the 116 real pulsar sky positions (the feather set's own `pos`)
  * their real published distance priors (tier `lit`, A2.load_real_priors)
  * the fringe ambiguity: dL[a] = min mode spacing over the population at the real sky
    geometry; K_counted unchanged. The comb is the real one.

THE ONLY TWO IDEALIZATIONS
--------------------------
  1. NOISE. White only, UNIFORM across the array, at {50, 20} ns. No red noise in the
     data; RN driven to numerically negligible amplitude in the model (basis shapes kept,
     so nothing downstream changes shape). No stochastic GWB in the data -- which is
     already true of this machinery by construction: the background IS the unresolved
     source sum.
  2. SOURCE KNOWLEDGE. The N_SEED=3 loudest members conditioned at tier T2
     (oracle_ignition.TIER_FREE verbatim): sky + orbital period + chirp mass pinned at
     drawn truth; cos_inc, log10_h, phase0, psi FREE (initialised from the declared
     analysis prior and fitted by the loop's own M-step).

DECLARED DEVIATIONS AND DISK-VS-BRIEF NOTES (also in GLACIER_LITE_prereg.md S7)
-------------------------------------------------------------------------------
 (a) POPULATION. The brief's "standard frozen 16-member population" is the house's
     trackB/B1 one -- trackB_b1_core.POP = dict(ncw=16, seed=3000,
     population=(3,-13.25,-14.25)) with N_LOUD=3 -- NOT GLACIER's own 256-member census
     (glacier_pop.N_POP_DEFAULT=256). It is re-used from its frozen seed, never redrawn.
     Its native venue is T=15 (span 22.15 yr), which is also the venue at which
     b1_L_gwb.npz is banked.
 (b) WHITE NOISE IS PUT IN BOTH THE DATA AND THE LIKELIHOOD. NoiseDrawer's `white_scale`
     scales the DATA only and leaves the likelihood covariance at its built value; using
     it alone would mis-specify the noise by the same factor and forfeit the entire
     sensitivity gain this campaign exists to test. The feather set already carries
     uniform toaerrs = 1e-6 s, and LOG10_EQUAD = -6.0, i.e. a uniform white rms of
     sqrt(2)*1e-6 s = 1414 ns. The lite venue sets toaerrs = equad = target/sqrt(2) on
     the DISCOVERY pulsars (the likelihood's own noise objects), so
     N_p = efac^2 (toaerr^2 + equad^2) = target^2 EXACTLY, and the drawer -- which reads
     sp.N_diag -- draws the same covariance the likelihood assumes. The ENTERPRISE
     pulsars and the enterprise-side equad are left at house values, so the injection
     (hence the population) is bit-identical to the frozen POP.
 (c) PROMOTION IS AT DRAWN TRUTH. PromoteLedger's frozen-census rule (a fitted-position
     promote raises CampaignStop = B0.2) means a member the frontier promotes enters the
     template at its TRUE parameters. The loop demo is therefore GENEROUS BY
     CONSTRUCTION on newly-recruited members; only the 3 conditioned seeds carry free
     axes. This is the machinery's own safety rule and is not modified here. It travels
     with every number as a scope line.
 (d) "ALL OTHER MEMBERS FULLY FREE" has no state in this machinery: a census member is
     either FED (present in the template) or CARRIED (H_ABSENT=-30, living in the fitted
     GP background). The 13 unconditioned members are CARRIED at lap 0 and become
     available only by promotion through frontier-v2 -- which is precisely the mechanism
     under test.
 (e) DRAIN BAND. BackgroundFit is band-matched to THIS population's generative band
     (-8.0,-7.5) = stagec_fisher.LOG10_FGW_RANGE, not GLACIER's remedy-A campaign band
     (-8.7,-7.5), which was cut for the 256-census. At T=15 this leaves 8 in-band GP
     frequencies, so the drain has a real basis.
 (f) M-STEP AXES. Conditioned seeds move on their T2-free axes (cos_inc, log10_h,
     phase0, psi); promoted members move on glacier_loop's own MSTEP_AXES
     (log10_fgw, log10_mc). Sweeps/step follow ORACLE-IGNITION's declared deviation
     (5 sweeps at step0=1.0) for the seeds, because their entry point is a prior draw up
     to a full box away, and the loop's 2/0.3 setting cannot reach it.

WHAT IS FROZEN AND UNTOUCHED
----------------------------
  * criterion v2.2, verbatim:
        cert = (dlnL_det > max(lnK + TRIALS_NAT, floor)) & (q_max > QBAR)
    TRIALS_NAT = 0.578, QBAR = 0.9.
  * floor adoption: recut_surface.adopt -- Gumbel only at zero-fraction <= 0.20, else
    emp-q95 + bootstrap sd. zero_fraction is a REQUIRED banked column.
  * frontier-v2 (ratio < 0.5 AND dlnL_feed > 0). v1 is not reachable.
  * THE KILL STEP RUNS SCORE-AND-LOG ONLY. Its false-negative gate has not passed
    (capstone S4.23.1). It is banked as columns (kill_*) and is NEVER read by the
    certification path -- structurally, `certify()` does not take it as an argument.

Env: cpus-per-task=8 enforced (check_affinity); x64; autotune off.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, glob, sys, time
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
for p in ("CW_lnL_check", "CW_transition", "hpc_harbor/forge", "hpc_harbor/ignite",
          "hpc_harbor/atlas", "hpc_harbor/spark", "hpc_harbor/glacier"):
    sys.path.insert(0, f"{HSYMT}/{p}")

from glacier_pop import (PromoteLedger, BackgroundFit, CampaignStop, bank_npz, lane_tag,
                         check_affinity, OUT, NP_SRC,
                         I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI)
from glacier_loop import (_stack, theta_with_absent, fisher_conditional, CertScoreboard,
                          QBAR, TRIALS_NAT, GATE_RATIO, H_ABSENT_GL, MSTEP_AXES,
                          localisation_area_deg2)
from oracle_ignition import (mstep_free, draw_free, TIER_FREE, TIER_DESC, AXIS_NAME,
                             FREE_PRIOR, adopt_floor)

# ============================================================
# PRE-LOGGED PREDICTION -- banked VERBATIM into every cell, before any scoring
# ============================================================
PREDICTION = (
    "Predicted: at 20-50 ns with 3 T2-conditioned seeds, lap 1 certifies "
    ">=1 true pulsar per realisation; certifications compound over laps "
    "2-4 and saturate when the readable sub-array is exhausted; residual "
    "power decreases monotonically with lap index; the anchor cell "
    "reproduces sub-pc single-shot distance widths for pulsars within "
    "~1 kpc. Failure of compounding (flat certifications after lap 1) "
    "falsifies the loop claim in the easiest available regime and is a "
    "headline result to be reported with equal prominence."
)

# SECONDARY ARM (Matt, 2026-07-27): the NG15-renormalised contrast. Banked separately,
# never merged with the primary arm's cells.
PREDICTION_RENORM = (
    "ignition here would bound the harvest's structure-limited conclusion to sub-SKA "
    "noise; non-ignition strengthens it."
)

# PART 3 (amendment, 2026-07-27): the completeness axis.
PREDICTION_COMPLETE = (
    "Predicted: C at saturation rises monotonically with prior tightening; "
    "at 1/30 widths with 20 ns noise, C > 0.8, with the uncertified "
    "remainder audibility-limited (low pulsar-term SNR from geometry), "
    "not ambiguity-limited. If C plateaus well below 1 even at 1/30 and "
    "the autopsy shows ambiguity-limited pulsars remaining, the "
    "completeness limit is NOT prior-purchasable and the closure claim "
    "must be scoped to the readable sub-array."
)

# ---- the venue (pre-registered) ----
T_LITE = 15                              # the 16-source POP's native venue
SPAN_LITE_YR = 22.15                     # asserted against the build, never assumed
YR_S = 365.25 * 86400.0
NCW_LITE = 16                            # trackB_b1_core.POP["ncw"]
SEED_POP_LITE = 3000                     # trackB_b1_core.POP["seed"]
POPULATION_LITE = (3, -13.25, -14.25)    # trackB_b1_core.POP["population"]
BAND_LITE = (-8.0, -7.5)                 # stagec_fisher.LOG10_FGW_RANGE (deviation (e))
RN_LOG10A_OFF = -20.0                    # RN driven negligible in the model
WHITE_LEVELS_NS = {"w50": 50.0, "w20": 20.0}
TIER_DIST = "lit"                        # REAL published distance priors -- kept real

# ---- the loop (pre-registered) ----
N_SEED = 3                               # trackB_b1_core.N_LOUD -- the 3 loudest
TIER_COND = "T2"
N_LAP = 8
REALS = (0, 1, 2, 3)
N_NULL = 100                             # >= 100 per the brief
SEED_NOISE_LITE = 6_100_000              # fresh disjoint block (house convention shape)
SEED_FREE_LITE = 6_200_000
SEED_NULL_LITE = 6_300_000
MSTEP_SWEEPS_SEED, MSTEP_STEP0_SEED = 5, 1.0     # ORACLE's declared deviation
MSTEP_SWEEPS_PROM, MSTEP_STEP0_PROM = 2, 0.3     # glacier_loop's own setting
BAR_R1 = 15.132                          # D2 kill step, SCORE-AND-LOG ONLY (S4.21)

# ---- the anchor cell (Part 2) ----
ANCHOR_N_SRC = 5
ANCHOR_WHITE = "w20"

# ---- SECONDARY ARM: the NG15 renormalisation (Matt, 2026-07-27) ----
# MEASURED: the frozen POP carries 572.8x the NG15 band power -- A_equivalent = -13.22,
# NOT "NG15-consistent". The phrase is struck from this campaign's outputs; every cell is
# labelled by its measured A_equivalent instead. The secondary arm rescales every member's
# strain by ONE common factor so the summed band power equals the NG15 target exactly.
# Geometry, frequencies, phases, sky and census rank order are untouched.
A_NG15 = -14.6
F_YR_HZ = 1.0 / (365.25 * 86400.0)

# ---- PART 3: the completeness axis (amendment, 2026-07-27) ----
PW_FACTORS = (1.0, 1.0 / 3.0, 1.0 / 10.0, 1.0 / 30.0)
PW_FACTORS_CUT = (1.0, 1.0 / 10.0, 1.0 / 30.0)   # pre-registered budget cut, 4 -> 3
PW_FRINGE_FLOOR = 3.0            # HARD floor: sigma_a >= 3 * dL[a]. Never sub-fringe.
COMPLETE_LAPS = 12
COMPLETE_REALS = (0, 1)
COMPLETE_WHITE = "w20"
DRY_LAPS_STOP = 2                # stop early after 2 consecutive laps adding 0 certs


def band_power_target(a_log10=A_NG15, band=BAND_LITE):
    """INT_band h_c^2 df/f for h_c = A (f/f_yr)^(-2/3) -- glacier_pop._band_power_target
    verbatim, reproduced here so the secondary arm's normalisation convention is
    auditable in one place and is banked with every cell that uses it."""
    A2 = 10.0 ** (2 * a_log10)
    flo, fhi = 10.0 ** np.array(band)
    x = lambda f: (f / F_YR_HZ) ** (-4.0 / 3.0)
    return A2 * 0.75 * (x(flo) - x(fhi))


TAG = ""          # '_smoke' on the smoke run: its banks must NEVER satisfy the
                    # production resume check (a 32-null floor is not a 100-null floor)


def _checksums():
    """Bank metadata: content hashes of the driver and of the pre-logged prediction, so a
    bank cannot be read against a driver it was not produced by, and the prediction
    cannot be edited after the fact without every bank disagreeing."""
    import hashlib
    h = lambda b: hashlib.sha256(b).hexdigest()[:16]
    with open(os.path.abspath(__file__), "rb") as f:
        drv = h(f.read())
    return dict(cksum_driver=drv, cksum_prediction=h(PREDICTION.encode()),
                driver_path=os.path.abspath(__file__))


def stem_loop(wkey, real, renorm=False):
    r = "_ng15" if renorm else ""
    return f"gl_lite_loop_{wkey}_r{real}_T{T_LITE}_{TIER_DIST}{r}{TAG}"


def stem_complete(pw_tag, real):
    return f"gl_lite_cmpl_{COMPLETE_WHITE}_pw{pw_tag}_r{real}_T{T_LITE}_{TIER_DIST}{TAG}"


def pw_tag_of(f):
    """Filename-safe tag for a prior-width factor: 1 -> '1', 1/3 -> '3', 1/30 -> '30'."""
    return "1" if abs(f - 1.0) < 1e-12 else f"{int(round(1.0 / f))}"


# ============================================================
# THE LITE VENUE
# ============================================================
def build_lite_venue(white_ns, renorm=False, verbose=True):
    """trackB_b1_core.build_b1_problem's construction at the FROZEN 16-source POP, with
    the white level put into BOTH the data and the likelihood (deviation (b)) and the
    model RN driven negligible (deviation (e) of the prereg).

    Returns a G dict shaped like build_glacier_problem's, plus `src` (the census in the
    I_* column convention) for CertScoreboard's fringe grid."""
    import jax
    jax.config.update("jax_enable_x64", True)
    jax.config.update("jax_compilation_cache_dir",
                      os.environ.get("HARBOR_JAXCACHE_IN",
                                     "/home/mattm/.cache/jax_stagec_cache"))
    jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
    import trackB_b1_core as C
    from cw_helpers import load_pulsars, build_enterprise_pta
    from stagec_fisher import make_geometry_injection, N_COMPONENTS, LOG10_EQUAD
    import discovery as ds

    t0 = time.time()
    ent_psrs, disco_psrs = load_pulsars(116)

    # ---- the frozen population: enterprise side untouched, so the draw is bit-identical
    pta, cwb, _ = build_enterprise_pta(ent_psrs, NCW_LITE, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)
    inj = make_geometry_injection(pta, ent_psrs, NCW_LITE, cwb, seed=SEED_POP_LITE,
                                  population=POPULATION_LITE)
    # frozen-population gate: the declared strain ladder must be exactly what POP says
    k_loud, h_loud, h_faint = POPULATION_LITE
    got = np.array([float(inj[f"{n}_log10_h"]) for n in cwb])
    want = np.array([h_loud] * k_loud + [h_faint] * (NCW_LITE - k_loud))
    if not np.allclose(got, want):
        raise CampaignStop(f"frozen-population gate FAILED: log10_h ladder {got} != "
                           f"declared POP {want}. The population is not the banked one.")

    # ---- SECONDARY ARM: renormalise the summed band power to NG15 (A = -14.6).
    # ONE common factor on every member's strain; geometry/frequency/phase untouched.
    h_all = np.array([float(inj[f"{n}_log10_h"]) for n in cwb])
    pop_power = float(np.sum(10.0 ** (2 * h_all)))
    tgt_power = float(band_power_target(A_NG15, BAND_LITE))
    a_equiv = A_NG15 + np.log10(np.sqrt(pop_power / tgt_power))
    renorm_dex = 0.0
    if renorm:
        renorm_dex = -np.log10(np.sqrt(pop_power / tgt_power))
        for n in cwb:
            inj[f"{n}_log10_h"] = float(inj[f"{n}_log10_h"]) + renorm_dex
        h_chk = np.array([float(inj[f"{n}_log10_h"]) for n in cwb])
        got = float(np.sum(10.0 ** (2 * h_chk)))
        if abs(got / tgt_power - 1.0) > 1e-9:
            raise CampaignStop(f"NG15 renormalisation gate FAILED: summed band power "
                               f"{got:.6e} != target {tgt_power:.6e}. STOP.")
        a_equiv = A_NG15
    if verbose:
        print(f"  loudness: sum h^2 = {pop_power:.4e}, NG15 target "
              f"({A_NG15}, band {BAND_LITE}) = {tgt_power:.4e} -> ratio "
              f"{pop_power/tgt_power:.1f}x power, "
              f"{np.sqrt(pop_power/tgt_power):.1f}x amplitude; "
              f"A_equivalent = {a_equiv:.3f}"
              + (f"  [RENORMALISED by {renorm_dex:+.3f} dex]" if renorm else
                 "  [FROZEN POP as banked -- NOT 'NG15-consistent']"), flush=True)

    # ---- idealization 1: no red noise in the model (basis shapes preserved)
    n_rn = 0
    for k in list(inj.keys()):
        if "rednoise_log10_A" in k:
            inj[k] = RN_LOG10A_OFF
            n_rn += 1

    # ---- idealization 1: uniform white at the target, in the LIKELIHOOD's own objects
    sig = white_ns * 1e-9
    per = sig / np.sqrt(2.0)                       # toaerr == equad == target/sqrt(2)
    for psr in disco_psrs:
        psr.toaerrs = np.full_like(psr.toas, per, dtype=np.float64)
    log10_equad_lite = float(np.log10(per))

    span0 = ds.getspan(disco_psrs)
    amo = C.build_b1_amortised(disco_psrs, NCW_LITE, inj, cwb,
                               components=N_COMPONENTS,
                               log10_equad=log10_equad_lite, rn_components=30)
    if verbose:
        print(f"  build_b1_amortised: {time.time()-t0:.1f}s (ncw={NCW_LITE}, "
              f"T={T_LITE}, white={white_ns:.0f} ns, RN off on {n_rn} keys)", flush=True)

    # the census in the I_* convention (CertScoreboard's fringe grid reads this)
    src = np.zeros((NCW_LITE, NP_SRC))
    for i, n in enumerate(cwb):
        src[i] = [inj[f"{n}_cos_gwtheta"], inj[f"{n}_gwphi"], inj[f"{n}_cos_inc"],
                  inj[f"{n}_log10_mc"], inj[f"{n}_log10_fgw"], inj[f"{n}_log10_h"],
                  inj[f"{n}_phase0"], inj[f"{n}_psi"]]
    order = np.argsort(src[:, I_H])[::-1]          # census rank: brightest first
    if not np.array_equal(order, np.arange(NCW_LITE)):
        # POP already puts the k_loud first; assert rather than silently reorder, because
        # theta/slot indices are shared with amo and must not be permuted.
        if not np.allclose(src[:, I_H][order], src[:, I_H]):
            raise CampaignStop("census is not brightest-first; slot indices are shared "
                               "with amo and must not be permuted. STOP.")

    span_s = float(ds.getspan(disco_psrs))
    if abs(span_s / (SPAN_LITE_YR * YR_S) - 1.0) > 0.02:
        raise CampaignStop(f"built span {span_s/YR_S:.3f} yr differs >2% from the "
                           f"declared {SPAN_LITE_YR} yr. STOP.")
    return dict(amo=amo, disco_psrs=disco_psrs, ent_psrs=ent_psrs, cwb=cwb, inj=inj,
                slots=src, src=src, T_label=T_LITE, span_s=span_s,
                white_ns=white_ns, white_sigma_s=sig, log10_equad=log10_equad_lite,
                n_rn_off=n_rn, renorm=renorm, renorm_dex=renorm_dex,
                a_equivalent=a_equiv, pop_band_power=pop_power,
                ng15_band_power=tgt_power, band_power_ratio=pop_power / tgt_power)


def apply_prior_width(sb, C, G, factor, verbose=True):
    """PART 3: rescale every pulsar's REAL published distance prior by `factor`, with the
    pre-registered HARD FLOOR sigma_a >= PW_FRINGE_FLOOR * dL[a] -- no prior may become
    sub-fringe, because a sub-fringe prior trivialises certification and would contaminate
    the completeness count.

    All three prior columns are scaled together, because build_EV sizes the +-3 sigma
    fringe-sampling grid from max(feather, script, lit): scaling only 'lit' would leave
    the sampled comb at its old width and the tightening would be cosmetic. EV and the
    FringeTables are then REBUILT, and the caller must re-draw truth through the same
    path -- the house's truth draw is prior-weighted over the sampled fringes, so a
    tighter prior necessarily concentrates the truth toward L0. That is the house
    convention applied consistently, and it is DECLARED: across prior-width factors the
    knowledge changes AND the drawn truth concentrates with it; the factors are not the
    same universe re-analysed.

    Returns a record of what the floor did, for the bank."""
    dL = sb.dL
    floor = PW_FRINGE_FLOOR * dL
    sig0 = np.asarray(sb.P["priors"][TIER_DIST], float).copy()
    n_floored = 0
    for col in ("feather", "script", "lit"):
        s = np.asarray(sb.P["priors"][col], float)
        scaled = np.maximum(s * factor, floor)
        if col == TIER_DIST:
            n_floored = int(np.sum(s * factor < floor))
        sb.P["priors"][col] = scaled
    sb.EV = C.build_EV(sb.P, dL)
    sb.FT = C.FringeTables(sb.P, sb.EV, dL, prior_key=TIER_DIST)
    sig1 = np.asarray(sb.P["priors"][TIER_DIST], float)
    K = np.asarray(sb.FT.K_counted, float)
    if verbose:
        print(f"  prior width x{factor:.4g}: sigma_lit median "
              f"{np.median(sig0)*1000:.1f} -> {np.median(sig1)*1000:.1f} pc; "
              f"floored (sigma < {PW_FRINGE_FLOOR:.0f}*dL) on {n_floored}/{sb.npsr} "
              f"pulsars", flush=True)
        print(f"  K_counted across {sb.npsr} pulsars: min {K.min():.0f} "
              f"median {np.median(K):.0f} max {K.max():.0f}", flush=True)
    return dict(pw_factor=factor, pw_floor_frac=PW_FRINGE_FLOOR,
                pw_n_floored=n_floored, sigma_lit_kpc=sig1, sigma_lit_kpc_pre=sig0,
                K_counted=K, K_min=float(K.min()), K_med=float(np.median(K)),
                K_max=float(K.max()))


def build_state(wkey, real=0, renorm=False, pw_factor=None, verbose=True):
    """The venue + scoreboard + data for one (white level, realisation)."""
    jax, jnp, C, B1, TE, IG, F, FL = _stack()
    white_ns = WHITE_LEVELS_NS[wkey]
    G = build_lite_venue(white_ns, renorm=renorm, verbose=verbose)
    amo = G["amo"]
    nd = amo["n_dist"]
    ones = jnp.ones(amo["npsr"])
    sb = CertScoreboard(G, amo, jnp, C, prior_key=TIER_DIST)
    pw_rec = None
    if pw_factor is not None:
        pw_rec = apply_prior_width(sb, C, G, pw_factor, verbose=verbose)

    # the white level as the LIKELIHOOD sees it -- gated, never assumed
    w_rms = float(np.median([np.sqrt(np.mean(v)) for v in sb.sp.N_diag]))
    if abs(w_rms / G["white_sigma_s"] - 1.0) > 1e-6:
        raise CampaignStop(f"white-noise gate FAILED: likelihood N_diag rms "
                           f"{w_rms*1e9:.3f} ns != target {white_ns:.3f} ns. The "
                           f"idealization did not reach the covariance. STOP.")
    if verbose:
        print(f"  white gate: likelihood N_diag rms = {w_rms*1e9:.3f} ns "
              f"(target {white_ns:.0f}) PASS", flush=True)

    ndraw = _lite_drawer(C, sb, amo, verbose)
    st = dict(jax=jax, jnp=jnp, C=C, TE=TE, IG=IG, FL=FL, G=G, amo=amo, nd=nd, ones=ones,
              sb=sb, ndraw=ndraw, wkey=wkey, white_ns=white_ns, renorm=renorm,
              pw_factor=pw_factor, pw_rec=pw_rec,
              scale=TE.phi_scale({"ncw": 1}), n_src=NCW_LITE,
              bf=BackgroundFit(amo, band_log10f=BAND_LITE))
    draw_realisation(st, real, verbose=verbose)
    return st


def _lite_drawer(C, sb, amo, verbose=True):
    """NoiseDrawer at the banked T=15 GWB square root. Phi_gwb depends only on the GWB
    prior and the Fourier basis (Tspan), NOT on the white level, so the banked artifact
    is still the right one -- loaded STRICT so a mismatch crashes rather than silently
    recomputing (the thread-count hazard). We never draw the gwb component anyway; the
    load is kept because it is the artifact's own provenance check."""
    import trackB_b1_core as Cm
    nd = C.NoiseDrawer(sb.sp, amo, lgwb_path=Cm.L_GWB_BANK, strict=True)
    if verbose:
        print(f"  L_gwb: {nd.lgwb_prov}", flush=True)
    return nd


def draw_realisation(st, real, verbose=True):
    """Noise + distance truth for realisation `real` on an already-built venue.
    WHITE ONLY -- components=("white",). No rn, no gwb."""
    jnp, sb, amo, ones = st["jnp"], st["sb"], st["amo"], st["ones"]
    noise_seed = SEED_NOISE_LITE + 1000 * real
    dist_seed = noise_seed + 10_000_000
    L_true, n_true = sb.draw_truth(st["IG"], dist_seed, tier=TIER_DIST)
    theta_true = np.asarray(amo["theta_truth"], float).copy()
    theta_true[sb.AI] = L_true
    clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
    noise = st["ndraw"].draw(noise_seed, components=("white",))
    st["theta_true"] = theta_true
    st["n_true"] = n_true
    st["data"] = tuple(jnp.asarray(np.asarray(c) + np.asarray(n))
                       for c, n in zip(clean, noise))
    st["noise_seed"], st["dist_seed"], st["real"] = noise_seed, dist_seed, real
    if verbose:
        nz = np.median([float(np.std(np.asarray(n))) for n in noise])
        sg = np.median([float(np.std(np.asarray(c))) for c in clean])
        print(f"  realisation r{real}: noise_seed {noise_seed} dist_seed {dist_seed}; "
              f"median noise rms {nz*1e9:.1f} ns, CW rms {sg*1e9:.1f} ns", flush=True)
    return st


# ============================================================
# CONDITIONING (idealization 2)
# ============================================================
def condition_seeds(st, n_seed=N_SEED, tier=TIER_COND, verbose=True):
    """theta_rec at lap 0: every member at drawn truth EXCEPT the n_seed loudest, whose
    tier-free axes are re-drawn from the DECLARED analysis prior (oracle_ignition's
    FREE_PRIOR verbatim). Returns (theta_rec, led, seed_idx, free_axes)."""
    theta_rec = st["theta_true"].copy()
    nd = st["nd"]
    free = TIER_FREE[tier]
    seed_idx = np.arange(min(n_seed, st["n_src"]))
    rng = np.random.default_rng(SEED_FREE_LITE + 1000 * st["real"])
    for k in seed_idx:
        for j in free:
            col = nd + NP_SRC * int(k) + j
            theta_rec[col] = draw_free(rng, j, float(st["theta_true"][col]))
    led = PromoteLedger(st["G"]["slots"])
    for k in seed_idx:
        led.promote(int(k), st["G"]["slots"][int(k)], iteration=0)
    if verbose:
        print(f"  conditioned {len(seed_idx)} seeds at {tier} ({TIER_DESC[tier]}); "
              f"free axes = {[AXIS_NAME[j] for j in free] or 'none'}", flush=True)
    return theta_rec, led, seed_idx, free


# ============================================================
# FLOORS (idealization-invalidated -> re-cut, mandatory, BEFORE any scoring)
# ============================================================
def null_offenders(st, theta_rec, led, n_null=N_NULL, p=print):
    """glacier_loop._null_offenders' statistic verbatim, on WHITE-ONLY null draws at THIS
    lap's joint state. The banked floors were cut at 1414 ns white + real RN and do not
    apply here -- this is the re-cut the brief requires."""
    sb, jnp = st["sb"], st["jnp"]
    ones = jnp.ones(sb.npsr)
    off = np.empty(n_null)
    t0 = time.time()
    for i in range(n_null):
        nz = st["ndraw"].draw(SEED_NULL_LITE + 1000 * st["real"] + i,
                              components=("white",))
        dnull = tuple(jnp.asarray(np.asarray(n)) for n in nz)
        dlnl, lnK, q_of, _ = sb.columns(theta_rec, led, dnull, ones,
                                        np.zeros(0, int), np.zeros(sb.npsr))
        o = np.where(q_of > QBAR, np.maximum(dlnl - np.maximum(lnK, 0.0), 0.0), 0.0)
        off[i] = float(np.max(o))
    p(f"      nulls: {n_null} white-only draws in {time.time()-t0:.0f}s")
    return off


def certify(dlnl, lnK, q_of, floor):
    """THE FROZEN GATE (criterion v2.2, glacier_loop.run_cell), verbatim. Do not modify.
    NOTE: the kill-step verdict is deliberately NOT an argument -- the kill step is
    SCORE-AND-LOG ONLY until its false-negative gate passes (capstone S4.23.1)."""
    return (dlnl > np.maximum(lnK + TRIALS_NAT, floor)) & (q_of > QBAR)


# ============================================================
# THE KILL STEP -- SCORE AND LOG ONLY. NEVER EXECUTES.
# ============================================================
def kill_score(st, theta_rec, led, cert_idx, q_of, verbose=True):
    """D2 rigidity (S4.21) R1/R2 at the certification's OWN state, scored and returned as
    columns. The caller banks it and NEVER acts on it."""
    if len(cert_idx) == 0:
        return dict(kill_scored=False, kill_note="no scorable claim (0 certs)")
    from spark import _adam_profile, H_ABSENT
    jax, jnp, sb, amo, nd = st["jax"], st["jnp"], st["sb"], st["amo"], st["nd"]
    fl_lnl, lb = amo["logL"], amo["logL_batch_theta"]
    data = st["data"]

    def lnl_at(th_base, kb, dat, pmask, free):
        th = th_base.at[kb + I_COSINC].set(free[0])
        th = th.at[kb + I_H].set(free[1])
        th = th.at[kb + I_PH0].set(free[2])
        th = th.at[kb + I_PSI].set(free[3])
        return fl_lnl(th, dat, pmask)
    prof = jax.jit(jax.vmap(_adam_profile(lnl_at), in_axes=(None, 0, None, None)))
    CH = 4

    fed = np.where(led.fed)[0]
    carried = np.setdiff1d(np.arange(len(st["G"]["slots"])), fed)
    th_state = theta_with_absent(theta_rec, nd, carried)
    sb.columns(theta_rec, led, data, st["ones"], np.zeros(0, int), np.zeros(sb.npsr))
    map_L = np.asarray(sb.FT.posterior(sb._last_LNL)["map_L"], float)

    def twoF(pmask, Ld):
        th = np.asarray(th_state, float).copy()
        th[sb.AI] = Ld
        kb_all = [nd + NP_SRC * int(k) for k in fed]
        thj, pmj = jnp.asarray(th), jnp.asarray(pmask)
        mx = []
        for c0 in range(0, len(kb_all), CH):
            kb = kb_all[c0:c0 + CH]
            KB = jnp.asarray(kb + [kb[0]] * (CH - len(kb)))
            m, _ = prof(thj, KB, data, pmj)
            mx.append(np.asarray(m)[:len(kb)])
        mx = np.concatenate(mx)
        offs = []
        for k in fed:
            t = th.copy(); t[nd + NP_SRC * int(k) + I_H] = H_ABSENT
            offs.append(t)
        ll_off = np.asarray(lb(jnp.asarray(np.stack(offs)), data, pmj))
        return 2.0 * (mx - ll_off)

    pm_coh = np.zeros(sb.npsr); pm_coh[cert_idx] = q_of[cert_idx]
    Ld_coh = sb.L0.copy(); Ld_coh[cert_idx] = map_L[cert_idx]
    tf_coh = twoF(pm_coh, Ld_coh)
    tf_s0 = twoF(np.zeros(sb.npsr), sb.L0)
    ks = int(np.argmax(tf_coh))
    c2F, d2F = float(tf_coh[ks]), float(tf_coh[ks] - tf_s0[ks])
    r1, r2 = bool(c2F >= BAR_R1), bool(d2F > 0.0)
    killed = not (r1 and r2)
    if verbose:
        print(f"      KILL STEP (score-and-log only): set n={len(cert_idx)} "
              f"2F_coh {c2F:.2f} (bar {BAR_R1}) Delta2F {d2F:+.2f} -> "
              f"would_be_{'KILLED' if killed else 'RETAINED'} "
              f"[NOT EXECUTED]", flush=True)
    return dict(kill_scored=True, kill_twoF_coh=tf_coh, kill_twoF_s0=tf_s0,
                kill_c2F=c2F, kill_d2F=d2F, kill_R1=r1, kill_R2=r2,
                kill_would_kill=killed, kill_argmax_member=ks,
                kill_mapL=map_L[cert_idx], kill_bar_r1=BAR_R1,
                kill_note="SCORE-AND-LOG ONLY; false-negative gate not passed "
                          "(capstone S4.23.1); never executed against any cert")


# ============================================================
# PART 3 -- THE PER-PULSAR AUTOPSY OF THE NEVER-CERTIFIED SET
# ============================================================
def pulsar_term_snr(st, k_loud=0):
    """Per-pulsar optimal SNR of the PULSAR TERM of census member `k_loud` (the loudest),
    against this venue's own white covariance:

        SNR_a^2 = (d_p,a | d_p,a) = sum_i (d_p,a[i]^2 / N_a[i])

    where d_p,a is the pulsar-term-only delay = (full two-term delay) - (Earth-term-only
    delay) for that member alone, everything else absent. This is the AUDIBILITY axis of
    the autopsy: how much matched power the array could ever harvest from that pulsar's
    own term, independent of whether the fringe can be identified."""
    jnp, amo, nd, sb = st["jnp"], st["amo"], st["nd"], st["sb"]
    slots = st["G"]["slots"]
    theta = st["theta_true"].copy()
    others = np.setdiff1d(np.arange(len(slots)), [k_loud])
    th = theta_with_absent(theta, nd, others)          # only k_loud present
    full = amo["inject_delay"](jnp.asarray(th), jnp.ones(sb.npsr))   # Earth + pulsar term
    earth = amo["inject_delay"](jnp.asarray(th), jnp.zeros(sb.npsr))  # Earth term only
    snr = np.zeros(sb.npsr)
    for a in range(sb.npsr):
        dp = np.asarray(full[a]) - np.asarray(earth[a])
        snr[a] = float(np.sqrt(np.sum(dp ** 2 / np.asarray(sb.sp.N_diag[a]))))
    return snr


def autopsy(st, ever_true, dlnl, lnK, q_of, snr_p):
    """PRE-STATED decomposition (fixed before scoring) of why a pulsar saturates
    UNCERTIFIED. The frozen criterion needs dlnL > max(lnK + 0.578, floor) AND q > 0.9.

      required   : dlnL_req  = lnK + TRIALS_NAT              (the trials/ambiguity bar)
      obtainable : dlnL_avail ~ SNR_p^2 / 2                  (the matched-filter ceiling
                   the pulsar term could deliver if its fringe were identified perfectly)

      AUDIBILITY-LIMITED : dlnL_avail <  dlnL_req  -- the pulsar term cannot clear its own
                           trials bar even with perfect fringe identification. More prior
                           does not help; only more signal or less noise does.
      AMBIGUITY-LIMITED  : dlnL_avail >= dlnL_req but the pulsar is uncertified -- the
                           power is there and the comb is what is unresolved. This is the
                           class a tighter prior CAN buy.
    Returns per-pulsar labels over the never-certified set."""
    dlnl_req = lnK + TRIALS_NAT
    dlnl_avail = 0.5 * snr_p ** 2
    aud = (~ever_true) & (dlnl_avail < dlnl_req)
    amb = (~ever_true) & (dlnl_avail >= dlnl_req)
    return dict(autopsy_dlnl_req=dlnl_req, autopsy_dlnl_avail=dlnl_avail,
                autopsy_snr_p=snr_p, autopsy_audibility_limited=aud,
                autopsy_ambiguity_limited=amb,
                n_audibility_limited=int(aud.sum()),
                n_ambiguity_limited=int(amb.sum()),
                autopsy_rule=("audibility-limited iff 0.5*SNR_p^2 < lnK + 0.578 (cannot "
                              "clear the trials bar even with perfect fringe ID); "
                              "ambiguity-limited otherwise. Stated before scoring."))


# ============================================================
# ONE LAP
# ============================================================
def run_loop(st, n_lap=N_LAP, n_null=N_NULL, stem=None, dry_stop=None,
             cond_tier=TIER_COND, verbose=True):
    """The closed loop. Per lap: Fisher -> frontier-v2 -> promote@truth -> drain ->
    E-step -> M-step -> floor re-cut -> criterion v2.2 -> kill score -> bank.

    `dry_stop`: stop early after this many CONSECUTIVE laps that add zero NEW true
    certifications (Part 3's saturation rule). None = run all n_lap laps.

    CUMULATIVE COUNT SEMANTICS: `ever_true` is the set of DISTINCT pulsars ever truly
    certified; a pulsar re-certified on a later lap is not a new certification. This is
    what makes the cumulative curve a compounding curve and what makes Part 3's
    completeness C = |ever_true| / npsr well defined."""
    p = print if verbose else (lambda *a, **k: None)
    sb, amo, nd, jnp = st["sb"], st["amo"], st["nd"], st["jnp"]
    ones, data = st["ones"], st["data"]
    slots = st["G"]["slots"]
    n_slot = len(slots)
    stem = stem or stem_loop(st["wkey"], st["real"], st.get("renorm", False))
    theta_rec, led, seed_idx, free_axes = condition_seeds(st, tier=cond_tier,
                                                          verbose=verbose)
    box_sigma = np.array([1.0, np.pi, 1.0, 0.5, 0.25, 0.25, np.pi, 0.5 * np.pi]) / np.sqrt(3.0)
    # The grid must bracket from ABOVE not just the population's A_equivalent but the
    # MISMATCH FIELD the background GP also absorbs when fed seeds are mis-specified.
    # -13 pegged at the CPU gates; -11 pegged again in smoke 12809208 (a_bg -11.000+-inf,
    # edge True) because prior-drawn seed axes leave a residual far larger than the
    # population itself. Widened to -8; `a_bg_edge` is banked every lap so a peg is
    # always visible rather than being read as a measurement.
    # ...and from BELOW: the T3 diagnostic (12810090) pegged at the -17 FLOOR, because
    # with all 16 members fed at truth the background is fully explained and the residual
    # GP amplitude runs off the bottom. The drain pegs in BOTH directions depending on
    # template quality, so the grid must be wide at both ends.
    a_grid = np.linspace(-20.0, -8.0, 241)
    prev_cert_idx, prev_q = np.zeros(0, int), np.zeros(sb.npsr)
    ever_true = np.zeros(sb.npsr, bool)
    dry = 0
    rows = []

    for lap in range(n_lap):
        ck = f"{OUT}/{stem}_lap{lap}__{lane_tag()}.npz"
        if os.path.exists(ck):
            z = np.load(ck, allow_pickle=True)
            led.fed = np.asarray(z["fed_mask"], bool).copy()
            theta_rec = np.asarray(z["theta_rec"], float).copy()
            prev_cert_idx = np.asarray(z["cert_idx"], int)
            prev_q = np.asarray(z["q_of_psr"], float)
            ever_true = np.asarray(z["ever_true_mask"], bool).copy()
            dry = int(z["dry_laps"])
            rows.append(dict(lap=lap, resumed=True))
            p(f"  lap {lap}: checkpoint exists -- resumed")
            continue
        t_lap = time.time()

        # (a) Fisher at the current joint state (own-term-live)
        carried = np.where(~led.fed)[0]
        sig_opt, sig_pes, F_ii = fisher_conditional(amo, jnp, theta_rec, carried, nd,
                                                    box_sigma, data, ones)
        ratio = np.max(sig_opt / sig_pes, axis=1)

        # (b) FRONTIER-v2 (mandatory): ratio < 0.5 AND dlnL_feed > 0
        new = [int(k) for k in np.where((ratio < GATE_RATIO) & ~led.fed)[0]]
        feed_dlnl = np.full(n_slot, np.nan)
        refused = []
        if new:
            th_off = theta_with_absent(theta_rec, nd, carried)
            for k in list(new):
                th_on = theta_with_absent(theta_rec, nd, np.setdiff1d(carried, [k]))
                ll = np.asarray(amo["logL_batch_theta"](
                    jnp.asarray(np.stack([th_on, th_off])), data, ones))
                feed_dlnl[k] = float(ll[0] - ll[1])
                if not feed_dlnl[k] > 0.0:
                    refused.append(k); new.remove(k)
            if refused:
                p(f"    frontier-v2 refused {refused} at the data-support term "
                  f"({[f'{feed_dlnl[k]:+.2f}' for k in refused]})")
        for k in new:
            led.promote(k, slots[k], iteration=lap)      # AT DRAWN TRUTH (deviation (c))
        fed_idx = np.where(led.fed)[0]
        n_new_resolved = len(new)

        # (c) THE DRAIN
        prof = st["bf"].profile(theta_rec, data, ones, led.fed.astype(float), a_grid)

        # (d)+(e) M-step: seeds on their T2-free axes, promoted on MSTEP_AXES
        n_eval = 0
        if len(fed_idx):
            carried_now = np.where(~led.fed)[0]
            th_eval = theta_with_absent(theta_rec, nd, carried_now)

            def marg_fn(srcs):
                ths = np.tile(th_eval, (len(srcs), 1))
                ths[:, nd:] = srcs
                return np.asarray(amo["logL_batch_theta"](jnp.asarray(ths), data, ones))

            src_cur = th_eval[nd:].copy()
            seeds_fed = np.array([k for k in fed_idx if k in set(seed_idx.tolist())], int)
            prom_fed = np.array([k for k in fed_idx if k not in set(seed_idx.tolist())], int)
            if len(seeds_fed) and len(free_axes):
                src_cur, w_s, ne = mstep_free(marg_fn, src_cur, seeds_fed, st["scale"],
                                              free_axes, n_sweep=MSTEP_SWEEPS_SEED,
                                              step0=MSTEP_STEP0_SEED)
                n_eval += ne
            if len(prom_fed):
                src_cur, w_p, ne = mstep_free(marg_fn, src_cur, prom_fed, st["scale"],
                                              MSTEP_AXES, n_sweep=MSTEP_SWEEPS_PROM,
                                              step0=MSTEP_STEP0_PROM)
                n_eval += ne
            theta_rec = theta_rec.copy()
            for k in fed_idx:
                theta_rec[nd + NP_SRC * k: nd + NP_SRC * (k + 1)] = \
                    src_cur[NP_SRC * k: NP_SRC * (k + 1)]

        # (f) FLOORS FIRST -- re-cut every lap (banked floors do not apply here)
        off = null_offenders(st, theta_rec, led, n_null=n_null, p=p)
        fl, err, est, zf = adopt_floor(off, st["FL"])

        # (g) criterion v2.2 -- PM cohered on the PREVIOUS lap's certified set (the loop)
        dlnl, lnK, q_of, on_true = sb.columns(theta_rec, led, data, ones,
                                              prev_cert_idx, prev_q)
        cert = certify(dlnl, lnK, q_of, fl)
        true_cert = cert & on_true
        wrong = cert & ~on_true
        cert_idx = np.where(cert)[0]
        fresh = true_cert & ~ever_true          # DISTINCT pulsars only
        new_true = int(fresh.sum())
        ever_true |= true_cert
        cum_true = int(ever_true.sum())
        dry = dry + 1 if new_true == 0 else 0

        # (h) THE KILL STEP -- scored, logged, NEVER executed
        kill = kill_score(st, theta_rec, led, cert_idx, q_of, verbose=verbose)

        areas = np.array([localisation_area_deg2(F_ii[k, [0, 1]]) for k in range(n_slot)])
        wall = time.time() - t_lap
        # PART 3 autopsy of the never-certified set (cheap; the SNR vector is static)
        if "snr_p" not in st:
            st["snr_p"] = pulsar_term_snr(st)
        aut = autopsy(st, ever_true, dlnl, lnK, q_of, st["snr_p"])
        pw = st.get("pw_rec") or {}
        bank_npz(stem + f"_lap{lap}",
                 prediction=PREDICTION,
                 prediction_renorm=PREDICTION_RENORM,
                 prediction_complete=PREDICTION_COMPLETE,
                 ever_true_mask=ever_true, dry_laps=dry,
                 completeness=float(cum_true) / sb.npsr,
                 renorm=bool(st.get("renorm", False)),
                 a_equivalent=st["G"]["a_equivalent"],
                 band_power_ratio=st["G"]["band_power_ratio"],
                 pop_band_power=st["G"]["pop_band_power"],
                 ng15_band_power=st["G"]["ng15_band_power"],
                 renorm_dex=st["G"]["renorm_dex"],
                 **{k: v for k, v in pw.items()},
                 **{k: v for k, v in aut.items()},
                 lap=lap, wkey=st["wkey"], white_ns=st["white_ns"],
                 real=st["real"], t_label=T_LITE, band_log10f=np.array(BAND_LITE),
                 n_cert=int(cert.sum()), n_true_cert=new_true,
                 n_true_cert_this_lap=int(true_cert.sum()),
                 cum_true_cert=cum_true, wrong_cert=int(wrong.sum()),
                 n_resolved=led.n_resolved(), n_new_resolved=n_new_resolved,
                 a_bg=prof["ahat"], a_bg_sig=prof["sig"], a_bg_grid=prof["grid"],
                 a_bg_lnl=prof["lnl"], a_bg_edge=prof.get("edge_hit", False),
                 dlnL_det=dlnl, lnK=lnK, qmax=q_of, on_true=on_true,
                 cert_idx=cert_idx, true_cert_idx=np.where(true_cert)[0],
                 q_of_psr=q_of, floor=fl, floor_err=err, floor_est=est,
                 zero_fraction=zf, null_offenders=off, n_null=n_null,
                 conc_ratio=ratio, sig_opt=sig_opt, sig_pes=sig_pes,
                 feed_dlnl=feed_dlnl, feed_refused=np.array(refused, int),
                 fed_mask=led.fed, promote_events=led.event_array(),
                 areas_deg2=areas, theta_rec=theta_rec, seed_idx=seed_idx,
                 tier_cond=cond_tier, tier_free=np.array([AXIS_NAME[j] for j in free_axes]),
                 n_mstep_eval=n_eval, wall_s=wall,
                 seed_noise=st["noise_seed"], seed_dist=st["dist_seed"],
                 seed_pop=SEED_POP_LITE, ncw=NCW_LITE,
                 population=np.array(POPULATION_LITE, float),
                 white_sigma_s=st["G"]["white_sigma_s"],
                 **_checksums(),
                 scope=("REPORT-ONLY. Promotion is AT DRAWN TRUTH (frozen-census rule): "
                        "recruited members enter the template at true parameters; only "
                        "the 3 seeds carry free axes. Kill step SCORE-AND-LOG ONLY."),
                 **{k: v for k, v in kill.items()})
        p(f"  lap {lap}: N_res {led.n_resolved()} (+{n_new_resolved}) "
          f"N_cert {int(cert.sum())} (true {new_true}, cum {cum_true}, "
          f"wrong {int(wrong.sum())}) A_bg {prof['ahat']:.3f}+-{prof['sig']:.3f} "
          f"floor {fl:.3f}({est}) zf {zf:.2f} [{wall:.0f}s]")
        rows.append(dict(lap=lap, n_res=led.n_resolved(), n_new=n_new_resolved,
                         n_cert=int(cert.sum()), n_true=new_true, cum_true=cum_true,
                         C=float(cum_true) / sb.npsr,
                         a_bg=float(prof["ahat"]), floor=fl, zf=zf, wall=wall))
        prev_cert_idx, prev_q = cert_idx, q_of
        if dry_stop is not None and dry >= dry_stop:
            p(f"  SATURATED: {dry} consecutive laps added 0 new true certifications "
              f"-- stopping at lap {lap} (C = {cum_true}/{sb.npsr} = "
              f"{cum_true/sb.npsr:.4f})")
            break
    return rows


# ============================================================
# PART 2 -- THE ANCHOR CELL (literature bridge)
# ============================================================
def anchor_cell(wkey=ANCHOR_WHITE, n_src=ANCHOR_N_SRC, real=0, verbose=True):
    """SINGLE-SHOT (no loop): 20 ns white, `n_src` sources conditioned at FULL truth (T3,
    every parameter), everything else carried. Reports per-pulsar distance posterior
    widths in BOTH house conventions, never merged:

      sigma_L_cond : within-fringe conditional width, 1/sqrt(-d2 logL / dL^2) at truth
                     with the CW parameters KNOWN -- the quantity the published
                     single-shot forecasts report, and the one comparable to their
                     sub-pc claims.
      sigma_L_marg : the fringe-comb width, sd of L over the fringe posterior -- the
                     honest number INCLUDING the real fringe ambiguity, which the
                     idealization does not remove.
    """
    p = print if verbose else (lambda *a, **k: None)
    check_affinity()
    st = build_state(wkey, real=real, verbose=verbose)
    sb, amo, nd, jnp = st["sb"], st["amo"], st["nd"], st["jnp"]
    ones, data = st["ones"], st["data"]
    slots = st["G"]["slots"]
    theta_rec = st["theta_true"].copy()          # T3: everything at truth
    led = PromoteLedger(slots)
    for k in range(min(n_src, len(slots))):
        led.promote(int(k), slots[int(k)], iteration=0)
    carried = np.where(~led.fed)[0]
    p(f"  anchor: {int(led.fed.sum())} sources fed at FULL truth (T3), "
      f"{len(carried)} carried; white {st['white_ns']:.0f} ns")

    # ---- sigma_L conditional: second differences of logL in each pulsar's distance
    lb = amo["logL_batch_theta"]
    th_base = theta_with_absent(theta_rec, nd, carried)
    dL = sb.dL
    sig_cond = {}
    for frac in (1e-2, 1e-3):
        step = frac * dL
        stack = []
        for a in range(sb.npsr):
            for d in (-1, 0, 1):
                t = th_base.copy(); t[sb.AI[a]] += d * step[a]
                stack.append(t)
        stack = np.stack(stack)
        lls = []
        for c0 in range(0, len(stack), 48):
            lls.append(np.asarray(lb(jnp.asarray(stack[c0:c0 + 48]), data, ones)))
        L = np.concatenate(lls).reshape(sb.npsr, 3)
        d2 = (L[:, 0] - 2.0 * L[:, 1] + L[:, 2]) / (step ** 2)
        F = -d2
        sig_cond[frac] = np.where(F > 0, 1.0 / np.sqrt(np.maximum(F, 1e-300)), np.inf)
    s_cond = sig_cond[1e-2]
    step_agree = np.abs(sig_cond[1e-2] / sig_cond[1e-3] - 1.0)

    # ---- the fringe posterior at the same state (the real ambiguity, kept real)
    dlnl, lnK, q_of, on_true = sb.columns(theta_rec, led, data, ones,
                                          np.zeros(0, int), np.zeros(sb.npsr))
    post = sb.FT.posterior(sb._last_LNL)
    map_L = np.asarray(post["map_L"], float)
    s_marg = np.zeros(sb.npsr)
    for a in range(sb.npsr):
        w = sb.FT.q[a]
        Lg = sb.L0[a] + sb.FT.uk[a] * dL[a]
        mu = float(np.sum(w * Lg))
        s_marg[a] = float(np.sqrt(max(np.sum(w * (Lg - mu) ** 2), 0.0)))

    L0 = sb.L0
    within = L0 < 1.0                              # the papers' "within ~1 kpc" cut
    KPC_PC = 1000.0
    p(f"\n  --- ANCHOR CELL: per-pulsar distance widths ({wkey}, {n_src} src at T3) ---")
    p(f"  pulsars within 1 kpc: {int(within.sum())}/{sb.npsr}")
    p(f"  sigma_L_cond (pc): median {np.median(s_cond)*KPC_PC:.4f}  "
      f"within-1kpc median {np.median(s_cond[within])*KPC_PC:.4f}  "
      f"sub-pc fraction {float((s_cond*KPC_PC < 1.0).mean()):.3f} "
      f"(within 1 kpc: {float((s_cond[within]*KPC_PC < 1.0).mean()):.3f})")
    p(f"  sigma_L_marg (pc): median {np.median(s_marg)*KPC_PC:.4f}  "
      f"within-1kpc median {np.median(s_marg[within])*KPC_PC:.4f}  "
      f"sub-pc fraction {float((s_marg*KPC_PC < 1.0).mean()):.3f}")
    p(f"  fringe: q_max median {np.median(q_of):.3f}; MAP == true fringe "
      f"{int(on_true.sum())}/{sb.npsr}; K_counted median "
      f"{np.median(sb.FT.K_counted):.0f}; dL median {np.median(dL)*KPC_PC:.2f} pc")
    p(f"  step-size stability |s(1e-2)/s(1e-3) - 1|: median {np.median(step_agree):.2e}, "
      f"max {np.max(step_agree):.2e}")

    path = bank_npz(f"gl_lite_anchor_{wkey}_n{n_src}_r{real}",
                    prediction=PREDICTION,
                    sigma_L_cond_kpc=s_cond, sigma_L_marg_kpc=s_marg,
                    sigma_L_cond_step1e3=sig_cond[1e-3], step_agree=step_agree,
                    L0_kpc=L0, dL_kpc=dL, within_1kpc=within,
                    q_max=q_of, on_true=on_true, map_L_kpc=map_L,
                    dlnL_det=dlnl, lnK=lnK, K_counted=sb.FT.K_counted,
                    psr_names=np.array([pe.name for pe in st["G"]["ent_psrs"]]),
                    n_src_fed=int(led.fed.sum()), tier_cond="T3_full_truth",
                    wkey=wkey, white_ns=st["white_ns"], real=real,
                    t_label=T_LITE, ncw=NCW_LITE, seed_pop=SEED_POP_LITE,
                    seed_noise=st["noise_seed"], seed_dist=st["dist_seed"],
                    **_checksums(),
                    note=("PART 2 ANCHOR CELL: single-shot, no loop. sigma_L_cond is the "
                          "within-fringe conditional width with CW params KNOWN (the "
                          "literature-comparable quantity); sigma_L_marg is the "
                          "fringe-comb sd, which retains the REAL fringe ambiguity. "
                          "Both banked, never merged."))
    p(f"  banked -> {path}")
    return 0


# ============================================================
# SMOKE (1 cell, mandatory before any array submit)
# ============================================================
def smoke(wkey="w50", cond_tier=TIER_COND, verbose=True):
    """One venue build + gates + ONE lap. Proves plumbing before anything wide runs."""
    global TAG
    TAG = "_smoke" if cond_tier == TIER_COND else f"_diag{cond_tier}"
    check_affinity()
    print(f"=== GLACIER-LITE SMOKE ({wkey}) ===", flush=True)
    print(f"PREDICTION (pre-logged, verbatim):\n  {PREDICTION}\n", flush=True)
    t0 = time.time()
    st = build_state(wkey, real=0, verbose=verbose)
    print(f"  venue+state built in {time.time()-t0:.0f}s", flush=True)
    rows = run_loop(st, n_lap=1, n_null=32, cond_tier=cond_tier, verbose=verbose)
    print(f"  SMOKE rows: {rows}", flush=True)
    print(f"  total {time.time()-t0:.0f}s", flush=True)
    return 0


def gates(wkey="w50", renorm=False, pw_factor=None, verbose=True):
    """BUILD-ONLY gates -- no E-step, no nulls, so this runs on the CPU lane (which needs
    no GPU entitlement) and fast-fails the riskiest new code before a GPU slot is spent:
      1. frozen-population gate   (log10_h ladder == the banked POP)
      2. white-noise gate         (likelihood N_diag rms == target: the idealization
                                   reached the COVARIANCE, not just the data)
      3. banked T=15 L_gwb loads STRICT against the rebuilt venue (the GWB block is
                                   untouched by the white-noise change)
      4. BackgroundFit constructs (>=1 in-band GP mode; both identity checks pass)
      5. the drain profile is finite at nothing-fed
    Deliberately does NOT touch estep_per_target (the XLA-CPU map_count hazard)."""
    print(f"=== GLACIER-LITE BUILD GATES ({wkey}) ===", flush=True)
    print(f"PREDICTION (pre-logged, verbatim):\n  {PREDICTION}\n", flush=True)
    t0 = time.time()
    st = build_state(wkey, real=0, renorm=renorm, pw_factor=pw_factor, verbose=verbose)
    amo, jnp, sb = st["amo"], st["jnp"], st["sb"]
    bf = st["bf"]
    if pw_factor is not None:
        pw = st["pw_rec"]
        K = pw["K_counted"]
        floor_ok = bool(np.all(np.asarray(pw["sigma_lit_kpc"])
                               >= PW_FRINGE_FLOOR * sb.dL - 1e-12))
        print(f"  PW GATE sigma >= {PW_FRINGE_FLOOR:.0f}*dL on all {sb.npsr} pulsars: "
              f"{'PASS' if floor_ok else 'FAIL'}  (floored on {pw['pw_n_floored']})",
              flush=True)
        print(f"  K distribution: min {K.min():.0f} median {np.median(K):.0f} "
              f"max {K.max():.0f}  (K at the floor would be "
              f"{2*int(3*PW_FRINGE_FLOOR)})", flush=True)
    snr = pulsar_term_snr(st)
    print(f"  pulsar-term SNR (loudest member): min {snr.min():.2f} "
          f"median {np.median(snr):.2f} max {snr.max():.2f}", flush=True)
    print(f"  in-band GP modes at band {BAND_LITE}: {int(bf.inband.sum())} of "
          f"{bf.ngp} (n_scaled {bf.n_scaled}) -> "
          f"{'PASS' if bf.inband.any() else 'FAIL'}", flush=True)
    ok = bf.gate_scaling_identity()
    print(f"  BackgroundFit scaling identity: {'PASS' if ok is None or ok else ok}",
          flush=True)
    theta = st["theta_true"]
    prof = bf.profile(theta, st["data"], st["ones"], np.zeros(len(st["G"]["slots"])),
                      np.linspace(-16.0, -13.0, 61))
    fin = bool(np.all(np.isfinite(prof["lnl"])))
    print(f"  drain profile @ nothing-fed: finite {fin}; ahat {prof['ahat']:.3f}"
          f"+-{prof['sig']:.3f} edge {prof.get('edge_hit')} -> "
          f"{'PASS' if fin else 'FAIL'}", flush=True)
    print(f"  fringe grid: dL median {np.median(sb.dL)*1000:.2f} pc, K_counted median "
          f"{np.median(sb.FT.K_counted):.0f}, npsr {sb.npsr}", flush=True)
    print(f"  ALL BUILD GATES DONE in {time.time()-t0:.0f}s", flush=True)
    return 0


def run_complete(pw_factor, real, n_lap=COMPLETE_LAPS, n_null=N_NULL, verbose=True):
    """PART 3 -- one completeness cell: 20 ns, Part-1 conditioning (3 loudest at T2),
    prior widths scaled by `pw_factor` with the 3*dL hard floor, run to saturation."""
    check_affinity()
    tag = pw_tag_of(pw_factor)
    print(f"=== GLACIER-LITE PART 3 (completeness) pw=1/{tag} r{real} ===", flush=True)
    print(f"PREDICTION (pre-logged, verbatim):\n  {PREDICTION_COMPLETE}\n", flush=True)
    st = build_state(COMPLETE_WHITE, real=real, pw_factor=pw_factor, verbose=verbose)
    rows = run_loop(st, n_lap=n_lap, n_null=n_null, stem=stem_complete(tag, real),
                    dry_stop=DRY_LAPS_STOP, verbose=verbose)
    if rows:
        last = [r for r in rows if "C" in r]
        if last:
            r = last[-1]
            print(f"  CELL DONE: saturation lap {r['lap']}, C = {r['C']:.4f} "
                  f"({r['cum_true']}/116)", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["smoke", "loop", "anchor", "gates", "complete"])
    ap.add_argument("--white", default="w50", choices=list(WHITE_LEVELS_NS))
    ap.add_argument("--real", type=int, default=0)
    ap.add_argument("--laps", type=int, default=N_LAP)
    ap.add_argument("--nnull", type=int, default=N_NULL)
    ap.add_argument("--nsrc", type=int, default=ANCHOR_N_SRC)
    ap.add_argument("--renorm", action="store_true",
                    help="SECONDARY ARM: renormalise the population's summed band power "
                         "to NG15 (A = -14.6) by one common dex shift on every strain.")
    ap.add_argument("--cond", default=TIER_COND, choices=("T1", "T2", "T3"),
                    help="conditioning tier for the seeds. T3 (all 8 axes at truth) is "
                         "the INSTRUMENT DIAGNOSTIC: it removes seed mis-specification "
                         "entirely, so a persistent on_true=0 under T3 is structural, "
                         "whereas recovery under T3 localises the effect to the "
                         "free-axis M-step.")
    ap.add_argument("--pw", type=float, default=None,
                    help="PART 3: prior-width factor (1, 0.3333, 0.1, 0.03333).")
    a = ap.parse_args()
    if a.mode == "gates":
        return gates(a.white, renorm=a.renorm, pw_factor=a.pw)
    if a.mode == "smoke":
        return smoke(a.white, cond_tier=a.cond)
    if a.mode == "anchor":
        return anchor_cell(a.white, n_src=a.nsrc, real=a.real)
    if a.mode == "complete":
        if a.pw is None:
            raise SystemExit("mode 'complete' requires --pw")
        return run_complete(a.pw, a.real, n_lap=a.laps, n_null=a.nnull)
    check_affinity()
    arm = "SECONDARY (NG15-renormalised)" if a.renorm else "PRIMARY (frozen POP)"
    print(f"=== GLACIER-LITE LOOP {a.white} r{a.real} -- {arm} ===", flush=True)
    print(f"PREDICTION (pre-logged, verbatim):\n  {PREDICTION}\n", flush=True)
    if a.renorm:
        print(f"SECONDARY-ARM PREDICTION (pre-logged, verbatim):\n"
              f"  {PREDICTION_RENORM}\n", flush=True)
    st = build_state(a.white, real=a.real, renorm=a.renorm, verbose=True)
    run_loop(st, n_lap=a.laps, n_null=a.nnull, verbose=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
