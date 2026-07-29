"""GLACIER Stage 1 -- THE LOOP (agent GLACIER, ACCRE). Does the array EAT THE BACKGROUND?

*** HELD: DO NOT SUBMIT THE FAN UNTIL FORGE-G2's runtime-SMASK LANDS AND RE-GATES ***
(the capstone g0 note: baked-SMASK pays 17-32 GPU-hr of pure recompile at campaign scale;
the re-gate expects all Bg5 PASS + bit-exact pattern flip max|diff|=0.0 + warm flip <10ms.
Mode `gate` here is runnable now -- it is one cheap zero-noise iteration, not the campaign.)

THE DESIGN (pre-registered; deviations from SPARK-3 SS5.3 are declared in the capstone SSL.1)
--------------------------------------------------------------------------------------------
Fixed-T primary: T = 30, lit priors, 8 sky draws. TWO IGNITERS (Matt's budget decision,
Option 1), PAIRED on census + noise seeds (the CHORUS C6 discipline -- the igniter contrast
is exact):
    e07  : the census's BRIGHTEST member at e = 0.7, embedded as its WEAVE-3.3 harmonic
           comb (CHORUS C1-C3 conventions via atlas_harmonics; 31 appended slots).
    none : the same census all-circular. The negative control.
SCOPE LINE (travels with every number): the drawn population's own eccentricities (the
LOTTERY realistic point, banked per source) are REPRESENTED CIRCULARLY in both data and
template for every member except the igniter -- eccentricity enters through the igniter
only (the dynamics experiment). 6 iterations to a fixed point; EXTENSION RIDER to 12 from
checkpoint for any cell whose DeltaN or drain still moves beyond its own error; ESCALATION
RIDER (e=0.5, e=0.3x2 arms) is Matt's call at first-sky readback, not before.

ONE SOFT-JOINT ITERATION (everything posterior-weighted, nothing hard-committed):
  (a) FISHER at the current joint state: per-census-member conditional widths (OPTIMISTIC
      bound, capped at the prior box; the PESSIMISTIC bound is the box itself -- both
      banked as columns, never merged).
  (b) THE FRONTIER: ModelQualityGate (ratio 0.5, DECLARED) on the optimistic width ->
      target fed set; new crossings promoted through the PromoteLedger AT DRAWN TRUTH
      (the frozen-census rule; a fitted-position promote raises CampaignStop = B0.2).
  (c) THE DRAIN: BackgroundFit profile of the background GP amplitude at the new smask --
      A_background(iter), the headline curve.
  (d) E-STEP: per-target fringe posteriors (spark3.estep_per_target; PM from the PREVIOUS
      iteration's certified set at their q -- uncertified DECOHERED, never pinned).
  (e) M-STEP: batched quadratic sweeps of lnL_marg over the fed members' (log10_fgw,
      log10_mc) scaled offsets (the TargetedMarg surface with DYNAMIC slots -- slot
      selection is numpy-side, so a frontier change recompiles nothing). Fed widths
      re-measured (Laplace) from the same curvatures.
  (f) END-OF-CYCLE SCOREBOARD (criterion v2.2, NEVER inside the loop): per-pulsar
      cert = (dlnL_det > max(lnK, floor)) & (q_max > 0.9); trials cost restored
      (+0.578 nat, saturating, SPARK-2 SS2) via lnK over the CARRIED census.
  (g) BANK per (sky, igniter, iter): N_cert, N_resolved, A_background +- sig, per-component
      concentration + both bounds, localisation areas (EPOCH crossings report-only),
      channel budget, wrong-certs (must stay 0), floor + zero-fraction + estimator label,
      host/thread metadata. Checkpoint -> resume.

FLOOR POLICY (Matt's decision on the held line): FULL per-iteration floors on the
verdict-carrying subset -- both scrambled-null cells + skies {0,1} per igniter (the
calibration pair). Remaining skies: floors REFIT at iterations {0,3,6}, carried between,
with the pre-registered DRIFT GATE: a refit floor moving beyond its bootstrap band from the
carried value escalates that sky to full refits from checkpoint. Floors per ANCHOR SS4
(recut_surface.adopt: Gumbel only at zf <= 0.20, else emp-q95 + bootstrap sd; zf a REQUIRED
column). Scrambled-counterpart null (forge_b1.scramble_source_sky on the igniter, data at
truth) through the FULL loop at 2 skies x >= 5 reals: the drain must not manufacture --
STOP + anatomy on any loop-grown wrong cert, any null-arm promote (spurious resolution),
or a null-arm drain fall beyond its own error.

VERDICTS (pre-registered per arm): COMPOUNDING (gains grow/hold >= 3 iterations; report the
saturator) / CONVERGENT (fixed point; height above ignition + shortfall to next recruit) /
INERT (nothing past iter 0; nearest-miss margins). Plus the DRAIN verdict: does
A_background measurably fall, with error. Scope lines everywhere.

SEEDS (house convention, disjoint): census 9_0xx_xxx (9_000_000 + sky); noise 9_5xx_xxx
(9_500_000 + 1000*sky + real); dist = noise + 10_000_000; scrambles 9_8xx_xxx.

Matt commits; never commit from here. Env: batch_gpu H200 lane, cpus-per-task=8 enforced,
autotune off, x64, SHARED warm campaign cache (TURBINE SS7).
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

from glacier_pop import (draw_population, draw_population_conditional, LADDER_RUNGS,
                         source_band_power, a_eff_projection,
                         PromoteLedger, BackgroundFit, CampaignStop, bank_npz, lane_tag,
                         check_affinity, build_glacier_problem, OUT,
                         N_POP_DEFAULT, SEED_POP_BASE, A_TARGET_LOG10,
                         I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI, NP_SRC)

# ---- campaign constants (pre-registered) ----
T_LABEL = 30
TIER = "lit"
# REMEDY A (authorized 2026-07-23): the generative band extends to ~2 nHz so the drain is a
# measurement, not an edge (GLACIER_capstone.md SSII/V.1 -- the incumbent 10-32 nHz box is
# outside the array's hearing; certified by the g2a-ii re-gate at this band). Every campaign
# draw and every BackgroundFit in this driver uses THIS band; banked per cell.
BAND_CAMPAIGN = (-8.7, -7.5)
# The canonical GWB square roots per venue (host-freedom fix; cut by gl_lgwb_t30.py at
# cpus=8 on hgx03, cross-host gated by gl_lgwb_xhost.py). strict=True -- a missing bank
# is a CRASH, never a silent thread/host-dependent recompute. T=40 added for the
# STAGE-2a ARRAY LADDER (2026-07-24 addendum).
LGWB_BANKS = {30: "/home/mattm/projects/HSYMT/CW_transition/b1_L_gwb_T30.npz",
              40: "/home/mattm/projects/HSYMT/CW_transition/b1_L_gwb_T40.npz"}
SKIES = list(range(8))
N_ITER, N_ITER_EXT = 6, 12
QBAR = 0.9
GATE_RATIO = 0.5                     # DECLARED (SPARK3_GATE_RATIO); the one re-derivable knob
TRIALS_NAT = 0.578                   # SPARK-2 SS2, saturating; restored to the bar
N_NULL_FULL = 100                    # nulls per full-floor refit (n_null=100 convention)
FULL_FLOOR_SKIES = (0, 1)            # the calibration pair, per igniter
CARRY_REFIT_ITERS = (0, 3, 6)        # carried-floor skies refit here; drift gate in between
EPOCH_AREA_DEG2 = 1.0                # host-unique localisation crossing, REPORT-ONLY, provisional
ARMS = ("e07", "none")
E_IGNITER = {"e07": 0.7, "none": 0.0}
FISHER_AXES = (I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI)
MSTEP_AXES = (I_FGW, I_MC)           # the TargetedMarg surface, per fed member
SEED_NOISE_BASE = 9_500_000
SEED_SCRAM_BASE = 9_800_000


def cell_stem(arm, sky, kind="cell"):
    return f"gl1_{kind}_{arm}_s{sky}_T{T_LABEL}_{TIER}"


# ============================================================
# the stack (jax config BEFORE heavy imports -- spark3._stack discipline)
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
    import trackB_b1 as B1
    import trackB_estimator as TE
    import ignite as IG
    import forge_b1 as F
    from recut_surface import gumbel_floor, emp_quantile, boot_sd, adopt, offender
    return jax, jnp, C, B1, TE, IG, F, dict(gumbel_floor=gumbel_floor,
                                            emp_quantile=emp_quantile, boot_sd=boot_sd,
                                            adopt=adopt, offender=offender)


# ============================================================
# igniter embedding (CHORUS C1-C3, via atlas_harmonics)
# ============================================================
TERM_FLOOR = 0.2                     # CHORUS/ATLAS merger guard: t_coal > (1-0.2)*span
TSUN_S = 4.925490947e-6              # discovery const.Tsun (s)
# venue spans, measured from the build logs (T=15: "span 22.15 yr"; T=30: "-> 37.14 yr").
# ASSERTED against the built problem's own getspan in run_cell/mode_gate -- a feather-set
# change fails loudly, never silently mis-guards the comb.
YR_S = 365.25 * 86400.0
VENUE_SPAN_S = {15: 22.15 * YR_S, 30: 37.14 * YR_S,
                40: 47.17 * YR_S}   # T=40 MEASURED (bank-cut job 12741097 build print:
                                    # "EXTENSION dT=25yr: span 22.15->47.17 yr")


def embed_igniter(pop, e_char, tmax_s):
    """Return (src_slots (n_slot,8), n_harm_appended, active_mask, chan_budget, n_clip).
    At e=0 the census is unchanged BIT-EXACTLY (n=2 harmonic IS the circular source --
    CHORUS C1's structural collapse). At e>0 the BRIGHTEST member's own slot carries n=2
    and its other 31 harmonics are APPENDED slots with the WEAVE-3.3 tie (mc_n, fgw_n,
    h_n; sky/inc/phase/psi shared), under the MERGER GUARD (chorus.tie_member verbatim:
    mc_n capped so every kernel's Earth term is evaluable over (1-TERM_FLOOR)*span --
    the WEAVE tie shortens t_coal by 2^(8/3)*F(e), ~173x at e=0.7, and an uncapped comb
    injects already-coalesced kernels whose waveform is NaN: SPARK-3 SS2's pathology in
    the DATA. Found by driver-gate NaN, job 12706069). n_clip banked per cell (CHORUS
    discipline: a fully-clipped comb is TOY-INVALID, flagged, never silently read)."""
    src = pop["src"].copy()
    if e_char <= 0.0:
        return src, 0, np.ones(len(src), bool), int(len(src)), 0
    import atlas_harmonics as H
    m = src[0]                                   # the brightest member (census rank 0)
    LOG2 = np.log10(2.0)
    f_orb = m[I_FGW] - LOG2
    ns_all = np.arange(1, 33)
    gvals = np.array([H.g_n(n, e_char) for n in ns_all])
    gmax = gvals.max()
    mc_all = np.array([H.tie_log10_mc_n(m[I_MC], e_char, n) for n in ns_all])
    fgw_all = np.array([H.tie_log10_fgw_n(f_orb, n) for n in ns_all])
    w_all = np.pi * 10.0 ** fgw_all
    cap = 0.6 * np.log10((1.0 - TERM_FLOOR) * 5.0
                         / (256.0 * w_all ** (8.0 / 3.0) * tmax_s)) - np.log10(TSUN_S)
    n_clip = int(np.sum(mc_all > cap))
    mc_all = np.minimum(mc_all, cap)
    harm, active = [], [True] * len(src)
    n_active_harm = 0
    for n in ns_all:
        g = max(gvals[n - 1], H.G_FLOOR)
        row = m.copy()
        row[I_MC] = mc_all[n - 1]
        row[I_FGW] = fgw_all[n - 1]
        row[I_H] = H.tie_log10_h_n(m[I_H], e_char, n, gvals=gvals[n - 1])
        is_active = g >= 1e-3 * gmax             # CHORUS C3's active-harmonic rule
        if n == 2:
            src[0] = row                          # the member's own slot carries n=2
        else:
            harm.append(row)
            active.append(bool(is_active))
        if is_active:
            n_active_harm += 1
    slots = np.vstack([src, np.array(harm)])
    chan = (len(src) - 1) + n_active_harm         # circular members + the comb's active lines
    return slots, len(harm), np.array(active, bool), int(chan), n_clip


# ============================================================
# per-iteration pieces
# ============================================================
def fisher_conditional(amo, jnp, theta_rec, carried, nd, box_sigma, data, pmask,
                       chunk=48):
    """Per-(census member, axis) conditional widths at the CURRENT joint state (central
    second differences on logL). OWN-TERM-LIVE (the G-d2a rung-mask principle, fixed
    2026-07-24): member k is measured on a template where K ITSELF IS PRESENT at its
    recovery params and the OTHER carried members are absent (they live in the fitted
    GP) -- the conditional question is 'how tight would k's widths be if k ENTERED the
    joint at this state'. The pre-fix code absented k too, so d2logL/dtheta_k^2 == 0
    identically and the frontier could never fire at ANY brightness -- caught by the
    r13p25 ladder cells (a 16x-background source reading ratio 1.000; driver gate G-d4
    is the standing positive control against regression). Capped at the prior box (the
    correct 'we know nothing'). Returns (n_slot, n_axis) sig_opt, sig_pes, F_ii."""
    lb = amo["logL_batch_theta"]
    n_slot = (len(theta_rec) - nd) // NP_SRC
    carried = set(int(k) for k in carried)
    box = np.array([box_sigma[j] for _ in range(n_slot) for j in FISHER_AXES])
    step = 1e-3 * box
    stack = []
    for k in range(n_slot):
        # member k's own evaluation template: k live, other carried absent
        others = np.array(sorted(carried - {k}), dtype=int)
        th_k = theta_with_absent(theta_rec, nd, others)
        for a, j in enumerate(FISHER_AXES):
            i = k * len(FISHER_AXES) + a
            col = nd + NP_SRC*k + j
            for d in (-1, 0, 1):
                t = th_k.copy(); t[col] += d * step[i]
                stack.append(t)
    stack = np.stack(stack)
    lls = []
    for c0 in range(0, len(stack), chunk):
        lls.append(np.asarray(lb(jnp.asarray(stack[c0:c0+chunk]), data, pmask)))
    L = np.concatenate(lls).reshape(n_slot * len(FISHER_AXES), 3)
    d2 = (L[:, 0] - 2.0*L[:, 1] + L[:, 2]) / (step**2)
    F_ii = -d2
    cond = np.where(F_ii > 0, 1.0/np.sqrt(np.maximum(F_ii, 1e-300)), np.inf)
    sig_opt = np.minimum(np.where(np.isfinite(cond), cond, box), box)
    return (sig_opt.reshape(n_slot, len(FISHER_AXES)),
            box.reshape(n_slot, len(FISHER_AXES)), F_ii.reshape(n_slot, len(FISHER_AXES)))


H_ABSENT_GL = -30.0
# GLACIER's absence strain is -30, NOT spark3's -18: G-d3 (job 12706109) MEASURED the -18
# remnant at 3.5e-4 relative to the population signal (max|d FtNmy| 1.4e4 of 4e7; d p0
# 9.1e-3 nat) -- harmless against every campaign margin (the smallest is the 0.578-nat
# trials cost) but NOT the bit-comparable identity the gate demands. At -30 the remnant is
# ~5e-16 relative: absence to fp. The waveform stays evaluable (h scales linearly).


def theta_with_absent(theta, nd, carried_idx):
    """Template-side statement of 'carried, not fed' on the raw logL / scoreboard paths:
    carried slots' strain to H_ABSENT_GL. The BackgroundFit path uses SMASK proper; the
    two statements are gated equal in driver gate G-d3."""
    t = theta.copy()
    for k in carried_idx:
        t[nd + NP_SRC*k + I_H] = H_ABSENT_GL
    return t


def mstep_quadratic(marg_fn, src0, fed_idx, scale, n_sweep=2, step0=0.3, bmax=64):
    """Batched coordinate-quadratic ascent of lnL_marg over the fed members' MSTEP_AXES.
    All slot selection is numpy-side (no recompile on frontier change). Returns
    (src_hat, widths (n_fed, n_axes), n_eval)."""
    src = src0.copy()
    slots = np.array([NP_SRC*k + j for k in fed_idx for j in MSTEP_AXES])
    sc = np.array([scale[j] for _ in fed_idx for j in MSTEP_AXES])
    n_eval = 0
    widths = np.full(len(slots), np.nan)
    for sweep in range(n_sweep):
        for i, sl in enumerate(slots):
            d = step0 * sc[i] * (0.5 ** sweep)
            trip = np.tile(src, (3, 1))
            trip[0, sl] -= d; trip[2, sl] += d
            v = np.asarray(marg_fn(trip)); n_eval += 3
            a, b, c = v
            den = a - 2*b + c
            if den < 0:                            # concave: quadratic step, clipped
                off = np.clip(0.5*(a - c)/den, -1.0, 1.0) * d
                src[sl] += off
                widths[i] = d / np.sqrt(-den)
            # non-concave: leave the slot (carried at entry value this sweep)
    return src, widths.reshape(len(fed_idx), len(MSTEP_AXES)), n_eval


def localisation_area_deg2(F_ii_sky):
    """Report-only 1-sigma sky-ellipse area from the two conditional sky curvatures
    (diagonal approximation, declared as such): pi * sig1 * sig2 in steradian-equivalent
    (cos_gwtheta x gwphi), converted to deg^2."""
    s1, s2 = F_ii_sky
    if s1 <= 0 or s2 <= 0:
        return np.inf
    return float(np.pi / np.sqrt(s1*s2) * (180.0/np.pi)**2)


# ============================================================
# THE CELL -- one (arm, sky), all iterations in-process
# ============================================================
def run_cell(arm, sky, n_src=N_POP_DEFAULT, scrambled=False, real=0, n_iter=N_ITER,
             full_floors=None, verbose=True, rung=None, t_label=None, wscale=1.0):
    """The campaign cell. Checkpointed per iteration; resume = skip-on-exist readback.
    scrambled=True: the counterpart-scrambled null arm (igniter sky scrambled in the
    RECOVERY template; data at truth) -- the manufacturing control, run through the FULL
    loop. STOPs (CampaignStop) on: null-arm promote, null-arm wrong cert, frozen-census
    violation.
    rung: a LADDER_RUNGS key -> STAGE-2a conditional-sky cell (Option C, authorized
    2026-07-24): the census is the conditional draw for (rung, sky); the igniter arm is
    a PROPERTY OF THE CONDITIONED LOUDEST MEMBER (e07 ties it into the eccentric comb --
    at the feasible rung this doubles as the probe of structure-assisted onset below the
    measured box); stems gl2_<rung>_*; the conditioning record is banked per iteration."""
    check_affinity()
    p = print if verbose else (lambda *a, **k: None)
    if full_floors is None:
        full_floors = scrambled or (sky in FULL_FLOOR_SKIES)
    jax, jnp, C, B1, TE, IG, F, FL = _stack()

    # STAGE-2a array ladder (2026-07-24 addendum): t_label + wscale are venue axes of
    # the ladder's ARRAY side. Defaults reproduce the incumbent stems bit-for-bit.
    T = T_LABEL if t_label is None else int(t_label)
    if T not in LGWB_BANKS:
        raise CampaignStop(f"no canonical L_gwb bank declared for T={T}")
    wtag = "" if wscale == 1.0 else "_w" + f"{wscale:g}".replace(".", "p")
    kind = f"null{real}" if scrambled else "cell"
    stem = cell_stem(arm, sky, kind) if rung is None else \
        f"gl2_{rung}{wtag}_{kind}_{arm}_s{sky}_T{T}_{TIER}"
    p(f"[GLACIER-{'1' if rung is None else '2'}] cell {stem} full_floors={full_floors} "
      f"wscale={wscale:g} lane={lane_tag()}")

    # ---- census + igniter embedding (paired arms share the census seed) ----
    if rung is None:
        pop = draw_population(SEED_POP_BASE + sky, n_src=n_src, band_log10f=BAND_CAMPAIGN)
        cond = None
    else:
        pop, cond = draw_population_conditional(sky, rung, n_src=n_src,
                                                band_log10f=BAND_CAMPAIGN)
        p(f"  conditioning [{cond['cond_mode']}]: brightest {cond['h_brightest']:+.3f} "
          f"(scan {cond['n_scanned']}), excess x{cond['excess_power_ratio']:.2f}, "
          f"A_eq {cond['a_equiv_log10']:+.3f}, tension {cond['tension_sigma']:+.1f}sigma")
    # amendment 2026-07-23: A_eff(drawn) is definitional (the drain's reference) -- banked
    a_eff_drawn, _ = a_eff_projection(pop["src"], band_log10f=BAND_CAMPAIGN)
    slots, n_harm, active, chan, n_clip = embed_igniter(pop, E_IGNITER[arm],
                                                        VENUE_SPAN_S[T])
    n_slot = len(slots)
    p(f"  census n={n_src} (+{n_harm} harmonic slots), channel budget {chan}, "
      f"merger-guard n_clip {n_clip}")

    # ---- problem + data ----
    pop_slots = dict(pop); pop_slots["src"] = slots; pop_slots["n_src"] = n_slot
    G = build_glacier_problem(T, pop_slots, verbose=verbose)
    if abs(G["span_s"] / VENUE_SPAN_S[T] - 1.0) > 0.02:
        raise CampaignStop(f"built span {G['span_s']:.4e}s differs >2% from the declared "
                           f"VENUE_SPAN_S[{T}] = {VENUE_SPAN_S[T]:.4e}s -- "
                           f"the merger guard was cut against a wrong span; STOP.")
    G["slots"] = slots
    amo = G["amo"]
    nd = amo["n_dist"]
    ones = jnp.ones(amo["npsr"])
    noise_seed = SEED_NOISE_BASE + 1000*sky + real
    dist_seed = noise_seed + 10_000_000            # house convention (module header)
    # G-d2 scoreboard machinery (also carries sp, needed by the noise drawer)
    sb = CertScoreboard(G, amo, jnp, C, prior_key=TIER)
    L_true, n_true = sb.draw_truth(IG, dist_seed, tier=TIER)
    theta_true = np.asarray(amo["theta_truth"], float).copy()
    theta_true[sb.AI] = L_true                     # data at the DRAWN distance truth
    clean = amo["inject_delay"](jnp.asarray(theta_true), ones)
    # NO stochastic GWB in the data -- the background IS the unresolved sum. Noise is
    # white+RN only (components below). The drawer loads the CANONICAL T=30 L_gwb bank,
    # STRICT (no silent recompute): with L frozen, draw(seed) is HOST-FREE up to BLAS
    # matmul reduction order (gated <1e-10 by gl_lgwb_xhost.py on the A100 lane) -- the
    # eigh basis freedom that forced host-pinning is gone (cut 2026-07-23, Matt's call:
    # the physics is hardware-invariant; freeze the artifact, not the host).
    ndraw = C.NoiseDrawer(sb.sp, amo, lgwb_path=LGWB_BANKS[T], strict=True)
    p(f"  L_gwb: {ndraw.lgwb_prov}")
    noise = ndraw.draw(noise_seed, components=("white", "rn"), white_scale=wscale)
    data = tuple(jnp.asarray(np.asarray(c) + np.asarray(n)) for c, n in zip(clean, noise))

    # ---- recovery template state (scramble applies HERE, never to the data) ----
    theta_rec = theta_true.copy()
    if scrambled:
        rng_s = np.random.default_rng(SEED_SCRAM_BASE + 1000*sky + real)
        theta_rec[nd + I_COSGT] = rng_s.uniform(-1, 1)     # igniter slot 0 sky, isotropic
        theta_rec[nd + I_GWPHI] = rng_s.uniform(0, 2*np.pi)

    led = PromoteLedger(slots)
    bf = BackgroundFit(amo, band_log10f=BAND_CAMPAIGN)
    box_sigma = np.array([1.0, np.pi, 1.0, 0.5, 0.25, 0.25, np.pi, 0.5*np.pi]) / np.sqrt(3.0)
    a_grid = np.linspace(A_TARGET_LOG10 - 1.0, A_TARGET_LOG10 + 1.0, 41)
    scale = TE.phi_scale({"ncw": 1})     # one 8-param block; mstep indexes by axis
    carried_floor = None
    prev_abg = prev_sig = prev_nfed = None   # null drain-fall STOP state (amended semantics)
    prev_cert_idx, prev_q = np.zeros(0, int), np.zeros(amo["npsr"])

    for it in range(n_iter):
        ck = f"{OUT}/{stem}_i{it}__{lane_tag()}.npz"
        if os.path.exists(ck):
            z = np.load(ck, allow_pickle=True)
            led.fed = np.asarray(z["fed_mask"], bool).copy()
            theta_rec = np.asarray(z["theta_rec"], float).copy()
            carried_floor = (float(z["floor"]), float(z["floor_err"]), str(z["floor_est"]))
            prev_cert_idx = np.asarray(z["cert_idx"], int)
            prev_q = np.asarray(z["q_of_psr"], float)
            p(f"  iter {it}: checkpoint exists -- resumed"); continue
        t_it = time.time()

        # (a) Fisher at the current joint state -- OWN-TERM-LIVE per member (fixed
        # 2026-07-24; the measured member is present, the OTHER carried are absent)
        carried = np.where(~led.fed)[0]
        sig_opt, sig_pes, F_ii = fisher_conditional(amo, jnp, theta_rec, carried, nd,
                                                    box_sigma, data, ones)
        # (b) the frontier: worst-axis ratio < GATE_RATIO -> fed (strict; SPARK-3's law)
        ratio = np.max(sig_opt / sig_pes, axis=1)
        want_fed = ratio < GATE_RATIO
        new = [int(k) for k in np.where(want_fed & ~led.fed)[0]]
        # FRONTIER-v2 (S4.20, wired 2026-07-26, Matt: "v1 must not be reachable going
        # forward"): feed requires BOTH the curvature condition above AND DATA SUPPORT
        # dlnL_feed > 0 -- the own-term-live present-vs-absent contrast (validated
        # S4.20.1: refused the wrong-sky copy at -6.5 nat, 0 false refusals). Refusals
        # are BANKED (feed_dlnl/feed_refused columns), never fed. Banked pre-2026-07-26
        # verdicts stand under the v1 scope line (S4.19/S4.20, no retro-apply); their
        # checkpoints resume PAST the feed decisions, so this path never re-decides them.
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
                p(f"  iter {it}: frontier-v2 refused {refused} at the data-support "
                  f"term (dlnL {[f'{feed_dlnl[k]:+.2f}' for k in refused]})")
        for k in new:
            # NULL SEMANTICS (amended 2026-07-25, Matt's S4.15 decision 1 -- the
            # pre-registration collision is recorded verbatim in the capstone as a
            # worked example): the census is SHARED, so a live frontier resolving REAL
            # members in the control arm is correct behaviour -- banked as the
            # control's own first-bite count (a baseline for the signal arm's).
            # STOP only on the manufacturing signal proper: a promote of the SCRAMBLED
            # member (slot 0 -- its template sky is wrong; resolving it from a wrong-sky
            # template is loop-grown structure). Wrong-cert and drain-fall STOPs below.
            if scrambled and k == 0:
                # S4.23 hook (additive): the STOP fired before the bank stage, so the
                # STOP-state anatomy was unreadable; bank it FIRST, then raise as before.
                bank_npz(f"{stem}_STOPANAT_i{it}",
                         stop_kind="scrambled_promote", iter=it,
                         theta_rec=theta_rec, fed_mask=led.fed.copy(),
                         promoted_k=k, conc_ratio=ratio,
                         arm=arm, sky=sky, real=real, t_label=T, wscale=wscale)
                raise CampaignStop(f"NULL-ARM SCRAMBLED-MEMBER PROMOTE at iter {it}: "
                                   f"the control resolved the sky-scrambled igniter -- "
                                   f"manufacturing. STOP + anatomy (banked at {ck}).")
            led.promote(k, slots[k], iteration=it)      # AT DRAWN TRUTH -- the frozen census
        fed_idx = np.where(led.fed)[0]

        # (c) THE DRAIN: refit the background at the new feed state
        smask = led.fed.astype(float)
        prof = bf.profile(theta_rec, data, ones, smask, a_grid)
        # null STOP (c), amended semantics: the control's drain must not FALL without a
        # new feed to explain it (a fall AT constant fed set = the M-step manufacturing
        # background removal from a template containing a wrong-sky member).
        if scrambled and prev_abg is not None and len(fed_idx) == prev_nfed and \
                prof["ahat"] < prev_abg - 2.0*np.sqrt(prof["sig"]**2 + prev_sig**2):
            raise CampaignStop(f"NULL-ARM DRAIN FALL at iter {it} without a new feed: "
                               f"{prof['ahat']:.3f} < {prev_abg:.3f} - 2sig -- "
                               f"manufacturing. STOP + anatomy.")
        prev_abg, prev_sig, prev_nfed = float(prof["ahat"]), float(prof["sig"]), len(fed_idx)

        # (d)+(e) E-step + M-step on the fed set (skipped while nothing is fed: the joint
        # template is empty and the fringe machinery has no live slot -- banked as such)
        n_eval = 0
        if len(fed_idx):
            # M-step surface at the SOFT-JOINT state: carried sources are ABSENT from the
            # template (they live in the fitted GP -- tiling theta_rec directly would
            # double-count them against the background). Only fed-slot values are written
            # back; carried slots in theta_rec keep their drawn-truth entries.
            # NB: recomputed AFTER the frontier promotes -- the (a)-stage `carried` is
            # stale here and would absent a just-promoted member from its own M-step.
            carried_now = np.where(~led.fed)[0]
            th_eval = theta_with_absent(theta_rec, nd, carried_now)
            def marg_fn(srcs):
                ths = np.tile(th_eval, (len(srcs), 1))
                ths[:, nd:] = srcs
                return np.asarray(amo["logL_batch_theta"](jnp.asarray(ths), data, ones))
            src_hat, widths, n_eval = mstep_quadratic(marg_fn, th_eval[nd:].copy(),
                                                      fed_idx, scale)
            theta_rec = theta_rec.copy()
            for k in fed_idx:
                theta_rec[nd + NP_SRC*k : nd + NP_SRC*(k+1)] = \
                    src_hat[NP_SRC*k : NP_SRC*(k+1)]

        # (f) scoreboard floors (policy) -- offender stat over null draws at THIS state
        zf = 0.0
        if full_floors or it in CARRY_REFIT_ITERS:
            off = _null_offenders(sb, ndraw, jnp, theta_rec, led, sky, it,
                                  prev_cert_idx, prev_q,
                                  n_null=N_NULL_FULL if full_floors else 32,
                                  wscale=wscale)
            zf = float((off == 0.0).mean())
            if float(np.ptp(off)) == 0.0:
                # Degenerate offender sample (typically ALL-ZERO at a faint joint state:
                # no null draw certifies anything). The Gumbel MLE is undefined here
                # (scipy's bracket search overflows on a constant sample), and ANCHOR
                # S4's high-zero-fraction branch is the empirical quantile anyway --
                # adopt it directly, labelled so the bank shows the estimator switch.
                fl, err, est = FL["emp_quantile"](off), 0.0, "emp_q95_degenerate"
            else:
                gu, mu, beta, gsd, nn = FL["gumbel_floor"](off)
                fl, err, est = FL["adopt"](gu, gsd, off, zf)
            if (carried_floor is not None and not full_floors
                    and abs(fl - carried_floor[0]) > 2.0 * max(err, carried_floor[1])):
                p(f"  DRIFT GATE FIRED at iter {it}: refit floor {fl:.2f} vs carried "
                  f"{carried_floor[0]:.2f} beyond band -> escalate this sky to full refits")
                full_floors = True
            carried_floor = (fl, err, est)
        fl, err, est = carried_floor

        # per-pulsar certification readout at the current joint state (G-d2 wiring:
        # PM rows cohere the PREVIOUS iteration's certified set at their q)
        dlnl, lnK, q_of, on_true = sb.columns(theta_rec, led, data, ones,
                                              prev_cert_idx, prev_q)
        cert = (dlnl > np.maximum(lnK + TRIALS_NAT, fl)) & (q_of > QBAR)
        wrong = cert & ~on_true
        if scrambled and wrong.any():
            # S4.23 hook (additive): bank the STOP-state columns first (the manufactured
            # certs' anatomy is the S4.23 pass's scoring material), then raise as before.
            bank_npz(f"{stem}_STOPANAT_i{it}",
                     stop_kind="wrong_cert", iter=it,
                     theta_rec=theta_rec, fed_mask=led.fed.copy(),
                     cert_idx=np.where(cert)[0], wrong_idx=np.where(wrong)[0],
                     dlnL_det=dlnl, lnK=lnK, q_of_psr=q_of, on_true=on_true,
                     floor=fl, floor_err=err, floor_est=est,
                     arm=arm, sky=sky, real=real, t_label=T, wscale=wscale)
            raise CampaignStop(f"NULL-ARM WRONG CERT at iter {it}: loop-grown "
                              f"manufacturing. STOP + anatomy.")
        cert_idx = np.argsort(np.where(cert, dlnl, -np.inf))[::-1][:int(cert.sum())]

        # (g) BANK
        areas = np.array([localisation_area_deg2(F_ii[k, [0, 1]]) for k in range(n_slot)])
        bank_npz(f"{stem}_i{it}",
                 n_cert=int(cert.sum()), n_resolved=led.n_resolved(),
                 a_bg=prof["ahat"], a_bg_sig=prof["sig"], a_bg_grid=prof["grid"],
                 a_bg_lnl=prof["lnl"],
                 conc_ratio=ratio, sig_opt=sig_opt, sig_pes=sig_pes,
                 feed_dlnl=feed_dlnl, feed_refused=np.array(refused, int),
                 fed_mask=led.fed, promote_events=led.event_array(),
                 areas_deg2=areas, epoch_cross=(areas < EPOCH_AREA_DEG2),
                 chan_budget=chan, n_clip=n_clip, band_log10f=np.array(BAND_CAMPAIGN),
                 a_eff_drawn=a_eff_drawn,
                 **({} if cond is None else {f"cond_{k}": v for k, v in cond.items()}),
                 t_label=T, wscale=wscale,
                 dlnL_det=dlnl, lnK=lnK, qmax=q_of, on_true=on_true,
                 wrong_cert=int(wrong.sum()), floor=fl, floor_err=err, floor_est=est,
                 zero_fraction=zf, full_floors=full_floors,
                 theta_rec=theta_rec, cert_idx=cert_idx, q_of_psr=q_of,
                 arm=arm, sky=sky, scrambled=scrambled, real=real, iter=it,
                 n_src=n_src, n_slot=n_slot, seed_census=SEED_POP_BASE + sky,
                 seed_noise=noise_seed, wall_s=time.time() - t_it,
                 orientation="row/col conventions: slot k == census rank k (brightest "
                             "first; harmonic slots appended after the census); psr "
                             "columns in P['names'] order; raw columns always banked")
        p(f"  iter {it}: N_res {led.n_resolved()} N_cert {int(cert.sum())} "
          f"A_bg {prof['ahat']:.3f}+-{prof['sig']:.3f} floor {fl:.2f}({est}) "
          f"wrong {int(wrong.sum())} [{time.time()-t_it:.0f}s]")
        prev_cert_idx, prev_q = cert_idx, q_of
    return 0


def _null_offenders(sb, ndraw, jnp, theta_rec, led, sky, it, prev_cert_idx, prev_q,
                    n_null, wscale=1.0):
    """Offender statistic over fresh no-CW null draws at THIS iteration's joint state
    (white+rn only -- consistent with the campaign's data model). The drawer is the
    CELL's drawer, built once (a NoiseDrawer construction refactors the GWB block --
    minutes at T=30 -- and must not be paid per refit).

    ================= DEFECT AND CORRECTION (E0), 2026-07-29, SIEVE-A =================
    WHAT WAS WRONG.  Until this date the statistic computed here was

        o[a] = q[a] > QBAR ? max(dlnL[a] - max(lnK[a], 0), 0) : 0        (WRONG)

    while the CANONICAL offender -- `recut_surface.offender:75`, and with it
    surface.py:334, chorus.py:605, ignite2.py:166, spark3.py:1063, anchor.py:322,
    kindle_d7_fall.py:67 and bank_surface_offenders.py:51, i.e. EVERY other campaign in
    this programme -- is

        off = max{ dlnL[a] : dlnL[a] > lnK[a] and q[a] > QBAR }, else 0   (CANONICAL)

    The floor returned here is compared against dlnL in criterion-v2.2
    (`run_cell`: cert = dlnL > max(lnK + TRIALS_NAT, floor)). So GLACIER measured its
    bar on the (dlnL - lnK) scale and applied it on the dlnL scale -- a scale mismatch
    in the PERMISSIVE direction. That this was unintended rather than a declared variant
    is visible two lines up: `_stack()` imports the canonical `offender` into FL and
    never calls it.

    HOW BIG.  Delta = floor_canonical - floor_glacier, measured on the banked null draws
    of 38 independent cells across three campaigns (GENERALISE arm A-SKY 32, CHORUS 4,
    IGNITE-2 2; SIEVE V1): median +6.97 nat, mean +7.42, range [+1.96, +15.37], and
    POSITIVE IN ALL 38. The sign is not empirical -- it is forced: lnK = log(max(K,1))
    >= 0, so the pulsar attaining the glacier max also satisfies dlnL > lnK and hence
    lies in the canonical mask, giving off_glacier <= off_canonical on every draw. For
    the same reason both statistics vanish on exactly the same draws, so `zero_fraction`
    and the degenerate-estimator branch below are UNAFFECTED by this correction.

    SCOPE.  Forward-correcting only. It changes the floor, and therefore only the
    criterion-v2.2 certification bar. It does NOT touch the D2 rigidity gate (R1
    2F_coh >= 15.132, R2 Delta-2F > 0), which never reads the floor -- so D2's
    manufactured-set kills stand unchanged. Banked verdicts are re-cut separately and
    reported as a before/after table (hpc_harbor/sieve/sv_v3_recut.py); they are NOT
    silently rewritten in place.
    ==================================================================================
    """
    off = []
    ones = jnp.ones(sb.npsr)
    n_nonfinite = 0
    for i in range(n_null):
        nz = ndraw.draw(9_000 + 100*sky + i, components=("white", "rn"),
                        white_scale=wscale)
        dnull = tuple(jnp.asarray(np.asarray(n)) for n in nz)
        dlnl, lnK, q_of, _ = sb.columns(theta_rec, led, dnull, ones,
                                        prev_cert_idx, prev_q)
        # Canonical statistic, with spark3.offender:1063's non-finite exclusion: a
        # pulsar with a single counted fringe (K = 1) has an infinite peak gap, and one
        # inf silently poisons the extreme-value fit. spark3 documents that guard as
        # part of the canonical form, so it rides along rather than being a variant.
        fin = np.isfinite(dlnl)
        n_nonfinite += int((~fin).sum())
        m = (dlnl > lnK) & (q_of > QBAR) & fin
        off.append(float(dlnl[m].max()) if m.any() else 0.0)
    if n_nonfinite:
        print(f"  [floor] excluded {n_nonfinite} non-finite dlnL rows across "
              f"{n_null} null draws", flush=True)
    return np.array(off)


class CertScoreboard:
    """G-d2 -- the per-target E-step scoreboard on the GLACIER venue.

    Reuses spark3's AUDITED wiring rather than reimplementing it: `estep_per_target`
    (SPARK-2 SS1's structural fix -- target p's own pulsar term live under ITS OWN mask, so
    every fringe row is live even at rung 0), `rung_masks` (certified coherent at their q,
    uncertified DECOHERED -- the hard-lock failure IGNITE-2 retired), and `score_from_LNL`
    (ignite.py:349-358 banked certification semantics, verbatim). Constructed ONCE per cell
    from GLACIER's own build.

    CARRIED sources are stated on this path via H_ABSENT slots in theta (exact to fp on a
    log-strain waveform; gated == the SMASK statement by driver gate G-d3) -- NOT via
    sp.set_smask, which on the pre-FORGE-G2 core invalidates all 116 compiled per-pulsar
    evaluators per flip (the recompile tax the fan is HELD on). The fringe machinery
    therefore never recompiles on a frontier change even on the baked core.

    THE FRINGE GRID: dL[a] = min mode spacing over ALL census+harmonic slots (the data's
    own comb -- more/finer sources sharpen the fringe pattern and grow K_counted, which is
    the correct trials pressure); EV/FringeTables/K_counted per the banked B_CERT=512
    convention (trackB_b1_core.build_EV / FringeTables, prior column = the declared tier)."""

    def __init__(self, G, amo, jnp, C, prior_key=TIER):
        import stagec_anchor_a2 as A2
        from cw_helpers import compute_mode_spacing
        import spark3 as SP3
        self.jnp, self.C, self.SP3 = jnp, C, SP3
        ent = G["ent_psrs"]
        names = [pe.name for pe in ent]
        npsr = len(ent)
        nd = amo["n_dist"]
        # dist keys are theta_keys[:nd] in disco==ent order; verify, never assume
        for i, pe in enumerate(ent):
            if amo["theta_keys"][i] != f"{pe.name}_cw_p_dist":
                raise CampaignStop(f"theta_keys[{i}] = {amo['theta_keys'][i]} != "
                                   f"{pe.name}_cw_p_dist -- distance-index map broken.")
        self.AI = np.arange(nd)
        self.L0 = np.array([pe.pdist[0] for pe in ent])
        self.priors = A2.load_real_priors(names)
        slots = G["slots"]
        self.dL = np.array([min(compute_mode_spacing(s[I_COSGT], s[I_GWPHI], s[I_FGW],
                                                     G["disco_psrs"][a].pos)
                                for s in slots) for a in range(npsr)])
        Pd = dict(npsr=npsr, L0=self.L0, priors=self.priors, ent_psrs=ent)
        self.P = Pd
        self.EV = C.build_EV(Pd, self.dL)
        self.FT = C.FringeTables(Pd, self.EV, self.dL, prior_key=prior_key)
        theta_truth = np.asarray(amo["theta_truth"], float)
        data_zero = amo["inject_delay"](jnp.asarray(theta_truth), jnp.ones(npsr))
        self.sp = C.B1Split(amo, theta_truth, data_zero, np.ones(npsr))
        self.npsr, self.nd = npsr, nd
        self.L_true = None                     # set per cell by draw_truth()

    def draw_truth(self, IG, dist_seed, tier=TIER):
        """Arm-B distance truth at the declared tier (ignite.draw_true_distances_tier --
        prior-weighted over SAMPLED fringes, so the truth is inside the search window;
        it returns (L_true, n_true), the within-fringe offset is folded into L_true)."""
        L_true, n_true = IG.draw_true_distances_tier(self.P, self.dL, self.EV,
                                                     seed=dist_seed, tier=tier)
        self.L_true = L_true
        return L_true, n_true

    def columns(self, theta_rec, led, data, pmask, prev_cert_idx, prev_q):
        """Per-pulsar (dlnL_det, lnK, q_max, on_true) at the current joint state.
        PM rows cohere the PREVIOUS iteration's certified set at their q (design line (d));
        carried census members are ABSENT from the scoring template (H_ABSENT statement)."""
        if self.L_true is None:
            raise CampaignStop("CertScoreboard.columns before draw_truth -- on_true would "
                               "be cut against an undrawn distance truth.")
        carried = np.where(~led.fed)[0]
        theta_base = theta_with_absent(np.asarray(theta_rec, float), self.nd, carried)
        theta_base[self.AI] = self.L0          # fringe offsets measured from L0 (ignite conv.)
        PM, take, pm_rung = self.SP3.rung_masks(self.npsr, np.asarray(prev_cert_idx, int),
                                                np.asarray(prev_q, float),
                                                k=len(prev_cert_idx))
        LNL = self.SP3.estep_per_target(self.sp, theta_base, self.EV, self.AI, data, PM,
                                        self.jnp)
        self._last_LNL = LNL                   # anatomy + the G-d2 liveness gate read this
        s = self.SP3.score_from_LNL(self.P, self.FT, LNL, self.C, self.L_true, self.dL)
        return (np.nan_to_num(np.asarray(s["dlnL"], float), posinf=1e30),
                np.asarray(s["lnK"], float), np.asarray(s["qmax"], float),
                np.asarray(s["on_true"], bool))


# ============================================================
# DRIVER GATES (mode gate) -- cheap, zero-noise, no campaign spend
# ============================================================
def mode_gate(verbose=True):
    """G-d1: one full soft iteration end-to-end at T=15, n=32, zero noise, no floors --
    plumbing only (census -> embed -> build -> Fisher -> frontier -> promote@truth ->
    drain refit -> M-step -> bank -> checkpoint resume identity).
    G-d2: the per-target E-step scoreboard (CertScoreboard) -- rung-mask structure,
    fringe-row LIVENESS with a fed member (the SPARK-2 0/116 failure class), rung-0 and
    rung-k paths finite, per-estep wall time printed for the budget.
    G-d3: theta_with_absent == SMASK-zero on the raw logL (the two statements of
    'carried' agree bit-comparably)."""
    check_affinity()
    print("=== GLACIER Stage-1 DRIVER GATES ===", flush=True)
    jax, jnp, C, B1, TE, IG, F, FL = _stack()
    pop = draw_population(SEED_POP_BASE, n_src=32, band_log10f=BAND_CAMPAIGN)
    slots, n_harm, active, chan, n_clip = embed_igniter(pop, 0.7, VENUE_SPAN_S[15])
    ok = True
    b = (len(slots) == 32 + 31) and (chan > 32)
    print(f"  embed: 32 -> {len(slots)} slots (+{n_harm}), chan {chan}, merger-guard "
          f"n_clip {n_clip} -> {'PASS' if b else 'FAIL'}")
    ok &= b
    # the guard itself: every slot's Earth term evaluable over (1-TERM_FLOOR)*span
    import spark3 as SP3g
    tcs = SP3g.t_coal(slots[:, I_MC], slots[:, I_FGW])
    bge = bool(np.all(tcs > (1.0 - TERM_FLOOR) * VENUE_SPAN_S[15] * (1 - 1e-9)))
    print(f"  merger guard: min t_coal/span = "
          f"{float(tcs.min())/VENUE_SPAN_S[15]:.3f} (> {1-TERM_FLOOR}) -> "
          f"{'PASS' if bge else 'FAIL'}")
    ok &= bge
    slots0, n0, _, chan0, nc0 = embed_igniter(pop, 0.0, VENUE_SPAN_S[15])
    b0 = (n0 == 0) and np.array_equal(slots0, pop["src"]) and chan0 == 32 and nc0 == 0
    print(f"  e=0 structural collapse (census unchanged, no appended slots): "
          f"{'PASS' if b0 else 'FAIL'}")
    ok &= b0
    pop_s = dict(pop); pop_s["src"] = slots; pop_s["n_src"] = len(slots)
    G = build_glacier_problem(15, pop_s, verbose=verbose)
    if abs(G["span_s"] / VENUE_SPAN_S[15] - 1.0) > 0.02:
        raise CampaignStop(f"built span {G['span_s']:.4e}s differs >2% from declared "
                           f"VENUE_SPAN_S[15] -- merger guard cut against a wrong span.")
    amo = G["amo"]
    nd = amo["n_dist"]
    theta = np.asarray(amo["theta_truth"], float)
    ones = jnp.ones(amo["npsr"])
    data = amo["inject_delay"](jnp.asarray(theta), ones)
    # G-d3: carried-absent (h=-18) vs smask-zero, on the background-fit residual terms
    bfit = BackgroundFit(amo, band_log10f=BAND_CAMPAIGN)
    smask0 = np.zeros(len(slots))
    p0_a, Ft_a = bfit._data_terms(theta, data, ones, smask0)
    th_abs = theta_with_absent(theta, nd, np.arange(len(slots)))
    p0_b, Ft_b = bfit._data_terms(th_abs, data, ones, np.ones(len(slots)))
    d_p0 = abs(p0_a - p0_b); d_ft = float(np.max(np.abs(Ft_a - Ft_b)))
    b3 = d_p0 < 1e-6 and d_ft < 1e-6
    print(f"  G-d3 carried-absent == smask-zero: |d p0| {d_p0:.3e}, max|d FtNmy| "
          f"{d_ft:.3e} (<1e-6) -> {'PASS' if b3 else 'FAIL'}")
    ok &= b3
    # drain profile at nothing-fed: gated on PLUMBING (finite everywhere); the edge/peak
    # QUESTION is g2a-ii's, under forensic (glacier_g2_forensic.py) -- at the realistic
    # point the T=15/n=32 profile legitimately pegs (the population sits below the in-band
    # noise), so an edge here is REPORTED, not failed.
    prof = bfit.profile(theta, data, ones, smask0, np.linspace(-15.6, -13.6, 21))
    b4 = bool(np.all(np.isfinite(prof["lnl"])))
    print(f"  drain profile: finite {b4}; ahat {prof['ahat']:.3f}+-{prof['sig']:.3f} "
          f"edge {prof['edge_hit']} (edge is REPORT-ONLY at this venue; see the g2a-ii "
          f"forensic) -> {'PASS' if b4 else 'FAIL'}")
    ok &= b4

    # ---- G-d2: the per-target E-step scoreboard ----
    G["slots"] = slots
    sb = CertScoreboard(G, amo, jnp, C, prior_key=TIER)
    L_true, n_true = sb.draw_truth(IG, SEED_NOISE_BASE + 10_000_000, tier=TIER)
    theta_t = theta.copy(); theta_t[sb.AI] = L_true
    data_t = amo["inject_delay"](jnp.asarray(theta_t), ones)
    # rung-mask structure (synthetic certified pair)
    import spark3 as SP3
    fake_q = np.zeros(sb.npsr); fake_q[[3, 7]] = (0.95, 0.99)
    PM, take, pm_rung = SP3.rung_masks(sb.npsr, np.array([3, 7]), fake_q, k=2)
    b5 = (np.allclose(np.diag(PM), 1.0) and np.allclose(PM[0, [3, 7]], (0.95, 0.99))
          and np.allclose(pm_rung[[3, 7]], (0.95, 0.99))
          and float(np.abs(pm_rung).sum()) == float(np.abs(pm_rung[[3, 7]]).sum()))
    print(f"  G-d2a rung-mask structure (own-term live, certified at q, rest decohered): "
          f"{'PASS' if b5 else 'FAIL'}")
    ok &= b5
    # rung-0 with a FED member: every fringe row must be LIVE (the SPARK-2 0/116 class)
    led = PromoteLedger(slots)
    led.promote(0, slots[0], iteration=0)
    t0 = time.time()
    dlnl, lnK, q_of, on_true = sb.columns(theta_t, led, data_t, ones,
                                          np.zeros(0, int), np.zeros(sb.npsr))
    t_estep = time.time() - t0
    LNL0 = sb._last_LNL.copy()
    row_rng = LNL0.max(axis=1) - LNL0.min(axis=1)
    live = float((row_rng > 0).mean())
    b6 = (np.isfinite(dlnl).all() and np.isfinite(lnK).all()
          and np.all((q_of >= 0) & (q_of <= 1.0 + 1e-12)) and live == 1.0)
    print(f"  G-d2b rung-0 columns, brightest fed: live rows {live*100:.1f}% "
          f"(need 100), dlnL finite {bool(np.isfinite(dlnl).all())}, q in [0,1] "
          f"{bool(np.all((q_of >= 0) & (q_of <= 1.0 + 1e-12)))}; estep wall "
          f"{t_estep:.1f}s -> {'PASS' if b6 else 'FAIL'}")
    ok &= b6
    # rung-k path (cohere a synthetic certified pair): runs, finite, and the coherence
    # actually changes the surface (structural, not a physics cut at this faint venue)
    dlnl2, lnK2, q2, _ = sb.columns(theta_t, led, data_t, ones, np.array([3, 7]), fake_q)
    moved = float(np.max(np.abs(sb._last_LNL - LNL0)))
    b7 = np.isfinite(dlnl2).all() and np.isfinite(q2).all()
    print(f"  G-d2c rung-2 coherent path finite: {'PASS' if b7 else 'FAIL'} "
          f"(max|LNL shift| from coherence {moved:.3e}, report-only at this faint venue)")
    ok &= b7

    # ---- G-d4: FRONTIER POSITIVE CONTROL (added 2026-07-24 after the dead-frontier
    # defect: the pre-fix fisher_conditional absented the measured member, so ratio was
    # 1.000 at ANY brightness and nothing could ever ignite -- caught by the r13p25
    # ladder cells, NOT by these gates, because no gate ever put a bright source across
    # the frontier. This one does, permanently.) ----
    box_sigma_g = np.array([1.0, np.pi, 1.0, 0.5, 0.25, 0.25, np.pi, 0.5*np.pi]) / np.sqrt(3.0)
    th_b = np.asarray(amo["theta_truth"], float).copy()
    th_b[nd + I_H] = -13.25                       # slot 0 at the onset class
    # ... AND at a fringe-CAPABLE frequency (25 nHz). The census's own brightest sits at
    # ~2 nHz (the f^-11/3 draw piles at the extended band's floor) where a source
    # completes ~1 cycle over the span -- pulsar-term fringe evidence is VENUE-starved
    # there at any brightness (measured: 0.06 nat at h=-13.25, first G-d4 run). This
    # control tests the MACHINERY, so the control source must live where fringes exist;
    # the venue-physics number at the census fgw is the ladder's own cert-curve result.
    th_b[nd + I_FGW] = -7.6
    data_b = amo["inject_delay"](jnp.asarray(th_b), ones)
    so, spn, _ = fisher_conditional(amo, jnp, th_b, np.arange(len(slots)), nd,
                                    box_sigma_g, data_b, ones)
    rat = np.max(so / spn, axis=1)
    b8a = rat[0] < GATE_RATIO
    b8b = float(np.median(rat[1:])) > 0.9
    print(f"  G-d4a frontier feeds a -13.25 member: ratio[0] = {rat[0]:.4f} "
          f"(< {GATE_RATIO}) -> {'PASS' if b8a else 'FAIL'}")
    print(f"  G-d4b faint census stays box-limited: median ratio = "
          f"{float(np.median(rat[1:])):.3f} (> 0.9) -> {'PASS' if b8b else 'FAIL'}")
    slots4 = slots.copy()
    slots4[0, I_H] = -13.25; slots4[0, I_FGW] = -7.6  # control census: slot 0 bright + fast
    led4 = PromoteLedger(slots4)
    led4.promote(0, slots4[0], iteration=0)           # at (control) drawn truth
    dl4, _, _, _ = sb.columns(th_b, led4, data_b, ones, np.zeros(0, int),
                              np.zeros(amo["npsr"]))
    # BAR CALIBRATION (2 measured rounds): the pre-fix structural zero was EXACTLY
    # 0.000000 (flat LNL over the sweep -- the regression this leg guards). Live values
    # at this cheap T=15/n=32 smoke venue: 0.06 nat at the census fgw (~2 nHz, fringe-
    # starved) and 0.20 nat at 25 nHz -- small because the NG15-level GWB GP legitimately
    # absorbs most of one source's in-band evidence HERE. Whether a fed source CERTIFIES
    # at the campaign venue (T=30/40) is the LADDER's own cert-curve measurement, not
    # this gate's. Bar = an order of magnitude above fp-noise on a flat sweep, i.e.
    # "the scoreboard RESPONDS to a fed source at all".
    b8c = float(np.max(dl4)) > 0.05
    print(f"  G-d4c fed bright member lights the scoreboard: max dlnL = "
          f"{float(np.max(dl4)):.2f} (> 0.05; machinery-live bar, see comment) -> "
          f"{'PASS' if b8c else 'FAIL'}")
    ok &= b8a and b8b and b8c

    print(f"\n=== DRIVER GATES: "
          f"{'PASS (fan HELD on FORGE-G2 runtime-SMASK + resume drill)' if ok else 'FAIL'} "
          f"===", flush=True)
    return 0 if ok else 1


def _redirect_out(sub):
    """Point BOTH modules' OUT at a drill subdirectory (bank_npz reads glacier_pop.OUT;
    run_cell's checkpoint path reads this module's OUT). Drill-only, declared."""
    import glacier_pop as GPM
    d = f"{GPM.OUT}/{sub}"
    os.makedirs(d, exist_ok=True)
    GPM.OUT = d
    globals()["OUT"] = d
    return d


def mode_drillcmp(sub_a="drill_resume", sub_b="drill_straight", n_iter=3):
    """THE RESUME DRILL VERDICT: the interrupted-and-resumed leg's final-iteration bank
    must equal the uninterrupted leg's, column for column (same lane, same seeds -- any
    drift means resume is not a continuation and the fan must not trust checkpoints)."""
    import glacier_pop as GPM
    it = n_iter - 1
    ok = True
    for kind in ("cell",):
        sa = glob.glob(f"{GPM.OUT}/{sub_a}/gl1_{kind}_*_i{it}__*.npz")
        sb_ = glob.glob(f"{GPM.OUT}/{sub_b}/gl1_{kind}_*_i{it}__*.npz")
        if len(sa) != 1 or len(sb_) != 1:
            print(f"  drillcmp: expected exactly one i{it} bank per leg, found "
                  f"{len(sa)}/{len(sb_)} -- run the drill legs first"); return 1
        za, zb = np.load(sa[0], allow_pickle=True), np.load(sb_[0], allow_pickle=True)
        for col in ("theta_rec", "fed_mask", "a_bg", "n_cert", "dlnL_det", "lnK", "qmax",
                    "floor", "zero_fraction", "cert_idx"):
            da = np.asarray(za[col]); db = np.asarray(zb[col])
            same = da.shape == db.shape and bool(np.array_equal(da, db))
            print(f"  drillcmp {col}: {'IDENTICAL' if same else 'DIFFERS'}"
                  + ("" if same or da.dtype.kind not in "fc" else
                     f" (max|d| {np.max(np.abs(da.astype(float)-db.astype(float))):.3e})"))
            ok &= same
    print(f"\n=== RESUME DRILL: {'PASS -- checkpoints are continuations' if ok else 'FAIL'}"
          f" ===", flush=True)
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "cell", "null", "cell2", "null2",
                                     "drillcell", "drillcmp"])
    ap.add_argument("--rung", choices=list(LADDER_RUNGS), default="r13p9",
                    help="STAGE-2a ladder rung (cell2/null2 only)")
    ap.add_argument("--wscale", type=float, default=1.0, choices=[1.0, 0.5, 0.25],
                    help="STAGE-2a ARRAY ladder: white-noise rms scale (cell2/null2)")
    ap.add_argument("--t", type=int, default=T_LABEL, choices=[30, 40],
                    help="STAGE-2a ARRAY ladder: venue T (cell2/null2)")
    ap.add_argument("--arm", choices=list(ARMS), default="e07")
    ap.add_argument("--sky", type=int, default=0)
    ap.add_argument("--real", type=int, default=0)
    ap.add_argument("--n", type=int, default=N_POP_DEFAULT)
    ap.add_argument("--iters", type=int, default=N_ITER)
    ap.add_argument("--sub", default="drill_resume",
                    help="drill output subdirectory (drillcell/drillcmp)")
    ap.add_argument("--sub-b", default="drill_straight")
    a = ap.parse_args()
    if a.mode == "gate":
        return mode_gate()
    if a.mode == "drillcell":
        # the RESUME DRILL leg: the full run_cell path at a small census (default n=32),
        # carried floors (32 nulls at refit iters), its own OUT subdir. Noisy T=30 --
        # run ONLY on the post-FORGE-G2 tree (the drill baseline must be the final tree).
        _redirect_out(a.sub)
        return run_cell(a.arm, a.sky, n_src=a.n, n_iter=a.iters, full_floors=False)
    if a.mode == "drillcmp":
        return mode_drillcmp(a.sub, a.sub_b, n_iter=a.iters)
    if a.mode in ("cell", "null"):
        # THE HOLD, AS CODE (amended 2026-07-23, all four items authorized): cell/null run
        # only when the holds-cleared marker exists. The marker is written BY HAND after the
        # three gate jobs are green on this tree -- (1) g2a-ii re-gate at BAND_CAMPAIGN,
        # (2) FORGE-G2 smask battery + Bg5 on this lane, (3) resume drill + driver re-gate --
        # and records their job ids. Deleting the marker re-arms the hold.
        import glacier_pop as GPM
        marker = f"{GPM.OUT}/HOLDS_CLEARED"
        if not os.path.exists(marker):
            print(f"REFUSED: no {marker}. The Stage-1 fan is HELD until the three gate "
                  "jobs are green on this tree: (1) g2a-ii re-gate at BAND_CAMPAIGN "
                  "(remedy A), (2) FORGE-G2 runtime-SMASK battery (SG1-SG6) + Bg5 on this "
                  "lane, (3) resume drill (drillcell x2 + drillcmp) + driver-gate re-run. "
                  "The marker records the green job ids. This refusal is the hold, as code.",
                  flush=True)
            return 3
        with open(marker) as fh:
            print(f"[GLACIER-1] holds cleared: {fh.read().strip()}", flush=True)
        return run_cell(a.arm, a.sky, scrambled=(a.mode == "null"), real=a.real)
    if a.mode in ("cell2", "null2"):
        # STAGE-2a conditioning ladder (Option C, Matt 2026-07-24). Same evidence gate
        # as stage 1 (the marker), plus the ladder's own CPU gates must have passed
        # (glacier_pop.py gate2 -- gl2_ladder_gates bank exists).
        import glacier_pop as GPM
        marker = f"{GPM.OUT}/HOLDS_CLEARED"
        gates2 = glob.glob(f"{GPM.OUT}/gl2_ladder_gates__*.npz")
        if not (os.path.exists(marker) and gates2):
            print(f"REFUSED: stage-2a needs {marker} AND a gl2_ladder_gates bank "
                  "(glacier_pop.py gate2).", flush=True)
            return 3
        return run_cell(a.arm, a.sky, scrambled=(a.mode == "null2"), real=a.real,
                        rung=a.rung, t_label=a.t, wscale=a.wscale)
    return 2


if __name__ == "__main__":
    sys.exit(main())
