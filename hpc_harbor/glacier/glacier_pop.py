"""GLACIER g2 -- the DRAWN POPULATION + the DRAINABLE BACKGROUND (agent GLACIER, ACCRE).

THE QUESTION THIS FEEDS. Does the array EAT THE BACKGROUND? The soft-joint loop on a
population where the background IS unresolved sources. This module is the population build
and its gates -- no campaign, no verdict.

WHAT g2 BUILDS
--------------
(1) POPULATION SYNTHESIS (draw_population): N ~ 200-500 binaries.
      f      : dN/df ~ f^(-11/3) (GW-driven residence time), inverse-CDF inside the
               GENERATIVE band LOG10_FGW_RANGE = (-8.0, -7.5) -- the estimator's own
               population box (stagec_fisher), NOT TE.phi_bounds (SPARK-3 SS2.2).
      h      : LOG-NORMAL strain hierarchy (sigma_logh dex, declared), NORMALISED so the
               drawn set's band-integrated characteristic-strain power equals the NG15-class
               target: sum_i h_i^2 = A^2 * INT_band (f/f_yr)^(-4/3) df/f with
               log10 A = -14.6, gamma = 13/3 (the repo's NG15 convention -- stagec_fisher
               GWB_LOG10_A/GWB_GAMMA; project_progress "GWB at NG15 -14.6").
      e      : LOTTERY's two-population mix, PARAMETERISED (f_ecc, e_char); the REALISTIC
               POINT default is (1.0, 0.3) -- the mix representation of the realistic
               astrophysical case, e ~ lnN(mode 0.3) (LOTTERY_tier1 SS"realistic"). Banked
               per source; the Stage-1 IGNITER ARMS OVERRIDE the loud members' e explicitly
               (e=0.7 single / all-circular NONE), so at Stage 1 the drawn e is report-only.
      mc, sky, cos_inc, phase0, psi : the generative box (uniform), per stagec_fisher.
(2) THE DRAINABLE BACKGROUND (BackgroundFit): the existing "gwb" global GP (HD, powerlaw
    gamma 13/3) with log10_A a FITTED parameter -- profiled by rebuilding
    Pinv(A) + blockdiag(c_truth) and re-factoring per grid point (the amortised build's
    cf_cached bypassed via its own internals; f64 factorisation mandate held). Everything
    below the model-quality frontier is carried COLLECTIVELY by this GP; resolved members
    are discrete templates. THE DRAIN CURVE is A_background(iter).
(3) THE FROZEN-CENSUS LEDGER (PromoteLedger): the campaign's whole safety case
    (GLACIER_capstone SSL.2). The census is drawn once per sky and never grows; a frontier
    update changes ONLY a source's feed weight w_s (carried -> fed); every component ENTERS
    at its TRUE DRAWN position, bit-exact. Any promote at a fitted position raises
    CampaignStop (B0.2 applies -- STORY S8.5.5: a self-found counterpart is a wrong
    counterpart essentially always).

THE GATES (blocking; artifact-readback discipline -- every verdict cut from the banked npz)
------------------------------------------------------------------------------------------
 g2a-i   DRAW NORMALISATION (CPU, numpy): the drawn set's projection onto the gamma=13/3
         powerlaw over the GP bins reproduces A_TARGET. |dlog10 A_eff| < 0.01 (the draw is
         normalised on the total band power; the projection residual is the drawn
         realisation's shape scatter, banked as anatomy, not gated).
 g2a-ii  THE FIT GATE (venue GPU/CPU job): at iteration 0 with NOTHING resolved (no source
         fed, smask all-zero), the FITTED background amplitude on the zero-noise Asimov
         data (all N sources injected as deterministic delays; no stochastic GWB in the
         data -- the background IS the unresolved sum) reproduces the drawn population's
         summed power: |log10 Ahat - log10 A_eff| < 0.15 dex, PRE-STATED. The tolerance
         absorbs (a) single-draw discreteness across bins, (b) the HD-vs-point-source ORF
         mismatch, (c) powerlaw-shape mismatch; the per-bin anatomy is banked beside the
         verdict either way. A failure is a FINDING, not a tuning knob.
 g2b     CONSERVATION (CPU, numpy): resolved + residual power is conserved across the
         frontier BY CONSTRUCTION -- sweep the frontier through every census rank and
         check P_fed + P_carried == P_total to 1e-12 relative, at every position.
 g2c     FROZEN CENSUS: the PromoteLedger accepts a promote ONLY at the drawn truth
         (bit-exact mean identity gated) and raises CampaignStop on a fitted-position
         promote (gated via a deliberate violation).
 g2d     PROVENANCE: every bank goes through bank_npz() -- lane tag (host_GPU) in the
         FILENAME (spark3 _lane_tag; the venueself sort-order accident is the 4th member
         of this hazard family), host/cpus/affinity/L_gwb-fingerprint/T metadata INSIDE,
         separate GLACIER_results/, and a REFUSAL to overwrite a bank carrying a
         different lane tag. Gated by write/readback/collision.

AFFINITY REFUSAL (adopted from LOTTERY_recut_band SS"metadata"): the driver REFUSES to run
if the observed CPU affinity disagrees with the declared convention (8) -- cpus-per-task is
part of the provenance (the NoiseDrawer thread-count hazard; here no GWB noise is ever
drawn, but the posture is held campaign-wide).

SEEDS (house convention, disjoint from every banked range): population draws 9_0xx_xxx;
Stage-1 noise 9_5xx_xxx (dist = noise + 10_000_000).

Matt commits; never commit from here. Env: batch_gpu H200 lane, autotune off, x64.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse, glob, socket, sys, time
import numpy as np

HSYMT = "/home/mattm/projects/HSYMT"
sys.path.insert(0, f"{HSYMT}/CW_lnL_check")
sys.path.insert(0, f"{HSYMT}/CW_transition")
sys.path.insert(0, f"{HSYMT}/hpc_harbor/ignite")
OUT = os.environ.get("GLACIER_OUT", f"{HSYMT}/GLACIER_results")

I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)
NP_SRC = 8
YR = 365.25 * 86400.0
F_YR = 1.0 / YR

# ---- the declared population (every number a parameter, every default recorded) ----
N_POP_DEFAULT = 256            # N ~ 200-500 per the brief; declared, parameterised
SIGMA_LOGH = 0.6               # log-normal strain hierarchy width, dex. DECLARED, provisional.
A_TARGET_LOG10 = -14.6         # the repo's NG15 convention (stagec_fisher.GWB_LOG10_A)
AEFF_SANITY_DEX = 0.16         # MANDATORY generator-bug tripwire on |log10 A_eff - target|:
                               # the observed max over the banked 200-seed ensemble at n=256
                               # (reports/glacier_aeff_ensemble.npz; amendment 2026-07-23)
GAMMA_BG = 13.0 / 3.0          # residual powerlaw index; h_c ~ f^(-2/3)
LOG10_FGW_RANGE = (-8.0, -7.5)   # the GENERATIVE band (stagec_fisher; == estimator pop box)
LOG10_MC_RANGE = (8.5, 9.5)
REALISTIC_POINT = dict(f_ecc=1.0, e_char=0.3)   # mix repr. of LOTTERY's lnN(mode 0.3) case
SEED_POP_BASE = 9_000_000
DECLARED_CPUS = 8


class CampaignStop(RuntimeError):
    """A pre-registered STOP condition fired. Never caught inside the campaign."""


def check_affinity(declared=DECLARED_CPUS):
    got = len(os.sched_getaffinity(0))
    if got != declared:
        raise CampaignStop(
            f"CPU affinity {got} != declared convention {declared}. cpus-per-task is part "
            f"of the provenance (thread-count hazard); refusing to run off-convention.")
    return got


# ============================================================
# provenance -- g2d. Lane in the FILENAME, basis in the METADATA.
# ============================================================
def lane_tag():
    """spark3._lane_tag verbatim in spirit: host_GPU. An artifact whose name does not
    distinguish the conditions it was produced under is not an artifact, it is a race."""
    h = socket.gethostname().split(".")[0]
    try:
        import subprocess
        g = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                           capture_output=True, text=True, timeout=20).stdout.split("\n")[0]
        g = g.strip().replace(" ", "") or "noGPU"
    except Exception:
        g = "noGPU"
    return f"{h}_{g}"


def bank_npz(stem, allow_lane_mismatch=False, **cols):
    """Write {OUT}/{stem}__{lane}.npz with host/thread provenance INSIDE. Refuses to
    overwrite any co-resident bank of the same stem from a DIFFERENT lane unless the caller
    says so explicitly -- the venueself sort-order accident (capstone g1.1b) must not get a
    fifth member. Returns the path."""
    os.makedirs(OUT, exist_ok=True)
    lt = lane_tag()
    others = [f for f in glob.glob(f"{OUT}/{stem}__*.npz")
              if not f.endswith(f"__{lt}.npz")]
    if others and not allow_lane_mismatch:
        raise CampaignStop(
            f"bank stem '{stem}' already exists under a different lane tag: {others}. "
            f"Two provenances must not share one stem (capstone g1 mechanism). Pass "
            f"allow_lane_mismatch=True only for a deliberate cross-lane study.")
    path = f"{OUT}/{stem}__{lt}.npz"
    np.savez(path, **cols,
             _lane=lt, _host=socket.gethostname(),
             _cpus_affinity=len(os.sched_getaffinity(0)),
             _declared_cpus=DECLARED_CPUS,
             _slurm_job=os.environ.get("SLURM_JOB_ID", "none"),
             _time=time.strftime("%Y-%m-%dT%H:%M:%S"))
    return path


# ============================================================
# (1) POPULATION SYNTHESIS
# ============================================================
def _band_power_target(a_log10=A_TARGET_LOG10, band_log10f=LOG10_FGW_RANGE):
    """INT_band h_c^2(f) df/f for h_c = A (f/f_yr)^(-2/3): the total the drawn strains sum
    to. Closed form: A^2 * (3/4) * [ (f/f_yr)^(-4/3) ]_hi^lo."""
    A2 = 10.0 ** (2 * a_log10)
    flo, fhi = (10.0 ** np.array(band_log10f))
    x = lambda f: (f / F_YR) ** (-4.0 / 3.0)
    return A2 * 0.75 * (x(flo) - x(fhi))


def draw_population(seed, n_src=N_POP_DEFAULT, sigma_logh=SIGMA_LOGH,
                    f_ecc=REALISTIC_POINT["f_ecc"], e_char=REALISTIC_POINT["e_char"],
                    a_target_log10=A_TARGET_LOG10, band_log10f=LOG10_FGW_RANGE):
    """The census, drawn ONCE per sky seed and frozen. Returns a dict of per-source columns
    plus the normalisation anatomy. Pure numpy -- bit-reproducible at any thread count.
    band_log10f: the GENERATIVE band. Default = the CW estimator pop box (the incumbent
    convention); the g2a-ii forensic showed an NG15-normalised population confined there
    sits below the array's in-band noise at every campaign venue, so the band is now a
    DECLARED PARAMETER (banked with every draw) rather than an inheritance."""
    rng = np.random.default_rng(seed)
    flo, fhi = 10.0 ** np.array(band_log10f)
    # f: dN/df ~ f^(-11/3)  ->  CDF^-1(u) = [flo^(-8/3) - u*(flo^(-8/3)-fhi^(-8/3))]^(-3/8)
    u = rng.uniform(size=n_src)
    p = -8.0 / 3.0
    fgw = (flo ** p - u * (flo ** p - fhi ** p)) ** (1.0 / p)
    # h: log-normal hierarchy, then normalise total band power to the target
    logh_raw = rng.normal(0.0, sigma_logh, size=n_src)          # shape only; centre fixed below
    h = 10.0 ** logh_raw
    total = float(np.sum(h ** 2))
    target = _band_power_target(a_target_log10, band_log10f)
    h *= np.sqrt(target / total)                                 # sum_i h_i^2 == target, exact
    # e: the two-population mix (LOTTERY draw_mix, per source)
    is_ecc = rng.random(n_src) < f_ecc
    ecc = np.where(is_ecc, e_char, 0.0)
    # the rest of the 8-block from the generative box
    src = np.zeros((n_src, NP_SRC))
    src[:, I_COSGT] = rng.uniform(-1, 1, n_src)
    src[:, I_GWPHI] = rng.uniform(0, 2 * np.pi, n_src)
    src[:, I_COSINC] = rng.uniform(-1, 1, n_src)
    src[:, I_MC] = rng.uniform(*LOG10_MC_RANGE, n_src)
    src[:, I_FGW] = np.log10(fgw)
    src[:, I_H] = np.log10(h)
    src[:, I_PH0] = rng.uniform(0, 2 * np.pi, n_src)
    src[:, I_PSI] = rng.uniform(0, np.pi, n_src)
    order = np.argsort(src[:, I_H])[::-1]                        # census rank: brightest first
    src = src[order]; ecc = ecc[order]; is_ecc = is_ecc[order]
    return dict(src=src, ecc=ecc, is_ecc=is_ecc, n_src=n_src, seed=seed,
                sigma_logh=sigma_logh, f_ecc=f_ecc, e_char=e_char,
                a_target_log10=a_target_log10, band_power_target=target,
                band_log10f=tuple(band_log10f))


# ---- STAGE-2a: the conditioning ladder (authorized 2026-07-24, Option C) ----
SIGMA_NG15_DEX = 0.05          # DECLARED NG15-class posterior width on log10 A -- the
                               # tension unit quoted per rung (the -13.25 rung's ~12-sigma
                               # / 16x-power number goes ON the figure axis, not hidden)
SEED_POP2_BASE = SEED_POP_BASE + 500     # stage-2 sky seeds: a fresh declared block
COND_SCAN_STRIDE = 1_000_000             # tail-mode seed stride per sky (no seed reuse)

# The ladder. Rung 1 is the NG15-CONSISTENT ceiling rung, realized as tail selection at
# threshold -13.9 (P_cond = 1.4e-2 from the banked 100k scan, gen_stage2_pcond_scan.npz;
# the hard ceiling itself is sqrt(target) = -13.867 and has measure ~3e-4). Rungs 2-3 are
# DECLARED super-NG15 skies: the NG15-consistent draw plus ONE exceptional source (the
# brightest member SET to the rung; the other 255 keep their NG15-normalised strains).
# P_cond(onset class) = 0 EXACTLY under power conservation -- that zero is a RESULT:
# "within the NG15-consistent class, no sky contains an onset-class source; the first
# resolved CW will itself be evidence of excess power or exceptional structure; the
# ladder measures how far above the median background the sky must be before the array
# begins to eat it." (framing banked verbatim per the 2026-07-24 readback)
LADDER_RUNGS = {
    "r13p9":  dict(rung_log10h=-13.9,  mode="tail"),
    "r13p5":  dict(rung_log10h=-13.5,  mode="set"),
    "r13p25": dict(rung_log10h=-13.25, mode="set"),
}


def draw_population_conditional(sky, rung_key, n_src=N_POP_DEFAULT,
                                band_log10f=LOG10_FGW_RANGE, max_scan=200_000):
    """One conditional-sky draw for ladder rung `rung_key`, sky index `sky`.

    TAIL mode (NG15-consistent): rejection-scan the unconditional ensemble from this
    sky's own seed block until the brightest member clears the rung threshold; the draw
    is an ORDINARY NG15 sky -- just a lucky one -- and the number of trials is banked.
    SET mode (super-NG15, declared): unconditional draw at the sky seed, then the
    brightest member's strain is SET to the rung. Total band power now exceeds the NG15
    target by the exceptional source's excess; the equivalent background amplitude and
    its tension vs the NG15 posterior (SIGMA_NG15_DEX) are computed here and banked.
    Returns (pop, cond) -- pop as draw_population, cond = the conditioning record."""
    r = LADDER_RUNGS[rung_key]
    rung, mode = r["rung_log10h"], r["mode"]
    base = SEED_POP2_BASE + sky * COND_SCAN_STRIDE
    if mode == "tail":
        for j in range(max_scan):
            pop = draw_population(base + j, n_src=n_src, band_log10f=band_log10f)
            if pop["src"][0, I_H] >= rung:
                break
        else:
            raise CampaignStop(f"tail rung {rung_key}: no qualifying draw in {max_scan} "
                               f"seeds from {base} -- conditioning mis-specified.")
        h0 = float(pop["src"][0, I_H])
        cond = dict(rung_key=rung_key, rung_log10h=rung, cond_mode="tail",
                    h_brightest=h0, h_unconditioned=h0, n_scanned=j + 1,
                    excess_power_ratio=1.0, a_equiv_log10=pop["a_target_log10"],
                    tension_sigma=0.0)
    elif mode == "set":
        pop = draw_population(base, n_src=n_src, band_log10f=band_log10f)
        h0 = float(pop["src"][0, I_H])
        pop["src"][0, I_H] = rung                      # the one exceptional source
        target = pop["band_power_target"]
        excess = (10.0 ** (2 * rung) - 10.0 ** (2 * h0))
        ratio = (target + excess) / target
        a_eq = pop["a_target_log10"] + 0.5 * np.log10(ratio)
        cond = dict(rung_key=rung_key, rung_log10h=rung, cond_mode="set",
                    h_brightest=rung, h_unconditioned=h0, n_scanned=1,
                    excess_power_ratio=float(ratio), a_equiv_log10=float(a_eq),
                    tension_sigma=float((a_eq - pop["a_target_log10"]) / SIGMA_NG15_DEX))
    else:
        raise CampaignStop(f"unknown conditioning mode {mode}")
    return pop, cond


def gate_ladder(n_sky=4, n_src=N_POP_DEFAULT, band_log10f=LOG10_FGW_RANGE, verbose=True):
    """CPU gates for the Stage-2a ladder draws (pure numpy). Per rung x sky:
    (i) the brightest member sits at/above the rung -- exact for set, >= for tail;
    (ii) power conservation: tail = NG15 target EXACT; set = target + declared excess,
         identity < 1e-12;
    (iii) the tension record is finite and monotone up the ladder.
    Banks the whole ladder census (gl2_ladder bank)."""
    p = print if verbose else (lambda *a, **k: None)
    ok = True
    rows = []
    p("\n== STAGE-2a LADDER GATES (conditional draws, CPU) ==")
    for rk in LADDER_RUNGS:
        for sky in range(n_sky):
            pop, cond = draw_population_conditional(sky, rk, n_src=n_src,
                                                    band_log10f=band_log10f)
            pw = float(np.sum(source_band_power(pop["src"])))
            want = pop["band_power_target"] * cond["excess_power_ratio"]
            d_pw = abs(pw / want - 1.0)
            b_h = pop["src"][0, I_H] >= cond["rung_log10h"] - 1e-12
            b_pw = d_pw < 1e-12
            ok &= b_h and b_pw
            rows.append((rk, sky, cond["h_brightest"], cond["n_scanned"],
                         cond["excess_power_ratio"], cond["a_equiv_log10"],
                         cond["tension_sigma"], d_pw))
            p(f"  [{rk:7s} sky {sky}] brightest {pop['src'][0, I_H]:+.3f} "
              f"(scan {cond['n_scanned']:>5d}) power/declared-1 = {d_pw:.2e} "
              f"A_eq {cond['a_equiv_log10']:+.3f} tension {cond['tension_sigma']:+.1f}s "
              f"-> {'PASS' if (b_h and b_pw) else 'FAIL'}")
    cols = np.array([(r[2], r[3], r[4], r[5], r[6], r[7]) for r in rows])
    bank_npz("gl2_ladder_gates",
             rungs=np.array([r[0] for r in rows]), skies=np.array([r[1] for r in rows]),
             h_brightest=cols[:, 0], n_scanned=cols[:, 1],
             excess_power_ratio=cols[:, 2], a_equiv_log10=cols[:, 3],
             tension_sigma=cols[:, 4], d_power_identity=cols[:, 5],
             sigma_ng15_dex=SIGMA_NG15_DEX, band_log10f=np.array(band_log10f),
             framing="within the NG15-consistent class, no sky contains an onset-class "
                     "source; the first resolved CW will itself be evidence of excess "
                     "power or exceptional structure; the ladder measures how far above "
                     "the median background the sky must be before the array begins to "
                     "eat it. P_cond(onset)=0 exactly; feasible-rung P_cond=1.4e-2 "
                     "(gen_stage2_pcond_scan.npz). Patience numbers travel beside it: "
                     "T~75yr median sky, ~53yr best sky.")
    p(f"== LADDER GATES: {'ALL PASS' if ok else 'FAIL'} ==")
    return 0 if ok else 1


def source_band_power(src):
    """Per-source contribution to the band-integrated h_c^2 df/f budget: h_i^2 exactly
    (see _band_power_target -- both sides of the ledger live in the same units, so the
    conservation identity g2b is arithmetic, not modelling)."""
    return 10.0 ** (2 * src[:, I_H])


def a_eff_projection(src, n_bins=24, band_log10f=LOG10_FGW_RANGE):
    """The drawn set's EFFECTIVE powerlaw amplitude: least-squares projection of the binned
    h_c^2(f) onto A^2 (f/f_yr)^(-4/3) over log-spaced band bins. Returns (log10_A_eff,
    per-bin anatomy). This is the number the g2a-ii FIT must reproduce."""
    f = 10.0 ** src[:, I_FGW]; h2 = 10.0 ** (2 * src[:, I_H])
    edges = np.logspace(*band_log10f, n_bins + 1)
    idx = np.clip(np.digitize(f, edges) - 1, 0, n_bins - 1)
    fc = np.sqrt(edges[:-1] * edges[1:])
    df = np.diff(edges)
    hc2 = np.zeros(n_bins)
    for i, k in enumerate(idx):
        hc2[k] += h2[i] * f[i] / df[k]                           # h_c^2 estimator per bin
    shape = (fc / F_YR) ** (-4.0 / 3.0)
    # power-weighted LSQ for A^2 in  hc2 = A^2 * shape  (weights: df/f, the band measure)
    w = df / fc
    A2 = float(np.sum(w * hc2 * shape) / np.sum(w * shape ** 2))
    return 0.5 * np.log10(A2), dict(f_bin=fc, hc2_bin=hc2, shape=shape,
                                    n_bin=np.bincount(idx, minlength=n_bins))


# ============================================================
# (3) THE FROZEN-CENSUS LEDGER -- g2c. The safety case as code.
# ============================================================
class PromoteLedger:
    """The census is drawn once and never grows. A promote (carried -> fed) is accepted ONLY
    at the source's TRUE DRAWN position, bit-exact; afterwards the M-step may move the fed
    template freely, but the ENTRY is always the truth (capstone L.2, the frozen-census
    rule). Every event is recorded for the bank."""

    def __init__(self, src_truth):
        self.truth = np.asarray(src_truth, float).reshape(-1, NP_SRC).copy()
        self.fed = np.zeros(len(self.truth), bool)
        self.events = []                     # (iter, idx) promote history

    def promote(self, idx, mean, iteration):
        if self.fed[idx]:
            raise CampaignStop(f"double promote of census member {idx} (iter {iteration}).")
        if not np.array_equal(np.asarray(mean, float), self.truth[idx]):
            raise CampaignStop(
                f"FROZEN-CENSUS VIOLATION at member {idx}, iter {iteration}: promote at a "
                f"FITTED position. B0.2 applies (search gap 3-4 orders, confident-wrong; "
                f"STORY S8.5.5) -- the campaign STOPs here, by pre-registration.")
        self.fed[idx] = True
        self.events.append((int(iteration), int(idx)))
        return self.truth[idx].copy()

    def n_resolved(self):
        return int(self.fed.sum())

    def event_array(self):
        return np.array(self.events, int).reshape(-1, 2)


# ============================================================
# (2) THE DRAINABLE BACKGROUND -- the fitted-amplitude GP
# ============================================================
class BackgroundFit:
    """Profile lnL over the background GP amplitude log10_A at a FIXED source template state.

    Uses the amortised build's own internals (kterms / P_var_inv / c_truth / frozen): for
    each A, Pinv(A) + blockdiag(c_truth) is re-factored (dsm.matrix_factor, f64 -- the
    Cholesky mandate stands) and lnL = p0 + 0.5*(FtNmy^T M^-1 FtNmy - ldP - logdet). The
    data-side terms (p0, FtNmy) are computed ONCE per (data, template) -- only the prior
    factorisation moves with A. c_truth is GWB-amplitude-independent (F^T N^-1 F against
    white+RN noise), which is exactly why this profile is affordable.

    BAND-MATCHED, and why (found by the g2a-ii gate itself, job 12645909): the population is
    band-limited to the generative box (10-32 nHz) while the 13/3 GP spans the full Fourier
    band down to 1/T. A full-band fit is dominated by the low-frequency modes the population
    does not populate and pegs at the grid floor -- the amplitude of a 13/3 powerlaw is not a
    well-defined summary of a band-limited sum. THE MODEL: the fitted amplitude scales ONLY
    the in-band GP modes; out-of-band modes stay FROZEN at the NG15 reference (comparability
    with every banked floor in the programme -- declared, not silent). Valid because the GP
    prior is mode-diagonal (Pinv = orf^-1 (x) diag(1/(A^2 S_k))): a diagonal rescale
    D Pinv D with D = A0/A on in-band rows is exactly the band-selective prior. Two identity
    checks run in-gate: (i) all-mode rescale == P_var_inv(A) directly (the scaling law);
    (ii) each pulsar block's Pinv diagonal is pairwise-monotone in the inferred mode
    frequency (the cos/sin-per-frequency ordering assumption fails loudly, not silently).
    """

    def __init__(self, amo, band_log10f=LOG10_FGW_RANGE):
        import jax.numpy as jnp
        from discovery import matrix as dsm
        self.jnp, self.dsm = jnp, dsm
        it = amo["internals"]
        self.kterms, self.Pvi = it["kterms"], it["P_var_inv"]
        self.c_truth, self.frozen = it["c_truth"], dict(it["frozen"])
        self.data_keys = it["data_keys"]
        self.truth_full = dict(it["truth_full"])
        self.theta_keys = amo["theta_keys"]
        self.npsr, self.ngp = it["npsr"], it["ngp_gwb"]
        self.akey = [k for k in self.truth_full if k.endswith("gwb_log10_A")]
        if len(self.akey) != 1:
            raise CampaignStop(f"expected exactly one gwb_log10_A key, found {self.akey}")
        self.akey = self.akey[0]
        self.a0 = float(self.truth_full[self.akey])
        # frozen reference Pinv + the in-band row mask (psr-major blocks of ngp modes;
        # mode m -> frequency (m//2 + 1)/T_span, cos/sin pairs)
        Pinv0, ldP0 = self.Pvi(self.truth_full)
        self.Pinv0 = np.asarray(Pinv0, float)
        self.ldP0_sum = float(np.sum(np.asarray(ldP0)))
        import discovery as ds
        T_span = float(ds.getspan(amo["disco_psrs"]))
        f_k = (np.arange(self.ngp) // 2 + 1) / T_span
        flo, fhi = 10.0 ** np.array(band_log10f)
        self.inband = (f_k >= flo) & (f_k <= fhi)
        if not self.inband.any():
            raise CampaignStop(f"no GP mode inside the population band "
                               f"[{flo:.2e}, {fhi:.2e}] Hz (modes reach "
                               f"{f_k.max():.2e}); the drain has no basis to live in.")
        self.mask_full = np.tile(self.inband, self.npsr)
        self.n_scaled = int(self.mask_full.sum())
        # identity check (ii): ordering -- each psr block's Pinv diag ~ orf^-1[p,p]/(A^2 S_k),
        # S ~ f^-13/3 -> diag increases with f; check pairwise (cos/sin pairs share f)
        d0 = np.diag(self.Pinv0).reshape(self.npsr, self.ngp)
        pair_mean = d0.reshape(self.npsr, self.ngp // 2, 2).mean(axis=2)
        if not np.all(np.diff(pair_mean, axis=1) > 0):
            raise CampaignStop("Pinv block diagonals are not monotone in the inferred mode "
                               "frequency: the cos/sin-per-frequency ordering assumption is "
                               "WRONG for this build -- band-matching would scale the wrong "
                               "modes. STOP (fix the mode map before any drain number).")

    def _pinv_band(self, a_log10):
        """Band-selective prior: in-band modes at amplitude a, out-of-band frozen at a0."""
        s = np.ones(self.npsr * self.ngp)
        s[self.mask_full] = 10.0 ** (self.a0 - a_log10)     # Pinv ~ 1/A^2 -> D = A0/A
        Pinv = self.Pinv0 * s[:, None] * s[None, :]
        ldP_sum = self.ldP0_sum + 2.0 * np.log(10.0) * (a_log10 - self.a0) * self.n_scaled
        return Pinv, ldP_sum

    def gate_scaling_identity(self, a_test=-14.2, tol=1e-8):
        """Identity check (i): the ALL-mode rescale must reproduce P_var_inv(a) directly."""
        pp = dict(self.truth_full); pp[self.akey] = a_test
        Pd, ld = self.Pvi(pp)
        Pd = np.asarray(Pd, float); ld_sum = float(np.sum(np.asarray(ld)))
        s = 10.0 ** (self.a0 - a_test)
        Ps = self.Pinv0 * s * s
        lds = self.ldP0_sum + 2.0 * np.log(10.0) * (a_test - self.a0) * self.npsr * self.ngp
        scale = max(np.max(np.abs(Pd)), 1e-300)
        d_p = float(np.max(np.abs(Ps - Pd))) / scale
        d_l = abs(lds - ld_sum) / max(abs(ld_sum), 1.0)
        return d_p < tol and d_l < tol, d_p, d_l

    def _data_terms(self, theta_arr, data_tuple, pmask, smask):
        """smask is REQUIRED, never defaulted: an ABSENT SMASK key means every source is in
        the template (the pre-FORGE-G identity), which for a background fit silently zeroes
        the residual. The drainable-background fit is always at an EXPLICIT feed state --
        zeros at iteration 0 (nothing resolved; the residual IS the full unresolved sum)."""
        jnp = self.jnp
        params = dict(self.frozen)
        for k, v in zip(self.theta_keys, theta_arr):
            params[k] = v
        for dk, d in zip(self.data_keys, data_tuple):
            params[dk] = d
        params["__pmask"] = pmask
        params["__smask"] = jnp.asarray(smask)
        terms = self.kterms(params)
        p0 = float(jnp.sum(terms[0]))
        FtNmy = np.asarray(terms[1]).reshape(self.npsr * self.ngp)
        return p0, FtNmy

    def _lnl_at(self, a_log10, p0, FtNmy):
        jnp, dsm = self.jnp, self.dsm
        Pinv, ldP_sum = self._pinv_band(a_log10)
        if not hasattr(self, "_C"):
            self._C = np.asarray(dsm.jsp.linalg.block_diag(*self.c_truth))
        cf = dsm.matrix_factor(jnp.asarray(Pinv + self._C))
        ld_cf = dsm.matrix_norm * float(jnp.sum(jnp.log(jnp.diag(cf[0]))))
        quad = float(FtNmy @ np.asarray(dsm.matrix_solve(cf, jnp.asarray(FtNmy))))
        return p0 + 0.5 * (quad - ldP_sum - ld_cf)

    def profile(self, theta_arr, data_tuple, pmask, smask, grid):
        """lnL over the log10_A grid + the quadratic-peak MAP and its curvature width."""
        p0, FtNmy = self._data_terms(theta_arr, data_tuple, pmask, smask)
        lnl = np.array([self._lnl_at(a, p0, FtNmy) for a in grid])
        k = int(np.argmax(lnl))
        if 0 < k < len(grid) - 1:
            a, b, c = lnl[k - 1], lnl[k], lnl[k + 1]
            d = grid[1] - grid[0]
            off = 0.5 * (a - c) / (a - 2 * b + c)
            ahat = grid[k] + off * d
            curv = (a - 2 * b + c) / d ** 2                     # d2lnL/dA2 < 0 at the peak
            sig = 1.0 / np.sqrt(max(-curv, 1e-300))
            edge = False
        else:
            ahat, sig, edge = grid[k], np.inf, True
        return dict(grid=np.asarray(grid), lnl=lnl, ahat=float(ahat), sig=float(sig),
                    edge_hit=bool(edge))


# ============================================================
# the GLACIER venue -- the T-extended problem at the drawn census
# ============================================================
def build_glacier_problem(T_label, pop, verbose=True):
    """build_ignite_problem's T-extension convention VERBATIM (extend_pulsar, span-scaled
    GP component counts), at GLACIER's OWN census: ncw = n_src, injection overridden
    per-source from the drawn population. RN frozen at the seed-0 reference draw exactly as
    make_geometry_injection does (fixed noise covariance across skies); the model GWB
    amplitude enters at the NG15 reference but the DATA never contains a stochastic GWB
    draw -- the background is the unresolved sum, which is the whole point."""
    import jax
    jax.config.update("jax_enable_x64", True)
    jax.config.update("jax_compilation_cache_dir",
                      os.environ.get("HARBOR_JAXCACHE_IN",
                                     "/home/mattm/.cache/jax_stagec_cache"))
    jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
    import trackB_b1_core as C
    import ignite as IG
    from cw_helpers import load_pulsars, build_enterprise_pta
    from stagec_fisher import (make_geometry_injection, N_COMPONENTS, LOG10_EQUAD)
    import discovery as ds

    n_src = pop["n_src"]
    t0 = time.time()
    ent_psrs, disco_psrs = load_pulsars(116)
    pta, cwb, _ = build_enterprise_pta(ent_psrs, n_src, components=N_COMPONENTS,
                                       log10_equad=LOG10_EQUAD)
    inj = make_geometry_injection(pta, ent_psrs, n_src, cwb, seed=pop["seed"])
    # override the per-source draw with the census (the drawn truth is THE truth)
    for i, name in enumerate(cwb):
        s = pop["src"][i]
        inj[f"{name}_cos_gwtheta"] = s[I_COSGT]
        inj[f"{name}_gwphi"] = s[I_GWPHI]
        inj[f"{name}_cos_inc"] = s[I_COSINC]
        inj[f"{name}_log10_mc"] = s[I_MC]
        inj[f"{name}_log10_fgw"] = s[I_FGW]
        inj[f"{name}_log10_h"] = s[I_H]
        inj[f"{name}_phase0"] = s[I_PH0]
        inj[f"{name}_psi"] = s[I_PSI]
    inj["gwb_log10_A"] = A_TARGET_LOG10
    inj["gwb_gamma"] = GAMMA_BG
    span0 = ds.getspan(disco_psrs)
    dT = float(T_label - 15)
    if dT > 0:
        diags = [IG.extend_pulsar(p, dT) for p in disco_psrs]
        span1 = ds.getspan(disco_psrs)
        rn_comp = int(round(30 * span1 / span0))
        gwb_comp = int(round(N_COMPONENTS * span1 / span0))
        if verbose:
            print(f"  EXTENSION dT={dT:.0f}yr: span {span0/YR:.2f}->{span1/YR:.2f} yr, "
                  f"rn 30->{rn_comp}, gwb {N_COMPONENTS}->{gwb_comp}", flush=True)
    else:
        rn_comp, gwb_comp = 30, N_COMPONENTS
    amo = C.build_b1_amortised(disco_psrs, n_src, inj, cwb,
                               components=gwb_comp, rn_components=rn_comp)
    if verbose:
        print(f"  build_b1_amortised: {time.time()-t0:.1f}s "
              f"(ncw={n_src}, T={T_label})", flush=True)
    return dict(amo=amo, disco_psrs=disco_psrs, ent_psrs=ent_psrs, cwb=cwb, inj=inj,
                T_label=T_label, rn_comp=rn_comp, gwb_comp=gwb_comp,
                span_s=float(ds.getspan(disco_psrs)))


# ============================================================
# GATES
# ============================================================
def gate_cpu(seed=SEED_POP_BASE, n_src=N_POP_DEFAULT, verbose=True,
             band_log10f=LOG10_FGW_RANGE):
    """g2a-i + g2b + g2c + g2d. Pure numpy + filesystem; runs anywhere at any thread count
    (nothing here consumes BLAS)."""
    p = print if verbose else (lambda *a, **k: None)
    ok = True

    p("\n-- g2a-i: draw normalisation (the drawn set sums to the NG15-class band power) --")
    pop = draw_population(seed, n_src=n_src, band_log10f=band_log10f)
    total = float(np.sum(source_band_power(pop["src"])))
    d_tot = abs(total / pop["band_power_target"] - 1.0)
    a_eff, anat = a_eff_projection(pop["src"], band_log10f=band_log10f)
    d_aeff = abs(a_eff - A_TARGET_LOG10)
    b_tot = d_tot < 1e-12
    # AMENDED 2026-07-23 (authorized; pre-registration trail): the former <0.01 projection
    # gate is RETIRED -- "the 0.01 tolerance was calibrated against a single draw (0.005)
    # that the 200-seed ensemble shows was 85th-percentile luck at the incumbent band --
    # a gate that most honest draws fail is not a gate." The projection is DEFINITIONAL
    # (it fixes A_eff, the g2a-ii fit gate's reference; that gate stays hard at 0.15 dex
    # vs A_eff-drawn) and is REPORT-ONLY + banked per draw. What remains gated:
    #   (i)  sum-power conservation, EXACT (<1e-12) -- the physical invariant, unchanged;
    #   (ii) the MANDATORY generator-bug tripwire |d log10 A_eff| < AEFF_SANITY_DEX = 0.16,
    #        the observed max over the banked 200-seed ensemble at n=256
    #        (reports/glacier_aeff_ensemble.npz, both bands) -- catches pathology, passes
    #        draw luck. Calibrated at campaign scale; at n < 200 the projection stays
    #        pure-report (small-N scatter is larger and the smoke rung is plumbing-only).
    gate_proj = n_src >= 200
    b_aeff = (d_aeff < AEFF_SANITY_DEX) if gate_proj else True
    p(f"   sum h_i^2 / target - 1        = {d_tot:.3e} (<1e-12, exact by construction)")
    p(f"   |log10 A_eff - log10 A_target| = {d_aeff:.4f} (REPORT-ONLY, banked; "
      f"{'tripwire <' + str(AEFF_SANITY_DEX) + ' dex, gated' if gate_proj else 'no tripwire at n<200'}; "
      f"projection over {len(anat['f_bin'])} bins)")
    p(f"   brightest member log10_h = {pop['src'][0, I_H]:.2f}; faintest = "
      f"{pop['src'][-1, I_H]:.2f}; n_ecc(drawn) = {int(pop['is_ecc'].sum())}/{n_src}")
    ok &= b_tot and b_aeff
    p(f"   [g2a-i] -> {'PASS' if (b_tot and b_aeff) else 'FAIL'}")

    p("\n-- g2b: resolved + residual conserved across the frontier, at EVERY position --")
    pw = source_band_power(pop["src"])
    tot = float(pw.sum())
    worst = 0.0
    for k in range(n_src + 1):                                  # frontier after census rank k
        worst = max(worst, abs((float(pw[:k].sum()) + float(pw[k:].sum())) / tot - 1.0))
    b_cons = worst < 1e-12
    p(f"   max |(P_fed + P_carried)/P_total - 1| over {n_src+1} frontier positions "
      f"= {worst:.3e} (<1e-12)")
    ok &= b_cons
    p(f"   [g2b] -> {'PASS' if b_cons else 'FAIL'}")

    p("\n-- g2c: the frozen-census ledger (promote at truth only; fitted position STOPs) --")
    led = PromoteLedger(pop["src"])
    m = led.promote(0, pop["src"][0], iteration=0)
    b_entry = np.array_equal(m, pop["src"][0]) and led.n_resolved() == 1
    try:
        led.promote(1, pop["src"][1] + 1e-9, iteration=1)       # a fitted position
        b_stop = False
    except CampaignStop:
        b_stop = True
    p(f"   promote@truth accepted, bit-exact entry: {b_entry}")
    p(f"   promote@fitted-position raises CampaignStop: {b_stop}")
    ok &= b_entry and b_stop
    p(f"   [g2c] -> {'PASS' if (b_entry and b_stop) else 'FAIL'}")

    p("\n-- g2d: provenance -- lane-tagged bank, metadata inside, cross-lane refusal --")
    stem = f"g2gate_pop_s{seed}_n{n_src}"
    path = bank_npz(stem, src=pop["src"], ecc=pop["ecc"], is_ecc=pop["is_ecc"],
                    seed=seed, n_src=n_src, sigma_logh=pop["sigma_logh"],
                    f_ecc=pop["f_ecc"], e_char=pop["e_char"],
                    a_target_log10=A_TARGET_LOG10, a_eff_log10=a_eff,
                    band_power_target=pop["band_power_target"],
                    f_bin=anat["f_bin"], hc2_bin=anat["hc2_bin"], n_bin=anat["n_bin"],
                    orientation="src row i == census rank i (brightest first); "
                                "ecc[i]/is_ecc[i] belong to src[i]; no permutation")
    z = np.load(path, allow_pickle=True)
    b_meta = all(k in z for k in ("_lane", "_host", "_cpus_affinity", "_slurm_job"))
    fake = path.replace(f"__{lane_tag()}.npz", "__otherhost_otherGPU.npz")
    np.savez(fake, dummy=1)
    try:
        bank_npz(stem, src=pop["src"])
        b_refuse = False
    except CampaignStop:
        b_refuse = True
    finally:
        os.remove(fake)
    b_read = np.array_equal(np.asarray(z["src"]), pop["src"])
    p(f"   metadata present in bank: {b_meta}; readback identical: {b_read}; "
      f"cross-lane overwrite refused: {b_refuse}")
    ok &= b_meta and b_read and b_refuse
    p(f"   [g2d] -> {'PASS' if (b_meta and b_read and b_refuse) else 'FAIL'}")
    return ok, pop, a_eff


def gate_fit(T_label, seed=SEED_POP_BASE, n_src=N_POP_DEFAULT,
             grid_lo=-15.6, grid_hi=-13.6, grid_n=41, tol=0.15, verbose=True,
             band_log10f=LOG10_FGW_RANGE):
    """g2a-ii: THE FIT GATE. Zero-noise Asimov data = the full drawn census injected as
    deterministic delays (pmask 1: pulsar terms in); NOTHING fed; the background GP
    amplitude profiled. PASS iff |log10 Ahat - log10 A_eff| < tol (0.15 dex, pre-stated)
    and the profile peaks inside the grid.

    band_log10f: generative AND fit band together (they are one declared band). At the
    incumbent default the gate CANNOT pass at any venue -- the forensic
    (glacier_g2_forensic.py; capstone session 3 SSII) showed the array carries no
    information on a 13/3 amplitude confined to 10-32 nHz (even a TRUE-A0 control pegs).
    The remedy-A band (-8.7, -7.5) is where the re-gate runs once authorised."""
    check_affinity()
    p = print if verbose else (lambda *a, **k: None)
    okc, pop, a_eff = gate_cpu(seed=seed, n_src=n_src, verbose=verbose,
                               band_log10f=band_log10f)
    if not okc:
        print("\n[g2] CPU-side gates FAILED -- not proceeding to the fit gate.", flush=True)
        return 1
    import jax.numpy as jnp
    G = build_glacier_problem(T_label, pop, verbose=verbose)
    amo = G["amo"]
    p(f"\n-- g2a-ii: the fit gate (T={T_label}, ncw={n_src}, zero-noise Asimov, "
      f"nothing resolved) --")
    theta = np.asarray(amo["theta_truth"], float)
    ones = jnp.ones(amo["npsr"])
    t0 = time.time()
    data = amo["inject_delay"](jnp.asarray(theta), ones)
    p(f"   Asimov injection: {time.time()-t0:.1f}s")
    bf = BackgroundFit(amo, band_log10f=band_log10f)
    ok_id, d_p, d_l = bf.gate_scaling_identity()
    p(f"   scaling identity (all-mode rescale == P_var_inv direct): max|dPinv|/scale "
      f"{d_p:.3e}, |d ldP| {d_l:.3e} (<1e-8) -> {'PASS' if ok_id else 'FAIL'}")
    p(f"   band-matched modes: {int(bf.inband.sum())}/{bf.ngp} per pulsar in "
      f"[{10**band_log10f[0]:.2e}, {10**band_log10f[1]:.2e}] Hz; out-of-band "
      f"frozen at the NG15 reference (declared)")
    grid = np.linspace(grid_lo, grid_hi, grid_n)
    smask0 = np.zeros(n_src)                      # NOTHING resolved: the gate's stated state
    t0 = time.time()
    prof = bf.profile(theta, data, ones, smask0, grid)
    p(f"   profile over {grid_n} amplitudes: {time.time()-t0:.1f}s "
      f"(rank {amo['npsr']*amo['internals']['ngp_gwb']})")
    d = abs(prof["ahat"] - a_eff)
    # THE GATE is the campaign venue (n >= 200): |Ahat - A_eff| < tol AND no edge. The
    # SMOKE rung (T=15/n=32) is PLUMBING-ONLY per the pre-registration (capstone SSV.2
    # item 2, authorized 2026-07-23): its job is the identity checks + a finite profile;
    # its tolerance line is REPORT-ONLY -- at n=32 the "background" is 32 discrete
    # sources and a 13/3 GP amplitude is a biased summary of that sum (measured
    # 2026-07-23: smoke |diff| 0.264 with finite sig 0.039 at the extended band -- the
    # venue HEARS now; the offset is small-N shape mismatch, not plumbing).
    gate_tol = n_src >= 200
    b_fit = ok_id and bool(np.all(np.isfinite(prof["lnl"]))) and \
        ((d < tol) and not prof["edge_hit"] if gate_tol else True)
    p(f"   log10 Ahat = {prof['ahat']:.4f} +- {prof['sig']:.4f} (profile curvature)")
    p(f"   log10 A_eff(drawn) = {a_eff:.4f};  |diff| = {d:.4f} "
      f"({'<' + str(tol) + ', PRE-STATED, gated' if gate_tol else 'REPORT-ONLY at n<200 (smoke = plumbing)'})"
      f"{'   [EDGE HIT]' if prof['edge_hit'] else ''}")
    bdef = tuple(band_log10f) == tuple(LOG10_FGW_RANGE)
    btag = "" if bdef else f"_flo{str(band_log10f[0]).replace('-', 'm').replace('.', 'p')}"
    stem = f"g2gate_fit_T{T_label}_s{seed}_n{n_src}{btag}"
    path = bank_npz(stem, grid=prof["grid"], lnl=prof["lnl"], ahat=prof["ahat"],
                    sig=prof["sig"], edge_hit=prof["edge_hit"], a_eff_log10=a_eff,
                    a_target_log10=A_TARGET_LOG10, tol=tol, T_label=T_label,
                    seed=seed, n_src=n_src, gwb_comp=G["gwb_comp"], rn_comp=G["rn_comp"],
                    band_log10f=np.array(band_log10f),
                    verdict=("PASS" if b_fit else "FAIL"),
                    orientation="lnl[i] is the profile at grid[i]; ahat is the quadratic "
                                "peak; the verdict is cut from these columns")
    p(f"   banked -> {path}")
    p(f"   [g2a-ii] -> {'PASS' if b_fit else 'FAIL'}")
    print(f"\n=== GLACIER g2 GATES: {'ALL PASS' if (okc and b_fit) else 'FAIL'} ===",
          flush=True)
    return 0 if (okc and b_fit) else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["gate", "fitgate", "gate2"])
    ap.add_argument("--t", type=int, default=30)
    ap.add_argument("--seed", type=int, default=SEED_POP_BASE)
    ap.add_argument("--n", type=int, default=N_POP_DEFAULT)
    ap.add_argument("--flo", type=float, default=LOG10_FGW_RANGE[0],
                    help="log10 f_lo of the declared band (remedy A: -8.7)")
    ap.add_argument("--fhi", type=float, default=LOG10_FGW_RANGE[1])
    ap.add_argument("--grid-lo", type=float, default=-15.6)
    ap.add_argument("--grid-hi", type=float, default=-13.6)
    a = ap.parse_args()
    band = (a.flo, a.fhi)
    if a.mode == "gate":
        ok, _, _ = gate_cpu(seed=a.seed, n_src=a.n, band_log10f=band)
        print(f"\n=== g2 CPU GATES: {'ALL PASS' if ok else 'FAIL'} ===", flush=True)
        return 0 if ok else 1
    if a.mode == "gate2":
        return gate_ladder(n_src=a.n, band_log10f=band)
    return gate_fit(a.t, seed=a.seed, n_src=a.n, band_log10f=band,
                    grid_lo=a.grid_lo, grid_hi=a.grid_hi)


if __name__ == "__main__":
    sys.exit(main())
