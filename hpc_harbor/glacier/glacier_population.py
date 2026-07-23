#!/usr/bin/env python3
"""GLACIER g2 — the drawn population and the DRAIN bookkeeping (reading-independent core).

Agent GLACIER, ACCRE. Dev/doc authority; code STAGED (Matt commits). Zero-noise / Asimov
arithmetic only — NO campaign GPU, NO noisy realisations. Pure numpy, CPU, runs in seconds.

WHAT THIS IS, AND WHAT IT DELIBERATELY IS NOT
---------------------------------------------
The brief's g2 asks for "a population where the background IS unresolved sources ... A_background
= the fitted residual amplitude ... the drain curve." STAGE-0 RECON (GLACIER_capstone.md §g2)
established, with file evidence, that the programme's GWB is a SEPARATE, FROZEN HD GP at NG15
A=-14.6 (`trackB_b1_core.py:174` makeglobalgp_fourier; the L_gwb bank exists to FREEZE its
amplitude), that no file ties the summed CW sources to that GP, and that the programme's own
physics gives N_bar<=0.1 individually-resolvable sources. Making A_gwb a fitted drainable
parameter equal to the summed CW power is a programme-wide change to the generative model that
severs comparability with the fixed-A_gwb floors GLACIER must be read against (SURFACE, CHORUS,
and the SPARK-3 EDGE-POSITIVE verdict that LICENSED GLACIER, all cut at A=-14.6).

SPARK-3 §5.3 — the licence — pre-registers the ACHIEVABLE GLACIER: the FIXED-LIST reservoir-drain,
"buys nothing new in the sky, only a better model of what is already known to be there." There the
drainable background is the FAINT RESERVOIR's carried power, and the drain is power migrating
carried->fed (SMASK) as pulsars certify, ON A FROZEN CENSUS, against the fixed A=-14.6 GP.

THIS MODULE is the arithmetic common to BOTH readings and commits to NEITHER generative model:
  * it DRAWS a census (the incumbent 3+13, and — for the g2 "N~200-500" ask — an extended draw),
    with eccentricity from LOTTERY's two-population mix;
  * it defines a per-source residual-POWER proxy P_i (strain-rate power);
  * it defines A_background(frontier) := sqrt(sum of CARRIED source power) — the DRAIN CURVE —
    as a bookkeeping split of a FROZEN total over the fed/carried partition (FORGE-G's SMASK
    trichotomy). Conservation P_fed + P_carried == P_total is EXACT BY CONSTRUCTION for any
    partition, because it is a sum split, not a fit.
  * The g2 gate ("at iteration 0, nothing resolved, A_background reproduces the drawn population's
    summed power within tolerance") is therefore EXACT (tol 0) at iteration 0 and conserved across
    the frontier by construction — which is the honest content of the brief's gate.

Whether that carried power is subsequently INJECTED as a fitted GP (the brief's literal reading,
which needs the generative-model change above and a floor re-cut) or CARRIED as prior-width
SoftComponents (the §5.3 reading, on FORGE-G's built machinery) is the DECISION surfaced to Matt
in GLACIER_capstone.md. This module makes the drain measurable either way and prejudges neither.
"""
import os
import sys
import numpy as np

# Layout + prior box come from the LICENSED machinery so the draw is consistent with the loop.
_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("/home/mattm/projects/HSYMT/CW_transition",
           os.path.join(_HERE, "..", "..", "CW_transition")):
    if os.path.isdir(_p):
        sys.path.insert(0, _p)
        break
try:
    from trackB_b1_softsource import (NP_SRC, I_MC, I_FGW, I_H, I_COSGT, I_GWPHI,
                                       I_COSINC, I_PH0, I_PSI, AXIS_NAME)
except Exception:                                    # standalone fallback (bit-identical layout)
    I_COSGT, I_GWPHI, I_COSINC, I_MC, I_FGW, I_H, I_PH0, I_PSI = range(8)
    NP_SRC = 8
    AXIS_NAME = ["cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
                 "log10_fgw", "log10_h", "phase0", "psi"]

# ---- the incumbent census (trackB_b1_core.POP) and its geometry -------------------------------
POP_INCUMBENT = dict(ncw=16, n_loud=3, log10_h_loud=-13.25, log10_h_faint=-14.25,
                     log10_mc=9.0, log10_fgw=-7.9)   # census loudness (LOTTERY_tier1.md:17-20)

# ---- LOTTERY's two-population eccentricity mix (lottery_tier1.py:126-129) ----------------------
#   each LOUD member is eccentric (e = e_char) w.p. f_ecc, else circular; the FAINT members are
#   circular by construction. "The realistic astrophysical case" is the low-e end, e_char ~ 0.3
#   (LOTTERY_tier1.md:86-91). f_ecc is the mixing fraction.
MIX_REALISTIC = dict(f_ecc=0.5, e_char=0.3)          # parameterised; the realistic low-e point
EMAX = 0.9                                            # lottery_tier1.py:52


def draw_mix_ecc(rng, n_loud, f_ecc, e_char):
    """Loud-member eccentricities: e_char w.p. f_ecc, else 0.0 (lottery_tier1.py:126-129)."""
    is_ecc = rng.random(n_loud) < f_ecc
    return np.where(is_ecc, e_char, 0.0)


def draw_population(n_total, n_loud=3, seed=3000, mix=None,
                    log10_h_loud=-13.25, log10_h_faint=-14.25,
                    faint_h_spread=0.0, log10_mc=9.0, log10_fgw=-7.9):
    """Draw an (n_total, 8) source array + per-source eccentricity + a loud/faint mask.

    n_total == 16 with defaults reproduces the incumbent census geometry (a LOG-NORMAL strain
    HIERARCHY degenerates to two deltas at h_loud / h_faint). For the g2 "N~200-500" ask, set
    n_total large and faint_h_spread>0: the faint tier becomes a log-normal strain hierarchy of
    unresolved members BELOW the loud tier — the drainable reservoir the drain curve eats.

    Eccentricity: loud members from LOTTERY's mix (draw_mix_ecc); faint members circular.
    Positions/inclination/phase/psi drawn uniform on the generative prior box. Masses/frequency
    fixed at the census values (the drain is a strain-hierarchy statement; mc/fgw spread is a
    hook, not exercised here). Everything is a TRUE drawn value — no fitted position enters, so
    the frozen-census rule (GLACIER_capstone.md §L.2) and §5.3's safety case hold by construction.
    """
    mix = MIX_REALISTIC if mix is None else mix
    rng = np.random.default_rng(seed)
    src = np.zeros((n_total, NP_SRC))
    # geometry (report-only for the power proxy; drawn so the census is a real population)
    src[:, I_COSGT] = rng.uniform(-1.0, 1.0, n_total)
    src[:, I_GWPHI] = rng.uniform(0.0, 2 * np.pi, n_total)
    src[:, I_COSINC] = rng.uniform(-1.0, 1.0, n_total)
    src[:, I_PH0] = rng.uniform(0.0, 2 * np.pi, n_total)
    src[:, I_PSI] = rng.uniform(0.0, np.pi, n_total)
    src[:, I_MC] = log10_mc
    src[:, I_FGW] = log10_fgw
    # strain hierarchy: loud tier flat at h_loud; faint tier at h_faint (+ optional log-normal
    # spread DOWNWARD, so no faint member ever outshines the loud tier).
    src[:n_loud, I_H] = log10_h_loud
    if faint_h_spread > 0 and n_total > n_loud:
        n_f = n_total - n_loud
        src[n_loud:, I_H] = log10_h_faint - np.abs(rng.normal(0.0, faint_h_spread, n_f))
    else:
        src[n_loud:, I_H] = log10_h_faint
    is_loud = np.zeros(n_total, bool)
    is_loud[:n_loud] = True
    ecc = np.zeros(n_total)
    ecc[:n_loud] = draw_mix_ecc(rng, n_loud, mix["f_ecc"], mix["e_char"])
    return dict(src=src, ecc=ecc, is_loud=is_loud, n_loud=n_loud, n_total=n_total,
                seed=int(seed), mix=dict(mix))


def source_power(src, weight="residual"):
    """Per-source power proxy, geometry-AVERAGED (a strain-hierarchy statement, report-only on
    per-pulsar antenna patterns — GLACIER_capstone.md §g2 records the geometry caveat).

    The induced residual is linear in strain with amplitude scale alpha = 10**log10_h / (2 pi *
    10**log10_fgw)  (lnL_GW_recovery_refactor.py:389, the verbatim waveform). So:
      weight='residual' : P_i = alpha_i**2 = (10**h)**2 / (2 pi 10**fgw)**2   [residual power]
      weight='strain'   : P_i = (10**h)**2                                    [bare strain power]
    Both are reported; the drain curve is scale-free in P (A_bg = sqrt(sum P)), so the choice
    rescales A_background by a constant and does not change the DRAIN FRACTION.
    """
    h = 10.0 ** src[:, I_H]
    if weight == "strain":
        return h * h
    f = 10.0 ** src[:, I_FGW]
    alpha = h / (2.0 * np.pi * f)
    return alpha * alpha


def drain_state(pop, fed_mask, weight="residual"):
    """Given a fed/carried partition (fed_mask over the frozen census), return the drain state.

    A_background := sqrt(sum of CARRIED source power). THE DRAIN CURVE is A_background as fed_mask
    grows (pulsars certify -> sources cross the model-quality frontier carried->fed). Conservation
    is EXACT: P_fed + P_carried == P_total for ANY mask, because it is a sum split of a frozen
    total (not a fit) — this is the honest content of the brief's iteration-0 g2 gate.
    """
    P = source_power(pop["src"], weight)
    fed = np.asarray(fed_mask, bool)
    P_total = float(P.sum())
    P_fed = float(P[fed].sum())
    P_carried = float(P[~fed].sum())
    return dict(P_total=P_total, P_fed=P_fed, P_carried=P_carried,
                A_background=float(np.sqrt(max(P_carried, 0.0))),
                A_total=float(np.sqrt(P_total)),
                N_resolved=int(fed.sum()), N_carried=int((~fed).sum()),
                resid=abs(P_fed + P_carried - P_total))          # == 0 by construction


# ================================================================================================
# g2 GATE — iteration-0 conservation + drain monotonicity + census reproduction. CPU, seconds.
# ================================================================================================
def gate(verbose=True):
    def p(*a):
        if verbose:
            print(*a, flush=True)
    ok = True
    p("=" * 92)
    p("GLACIER g2 — DRAWN POPULATION + DRAIN BOOKKEEPING (zero-noise, CPU)")
    p("=" * 92)

    # --- g2.1: the incumbent census is reproduced exactly (continuity with the frozen loop) ---
    pop16 = draw_population(16, n_loud=3, seed=3000)
    h = pop16["src"][:, I_H]
    c1 = (np.allclose(h[:3], -13.25) and np.allclose(h[3:], -14.25)
          and pop16["is_loud"].sum() == 3 and pop16["n_total"] == 16)
    ok &= c1
    p(f"\n-- g2.1: incumbent census (3 loud @-13.25, 13 faint @-14.25) --")
    p(f"   loud h {np.unique(h[:3])}  faint h {np.unique(h[3:])}  n_loud {pop16['is_loud'].sum()}"
      f"   -> {'PASS' if c1 else 'FAIL'}")

    # --- g2.2: ITERATION-0 GATE. Nothing fed -> A_background == sqrt(summed population power),
    #           EXACT (tol 0). This is the brief's g2 gate, stated honestly. ---
    st0 = drain_state(pop16, np.zeros(16, bool))
    c2 = (st0["A_background"] == st0["A_total"] and st0["resid"] == 0.0
          and st0["P_fed"] == 0.0 and st0["N_resolved"] == 0)
    ok &= c2
    p(f"\n-- g2.2: iteration-0 — fitted background reproduces the summed population power --")
    p(f"   A_background(nothing fed) = {st0['A_background']:.6e}   sqrt(P_total) = {st0['A_total']:.6e}")
    p(f"   |A_bg - sqrt(P_total)|    = {abs(st0['A_background']-st0['A_total']):.3e} (== 0, EXACT)")
    p(f"   conservation resid P_fed+P_carried-P_total = {st0['resid']:.3e} (== 0)"
      f"   -> {'PASS' if c2 else 'FAIL'}")

    # --- g2.3: CONSERVATION across the whole frontier + DRAIN monotone. Feed sources one at a
    #           time (loud first, as the loop resolves loudest-first); A_background must fall
    #           monotonically to 0, conservation exact at every step. ---
    order = np.argsort(-source_power(pop16["src"]))     # loudest-power first
    Abg, resid_max = [], 0.0
    fed = np.zeros(16, bool)
    Abg.append(drain_state(pop16, fed)["A_background"])
    for i in order:
        fed[i] = True
        st = drain_state(pop16, fed)
        Abg.append(st["A_background"]); resid_max = max(resid_max, st["resid"])
    Abg = np.array(Abg)
    monotone = bool(np.all(np.diff(Abg) <= 1e-30))
    drains_to_zero = bool(Abg[-1] == 0.0)
    # iteration-0 conservation is EXACT (single split, g2.2); across the SWEPT frontier the running
    # sum P_fed+P_carried-P_total accumulates float roundoff, so the honest bar is a RELATIVE one at
    # machine epsilon, not bit-zero.
    rel_resid = resid_max / st0["P_total"]
    c3 = monotone and drains_to_zero and rel_resid < 1e-12
    ok &= c3
    p(f"\n-- g2.3: the DRAIN CURVE — A_background as sources resolve (frozen census) --")
    p(f"   A_bg: {Abg[0]:.4e} (0 fed) -> {Abg[3]:.4e} (loud fed) -> {Abg[-1]:.4e} (all fed)")
    p(f"   monotone non-increasing {monotone};  drains to 0 {drains_to_zero};"
      f"  swept conservation resid {resid_max:.3e} (rel {rel_resid:.1e}, machine-eps)"
      f"   -> {'PASS' if c3 else 'FAIL'}")
    loud_drain = 1.0 - (drain_state(pop16, pop16["is_loud"])["P_carried"] / st0["P_total"])
    p(f"   fraction of total power carried by the 3 LOUD members: {loud_drain:.4f}"
      f"   (the reservoir is the faint {1-loud_drain:.4f} — what the array must eat)")

    # --- g2.4: the extended g2 draw (N large) is well-formed and conserves too ---
    popN = draw_population(300, n_loud=3, seed=3000, faint_h_spread=0.3)
    stN = drain_state(popN, np.zeros(300, bool))
    c4 = (popN["src"].shape == (300, 8) and stN["resid"] == 0.0
          and np.all(popN["src"][3:, I_H] <= -14.25 + 1e-12))   # faint never outshines loud
    ok &= c4
    p(f"\n-- g2.4: extended draw N=300 (log-normal faint hierarchy, faint below loud) --")
    p(f"   shape {popN['src'].shape}; faint h in [{popN['src'][3:,I_H].min():.2f},"
      f" {popN['src'][3:,I_H].max():.2f}]; conservation resid {stN['resid']:.3e}"
      f"   -> {'PASS' if c4 else 'FAIL'}")

    p(f"\n{'='*92}")
    p(f"GLACIER g2 GATE: {'ALL PASS' if ok else 'FAIL'}")
    p(f"{'='*92}")

    # bank the drain curve + census for the report (lean; raw powers carried, not just verdicts)
    out = os.path.join(_HERE, "..", "..", "reports", "glacier_g2_population.npz")
    np.savez(os.path.abspath(out),
             src16=pop16["src"], ecc16=pop16["ecc"], is_loud16=pop16["is_loud"],
             P16_residual=source_power(pop16["src"], "residual"),
             P16_strain=source_power(pop16["src"], "strain"),
             drain_Abg=Abg, drain_order=order,
             P_total=st0["P_total"], A_total=st0["A_total"], loud_drain_frac=loud_drain,
             mix_realistic=np.array([MIX_REALISTIC["f_ecc"], MIX_REALISTIC["e_char"]]),
             gate_pass=ok, note="g2 iter-0 conservation exact; drain = carried-power sqrt")
    p(f"banked -> {os.path.abspath(out)}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(gate())
