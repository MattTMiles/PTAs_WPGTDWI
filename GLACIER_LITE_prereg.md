# GLACIER-LITE — PRE-REGISTRATION

Report-only campaign. Branch `glacier_lite`. Budget cap 25 GPU-hr.
Written 2026-07-27, **before any cell of this campaign was scored**.
Nothing in this campaign arms any protocol step or enters closure claims.

---

## 1. THE PRE-LOGGED PREDICTION BLOCK (verbatim, as supplied in the brief)

> "Predicted: at 20-50 ns with 3 T2-conditioned seeds, lap 1 certifies
> >=1 true pulsar per realisation; certifications compound over laps
> 2-4 and saturate when the readable sub-array is exhausted; residual
> power decreases monotonically with lap index; the anchor cell
> reproduces sub-pc single-shot distance widths for pulsars within
> ~1 kpc. Failure of compounding (flat certifications after lap 1)
> falsifies the loop claim in the easiest available regime and is a
> headline result to be reported with equal prominence."

This block is reproduced byte-identically as `PREDICTION` in
`hpc_harbor/glacier/gl_lite.py` and is banked into **every** `gl_lite_*` npz
produced by this campaign, so no bank can be read without it.

## 1b. THE SECONDARY-ARM PREDICTION (verbatim, Matt 2026-07-27)

> "ignition here would bound the harvest's structure-limited conclusion to sub-SKA
> noise; non-ignition strengthens it."

Banked as `PREDICTION_RENORM` into every cell.

## 1c. THE PART-3 PREDICTION (verbatim, amendment 2026-07-27)

> "Predicted: C at saturation rises monotonically with prior tightening;
> at 1/30 widths with 20 ns noise, C > 0.8, with the uncertified
> remainder audibility-limited (low pulsar-term SNR from geometry),
> not ambiguity-limited. If C plateaus well below 1 even at 1/30 and
> the autopsy shows ambiguity-limited pulsars remaining, the
> completeness limit is NOT prior-purchasable and the closure claim
> must be scoped to the readable sub-array."

Banked as `PREDICTION_COMPLETE` into every cell.

## 1d. THE CAMPAIGN AS AMENDED (budget cap now 35 GPU-hr)

| arm | part | noise | reals | laps | notes |
|---|---|---|---|---|---|
| **PRIMARY** (frozen POP, A_eq **-13.22**) | 1 loop | 50, 20 ns | 4 | 8 | the loop demonstration |
| PRIMARY | 2 anchor | 20 ns | 1 | single-shot | 5 sources at full truth |
| PRIMARY | 3 completeness | 20 ns | 2 | <=12, saturation | prior widths {1, 1/3, 1/10, 1/30} |
| **SECONDARY** (NG15-renormalised, A_eq **-14.60**) | 1 loop only | 20 ns | 2 | 8 | no anchor, no Part 3 |

The primary arm stays **population-identical to the realistic campaign banks** so that
only noise and knowledge change (single-variable control). The secondary arm is banked
under `*_ng15` stems and is **never merged** with the primary arm's cells.

**Budget-cut order (amended):** cut the SECONDARY arm's realisations 2 -> 1 before
touching the primary. Within Part 3, cut the prior-width factors from four to three
`{1, 1/10, 1/30}` before cutting realisations or laps. Within Part 1, realisations
4 -> 2 before laps, and drop 50 ns before 20 ns. Every cut reported.

## 1e. PART-3 CONVENTIONS FIXED BEFORE SCORING

- **The prior-width floor.** `sigma_a := max(factor * sigma_a_real, 3 * dL[a])`. All
  three prior columns (`feather`, `script`, `lit`) are scaled together, because
  `build_EV` sizes the fringe-sampling grid from `max` over the three — scaling only
  `lit` would leave the sampled comb at its old width and the tightening would be
  cosmetic.
  *Reconciliation of the floor to K:* the brief glosses the `3*dL` floor as "never
  sub-fringe (K < ~3)". Under the house's own definition
  `K_counted = 2*floor(3*sigma/dL)`, a prior at the floor gives **K_counted = 18**, not
  3; the informal count `sigma/dL` gives 3. The operative rule implemented is the
  literal one — `sigma >= 3*dL` — and the resulting K distribution (min/median/max over
  the 116 pulsars) is banked and reported per cell, so both counts are visible and
  nothing is silently reinterpreted.
- **Truth moves with the prior (declared).** The house truth draw
  (`ignite.draw_true_distances_tier`) is prior-weighted over the sampled fringes, so a
  tighter prior necessarily concentrates the drawn truth toward `L0`. Across prior-width
  factors the knowledge changes *and* the drawn truth concentrates with it; the factors
  are **not** the same universe re-analysed. The alternative — holding truth fixed while
  shrinking sigma — would place the truth many new-sigma outside its own prior and is
  unphysical.
- **The autopsy rule (stated before scoring).** For a pulsar that saturates
  uncertified, with `dlnL_req = lnK + 0.578` and `dlnL_avail = 0.5 * SNR_p^2` where
  `SNR_p` is the optimal SNR of the *pulsar term* of the loudest census member against
  this venue's own white covariance:
  - **audibility-limited** iff `dlnL_avail < dlnL_req` — the pulsar term cannot clear
    its own trials bar even with perfect fringe identification. More prior cannot help.
  - **ambiguity-limited** otherwise — the matched power is there and the comb is what is
    unresolved. This is the class a tighter prior *can* buy.
- **Saturation** = 2 consecutive laps adding zero *new* true certifications.
- **Cumulative counts are over DISTINCT pulsars.** A pulsar re-certified on a later lap
  is not a new certification; `C = |ever-truly-certified| / 116`.

## 2. WHAT IS FROZEN AND UNTOUCHED

- **Criterion v2.2**, verbatim from `glacier_loop.run_cell`:
  `cert = (dlnL_det > max(lnK + TRIALS_NAT, floor)) & (q_max > QBAR)`
  with `TRIALS_NAT = 0.578`, `QBAR = 0.9`.
- **Floor adoption**: `recut_surface.adopt` — Gumbel only where the null
  zero-fraction <= 0.20, else empirical q95 + bootstrap sd. Zero-fraction is a
  required banked column.
- **Frontier-v2** (`ratio < 0.5` AND `dlnL_feed > 0`). The v1 path is not reachable.
- **The kill step (D2 rigidity, R1/R2) runs SCORE-AND-LOG ONLY.** Its
  false-negative gate has not passed (capstone S4.23.1). It must not execute
  against any certification. Enforced structurally: the kill verdict is a banked
  column and is never read by the certification path.
- All existing banks, populations and the frozen gate are read-only.

## 3. WHAT IS IDEALIZED (the only two axes)

1. **Noise.** White only, uniform across the array, two levels {50 ns, 20 ns}.
   No red noise, no GWB in the data.
2. **Source knowledge.** The 3 loudest members conditioned at tier **T2**
   (`oracle_ignition.TIER_FREE`): sky (`cos_gwtheta`, `gwphi`) + orbital period
   (`log10_fgw`) + chirp mass (`log10_mc`) pinned at drawn truth; `cos_inc`,
   `log10_h`, `phase0`, `psi` free. All other members are not conditioned.

## 4. WHAT IS KEPT REAL (explicitly not idealized)

- The **116 real pulsar sky positions** (the feather set's own `pos`).
- Their **real published distance priors** (tier `lit`, `A2.load_real_priors`).
- The **fringe ambiguity is the real one** — `dL[a]` is the minimum mode spacing
  over the population at the real sky geometry; `K_counted` is unchanged.

## 5. PRE-DECLARED STOP CONDITIONS (from the brief)

Report, do not improvise, if any of these fire:
- a floor re-cut that comes back degenerate;
- disagreement between disk and the brief on naming, machinery, population
  structure;
- the 25 GPU-hr cap.

## 6. PRE-DECLARED BUDGET-CUT ORDER (from the brief)

If the design does not fit 25 GPU-hr:
1. realisations per noise level 4 -> 2 (**before** cutting laps — lap depth is the point);
2. drop the 50 ns level **before** touching the 20 ns level.
Every cut is reported. No silent truncation.

## 7. DEVIATIONS AND DISK-VS-BRIEF NOTES DECLARED IN ADVANCE

Recorded here before scoring so they cannot be read as post-hoc:

- **(a) Branch name.** The brief says branch `MM_playground`; the checkout was on
  `oracle_ignition`. Both are the *same commit* `ed78a1c` (`git rev-list --count`
  both directions = 0). No material difference. New code goes on `glacier_lite`
  as instructed.
- **(b) "16-member population".** GLACIER's own census is 256 members
  (`glacier_pop.N_POP_DEFAULT = 256`; every banked `gl1_*` cell has `n_slot`
  256 or 287). The **frozen 16-member population is the house's trackB/B1 one** —
  `trackB_b1_core.POP = dict(ncw=16, seed=3000, population=(3, -13.25, -14.25))`
  with `N_LOUD = 3` — which is the population SPARK, CHORUS, GEO, SIREN and FORGE
  all ride, and whose "3 loudest" is a declared house constant. This campaign uses
  **that** population, re-used from its frozen seed, not redrawn, and not GLACIER's
  256-census. Its native venue is T=15 (span 22.15 yr), which is also the only
  venue at which `b1_L_gwb` is banked.
- **(b2) THE POPULATION IS NOT "NG15-CONSISTENT LOUDNESS" — measured, before scoring.**
  The brief describes the frozen 16-member population as being "at NG15-consistent
  loudness". It is not. Its declared strain ladder is 3 sources at `log10_h = -13.25`
  and 13 at `-14.25`, giving `sum h_i^2 = 9.898e-27` against an NG15 band-power target
  (`A = -14.6`, band 10–31.6 nHz) of `1.728e-29`:

  | quantity | value |
  |---|---|
  | band power ratio | **572.8 x** |
  | amplitude ratio | **23.9 x  (+1.38 dex)** |
  | A_equivalent of the population | **-13.22** (vs NG15 -14.6) |

  For scale, GLACIER's own ladder called a *single* member at `-13.25` a "declared
  super-NG15" rung; this population has three of them. The loud members induce a
  residual amplitude `h/(2 pi f)` of **895 ns at 10 nHz** (283 ns at 31.6 nHz) against
  the 20–50 ns noise of this campaign.

  This is a genuine contradiction between the brief and disk, and it is recorded here
  *before* any cell is scored so it cannot be read as post-hoc. It is a **labelling**
  contradiction, not an operational one: there is exactly one standard frozen 16-member
  population, and the brief's own instruction is to *"reuse a banked population; do not
  redraw"*. Renormalising it to NG15 would be redrawing it. This campaign therefore uses
  the frozen population as banked, and **every headline number carries the scope line:
  the loop is being run at ~24x NG15 band power, not at NG15.** If the loop ignites here,
  that is ignition in a loud, idealized corner — not a statement about an NG15 sky.

- **(c) "All other members fully free".** The machinery has no "free but present"
  state for an unconditioned member: a census member is either FED (present in the
  template) or CARRIED (`H_ABSENT = -30`, living in the fitted GP background). The
  13 unconditioned members are therefore CARRIED at lap 0 and become available to
  the loop only by being promoted through frontier-v2 — which is precisely the
  mechanism the loop demo is measuring. Read "fully free" as "not conditioned, and
  available for the frontier to resolve".
- **(d) White-noise implementation.** `white_scale` in `NoiseDrawer.draw` scales the
  **data** only and leaves the likelihood covariance at its built value. Using it
  alone would mis-specify the noise by the same factor and forfeit the entire
  sensitivity gain the idealization is supposed to buy. This campaign therefore
  **rebuilds the venue** so that the target white level is in *both* the data and
  the likelihood: the feather set already carries uniform `toaerrs = 1e-6 s` and
  `LOG10_EQUAD = -6.0`, i.e. a uniform white rms of `sqrt(2) * 1e-6 s = 1414 ns`;
  the lite venue sets `toaerrs = equad = target/sqrt(2)` so that
  `N_p = efac^2 (toaerr^2 + equad^2) = target^2` exactly.
- **(e) Red noise.** Removed from the data (components = `("white",)`) and driven
  to numerically negligible amplitude in the model, with the RN basis shapes kept
  intact so no downstream factorisation changes shape.
- **(f) The GWB model block is retained** — it is the drain's instrument
  (`BackgroundFit` profiles its amplitude, which is the residual-power readout).
  The data never contains a stochastic GWB draw in this machinery; the background
  IS the unresolved source sum. So "no GWB" is already true of the data by
  construction.

## 7b. BUILD GATES LANDED (CPU lane, jobs 12803564 + 12804324; pre-scoring)

Run on the free CPU lane because the A100 entitlement was saturated. Deliberately does
NOT call `estep_per_target` (the XLA-CPU `vm.max_map_count` hazard), so the E-step, the
floor re-cut, frontier-v2 promotion and the kill scorer remain unexercised until the GPU
smoke.

| gate | primary (50 ns) | secondary (20 ns, renorm, pw 1/30) |
|---|---|---|
| white level reached the LIKELIHOOD covariance | **50.000 ns** PASS | **20.000 ns** PASS |
| frozen-population strain ladder == banked POP | PASS | PASS |
| NG15 renormalisation (one common dex shift) | n/a | A_eq **-14.600** (shift -1.379) PASS |
| banked `b1_L_gwb.npz` loads STRICT (`fp=71677a810cbc7187`) | PASS | PASS |
| prior floor `sigma >= 3*dL` on all 116 | n/a | PASS (floored on 8/116) |
| in-band GP modes at band (-8.0,-7.5) | 16 of 28 PASS | 16 of 28 PASS |
| BackgroundFit scaling identity | PASS | PASS |
| drain profile finite, not at grid edge | PASS (after fix) | PASS |

**A defect the gates caught, fixed before any production run.** The first pass returned
`ahat -13.000 +- inf, edge_hit True`: the drain's amplitude grid topped out at -13.0,
*below* the frozen population's own A_equivalent of -13.22, so at nothing-fed (where the
entire census is unresolved background) the profile pegged. Grid widened to
`linspace(-17, -11, 121)`.

**DRAIN ZERO-POINT (carry with every residual-vs-lap number).** On the renormalised
venue, with nothing fed, the drain reads **-14.501 +- 0.013** against a population truth
of **-14.600** -- biased HIGH by ~0.10 dex, which is ~7.7x its own quoted error. This is
expected (a band-matched Gaussian GP standing in for 16 discrete sinusoids is not an
unbiased estimator of their summed amplitude) and is NOT a failure, but it means the
residual-power curve is a **relative** measurement against this zero-point. The absolute
`a_bg` value must not be quoted as the background amplitude.

**Measured scales, for reading Part 3's autopsy.**

| quantity | value |
|---|---|
| `dL` (fringe spacing) median | 0.22 pc |
| `K_counted` median, real `lit` priors | **6381** |
| `K_counted` min/median/max at pw 1/30 | **18 / 212 / 1882** |
| `sigma_lit` median, 1 -> 1/30 | 226.8 -> 7.6 pc |
| pulsar-term SNR of loudest member (renormalised venue, 20 ns) | min 0.04, median **7.38**, max 19.87 |

The `K_min = 18` confirms the S1e reconciliation empirically: a prior at the `3*dL` floor
gives `K_counted = 18` under the house definition, not the "~3" the brief's gloss implies.

**Honest caveat on the autopsy statistic.** `dlnL_avail = 0.5*SNR_p^2` is a matched-filter
CEILING assuming perfect fringe identification. The criterion's actual statistic
`dlnL_det` is the *gap between the best and second-best fringe*, which is a different and
generally much smaller quantity -- a large SNR does not guarantee a large gap when
fringes are near-degenerate. The autopsy rule is therefore an ATTRIBUTION device ("the
matched power is present; whether the comb resolves is the open question"), not a
forecast of certification. Only the measured `dlnL_det` settles it.

## 8. RELEVANT PRIOR RESULT (context, not a prediction)

`ORACLE_IGNITION.md` (this same branch lineage, 2026-07-27) measured **276 cells,
zero certifications** across N in {1,2,3,5,8} x tier {T1,T2,T3}, and concluded that
*knowledge alone cannot buy ignition* — the binding levers being in-band frequency
and prior width, at **realistic** noise. GLACIER-LITE moves the one lever ORACLE
held fixed: the noise. It is therefore not a repeat, and ORACLE's verdict does not
predetermine it.
