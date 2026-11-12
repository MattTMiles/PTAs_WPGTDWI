# wheel_demo_v4 — build report

**Built:** 2026-08-07. **Source:** `wheel_demo_v3.ipynb`, which is **untouched** —
verified by md5 before and after every build step (`d897b1e2802692d30b1deb074a3ac94c`).
v3 remains the frozen presentation fallback.

**Deliverable:** `wheel_demo_v4.ipynb` — 24 cells, 12 code cells, **0 errors**, executes
top to bottom in **3 min 18 s** (1 min 33 s before the seed sweep was added).

---

## 1. Headline numbers, old → new

| quantity | v3 | v4 (hybrid) |
|---|---|---|
| population | 340,000 drawn, treated as the whole sky | **5.5×10⁷ modelled**: 14,793 discrete + analytic sea |
| realized A / NG15 | 6.24e-15 (0.98×) | **6.77e-15 (1.06×)** |
| λ_fore | 0.28 | **0.31** this draw; **0.16 median** of 10 |
| counts slope γ_eff | 1.9 (pooled fit) | **1.57 (per bin, V24's definition)** |
| occupancy at 10 nHz | 170/bin — **271× low** | **2.7×10⁴/bin — 0.59× audit** |
| mass function | hard ceiling at 10¹⁰ M⊙ | **Schechter cutoff** at 10¹⁰ M⊙ |
| CW upper limit | 2 drawn sources **violate** it | **0 violate** (2 rejected and redrawn) |
| faint sea | did not exist | **12 %** of the in-band background, unsubtractable |
| **K ladder** (60/80/100/120/150/200) | 9 / 13 / 16 / 24 / 28 / 36 | **10 / 13 / 20 / 25 / 30 / 40** |
| **ignition flip point** | between **80 and 100** psr | **between 80 and 100 psr — unchanged** |
| 100-psr wheel | 16 → 36 → 44, 80 % cert turn 0 | **20 → 44 → 44, 97 % cert turn 0** |
| 120-psr wheel | 24 → 54 | **25 → 51** |
| wheel multiplier | ~2.2× | **~2.0×** (2.04 / 2.03 / 1.98 at 120/150/200) |
| SNR ≥ 8 bracket | 3 / 5 / 6 / 6 / 6 / 7 | **0 / 0 / 1 / 1 / 1 / 3** |
| V4c real per-pulsar priors | 24 → 30 → 49 → 54, 42 %→100 % | **25 → 40 → 51, 56 %→100 %** |
| V4d strict gate | K = 30, 37 % certified | **K = 49, 92 %** |
| V4e calibrated gate | K = 41, 78 % | **K = 51, 100 %** |
| V6 computed bound | 61 % tighter | **79 % tighter** (banked V26: 38.9 % at R = 22) |

**Direction check against the brief's expectation.** The UL cut and the softened slope
were expected to push K *down*; occupancy realism and cleaning to push it back *up*. Net:
K rose modestly (16 → 20 at the working array) because the Schechter cutoff removed the
few pathological loud sources while the realistic occupancy added many more moderate ones.
The two effects very nearly cancelled at the flip point, which is why it did not move.

---

## 2. Population scorecard — five checks

| # | check | measured | target | verdict |
|---|---|---|---|---|
| 1 | realized A | 6.769e-15, ratio **1.058** | 6.4e-15 | **PASS** |
| 2 | λ_fore | **0.313** (this draw) | 0.31 audited | **PASS — but see caveat** |
| 3 | loudest discrete vs CW UL | 1.14e-14 vs 3.26e-14 at 1.8 nHz, ratio **0.348** | < 1 | **PASS** |
| 4 | counts slope γ_eff | **1.569** per bin | 1.56 audited | **PASS** |
| 5 | occupancy at 10 nHz | 2.71e4 vs 4.60e4, factor **0.59** | within a few × | **PASS** |

**Caveat on check 2, printed in the notebook itself.** Over 10 draws λ_fore has median
**0.162**, range 0.088–0.313. Seed 3 — the notebook's — sits at the very **top** of that
range. The PASS is a favourable draw; the typical draw gives 0.16 against the audited 0.31.
**λ_fore is the one calibration target v4 does not close**, and the scorecard cell prints
the spread and says so explicitly so it cannot be read as representative.

Why it misses: our population is built from physical mass and distance distributions, so
the in-bin strain function steepens below the knee. V24's 0.31 comes from drawing the
*entire* bin from a single Pareto. This is a genuine model difference, not a tuning
residual, and it is worth the monitoring agent's attention — it suggests the audited
λ_fore = 0.31 may be partly a property of the Pareto parametrisation.

A is also consistently ~5 % high (median 1.055, range 0.939–1.118).

---

## 3. The three K ladders (Part B, cleaning, at the 100-psr working array)

| wheel | K ladder | final K | stalls? |
|---|---|---|---|
| coherence only | 20 → 44 → 44 | 44 | **yes** |
| + cleaning, banked η [V26] | 20 → 47 → 49 → 49 | 49 | **yes** |
| + cleaning, η = 1 (perfect-subtraction bound) | 20 → 58 → 62 → 62 | 62 | **yes** |

η reached 0.293 → 0.389 (banked, R-dependent). A fell 6.77e-15 → 6.28e-15 (banked η) and
→ 5.23e-15 (η = 1). Total curve depth went 0.190 dex (coherence alone) → 0.222 → 0.300 dex.
The two channels multiply, as specified.

---

## 4. Does ignition still land between the 80- and 100-pulsar eras?

**On the presentation draw, yes.** K_rec = 6 at 80 psr (below combs\* = 7.8) and 10 at
100 psr (above), unchanged from v3 despite the population being rebuilt from scratch.

**But this is a single-draw statement and it does not survive a seed sweep — see the
follow-up section below.** Across 13 draws the flip appears in only 4; seed 3 is the
kindest of the 13. The claim as originally written here was over-confident.

---

## 5. Verdicts on the two pre-registered predictions

**(a) "Cleaning extends the wheel but does not prevent the stall — the faint sea is
unsubtractable and is now structural." → HELD.** Cleaning extends final K from 44 to 49
(banked η) and to 62 (η = 1), and **all three wheels still stall**. The faint sea carries
12 % of the in-band background and no amount of discrete-source subtraction touches it.

Worth flagging: at η = 1 the wheel extends by 41 % in K, so the cleaning channel is not
negligible — it is the largest single upgrade in v4 after the population itself.

**(b) "h_split robustness: K at the working array changes by ≤ 1 across ±0.3 dex." → HELD,
ΔK = 0.** The paired test gives K = 16 at all three splits (−0.3 / 0.0 / +0.3 dex) while
the discrete count varies 44,815 / 15,153 / 4,873. The split is bookkeeping, not a hidden
knob. **No raise of the split was needed.**

*This required fixing the test itself* — see decision 5 below.

---

## 6. Autonomous decisions taken, with rationale

1. **NG15 CW upper limit: constructed shape, published normalisation.** Only one point is
   quoted in the accessible paper text ("8×10⁻¹⁵ at 6 nHz", Agazie+2023 ApJL 951 L50,
   confirmed by web search). Rather than invent a published curve, the *shape* is taken
   from a hasasia NG15-like array (67 psr, 15 yr, 100 ns — real code) and normalised
   through that anchor. Labelled in the code and honesty box as **not** the published curve.
2. **γ_eff measured per frequency bin, not pooled.** V24 defines it as the Pareto index of
   the *in-bin* strain function. My first pooled fit gave 2.6–2.9 and looked like a hard
   miss; measured per bin the same population gives 1.52–1.60. Pooling mixes the
   f-dependence of h and biases the slope high. **This single correction is what let target
   (ii) close** — without it I would have wrongly driven MC_CUT to an unphysical 3×10¹⁰ M⊙.
3. **MC_CUT = 10¹⁰ M⊙, D_MAX = 4000 Mpc, MC_MIN = 10⁸ M⊙.** Chosen by scanning against the
   priority order. 10¹⁰ M⊙ is a defensible Schechter *scale* (exponential suppression, not
   a hard cap; component BHs ~1.15×10¹⁰, within the observed range). D_MAX = 4000 gave
   γ_eff = 1.556 ± 0.050 — exactly on target — at the cost of λ_fore; per the stated
   priority (γ_eff exact **before** λ_fore close) that is the correct trade.
4. **N_TOT_BOOST = 1.44.** The ensemble-mean normalisation put the *realized* A at 0.83 of
   target, because bin power is heavy-tailed (γ_eff < 2) so a median draw falls below the
   mean. Target (i) says *realized* A, so I measured the median over 10 seeds and corrected.
   Stated as a named constant with its reason in the code.
5. **The h_split robustness test was re-designed.** My first version redrew the population
   at each split, which changes the RNG stream entirely — it reported ΔK = 2 and "FAIL"
   while actually measuring Poisson noise between independent realisations. The test now
   draws **once**, deep (−0.6 dex), and **re-cuts the same sources** at each split. That is
   the only valid comparison, and it gives ΔK = 0.
6. **Kept V4c / V4d / V4e and the prior experiment.** My first build replaced v3's cell 10
   wholesale, which would have silently deleted the real-priors panel, both gated wheels and
   the prior-sensitivity experiment, orphaning ~150 lines of markdown. I patched `run_wheel`
   in place instead. The brief asked to *add* cleaning, not to remove the gate work.
7. **Seed-sensitivity print added to the scorecard.** Seed 3 gives λ_fore = 0.313, which
   reads as a perfect match to the audit. The median draw gives 0.16. Leaving that PASS
   unqualified would have overstated the calibration, so the cell now prints the spread and
   an explicit NOTE. This was not requested; it is required for the scorecard to be honest.
8. **A-grid interpolation for the cleaning channel.** Cleaning lowers A, which lowers the
   hasasia curve; rebuilding hasasia every iteration would dominate runtime. Six curves per
   era on a log-A grid, interpolated in log-log. The curve is *not* a simple rescaling of A
   (GWB-dominated at low f, white-noise dominated at high f), which is why a grid is needed
   rather than a scalar factor.
9. **Candidate set widened.** With cleaning the curve reaches deeper than
   `hc0 × 10^(−DEPTH_DEX_MAX)`, so the original geometry-candidate set raised a `KeyError`.
   It now covers the deepest curve reachable *with* cleaning.
10. **Faint sea excluded from detection/combs.** Per the brief: the sea enters the realized
    background and the comb-jitter scale, but detection and the V3/V5 machinery see discrete
    sources only.

---

## 7. Wall clock and compute

- **Notebook execution:** 1 min 33 s (14 min CPU — hasasia is threaded).
- **Part A calibration:** ~35 min of scripted scans (well inside the 2-hour time-box).
- **Total build:** ~2.5 h including the two rebuilds forced by decisions 5 and 6.
- **GPU: not used.** As the brief anticipated, essentially all of this is numpy quadrature
  and hasasia CPU work. The faint sea is closed-form (no MC at all), and the calibration
  loop is ~10 draws of ~15,000 sources — milliseconds each. There was no point in the
  session where a GPU would have helped, and I did not restructure code to manufacture one.

---

## 8. What a reviewer should look at first

1. **λ_fore = 0.16 median vs 0.31 audited** — the one open calibration gap, and possibly a
   statement about the audit's Pareto parametrisation rather than about our population.
2. **The per-bin vs pooled γ_eff distinction** (decision 2) — it changes the answer by 1.0
   in γ and would silently corrupt any future recalibration.
3. **The cleaning channel is large** (K 44 → 62 at η = 1). If the banked V26 η is the right
   efficiency to reuse here, the wheel is meaningfully stronger than v3 suggested; if not,
   that is the number to challenge.
4. **Still absent:** intrinsic pulsar red noise (spin noise, DM, jitter). The only red noise
   is the GWB, so every array still looks better than its real counterpart. This is now the
   largest un-modelled effect in the notebook.

---

# Follow-up: seed robustness of the ignition verdict

**Added 2026-08-07.** The v4 headline — *"ignition lands between 80 and 100 pulsars"* —
was measured on seed 3, the draw where λ_fore = 0.313. The scorecard shows the model's
typical sky is λ_fore ≈ 0.16. This section asks whether the verdict survives typical
skies, and whether it leans on the disputed census number.

**Method.** Every calibrated knob frozen exactly as in v4 (normalisation, MC_CUT,
MC_SLOPE, D_RANGE, h_split rule, CW-UL cut); only the population seed varies. Per-draw A
refit is the v4 pipeline's own. **13 seeds** (seed 3 included as a labelled row), all six
arrays, counting-rule wheel **with cleaning at banked η**. New notebook cell sits after
the V2 ladder; new figure `fig_seed_robustness.png`.

## The table

| seed | A/NG15 | λ_fore | K ladder (60/80/100/120/150/200) | rung | K₈₀ | K₁₀₀ | final K |
|---|---|---|---|---|---|---|---|
| **3** | 1.058 | **0.313** | 10 / 13 / 20 / 25 / 30 / 40 | 100 | 13 | 20 | 49 |
| 4 | 1.082 | 0.162 | 9 / 10 / 12 / 16 / 19 / 35 | 120 | 10 | 12 | 15 |
| 5 | 1.051 | 0.184 | 10 / 14 / 16 / 19 / 24 / 37 | 100 | 14 | 16 | 52 |
| 6 | 1.082 | 0.163 | 8 / 9 / 11 / 16 / 20 / 28 | 120 | 9 | 11 | 11 |
| 7 | 0.944 | 0.088 | 8 / 10 / 14 / 22 / 30 / 35 | 120 | 10 | 14 | 14 |
| 8 | 0.996 | 0.111 | 5 / 6 / 13 / 16 / 23 / 27 | 120 | 6 | 13 | 13 |
| 9 | 1.118 | 0.216 | 9 / 13 / 14 / 22 / 24 / 32 | 100\* | 13 | 14 | 41 |
| 10 | 1.026 | 0.217 | 6 / 11 / 14 / 18 / 26 / 34 | 120 | 11 | 14 | 15 |
| 11 | 0.939 | 0.106 | 11 / 20 / 24 / 27 / 34 / 43 | 80 | 20 | 24 | 59 |
| 12 | 1.083 | 0.116 | 14 / 14 / 14 / 19 / 25 / 36 | 120 | 14 | 14 | 14 |
| 13 | 1.047 | 0.154 | 15 / 18 / 24 / 28 / 33 / 40 | 80 | 18 | 24 | 53 |
| 14 | 0.948 | 0.114 | 11 / 16 / 19 / 23 / 30 / 37 | 80 | 16 | 19 | 46 |
| 15 | 0.970 | 0.093 | 8 / 9 / 16 / 22 / 28 / 40 | 100 | 9 | 16 | 50 |

\* seed 9 fails on turn 0 at 100 psr but certifies later — see the definition note below.

## Headline fractions

- **ignite by the 100-psr era: 7 / 13 = 54 %**
- **flip in (80, 100] — the v3/v4 claim: 4 / 13 = 31 %**
- rung distribution: **80 psr ×3, 100 psr ×4, 120 psr ×6**
- median rung **100 psr**; on turn-0 firing alone the median is **120 psr**
- seed 3's λ_fore ranks **13 / 13** — the single kindest draw in the sweep

## The correlation that decides the ticket

| relation | Spearman ρ | p |
|---|---|---|
| λ_fore vs ignition rung | **+0.065** | **0.833** |
| λ_fore vs K at 100 psr | −0.084 | 0.785 |
| **K at 100 psr vs ignition rung** | **−0.88** | **0.0001** |
| λ_fore vs A ratio | +0.670 | 0.012 |

Mann-Whitney on λ_fore between skies that ignite by 100 psr and those that do not:
**p = 0.836**. Median λ_fore of the igniting group is 0.154 against 0.139 for the rest —
*higher*, i.e. the opposite direction to the hypothesis, and not significant.

## VERDICT

**Neither (a) nor (b) as pre-written — a third outcome, and I am reporting it rather than
forcing a fit.**

> **Ignition is NOT seed-robust, but it does NOT track λ_fore.** Only 7/13 skies certify
> by the 100-psr era and only 4/13 show the 80–100 flip; the median sky needs 100 psr and
> the turn-0 flip needs 120. Yet ignition is uncorrelated with λ_fore (ρ = +0.07,
> p = 0.83; Mann-Whitney p = 0.84). The calendar leans on **realization scatter in the
> detectable source count**, not on the disputed census number — so **the λ_fore ticket is
> NOT load-bearing for the ignition claim.**

**Consequences for the two open tickets, which now separate cleanly:**

- **λ_fore ticket — downgraded for this purpose.** Whatever P(D) delivers on λ_fore, it
  will not move the ignition rung. λ_fore still matters for the population's fidelity to
  the census (it is v4's one missed calibration target), but the calendar does not rest on it.
- **A new ticket is created, and it is the load-bearing one:** the ignition rung is set by
  **K at the candidate array**, which needs **K ≥ 16** (rec ≥ combs\* = 7.8) and is achieved
  by only 6/13 draws at 100 psr. What a measurement must deliver is not a λ_fore value but
  a tighter constraint on **source counts at the detection bar**. Because there is no
  λ_fore dependence, "what λ_fore does the flip need" has no answer — the requested number
  does not exist.

**The wheel is bimodal, which is why the rung is sharp.** Final K at the 100-psr array is
**41–59 when it certifies and 11–15 when it does not** — nothing in between. The array
either catches enough sources to start, or it never moves.

## Judgment calls in this follow-up

1. **Two ignition definitions, both reported.** "Certifies on turn 0" (the v4 *FIRES*
   label) and "certifies at any iteration". They differ because **cleaning does not require
   certification** — you subtract a recovered source from the background whether or not its
   distance is pinned. So a sky can fail turn 0, clean its way deeper, and certify later.
   Exactly one seed (9) is rescued that way. I used *ever-certifies* for the headline as the
   physically meaningful definition, and report the turn-0 median (120 psr) alongside.
   This was not anticipated in the brief and is a genuine finding about the cleaning channel:
   **cleaning can ignite the wheel, not merely extend it.**
2. **Cell placement moved.** The brief asked for the cell after the scorecard; it needs
   `ERAS`/`WORK`, which V2 defines *later*. Rather than duplicate the era list (drift risk),
   I placed it immediately after the V2 ladder — narratively still correct: you see the
   ladder, then ask whether it is robust.
3. **One absolute A-grid shared across seeds.** hasasia curves depend on the population only
   through A, so a single 6-era × 8-amplitude grid (48 builds) serves all 13 seeds instead
   of ~600 rebuilds. Interpolation is log-log in A, the same machinery v4 already uses for
   cleaning.
4. **13 seeds, not 10.** Marginal cost, and n = 13 gives the Spearman test enough power that
   a null result is worth reporting rather than merely inconclusive.
5. **Gated wheel not swept**, as the brief allowed. The gate is strictly harder than the
   counting rule, so **every rung here is a lower bound on the array the gated wheel needs** —
   the real picture is at least this pessimistic.

## Runtime

Sweep adds ~1 min 45 s to the notebook (48 hasasia builds + 13 × 6 wheel runs); full
notebook now executes in **3 min 18 s**, 0 errors, 24 cells. No GPU — as before, this is
numpy quadrature and hasasia CPU work.

## Markdown honesty updates applied

- **Honesty box**: the ignition claim now carries "4 of 13 draws", names seed 3 as the
  kindest, and records the λ_fore null.
- **V4 figure suptitle**: now reads "SEED 3 IS THE KINDEST OF 13 DRAWS: the 80-100 psr flip
  shown here appears in only 4/13 skies; the rung is 120 psr in 6/13."
- **V7 summary** and **appendix**: both now state the rung distribution and seed status
  rather than the single-draw claim.
