# ANCHOR — the REALISM anchor. Floors against noise the analysis did not specify.

Agent ANCHOR, ACCRE, 2026-07-12/13. Tag `criterion-v2.1` @ `6bec3d6`.
Untracked working report. **Nothing committed.** Results: `ANCHOR_results/`. Code: `hpc_harbor/anchor/`.

**Two verdicts, and they are different in kind:**

1. **THE LADDER: CONSISTENT AT EVERY RUNG.** Noise mis-specification at realistic amplitudes does
   **not** move the detection floor. **The onset surface does not shift.** 18/40 cells "inflated"
   (a coin flip), median Δq95 = **−0.18 nat**, 1/40 cells significant (2 expected by chance) —
   *while the realised noise rms rises by up to 32 %.*
2. **A DEFECT IN THE CRITERION'S OWN FLOOR ESTIMATOR, found on the way.** The D2 Gumbel fit is
   **biased LOW — permissive — wherever the offender distribution is dominated by its zero point
   mass**, by up to **2.8×** (0.85 nat fitted vs 2.40 nat empirical at h = −13.25 lit, where
   **93 % of nulls have no offender at all**). This is S6.2.2's caveat, at three times the
   severity, and it is not a caveat: it is a bias, and it errs in the dangerous direction.

---

## 0. WHAT CHANGED SINCE THE STOP

The first ANCHOR brief was refuted at Task 1: the repo has no real residuals (the 116-psr feather
set is a mock — telescope `AXIS`, one 1440 MHz channel — and its residual column *is* the injected
CW+CURN realisation `b20_cw_curn_r0`). That forensics stands and is banked
(`anchor_data_forensics.npz`); §7 carries the doc-correction drafts it generated.

The amended premise: **the anchor is REALISM, not provenance.** Regenerate noise-only nulls on the
*same* array with noise the frozen analysis does **not** model, and ask what that does to the floor.

**The seam is exact.** `NoiseDrawer.draw` is the *only* thing the ladder replaces. `B1Split` (`sp`),
the `FringeTables`, `ln K_counted`, the EM distance priors, the tiers, the geometry ensembles are
the banked frozen objects, untouched. **The mismatch IS the experiment.**

---

## 1. GATES — BOTH PASS

**g1 — R0 reproduces the banked modelled-noise behaviour. PASS, BIT-FOR-BIT.**
All **80** banked `ig_nullN_*_T15_*` realisations replayed through the ANCHOR harness under R0:

| column | `dlnL_det` | `lnK` | `qmax` | `mapk` | `ptrue` | `on_true` | **offender** |
|---|---|---|---|---|---|---|---|
| max abs diff | **0.0** | **0.0** | **0.0** | **0.0** | **0.0** | **0.0** | **0.0** |

80/80 bit-identical (`anchor_g1.npz`). R0 is a literal passthrough of `P["nd"].draw(seed)` — not a
reimplementation — so this gates the *harness*, not the arithmetic. (It also silently confirms
[[noise-draw-thread-count-hazard]]: all jobs ran at `cpus-per-task=8`.)

**g2 — resume drill. PASS, AND STRONGER THAN THE HOUSE DRILL.**
A production array task was `scancel`'d mid-flight. Then: (i) **all 1800 banked npz verified intact**
— no truncated write from the kill; (ii) 12 npz deleted with their SHA-256 recorded; (iii) the
resubmitted shard logged **`already banked: 1788; to run: 12`** and re-ran exactly 12; (iv) all
**12/12 refilled files are BYTE-IDENTICAL** to the originals. Resume is not merely *skipping* — it is
**deterministic**. Final bank verified: **7200/7200 npz, 0 corrupt**.

**A free third gate.** R0's floor reproduces the programme's own `h`-law without being asked to:
empirical q95 (lit) = 2.40 / 10.42 / 30.05 / 78.80 nat at h = −13.25 / −13.0 / −12.75 / −12.5.
Between the two loudest cells that is **`floor ∝ h^1.67`** — against S5.3.1's independently measured
**`floor_fN ∝ h^1.66`**. The harness reproduces the floor mechanism from scratch, at a baseline
(T = 15) where that law was never fitted.

---

## 2. THE LADDER AS BUILT

`T = 15` — the **native, unextended array** (`ignite.py:207`, `dT = T_label − 15`; `extend_pulsar` is a
strict no-op at `dT ≤ 0`). 8 cells = h {−13.25, −13.0, −12.75, −12.5} × tier {lit, vlbi}.
**N = 150 nullN per (rung, cell)** (D2's onset-cell sizing), skies cycled over `SKIES = [−1,0,1,2,3]`,
seeds `7_xxx_xxx` (disjoint from IGNITE 1–4.xM and IGNITE-2 5–6.xM). **7 200 realisations**, ~0.5 s
each on an A100-40.

| rung | what is added to the DRAW (the analysis never changes) | realised noise rms | ×R0 |
|---|---|---|---|
| **R0** | control — the modelled spectra exactly, `nd.draw(seed)` verbatim | **1.993 µs** | 1.000 |
| **R1** | +per-pulsar RN (log10 A, γ) mis-specification; GWB amplitude × U(0.5, 1.5) | 2.030 µs | 1.018 |
| **R2** | +unmodelled DM power-law GP per pulsar | 2.207 µs | 1.107 |
| **R3** | +non-Gaussian tails (excess white + 2 % shot outliers at 6σ) | **2.634 µs** | **1.321** |
| R2o | *isolation:* R0 + DM only | 2.169 µs | 1.088 |
| R3o | *isolation:* R0 + tails only | 2.450 µs | 1.229 |

The mis-specification is **redrawn per null**, so each floor is *marginalised over the
measurement-uncertainty ensemble*: "the floor you get when the true noise differs from the model by
as much as the measurements cannot exclude."

**R1's amplitudes.** σ(log10 A_red) = 0.25 dex, σ(γ_red) = 0.60, ±3σ clipped — representative NG15-era
marginal widths (Agazie et al. 2023, ApJL 951 L10, the noise-budget paper). **Stated plainly: these
are representative published *widths*, not the published *per-pulsar posteriors*** — that table is not
in-repo. R1's median rms barely moves (×1.018) **by construction**, because a symmetric log
perturbation is median-preserving; the *per-pulsar* amplitudes swing by up to ±0.75 dex (a factor of
**5.6 in power**), which is the mis-specification that matters. The array's modelled noise is
red-dominated (RN 1.56 µs vs white 1.41 µs, fitted log10 A_rn = −13.85, γ = 3.39), so R1 perturbs the
*dominant* channel.

**R2's amplitude, and why it is not at the literature's loud edge.** log10 A_dm ~ N(−13.5, 0.3),
γ_dm ~ N(2.0, 0.5), clipped to the NG15 DMGP support. At the loud edge (log10 A_dm = −12.9) the DM
channel alone would carry **~2.1 µs** — as loud as the *entire* noise model — and R2 would have been
measuring our amplitude choice, not realism. At −13.5 it carries **0.52 µs** (26 % of the modelled
noise): real, subdominant, defensible. The achieved rms is **computed at build time and banked**, not
asserted. *Sensitivity, stated not measured:* DM power ∝ A², so a 0.6-dex louder DM would roughly
triple the DM rms.

### 2.1 THE LADDER'S STATED CEILING — R2 is not a test of chromaticity

**The array is single-frequency.** All 30 225 TOAs sit at 1440.0 MHz, so the chromatic factor
`(f_ref/f_obs)²` is **a constant** — measured, 0.9450 to 0.9454. Chromaticity is therefore
**unidentifiable by construction**: a DM GP enters the data as an extra *red* process carrying the DM
power law's shape, and no analysis could separate it from an achromatic one. **R2 tests UNMODELLED
CHROMATIC-BAND POWER, not chromaticity.** Splitting the channel would mean rebuilding the array
(new TOAs, new design matrix, real DM columns), which the brief excluded. This is the ceiling the
redirect anticipated, and it is the ladder's binding limitation. *A real multi-frequency array would
lift it* — see the REAL-ARRAY PENDING in §7.

*(ECORR is likewise degenerate with EQUAD here: one TOA per epoch. R3 folds it into the excess-white
channel and says so.)*

---

## 3. THE COMPARISON ROW — THE LADDER TABLE

Floors on the offender statistic (max over pulsars of the fringe `dlnL` gap, among cells passing
layer 1 ⊕ layer 3; 0.0 if none). **Gumbel = the ADOPTED D2 estimator** (α = 0.05, fit error
2.80·β̂/√N). **EMP q95 = the empirical (1−α) quantile**, bootstrap error — *valid regardless of the
zero point mass*. `*` marks cells where the zero-fraction exceeds 20 % and **the Gumbel is fitting a
point mass, not a tail**; there the verdict is issued on the empirical quantile (see §4).

| rung | cell | zero % | Gumbel (D2) | ±fit | Δ | **EMP q95** | ±bs | **Δ vs R0** | P(off>0) | MW p | **VERDICT** |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 | −13.25 lit `*` | 93.3 | 0.845 | 0.064 | — | **2.395** | 1.342 | — | 0.067 | — | CONTROL |
| R1 | −13.25 lit `*` | 88.7 | 1.162 | 0.086 | +0.32 | 2.695 | 0.857 | **+0.30** | 0.113 | 0.180 | **CONSISTENT** |
| R2 | −13.25 lit `*` | 93.3 | 0.740 | 0.056 | −0.10 | 1.982 | 0.942 | **−0.41** | 0.067 | 0.967 | **CONSISTENT** |
| R3 | −13.25 lit `*` | 88.7 | 1.261 | 0.093 | +0.42 | 2.577 | 0.995 | **+0.18** | 0.113 | 0.175 | **CONSISTENT** |
| R0 | −13.25 vlbi `*` | 80.0 | 2.669 | 0.191 | — | **4.982** | 0.361 | — | 0.200 | — | CONTROL |
| R1 | −13.25 vlbi `*` | 78.7 | 2.349 | 0.168 | −0.32 | 3.927 | 0.329 | **−1.06** | 0.213 | 0.947 | DEFLATED −1.06 |
| R2 | −13.25 vlbi `*` | 72.0 | 3.316 | 0.231 | +0.65 | 4.862 | 0.370 | **−0.12** | 0.280 | 0.196 | **CONSISTENT** |
| R3 | −13.25 vlbi `*` | 54.7 | 5.583 | 0.365 | +2.91 | 5.338 | 0.292 | **+0.36** | **0.453** | **0.0000** | **CONSISTENT** (but see §5) |
| R0 | −13.00 lit `*` | 57.3 | 7.273 | 0.487 | — | **10.416** | 0.837 | — | 0.427 | — | CONTROL |
| R1 | −13.00 lit `*` | 58.7 | 6.953 | 0.468 | −0.32 | 10.661 | 0.585 | **+0.25** | 0.413 | 0.762 | **CONSISTENT** |
| R2 | −13.00 lit `*` | 48.7 | 9.256 | 0.600 | +1.98 | 10.177 | 0.400 | **−0.24** | 0.513 | 0.079 | **CONSISTENT** |
| R3 | −13.00 lit `*` | 54.7 | 8.635 | 0.573 | +1.36 | 10.796 | 0.740 | **+0.38** | 0.453 | 0.432 | **CONSISTENT** |
| R0 | −13.00 vlbi | 19.3 | 11.016 | 0.615 | — | 9.297 | 0.850 | — | 0.807 | — | CONTROL |
| R1 | −13.00 vlbi `*` | 28.7 | 10.651 | 0.628 | −0.36 | 9.287 | 0.494 | −0.01 | 0.713 | 0.279 | **CONSISTENT** |
| R2 | −13.00 vlbi `*` | 24.0 | 12.393 | 0.718 | +1.38 | 11.169 | 0.813 | +1.87 | 0.760 | 0.566 | INFLATED +1.87 |
| R3 | −13.00 vlbi | 18.0 | 12.467 | 0.694 | **+1.45** | 11.650 | 0.895 | +2.35 | 0.820 | 0.099 | INFLATED +1.45 |
| R0 | −12.75 lit | 5.3 | 31.493 | 1.586 | — | 30.046 | 1.160 | — | 0.947 | — | CONTROL |
| R1 | −12.75 lit | 5.3 | 29.849 | 1.530 | −1.64 | 25.964 | 1.738 | −4.08 | 0.947 | 0.285 | **CONSISTENT** |
| R2 | −12.75 lit | 3.3 | 30.888 | 1.513 | −0.61 | 29.804 | 1.045 | −0.24 | 0.967 | 0.823 | **CONSISTENT** |
| R3 | −12.75 lit | 4.7 | 29.526 | 1.464 | −1.97 | 29.587 | 2.289 | −0.46 | 0.953 | 0.488 | **CONSISTENT** |
| R0 | −12.75 vlbi | 0.0 | 28.032 | 1.260 | — | 30.687 | 2.716 | — | 1.000 | — | CONTROL |
| R1 | −12.75 vlbi | 0.0 | 26.649 | 1.229 | −1.38 | 29.081 | 1.088 | −1.61 | 1.000 | 0.165 | **CONSISTENT** |
| R2 | −12.75 vlbi | 0.7 | 27.364 | 1.241 | −0.67 | 29.155 | 1.334 | −1.53 | 0.993 | 0.452 | **CONSISTENT** |
| R3 | −12.75 vlbi | 0.0 | 26.882 | 1.191 | −1.15 | 28.158 | 2.337 | −2.53 | 1.000 | 0.638 | **CONSISTENT** |
| R0 | −12.50 lit | 0.0 | 74.728 | 3.190 | — | 78.796 | 4.344 | — | 1.000 | — | CONTROL |
| R1 | −12.50 lit | 0.0 | 82.499 | 3.720 | +7.77 | 72.673 | 4.253 | −6.12 | 1.000 | 0.395 | INFLATED +7.77 (Gumbel) |
| R2 | −12.50 lit | 0.0 | 72.248 | 3.079 | −2.48 | 74.516 | 4.354 | −4.28 | 1.000 | 0.557 | **CONSISTENT** |
| R3 | −12.50 lit | 0.0 | 75.320 | 3.300 | +0.59 | 77.668 | 2.850 | −1.13 | 1.000 | 0.637 | **CONSISTENT** |
| R0 | −12.50 vlbi | 0.0 | 77.681 | 3.430 | — | 78.625 | 3.395 | — | 1.000 | — | CONTROL |
| R1 | −12.50 vlbi | 0.0 | 80.532 | 3.604 | +2.85 | 78.909 | 3.838 | +0.28 | 1.000 | 0.586 | **CONSISTENT** |
| R2 | −12.50 vlbi | 0.0 | 76.791 | 3.363 | −0.89 | 81.677 | 4.262 | +3.05 | 1.000 | 0.905 | **CONSISTENT** |
| R3 | −12.50 vlbi | 0.0 | 76.273 | 3.368 | −1.41 | 77.509 | 3.198 | −1.12 | 1.000 | 0.723 | **CONSISTENT** |

*(Isolation rungs R2o/R3o are in `anchor_ladder.npz` and in §6; they carry the same message.)*

### THE PRE-REGISTERED VERDICT, PER RUNG

**R1 (spectral mis-specification): CONSISTENT.**
**R2 (unmodelled chromatic-band power): CONSISTENT.**
**R3 (non-Gaussian tails): CONSISTENT.**

**The onset surface does not shift. Not by any amount worth quoting.** Pooled over all 40
(rung × cell) comparisons:

| statistic | value | what a real inflation would look like |
|---|---|---|
| cells inflated (Δq95 > 0) | **18 / 40** | ≫ 20/40 |
| median Δq95 | **−0.179 nat** | ≫ 0 |
| mean Δq95 | **−0.658 nat** | ≫ 0 |
| Mann-Whitney p < 0.05 | **1 / 40** (2 expected by chance) | many |

The per-cell "INFLATED"/"DEFLATED" entries above are **not a signal**: they change sign between
neighbouring cells of the same rung (R1 is +7.77 at −12.50 lit and −1.64 at −12.75 lit), and they sit
at ~1–2σ against fit errors that are themselves 3–5 nat at the loud cells. **They are the calibration
noise S5.3.7 already told us to expect** — and S5.3.7's rule applies to them exactly as written.

---

## 4. THE SECOND FINDING — THE D2 ESTIMATOR IS BIASED LOW WHERE THE NULL IS MOSTLY SILENT

This was not in the brief. It fell out of the control arm, and it matters more than the ladder does.

The offender statistic is **0.0 whenever a realisation has no cell passing layer 1 ⊕ layer 3**. At
faint `h` that is *most* realisations. Measured zero-fractions in **R0** (the modelled noise, the
banked convention):

| cell | zero-fraction | Gumbel floor (D2, adopted) | **empirical q95** | Gumbel / empirical |
|---|---|---|---|---|
| h = −13.25 lit | **93.3 %** | **0.845 ± 0.064** | **2.395** | **0.35 ×** |
| h = −13.25 vlbi | 80.0 % | 2.669 ± 0.191 | 4.982 | 0.54 × |
| h = −13.00 lit | 57.3 % | 7.273 ± 0.487 | 10.416 | 0.70 × |
| h = −13.00 vlbi | 19.3 % | 11.016 ± 0.615 | 9.297 | 1.18 × |
| h = −12.75 lit | 5.3 % | 31.493 ± 1.586 | 30.046 | 1.05 × |
| h = −12.75 vlbi | 0.0 % | 28.032 ± 1.260 | 30.687 | 0.91 × |
| h = −12.50 lit | 0.0 % | 74.728 ± 3.190 | 78.796 | 0.95 × |
| h = −12.50 vlbi | 0.0 % | 77.681 ± 3.430 | 78.625 | 0.99 × |

**Where the zero-fraction is small (≤ 20 %), the Gumbel and the empirical quantile agree to within a
few per cent — the estimator is sound and the fit errors are honest.** Where the point mass dominates,
the Gumbel fit is pulled down toward it and **understates the α = 0.05 bar by up to 2.8×**, by
**1.55 nat at (−13.25, lit)** and **2.31 nat at (−13.25, vlbi)** — gaps of **24σ and 12σ against its own
quoted fit error.** The fit error is not just wrong, it is *confidently* wrong: a Gumbel fitted to a
93 % point mass reports ±0.064 nat.

**This errs in the dangerous direction.** Detection is `dlnL > max(ln K, floor)`. A floor that is too
low is **too permissive** — it lets pure-noise offenders through. Every faint-`h` cell in the onset map
is calibrated by this estimator.

**S6.2.2 already saw the edge of this** ("at that cell 45 % of nullN realisations have NO offender — a
point mass at zero the Gumbel does not model… the fit is serviceable, but the zero-fraction travels
with the number"). At 45 % it *was* serviceable. **At 57 %, 80 % and 93 % it is not, and "serviceable"
is not a property that extrapolates.** The honest statement is now:

> **CONVENTION (proposed): the D2 Gumbel floor is valid only where the nullN zero-fraction is
> ≲ 20 %. Above that, quote the empirical (1−α) quantile with a bootstrap error, and bank the
> zero-fraction beside it. The zero-fraction is a REQUIRED column, not a caveat.**

This does not retract anything. **IGNITE-2's two calibrated cells are safe**: (−12.75, 30, lit) has
zero-fraction 0.00, and (−13.25, 30, vlbi) has 0.45 with a measured Gumbel-vs-q95 agreement of 0.5 nat.
The exposure is the **uncalibrated faint-h cells** of the SURFACE grid, where nobody has ever looked at
the zero-fraction — and where, on this evidence, the floors are **too low**.

---

## 5. THE ONE THING REALISM *DOES* MOVE — THE BODY, NOT THE TAIL

Exactly **1 of 40** cells is significant under Mann-Whitney — and it is not marginal:
**R3 at (h = −13.25, vlbi), p = 0.0000.** Its floor barely moves (**Δq95 = +0.36 ± 0.29**, CONSISTENT),
but the **rate at which the null produces a candidate at all more than doubles: P(offender > 0) =
0.200 → 0.453.** The same cell's zero-fraction falls 80.0 % → 54.7 %. R2o and R3o push it the same way
(0.293, 0.280), so it is a channel effect, not a fluke.

**Unmodelled noise makes the null throw up MORE candidates, but not LOUDER ones.** The floor is an
*upper-tail* quantile, and the upper tail is set by the **template**, not by the data: the E-step's
matched-filter cross term is linear in the *model* amplitude, so the loudest fringe fluctuations scale
with the hypothesis `h` (that is the `floor ∝ h^1.66` mechanism of S5.3.1, which this campaign
independently reproduces at 1.67). **The tail is template-dominated; the body is noise-dominated.**

*This is the mechanistic reason the ladder came back CONSISTENT*, and it is a stronger statement than
the null result itself: **the criterion's floor is robust to noise mis-specification BECAUSE it is
loudness-relative.** The property S5.3.2 recorded as a *cost* (the bar rises almost as fast as the
signal, so louder alone buys nothing) is the same property that makes the floor **immune to getting
the noise model wrong.** The one-wall coincidence has a sibling.

**Operationally the criterion is safe**: detection requires `dlnL > max(ln K, floor)`, the floor is
unmoved, so those extra sub-threshold candidates are vetoed exactly as intended. **The margin between
"more candidates" and "more detections" is the floor doing its job** — and it only does that job if the
floor is right, which §4 says it is not at faint `h`.

---

## 6. ANATOMY — no new offenders, and no chromatic suspects

Pooled over all 8 cells (1200 nulls per rung), the pulsars carrying the offender:

| pulsar | R0 | R1 | R2 | R3 | R2o | R3o | ΔR3−R0 |
|---|---|---|---|---|---|---|---|
| J0711−6830 | 163 | 149 | 181 | 170 | 181 | 178 | +7 |
| J1603−7202 | 145 | 141 | 114 | 117 | 111 | 127 | **−28** |
| J0437−4715 | 108 | 128 | 121 | 139 | 130 | 122 | **+31** |
| J1713+0747 | 93 | 77 | 74 | 86 | 88 | 76 | −7 |
| J1909−3744 | 90 | 75 | 81 | 89 | 82 | 78 | −1 |
| J1045−4509 | 68 | 63 | 72 | 54 | 54 | 62 | −14 |
| J1824−2452A | 6 | 6 | 10 | 13 | 8 | 7 | +7 |
| **TOTAL** | **817** | **810** | **837** | **869** | **842** | **833** | **+52 (+6 %)** |

**The offender population does not change.** The same five pulsars — J0711, J1603, J0437, J1713,
J1909 — carry ~75 % of offenders in *every* rung, including the control. Total offender count moves
+6 % from R0 to R3, against a 32 % rise in noise rms.

**The suspects the brief expected did not show up.** **J1824−2452A — the array's reddest pulsar by far
(feather χ²/N = 3243) — contributes 6 offenders in R0 and 13 in R3.** It is not the story. Neither is
B1937+21 (χ²/N = 201), which never appears. **Unmodelled DM and unmodelled red power do NOT route the
false alarms through the chromatically-dirty pulsars.**

What *does* move, modestly, is **J0437−4715: +31 under R3, the largest single shift in the table** —
and J0437 is precisely **the array's smallest trials factor (`ln K` = 1.39, the array minimum, S5.2.3)**.
It is the pulsar whose bar is so low that a noise fluctuation clears it, and it is therefore the pulsar
most exposed when you add noise the model does not know about. **The offender anatomy is set by the
trials factor, not by the noise budget.** That is the "J0437 double edge" (S5.2.4 — *"robustness to
source error and vulnerability to noise are THE SAME PROPERTY"*), showing up in a third independent
place. The counter-moving J1603−7202 (−28) is the same mechanism in reverse and is within Poisson
scatter on 145 counts (±12).

**One paragraph, no rabbit hole, as briefed:** *there is no inflated rung to anatomise.* The excess is
+6 % of offenders, spread across the existing offender set, concentrated on the smallest-K pulsar and
not on the reddest ones.

---

## 7. THE VLBI PRICE ROW

`floor(vlbi) − floor(lit)` at the same `h`, empirical q95, bootstrap errors. **POSITIVE = the VLBI tier's
null floor is HIGHER** — tightening the distance prior *raises* the bar the null must clear.

| rung | h = −13.25 | h = −13.00 | h = −12.75 | h = −12.50 |
|---|---|---|---|---|
| **R0** (control) | **+2.59 ± 1.41** | −1.12 ± 1.22 | +0.64 ± 2.96 | −0.17 ± 5.44 |
| **R2** (the rung the brief names) | **+2.88 ± 1.01** | +0.99 ± 0.89 | −0.65 ± 1.75 | +7.16 ± 6.16 |
| **R3** (the *measured* worst rung, ×1.32 noise) | **+2.76 ± 1.02** | +0.85 ± 1.17 | −1.43 ± 3.29 | −0.16 ± 4.16 |

**THE PRICE, STATED: a VLBI campaign costs ≈ +2.9 ± 1.0 nat of null floor at h = −13.25, and nothing
measurable anywhere else in the box.** At h = −13.0, −12.75 and −12.50 the tier difference is
consistent with zero (the loud cells' errors are 3–6 nat — they cannot resolve a price of this size,
and I will not quote one).

**Two things must travel with that number.**

1. **It is not a realism effect. It is a tier effect.** The control (R0) already carries
   **+2.59 ± 1.41**, statistically indistinguishable from R2's +2.88 and R3's +2.76. **Realism does not
   change the price of VLBI.** The brief asked for the price "on the worst realistic rung"; the answer
   is that the rung does not matter, which is a cleaner result than the one requested. *(The brief
   named R2 as the worst rung; measured, R3 is — 2.63 µs vs 2.21 µs. Both give the same price, so
   nothing turns on it.)*
2. **The mechanism, and the sign is not a paradox.** VLBI shrinks `σ_d` → fewer fringes in the prior
   window → smaller `K_counted` → a **lower trials bar** `ln K` → a pure-noise fluctuation clears
   layer 1 more easily → more and louder offenders → **a higher absolute floor**. This is the J0437
   double edge again, now measured as a *tier-level* quantity. **VLBI buys detections on the signal
   side (ΣK 88 454 → 470, S6.1.3) and pays for them on the null side.** ANCHOR prices only the null
   side. **A VLBI campaign is not free in onset-map units, and this is the first number for what it
   costs.**

---

## 8. HOUSEKEEPING — doc corrections drafted, NOT applied

Per the brief I did **not** edit `MANIFEST.md` or `STORY.md`. Draft text for the next cronus doc
session:

### 8.1 MANIFEST §A / §E — the b20 pickle, stated positively

**REMOVE from §E ("UNKNOWN — can't classify from headers alone")** the entire
`data_products/b20_cw_curn_r0.pkl` bullet. It is not unknown, and it is not "a benign data artifact"
sitting *beside* the canonical data.

**ADD to §A (CANONICAL — load-bearing):**

> - `data_products/*.feather` (116 files, 20.5 MB, git-tracked) — **THE ARRAY, AND IT IS A MOCK.**
>   116 pulsars, 30 225 TOAs, span 22.1549 yr (MJD 52467.38–60559.45). Real pulsar names, sky
>   positions, distance priors and timing-model structure; **simulated TOAs** — telescope `AXIS`, a
>   **single 1440.0 MHz channel** across every TOA, median cadence ≈ 23 d. Supplies TOAs, TOA errors,
>   the timing design matrix, sky positions and `pdist` to `cw_helpers.load_pulsars`,
>   `ignite.build_ignite_problem`, `forge_b1`. **Its `residuals` column is an injected CW + CURN
>   realisation** — bit-identical (max|diff| = 0.0, all 116) to `b20_cw_curn_r0.pkl`. **The
>   certification chain never reads it** (`data = inject_delay(θ_true) + NoiseDrawer.draw(seed)`), so
>   no banked result depends on it — *but any future task that treats these residuals as "the data"
>   is measuring an injected CW.* [ANCHOR 2026-07-12, `ANCHOR_results/anchor_data_forensics.npz`]
> - `data_products/b20_cw_curn_r0.pkl`, `..._w_flags.pkl` — the same array as the feathers, as
>   `enterprise.pulsar_edited.Tempo2Pulsar` objects (needs the cronus fork; `_w_flags` also needs
>   `dill`). **Not a stray artifact — it IS the feathers' residual column**, in enterprise form.

### 8.2 STORY §S6.4 — new PENDING candidate

> [PENDING: **REAL-ARRAY** — port the criterion stack to a real PTA. **Why it is a campaign, not a
> task:** every prior in the criterion is keyed to the 116-pulsar mock — `best_distances.txt`, the
> per-pulsar `ln K_counted`, `ΣK`, the lit/vlbi tiers, the geometry ensembles, the census, and the
> `NoiseDrawer`'s hyperparameters. A real-array anchor is a **re-derivation of the prior stack on a
> different array**, not a substitution of a residual vector.
> **On disk (verified, ANCHOR 2026-07-12):** NG 15 yr `~/projects/NANOGRAV/15yr_dataset/psr_feathers`
> (66 psr, 615 294 TOAs, 16.03 yr, 1705 radio frequencies, loads in `discovery` unmodified);
> NG 20 yr `~/projects/NANOGRAV/20yr_dataset/feather_files` (77 psr, 1 131 412 TOAs, 20.62 yr);
> MPTA DR3 `/data/taylor_group/matt_miles/MPTADR3/partim/fifth_pass` (83 psr, par/tim only).
> **Deliverables:** (i) the prior stack re-derived on the target array; (ii) the floor re-fit against
> that array's own noise; (iii) the programme's first certification numbers that touch a real TOA;
> (iv) **it lifts ANCHOR's ceiling** — a multi-frequency array makes chromaticity *identifiable*, so
> the DM channel can be tested as chromaticity rather than as unmodelled red power (§2.1).]

### 8.3 A convention for the criterion (§4)

The zero-fraction gate on the D2 estimator. Proposed text in §4 above. **This one is not
housekeeping — it is a correction to a live estimator, and it errs permissive.**

---

## 9. BANKED

| artifact | contents |
|---|---|
| `ANCHOR_results/an_{rung}_{cell}_g{sky}_n{seed}.npz` | **7 200 files.** Per null: `dlnL_det` (116), `lnK` (116), `qmax` (116), `mapk`, `ptrue`, `on_true`, `n_true`, `offender`, `resid_rms`, `names`. **Raw per-pulsar `dlnL`/`lnK`/`qmax`, per [[lean-npz-discipline-amended]]** — every floor is re-cuttable from the bank without rerunning a GPU. |
| `anchor_ladder.npz` | the §3 table: per (rung, cell) Gumbel floor + fit error, empirical q95 + bootstrap error, zero-fraction, P(off>0), Mann-Whitney p, verdict, μ, β; the full 150-offender vector for all 48 (rung, cell); per-pulsar offender counts. |
| `anchor_floors.npz` | the Gumbel-only cut (D2 estimator as adopted). |
| `anchor_g1.npz` | the 80-realisation bit-identity replay. |
| `anchor_data_forensics.npz` | the STOP's forensics (per-pulsar χ²/N, telescope, frequency, set hash). |
| `hpc_harbor/anchor/` | `anchor.py` (ladder + floors + g1), `anchor_analysis.py`, `test_ladder.py` (CPU unit test, 20/20 PASS), `g2_check.py`, `probe_*.py`, sbatch. |

**Compute:** 7 200 realisations ≈ 0.5 s each; 4 shards on A100-40 (the shared `a100-80gb` entitlement
was exhausted by concurrent agents). No failures in any shard. **Nothing committed.**

---

## 10. WHAT I WOULD DO NEXT, IN ORDER

1. **Fix the floor estimator before anything else uses it (§4).** Re-cut every faint-`h` cell in the
   SURFACE grid with the empirical quantile + zero-fraction. On this evidence the faint-`h` floors are
   **too low**, which means the uncalibrated cells' counts are **over-stated**, not under-stated. This
   is cheap — the banks exist — and it moves in the direction that *closes* the box further.
2. **Do not spend GPU-hours hardening the floor against noise realism.** The ladder says the floor
   does not care, and §5 says *why* it cannot care: the tail is template-dominated. That is a settled
   question at these amplitudes, and the reason it is settled is a mechanism, not a null result.
3. **The VLBI price (+2.9 ± 1.0 nat at h = −13.25) belongs in the payoff arithmetic**, on the cost side,
   next to the ΣK benefit. It is small, it is real, and it is the first time the null side of a VLBI
   campaign has been priced.
4. **REAL-ARRAY (§8.2)** remains the only way to answer the question the original brief actually asked.

**STOP.**
