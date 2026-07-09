# PENCIL — Crossing time of the registration wall

Agent PENCIL, 2026-07-07. CPU-only pencil calculation (no simulation, no likelihood
calls). Files: `PENCIL_t_crossing.py`, `PENCIL_t_crossing.npz`, `PENCIL_t_crossing.png`.

## Question

The F2 lever-arm ladder's **loosest** pulsar-term registration baseline tolerates a
**1.85e-3 scaled** sky error; blind Earth-term (F-stat) source localisation floors at
**0.05 scaled (6.4°)** on the current 15-yr dataset. The wall closes from *above* as data
accumulate (σ_sky shrinks with SNR). When does the float cross onto the first rung?

## Anchor (measured, from PROJECT_PROGRESS)

| quantity | value | source |
|---|---|---|
| span T₀ | 15 yr | current dataset |
| blind sky float σ₀ | 0.05 scaled = 6.4° | F-stat seeder floor (B0.2 / F1c) |
| per-source SNR₀ | 10.7 | D4 equal-strain, h = 10⁻¹³·⁷⁵ |
| loosest rung | 1.85e-3 scaled | F2 lever-arm ladder (top) |
| L2c pull-in | 1e-4 scaled | conditional-re-solve radius (banked) |

## Assumptions (named)

- **A1** Persistent monochromatic CW, always on — coherent matched filter over the full span.
- **A2** White-noise-dominated band, **fixed cadence, no new pulsars, no noise re-weighting**.
  For a sinusoid in white noise SNR² = h²T/(2Sₙ) is **linear in T** ⇒ SNR ∝ T^½.
- **A3** Blind sky precision is SNR-limited via the geometric array response:
  σ_sky ∝ 1/SNR × (fixed baseline geometry). Geometry frozen under A2 ⇒ only the SNR
  factor moves: **σ_sky(T) = 0.05 · (T/15)^(−½)**.
- **A4** Frequency resolution improves *separately* and *faster*, σ_f ∝ 1/(SNR·T) ∝ T^(−3/2),
  but is **irrelevant**: F2 showed **sky binds** — a 0.05 scaled sky error alone wraps the
  loosest pulsar-term phase by 0.05/1.85e-3 ≈ 27 rad regardless of frequency precision. The
  binding axis therefore scales as T^(−½).

General solution: **T = T₀ · (σ₀/σ_target)^(1/p)**, p = ½ baseline.

## (2) Crossing times — baseline SNR² ~ T (σ_sky ~ T^−½)

Shrink factor = σ₀/σ_target; T/T₀ = (factor)² :

| target rung | shrink factor | T/T₀ | **crossing T** |
|---|---|---|---|
| loosest 1.85e-3 | 27.03 | 730.5 | **≈ 10,960 yr (≈ 11 kyr)** |
| L2c 1e-4 | 500 | 250,000 | **≈ 3,750,000 yr (≈ 3.75 Myr)** |

## (3) vs STRAIN instead of T (fix T = 15 yr)

σ_sky ∝ 1/SNR and SNR ∝ h at fixed span. Reach a rung by loudness instead of waiting:
SNR_thresh = SNR₀·(σ₀/σ_target); log₁₀h = −13.75 + log₁₀(factor).

| target rung | **SNR threshold** | log₁₀h | h |
|---|---|---|---|
| loosest 1.85e-3 | **≈ 289** | −12.32 | 4.8e-13 |
| L2c 1e-4 | ≈ 5,350 | −11.05 | 8.9e-12 |

**First rung needs per-source SNR ≈ 289** (h ≈ 10⁻¹²·³²) — ~27× the current loudest-source
SNR, a strain no realistic single CW in the nHz band provides.

## (4) Sensitivity — SNR² ~ T³ (σ_sky ~ T^−3/2)

**Justify/reject:** for a sinusoid in white noise SNR² is strictly linear in T (A2). T³ needs
the *per-sample* sensitivity to also grow ~T — i.e. red-noise / spin-down-limited spectra
where the narrowing 1/T frequency bin suppresses in-band noise. Stage C freezes noise at
truth in a white-dominated 3–20 nHz band, so **T³ is rejected as the baseline** and reported
only as an optimistic bracket.

| target rung | T = T₀·(σ₀/σ_target)^(2/3) |
|---|---|
| loosest 1.85e-3 | ≈ 135 yr |
| L2c 1e-4 | ≈ 945 yr |

Even this optimistic bound puts the first crossing a century out and full pull-in ~10
centuries — beyond any mission horizon.

## (5) Figure

`PENCIL_t_crossing.png` — σ_sky vs T (log-log), baseline (T^−½, solid) and optimistic
(T^−3/2, dotted) curves, both rungs as horizontal dashed lines, the 15-yr anchor and the
achievable-float band marked, crossing points as markers.

## Scope — what these crossings price (and what they do NOT)

These times price the **blind, self-contained float only**: same 116-pulsar array, noise
and geometry frozen, **no external information**, precision bought purely by integrating
in-band SNR. Explicitly **not priced**:

- **(a) EM counterpart** localizing a loud source (the anchor-paper regime). An external
  host-galaxy-quality sky position supplies the sky precision the blind float cannot reach.
  In this plot's currency a host position is worth **~10⁴ yr of blind integration** — it
  jumps the wall rather than climbing it. Different information channel, not on this axis.
- **(b) New pulsars.** Adding baselines adds *rungs* (more, and possibly looser, registration
  lanes) and lifts the array SNR — this is **lever i (more/looser baselines)**, not lever iii
  (integrate the existing float). Changes the ladder, not priced here.
- **(c) Inference-paradigm routes (deliverable R).** Marginalizing the fringe/source posterior
  *integrates over* the wall instead of localizing across it — a different estimator class, not
  a longer blind integration. Out of scope for this pencil. **RESULT (R, 2026-07-08):** it does NOT
  rescue certification — the fringe-marginalised posterior smears (f = Z_needle/(Z_needle+Z_plateau)
  = 6.9e-7; break-even sky prior 0.188°/loud source). That 0.188° is the INFERENCE-side counterpart
  of this pencil's point-estimation wall: both say the 15-yr array cannot self-localise the sky
  finely enough — by optimisation (this pencil) OR by marginalisation (R).

This calculation bounds lever iii (blind integration of the current float) only.

## Consistency — three independent derivations agree

Internal cross-check, all landing on **27×** at the 0.05 scaled ceiling:
- **Strain axis (step 3):** SNR_thresh / SNR₀ = 289 / 10.7 ≈ **27**.
- **Time axis (step 2):** √(T_cross/T₀) = √(10,960/15) ≈ **27**.
- **F2 wall (geometry):** blind float 0.05 / loosest rung 1.85e-3 ≈ **27**.

σ_sky ∝ 1/SNR ∝ T^−½ ties them: the same 27× gap appears as a strain factor, as √(time ratio),
and as the raw geometric ladder wall. Three routes, one number.

## Bottom line

The registration wall does close from above, but on the baseline white-noise scaling the
float reaches even the **loosest** rung only at **T ≈ 11 kyr** and the L2c pull-in at
**T ≈ 3.75 Myr**. Buying it with strain instead needs **per-source SNR ≈ 289** (h ≈ 10⁻¹²·³).
The optimistic T³ bound (rejected physically) only softens 11 kyr → 135 yr. This is the same
≥1.4-dex separation F2/L2c found, now expressed as a timescale: the wall does not close on
any observationally meaningful horizon. Consistent with the standing **NO-GO** — do not start B1.
