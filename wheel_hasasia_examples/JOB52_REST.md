# JOB 52 — Parts 3 to 8

HEAD at the time of writing: `74cdf93` (branch `MM_playground`, in sync with `origin`). `wheel_hasasia_examples/` **is committed**. All Job 52 files are uncommitted, written 2026-11-11/12.

## Five-line summary

1. **A correct distance cannot be told from a wrong one when only the detectable sources are known.** Census priors, N = 7 or 12, reference configuration: no catalogue reaches 0.90 purity (best 0.74 and 0.83). With every candidate known a catalogue exists: 12 pulsars per sky at 0.87 purity.
2. **The bins behave exactly as the exact (Rice) phase-error law says** (8 of 8 bins within 2 standard errors). The typical known pulsar term holds 5.6% of its bin after the known set is subtracted; without retardation it would hold 20%.
3. **What sets the ceiling is the discrete sources outside the known set.** Removing them lifts N = 12 from 5% to 33%; removing the sea lifts it to 6%; removing timing noise or changing the amplitude does nothing. Retardation costs a factor 2.7 at N = 12.
4. **Information per bin does predict decoding, but arms do not fall on one curve.** Pooled, 50% is reached at C − ln K = −1.8 nats (finest tooth) and 90% at +0.9; individual arms sit up to 2 nats away, and Arm D is offset by +1.7.
5. **Part 7 stopped: the two codes cannot share a sky.** The machine draws the angle between each pulsar and each source independently, so its skies correspond to no real arrangement of pulsars and sources. **Part 8:** the Fisher code rerun today (36 minutes) reproduces the June file to 6e-7 or better on every curve.

## Where the prompt and the code disagree

- **Part 7's premise fails.** `job50_common.py:181`: `mu_all = np.arccos(rx.uniform(-0.9, 0.9, len(f2)))`, drawn per pulsar from a pulsar-specific stream (line 178: `rx = np.random.default_rng([50, seed0, j, i])`). Test: for real positions the matrix of cos μ (pulsars × sources) has three non-zero singular values. For sky 100 (40 pulsars × 2,000 sources) the first six are 1.00, 0.99, 0.98, 0.97, 0.96, 0.95 and the first three carry 9% of the variance.
- **"Amplitude 3.7 and 10.6 × 10⁻¹⁵" cannot change a joint fit with a fixed candidate list**, except through timing noise: every gravitational-wave term scales together. The rows are run that way (timing noise rescaled by the inverse amplitude ratio) and show no change. The amplitude dependence in earlier jobs came from the "counts only" scaling and from the candidate list.
- **Prediction P11's numbers are off by about 2 nats** and the arms do not collapse onto one curve (Part 6).

## Provenance

| Numbers | Produced by | Raw output |
|---|---|---|
| Parts 3, 6 | margins and C saved by `job52_common.Joint52.psr` in the Part 1/2 runs (`job52_data/run_main_*.pkl`, `run_mainx_*.pkl`); tables by `job52_rest.py` | `job52_logs/analysis_rest.log`, `job52_data/analysis_rest.json` |
| Parts 4, 5 | `job52_run.py` runs `p5_c` (retardation on) and `p5off_c` (off): `Joint52(mode="p5", part4=True)` | `job52_data/run_p5_c.pkl`, `run_p5off_c.pkl`; same log |
| Part 7 test | inline script recorded in this report's text | — |
| Part 8 | copy of `CW_transition/prong2_transition.py` with only its output path changed, run in env `sculpt` | scratch directory (not in the repository) |
| Figures | `job52_fig2_phase_error.png`, `job52_fig3_information.png` | `job52_rest.py` |

All code uncommitted, modified 2026-11-11/12. Reference configuration as in `JOB52_MAIN.md`. Census priors: skies 100–119 (1,580 pulsars). Harness: skies 100–199.

---

## Part 3 — Can a correct distance be told from a wrong one?

**Score:** the joint-fit margin = (statistic plus prior at the best grid point) − (largest value of statistic plus prior at any grid point more than 0.25 pc from the best). **Gate:** train on even-numbered skies, apply blind to odd-numbered; purity of the top g% by score for g = 0.5% to 100% in 0.5% steps; 7-point mean with the window shrinking at the edges (no zero padding); the last g at which smoothed train purity is at least 0.90. No fallback. Timing ×1.

**Census priors (10 train skies, 10 held-out):**

| Configuration | N known | Result on held-out skies | Catalogue K < 1 / K < 3 / K ≥ 10 / K ≥ 60 | Right per sky before gating |
|---|---|---|---|---|
| reference | 7 | **no 0.90 catalogue** (best smoothed train purity 0.740) | — | 3.3 |
| reference | 12 | **no 0.90 catalogue** (best 0.828) | — | 4.2 |
| reference | every candidate | 12.0 [9.6, 15.3] per sky at purity 0.867 [0.765, 0.950] | 0.04 / 0.12 / 0.80 / 0.48 | 21.9 |
| Arm D | 7 | 1.0 [0.4, 1.6] per sky at purity 1.000 (10 pulsars in all) | 0.30 / 0.90 / 0.00 / 0.00 | 5.7 |
| Arm D | 12 | 2.1 [1.2, 2.9] per sky at purity 1.000 (21 pulsars in all) | 0.19 / 0.57 / 0.38 / 0.14 | 7.7 |
| Arm D | every candidate | 31.0 [26.8, 35.9] per sky at purity 0.897 [0.857, 0.932] | 0.02 / 0.05 / 0.91 / 0.65 | 40.6 |

**Harness priors (50 / 50 skies), for comparison:**

| Configuration | N known | Result on held-out skies | Catalogue K < 1 / K < 3 / K ≥ 10 / K ≥ 60 |
|---|---|---|---|
| reference | 7 | 26.3 [23.9, 28.5] per sky at 0.909 | 0.40 / 0.97 / 0.00 / 0.00 |
| reference | 12 | 33.1 [31.3, 34.7] per sky at 0.913 | 0.41 / 0.96 / 0.01 / 0.00 |
| reference | every candidate | 68.9 [68.0, 69.6] per sky at 0.901 | 0.35 / 0.79 / 0.18 / 0.11 |
| Arm D | 7 | 34.2 [32.1, 36.2] per sky at 0.891 | 0.41 / 0.97 / 0.00 / 0.00 |
| Arm D | 12 | 42.5 [41.0, 44.0] per sky at 0.906 | 0.42 / 0.96 / 0.01 / 0.01 |
| Arm D | every candidate | 77.8 [77.5, 78.1] per sky at 0.897 | 0.31 / 0.70 / 0.27 / 0.19 |

Reading: with the detectable sources known and real priors, nothing can be certified in the reference configuration, and 1–2 pulsars per sky in Arm D, mostly ones the prior nearly decides. The harness catalogues at N = 7 and 12 are 96–97% pulsars with K < 3: the prior does the work. Only when every candidate is granted does a catalogue of wide-prior pulsars appear. Census catalogues rest on 10 held-out skies.

---

## Part 4 — The jitter itself, measured in the bins

N = 12 known, census priors, reference configuration; 18,960 (pulsar, source) pairs (12 per pulsar).

A pulsar term spreads over several bins, so each quantity uses its leakage profile as weights:

- measured phase = phase of Σ over bins of (datum − known Earth terms − other known pulsar terms at their true distances) × conj(the pair's true phasor in that bin);
- β_unmod = own power ÷ (own power + unmodelled power weighted by the pair's leakage profile);
- β_bin = own power ÷ (all gravitational-wave power under the same weights, sea included).

**Distribution over pairs (10th / 50th / 90th percentile):**

| | β_bin | β_unmod (×1) | β_unmod (timing ÷ 10⁶) |
|---|---|---|---|
| retardation on | 0.007 / 0.046 / 0.211 | 0.009 / 0.056 / 0.229 | 0.009 / 0.057 / 0.232 |
| retardation off | 0.051 / 0.147 / 0.418 | 0.072 / 0.202 / 0.473 | 0.073 / 0.206 / 0.474 |

**Phase error in bins of β_unmod, retardation on, timing ×1:**

| β_unmod | Pairs | Measured rms δφ (rad) | Exact rms | Small-angle √((1−β)/2β) | Measured ⟨cos δφ⟩ ± s.e. | Exact ⟨cos δφ⟩ |
|---|---|---|---|---|---|---|
| 0.003–0.006 | 834 | 1.716 | 1.747 | 10.28 | 0.094 ± 0.025 | 0.060 |
| 0.006–0.013 | 1,468 | 1.688 | 1.717 | 7.20 | 0.109 ± 0.018 | 0.086 |
| 0.013–0.026 | 2,510 | 1.678 | 1.676 | 5.07 | 0.126 ± 0.014 | 0.123 |
| 0.026–0.054 | 3,906 | 1.601 | 1.614 | 3.51 | 0.187 ± 0.011 | 0.177 |
| 0.054–0.113 | 4,677 | 1.529 | 1.523 | 2.43 | 0.245 ± 0.010 | 0.254 |
| 0.113–0.232 | 3,150 | 1.396 | 1.387 | 1.66 | 0.354 ± 0.011 | 0.366 |
| 0.232–0.480 | 1,825 | 1.100 | 1.110 | 0.99 | 0.570 ± 0.013 | 0.565 |
| 0.480–0.990 | 27 | 0.728 | 0.882 | 0.72 | 0.799 ± 0.071 | 0.712 |

3.0% of pairs fall below 0.003; none above 0.99. Exact values: phasor 1 + n with complex Gaussian n, E|n|² = (1 − β)/β; ⟨cos δφ⟩ in closed form (Rice), rms by 200,000 draws at the bin's median β (⟨cos⟩ averaged over the β values in the bin).

**Same table with timing noise ÷ 10⁶:** measured ⟨cos δφ⟩ 0.095, 0.109, 0.125, 0.188, 0.249, 0.350, 0.569, 0.789 (pairs 833, 1,466, 2,505, 3,895, 4,642, 3,166, 1,852, 38); measured rms 1.716, 1.687, 1.679, 1.600, 1.524, 1.402, 1.101, 0.762. Unchanged within errors: timing noise is not what jitters the phase.

**Figure 2:** `job52_fig2_phase_error.png`.

**One sentence:** the machine's actual bins follow the exact relation in every bin (largest departure 1.4 standard errors), the small-angle formula is wrong by a factor 1.2 to 6 everywhere below β = 0.23, and the typical pair sits at β_unmod = 0.056, where the rms phase error is 1.5 rad and ⟨cos δφ⟩ is 0.25.

---

## Part 5 — What sets the ceiling

Census priors, skies 100–119, ×1, reference configuration, one change at a time (all pulsars / K ≥ 60 / K ≥ 600):

| Change | N = 12 | Every candidate |
|---|---|---|
| reference | 0.053 / 0.021 / 0.000 | 0.271 / 0.188 / 0.080 |
| retardation off | 0.141 / 0.104 / 0.064 | 0.535 / 0.484 / 0.325 |
| sea removed | 0.061 / 0.024 / 0.000 | 0.472 / 0.404 / 0.250 |
| discrete sources outside the known set removed | 0.330 / 0.246 / 0.084 | 0.766 / 0.725 / 0.498 |
| timing noise removed | 0.054 / 0.022 / 0.000 | 0.275 / 0.193 / 0.082 |
| amplitude 3.7 × 10⁻¹⁵ | 0.054 / 0.021 / 0.000 | 0.270 / 0.188 / 0.075 |
| amplitude 10.6 × 10⁻¹⁵ | 0.053 / 0.022 / 0.000 | 0.275 / 0.190 / 0.080 |
| known Earth terms left in | 0.042 / 0.015 / 0.000 | 0.162 / 0.097 / 0.022 |

Intervals are in `job52_logs/analysis_rest.log`; all-pulsar half-widths are 0.01 at N = 12 and 0.03–0.04 at every candidate.

- **Retardation off:** 2.7 times more distances at N = 12 and 2.0 times at every candidate; retardation is the second-largest loss.
- **Sea removed:** nothing at N = 12 (0.053 → 0.061), a large gain at every candidate (0.27 → 0.47); the sea matters only once the discrete sources are known.
- **Discrete sources outside the known set removed:** the largest gain in both columns (0.05 → 0.33; 0.27 → 0.77); these sources are the ceiling.
- **Timing noise removed:** no change.
- **Amplitude 3.7 and 10.6:** no change (see top: with a fixed source list the amplitude only rescales timing noise).
- **Known Earth terms left in:** a loss of a fifth at N = 12 and of 40% at every candidate; subtracting the known Earth terms is worth having but is not the main effect.

Even with the sea alone as noise (discrete sources removed), half of the K ≥ 600 pulsars still fail at every candidate.

---

## Part 6 — Information per bin

C = Σ over bins of ln(1 + P_known / P_unmodelled), per pulsar, per arm (104 arms: 52 per prior set; 492,960 decodes).

**Distribution of C, census priors (10th / 50th / 90th percentile, nats):**

| Arm | N = 7 | N = 12 | N = 30 | Every candidate |
|---|---|---|---|---|
| reference ×1 | 0.36 / 1.00 / 1.73 | 0.57 / 1.26 / 2.16 | 1.00 / 1.82 / 2.95 | 2.85 / 4.05 / 5.43 |
| reference, zero noise | 0.37 / 1.01 / 1.78 | 0.58 / 1.28 / 2.20 | 1.01 / 1.85 / 3.01 | 2.89 / 4.09 / 5.50 |
| A ×1 | 0.77 / 1.91 / 3.24 | 1.23 / 2.43 / 4.00 | 2.02 / 3.52 / 5.34 | 5.14 / 7.06 / 9.26 |
| B ×1 | 0.36 / 0.99 / 1.78 | 0.61 / 1.29 / 2.24 | 1.04 / 1.88 / 3.03 | 2.97 / 4.17 / 5.63 |
| C ×1 | as reference | as reference | as reference | as reference |

Other N (reference ×1): N = 1: 0.01 / 0.08 / 0.37; N = 100: 1.61 / 2.66 / 3.96. All arms, both prior sets: `job52_logs/analysis_rest.log`.

**Fraction decoded against C − ln K, all arms pooled:**

| K taken from | 50% at | 90% at | Logistic width |
|---|---|---|---|
| finest tooth of the known set | −1.81 nats | +0.85 nats | 1.21 |
| coarsest tooth of the known set | +0.67 nats | +4.14 nats | 1.58 |

Binned, finest tooth: 0.02 at −7.5 to −6; 0.14 at −4.5 to −3; 0.43 at −3 to −1.5; 0.71 at −1.5 to 0; 0.84 at 0 to 1.5; then a plateau at 0.89–0.90 from 1.5 to 7.5 nats, and 0.94–0.98 above 7.5. The plateau means a tenth of pulsars with apparently ample information still fail.

**Do the arms fall on one curve? No.** Fitting each arm's 50% point with the pooled width (finest tooth):

- spread across arms: minimum −2.37, median −1.61, maximum +5.63 nats;
- by prior set: census −1.33, harness −2.04;
- by arm: reference −1.70, A −0.92, B −1.62, C −1.94, D +0.04;
- by known-set size: N ≤ 12 −1.60, N = 20–300 −1.58, every candidate −1.63;
- by timing: ×1 −1.62, zero noise −1.38.

Known-set size and timing collapse well (within 0.3 nats). Prior set shifts the curve by 0.7 nats, Arm A by 0.8, Arm D by 1.7. With K from the coarsest tooth the collapse is worse (census 2.03 against harness 0.29; N ≤ 12 at 1.08 against every candidate at 2.81). C is a useful ordering variable within one configuration but not a universal one.

**Figure 3:** `job52_fig3_information.png` (one line per arm, both definitions of K, pooled logistic dashed).

---

## Part 7 — Independent check in the real-likelihood code: stopped

The two codes cannot share a sky without major work.

- **Why.** The machine has no source positions and no pulsar positions. It draws the angle between each pulsar and each source independently (quoted at the top; singular-value test at the top). A real likelihood needs positions, and no set of positions reproduces the machine's angles.
- **What the other side offers.** `CW_transition/prong2_transition.py` has a full time-domain residual model with positions (`single_source_residual(toas, p_hat, L_kpc, src)`, line 81; `total_residual`, line 105), so a grid likelihood in one pulsar's distance with all source parameters fixed is straightforward there.
- **What is needed.**
  1. A geometry mode in `collect50`: give every source and pulsar a sky position and derive μ from them, replacing the independent draws. This changes every stream downstream of line 163, so Jobs 49–52 would not reproduce in that mode.
  2. Per-source inclination, polarisation and phase, and real antenna patterns in place of the machine's fixed factor `np.sqrt(0.2*120)`.
  3. Injection of about 18,700 sources (Earth and pulsar terms, chirping) plus the sea as a Gaussian process into time-domain residuals at 14-day cadence over 45 years, and a check that the machine's bin powers match the injected ones.
- **Estimate:** one to two days, most of it item 3 and re-validating the machine in geometry mode.

---

## Part 8 — Fisher figures

**Diff.** `prong2_results.npz` was written 2026-06-25 20:06 and committed in `f0e7827` (20:08). Between that commit and HEAD, `prong2_transition.py` changed by 576 added lines and 1 removed:

- one 7-line block defining a constant, `KNEE_COEF_PER_PSR = 52.75 / 116.0`, with its comment (after `confusion_onset`, line 348);
- 569 lines of new functions (coherence experiment and N × t_c map: `coh_single_residual` … `run_ntc_map`), inserted before the main block;
- the main block gained two early branches (`coherence`, `ntcmap`) that are skipped when the script is run without arguments.

No line of any function used by the main sweep was modified.

**Rerun.** First attempt was stopped at the 30-minute cap (4 of 5 sweeps done); Matt then lifted the cap and the rerun was repeated to completion on 2026-11-12. Cost: 2,179 s (36 minutes) on the GPU, env `sculpt`: fixed-per-source sweep 392 s, fixed-total 312 s, population 312 s, distance-error sweep 624 s, knee diagnostic 538 s.

Comparison with the June file `CW_transition/prong2_results.npz` (commit f0e7827), key by key: same 103 numeric arrays in both, none missing, none extra, no shape or label difference, no NaN mismatch. 66 arrays are bit-identical.

| curve | max relative difference |
|---|---|
| fixed per-source: conditional, marginal, 16–84% bands | ≤ 1.2e-10 |
| fixed total: same | ≤ 4.0e-11 |
| population: same | ≤ 1.7e-10 |
| marginal curves at the three pseudo-inverse thresholds | ≤ 5.9e-11 |
| distance-error curves (`j_*`) | ≤ 2.2e-11 |
| knee diagnostic ratios, five settings | ≤ 5.8e-7 (T = 15 yr, 3–12 nHz); others ≤ 4.3e-9 |
| knees, N*, signal-to-noise at the knee | identical |

One scalar differs by a large ratio and means nothing: `mode_indep_maxdev`, the largest difference between two curves that should be equal, is 2.8e-12 in June and 5.5e-12 now; both are rounding noise. **The Fisher code as it stands reproduces the June file to 6e-7 or better on every curve**, so the 576 lines added since June did not change the main sweep. The June file and figures are untouched (file date still 25 June; the rerun wrote only to the session scratch directory). Source: scratch copy of `CW_transition/prong2_transition.py` (uncommitted changes since f0e7827; output path changed, nothing else).

## Verdicts on the predictions

| | Prediction | Result | Verdict |
|---|---|---|---|
| P8 | at N = 7 and 12 the catalogue, if one exists, is more than 90% K < 3 | census reference: no catalogue at either N. Census Arm D: 90% (N = 7, 10 pulsars) and 57% (N = 12, 21 pulsars). Harness: 96–97% in all four cases | **Holds** where the catalogue is of useful size (harness) and trivially in the census reference; **fails** for census Arm D at N = 12 on 21 pulsars |
| P9 | measured ⟨cos δφ⟩ follows the exact curve within sampling error | 8 of 8 bins within 1.4 standard errors, at ×1 and at zero noise | **Holds** |
| P10 | retardation off raises the N = 12 all-pulsar fraction by at least 1.5× | 0.053 → 0.141, a factor 2.7 | **Holds** |
| P11 | finest tooth: all arms on one curve of C − ln K, 50% near 0, 90% near +3 | pooled 50% at −1.8, 90% at +0.9; arms spread over several nats; Arm D offset by +1.7 | **Wrong** in both the numbers and the collapse |

## Limits

- Census catalogues in Part 3 rest on 10 held-out skies; the Arm D catalogues at N = 7 and 12 hold 10 and 21 pulsars.
- Part 4's weighted definitions of β reduce to the single-bin ones only when a pulsar term sits in one bin.
- Part 5's "amplitude" rows keep the ×1 candidate list.
- Arm D inherits the caveats of `JOB52_MAIN.md` (perfect Earth-term removal; lifted candidate cap).

---

## Addendum (2026-11-12): two questions on the figures

**1. Why β_unmod stops near 0.5 (Figure 2).** It is the machine's population, mainly a built-in construction, with the upper-limit rejection as a second, weaker cap.

- **The machine's skies are not the raw draw.** `job50_common.py:115`: `sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)`. `_pop18` (notebook cell 18:11-15) calls `build_pop15(_f2, _h2, _sea50v[j], s_to_s30(s15), 3, seed=_sd)`, and `build_pop15` (cell 14:4179-4186) takes the three bins in 2–10 nHz with the loudest sources and sets `Ptop = max(Pneed, s_tgt*Pel0)`, rescaling everything else in the bin to `(1 - s_tgt)*Pel`. With `s_to_s30(0.35)` = 0.5185, three sources per sky are forced to hold 52% of their 30-year bin (strain raised by ×1.0 to ×4.0, median ×1.6).
- **Predicted edge.** At 45 years that share becomes 0.618 of the Earth-term power; a pulsar term then has its own power against the unknown Earth and pulsar power, so β_unmod = 1/(1 + 2 × 0.382/0.618) = 0.447 before the other known sources are subtracted. Measured, retardation off, one value per constructed source (skies 104, 109, 118): 0.33 to 0.51, median 0.476.
- **Test** (`job52_diag_beta.py`, skies 104, 109, 118, 2,844 pairs): of the 210 pairs with β_unmod in 0.40–0.48, 182 belong to the constructed sources (711 pairs) and 28 to natural ones (2,133 pairs). The constructed sources rank 1st to 5th in the known set. What limits every pair is the unknown discrete sources in its bins (at least 0.48 of its own power in every pair); the sea is as low as 0.025, timing 0.1 or less.
- **The few pairs above 0.48 are natural sources at 18–20 nHz** (9 of 13 in these skies; the largest, 0.62, is sky 109's loudest at 19.7 nHz). The upper limit bounds these: the raw draw rejected 30 sources in 20 skies (0 to 4 per sky), the loudest surviving source sits at a median 0.87 of the limit, and 13 of 20 skies have one within 20% of it (`job52_diag_pop.py` → `job52_logs/diag_pop.log`). After the construction no source exceeds the limit (largest 0.998).
- **Consequence.** The Part 4 histogram's upper edge, and the Part 1 share-of-power column for N ≤ 3, partly reflect a designed dominance of 0.35 at 15 years, not the mass function alone. This applies to every Job 49–52 result from this machine.

**2. The lines offset to the right in Figure 3 are the Arm D zero-noise runs.** Per-arm 50% point, finest tooth (pooled −1.8): census D every candidate, zero noise +5.6; harness the same +5.1; census D N = 30 / 12 / 7 zero noise +3.0 / +2.3 / +1.8; harness +1.4 / +1.0 / +0.8. The next arm is the reference at N = 1 (+0.3); all others lie between −2.4 and −0.4. Cause as in `JOB52_MAIN.md`: with bins to 100 nHz and no timing noise, C counts about 100 correlated bins as independent (median C 12.1 against 7.3 at ×1 for census D every candidate) while the decode fraction falls. In the coarsest-tooth panel the same eight arms lead (+3.1 to +10.0), followed by census Arm A and Arm B every candidate (+3.6 to +3.8).
