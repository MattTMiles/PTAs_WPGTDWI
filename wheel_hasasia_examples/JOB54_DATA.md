# JOB 54 — Part 1 (data for the talk), Part 3 (per-bin model), Part 4 (the pile-up split)

HEAD: `74cdf93` (branch `MM_playground`). All Job 52–54 files are uncommitted (staged), written 2026-11-11/12. Part 2 is in `JOB54_LOCALISATION.md`.

## Five-line summary

1. **`talk_data/` is ready to copy: 22 files, 7.0 MB**, with `README.md` and `check.py` (which passes). Every sky-based file uses natural skies.
2. **A defect in the simulated skies was found and corrected.** Unknown sources' pulsar terms had no random phase, so for tight priors they lined up with the true distance and helped the fit. Census-prior results drop a little (N = 12: 0.046 → 0.041; every candidate 0.260 → 0.238). **Job 53's prior-width slice was too optimistic** (N = 12 at 1 pc: 0.72 → 0.50; at 3 pc: 0.39 → 0.29). All talk files use the corrected machine.
3. **With that corrected, the one-source-per-bin model fed the real shares reproduces the skies within 0.007 at every prior width** (Part 3). The paper can say the per-bin model explains the sky result.
4. **The ladder (N = 12 known, all pulsars):** quiet sky 0.999 → add the sea 0.744 → add every other binary 0.071 → add retardation 0.041 → better timing 0.042 → 25-year span 0.023 → 15-year span 0.015. Out of 80 pulsars: 80, 60, 6, 3, 3, 2, 1.
5. **The loss from retardation is mostly "fainter" (Part 4):** giving the known pulsar terms their retarded amplitude alone accounts for 80–100% of the loss; moving them alone, 45–55%. The destination bins are not made louder by the migration (pulsar-term power below 5 nHz: 0.97 of its unretarded value). Pulsar terms are 47% of what a known pulsar term competes with, Earth terms 52%.

## Where the prompt and the code disagree

- **The phase defect (summary point 2).** `job53_common.py:119` (same in `job51_common.py:68`, `job52_common.py:110`): `eP = np.exp(1j*kq*Lt)[:, None]*D["cP"]`. Every candidate's pulsar term is given the phase κ × (true offset from the prior centre) and nothing else. For a known source that is correct (its fixed phase is part of the granted template). For the roughly 1,290 candidates *outside* the known set it means their phases are all κ_j × Lt: random when Lt is tens of parsecs, nearly equal and pointing at the truth when Lt is a parsec or less. Non-candidates and all Earth terms did have random phases. Test: Part 3 below (the gap of 0.10–0.22 at 0.3–3 pc closes to ±0.007 when each candidate is given its own random phase, `job54_common.phasors`, `B["phi"]`).
  - **Affected:** Job 53 Part 3b and 3c (prior-width slice and certification at 1 and 3 pc); the harness-prior results of Jobs 51 and 52 (those priors are about 1 pc for 70% of pulsars); census results by 0.005 to 0.02; the Job 53 span results (0.022 → 0.015 at 15 yr).
  - **Not affected:** the heatmap (no skies), the localisation study, anything from Jobs 49–50 (different data model).
  - The brief's "with 12 binaries known, 90% are right at 0.4 pc, 72% at 1 pc, 39% at 3 pc" becomes **90% at 0.16 pc, 50% at 1 pc, 29% at 3 pc**.
- **The "detectable" known set misses most detectable sources.** The candidate list requires single-pulsar signal-to-noise 3 (`job50_common.py`: `cand = np.where((arr_all >= 3.0) & (f2 >= 2e-9) & (f2 <= fmax_cand) & (rp_all >= 3.0))[0]`). On natural skies that keeps only low-frequency sources (median 2.75 nHz, none above 14 nHz on sky 100). Admitting the 300 in-band sources loudest against the background (`collect50(cand_bg=300)`, natural skies 100–139): sources detectable against the background rise from 4.6 to **12.6** per sky (`arr_bg` ≥ 3) and from 10.1 to **30.5** (`arr_bg_opt` ≥ 3); 7.8 of the 12 top-ranked are then above 15 nHz. **The fraction right does not change for N = 7 or 12** (0.026 → 0.028; 0.045 → 0.041), and rises at N = 30 (0.078 → 0.103). So the result stands, but the grey "detectable" band on the curve belongs at 13 to 30 sources, where the curve reads 0.04 to 0.07 (reference). (Run with the uncorrected phases; census priors, so unaffected to 0.005.) Source: `job54_bgknown.py` → `job54_logs/bgknown.log`.
- **1A selection.** No pulsar of the first natural sky with a 1–3 pc prior is decoded correctly at N = 12 once the phases are corrected; the "right" example is from the second natural sky (101), as the rule allows ("first sky that has one").
- **1A, the clean N = 3 statistic of the "wrong" pulsar peaks 21.9 pc from the truth.** With only three combs and a 25 pc prior the comb pattern itself repeats; three clean combs do not single out the truth. The N = 12 clean statistic does (−0.006 pc).
- **1B, the chosen bin is not typical.** For that pulsar the loudest known source holds 53% of what it competes with (β_unmod 0.535); the median over all pairs is 0.055. The rule picked it; say so in the talk or use `spectrum_example.json` to show a typical one.
- **1E, "do not skip skies":** done (`collect50(noskip=True)`). The skipped skies are poorer at 15 years (0.010 against 0.017), as Job 53 inferred.
- **1G:** the Fisher curves' information is in pc⁻² (the plotting code converts with 10³/√I to milliparsec).

---

## Part 1 — `talk_data/`

Size 7.0 MB. `python talk_data/check.py` prints every file's shape and one sanity number (`job54_logs/check.log`). Full column lists and units: `talk_data/README.md`.

| File | Exported or rerun | Sanity number |
|---|---|---|
| `combs_example.json` | rerun, one pulsar each (`job54_export.py`) | "wrong" pulsar (sky 100, pulsar 103, prior 25.0 pc, K 270): real N = 12 winner +43.3 pc from truth; clean winner −0.006 pc. "Right" pulsar (sky 101, pulsar 97, prior 1.00 pc, K 9): −0.08 pc. Direct evaluation equals the fit's FFT statistic to 10⁻⁶ |
| `phasor_example.json` | rerun | terms sum to the bin datum: residual 2 × 10⁻¹⁶; minus the known Earth terms equals the datum the fit uses: 2 × 10⁻¹⁵ |
| `spectrum_example.json` | rerun | parts sum to the total in every bin (difference 0) |
| `curve.csv` | **rerun** (corrected machine, 200 natural populations) | reference, every candidate 0.238 (Job 53: 0.260); N = 12 0.041 (0.046) |
| `prior_slice.csv` | **rerun** (corrected; 100 natural populations, skies 100–199) | N = 12 at 1 pc: 0.503 |
| `prior_slice_job53.csv` | exported from Job 53 (kept for the record; superseded) | 0.72 at 1 pc |
| `heatmap_Ns{2,4,8,12}.csv`, `heatmap_contours.csv` | exported from `job53_data/heat_grid.npz` | 400 cells each; contours at 0.5 and 0.9 |
| `jitter_points.csv`, `beta_hist.csv` | rerun (Job 53 did not save phases) | 189,600 pairs; measured ⟨cos δφ⟩ within 2.8 standard errors of the exact value in all 18 bins |
| `ladder.csv` | rerun | see below |
| `answers_example.csv` | rerun | sky 100: 3 of 79 right; 20 marked (every fourth row) |
| `fisher_curves.csv` | exported from `CW_transition/prong2_results.npz` and `PAPER_figures/job48b_turnover_per_bin.csv` | 672 rows |
| `localisation.csv`, `skymaps.npz` | Part 2 (`job54_loc_*.py`) | 1,830 rows; 11 maps of 256 × 256 plus 5 zoomed |
| `sensitivity_vs_nph.csv`, `sources_strain.csv` | Part 2 (`job54_sens.py`) | 47,520 rows; 3,000 binaries |

Tests behind the reruns: `job54_common.Joint54` returns Job 53's flags for 79 of 79 pulsars in 16 arms on natural sky 100 (`job54_logs/test.log`), before the phase correction is switched on.

### The ladder (`ladder.csv`; natural skies, census priors, corrected machine, 200 populations, no sky skipped)

Fraction right, all pulsars [95% interval] / K ≥ 60 (booked list) / count out of 80:

| Step | N = 12 known | N = 7 known |
|---|---|---|
| 1. known binaries + timing noise only | 0.999 [0.999, 1.000] / 0.999 / 79.9 | 0.954 [0.947, 0.960] / 0.947 / 76.3 |
| 2. add the sea | 0.744 [0.722, 0.767] / 0.712 / 59.5 | 0.486 [0.465, 0.509] / 0.428 / 38.9 |
| 3. add every other binary (no retardation anywhere) | 0.071 [0.066, 0.077] / 0.040 / 5.7 | 0.044 [0.040, 0.048] / 0.020 / 3.5 |
| 4. add retardation (the reference) | **0.041 [0.037, 0.044] / 0.014 / 3.2** | **0.028 [0.025, 0.030] / 0.008 / 2.2** |
| 5. step 4, timing ×10 better | 0.042 / 0.015 / 3.4 | 0.027 / 0.007 / 2.2 |
| 6. step 4, timing noise ÷ 10⁶ | 0.042 / 0.015 / 3.4 | 0.027 / 0.007 / 2.2 |
| 7. step 4, 25-year span (all 200 skies) | 0.023 [0.021, 0.026] / 0.008 / 1.9 | 0.018 [0.015, 0.020] / 0.006 / 1.4 |
| 7. step 4, 15-year span (all 200 skies) | 0.015 [0.013, 0.017] / 0.005 / 1.2 | 0.014 [0.012, 0.016] / 0.003 / 1.1 |

- **Steps 2 and 3 swapped:** known binaries → add every other binary (no sea): 0.081 [0.075, 0.088] (N = 12), 0.047 (N = 7) → add the sea: 0.071, 0.044. The other binaries do nearly all the damage in either order; the sea alone costs a quarter (N = 12) to a half (N = 7).
- **Span, by sky group, N = 12:** 15 yr: 0.017 on the 157 skies Job 53 used, 0.010 on the 43 it skipped. 25 yr: 0.023 on 197, 0.013 on 3.
- Prior only: 0.006. The file also gives K ≥ 60 by the known set's coarsest tooth, which is defined on every sky.

### The curve (`curve.csv`; corrected machine)

| N known | Reference: all / K ≥ 60 / K ≥ 600 | Right per sky (of 79) | Earth terms removed: all / K ≥ 60 / K ≥ 600 |
|---|---|---|---|
| 1 | 0.009 / 0.001 / 0.000 | 0.7 | 0.010 / 0.001 / 0.001 |
| 3 | 0.015 / 0.004 / 0.000 | 1.2 | 0.020 / 0.003 / 0.001 |
| 7 | 0.028 [0.025, 0.030] / 0.008 / 0.001 | 2.2 | 0.047 / 0.018 / 0.005 |
| 12 | 0.041 [0.037, 0.044] / 0.014 / 0.002 | 3.2 | 0.079 / 0.041 / 0.010 |
| 30 | 0.072 [0.068, 0.077] / 0.036 / 0.010 | 5.7 | 0.153 / 0.101 / 0.036 |
| 100 | 0.119 / 0.074 / 0.026 | 9.4 | 0.271 / 0.210 / 0.096 |
| 300 | 0.174 / 0.122 / 0.048 | 13.7 | 0.382 / 0.319 / 0.166 |
| every candidate | 0.238 [0.229, 0.247] / 0.180 / 0.082 | 18.8 | 0.492 [0.480, 0.506] / 0.434 / 0.249 |

Timing ×10 and zero noise change nothing (N = 12: 0.042, 0.042; every candidate 0.241, 0.241).

### The prior-width slice, corrected (`prior_slice.csv`; 100 natural populations)

Reference, fraction right (Job 53's uncorrected value in brackets):

| N known | 0.1 pc | 0.3 | 1 | 3 | 10 | 30 | 100 | 300 | crosses 0.9 / 0.5 / 0.1 at (pc) |
|---|---|---|---|---|---|---|---|---|---|
| prior only | 0.988 | 0.595 | 0.197 | 0.066 | 0.020 | 0.007 | 0.002 | 0.001 | 0.15 / 0.37 / 2.0 |
| 7 | 0.988 | 0.745 (0.92) | 0.435 (0.63) | 0.227 (0.30) | 0.102 | 0.039 | 0.018 | 0.007 | 0.15 / 0.78 / 10 |
| **12** | 0.990 | **0.785 (0.96)** | **0.503 (0.72)** | **0.290 (0.39)** | 0.147 (0.18) | 0.068 | 0.033 | 0.011 | **0.16 / 1.0 / 19** |
| 30 | 0.993 | 0.838 (0.99) | 0.608 (0.85) | 0.409 (0.54) | 0.232 | 0.126 | 0.064 | 0.035 | 0.19 / 1.8 / 49 |
| every candidate | 0.996 | 0.937 (1.00) | 0.856 (1.00) | 0.733 (0.95) | 0.567 (0.73) | 0.413 (0.52) | 0.280 | 0.189 | 0.52 / 16 / above 300 |

Earth terms removed: N = 12: 0.871 / 0.659 / 0.448 / 0.267 at 0.3 / 1 / 3 / 10 pc (crosses 0.9 / 0.5 / 0.1 at 0.23 / 2.3 / 67 pc); every candidate: 0.979 / 0.957 / 0.915 / 0.840 (0.9 at 3.8 pc). All rows: `job54_logs/tables2.log`.

**Certification, corrected** (Job 52 Part 3 gate, unchanged; pulsars per sky of 79 at held-out purity; Job 53's value in brackets):

| Priors | Arm | N = 7 | N = 12 | every candidate |
|---|---|---|---|---|
| census | reference | none (none) | **none** (0.34) | 5.4 at 0.889 (9.1) |
| census | Earth terms removed | none (0.39) | 0.76 at 0.895 (1.5) | 26.3 at 0.881 (26.8) |
| σ_d = 1 pc | reference | 1.6 at 0.899 (25.4) | **3.9 at 0.907 (45.7)** | 70.8 at 0.899 (79.0) |
| σ_d = 1 pc | Earth terms removed | 13.1 at 0.872 (44.6) | 25.5 at 0.889 (62.0) | 79.0 at 0.957 |
| σ_d = 3 pc | reference | none (3.6) | **none** (10.5) | 45.5 at 0.902 (79.0) |
| σ_d = 3 pc | Earth terms removed | 1.1 at 0.889 (9.2) | 4.2 at 0.938 (20.1) | 78.9 at 0.917 |

Job 53's reading "with priors of 1 to 3 pc a useful certified catalogue exists (10 to 46 pulsars per sky at N = 12)" is withdrawn: it is 4 per sky at 1 pc and none at 3 pc.

---

## Part 3 — Does the per-bin model reproduce the skies?

**Set-up** (`job54_parts34.py`, `one_psr`). For each of 6,000 (pulsar, N = 12 known set) cases drawn from the 15,800 of the natural skies: the 12 sources keep their own tooth constants and their own shares β_unmod (retardation on); each is alone in a bin with a complex Gaussian holding the rest; joint fit on the fine grid (the `job53_heatmap.py` decoder with weights 2/P_k per source). Same true offsets as the sky run.

| Prior | Skies, Job 53 (phase defect) | **Skies, corrected** | **Per-bin model, real shares** | Model − corrected skies |
|---|---|---|---|---|
| census | 0.046 | 0.041 | 0.040 ± 0.003 | −0.001 |
| 0.1 pc | 0.999 | 0.990 | 0.990 | 0.000 |
| 0.3 pc | 0.963 | 0.785 | 0.784 | −0.001 |
| 1 pc | 0.723 | 0.503 | 0.505 | +0.002 |
| 3 pc | 0.388 | 0.290 | 0.292 | +0.002 |
| 10 pc | 0.173 | 0.147 | 0.144 | −0.003 |
| 30 pc | 0.077 | 0.068 | 0.069 | 0.000 |
| 100 pc | 0.032 | 0.033 | 0.025 | −0.007 |
| 300 pc | 0.013 | 0.011 | 0.012 | 0.000 |

- **The two agree within 0.007 everywhere.** The per-bin model with the real spread of shares explains the sky result; Job 53's attribution of the heatmap gap (0.038 against 0.077 at 30 pc) to the spread of shares is confirmed (0.069 with real shares).
- **How the defect was found.** Against Job 53's skies the model was low by 0.18 at 0.3 pc and 0.22 at 1 pc. Two refinements of the model did not close it: letting known sources that land in the same bin share one noise draw (3.8 of 12 do; moves the result by +0.03 at most), and combining the bins a pulsar term leaks into by inverse variance (raises the effective signal-to-noise by 5%; +0.01). The gap was in the skies, not the model.
- Source: `job54_logs/parts34.log`, `job54_logs/part3b.log`, `job54_logs/tables2.log`.

## Part 4 — The pile-up, split into "fainter" and "more crowded"

Corrected machine, natural skies, census priors, 200 populations; same known sets and draws in every arm (`job54_parts34.py`).

### 4.1 Two more arms

In both, every *other* source is retarded as in the reference. Fraction right, all pulsars / K ≥ 60:

| N known | off (nothing retarded) | others only | **frequency only** | **amplitude only** | on (reference) |
|---|---|---|---|---|---|
| 7 | 0.044 / 0.020 | 0.046 / 0.023 | 0.037 [0.034, 0.040] / 0.014 | 0.027 [0.024, 0.030] / 0.009 | 0.028 / 0.008 |
| 12 | 0.071 / 0.040 | 0.078 / 0.049 | 0.057 [0.054, 0.062] / 0.026 | 0.044 [0.041, 0.048] / 0.019 | 0.041 / 0.014 |
| every candidate | 0.389 / 0.339 | 0.398 / 0.348 | 0.321 [0.309, 0.333] / 0.264 | 0.267 [0.258, 0.277] / 0.209 | 0.238 / 0.180 |

Share of the loss ("others only" → "on") that each arm accounts for on its own:

| N known | frequency only (moved, not dimmed) | amplitude only (dimmed, not moved) |
|---|---|---|
| 7 | 0.48 | 1.03 |
| 12 | 0.55 | 0.90 |
| every candidate | 0.48 | 0.82 |

- The two do not add to 1: either change alone is enough to lose most of the pulsars that the other would also lose.
- The typical known pulsar term has 0.57 of its unretarded power (10th–90th percentile 0.15–0.94) and moves 2.7 bins (0.3–9.2).

### 4.2 Earth terms against pulsar terms (N = 12; share of the unmodelled power a known pulsar term competes with; 10th / median / 90th percentile over 189,600 pairs)

| Bin where it lands | Pairs | Pulsar terms of other binaries (discrete + half the sea) | Earth terms (discrete + half the sea) | Timing noise |
|---|---|---|---|---|
| **retardation on, all** | 189,600 | 0.424 / **0.473** / 0.502 | 0.495 / **0.523** / 0.570 | 0.000 / 0.002 / 0.012 |
| below 5 nHz | 49,465 | 0.451 / 0.487 / 0.504 | 0.496 / 0.513 / 0.548 | 0 |
| 5–10 nHz | 111,121 | 0.426 / 0.469 / 0.499 | 0.499 / 0.528 / 0.571 | 0.001 / 0.002 / 0.006 |
| 10–20 nHz | 29,010 | 0.386 / 0.457 / 0.508 | 0.472 / 0.525 / 0.590 | 0.009 / 0.015 / 0.038 |
| above 20 nHz | 4 | 0.307 / 0.404 / 0.572 | 0.302 / 0.479 / 0.554 | 0.116 / 0.124 / 0.144 |
| **retardation off, all** | 189,600 | 0.481 / **0.496** / 0.500 | 0.481 / 0.496 / 0.500 | 0.000 / 0.007 / 0.037 |

With retardation off the two are equal by construction (each source's pulsar term copies its Earth term). With it on, pulsar terms are slightly the smaller half. The sea is 0.14 of the total (0.07 in each half).

### 4.3 Is the destination bin louder because of the migration?

Per bin of the fit, all discrete sources, retardation on ÷ off, median over 15,800 pulsars:

| Bins | Pulsar-term power, on ÷ off | Range over the bins | Count of candidate pulsar terms, on ÷ off (mean count per pulsar) | All gravitational-wave power, on ÷ off |
|---|---|---|---|---|
| below 5 nHz (7 bins) | **0.97** | 0.90–1.02 | 1.015 (1,256 against 1,236) | 0.99 |
| 5–10 nHz (7 bins) | 0.82 | 0.73–0.90 | 0.75 (42 against 58) | 0.92 |
| 10–20 nHz (14 bins) | 0.61 | 0.49–0.70 | 0.25 (2 against 6) | 0.81 |
| above 20 nHz (19 bins) | 0.36 | 0.25–0.49 | 0 | 0.69 |

**No.** Retardation takes pulsar-term power out of every band above 5 nHz and does not add it below: 20 more candidate pulsar terms arrive below 5 nHz (of 1,240), each dimmer than when it left. The low-frequency end is crowded because the population is, with or without retardation.

**Reading for the talk.** The pulsar term does slide to lower frequency, and that alone costs about half the loss; but dimming alone costs nearly all of it, and the place it slides into is not made louder by other pulsar terms sliding with it. "An older, fainter copy that has to compete with the same crowd from a weaker position" is closer than "pulsar terms pile up".

---

## Verdicts

| | Prediction | Result | Verdict |
|---|---|---|---|
| W7 | Per-bin model with real shares agrees with the skies within 0.03 at every prior width | Against the corrected skies: within 0.007 at all nine. Against Job 53's skies: off by up to 0.22, which exposed the phase defect | **Holds**, after correcting the skies |
| W8 | "Amplitude only" 40–70% of the loss, "frequency only" 30–60%; pulsar terms 40–55% of what a known term competes with; per-bin pulsar-term power on within 20% of off below 20 nHz | amplitude only 82–103%; frequency only 48–55%; pulsar terms 47%; power on ÷ off 0.97 below 5 nHz, 0.82 at 5–10, 0.61 at 10–20 | **Two of four hold.** Amplitude matters more than predicted; pulsar-term power above 10 nHz falls by 40%, not within 20% |

## Limits

- The corrected prior-width slice rests on 100 populations (skies 100–199); everything else on 200.
- The enlarged-pool test (detectable sources above 15 nHz in the known set) used 40 skies and the uncorrected phases; it should be repeated corrected before the paper quotes N = 13 to 30 as "the detectable set".
- Jobs 51–53 reports are not edited; this report supersedes their tight-prior numbers.

## What was rerun, and cost

| Run (`job54_run.py`) | What | Wall time, 18–20 workers |
|---|---|---|
| `stat15`, `stat25` | ranking at 15 and 25 yr, no sky skipped | 3 min |
| `main`, `span15`, `span25` | ladder, split arms, diagnostics (uncorrected phases; used for tests and Part 3 inputs) | 46 min |
| `sigr` | corrected prior-width slice (100 skies) | 37 min |
| `span15r`, `span25r`, `mainr` | corrected ladder, curve, split arms (200 skies) | 62 min |
| `statbg`, `mainbg`, `inf` | side tests (40 and 60 skies) | 20 min |
