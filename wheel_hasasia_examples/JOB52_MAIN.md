# JOB 52 — The best case: how many pulsar distances does the sky allow? (Parts 1 and 2)

HEAD at the time of writing: `74cdf93` (branch `MM_playground`, in sync with `origin`). `wheel_hasasia_examples/` **is now committed** (420 files in that commit). All Job 52 files are uncommitted, written 2026-11-11.

## Five-line summary

1. **With the sources an array can actually detect granted as perfectly known (7 to 12 per sky), 4–5% of pulsars with real (census) priors get the right distance** (prior alone: 0.7%), and 1.4–2.1% of those with K ≥ 60. Better timing changes nothing.
2. **The curve rises steadily with the number of sources granted**: 1% at N = 1, 5% at 12, 8% at 30, 15% at 100, 27% at every candidate (about 1,170 per sky). There is no plateau and no threshold.
3. **Removing the unknown sources' Earth terms is the one thing that matters among the three missing ingredients.** It roughly doubles every figure (every candidate, census, K ≥ 60: 0.19 → 0.44). Sources above 30 nHz and letting frequency follow the trial distance each change K ≥ 60 by 0.02 or less.
4. **Best case with everything granted (Arm D, every candidate, census priors): 54% of all pulsars, 47% of K ≥ 60, 29% of K ≥ 600.** With only the detectable sources granted, Arm D gives 7–11% of all pulsars.
5. **Arm D gets worse, not better, when timing noise is removed** (census K ≥ 60: 0.47 → 0.32). The cause is the bins above 33 nHz: with those bins ignored the zero-noise result returns to the ×1 value. Quote the ×1 figure.

## Where the prompt and the code disagree

- **Almost no source above 30 nHz passes the detection cut.** 0.04–0.05 sources per sky (at most 1 in any sky). Arm B's premise, that fast-chirping sources above 30 nHz belong in the known list, is nearly empty in this population. The known set holds about 1% of the power in the 15–40 nHz bins with or without them.
- **The candidate list is capped at 1,300 by `job50_common.py:133`** (`if len(cand) > cand_cap: cand = cand[np.argsort(rp_all[cand])[::-1][:cand_cap]]`). Lifting the cap for Arm B adds 45–48 in-band sources per sky on average, far more than the frequency extension does. Arm B and Arm D therefore include that change too.
- **"Every candidate" is 1,161 (harness skies) or 1,175 (census skies) per sky at fixed amplitude**, not 1,127 (that was the per-pulsar median in Job 51).
- **Detectability marks at fixed amplitude are 6.5–6.8 (`arr_bg` ≥ 3) and 10.8–11.1 (`arr_bg_opt` ≥ 3)**, not 7.4 and 12.8 (those used the marginalised amplitude).
- **No pulsar has K < 3 from the finest tooth** of the booked source list, under either prior set. The two lowest classes of the class table are empty.
- **Prediction P7 is wrong in an unexpected direction** (summary point 5).

## Method and provenance

| Item | Where |
|---|---|
| Fit | `job52_common.py`, `Joint52.psr` (uncommitted, 2026-11-11), called from `job50_common.collect50` (modified 2026-11-11: three optional arguments, defaults unchanged) |
| Runs | `job52_run.py`: `stat_*`, `statx_*` (ranking), `main_*` (reference, A, C), `mainx_*` (B, D) → `job52_data/run_*.pkl` |
| Tables, figure | `job52_analysis.py` → `job52_logs/analysis_main.log`, `job52_data/analysis_main.json`, `job52_fig1_fraction_vs_N.png` |
| Checks | `job52_test.py` → `job52_logs/test.log`; `job52_diag_bins.py` → `job52_logs/diag_bins.log` |

Reference configuration as in the brief: retardation on, redshifted chirp masses, joint fit, fine grid, fixed amplitude 6.4 × 10⁻¹⁵, 79 red-noise-clean pulsars per sky.

- **Skies.** Harness priors: 100–199 (7,900 pulsars; 1,516 with K ≥ 60, 223 with K ≥ 600). Census priors: 100–119 only, for cost (1,580 pulsars; 1,318 with K ≥ 60, 452 with K ≥ 600). The same skies serve every arm within a prior set.
- **Same draws across arms.** Every arm of a sky uses the Job 51 random stream; Arms B and D rebuild the Job 51 draws for the sources and bins they share with it.
- **Timing levels** rescale the timing noise only; the candidate list stays the ×1 list (Job 51 let it grow at ×3 and ×10).
- **Rule 2 (reproduce first).** Every candidate, ×1: harness 0.238 / 0.076 / 0.803; census 0.188 / 0.080 / 0.271. Exact. On sky 100 the per-pulsar flags equal Job 51's for 79 of 79 pulsars under both prior sets.
- **Fit algebra.** The statistic is evaluated with one fast Fourier transform per arm instead of one per bin; it is the same statistic (the 79/79 check above would fail otherwise).
- **Grid.** Step ≤ 0.05 pc. The size cap binds for 5 of 7,900 harness pulsars (step up to 41 pc) and 0 of 1,580 census pulsars.
- **Intervals** are 95%, from resampling skies. K in the decode tables is from the coarsest tooth of the booked list.

**Ranking.** Sources are ranked by `arr_bg_opt` computed at fixed amplitude. Of the booked sources, 10.5 of 11.5 per sky (harness) and 10.3 of 11.6 (census) are in the top 12 of this ranking.

---

## Part 1 — The curve

### Census priors (main result), reference configuration, timing ×1

Prior only: all 0.007, K ≥ 60 0.003, K ≥ 600 0.000.

| N known | All pulsars | K ≥ 60 | K ≥ 600 | Right per sky (of which K ≥ 3 / K ≥ 10) | Share of band power |
|---|---|---|---|---|---|
| 1 | 0.011 [0.006, 0.016] | 0.002 [0.000, 0.005] | 0.000 | 0.8 (0.4 / 0.3) | 0.017 |
| 2 | 0.019 [0.013, 0.026] | 0.005 [0.002, 0.008] | 0.000 | 1.5 (0.8 / 0.7) | 0.043 |
| 3 | 0.026 [0.021, 0.032] | 0.008 [0.004, 0.011] | 0.004 [0.000, 0.011] | 2.0 (1.2 / 1.1) | 0.080 |
| 5 | 0.033 [0.025, 0.041] | 0.009 [0.004, 0.015] | 0.000 | 2.6 (1.6 / 1.4) | 0.138 |
| **7** | **0.042 [0.033, 0.053]** | **0.014 [0.008, 0.023]** | 0.000 | 3.3 (2.2 / 2.0) | 0.175 |
| **12** | **0.053 [0.045, 0.060]** | **0.021 [0.013, 0.030]** | 0.000 | 4.2 (2.9 / 2.5) | 0.195 |
| 20 | 0.068 [0.056, 0.079] | 0.032 [0.023, 0.041] | 0.007 [0.000, 0.014] | 5.3 (4.0 / 3.6) | 0.217 |
| 30 | 0.080 [0.066, 0.097] | 0.040 [0.032, 0.050] | 0.011 [0.002, 0.020] | 6.3 (5.0 / 4.5) | 0.242 |
| 60 | 0.117 [0.097, 0.139] | 0.065 [0.051, 0.082] | 0.013 [0.004, 0.024] | 9.2 (7.8 / 7.1) | 0.284 |
| 100 | 0.149 [0.128, 0.173] | 0.093 [0.074, 0.113] | 0.029 [0.012, 0.047] | 11.8 (10.3 / 9.4) | 0.315 |
| 300 | 0.206 [0.183, 0.231] | 0.136 [0.118, 0.156] | 0.042 [0.026, 0.059] | 16.2 (14.8 / 13.9) | 0.383 |
| every candidate (1,175) | 0.271 [0.241, 0.308] | 0.188 [0.162, 0.218] | 0.080 [0.049, 0.110] | 21.4 (19.9 / 18.8) | 0.447 |

Detectability marks (census skies): 6.5 sources per sky with `arr_bg` ≥ 3; 10.8 with `arr_bg_opt` ≥ 3.

**Timing, census priors (all pulsars / K ≥ 60 / K ≥ 600):**

| N known | ×1 | ×10 | noise ÷ 10⁶ |
|---|---|---|---|
| 7 | 0.042 / 0.014 / 0.000 | 0.039 / 0.014 / 0.000 | 0.039 / 0.014 / 0.000 |
| 12 | 0.053 / 0.021 / 0.000 | 0.055 / 0.024 / 0.000 | 0.054 / 0.022 / 0.000 |
| 30 | 0.080 / 0.040 / 0.011 | 0.085 / 0.046 / 0.015 | 0.082 / 0.044 / 0.013 |
| every candidate | 0.271 / 0.188 / 0.080 | 0.277 / 0.195 / 0.084 | 0.275 / 0.193 / 0.082 |

### Harness priors (for comparison), timing ×1

Prior only: all 0.130, K ≥ 60 0.003, K ≥ 600 0.000.

| N known | All pulsars | K ≥ 60 | K ≥ 600 | Right per sky (K ≥ 3 / K ≥ 10) | Share of band power |
|---|---|---|---|---|---|
| 1 | 0.275 [0.259, 0.292] | 0.003 [0.001, 0.006] | 0.000 | 21.7 (0.7 / 0.1) | 0.017 |
| 2 | 0.348 [0.334, 0.365] | 0.005 [0.002, 0.009] | 0.000 | 27.5 (0.9 / 0.2) | 0.045 |
| 3 | 0.403 [0.389, 0.420] | 0.007 [0.004, 0.012] | 0.000 | 31.9 (1.1 / 0.4) | 0.082 |
| 5 | 0.462 [0.450, 0.474] | 0.007 [0.003, 0.011] | 0.000 | 36.5 (1.3 / 0.6) | 0.145 |
| 7 | 0.495 [0.484, 0.506] | 0.018 [0.011, 0.025] | 0.000 | 39.1 (1.7 / 0.8) | 0.169 |
| 12 | 0.535 [0.524, 0.546] | 0.033 [0.024, 0.044] | 0.004 [0.000, 0.017] | 42.2 (2.3 / 1.3) | 0.194 |
| 20 | 0.571 [0.560, 0.582] | 0.032 [0.024, 0.041] | 0.004 [0.000, 0.017] | 45.1 (2.7 / 1.5) | 0.218 |
| 30 | 0.612 [0.601, 0.622] | 0.055 [0.043, 0.067] | 0.013 [0.000, 0.031] | 48.3 (3.3 / 2.0) | 0.242 |
| 60 | 0.662 [0.651, 0.672] | 0.078 [0.066, 0.089] | 0.027 [0.005, 0.053] | 52.3 (4.5 / 2.9) | 0.284 |
| 100 | 0.697 [0.686, 0.707] | 0.110 [0.093, 0.128] | 0.036 [0.014, 0.063] | 55.0 (5.4 / 3.7) | 0.316 |
| 300 | 0.766 [0.757, 0.775] | 0.175 [0.155, 0.195] | 0.049 [0.023, 0.079] | 60.5 (7.3 / 5.3) | 0.381 |
| every candidate (1,161) | 0.803 [0.794, 0.812] | 0.238 [0.216, 0.260] | 0.076 [0.037, 0.122] | 63.5 (9.3 / 7.1) | 0.444 |

Timing (all / K ≥ 60), harness: N = 7: 0.495 / 0.018 at ×1, 0.496 / 0.016 at zero noise; N = 12: 0.535 / 0.033 and 0.536 / 0.030; N = 30: 0.612 / 0.055 and 0.613 / 0.055; every candidate: 0.803 / 0.238 and 0.804 / 0.242. Detectability marks: 6.8 and 11.1.

**Figure 1:** `job52_fig1_fraction_vs_N.png` (census priors; all pulsars and K ≥ 60 against N; three timing levels; prior-only dotted; grey band between the two detectability marks).

Reading: the harness "all pulsars" column is carried by pulsars that are nearly decided by their prior (right per sky with K ≥ 10: 1.3 of 42 at N = 12). Under census priors almost every correct decode is a K ≥ 10 pulsar.

---

## Part 2 — Is the best case a true upper bound?

**Arms.**

- **A — unknown Earth terms removed.** The Earth terms of every source outside the known set are deleted from the data and from the unmodelled variance. The sea is split in two equal halves (the machine books it as `2*bsea`, Earth plus pulsar: `job50_common.py`, `X_on = np.maximum(gt-own+2*bsea(fPc), 0)*rt/Pt_on`); the Earth half is deleted, the pulsar half kept as the same random draw at half the power.
- **B — pool to 100 nHz in Earth frequency, bins to 100 nHz** (142 bins instead of 47), candidate cap lifted (see top).
- **C — retarded frequency, amplitude, leakage and phase follow the trial distance.** Phase is the exact integral of the retarded frequency along L. Amplitude and leakage are recomputed in steps of 0.01 bin of the fastest-moving source (up to 4,000 steps across ±3σ_d; the cap never binds under census priors and binds for 1 of 7,900 harness pulsars); phase is exact at every grid point. Trial distances below zero are clipped to zero.
- **D — A, B and C together.**

### Census priors (all pulsars / K ≥ 60 / K ≥ 600)

| N known | Timing | Reference | A | B | C | **D** |
|---|---|---|---|---|---|---|
| 7 | ×1 | 0.042 / 0.014 / 0.000 | 0.066 / 0.031 / 0.007 | 0.042 / 0.017 / 0.000 | 0.041 / 0.017 / 0.004 | **0.073 / 0.035 / 0.013** |
| 7 | zero | 0.039 / 0.014 / 0.000 | 0.067 / 0.028 / 0.002 | 0.037 / 0.010 / 0.000 | 0.043 / 0.020 / 0.004 | 0.053 / 0.024 / 0.004 |
| 12 | ×1 | 0.053 / 0.021 / 0.000 | 0.098 / 0.054 / 0.011 | 0.054 / 0.024 / 0.002 | 0.058 / 0.024 / 0.002 | **0.106 / 0.059 / 0.027** |
| 12 | zero | 0.054 / 0.022 / 0.000 | 0.103 / 0.054 / 0.013 | 0.046 / 0.014 / 0.000 | 0.055 / 0.022 / 0.002 | 0.065 / 0.032 / 0.015 |
| 30 | ×1 | 0.080 / 0.040 / 0.011 | 0.177 / 0.115 / 0.040 | 0.085 / 0.043 / 0.011 | 0.085 / 0.042 / 0.004 | **0.191 / 0.136 / 0.060** |
| 30 | zero | 0.082 / 0.044 / 0.013 | 0.185 / 0.124 / 0.038 | 0.069 / 0.031 / 0.011 | 0.085 / 0.046 / 0.009 | 0.117 / 0.071 / 0.024 |
| every candidate | ×1 | 0.271 / 0.188 / 0.080 | 0.509 / 0.436 / 0.246 | 0.277 / 0.197 / 0.091 | 0.286 / 0.209 / 0.091 | **0.540 / 0.472 / 0.290** |
| every candidate | zero | 0.275 / 0.193 / 0.082 | 0.516 / 0.445 / 0.265 | 0.261 / 0.184 / 0.071 | 0.289 / 0.212 / 0.088 | 0.385 / 0.320 / 0.181 |

95% intervals for Arm D, ×1: N = 7: [0.058, 0.087] / [0.025, 0.045] / [0.004, 0.024]; N = 12: [0.085, 0.127] / [0.041, 0.080] / [0.011, 0.043]; N = 30: [0.163, 0.223] / [0.107, 0.168] / [0.030, 0.090]; every candidate: [0.502, 0.582] / [0.428, 0.519] / [0.234, 0.346]. All other intervals are in `job52_logs/analysis_main.log`.

### Harness priors (all pulsars / K ≥ 60 / K ≥ 600)

| N known | Timing | Reference | A | B | C | **D** |
|---|---|---|---|---|---|---|
| 7 | ×1 | 0.495 / 0.018 / 0.000 | 0.540 / 0.033 / 0.009 | 0.495 / 0.020 / 0.004 | 0.494 / 0.020 / 0.000 | **0.540 / 0.034 / 0.027** |
| 7 | zero | 0.496 / 0.016 / 0.004 | 0.542 / 0.032 / 0.009 | 0.465 / 0.015 / 0.004 | 0.495 / 0.016 / 0.004 | 0.385 / 0.015 / 0.013 |
| 12 | ×1 | 0.535 / 0.033 / 0.004 | 0.594 / 0.063 / 0.018 | 0.533 / 0.032 / 0.004 | 0.534 / 0.031 / 0.000 | **0.596 / 0.065 / 0.018** |
| 12 | zero | 0.536 / 0.030 / 0.004 | 0.598 / 0.061 / 0.009 | 0.503 / 0.018 / 0.000 | 0.535 / 0.026 / 0.000 | 0.438 / 0.028 / 0.009 |
| 30 | ×1 | 0.612 / 0.055 / 0.013 | 0.683 / 0.123 / 0.027 | 0.608 / 0.050 / 0.013 | 0.611 / 0.052 / 0.022 | **0.683 / 0.135 / 0.072** |
| 30 | zero | 0.613 / 0.055 / 0.013 | 0.687 / 0.135 / 0.045 | 0.591 / 0.045 / 0.018 | 0.613 / 0.053 / 0.018 | 0.559 / 0.071 / 0.045 |
| every candidate | ×1 | 0.803 / 0.238 / 0.076 | 0.880 / 0.491 / 0.296 | 0.805 / 0.243 / 0.090 | 0.803 / 0.237 / 0.085 | **0.889 / 0.532 / 0.395** |
| every candidate | zero | 0.804 / 0.242 / 0.081 | 0.881 / 0.499 / 0.278 | 0.796 / 0.222 / 0.094 | 0.804 / 0.240 / 0.094 | 0.829 / 0.338 / 0.215 |

### What each arm shows

**Arm A.** Earth-term power is half or more of what the reference treats as noise. Share of the unmodelled gravitational-wave power that was Earth-term (median over bins and pulsars; every candidate known; census / harness): 0.50 / 0.50 in bins at or below 10 nHz; 0.68 / 0.71 in bins 15–33 nHz; 0.70 / 0.73 in bins 15–40 nHz (extended configuration). For N = 7 to 30 the low-frequency figure is 0.51–0.53. Arm A assumes those Earth terms are removed perfectly, with no error and no loss of signal; it is a bound, not an analysis.

**Arm B.** Sources per sky with Earth frequency above 30 nHz: (i) passing the booked detection cut against timing noise 0.04 (harness) and 0.05 (census), never more than 1 in a sky; (ii) `arr_bg` ≥ 3: the same 0.04 / 0.05; (iii) `arr_bg_opt` ≥ 3: the same. They sit at 30.7–32.4 nHz, and for the first pulsar of their sky the pulsar term lands at 6–8.5 nHz. Share of each 15–40 nHz bin's power held by the known set (median over bins and pulsars): 0.006 at N = 7 and 12, 0.007 at N = 30, 0.012 at every candidate; the reference configuration gives the same to within 0.001. So the pile-up band stays unmodelled whatever is granted.

**Arm C.** Fraction of (pulsar, source) pairs whose retarded frequency moves by more than 0.1 bin across ±3σ_d: census 0.34 of every-candidate pairs, 0.90 of top-30 pairs, 0.96 of top-12 pairs (median largest move per pulsar 4.7 bins); harness 0.05, 0.25 and 0.28 (median largest move 0.02 bins). Despite that, the decode fraction barely changes: the template and the data move together.

**Arm D at zero timing noise is lower than at ×1.** Test (`job52_diag_bins.py`, harness skies 100–103, 316 pulsars, zero timing noise): with bins to 100 nHz, D every candidate gives 0.835 and D N = 12 gives 0.443; with the bins above 33 nHz ignored the same runs give 0.905 and 0.585, equal to the ×1 values (0.902, 0.579). Arm B behaves the same way (0.801 → 0.816; 0.491 → 0.554). The fit treats bins as independent. Above 33 nHz, once timing noise is gone, what is left is the leakage of low-frequency sources, which is correlated from bin to bin (measured neighbour correlation 0.19 against 0.11 expected for independent bins), so 100 extra bins are counted as 100 independent measurements. This is a limit of the per-bin noise model, not physics. **The ×1 column is the one to quote.**

### By class of K (finest tooth of the booked source list)

No pulsar falls below K = 3 in either prior set. Fraction decoded correctly [95% interval]:

**Census priors** (pulsars per class: 2 / 34 / 495 / 1,049)

| Arm | K 3–10 (n = 2) | K 10–60 (n = 34) | K 60–600 (n = 495) | K ≥ 600 (n = 1,049) |
|---|---|---|---|---|
| prior only | 0.000 | 0.176 | 0.008 | 0.001 |
| reference, N = 7 | 0.500 | 0.676 [0.47, 0.85] | 0.065 [0.05, 0.08] | 0.010 [0.00, 0.02] |
| D, N = 7 | 0.500 | 0.794 [0.63, 0.94] | 0.127 [0.09, 0.16] | 0.023 [0.01, 0.03] |
| reference, N = 12 | 0.500 | 0.824 [0.66, 0.97] | 0.079 [0.06, 0.10] | 0.014 [0.01, 0.02] |
| D, N = 12 | 0.500 | 0.824 [0.67, 0.97] | 0.192 [0.16, 0.23] | 0.042 [0.03, 0.06] |
| reference, every candidate | 1.000 | 0.941 [0.88, 1.00] | 0.457 [0.40, 0.53] | 0.160 [0.13, 0.19] |
| **D, every candidate** | 1.000 | 1.000 | **0.760 [0.71, 0.81]** | **0.420 [0.38, 0.47]** |

**Harness priors** (pulsars per class: 285 / 5,270 / 1,466 / 879)

| Arm | K 3–10 | K 10–60 | K 60–600 | K ≥ 600 |
|---|---|---|---|---|
| prior only | 0.161 | 0.184 | 0.008 | 0.003 |
| reference, N = 7 | 0.719 [0.67, 0.77] | 0.684 [0.67, 0.70] | 0.061 [0.05, 0.07] | 0.013 [0.01, 0.02] |
| D, N = 7 | 0.772 [0.73, 0.82] | 0.736 [0.72, 0.75] | 0.100 [0.09, 0.12] | 0.024 [0.01, 0.03] |
| reference, N = 12 | 0.740 [0.69, 0.79] | 0.733 [0.72, 0.75] | 0.092 [0.08, 0.11] | 0.018 [0.01, 0.03] |
| D, N = 12 | 0.821 [0.78, 0.86] | 0.796 [0.78, 0.81] | 0.163 [0.14, 0.18] | 0.042 [0.03, 0.05] |
| reference, every candidate | 0.996 | 0.998 | 0.443 [0.41, 0.47] | 0.174 [0.15, 0.20] |
| **D, every candidate** | 1.000 | 1.000 | **0.705 [0.67, 0.73]** | **0.491 [0.45, 0.53]** |

By class the two prior sets agree: the difference in the all-pulsar figures (0.27 against 0.80) is the class mix, 98% of census pulsars against 30% of harness pulsars having K ≥ 60 on the finest tooth.

### Where the harness priors come from

`job50_common.py:118-119` (same in `job43c_rebooking.py`, `collect36m`):

```
px = 2.34/np.cos(beta)**2*Dk**2*10.0*np.sqrt(30.0/T); pb = 3.0*(T/30.0)**(-2.5)
pri = np.maximum(np.where(rng.random(Np) < 0.7, np.minimum(px, pb), px), 0.05)
```

- `px` is a timing-parallax distance error: 2.34 pc × (D/kpc)² / cos²β at 10 ns, scaled to 100 ns and by √(30 yr / T). At T = 45 yr and D = 1 kpc it is about 19 pc / cos²β. Distances are uniform in 0.5–3 kpc, ecliptic latitudes isotropic.
- `pb` is an orbital-period-derivative distance error: 3 pc at 30 yr, falling as T^−2.5, so 1.09 pc at 45 yr. 70% of pulsars are given it.
- The report labels the span law, the binary fraction and the period-derivative law as demonstration values (`ACTIVE_WHEEL_REPORT.md:717-718`: "parallax ∝ 1/√N_i, Pbdot ∝ T_i^−2.5, binary fraction 0.7, all [DEMO]").

Consistency with the simulated array: the parallax term uses the array's own 100 ns and 45 yr, so it is consistent. The period-derivative term is not derived from the array at all: it asserts that 70% of the 120 pulsars are binaries with distances good to about 1 pc. The census of real pulsars has a median distance error of 242 pc (10th percentile 31 pc). The harness priors are therefore the simulated array's parallax plus an optimistic assumption about binaries, and the census is the realistic set.

---

## Verdicts on the predictions

| | Prediction | Result | Verdict |
|---|---|---|---|
| P1 | N = 12, census, ×1: all between 0.08 and 0.20; K ≥ 60 below 0.10 | all 0.053 [0.045, 0.060]; K ≥ 60 0.021 | **Half wrong.** All-pulsar fraction is below the predicted range; K ≥ 60 holds |
| P2 | N ≤ 30: all-pulsar change below 0.02 from ×1 to zero noise | largest change 0.003 (census), 0.001 (harness) | **Holds** |
| P3 | rises with N at every step; reaches Job 51's values | all-pulsar fraction rises at every step in both prior sets and ends on Job 51's values exactly. Harness K ≥ 60 is flat from N = 12 to 20 (0.033, 0.032) | **Holds** for all pulsars; one flat step at K ≥ 60 within noise |
| P4 | Arm A raises K ≥ 60 (every candidate, harness) from 0.24 to 0.35–0.60 | 0.491 [0.463, 0.517] | **Holds** |
| P5 | Arm B changes K ≥ 60 by less than 0.05 | +0.005 (harness), +0.009 (census) at ×1; −0.020, −0.009 at zero noise | **Holds** |
| P6 | Arm C changes K ≥ 60 by less than 0.03, K ≥ 600 by less than 0.10 | K ≥ 60: −0.001 (harness), +0.021 (census); K ≥ 600: +0.009, +0.011 | **Holds** |
| P7 | Arm D flat with timing (change below 0.05) | every candidate K ≥ 60: 0.532 → 0.338 (harness), 0.472 → 0.320 (census) | **Wrong.** It falls. Cause: bins above 33 nHz counted as independent (test above) |

## Limits

- Census results rest on 20 skies; harness on 100.
- Arm A removes unknown Earth terms perfectly. A real array removes them imperfectly and loses some pulsar-term signal doing so.
- The fit treats bins as independent, which fails above 33 nHz once timing noise is removed.
- Arms B and D also lift the 1,300-source cap on the candidate list.
- "Detectable" means array ratio ≥ 3 against the background under the Job 51 approximation (independent background in each pulsar).

## The two closing sentences

- **With perfect knowledge of the sources an array can detect (N = 7 to 12), 4.2–5.3% of pulsars with real priors get the right distance (prior alone 0.7%; 7.3–10.6% if the unknown Earth terms are also removed perfectly), and better timing does not change it (largest change 0.003 down to zero timing noise).**
- **With everything a real analysis could do granted (Arm D, every candidate known), 47% [43, 52] of wide-prior (K ≥ 60) pulsars with census priors get the right distance, and 29% [23, 35] of those with K ≥ 600; with only the detectable sources known, Arm D gives 3.5–5.9% of K ≥ 60.**
