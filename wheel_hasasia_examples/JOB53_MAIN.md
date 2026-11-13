# JOB 53 — Parts 1 to 3: the pile-up split, the heatmap, the final curve

HEAD: `74cdf93` (branch `MM_playground`). All Job 52 and Job 53 files are uncommitted (staged), written 2026-11-11/12. Part 0 is in `JOB53_SKIES.md`.

**Natural skies here are 200 distinct populations** (skies 100–299, seeds 5100–5299; 15,800 pulsars; 13,255 with K ≥ 60; 4,638 with K ≥ 600). Constructed skies are Job 52's 20 (skies 100–119). Census priors unless stated. Intervals are 95%, from resampling populations.

## Five-line summary

1. **Final numbers on natural skies.** With the 7 to 12 loudest sources known perfectly, **2.9% to 4.6% of pulsars get the right distance** (prior alone 0.6%); 0.7% to 1.6% of those with K ≥ 60. With every candidate known (1,300 per sky): 26.0%. Better timing changes nothing (largest change 0.001).
2. **Retardation costs a factor 1.6 to 1.7 on natural skies, and all of it is the known source's own pulsar term.** Retarding only the known sources reproduces the full loss; retarding only the other sources costs nothing (it helps slightly). "Pulsar terms pile up on the known source" is the wrong description.
3. **What hurts a known source is that its own pulsar term arrives weaker (median ×0.57 in power) at a lower frequency where the bin is louder (×1.5).** Other sources' pulsar terms arriving from more than one bin higher are 8% of what it competes with. Earth terms of unknown sources are 45%.
4. **The heatmap with the joint fit on a fine grid.** For 4 sources the 90% contour is at 2.6 pc for a share of 0.8 and 0.34 pc for a share of 0.5. The typical known pulsar term on natural skies has a share of 0.055 (90th percentile 0.18): there, 90% is reached only where the prior alone nearly decides (σ_d below 0.2 to 0.3 pc), for 4 sources or 12.
5. **Against prior width, natural skies, N = 12:** 96% right at σ_d = 0.3 pc, 72% at 1 pc, 39% at 3 pc, 17% at 10 pc, 8% at 30 pc. A 0.90-pure catalogue exists for σ_d = 1 pc (46 of 79 pulsars per sky) and 3 pc (10 per sky), and essentially not for census priors (0.3 per sky, all with K < 3).

## Where the prompt and the code disagree

- **"It costs a factor 2.7 at N = 12"** is the constructed-sky figure, and part of it was an artefact: the three raised sources were given too large a chirp mass (`JOB53_SKIES.md`), so their pulsar terms were over-retarded. On natural skies with each source's own mass the factor is **1.68**.
- **"As Job 47c" for the heatmap sources** means frequencies taken from the booked source lists of the 50 *constructed* populations. I kept that (only the frequencies are used; shares are set by hand), and say so.
- **Job 49's D-mf panel did not use a complex Gaussian** for the rest of the bin; it used each comb's actual mix of contributors as explicit phasors. The new heatmap uses the Gaussian, as the brief asks. The reproduction below is of Job 49 as it was.
- **Natural skies give 10.2 sources per sky detectable with optimal combination**, not fewer than 8 (V1).

## Method and provenance

| Item | Where (all uncommitted, 2026-11-12) |
|---|---|
| Sky switch | `job50_common.collect50(sky=…, mass=…)`, `natural_pop` |
| Fit | `job53_common.Joint53` (Job 52's joint fit with the retardation of the known set and of the other sources set separately; reproduces Job 52 on 79/79 pulsars in 26 arms, `job53_logs/test.log`) |
| Runs | `job53_run.py`: `stat_nat`, `main_nat`, `sig_nat` (natural); `main_con` (constructed); `p1x_nat20` |
| Parts 1 and 3 | `job53_main.py` → `job53_logs/main.log`, `job53_data/main.json`, Figures 3 and 4 |
| Part 2 | `job53_heatmap.py` → `job53_logs/heat_*.log`, `job53_data/heat_*.{npz,json}`, Figures 1 and 2 |

Reproduction (rule 2), constructed census skies 100–119: reference N = 12: 0.053 / 0.021 / 0.000; every candidate: 0.271 / 0.188 / 0.080; Arm A, N = 12: 0.098 / 0.054 / 0.011. Exact.

Natural skies, 200 populations: 4.7 sources per sky detectable against the background (`arr_bg` ≥ 3; 10th to 90th percentile 2 to 7); 10.2 with optimal combination (6 to 14); 1.7 sources per sky rejected by the upper limit in the draw (at most 6); 1,300 candidates per sky.

---

## Part 1 — What retardation costs, and why

### The four arms (natural skies; all pulsars / K ≥ 60 / K ≥ 600 | right per sky of 79)

| N known | off | known only | others only | on (reference) |
|---|---|---|---|---|
| 7 | 0.046 [0.042, 0.050] / 0.019 / 0.005 \| 3.6 | 0.030 [0.027, 0.033] / 0.008 / 0.002 \| 2.4 | 0.050 [0.046, 0.054] / 0.024 / 0.008 \| 3.9 | **0.029 [0.026, 0.032] / 0.007 / 0.001 \| 2.3** |
| 12 | 0.077 [0.072, 0.084] / 0.042 / 0.017 \| 6.1 | 0.045 [0.042, 0.049] / 0.016 / 0.004 \| 3.6 | 0.085 [0.079, 0.092] / 0.049 / 0.024 \| 6.7 | **0.046 [0.043, 0.049] / 0.016 / 0.003 \| 3.6** |
| every candidate | 0.407 [0.394, 0.420] / 0.349 / 0.230 \| 32.1 | 0.261 [0.252, 0.270] / 0.195 / 0.083 \| 20.6 | 0.418 [0.404, 0.433] / 0.361 / 0.253 \| 33.0 | **0.260 [0.252, 0.269] / 0.194 / 0.080 \| 20.6** |

Prior only: 0.006 / 0.001 / 0.000. All intervals are in `job53_logs/main.log`.

**Constructed skies, N = 12:** off 0.141 / 0.104 / 0.064; known only 0.059 / 0.023 / 0.002; others only 0.166 / 0.126 / 0.069; on 0.053 / 0.021 / 0.000. Same pattern; the loss is larger (factor 2.7) for the reason given at the top.

- "known only": the known sources' pulsar terms sit at their retarded frequency and amplitude; every other source's pulsar term is left at its Earth frequency and amplitude.
- "others only": the reverse.
- The known set (ranking), pulsars, distances, angles and noise draws are the same in all four arms.

**Attribution** (all pulsars, natural):

| N known | off ÷ on | share of the loss from retarding the known set alone | share from retarding the others alone |
|---|---|---|---|
| 7 | 1.57 | 0.93 | −0.23 |
| 12 | 1.68 | 1.03 | −0.25 |
| every candidate | 1.56 | 0.99 | −0.08 |

A negative share means that arm does better than "off": a retarded pulsar term is weaker than an unretarded one, so retarding the *other* sources slightly lowers what the known source competes with.

### Per-pair bookkeeping (N = 12; every pulsar and known source; "on" ÷ "off")

| Quantity | Natural (189,600 pairs), 10th / 50th / 90th | Constructed (18,960 pairs) |
|---|---|---|
| Own pulsar-term power | 0.148 / **0.572** / 0.940 | 0.107 / 0.454 / 0.887 |
| Unmodelled power it competes with (leakage-weighted) | 0.969 / **1.407** / 3.536 | 0.980 / 1.637 / 4.030 |
| Share β_unmod | 0.054 / **0.457** / 0.963 | 0.041 / 0.339 / 0.903 |
| Bins the pulsar term moves | 0.27 / **2.7** / 9.2 | 0.57 / 3.5 / 10.4 |
| β_unmod itself, off | 0.065 / 0.137 / 0.363 | 0.072 / 0.202 / 0.473 |
| β_unmod itself, on | 0.008 / 0.055 / 0.180 | 0.009 / 0.056 / 0.229 |

74% of natural pairs move by more than one bin; 96% by more than a tenth of a bin.

**Splitting the rise in competing power** (`p1x_nat20`, natural skies 100–119, 18,960 pairs):

| Step | Competing power, 10th / 50th / 90th |
|---|---|
| Move the known pulsar term to its retarded bin, leave every other pulsar term where it was | × 1.01 / **1.50** / 3.60 |
| Then retard every other pulsar term too | × 0.87 / **0.95** / 1.01 |
| Both (on ÷ off) | × 0.97 / 1.40 / 3.36 |

So the rise is the move to a lower-frequency bin, where the background is louder. The other sources' retardation lowers it by 5%.

### What the unmodelled power is made of where a known pulsar term lands (retardation on, N = 12, natural; median share over pairs)

| Bin frequency | Pairs | Earth terms of outside sources | Outside pulsar terms from more than one bin higher | Outside pulsar terms moved one bin or less | Sea | Timing noise |
|---|---|---|---|---|---|---|
| all | 189,600 | **0.450** | **0.077** | 0.308 | 0.144 | 0.002 |
| below 5 nHz | 49,465 | 0.396 | 0.025 | 0.328 | 0.240 | 0.000 |
| 5–10 nHz | 111,121 | 0.459 | 0.087 | 0.312 | 0.136 | 0.002 |
| 10–20 nHz | 29,010 | 0.483 | 0.183 | 0.226 | 0.082 | 0.015 |
| above 20 nHz | 4 | 0.468 | 0.208 | 0.151 | 0.022 | 0.124 |

Constructed skies give the same within 0.04 in every cell (`job53_logs/main.log`). Only 4 of 189,600 known pulsar terms land above 20 nHz.

### The two sentences

- **Of the loss from retardation, all of it (93% to 103%) comes from the known source's own pulsar term arriving weaker (median 0.57 of its unretarded power) and in a lower, louder bin (competing power × 1.5); none comes from other sources' pulsar terms arriving, which if anything help by 5%.**
- **"Pulsar terms pile up" is not the right description of what hurts a detectable source: arrivals from more than one bin away are 8% of what its pulsar term competes with; the right description is that the source's own pulsar term is an older, weaker, lower-frequency copy of the signal that lands where unknown Earth terms (45%), unknown pulsar terms that barely moved (31%) and the sea (14%) are already louder than it is.**

---

## Part 2 — The heatmap

**Set-up** (`job53_heatmap.py`):

- **Data.** N_s sources, each alone in its own bin: the source's pulsar-term phasor plus a complex Gaussian holding the rest of the bin, so the source's share is exactly β. No timing noise, no template error.
- **Decoder.** The joint fit. With one unit-amplitude source per bin its statistic reduces exactly to the matched filter of Job 49 ("D-mf"); the difference from Job 49 is the grid: step ≤ 0.05 pc and ≤ one sixth of the shortest tooth spacing, instead of 601 points across ±3σ_d.
- **Sources, as Job 47c.** Frequencies: the N_s highest-ranked sources (by information, `nats`) of the machine's booked list on population t mod 50, for realisation t (constructed populations; if a list is shorter than N_s the next population with enough is used: 38 of 50 have 12). Median frequency 6.1 to 10.3 nHz. Angle between pulsar and source: cos μ uniform in (−0.9, 0.9), drawn per source per realisation. Tooth spacing from κ = 2π f (1 − cos μ)/c.
- **Axes.** σ_d: 20 log steps from 0.1 to 300 pc. β: 20 values from 0.02 to 0.98 in equal steps of ln[β/(1 − β)] (0.02, 0.03, 0.044, 0.065, 0.095, 0.137, 0.193, 0.264, 0.351, 0.449 and their mirror images), so both ends are resolved.
- **Statistics.** 500 realisations per cell (standard error at most 0.022); cells drawn raw; contours from the 3 × 3 mean. Crossing tables from dedicated rows: 40 log steps in σ_d, 1,000 realisations per cell.
- **K per parsec of σ_d** (median): N_s = 2: 1.9 (coarsest tooth), 5.8 (finest); N_s = 4: 1.4, 10.0; N_s = 8: 1.1, 13.4; N_s = 12: 1.0, 14.2.

### Tie to Job 49

- **Reproduction.** Job 49's D-mf panel (`job49_data/D1_grid.npz`, 250 realisations, 601-point grid, 4 combs): largest difference over 400 cells **0** (`job53_logs/heat_repro.log`).
- **What the fine grid changes, same data, cell by cell:** 307 of 400 cells change; mean +0.035; largest gain +0.756; largest loss −0.028.

| σ_d | Coarse-grid step | Mean change | Largest change |
|---|---|---|---|
| below 0.6 pc | 0.001–0.005 pc | −0.001 | +0.012 |
| 0.6–3 pc | 0.008–0.029 pc | 0.000 | +0.020 |
| 3–20 pc | 0.044–0.16 pc | +0.014 | +0.268 |
| above 20 pc | 0.24–3.0 pc | **+0.091** | **+0.756** |

The coarse grid was adequate below about 5 pc and lost real decodes above 10 pc, where its step exceeded the 0.25 pc tolerance. The 90% crossings Job 49 quoted (all below 3 pc) are unchanged: β = 0.79: 2.33 → 2.22 pc; 0.51: 0.31 → 0.33; 0.28: 0.18 → 0.18; 0.19: 0.16 → 0.16.

### σ_d at which the 90% contour crosses each β (pc; in brackets the 50% crossing)

| N_s | β = 0.95 | 0.8 | 0.5 | 0.3 | 0.2 | 0.1 |
|---|---|---|---|---|---|---|
| 2 | 0.93 (5.8) | 0.41 (2.9) | 0.18 (1.0) | 0.15 (0.6) | 0.14 (0.5) | 0.15 (0.5) |
| **4** | **23.9 (275)** | **2.64 (34)** | **0.34 (3.5)** | 0.17 (1.3) | 0.16 (0.8) | 0.15 (0.5) |
| 8 | above 300 (1.00 at 300 pc) | 255 (above 300) | 1.64 (41) | 0.32 (3.9) | 0.21 (1.7) | 0.15 (0.7) |
| 12 | above 300 (1.00) | above 300 (1.00) | 11.0 (above 300) | 0.68 (13.9) | 0.28 (3.6) | 0.17 (1.1) |

The prior alone gives 90% at σ_d = 0.152 pc. So for β ≤ 0.2 the 90% contour sits on the prior-alone boundary for any N_s up to 8, and within a factor 2 of it for 12.

### Fraction decoded along the two natural-sky values of β (4,000 realisations per cell)

| N_s | β | σ_d = 0.3 pc | 1 | 3 | 10 | 30 | 100 | 300 |
|---|---|---|---|---|---|---|---|---|
| 4 | 0.055 (median known pulsar term) | 0.638 | 0.276 | 0.127 | 0.043 | 0.018 | 0.004 | 0.002 |
| 4 | 0.180 (90th percentile) | 0.726 | 0.426 | 0.216 | 0.090 | 0.043 | 0.014 | 0.004 |
| 12 | 0.055 | 0.708 | 0.384 | 0.198 | 0.091 | 0.038 | 0.012 | 0.005 |
| 12 | 0.180 | 0.875 | 0.687 | 0.482 | 0.291 | 0.176 | 0.083 | 0.044 |

Prior only: 0.595, 0.197, 0.066, 0.020, 0.007, 0.002, 0.001.

Check against the skies: 12 sources each at the median share give 0.038 at 30 pc and 0.012 at 100 pc; the sky result for N = 12 with every prior at 30 pc is 0.077 and at 100 pc 0.032 (Part 3b). The skies do a little better because shares are spread (a few pairs have large shares), which the equal-share heatmap does not represent.

**Figure 1:** `job53_fig1_heatmaps.png` (four panels, one colour scale; 50% and 90% contours; top axis K from the coarsest tooth, with the factor to the finest tooth; on the N_s = 4 and 12 panels: red lines at β = 0.055 and 0.180, orange marks at 31 and 242 pc, green dotted line where the prior alone reaches 90%).
**Figure 2:** `job53_fig2_heatmap_Ns4.png` (N_s = 4 alone).

---

## Part 3 — The curve, the prior-width slice, certification

### 3a. The curve (natural skies, census priors, 200 populations)

Fraction right [95% interval], reference and Arm A (unknown Earth terms removed perfectly), timing ×1:

| N known | Reference: all pulsars | K ≥ 60 | K ≥ 600 | Right per sky | Arm A: all pulsars | K ≥ 60 | K ≥ 600 | Right per sky | A ÷ ref |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.009 [0.008, 0.011] | 0.002 | 0.000 | 0.7 | 0.010 [0.009, 0.012] | 0.001 | 0.001 | 0.8 | 1.12 |
| 2 | 0.012 [0.010, 0.014] | 0.001 | 0.001 | 1.0 | 0.017 [0.015, 0.020] | 0.003 | 0.000 | 1.4 | 1.41 |
| 3 | 0.017 [0.015, 0.019] | 0.003 | 0.001 | 1.4 | 0.023 [0.021, 0.026] | 0.005 | 0.001 | 1.8 | 1.34 |
| 5 | 0.024 [0.021, 0.026] | 0.005 | 0.001 | 1.9 | 0.035 [0.032, 0.038] | 0.011 | 0.002 | 2.8 | 1.48 |
| **7** | **0.029 [0.026, 0.032]** | 0.007 [0.005, 0.008] | 0.001 | 2.3 | 0.050 [0.046, 0.054] | 0.018 | 0.004 | 3.9 | 1.71 |
| **12** | **0.046 [0.043, 0.049]** | 0.016 [0.014, 0.018] | 0.003 | 3.6 | 0.083 [0.077, 0.089] | 0.041 | 0.012 | 6.6 | 1.81 |
| 20 | 0.063 [0.059, 0.066] | 0.027 | 0.008 | 5.0 | 0.119 [0.112, 0.127] | 0.068 | 0.020 | 9.4 | 1.91 |
| 30 | 0.080 [0.075, 0.085] | 0.039 | 0.012 | 6.3 | 0.158 [0.149, 0.167] | 0.101 | 0.034 | 12.5 | 1.97 |
| 60 | 0.109 [0.103, 0.115] | 0.062 | 0.019 | 8.6 | 0.223 [0.213, 0.234] | 0.164 | 0.062 | 17.7 | 2.06 |
| 100 | 0.132 [0.126, 0.139] | 0.080 | 0.027 | 10.5 | 0.275 [0.263, 0.287] | 0.211 | 0.087 | 21.7 | 2.08 |
| 300 | 0.190 [0.182, 0.198] | 0.128 | 0.047 | 15.0 | 0.386 [0.374, 0.399] | 0.322 | 0.164 | 30.5 | 2.03 |
| every candidate (1,300) | **0.260 [0.252, 0.269]** | 0.194 [0.186, 0.203] | 0.080 [0.072, 0.087] | 20.6 | **0.502 [0.490, 0.514]** | 0.440 [0.427, 0.453] | 0.249 [0.234, 0.262] | 39.6 | 1.93 |

Prior only: 0.006 / 0.001 / 0.000. Intervals for every cell are in `job53_logs/main.log`.

**Timing** (reference; all pulsars / K ≥ 60 / K ≥ 600):

| N known | ×1 | ×10 | timing noise ÷ 10⁶ |
|---|---|---|---|
| 7 | 0.029 / 0.007 / 0.001 | 0.029 / 0.006 / 0.001 | 0.029 / 0.006 / 0.001 |
| 12 | 0.046 / 0.016 / 0.003 | 0.046 / 0.016 / 0.002 | 0.046 / 0.016 / 0.002 |
| 30 | 0.080 / 0.039 / 0.012 | 0.079 / 0.038 / 0.013 | 0.079 / 0.039 / 0.013 |
| every candidate | 0.260 / 0.194 / 0.080 | 0.262 / 0.195 / 0.081 | 0.261 / 0.195 / 0.082 |

**Constructed skies, reference, ×1, for comparison** (all pulsars): N = 1: 0.011; 2: 0.019; 3: 0.026; 5: 0.033; 7: 0.042; 12: 0.053; 20: 0.068; 30: 0.080; 60: 0.117; 100: 0.149; 300: 0.206; every candidate 0.271. The two curves differ for N ≤ 12 (where the three raised sources are most of the known set) and merge from N = 30.

**Figure 3:** `job53_fig3_curve.png` (top: all pulsars; bottom: K ≥ 60; reference and Arm A with the region between shaded; timing levels as open triangles; prior-only dotted; grey band from 4.7 to 10.2 sources; constructed-sky curve thin dashed; right axis pulsars per sky).

### 3b. Fraction right against prior width (natural skies; every pulsar given the same σ_d)

True distances as in the census runs; the true offset of each pulsar from its prior centre is its census draw rescaled to the new σ_d (same random number). Reference:

| N known | 0.1 pc | 0.3 | 1 | 3 | 10 | 30 | 100 | 300 | crosses 0.9 / 0.5 / 0.1 at (pc) |
|---|---|---|---|---|---|---|---|---|---|
| prior only | 0.988 | 0.599 | 0.198 | 0.068 | 0.020 | 0.006 | 0.002 | 0.001 | 0.15 / 0.37 / 2.0 |
| 7 | 0.998 | 0.920 | 0.624 | 0.305 | 0.120 | 0.048 | 0.018 | 0.006 | 0.33 / 1.5 / 13.6 |
| **12** | 0.999 | **0.963** | **0.723** | **0.388** | **0.173** | 0.077 | 0.032 | 0.013 | **0.41 / 2.1 / 23** |
| 30 | 1.000 | 0.992 | 0.842 | 0.538 | 0.279 | 0.143 | 0.067 | 0.032 | 0.63 / 3.6 / 59 |
| every candidate | 1.000 | 1.000 | 1.000 | 0.951 | 0.725 | 0.501 | 0.301 | 0.193 | 3.9 / 30 / above 300 |

Arm A:

| N known | 0.1 pc | 0.3 | 1 | 3 | 10 | 30 | 100 | 300 | crosses 0.9 / 0.5 / 0.1 at (pc) |
|---|---|---|---|---|---|---|---|---|---|
| 7 | 1.000 | 0.966 | 0.720 | 0.408 | 0.187 | 0.089 | 0.039 | 0.017 | 0.41 / 2.2 / 27 |
| 12 | 1.000 | 0.987 | 0.808 | 0.514 | 0.282 | 0.158 | 0.075 | 0.038 | 0.54 / 3.2 / 70 |
| 30 | 1.000 | 0.997 | 0.902 | 0.674 | 0.447 | 0.293 | 0.172 | 0.102 | 1.0 / 7.6 / above 300 |
| every candidate | 1.000 | 1.000 | 1.000 | 0.989 | 0.905 | 0.775 | 0.607 | 0.468 | 10.4 / 233 / above 300 |

Intervals are ±0.01 or smaller everywhere (`job53_logs/main.log`). Median K of an N = 12 known set (coarsest tooth): 0.9 at 1 pc, 2.7 at 3 pc, 8.9 at 10 pc, 27 at 30 pc, 89 at 100 pc.

**Figure 4:** `job53_fig4_prior_width.png` (left reference, right Arm A; one line per N; prior-only dotted).

Reading: with the detectable sources known, the method is worth having only where the distance is already known to about a parsec. At the census 10th percentile (31 pc) N = 12 gives 8%; at the census median (242 pc) about 1.5%.

### 3c. Can anything be certified? (Job 52 Part 3 gate, unchanged: train on even skies, apply blind to odd skies)

Pulsars per sky (of 79) in the held-out catalogue, with its measured purity:

| Priors | Arm | N = 7 | N = 12 | every candidate |
|---|---|---|---|---|
| census | reference | **no 0.90 catalogue** (best 0.84) | 0.34 [0.22, 0.47] at 1.00; 97% of them K < 3 | 9.1 [8.5, 9.7] at 0.896; 42% K ≥ 60 |
| census | Arm A | 0.39 [0.26, 0.53] at 1.00; all K < 3 | 1.5 [1.3, 1.8] at 0.987; 75% K < 3 | 26.8 [25.3, 28.1] at 0.902; 65% K ≥ 60 |
| σ_d = 1 pc | reference | 25.4 [23.4, 27.3] at 0.905 | 45.7 [44.2, 47.2] at 0.904 | 79.0 at 1.000 |
| σ_d = 1 pc | Arm A | 44.6 [42.7, 46.2] at 0.902 | 62.0 [60.9, 63.0] at 0.899 | 79.0 at 1.000 |
| σ_d = 3 pc | reference | 3.6 [3.0, 4.2] at 0.906 | 10.5 [9.6, 11.3] at 0.920 | 79.0 at 0.950 |
| σ_d = 3 pc | Arm A | 9.2 [8.1, 10.1] at 0.912 | 20.1 [18.7, 21.4] at 0.896 | 79.0 at 0.990 |

(K in the census rows is from the coarsest tooth of the booked list; in the σ_d rows every pulsar has K below 3 by construction at 1 pc and a median 2.7 at 3 pc.)

Reading: with the detectable sources known and real priors, almost nothing can be certified, and what can is decided by its prior. With priors of 1 to 3 pc a useful certified catalogue exists (10 to 46 pulsars per sky at N = 12).

---

## Verdicts on the predictions

| | Prediction | Result | Verdict |
|---|---|---|---|
| V1 | Natural, N = 12: all-pulsar fraction below constructed (0.053), between 0.015 and 0.045; `arr_bg_opt` ≥ 3 for fewer than 8 per sky | 0.046 [0.043, 0.049], below 0.053 but just above the range; 10.2 per sky | **Half wrong.** Direction right, fraction 0.001 above the range; the detectable count is wrong (10.2, not below 8). `arr_bg` ≥ 3 does fall, 6.5 → 4.7 |
| V2 | Natural, every candidate within 0.05 of 0.271 | 0.260 | **Holds** |
| V3 | "known only" closer to "on"; "others only" closer to "off" | N = 12: 0.045 against on 0.046 and off 0.077; 0.085 against off 0.077 | **Holds**, and more strongly than predicted: all of the loss, not most |
| V4 | Pulsar terms from more than one bin higher below 30% of the unmodelled power (median) | 7.7% | **Holds** |
| V5 | N_s = 4: 90% contour crosses β = 0.8 at 2–6 pc, β = 0.5 at 0.3–0.8 pc, and lies beyond 30 pc at β = 0.95 | 2.64 pc; 0.34 pc; 23.9 pc | **Two of three hold.** At β = 0.95 the crossing is 24 pc (0.88 at 31 pc) |
| V6 | Flat with timing (below 0.01 at N ≤ 30); Arm A ÷ reference 1.5–2.2 for N ≥ 7 | largest change 0.0011; ratios 1.71 to 2.08 | **Holds** |
| V7 | N = 12, reference: above 0.5 only for σ_d ≤ 1 pc; below 0.1 for σ_d ≥ 10 pc | 0.72 at 1 pc, 0.39 at 3 pc (crosses 0.5 at 2.1 pc); 0.17 at 10 pc, 0.08 at 30 pc (crosses 0.1 at 23 pc) | **Half wrong.** The first half holds on the grid; the method keeps more than predicted at 10 pc |

## Limits

- "Known perfectly" grants frequency, sky position, amplitude and chirp mass of each known source with no error; Arm A removes unknown Earth terms with no error and no signal loss. Both are bounds.
- The fit treats bins as independent (adequate below 33 nHz; `JOB52_MAIN.md`).
- Each pulsar–source angle is drawn independently (`JOB52_REST.md`, Part 7), so the skies do not correspond to one real arrangement of pulsars.
- Heatmap sources have equal shares; real known sets do not.
- Constructed-sky comparisons rest on 20 skies.
