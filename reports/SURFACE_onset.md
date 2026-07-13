# SURFACE — the GENERAL onset surface, calibrated everywhere

**Agent:** SURFACE · ACCRE · tag `criterion-v2.1` (`git rev-parse HEAD` → `6bec3d6` ✓;
discovery `136b270f`, ee `f73b8e0`) · **Date:** 2026-07-12 · **PURE EXECUTION** (no commits,
no tracked-file edits).

**Scratch paths (host):** code `hpc_harbor/surface/` (`surface.py`, `surface_analysis.py`,
`surface_addendum.py`, `surface_figures.py`, 3 sbatch), results `SURFACE_results/` (24 840 lean
keyed npz + 4 summary npz + 3 figures; raw statistics always banked, **prior-referenced frame,
oracle columns beside and labelled**), logs `hpc_harbor/logs/sf_*`. Lane
`-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1`.

---

## 0. THE ANSWER IN ONE PARAGRAPH

**The onset exists. IGNITE-2's "no properly-calibrated onset exists anywhere in the modelled box"
was a statement about two cells, generalised to twenty-four — and when the D2 sizing is actually
paid at all 108 cells of the extended box, 59 of them clear onset.** The specific retraction
survives (IGNITE's `h*` = −12.75 lit / −13.25 vlbi remain **below** onset under honest floors:
0.87 and 0.47 correct certs/realisation here, reproducing IGNITE-2's 0.92 / 0.54 within the sky
error). What does **not** survive is the generalisation: STORY S6.1.5's *"nothing measured
contradicts the expectation that paying the sizing everywhere would close the rest of the box"* is
**refuted — paying it OPENS the box.** In the census structure at the census baseline the onset
sits at **h\* = −12.50** (1.13 [1.10, 1.23] correct certs/real, floor 106.04 ± 4.62 nat, N = 100) —
one grid step **louder** than IGNITE's retracted claim, not absent. And the axis IGNITE never
had is the one that matters most: **the population's loudness STRUCTURE is a stronger lever than
loudness itself.** Promoting two of the sixteen sources from faint to loud (3+13 → 5+11) raises
the count **6×** at fixed (h, T, tier) — *super-linearly*, and against a floor that itself nearly
doubles — while demoting one (3+13 → 2+14) all but extinguishes certification. In the 5+11
structure the map **ignites under `fALL` for the first time in the programme** (21 cells above
onset on the wrong-counterpart-robust column, best 2.57/real, against IGNITE's *"the map never
ignites under fALL, best 0.24"*). **Certification is a property of the POPULATION, not of the
source** — which is CHORUS's central premise, and it is now measured, not assumed. Two further
structural results: the certified count is **non-monotone in T**, with an h-dependent optimum
(T = 30 yr is optimal in **0 of 36** columns; loud cells peak at **40 yr** and *lose* at 50 as the
floor's data-volume growth overtakes the `T^{5/2}` leverage; faint cells keep gaining to 50) — so
*"T is the strongest lever"* needs its own ceiling attached; and **J1909-first is dead**: the
anchor **J0437-4715** leads at **19 of 22** near-onset cells in the census structure, and the
"J1909 first" of the census was itself a floor artifact (in the zero-noise re-cut below, J0437
leads 25/40 to 21/40 the moment the floor is removed). Both disputes are closed: **D-3** — the
zero-noise ceiling is **1.350 ± 0.82/draw** under the flat gate and `0.275` is retired as a
floor-concept category error; **D-4** — of the four uncalibrated cells, **two retract, one
confirms, one stays marginal.** The wrong-counterpart hole (D1) is untouched and travels beside
every count below.

---

## 1. GATES — both passed before the wide launch (STOP was not reached)

### g1 — reproduce IGNITE-2's two D2-sized cells bit-comparably (job **12507751**, ALL PASS, 681 s)

| gate | check | result |
|---|---|---|
| **S1** | banked `ig_sig_h1275_T30_lit_g0_n1240100` re-run through the SURFACE path | `dlnL`, `lnK`, `qmax`, `mapk`, `on_true`, `n_true_grid` — **all six columns max\|diff\| = 0.0** |
| **S2** | IGNITE-2's 540 banked nulls refit through SURFACE's own Gumbel estimator | **38.855 ± 1.465** (banked 38.86 ± 1.47) · **7.591 ± 0.481** (banked 7.59 ± 0.48); \|d\| ≤ 5e-4 |
| **S3** | IGNITE's banked signal realisations re-cut under those floors | **0.92** and **0.54** correct certs/real; **14/50** and **3/50** wrong — all four numbers exact |
| **S4** | structure axis at k_loud = 3 ≡ IGNITE | `cell_geometry_s` vs `IG.cell_geometry`: max\|dθ\| = 0.0, max\|dK\| = 0.0; `scramble_s(3\|16)` vs `F.scramble_source_sky(loud_only=True\|False)`: max\|diff\| = 0.0 |
| **S5** | T-extension is a strict no-op at dT = 0 | `lnL(truth\|zero-noise)` = **405413.512739**, bit-equal to the banked value |

S4 is what licenses the whole structure axis: the generalisation from IGNITE's hard-wired 3-loud
population to a `k_loud` axis is **bit-for-bit identical to IGNITE at k = 3**, so the census
column of this report is continuous with every banked number in the campaign.

### g2 — resume drill on the production path (jobs **12507894_0** → **12507979_0**)

Production array task `12507894_0` `scancel`'d **deliberately** at **6 m 51 s** mid-shard (41 of
its 828 realisations banked). Resubmission of the identical sbatch logged:

```
T=30: 8280 entries, shard 0/10 -> 828 mine; already banked: 41; to run: 787
```

and completed. **Skip-on-exist resume proven on the production SLURM path before the wide
launch**, per the standing convention. No banked artifact was touched.

---

## 2. WHAT WAS RUN, AND THE CONVENTIONS

**The grid — 108 cells.** h ∈ {−13.25, −13.0, −12.75, −12.5, −12.25, −12.0} × T ∈ {30, 40, 50} yr
× tier ∈ {lit, vlbi} × structure ∈ {3+13, 5+11, 2+14}. T = 15/20 are measured dead and were not
respent. **Per cell:** 5 GEO-protocol skies × 6 noise realisations = **30 Arm-B signal
realisations** (the honest arm is the only arm) + **100 nullN + 50 nullA + 50 nullL = 200 nulls**.
**24 840 realisations**, all banked.

**Floors — the campaign's credibility, and most of its array volume.** `floor(h, T, tier, struct)`
= **Gumbel-MLE (1 − α) quantile at α = 0.05** (z = 2.9702), fit error `2.80·β̂/√N` (G7's measured
constant). **Never max-of-N. Refit per cell, never inherited** — and the structure axis *forces*
this, because promoting a source changes the recovery model's amplitude and therefore the null.
`fN` (operative, D1) from N = 100 counterpart-matched nullN; **`fALL` (blind-robust column, D1)
from N = 200** scrambles+nulls, carried beside every number. **Every count in this report is
quoted three times** — at the floor, and at floor ± its own fit error — and an onset claim is
**CONFIRMED only if it survives floor + sd**. A count that dies when its floor moves by one sigma
is not an onset; it is a calibration artifact, and IGNITE's was one.

**Cadence extension (stated, per the mission).** Inherited **verbatim** from IGNITE §2 and
extended to T = 40/50: append dT = T − 15 yr of TOAs past each pulsar's own last TOA at its own
median cadence (median frequency, modal backend flag, median toaerr); Mmat columns extrapolated by
per-column least squares on a smooth basis; columns below R² = 0.99 — the binary-orbital class,
214/1141 — **zero-extended**, so their timing-fit leverage stays on the observed span; RN/GWB
Fourier counts scale with the span so the noise-model f_max is span-invariant; `fdot` leverage
∼ T^{5/2} rides free. Measured: T = 40 → +71 394 TOAs (span 47.17 yr, rn 30→64, gwb 14→30);
T = 50 → +99 955 TOAs (span 57.16 yr, rn 30→77, gwb 14→36). **This is a stated convention, not a
forecast of real future data** (no new backends, no DM events, no solar-wind realism), and at
T = 50 it is extrapolating 35 yr past the last real TOA — see §9.

**Loudness structure (the generality axis).** The 16-source population is frozen (seed 3000; sky
redrawn per GEO draw). Structure `k_loud` promotes the **first k_loud** sources to `log10_h = h`
and demotes the rest to −14.25, in **both** injection and recovery. The structures are therefore
**nested (2 ⊂ 3 ⊂ 5) on one frozen population** — a controlled comparison, not three different
universes. `dL`/`EV`/`K_counted` are strain-independent, so **the fringe grid and the trials
factor do not move with structure**; only the injected signal and the model amplitude do. k = 3 is
the census structure (gate S4).

---

## 3. DELIVERABLE 1 — THE SURFACE

Correct certifications per realisation. **Bold** = count at the fitted floor; `[lo, hi]` = count at
floor ± its fit error; `w` = wrong certs/real (**reported, never hidden**); `fA` = count on the
`fALL` blind-robust column. Floors below each block, `± 2.80·β̂/√N`, N = 100. Figure:
`SURFACE_results/surface_onset.png` (green box = onset surviving the floor's own calibration
error; green contour = count ≥ 1; blue contour = count ≥ 3).

### structure 3+13 (the census structure) — tier lit

| T | −13.25 | −13.00 | −12.75 | −12.50 | −12.25 | −12.00 |
|---|---|---|---|---|---|---|
| **30** | 0.67 [0.67,0.70] w0.03 | 0.37 [0.23,0.47] w0.03 | 0.87 [0.80,1.00] w0.23 | **1.13** [1.10,1.23] w0.23 | **1.20** [1.17,1.23] w0.40 | **1.13** [1.03,1.30] w0.27 |
| **40** | **1.37** [1.20,1.57] w0.10 | 0.47 [0.43,0.50] w0.03 | **1.20** [1.07,1.37] w0.33 | **2.33** [2.27,2.47] w0.50 | **2.47** [2.20,2.80] w1.13 | **2.73** [2.63,2.93] w0.70 |
| **50** | 1.10 [0.87,1.47] w0.00 | **1.13** [1.07,1.40] w0.00 | **1.33** [1.23,1.50] w0.17 | **1.57** [1.47,1.77] w0.33 | **1.83** [1.57,2.03] w0.20 | **1.83** [1.53,1.87] w0.43 |
| floor T=30 | 4.76±0.40 | 19.46±1.40 | 37.68±1.63 | 106.04±4.62 | 334.40±14.24 | 1064.59±45.07 |
| floor T=40 | 6.81±0.56 | 24.86±1.43 | 46.91±2.03 | 119.97±4.78 | 344.17±13.62 | 1088.04±48.94 |
| floor T=50 | 12.06±0.93 | 27.63±1.43 | 68.12±3.53 | 203.96±11.47 | 640.79±37.22 | 1977.69±110.44 |

### structure 3+13 — tier vlbi

| T | −13.25 | −13.00 | −12.75 | −12.50 | −12.25 | −12.00 |
|---|---|---|---|---|---|---|
| **30** | 0.47 [0.43,0.50] w0.03 | 0.47 [0.43,0.53] w0.03 | 0.73 [0.63,0.80] w0.03 | 1.17 [1.00,1.30] w0.27 † | **1.20** [1.17,1.43] w0.27 | 1.03 [0.97,1.27] w0.20 † |
| **40** | 0.63 [0.60,0.87] w0.03 | 1.00 [0.73,1.13] w0.00 | **1.63** [1.53,1.83] w0.27 | **2.50** [2.23,2.77] w0.57 | **2.70** [2.47,2.87] w0.70 | **2.70** [2.40,2.87] w0.80 |
| **50** | 0.83 [0.57,1.00] w0.00 | 0.97 [0.83,1.17] w0.00 | **1.70** [1.53,1.87] w0.30 | **1.30** [1.13,1.37] w0.37 | **1.43** [1.30,1.60] w0.27 | **1.37** [1.17,1.40] w0.47 |
| floor T=30 | 8.16±0.63 | 17.40±1.01 | 38.29±1.70 | 103.99±4.46 | 331.66±13.71 | 1074.03±45.49 |
| floor T=40 | 10.24±0.73 | 20.14±1.06 | 46.21±1.88 | 118.84±4.64 | 347.34±14.59 | 1064.78±44.94 |
| floor T=50 | 13.23±0.91 | 26.87±1.31 | 60.41±2.98 | 204.96±11.50 | 620.29±34.65 | 2006.73±113.36 |

† MARGINAL: count > 1 at the fitted floor but **not** at floor + sd. Not an onset.

### structure 5+11 — tier lit / tier vlbi

| T | −13.25 | −13.00 | −12.75 | −12.50 | −12.25 | −12.00 |
|---|---|---|---|---|---|---|
| **lit 30** | 0.60 [0.47,0.67] | **2.13** [1.93,2.50] | **2.47** [2.30,2.83] | **2.63** [2.33,2.90] | **3.03** [2.73,3.40] | **3.03** [2.53,3.53] |
| **lit 40** | **1.17** [1.03,1.47] | **2.87** [2.53,3.33] | **4.07** [3.63,4.30] | **6.53** [6.10,7.13] | **7.37** [7.07,7.87] | **7.80** [7.50,8.03] |
| **lit 50** | **1.30** [1.03,1.50] | **1.83** [1.63,2.03] | **2.50** [2.27,2.93] | **2.03** [1.83,2.57] | **2.60** [2.33,2.90] | **2.43** [2.17,2.60] |
| **vlbi 30** | **1.27** [1.03,1.47] | **2.33** [1.97,2.63] | **1.90** [1.70,2.27] | **2.70** [2.47,2.93] | **3.20** [2.90,3.60] | **2.77** [2.63,3.23] |
| **vlbi 40** | **1.63** [1.57,1.80] | **3.57** [3.33,3.87] | **4.07** [3.77,4.67] | **5.83** [5.50,6.23] | **7.07** [6.70,7.57] | **7.93** [7.73,8.47] |
| **vlbi 50** | **1.90** [1.57,2.10] | **2.40** [2.10,3.07] | **2.20** [1.97,2.63] | **2.50** [2.27,2.63] | **2.97** [2.60,3.40] | **2.30** [2.17,2.63] |

Floors, 5+11 lit: T=30 {20.17±1.41, 36.19±1.65, 105.80±5.08, 327.40±16.52, 962.30±46.41,
2985.07±136.63}; T=40 {23.41±1.49, 47.23±2.11, 129.07±5.29, 321.81±11.17, 965.94±31.65,
2829.50±80.24}; T=50 {29.38±1.70, 72.78±3.91, 218.26±11.96, 694.79±38.83, 2095.57±109.81,
6690.86±343.25}. (vlbi within ~10 % of these; full table in `surface_floors.npz`.)

### structure 2+14 — the axis in the other direction

| T | −13.25 | −13.00 | −12.75 | −12.50 | −12.25 | −12.00 |
|---|---|---|---|---|---|---|
| **lit 30** | 0.40 | 0.13 | 0.13 | 0.20 | 0.17 | 0.30 |
| **lit 40** | 0.63 | 0.20 | 0.23 | 0.40 | 0.33 | 0.43 |
| **lit 50** | **1.43** [1.37,1.57] | 0.23 | 0.47 | 0.73 | **1.03** [1.03,1.07] w0.63 | 1.03 [0.93,1.10] † |
| **vlbi (all T)** | ≤ 0.90 | ≤ 0.33 | ≤ 0.63 | ≤ 0.63 | ≤ 0.60 | ≤ 0.90 |

**Removing ONE loud source from the census population closes the box almost completely** — the
entire vlbi half of 2+14 posts **no onset at any (h, T) in the grid**.

### The tally

| | count |
|---|---|
| cells with corr > 1 at the fitted floor | **63 / 108** |
| **ONSET** (corr > 1 surviving floor + sd) | **59** |
| MARGINAL (corr > 1, dies at floor + sd) | 4 |
| **count ≥ 3 contour** (corr_lo > 3) | **9** (all 5+11, T = 40, h ≥ −13.0) |
| cells above onset on the **`fALL`** column | **21** (all 5+11) |

---

## 4. DELIVERABLE 2 — WHERE ONSET EXISTS: THE FRONTIER

The **faintest** h clearing onset (corr_lo > 1), per (structure, tier, T):

| structure | tier | T = 30 | T = 40 | T = 50 |
|---|---|---|---|---|
| 3+13 | lit | **−12.50** | **−13.25** ‡ | −13.00 |
| 3+13 | vlbi | −12.25 | −12.75 | −12.75 |
| 5+11 | lit | −13.00 | **−13.25** ‡ | **−13.25** ‡ |
| 5+11 | vlbi | **−13.25** ‡ | **−13.25** ‡ | **−13.25** ‡ |
| 2+14 | lit | none | none | **−13.25** ‡ |
| 2+14 | vlbi | **none — the extended box contains no calibrated onset** | none | none |

‡ **The frontier sits ON the faint edge of the grid in 7 of the 18 columns.** For those, `h*` is
**not bounded below by this box**: SURFACE locates an onset but cannot locate its faint boundary.
This is the honest limit of the result and the obvious next grid.

**The headline comparison, stated exactly.**

| claim | status |
|---|---|
| IGNITE: `h*` = −12.75 (30 yr, lit) | **STILL RETRACTED** — 0.87 [0.80, 1.00] here, below the bar |
| IGNITE: `h*` = −13.25 (30 yr, vlbi) | **STILL RETRACTED** — 0.47 [0.43, 0.50] |
| IGNITE-2 / STORY S6.1.2: *"NO properly-calibrated onset exists anywhere in the modelled box"* | **SUPERSEDED.** It was inferred from 2 calibrated cells. With all 108 calibrated: **59 onsets**, one of them (**−12.50, 30 yr, lit**) *inside IGNITE's own box* |
| STORY S6.1.5: *"nothing contradicts the expectation that paying the sizing everywhere would close the rest of the box"* | **REFUTED.** Paying it **opens** the box |

The onset was never absent. It sits **one grid step louder** than IGNITE claimed (h\* = −12.50, not
−12.75, in the census structure at 30 yr), and IGNITE's max-of-10 floors were drawn low enough to
*move it into the wrong cell* — which is a subtler and more instructive failure than "there is no
onset".

---

## 5. DELIVERABLE 4 — LOOP CELLS (reserved; the loop was NOT run)

The mission asks for the ≥2 cells most comfortably above onset, as the reservation for the
soft-loop above-onset test and PIPELINE. Two pairs, because *comfort* and *purity* do not pick the
same cells and the choice is a design decision, not mine to make:

**Pair A — maximum margin** (deepest above onset; the loop gets the richest seed set):

| cell | corr/real | margin above the bar at floor + sd | wrong/real | floor (N = 100) |
|---|---|---|---|---|
| **h = −12.00, T = 40, vlbi, 5+11** | **7.93** [7.73, 8.47] | **+6.73** | 0.50 | 2820.93 ± 80.37 |
| **h = −12.00, T = 40, lit, 5+11** | **7.80** [7.50, 8.03] | **+6.50** | 0.37 | 2829.50 ± 80.24 |

**Pair B — purity-preferred** (still ≥ 3× the bar, but with the cleanest seeds; the wrong-cert rate
is what poisoned every previous above-onset test — D3's reference and D4's protected population
both died of it):

| cell | corr/real | margin | wrong/real | floor (N = 100) |
|---|---|---|---|---|
| **h = −12.75, T = 40, vlbi, 5+11** | **4.07** [3.77, 4.67] | +2.77 | **0.13** | 122.46 ± 4.90 |
| **h = −13.00, T = 40, vlbi, 5+11** | **3.57** [3.33, 3.87] | +2.33 | **0.07** | 44.40 ± 1.94 |

**Recommendation, recorded not executed:** run the soft loop on **Pair B**. Pair A's cells carry
0.37–0.50 wrong certifications per realisation, and IGNITE-2 §3.2 showed that *every* soft-loop
failure was inherited from impure seeds, never generated by the loop. Pair B gives the loop a
genuinely above-onset seed set at **1/4 the impurity** — which is the first time in this programme
that such a cell has existed. **The loop was not run here, per the mission.**

---

## 6. DELIVERABLE 3 — WHO CERTIFIES FIRST: **J1909-first is dead**

Correct certifications pooled over all cells of a structure, and the leader at **near-onset** cells
(0.5 ≤ corr ≤ 1.6 — the cells actually *crossing* the bar, which is what the question asks):

| structure | pooled leader | J0437 rank | J1909 rank | leads at near-onset cells |
|---|---|---|---|---|
| **3+13** (census) | **J0437-4715** (331) | **1** | 4 (105) | **J0437 in 19 of 22**; J1713 2; J1909 **1** |
| **5+11** | J1713+0747 (337) | 3 (222) | 2 (278) | J1909 2, J1713 2 |
| **2+14** | **J0437-4715** (137) | **1** | 5 (36) | **J0437 in 5 of 10**; J0711 3; J1713 2 |

**J1909-first does not hold, in any structure, at the boundary where cells cross onset.** The
array's prior-pinned anchor **J0437-4715** — smallest `K_counted`, so the smallest trials bar —
arrives first wherever the floor is the binding constraint, which is exactly what onset means.
This confirms IGNITE §4.4's "above onset the anchor leads" and extends it across the structure
axis. J1909's primacy is a property of the *census loudness at zero noise*, and §8 shows it was
**itself a floor artifact**: the moment the floor is removed from the zero-noise scoring, J0437
overtakes it there too (25/40 vs 21/40). The one place J1909 leads is the 5+11 structure, where a
richer loud population restores the *data-driven* channel that J1909 was always the exemplar of —
which is a pleasing consistency, not an exception.

---

## 7. DELIVERABLE 7 — THE FOUR D-4 CELLS (STORY Appendix B): **uncalibrated → resolved**

Each rested on an under-sized 10-null floor (±2–18 nat) and was, per STORY, *"neither retracted nor
confirmed — UNCALIBRATED"*. N = 100 Gumbel refits resolve all four:

| cell | 10-null floor → corr | **N = 100 floor** | **→ corr [lo, hi]** | **verdict** |
|---|---|---|---|---|
| (−13.00, 30, lit) | 13.48 ± 3.21 → 1.10 | **19.46 ± 1.40** | **0.37** [0.23, 0.47] | **RETRACTED** |
| (−12.75, 30, vlbi) | 35.38 ± 5.03 → 1.24 | **38.29 ± 1.70** | **0.73** [0.63, 0.80] | **RETRACTED** |
| (−12.50, 30, lit) | 119.27 ± 18.30 → 1.06 | **106.04 ± 4.62** | **1.13** [1.10, 1.23] | **CONFIRMED** |
| (−12.50, 30, vlbi) | 104.93 ± 12.53 → 1.32 | **103.99 ± 4.46** | **1.17** [1.00, 1.30] | **MARGINAL** — dies at floor + sd |

**Two retract, one confirms, one stays marginal.** Note the mechanism is *not* uniformly "10-null
floors are biased low", which is what IGNITE-2's single datum suggested: at (−12.50, lit) the
properly-sized floor came out **11 % lower** than the 10-null value and the cell **survived**. A
max-of-N floor is not biased in a fixed direction — it is an order statistic with ±1.283β of
scatter and no fixed false-alarm rate, so it lands wherever its 10 draws put it. That is exactly
D2.2's argument, and this is the first four-cell test of it. **(−12.50, 30, lit) is the
programme's first confirmed onset cell.**

---

## 8. DELIVERABLE 6 — THE GEO ZERO-NOISE RE-CUT (STORY dispute D-3): **0.275 RETIRED**

CPU-only, from the banked 40 draws (`reports/geo_dlnl_bank.npz`); no new realisations.

**The adjudication.** At zero noise there are **no fluctuations for a noise-floor to gate**. GEO's
data are Asimov-at-truth; the 9.01-nat floor was fitted to *noisy* nulls. Applying it to a
noiseless statistic is not a mis-sized number, it is a **category error** — so the zero-noise
ceiling is quoted under the **flat gate alone** (criterion-v2.1 layers 1 + 3: `dlnL > ln K_counted`
and `q > 0.9`, no floor term).

| stage | per draw | sd | range |
|---|---|---|---|
| Bayes P_true > 0.9 (GEO's headline) | 4.50 | 1.47 | 1–9 |
| flat, layer 1 only (`dlnL > lnK`) | 1.38 | 0.89 | 0–4 |
| **v2.1 ZERO-NOISE CEILING (layers 1+3)** | **1.350** | **0.82** | **0–3** |
| strict (q > 0.99) | 1.025 | 0.82 | 0–3 |
| ~~criterion-v1 three-layer (floor 9.01)~~ | ~~0.275~~ | — | **RETIRED** |

**1.350 ± 0.82/draw** — the 1.35/draw class exactly, and it reproduces FORGE §9.2's independent
two-layer number (1.35 ± 0.82) to three digits. **Trail recorded: `0.275` is retired as a
floor-concept category error, not as a wrong measurement.** Its *sign* was always safe; its
*value* was never meaningful.

**And a caveat is closed, not inherited.** GEO's certification column is `P_true` — the posterior
at the **true** fringe, an **ORACLE** statistic under the criterion-v2.1 convention. Checked rather
than assumed: **0 of 4640 (draw, pulsar) cells have `dlnL < 0`**, i.e. at zero noise the true
fringe is the likelihood maximum in **every** cell, so the MAP fringe **is** the true fringe and
`q_max ≡ P_true` identically. **GEO's count is implementable, not an oracle** — the D4
oracle/implementable gap does not bite at zero noise. Per-pulsar at the adopted gate: **J0437-4715
25/40, J1909-3744 21/40**, J1713 4/40 — where the *retired* floor gave J1909 9/40 and J0437 1/40.
**The floor was manufacturing J1909-first.** (§6.)

---

## 9. THE STRUCTURAL RESULTS

### 9.1 The structure lever — certification is a property of the POPULATION

Correct certs/real vs `k_loud` at fixed (h, T, tier), with per-loud-source yield in brackets:

| cell | k = 2 | k = 3 | k = 5 |
|---|---|---|---|
| h = −13.00, T = 30, lit | 0.13 [0.07/src] | 0.37 [0.12/src] | **2.13 [0.43/src]** |
| h = −13.00, T = 40, lit | 0.20 [0.10/src] | 0.47 [0.16/src] | **2.87 [0.57/src]** |
| h = −12.50, T = 40, lit | 0.40 [0.20/src] | 2.33 [0.78/src] | **6.53 [1.31/src]** |

**The gain is super-linear in the number of loud sources, and it is measured against a floor that
rises with them** (the 5+11 floors are 2–3× the 3+13 floors at the same h, because the recovery
model carries more amplitude — the same `h^1.7` mechanism, now driven by *count* instead of
*strain*). Two loud sources is a dead population; five is a factory. **This is the
population-vs-source premise, and it is the first measurement of it in the repo.** It is *not* the
clock-sharing mechanism CHORUS pre-registered (that is about an eccentric comb pinning `f_gw`/`mc`
for the array): here the extra loud members are circular, and each simply adds fringe-breaking
evidence to every pulsar within its lag. But it establishes the premise CHORUS rests on —
**certification arithmetic is not a per-source property** — and it means **every single-source
no-go in this repo is scoped to a population structure nature does not have to supply.**

### 9.2 T is not monotone — the strongest lever has a ceiling

Which T maximises the count, over all 36 (h, tier, structure) columns:

| T | columns won |
|---|---|
| 30 yr | **0 / 36** |
| 40 yr | 18 / 36 |
| 50 yr | 18 / 36 |

And the split is **h-dependent**, not random:

| cell | T = 30 | T = 40 | T = 50 |
|---|---|---|---|
| 5+11, lit, h = −12.00 | 3.03 (floor 2985) | **7.80** (floor 2830) | 2.43 (floor **6691**) |
| 5+11, lit, h = −12.50 | 2.63 (floor 327) | **6.53** (floor 322) | 2.03 (floor **695**) |
| 3+13, lit, h = −13.00 | 0.37 (floor 19) | 0.47 (floor 25) | **1.13** (floor 28) |

**Loud cells peak at T = 40 yr and LOSE at 50; faint cells keep gaining out to 50.** The mechanism
is visible in the floors: the counterpart-matched floor grows with data volume as well as with
loudness (the matched-filter cross term integrates more data), and between 40 and 50 yr that growth
**overtakes the `T^{5/2}` fdot/coherence leverage** at loud h — while at faint h the floor is still
small enough that the leverage wins. IGNITE's *"onset is baseline-driven, not loudness-driven;
T^{5/2} beats the h^{1.66} floor race"* is **true up to a ceiling, and the ceiling is inside the
box**: past ~40 yr the floor race resumes and wins. **T = 30 yr is optimal nowhere.**

*Convention caveat, load-bearing:* the T = 50 cells extrapolate the timing model **35 yr** past the
last real TOA under the §2 convention. The *sign* of the T = 40 → 50 fall is robust where the floor race
bites — it appears in **12 of 12** loud columns (h ≥ −12.5) of the 3+13 and 5+11 structures, and
the floors' fit errors there are ~5 % — but in the 2+14 structure the count is small enough that
the leverage still wins and the count **rises** to 50 yr in all 6 loud columns. The *magnitude* is
a property of the extension convention and must not be quoted as a forecast.

### 9.3 The `fALL` column ignites — for the first time

IGNITE: *"under the fALL calibration the map never ignites anywhere (best cell 0.24/real)"*;
IGNITE-2: *"under fALL the surface never exceeds 0.36 anywhere"*. **SURFACE: 21 cells clear onset
on the fALL column, best 2.57/real** — and **all 21 are 5+11**. (3+13: 0 cells. 2+14: 0 cells.)

| cell | fALL corr | fALL floor (N = 200) | fN corr |
|---|---|---|---|
| h = −12.00, T = 40, lit, 5+11 | **2.57** | 5107.6 ± 176.5 | 7.80 |
| h = −12.00, T = 40, vlbi, 5+11 | **2.33** | 5314.9 ± 189.8 | 7.93 |
| h = −12.25, T = 40, vlbi, 5+11 | **1.93** | 1758.5 ± 64.6 | 7.07 |

**The D1 fork is no longer a window-closing choice — in a richer loud population, the
wrong-counterpart-robust criterion also turns on.** D1's price (*"the targeted programme's onset
exists because of the null family it is calibrated against"*) is real at the census structure and
**vanishes at 5+11**. This does not close the D1 wrong-counterpart hole (a floor tall enough to
outrun the noise-lock is not the same as a defence against a wrong counterpart, and D4 proved no
co-registration statistic can supply one) — but it removes the criterion's most uncomfortable
dependence on a design decision.

### 9.4 The loudness law, refit against a converging estimator

| | SURFACE (Gumbel, N ≥ 100, 18 independent per-(T,tier,struct) fits) | IGNITE (max-of-N) |
|---|---|---|
| `floor_fN` | **∝ h^1.72** (span 1.53–1.91) | h^1.66 |
| `floor_fALL` | **∝ h^1.80** (span 1.63–1.93) | h^1.88 |

**The h-law survives the estimator change intact.** It was never an artifact of max-of-N — only
its *normalisation* was, and that is what moved the onset cell.

---

## 10. DELIVERABLE 5 — THE JOIN: surface × ATLAS corner × SCOUT clock

Figure: `SURFACE_results/surface_join.png`. Strain convention verified against SIREN's banked DL
table (Mc = 10⁹, f_gw = 2×10⁻⁸ Hz at 7.77 Mpc → log₁₀h = −13.25).

Maximum `D_L` at which the first source can sit and still clear onset:

| f_orb (Hz) | log₁₀Mc | ATLAS min-e (self-clock > 20×) | D_L ≤ … at **h\* = −12.50** (census 3+13, 30 yr, lit) | D_L ≤ … at **h\* = −13.25** (5+11, T ≥ 40) |
|---|---|---|---|---|
| 10⁻⁸ | 9.5 | 0.59 | **9.4 Mpc** | **53.1 Mpc** |
| 10⁻⁸ | 9.0 | 0.58 | 1.4 Mpc | 7.8 Mpc |
| 10⁻⁸ | 8.5 | 0.70 | 0.2 Mpc | 1.1 Mpc |
| 10⁻⁸·⁵ | 9.5 | 0.66 | 4.4 Mpc | 24.6 Mpc |
| 10⁻⁸·⁵ | 9.0 | 0.77 | 0.6 Mpc | 3.6 Mpc |
| 10⁻⁸·⁵ | 8.5 | 0.84 | 0.1 Mpc | 0.5 Mpc |

**The join, read honestly, moves in BOTH directions at once.**

- **Against the programme:** in the census structure the onset is **louder** than IGNITE's
  retracted `h*` (−12.50 vs −12.75), so the (10⁹·⁵, 10⁻⁸) corner must sit inside **9.4 Mpc** —
  *inside* the Virgo distance, not at 17 Mpc. SCOUT's population clock (N̄ ≲ 0.01–0.1 at
  h ≈ −13.75) prices a source **1.25 dex quieter** than that. On the census population,
  **certification remains a variance play on the loud-nearby tail, and §9.2 caps the baseline
  lever at ~40 yr.**
- **For the programme, and this is the new thing:** the frontier is **not a property of the
  source alone**. At 5+11 the onset reaches **h\* ≤ −13.25 (unbounded below by this grid)**, which
  puts the same (10⁹·⁵, 10⁻⁸) corner at **53 Mpc** — *outside* Virgo, a ~180× larger volume than
  the census-structure onset allows. **The single most valuable unknown in the forecast is
  therefore no longer "how loud is the loudest source" but "how many loud sources are there".** The
  self-clocking corner (ATLAS: e ≳ 0.58, Mc ≳ 10⁹, top of band) is unchanged and still required for
  the *eccentric* road; what SURFACE adds is that the *circular* road's distance ceiling is set by
  a population parameter nobody has measured.

---

## 11. WHAT TRAVELS — CAVEATS, STATED

1. **The wrong-counterpart hole (D1) is untouched and OPEN.** Wrong certifications travel beside
   every count in §3 (up to **1.13/real** at (−12.25, 40, lit, 3+13); 0.37–0.53 at the Pair-A loop
   cells; 0.07–0.13 at Pair B). D3 and D4 are both rejected; no purity layer exists; and D4 proved
   the hole is **structurally un-closable by co-registration**. A count above onset is not a count
   of *correct* certifications you can defend against a mis-specified counterpart.
2. **30 realisations/cell** (5 skies × 6 weathers). GEO says the sky draw dominates yield variance;
   the quoted `± se` on each count (0.17–1.27) is the realisation scatter on the mean and is
   **large at the loud 5+11 cells** (±1.2 on 7.9). The *ordering* results (structure lever, T
   non-monotonicity, J0437-first) are consistent across all 36/18/22 columns respectively and do
   not rest on any single cell.
3. **The T = 50 extension convention** extrapolates 35 yr past the last real TOA (§9.2).
4. **The structures are nested on ONE frozen population** (seed 3000). The 5+11 result is a
   statement about promoting *these two* sources; a redraw of the population (different `f_gw`,
   `mc` for the promoted members) is untested, and the per-source yield would move with it.
5. **`h*` is not bounded below** in 7 of 18 frontier columns (§4). The faint edge of the onset
   surface is outside this grid.
6. **Nothing here touches a real TOA** — simulated array, simulated noise, injected sources,
   a stated cadence-extension convention.
7. **The 2 cells IGNITE-2 calibrated reproduce here** (0.87 vs 0.92; 0.47 vs 0.54, both within the
   ±0.2-class sky error), on **independent** signal seeds and **independent** N = 100 null banks.
   The agreement is a cross-check, not a re-use.

---

## 12. EXECUTION RECORD

| item | value |
|---|---|
| lane | `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1` |
| **g1 gate** | **12507751** — ALL PASS, 681 s (job wall 11 m 51 s) (S1 bit-identical, 6 columns max\|diff\| = 0.0; S2/S3 exact; S4 bit-for-bit; S5 no-op) |
| warm T=40 / T=50 | **12507752** (10 m 34 s), **12507753** (12 m 31 s) — extension diagnostics + steady-state timing |
| **g2 resume drill** | **12507894_0** `scancel`'d deliberately at 6 m 51 s (41 npz banked) → **12507979_0** logged `already banked: 41; to run: 787`, COMPLETED |
| production arrays | **12507894** (T=30, 10 shards), **12508061** (T=40, 10), **12508062** (T=50, 12) — **31/31 COMPLETED**, 16–22 min/task |
| batching | 828 (T=30/40) / 690 (T=50) realisations **per process**, one problem build per task — the ~230 s per-process materialisation is amortised ~800× |
| per-realisation wall (steady) | T=30 **0.92 s** · T=40 **1.04 s** · T=50 **1.40 s** |
| per-task fixed cost | build 208 s (T=30) / 359 s (T=40) / 463 s (T=50) + ~150–250 s first-call graph materialisation |
| realisations | **24 840** = 108 cells × (30 signal + 100 nullN + 50 nullA + 50 nullL) — all banked, none re-run |
| GPU cost | ≈ **11 GPU-hours** total (≈ 7 h science + ≈ 4 h per-process fixed cost) |
| disk | 24 840 lean npz + 4 summary npz + 3 figures; **1.8 GB on-disk** at panfs block granularity (~250 MB content) — 3 orders below the 3.1 TiB headroom |
| device log | every task logs GPU UUID + memory + foreign-process line as its first `[SURFACE]` line; **the squat lottery drew clean on every task** |

**Nothing was committed. Nothing was pushed. No tracked file was edited.** `SURFACE_onset.md`,
`SURFACE_results/`, `hpc_harbor/surface/` are untracked.

---

## 13. WHAT SURFACE SETTLES, AND WHAT IT OPENS

**Settled.**
1. **The onset exists** — 59 calibrated cells, including one inside IGNITE's own box. IGNITE's
   specific `h*` values remain retracted; the *generalisation* ("no onset anywhere") is superseded.
   The onset in the census structure is **h\* = −12.50 at T = 30 yr (lit)**, with a floor of
   **106.04 ± 4.62 nat (N = 100, α = 0.05)**.
2. **The floor is a converging function and its law is real** — `fN ∝ h^1.72`, `fALL ∝ h^1.80`,
   refit on 18 independent Gumbel fits. Only the *normalisation* was ever a max-of-N artifact.
3. **Certification is a property of the POPULATION.** The `k_loud` axis moves the count
   super-linearly and moves the *frontier* by ≥ 0.75 dex — a bigger lever than loudness, tier, or
   baseline. **`fALL` ignites at 5+11 and nowhere else.**
4. **T has a ceiling** (~40 yr at loud h), and **T = 30 yr is optimal nowhere.**
5. **J0437-first, not J1909-first**, at every onset boundary in the census structure — and
   J1909-first at zero noise was itself a floor artifact.
6. **D-3 closed** (zero-noise ceiling = 1.350 ± 0.82/draw; 0.275 retired, trail recorded; GEO's
   column proven implementable, not oracular). **D-4 closed** (2 retracted, 1 confirmed, 1
   marginal).

**Open.**
1. **The faint frontier.** `h*` is unbounded below in 7 of 18 columns. The next grid is fainter,
   not louder.
2. **The population axis proper.** SURFACE moved `k_loud` on a frozen 16-source population. The
   real question — the mix, the eccentric minority, the clock-sharing test — is **CHORUS's**, and
   SURFACE has now established the premise it rests on.
3. **The loop above onset.** Pair B (§5) is the first genuinely above-onset, low-impurity seed set
   the programme has ever had. **Not run here.**
4. **The purity hole travels unchanged.** Every count above onset carries wrong certifications, and
   there is no layer left to remove them.
