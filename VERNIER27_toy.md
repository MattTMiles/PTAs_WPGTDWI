# VERNIER-27 — firm the K\* number; and the number turns out to be two numbers

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, **243 min** — the
certification pass was the point. Clean room: imports only `vernier27_lib` / `vernier25` /
`vernier24` (and through them `vernier_toy` / `vernier6` / `vernier7` / `vernier16`). None of
those files is edited.

Code `hpc_harbor/vernier/vernier27.py` + `vernier27_lib.py` · data `reports/vernier27.npz` ·
log `CW_transition/logs/vernier27.log`. Staged, not committed.

**2 arms × 6 K × 20 realisations × 6–8 chains = 1 680 chains**, plus a 480-chain comb-cap arm
and a 36-chain safety kit.

---

## Lead — K\*\_50 = 13–15 and K\*\_90 = 16–19, but only if the sources are in different bins

| | K\*\_50 | K\*\_90 |
|---|---|---|
| **Arm W** (one source per 1/T bin) — contrast | **13.1 ± 0.6**  [12.0, 14.0] | **16.0 ± 1.5**  [15.3, 21.3] |
| **Arm W** — certification | **14.0 ± 0.7**  [12.7, 15.3] | **19.0 ± 1.3**  [15.8, 21.3] |
| **Arm P** (all sources inside 10–25 nHz) | **does not exist** | **does not exist** |

Errors are cluster bootstrap over realisations (4 000 resamples), interpolated between grid
points; brackets are 95 %.

**The deliverable is not one number with an error bar, and it is not a factor-of-two range.
It is a fork.** The variable that decides whether the wheel ever turns is **sources per
frequency bin**, not K. Spread K sources one per bin and the mechanism works, comfortably,
by K ≈ 14–19. Pack the same K sources into the ten 1/T bins of 10–25 nHz and **contrast peaks
at K = 16 and then falls**; the crossing fraction plateaus near 0.4 and declines; there is no
K at which it certifies reliably, out to K = 32.

**Slope verdict: no saturation in Arm W.** Contrast is still climbing at **+0.8 to +1.1 nat
per comb** between K = 28 and K = 32, on the same ~1–3 nat/comb curve it started on, and by
K = 32 it stands at **29.1 nat against a bar of 6.2 — 4.7× the bar.** Nothing in this arm needs
extrapolating. **In Arm P the slope is not small, it is undefined**: the accumulated comb count
is pinned at 7–8 at every K from 12 to 32, so there is no d(contrast)/d(comb) to read. What
moves there is d(contrast)/d**K**, and beyond K = 16 it is **negative**.

**FLAG, as the brief asked: yes, things cross below K = 16, in both arms.** K\*\_50 lands at
13–15. Before that is used to move the paper from SKA-era to MeerKAT-era, read
[§ What this does not settle](#what-this-does-not-settle) item 1 — the certification gate is
**not null-calibrated**, so every onset threshold here is a **lower** bound on the true K\*.

---

## The V25 defect this run was really about

**VERNIER-25's Arm K could not place more than ten seeds.** `vernier25.chain()` seeds bin `b`
iff `b < n_seed`, `vernier25.build_bin()` prepends exactly **one** seed, and `vernier25.BAND`
holds **ten** 1/T bins:

| K | 2 | 3 | 5 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|
| **seeds actually placed** | 2 | 3 | 5 | 8 | **10** | **10** |

So V25's K = 12 and K = 16 rows are **the same cell**, differing only by chain seed. The
contrast-vs-combs curve that produced **K\* ≈ 20–25 was extrapolated off a tail that had
stopped moving for a construction reason.** Every V25 result at K ≤ 10 is unaffected, and
Arms G and G-FP (no seeds) are untouched.

**Fixing it has no unique answer, and that is the finding.** At T = 20 yr the band 10–25 nHz
holds exactly ten 1/T bins, so K = 32 sources at distinct frequencies inside that band **must**
share bins. Two constructions, both defensible, run as co-primaries:

- **Arm P — PACKED.** The brief read literally: K seeds spread over the ten bins of
  10–25 nHz, ⌈K/10⌉ per bin, each on its own sub-bin frequency (regular grid, ±0.35/T about
  the bin centre; sub-bin frequency is resolvable at these SNRs, σ_f/Δf = 1/2πρ = 1.3 % at
  ρ = 12, and it is exactly the quantity the comb spacing dL = c/(f(1−cos μ)) uses).
- **Arm W — WIDE.** V25 Arm K's own construction continued: one seed per bin, band widened to
  K bins and centred on the sensitive decade (K = 32 → 3.0–52.1 nHz).

Arm P is the crowded, pessimistic case; Arm W the spread, optimistic one. They are reported as
a bracket and never averaged.

## Three other changes, each of which moved the answer

3. **`COMB_CAP` 8 → 32.** V25's own caveat 1 called this "the single most leveraged assumption
   in the note". It was right — see [Arm C](#arm-c--what-v25s-comb_cap--8-cost).
4. **The GENUINE pool made K-invariant** — seeds + top-10 sky, not a global top-10 that the
   seeds themselves crowd out as K rises. Chance-match floor rises 0.060 → 0.076 with the
   enlarged pool and is reported.
5. **The detection-floor RNG isolated from the data RNG.** V25 drew its recut floors from the
   same stream as the data, and `floor_for` is memoised — a cache *hit* consumes no randomness
   and a *miss* does, so a V25 chain's data depended on how many chains had run before it.
   Chains here are exactly reproducible from (`seed`, `noise_seed`). This is a latent V25
   reproducibility hazard, not a design choice.

Also `max_turn` 16 → 2K + 16. Checked: chains stop on their own well inside it.

---

## Arm P — packed into 10–25 nHz. 20 realisations × 8 chains per K

| K | A_eq | K_rec | combs | S_end | ±SE | IQR | bar | S/bar | cross | ±SE | CERT/run | true | WRONG | null |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 | 1.71 | 7.0 | 7.0 | 1.88 | 0.29 | 0.94–3.72 | 5.29 | 0.36 | 0.50 | 0.11 | 0.07 | 0.08 | 0.01 | 0.00 |
| **16** | 1.59 | 8.0 | 8.0 | **2.54** | 0.49 | 1.20–5.16 | 5.29 | **0.48** | 0.55 | 0.11 | 0.18 | 0.11 | 0.09 | 0.01 |
| 20 | 1.83 | 8.0 | 8.0 | 1.76 | 0.34 | 0.87–3.19 | 5.26 | 0.34 | 0.45 | 0.11 | 0.18 | 0.18 | 0.02 | 0.00 |
| 24 | 1.77 | 8.0 | 8.0 | 1.46 | 0.35 | 0.64–3.15 | 5.26 | 0.28 | 0.40 | 0.11 | 0.09 | 0.11 | 0.02 | 0.01 |
| 28 | 1.87 | 7.0 | 7.0 | 1.38 | 0.35 | 0.49–2.66 | 5.21 | 0.27 | 0.30 | 0.10 | 0.11 | 0.13 | 0.00 | 0.00 |
| 32 | 1.49 | 7.0 | 7.0 | 1.15 | 0.23 | 0.50–2.53 | 5.26 | 0.22 | 0.35 | 0.10 | 0.07 | 0.09 | 0.01 | 0.00 |

`cross` = fraction of realisations whose **best chain** reaches its own trials bar.

**Recovered multiplicity is capped at ~8 no matter how many sources exist.** `K_rec` reads
7, 8, 8, 8, 7, 7 while K runs 12 → 32. Seed recovery efficiency collapses from 0.50 to 0.19
and the blend fraction climbs 0.59 → 0.75. **The extra sources are not merely wasted — they
are actively harmful**, because every unrecovered seed is a loud interloper in the
λ_eff-consistent confusion budget of the ones the chain is trying to work with. That is why
contrast *declines*: dS/dK is +0.17 from K = 12 → 16, then −0.19, −0.08, −0.02, −0.06.

**Arm P's "K\*\_50 = 12.0 ± 2.7" is an artefact and must not be quoted.** The crossing fraction
sits at 0.50 at the first grid point purely by scatter and then goes *down*; the interpolator
reports the first K reaching 0.5 without knowing the curve is not monotone. **The correct
statement is that Arm P has no K\*: crossing plateaus at 0.4 ± 0.1 and never improves.**

## Arm W — one source per 1/T bin. 20 realisations × 6 chains per K

| K | bins | A_eq | K_rec | combs | S_end | ±SE | IQR | bar | S/bar | cross | ±SE | CERT/run | true | WRONG | null |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 12 | 8.3–25.7 nHz | 1.32 | 7.5 | 7.5 | 2.87 | 0.33 | 1.5–3.8 | 5.32 | 0.54 | 0.35 | 0.11 | 0.07 | 0.07 | 0.03 | 0.01 |
| **16** | 5.1–28.9 | 0.98 | 11.0 | 11.0 | **6.17** | 0.96 | 4.0–9.9 | 5.52 | **1.12** | **0.90** | 0.07 | 0.57 | 0.77 | 0.06 | 0.03 |
| 20 | 3.0–33.1 | 1.04 | 13.0 | 13.0 | 11.05 | 1.38 | 6.7–16.9 | 5.66 | 1.95 | 0.95 | 0.05 | 1.53 | 2.02 | 0.32 | 0.16 |
| 24 | 3.0–39.4 | 1.18 | 16.5 | 16.5 | 21.27 | 2.65 | 10.9–29.7 | 5.87 | 3.62 | 1.00 | 0.00 | 3.42 | 7.65 | 0.59 | 0.39 |
| 28 | 3.0–45.8 | 1.22 | 21.0 | 21.0 | 24.91 | 3.91 | 15.2–40.1 | 6.04 | 4.13 | 1.00 | 0.00 | 4.07 | 7.17 | 0.82 | 0.46 |
| 32 | 3.0–52.1 | 1.50 | 25.0 | 25.0 | 29.12 | 5.22 | 20.1–51.8 | 6.20 | 4.70 | 1.00 | 0.00 | 5.55 | 10.36 | 1.12 | 0.68 |

**Recovery efficiency *rises* with K here — 0.50, 0.53, 0.55, 0.58, 0.64, 0.66 seeds hit per
seed offered — and the blend fraction *falls*, 0.62 → 0.49.** The exact opposite of Arm P, from
the same sky, the same array, the same search and the same certification statistic. One source
per bin is the whole difference.

**Ignition, for the first time in the series.** V25 reported an ignition fraction of 0.00 in
all eighteen of its cells. Arm W ignites (two consecutive deepening turns, scrambled-null
clean by construction) in **0.03 → 0.06 → 0.28 → 0.27 → 0.39** of chains across K = 16 → 32,
with deepening turns per chain rising 0.07 → 4.22. **The detection channel does acquire positive
feedback — it just needs ~11 recovered combs before it does.** Arm P never gets there: its
ignition fraction is 0.00–0.01 at every K.

### Crossing fraction, the curve the deliverable is read off

| K | 12 | 16 | 20 | 24 | 28 | 32 |
|---|---|---|---|---|---|---|
| Arm W, contrast | 0.35±0.11 | **0.90±0.07** | 0.95±0.05 | 1.00 | 1.00 | 1.00 |
| Arm W, certification | 0.25±0.10 | 0.75±0.10 | **0.95±0.05** | 1.00 | 1.00 | 1.00 |
| Arm P, contrast | 0.50±0.11 | 0.55±0.11 | 0.45±0.11 | 0.40±0.11 | 0.30±0.10 | 0.35±0.10 |
| Arm P, certification | 0.30±0.10 | 0.55±0.11 | 0.45±0.11 | 0.40±0.11 | 0.30±0.10 | 0.35±0.10 |

### Shape — contrast vs accumulated combs

| K | combs | S_end | ±SE | nat/comb (local) | nat/comb (within-K) | gap to bar |
|---|---|---|---|---|---|---|
| 12 | 7.5 | 2.87 | 0.33 | 0.94 | 0.72 | +2.45 |
| 16 | 11.0 | 6.17 | 0.96 | 2.44 | 0.82 | −0.65 |
| 20 | 13.0 | 11.05 | 1.38 | 2.92 | 0.95 | −5.39 |
| 24 | 16.5 | 21.27 | 2.65 | 0.81 | 1.81 | −15.40 |
| 28 | 21.0 | 24.91 | 3.91 | 1.05 | 2.07 | −18.88 |
| 32 | 25.0 | 29.12 | 5.22 | — | 0.82 | −22.92 |

**Not saturating.** The local slope averages +1.22 nat/comb over K = 12–16 and +0.99 over
K = 28–32 — statistically the same, and both inside V25's ~2 nat/comb read. The within-K
regression of contrast on comb count is flat-to-rising (0.72 → 2.07 → 0.82). **The SKA-era
claim is not killed by curve shape; if anything the curve is more generous than V25's, because
it crosses the bar at ~10–11 combs exactly as V25 predicted and then keeps going.**

## Arm C — what V25's `COMB_CAP` = 8 cost

Same universes, same noise, cap 8 vs cap 32. 10 realisations × 4 chains.

| arm | K | S_end cap 8 | S_end cap 32 | ratio | cross cap 8 | cross cap 32 | CERT cap 8 | CERT cap 32 |
|---|---|---|---|---|---|---|---|---|
| P | 16 | 2.23 | 2.14 | 0.96 | 0.30 | 0.30 | 0.12 | 0.12 |
| P | 24 | 2.11 | 2.25 | 1.06 | 0.10 | 0.20 | 0.03 | 0.03 |
| P | 32 | 1.51 | 1.26 | 0.83 | 0.30 | 0.30 | 0.07 | 0.10 |
| W | 16 | 3.81 | 5.11 | 1.34 | 0.40 | 0.60 | 0.07 | 0.45 |
| **W** | **24** | **9.47** | **21.53** | **2.27** | 0.80 | 1.00 | 0.90 | **4.65** |
| **W** | **32** | **9.47** | **22.30** | **2.36** | 0.90 | 1.00 | 1.12 | **6.53** |

**V25's caveat 1 was correct and the cap was doing real damage.** Look at the cap-8 column in
Arm W: **9.47 at K = 24 and 9.465 at K = 32 — the same number.** That is the cap, not the sky.
Certifications per run are suppressed **5–6×**. Arm P is untouched because it never accumulates
8 combs in the first place, which is why V25 — whose Arm K was effectively packed at ten
seeds — saw so little.

## The safety kit

- **M = 1 control** (single bin, seeded, sea removed, 12 realisations): detected 6/12,
  **0 certifications**, ρ_P max 3.06, fringe contrast **0.000 nat**. One comb cannot break its
  own degeneracy at any SNR — V24 and V25 reproduced.
- **Pure-noise null through the whole chain**, K = 32 geometry, no signal, 12 realisations:
  Arm P fired **3/12 = 0.25** against an expected 0.40 for 10 bins at α = 0.05; Arm W fired
  **11/12 = 0.92** against an expected 0.81 for 32 bins. Both consistent. **0 certifications
  and 0 ignitions in both.**
- **Per-turn source-level scrambled null.** Null certifications per genuine+wrong
  certification: **0.00–0.10 in Arm P, 0.03–0.08 in Arm W.** A scrambled-null certification
  vetoes the deepening flag by construction, so no ignition reported above is scrambled-dirty.
- **Wrong-certification ledger.** Wrong fraction is **0.07–0.14 in Arm W at every K**, and
  **0.44 at Arm P K = 16** — that single cell is the worst number in this note. It sits on
  0.11 true and 0.09 wrong per run, so it is 18 wrong certifications out of 41 across 160
  chains, not a small-count mirage. It is also the cell Arm P's spurious "K\*\_50 = 12" is
  interpolating through.
- **V18 accounting.** Splits are **0.0 everywhere** — the chain never re-finds a source it
  already subtracted. Blends run 0.49–0.77 of turns.
- **Chance-match floors** 0.060 (Arm W, 1 seed + top-10 sky) to 0.076 (Arm P at K = 32,
  4 seeds + top-10).

### The sampling requirement, honestly

The brief asked for enough chains that SE(median S_end) < half the gap between adjacent K
medians. **Met in Arm W at K = 12, 16, 20 — the three cells that actually set K\*** (0.33 vs
1.65; 0.96 vs 1.65; 1.38 vs 2.44). **Not met at Arm W K ≥ 24**, where the contrast distribution
develops a heavy tail (IQR 20–52 nat at K = 32) — but crossing there is already 1.00 ± 0.00, so
the criterion is moot. **Not met anywhere in Arm P**, where the K-to-K gaps are 0.04–0.33 nat
against SE ≈ 0.23–0.49. **That is not a sampling failure, it is the result**: in Arm P the
K-dependence of contrast is smaller than the realisation scatter, which is what "no K\*" means
quantitatively. Meeting the literal criterion there would need ~600 realisations per cell to
resolve gaps the mechanism is not producing.

---

## What this does not settle

1. **The certification gate is not null-calibrated, so every K\* here is a LOWER bound.**
   `vernier24.certify` sets `bar = max(ln K_i, floor_i + err_i)`. Measured directly in this
   run's geometry: **ln K_i = 5.78 nat against an incoherent-comb floor of 0.003–0.45 nat at
   every comb count from 1 to 19.** The floor term never binds, so CRITERION-v2.2's `recut_floor`
   and its α = 0.05 are dead code and the gate is a bare trials term whose realised false-alarm
   rate is calibrated to nothing. Contamination bites hardest where the signal rate is lowest —
   i.e. at the onset — so **K\*\_50 = 13–15 and K\*\_90 = 16–19 should be read as lower bounds,
   and the sub-16 crossing flag should NOT on its own move the paper's framing to MeerKAT-era.**
   Fixing this means calibrating `recut_floor` at α/K_i, which moves banked columns in five
   instalments and is not done here.
2. **Arm W leaves the brief's band.** At K = 32 it runs 3.0–52.1 nHz. The low bins carry far more
   confusion and the high bins far less; both effects are in the numbers but they are not the
   10–25 nHz the brief specified. Arm P is the arm that honours the band.
3. **Sub-bin frequencies cost Arm P no detection SNR**, because the toy's F-statistic searches
   sky position at a fixed bin frequency. A real search would pay a mismatch penalty, so Arm P
   as built is if anything *optimistic* about packed sources — and it still turns over.
4. **All seeds are at ρ = 12 and A_eq ≈ 1.0–1.9.** A realistic K = 32 population is a luminosity
   function, not 32 copies of one source, and the loudest few would dominate the comb set.
   Untested here.
5. **Detection is Earth-term only** (V16/V22's standing limitation), so the blend rate that
   caps `K_rec` at 8 in Arm P is measured with the weaker search.
6. **`K_rec` and `combs` are the same column by construction** (the comb set is the genuine
   detections). Arm P's cap at 8 combs is therefore a cap on recovery, not on the certifier.
7. **20 realisations bounds a measured 0.00 at ≲ 14 % (95 %)**; the 1.00 cells in Arm W bound
   the failure rate at ≲ 14 %, not at zero.
8. **`S_end` is max-over-turns of a per-pulsar max** — the most optimistic honest summary of the
   certification channel, inherited from V25 §6 so the columns are comparable.
9. **Arm P's `K*_50` values in the banked npz (12.0 and 15.2) are interpolator artefacts on a
   non-monotone curve** and are contradicted in this note's text. The npz keys are kept for
   provenance; do not read them without this paragraph.

## One-line verdict

**V25's K\* ≈ 20–25 is superseded and its Arm K was broken — `BAND` holds ten bins and
`chain()` placed one seed per bin, so K = 12 and K = 16 were the same cell and the
contrast-vs-combs extrapolation ran off a dead tail (K ≤ 10 unaffected) — and correcting it
does not give one number but a fork, because at T = 20 yr the band 10–25 nHz holds exactly ten
1/T bins and K > 10 distinct frequencies inside it must share them: spread one source per bin
(Arm W) and the wheel works, K\*\_50 = 13.1 ± 0.6 on contrast and 14.0 ± 0.7 on certification,
K\*\_90 = 16.0 ± 1.5 and 19.0 ± 1.3, with contrast still climbing at ~1 nat/comb at K = 32 where
it stands 4.7× above the bar, recovery efficiency *rising* 0.50 → 0.66, blend fraction *falling*
0.62 → 0.49, and the series' first non-zero ignition (0.39 of chains at K = 32) — while packing
the same K sources into 10–25 nHz (Arm P) caps recovered multiplicity at 7–8 combs at *every* K
from 12 to 32, turns contrast over after K = 16 (dS/dK negative thereafter), drives seed
recovery from 0.50 to 0.19 and the blend fraction from 0.59 to 0.75, and leaves the crossing
fraction plateaued at 0.4 ± 0.1 with no K\* at all; V25's own caveat-1 comb cap is confirmed as
the second big lever (cap 8 pins Arm W's contrast at 9.47 nat at both K = 24 and K = 32 and
suppresses certifications 5–6×, while doing nothing in Arm P, which never reaches 8 combs —
explaining why V25 saw so little); the safety kit is clean (M = 1 gives 0.000 nat and zero
certifications, the pure-noise null gives zero certifications and zero ignitions in both arms
at the K = 32 geometry, scrambled-null fractions run 0.00–0.10, splits are zero everywhere)
with one bad cell, Arm P K = 16, where 18 of 41 certifications are wrong; **but the whole
certification gate is `max(ln K_i, floor_i)` with ln K_i = 5.78 nat against a floor of
0.003–0.45 nat, so the α = 0.05 recut floor never binds, the false-alarm rate is calibrated to
nothing, and every K\* above — including the sub-16 crossings that would otherwise move the
claim from SKA-era to MeerKAT-era — must be read as a LOWER bound until that gate is
recalibrated at α/K_i.**
