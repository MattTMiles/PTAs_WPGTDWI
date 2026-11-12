# BASELINE — the field-baseline comparison the paper lacks

**Agent BASELINE, ACCRE, 2026-07-29.** Working report. **REPORT-ONLY: nothing here arms a
protocol step, moves a banked verdict, or enters a closure claim.** I stage; Matt commits.

---

## HEADER — the campaign's declared conditions

| | |
|---|---|
| **Spec** | `CW_transition/trackB_estimator_spec.md` §CRITERION-v2.1, extended to **criterion-v2.2** (spec text authoritative; tags `criterion-v2`, `v2.1`, `v2.2` all exist) |
| **Floor convention** | `reports/RECUT_floors.md` §0 — Gumbel MLE (1−α) quantile at α = 0.05, fit error `2.80·β/√N`; **above a 20 % nullN zero-fraction the Gumbel is INVALID** and the floor is the empirical q95 with a bootstrap error. **The zero-fraction is a required column, not a caveat.** |
| **F-stat machinery** | `hpc_harbor/spark/spark.py` — `build_detectors` (`earth1`, EarthDelay), `make_fstat_earth` (**gate g0a: == `TE.make_fstat` verbatim**), `scan_grid`, `flat_grid`, `ll0_earth` |
| **Criterion recap** | STORY §S5 — DETECTION `dlnL_a > max(ln K_a, floor)`; CERTIFICATION `q_max,a > 0.9` **within detections only**; PURITY **NONE** (D3 and D4 both tested, both rejected) |
| **Venue** | T = 30 yr cadence-extension convention, 116-pulsar MOCK array (AXIS, 1440 MHz), `lit` distance tier |
| **Lane** | `interactive_gpu` **A100-40GB on dgx01**, `-A dsi_dgx_iacc -q dgx_iacc`, `--cpus-per-task=8`, array cap `%4`. See the lane claim in `GLACIER_capstone.md` (2026-07-29). **PHOENIX (`frzgate`) holds the A100-80GB share on dgx03; BASELINE does not enter it.** `batch_gpu` H200 not entered (its `--test-only` start date is 2026-11-04). |
| **Budget class** | few GPU-hr. Header estimate and the measured spend are in §6. |
| **Banks** | `BASELINE_results/` only. `SURFACE_results/`, `GENERALISE_results/`, `GLACIER*/`, `SPARK*/` are **read-only inputs**; nothing in them is written or moved. |

### Scope lines that ride with every number below

1. **The array is a MOCK.** 116 real pulsar sky positions, real TOA uncertainties, real
   published distance priors — but the residuals **are** the injected CW plus drawn noise.
   No real TOA is touched. (`data-spine-is-simulated`; SPARK "SCOPE OF INFERENCE".)
2. **T = 30 is a convention, not a forecast.** Cadence extension per IGNITE §2 / SURFACE:
   no new backends, no DM events, no solar-wind realism.
3. **These cells are FAR above NG15.** The frozen 16-source POP with 3 (or 5) members at
   `log10 h = −12.75` sits well above the +1.38 dex that the GLACIER-LITE addendum already
   measured for the `−13.25` version of the same population. Per that addendum's binding
   labelling rule, **the phrase "NG15-consistent" is struck**; cells are labelled by their
   loudness parameters only.
4. **The F_e grid is the incumbent seeder's** (healpix `nside = 4`, ~14.7°, × `1/(3T)` in
   frequency), **followed by a full 8-parameter polish** (§1.3–1.4). A finer grid raises both
   the signal statistic and the null threshold.
5. **`q_max > 0.9` is the criterion's own certification bar.** Applying it to the field's
   posterior is the whole point of B2 — it is not a bar the field would state, and B2 is
   explicit that the field's confidence number is *the same kind of object* (posterior mass
   on the MAP fringe), computed on a different comb.

---

## 0. PRE-REGISTRATION — written and staged BEFORE the fan was submitted

Everything in this section is fixed in advance. Section 4 reports against it and nothing
else.

### 0.1 The two cells and the null set

| id | cell | source campaign | realisations |
|---|---|---|---|
| **C1** | census `3+13`, `h = −12.75`, `T = 30`, `lit` | SURFACE (`sf_sig_h1275_T30_lit_k3_*`) — the D2-sized `G1_CELLS` cell, and D4's onset cell | 5 skies × 6 noise = **30** |
| **C2** | A-SKY survivor: `5+11`, `e = 0.3`, `n_ecc = 1`, `h = −12.75`, `T = 30`, `lit` | GENERALISE Arm A-SKY (`gen_sig_AS_e03_h1275_k5_s*`) — **the one boundary cell of four that survives sky marginalisation** (6/8 skies, sky median 1.77) | 8 skies × 15 = **120** |
| **N** | pure noise, no CW, `T = 30` | SURFACE `sf_nullN_h1275_T30_lit_k3_*` seeds | **100** |

**One null set, and it is licensed by a gate, not by assumption.** A `nullN` realisation has
`clean = 0`: its data is `NoiseDrawer.draw(seed)` and nothing else. The noise objects
(`N_diag`, `Fs_rn`, `Ts_gwb`, `Phi_gwb`) do not depend on how many CW slots the template
carries, so the same seed must produce the *same bytes* in the `ncw = 16` (C1) and
`ncw = 47` (C2) builds. **Gate G-B6** re-runs 8 of the null seeds inside the C2 build and
requires `max|Δ2F| < 1e−6` and `max|Δ dlnL_field| < 1e−6`. If G-B6 fails, the pooling is
withdrawn and each cell is re-run against its own banked nulls.

### 0.2 B1 — DETECTION, at a matched false-alarm probability

**The field's statistic.** The Earth-term F_e-statistic: a single-source matched filter,
maximised over sky and frequency and over the four amplitude parameters
(`cos_inc, log10 h, phase0, psi`). `2F = 2(lnL_max − lnL_0)`, with `lnL_0` at
`log10 h = −18` (`TE.seed_scan`'s own no-signal convention — adequate for a ONE-slot
reference, which is all this is). The maximisation is the incumbent grid
(`spark.flat_grid` at the problem's **real** span — `TE.seed_scan` hard-codes the T = 15
span, which would be the wrong `df` here) followed by the 8-parameter polish of §1.3–1.4.

**The threshold.** `t95` = the empirical 0.95 quantile of `2F` over the **same 100
pure-noise nulls**, with a bootstrap standard error, and with the zero-fraction column
reported beside it exactly as criterion-v2.2 requires of a floor. The Gumbel fit is
reported alongside; **the empirical quantile is the adopted number** (a max-over-grid
statistic is in the Gumbel domain, so the two should agree — if they disagree by more than
their errors that is itself a reportable fact).

**The criterion's matched number.** The criterion's floor is the 0.95 point of the **D2
offender statistic** over those same 100 nulls — `offender = max over pulsars of dlnL among
pulsars passing (dlnL > lnK) and (q_max > 0.9)`, 0 if none — adopted under criterion-v2.2
(Gumbel if zero-fraction ≤ 0.20, else empirical q95 + bootstrap). The **realisation-level**
detection fraction is then `fraction of realisations with ≥ 1 pulsar at
dlnL > max(lnK, floor)`. That is the object matched to the F_e fraction: both are
"fraction of realisations that clear a FAP = 0.05 bar set on the same nulls."

> **The asymmetry is deliberate and is the measurement.** The criterion carries a
> *per-pulsar trials factor* `ln K` that the field's statistic has no analogue of, and the
> field's statistic carries a *search over sky and frequency* that the criterion's
> conditional pipeline does not pay. Both costs are inside their own null calibration.
> Neither is subtracted.

### 0.3 B2 — the DISTANCE question at baseline

For **F_e-detected** realisations only, the field's standard per-pulsar pulsar-term
inference is run, and scored with the criterion's own vocabulary:

1. **The field's comb.** `dL_a` is the mode spacing at the **ESTIMATED** sky and frequency,
   for a **single** source — the field does not know there are sixteen, so it cannot take
   the min over a population. `EV`, `FringeTables`, `K_counted` all follow from that comb.
2. **The field's template.** Source slot 0 at the F_e estimate; **every** other slot at
   `H_ABSENT` (gate G-B4 proves they are numerically absent — the inherited value of
   −18 FAILED that gate and was replaced by −30; see §1.2).
3. **The field's posterior.** The incumbent `B1Split.estep` + `FringeTables.posterior` —
   the same E-step the criterion uses, on the field's comb and template. Out come
   `q_max`, `map_k`, `dlnL_det`, `lnK` by the field's method.
4. **The scoring.** *A MAP fringe is right or wrong only relative to the comb it was picked
   from*, so `n_true` is **re-derived on the field's own comb**
   (`n_true_field = round((L_true − L0)/dL_field)`). `correct = (map_k == n_true_field)`.
   **purity = correct / confident**, confident = `q_max > 0.9`.

Three purity columns, all at the same bar:

| column | gate applied to the field's picks |
|---|---|
| **FIELD (as practised)** | `q_max > 0.9`. No absolute floor, no trials factor — *this is what the literature reports.* |
| **FIELD + D1 gate** | additionally `dlnL > max(lnK, floor_field)`, where `floor_field` is the D2 offender floor **refit on the field's own columns over the same nulls**, with its own zero-fraction column. Isolates what the criterion's first two layers buy *when applied to the field's method*. |
| **CRITERION v2.2** | the programme's own number, on the same realisations. |

### 0.4 The declared `mc` fork, and its a-fortiori direction

SPARK's s2 note (`spark.py:738`) records that a fixed reference chirp mass is "nearly
harmless" for the **Earth** term and "fatal" for the **pulsar** term, whose look-back is
~kyr. The field's fringe posterior is a pulsar-term object, so this is live.

- **HEADLINE (matched, no oracle).** `log10 mc` is a **free parameter of the polish**.
  Signal and null are treated identically, so the refit floor is matched.
- **BOUND (oracle-generous, SIGNAL ONLY).** The pulsar-term template is re-scored at the
  **true** chirp mass of the injected member the F_e statistic actually locked onto. This
  has no matched null — a pure-noise realisation has no true `mc` — so it is scored **for
  purity only**, which is conditional on confident picks and needs no floor. It bounds how
  much of any field wrongness is `mc`-ignorance rather than the mechanism below.

**Every convention that is free to move is moved to FAVOUR the field**, so a "the field
certifies confidently and wrongly" verdict is a floor on the effect, not a ceiling.

### 0.4a THE ORACLE-SOURCE BOUND — added before the fan, and it is the load-bearing control

A third column, strictly more generous than either above, and declared for the same
a-fortiori reason: **the field's comb and pulsar-term template placed at the TRUE sky, TRUE
frequency and TRUE chirp mass of the member the F_e statistic locked onto**, with only the
amplitude parameters left at the F_e estimate. *This is the field's method with its search
error set to zero.* Purity-only (a pure-noise realisation has no true source), so it owes no
matched floor.

It exists for a specific methodological reason, and the reason is worth stating plainly:
**the headline field number depends on how well my polish converges, and a maximiser that
under-converges would understate the field's sky accuracy and thereby manufacture wrong
picks** — the exact artefact §1.4 warns about, pointing the wrong way. The oracle-source
column removes that dependence entirely. If confident picks are still wrong *there*, the
wrongness cannot be attributed to sky or frequency estimation at all, however good or bad
the maximiser is. Costs one extra E-step (~0.5 s) per signal realisation.

### 0.4b The ELIGIBILITY rule — added during the build, before any realisation was scored

Found by reading the comb construction, not by reading an output. The mode spacing goes as
`dL_a ~ 1/(1 − cos μ_a)`, so it **diverges** for a pulsar lying near the source direction.
On the criterion's comb this never bites (the comb is a min over the whole population, and
`K_counted ≥ 2` for all 116 pulsars). On the **field's** comb — one source, at an
*estimated* sky — it can: a pulsar's ±3σ window can hold **one** fringe. Such a pulsar has
no fringe ambiguity at all, so its `dlnL` (a gap between the best and second-best fringe) is
`+inf` by construction and its `q_max` is trivially `1.0`.

> **RULE (pre-registered, applied to BOTH arms identically): a pulsar whose window holds
> fewer than two sampled fringes is INELIGIBLE** — it is dropped from the confident-pick
> pool and from the offender statistic, in the field arm and the criterion arm alike. The
> count dropped is reported per arm.

The direction matters and is stated: scoring a single-fringe pulsar as an *infinitely
significant, perfectly confident* certification would be an artefact that **flatters the
field**, so removing it is the conservative choice against this report's own expectation.
`n_fringe` is banked as a column for both arms; the incumbent per-arm code paths are left
byte-identical and the exclusion lives only in the scorer.

### 0.5 The pre-registered expectation, and what would refute it

> **EXPECTATION (P14's mechanism).** Above onset the field's confident picks are
> **wrong at a rate its own machinery cannot see**: `q_max` is posterior mass on the MAP
> fringe *given the assumed source*, so a comb built at a mis-estimated sky/frequency
> concentrates mass on a fringe that is not the true one, and concentrates it *confidently*.
> The criterion's absolute floor is the only layer that can see this, because it is the only
> layer calibrated on data containing no source at all.

Pre-registered readouts, stated **before** any number was scored:

| # | readout | what would REFUTE the expectation |
|---|---|---|
| **R1** | FIELD purity on F_e-detected reals at C1 and C2 | purity ≥ 0.9 with a Wilson interval excluding 0.9 from below — the field's confident picks are then simply right, and the criterion's purity layer buys nothing here |
| **R2** | FIELD-vs-CRITERION purity gap | a gap consistent with zero at both cells |
| **R3** | what the D1 gate recovers (FIELD → FIELD + D1) | if D1 restores purity to the criterion's level, then the *floor* — not the *comb* — is the whole story, and the field's method is repairable by a threshold rather than by conditioning |
| **R4** | detection fractions at matched FAP | if F_e detection fraction ≥ the criterion's at both cells, the criterion costs sensitivity and buys purity; **if it is lower, the criterion costs nothing on detection** — either is a result, and neither was assumed |
| **R5** | the `mc` BOUND | if oracle-`mc` purity ≈ 1, the field's wrongness is `mc`-ignorance, **not** the comb mechanism, and R1's headline must be re-attributed |

**A null result is publishable here.** If the field's baseline matches the criterion on both
axes, the paper's contribution is a *calibration* claim, not a *correction* claim, and this
report says so.

### 0.5b The BUDGET RULE — stated before the per-realisation wall was measured

The fan as specified is **258 realisations**: C1 30 signal + 100 null, C2 120 signal +
8 cross-build-gate nulls. Each carries one F_e grid scan (14 976 points) plus **three**
E-steps for a signal realisation (criterion, field, field-at-oracle-`mc`) or **two** for a
null.

> **CAP: 12 GPU-hr for the whole campaign**, smokes and venue builds included. If the
> measured per-realisation wall puts the full fan over that, the trim order is fixed **in
> advance** and is applied in this order until it fits:
> 1. the oracle-`mc` BOUND drops to **C1 only** (it is a bound, not the headline);
> 2. **C2 signal** is subsampled to the largest multiple of 8 that fits, **keeping all
>    8 skies equally represented** — sky coverage is the entire reason C2 is the cell that
>    survived, so the sky axis is never the thing that gets thinned;
> 3. nothing else. In particular **the 100-realisation null set is never cut**: it is what
>    sizes both thresholds, and criterion-v2.2's own D2 rule requires N ≥ 100.
>
> Whatever is dropped is **stated in the results section with the number dropped**, never
> silently.

### 0.6 The gate chain (all must pass before the fan)

| gate | what it proves |
|---|---|
| **G-B0** | `L_gwb` is **BANKED** (`GENERALISE_results/gen_L_gwb_n5336.npz`), so every noise draw is thread- and host-independent. BASELINE **refuses to recompute** it. |
| **G-B1** | the rebuilt comb reproduces the source campaign's **noise-free** columns exactly (`max|Δ lnK| = 0`) — i.e. the geometry path here *is* SURFACE's / GENERALISE's |
| **G-B2** | the paired reproduction residual against the source bank. **C2 must be bit-identical** (it was drawn under the banked `L_gwb`). **C1 cannot be**: SURFACE's T = 30 realisations predate the bank and were drawn under the `eigh` recompute path, so they carry a different — equally valid — GWB rotation. C1 is therefore **re-drawn, and BOTH arms are recomputed on the re-draw**, making the comparison paired within this campaign and only *distributional* against SURFACE's bank. The residual is measured and reported, never hidden. |
| **G-B3** | the 8-parameter polish never *lowers* `2F` |
| **G-B3b** | *(added after G-B3 failed — §1.4)* the polished sky is ≥ 2× closer to the injected source than the grid pixel. A gate that only checks "did not get worse" is not enough when the quantity is an **input to the next stage**. |
| **G-B4** | `H_ABSENT` slots are numerically absent (−6 dex moves `lnL` by < 1e−6) |
| **G-B5** | `field_comb` reproduces the incumbent single-source mode spacing at true parameters |
| **G-B6** | the pure-noise null set is venue-only (`ncw = 16` build == `ncw = 47` build) — licenses one null set for both cells |

---

## 1. THE GATE CHAIN — and the six defects it caught before a single realisation was banked

**Smoke #1: job `12833821`, dgx01, A100-40GB, `cpus-per-task=8`, arm C1.**

### 1.1 What passed immediately

| gate | verdict | number |
|---|---|---|
| **G-B0** | **PASS** | `BANKED gen_L_gwb_n5336.npz fp=f92c9e36b460d6f5 cpus=8`. **This is the same fingerprint the GENERALISE C2 realisations were drawn under on `hgx03`** — so the noise draws here, on a different host and a different GPU, are the same bytes. |
| **G-B1** | **PASS** | `max|Δ lnK| = 0.000e+00` against `sf_sig_h1275_T30_lit_k3_g-1_n20020000.npz`. The geometry path used here **is** SURFACE's, exactly. |
| **G-B5** | **PASS** | `max|Δ dL| = 0.000e+00`. `field_comb` at true parameters reproduces the incumbent single-source mode spacing. |

The venue is also confirmed as the right one: `EXTENSION dT=15yr: +42835 TOAs, span
22.15 → 37.14 yr, rn 50, gwb 23`, giving `npsr × ngp = 116 × 46 = 5336` — the shape of the
banked `L_gwb`.

### 1.2 DEFECT 1 — the inherited "absent source" strain is not absent (G-B4 FAIL)

    *** FAIL  G-B4 H_ABSENT slots numerically absent (-6 dex moves lnL < 1e-6)  |dlnL| = 1.003e-02

`H_ABSENT = −18` is `TE.seed_scan`'s own no-signal convention and is what SPARK adopted.
At this venue, **fifteen** slots parked there still carry **1.003e−02 nat** of spurious
source. The convention is fine for the ONE-slot no-signal reference TE actually used it for;
it does not survive being used for fifteen. **Fixed: `H_ABSENT = −30`** (GLACIER-LITE's
value, inside chorus's `_PHI_LO` box of −32).

*This is reported, not patched upstream:* `TE.seed_scan` and `spark.py` are untouched. The
value is BASELINE-local.

### 1.3 DEFECT 2 — my polish DESTROYED the statistic (G-B3 FAIL)

    [Fe] grid 14976 pts 1369.4s  2F_grid=    25.4 -> 2F_polish= 8.0     (null)
    [Fe] grid 14976 pts  103.6s  2F_grid= 13131.4 -> 2F_polish= 2.5     (signal)
    [Fe] grid 14976 pts  103.2s  2F_grid= 13362.3 -> 2F_polish=-0.6     (signal)
    *** FAIL  G-B3 polish never lowers 2F   d2F = -17.38, -13128.92, -13362.83

Not a fluke — a design error, and the gate caught it on the first three draws. Adam is a
fixed-step method, not a line search: `mh/√vh ≈ ±1`, so the step is `≈ lr·scale`
**regardless of the gradient**. At the incumbent's `lr = 0.05` that is `0.05 × 3.14 = 0.157`
rad `= 9°` in `gwphi` — **coarser than the healpix cell the polish exists to refine**. It
walks off the peak the grid already found and 200 steps need not walk back. Returning the
last iterate then throws away a maximum that had already been evaluated.

**Fixed, two changes** (`make_polish`):
1. **Return the best point the trajectory ever visited, seeded with the grid argmax.** A
   detection statistic is a maximum over evaluated points; returning less than a point
   already evaluated is simply a bug in the maximiser. Free — `value_and_grad` already
   computes the value the gradient is taken at.
2. **Anneal the step** by `1e−4` over the 200 steps, so the refinement converges instead of
   orbiting at radius `lr·scale`. Final step ≈ `1e−3` degrees.

### 1.4 DEFECT 3 — the fix for defect 2 exposed a *direction* problem (new gate G-B3b)

Change (1) alone would have satisfied G-B3 while leaving the polish **inert**: with a fixed
step it never improves, so `p_hat` stays at the grid argmax and the field's sky estimate is
stuck at the healpix half-width, ~7°. B2 builds the field's fringe comb **at the estimated
sky** — so an artificially coarse sky error would have **manufactured** exactly the
confidently-wrong picks this campaign is testing for. That is the wrong direction: every
free convention here must favour the field. Hence change (2), and a new gate:

> **G-B3b: the polished sky must be ≥ 2× closer to the injected source than the grid pixel**,
> on every signal draw.

*A gate that only checks "did not get worse" is not enough when the quantity is an input to
the next stage.* — **and the gate immediately proved its own point.** Smoke #2 (job
`12834517`), with best-point tracking in place and the anneal started at the incumbent's
`lr = 0.05`:

    [Fe] grid 14976 pts  2F_grid= 13131.4 -> 2F_polish= 13131.4
    [Fe] grid 14976 pts  2F_grid= 13362.3 -> 2F_polish= 13362.3
    *** FAIL  G-B3b  grid 6.37 deg -> polish 6.3724 deg   (improved on the grid in 0/3)

**Safe, and completely inert.** `6.37°` is the healpix half-width, exactly as predicted:
the polish still leaves the peak on step 1 and the annealed remainder wanders without
returning. Under the old code this would have passed a "never lowers" gate while silently
handing B2 a 6.4°-wrong sky.

**And the C2 smoke explains the whole thing.** Same code, same annealed `lr₀ = 0.05`,
different geometry:

| arm | grid argmax distance from the source | polish |
|---|---|---|
| **C1** | **6.37°** | `6.37° → 6.3724°` — **inert** (`d2F = −0.00, −0.00, −0.00`) |
| **C2** | **11.38°** | `11.38° → 2.51°, 2.39°` — **PASS**, 4.5× (`d2F = +6486, +6341`) |

*A 9° first step helps when the argmax is 11° away and jumps clean past when it is 6° away.*
The behaviour was never about the venue or the source; it was the ratio of the step size to
the initial distance, and a campaign that had only ever smoked C2 would have shipped a
maximiser that silently fails on half its cells.

**The schedule is now DERIVED from the two lengths that matter, not tuned:**

| | |
|---|---|
| how far it must **travel** | up to a healpix half-width, ~7° |
| how finely it must **land** | the peak of a `2F ≈ 1e4` source is far narrower than a grid cell |

`lr₀ = 0.004` ⇒ first step **0.72°** (resolves the peak, does not jump the cell);
400 annealed steps give **31.6°** of cumulative travel (≫ the 7° ever needed) and a final
step of **7.4e−05°**. Because the best point is tracked, *a mis-set schedule can only ever
cost inertness, never damage* — which is why this failure showed up as an **equality**
rather than a loss, and why it would have been easy to miss without G-B3b.

### 1.4b DEFECT 4 — G-B3's own bar was absolute where it had to be relative

    *** FAIL  G-B3  d2F = -0.00, -0.00, -0.00

The polish returned the grid value — and was flagged anyway. The grid statistic comes out of
the **vmapped** profile and the polish's seed value out of a **non-vmapped** call to the same
likelihood; XLA fuses the two differently, so they agree only to ~`1e−6` **relative**. At
`2F ≈ 1.3e4` that is ~`1e−2` absolute, and a `−1e−9` absolute bar fails on *arithmetic*, not
on anything about the maximiser. Re-cut to EMBER 2.2(b)'s continuous-column convention
(`≥ −1e−6 × 2F`), in both the gate and the scorer's `WORST change` diagnostic.

*Worth stating because the failure mode is generic: a bar copied from a discrete-column gate
onto a continuous column of a different magnitude will fire on floating point and look like
physics.*

### 1.5 DEFECT 5 — my own reasoning about the JAX compile-cache race was wrong

I switched the fan from the harbor-default per-task compile cache to a shared one, arguing
the documented race (array tasks colliding on
`xla_gpu_per_fusion_autotune_cache_dir/tmp/*.textproto`) cannot fire because `harbor_env`
pins `--xla_gpu_autotune_level=0`. **The smoke created and kept writing that very directory
under exactly that flag.** Reverted to per-task caches — and then made isolation free by
**seeding each task's cache from the arm's warm smoke cache**, so a private directory does
not mean a private recompile.

A fourth scan, on a real C1 signal realisation, makes the diagnosis airtight:

    [Fe] grid 14976 pts 103.6s  2F_grid= 15414.4 -> 2F_polish= 17962.2     (+2548, an IMPROVEMENT)

**The un-annealed polish sometimes helped and sometimes destroyed the statistic** — on the
same code, the same venue, four consecutive draws. That is not a maximiser with a tuning
problem; it is a random walk that happens to be initialised at a good point. The C2 smoke
adds the fourth outcome on a *different* geometry with the *annealed* version — an
improvement of **+6486** (`2F_grid = 21268.6 → 27754.7`) where the identical code on C1 was
inert. **Destroys / improves / inert / improves-a-lot, across two venues.** A statistic
whose value depends on which of those the walk happens to land in is not a statistic. Under the fixed
version the improvement can only be ≥ 0 by construction, and the interesting number becomes
*how often and by how much* it improves — which is banked (`twoF_grid` beside `twoF`) and
reported in §2.

### 1.5b G-B2 — the C1 residual, and why it settles a methodological question

    PASS  G-B2  max|dq| = 5.056e-01   max|d dlnL| = 4.091 nat     [ROTATED, not bit-identical]

As pre-registered: SURFACE's T = 30 realisations predate the shape-keyed `L_gwb` bank and
were drawn under the `eigh` recompute path, so at the same seed this run draws a **different
— equally valid — GWB rotation**. The size of the residual is the point: **`q_max` moves by
up to 0.51 and `dlnL` by up to 4.1 nat on a single pulsar.** Those are not rounding errors;
they are a different realisation of the same statistical object.

> **This retires an option that looked cheaper.** Reading SURFACE's *banked* criterion
> columns and putting the field's freshly-computed columns beside them would have compared
> **two different realisations** and called the difference a method effect. Recomputing BOTH
> arms on the same re-drawn data — which is what this campaign does — is not belt-and-braces;
> at this residual size it is the only way the comparison means anything.

C2 is the control on that statement: it *was* drawn under the banked `L_gwb`, so its G-B2
must reproduce. **MEASURED (job `12834510`):**

    G-B2  C2 vs the GENERALISE bank:  max|dq| = 1.171e-09   max|d dlnL| = 6.519e-09

> **The banked-`L_gwb` host-independence claim holds across hosts.** GENERALISE drew those
> realisations on **`hgx03` (H200)**; this run redrew them on **`dgx01` (A100-40GB)**, a
> different host and a different GPU, and the criterion's `q_max` and `dlnL` come back
> agreeing to **~1e−9**. Set against C1's `0.51` / `4.09` for the *unbanked* path, the two
> numbers together are a clean demonstration of exactly what banking the GWB square root
> buys: it converts a 0.5-in-`q_max` host/thread systematic into floating-point noise.

**DEFECT 6 — and it is the same mistake as DEFECT 4.** I wrote that gate as **exact
bit-identity** (`dq == 0.0 and dd == 0.0`), so it reported `FAIL` on a `1e−9` agreement.
`q_max` and `dlnL` are *continuous* columns and this driver reaches them by a different
summation order over the `ncw` slots than `generalise.py` does, so `~1e−9` is the correct
expectation and exact equality was never the right bar. Re-cut to EMBER 2.2(b)
(*discrete exact, continuous `< 1e−6`*).

*Twice now I have put a discrete-column bar on a continuous column and had it fire on
arithmetic. Both times the underlying measurement was fine and the gate was wrong — which is
the failure mode that most easily gets "fixed" by loosening a threshold until it passes, so
it is recorded here with the reasoning rather than the number.*

### 1.6 The measured cost model (this is what sizes the fan)

| stage | wall | note |
|---|---|---|
| venue build (`build_chorus_problem`, `n_ecc=0`, `ncw=16`, T=30) | **480 s** | once per task |
| F_e detector + grid (`build_fisher_amortised` ncw=1 + `flat_grid`) | **105 s** | once per task; grid = **14 976 pts** (78 freq × 192 healpix), span 37.14 yr |
| F_e scan, **cold** (first call, compiles the vmapped 40-step Adam graph) | **1 369 s** | once per process — **the reason the fan seeds a warm cache** |
| F_e scan, **warm** | **103.6 s / 103.2 s** | the per-realisation cost that matters |
| first E-step in a process | **~6–8 min** — compiles 116 per-pulsar graphs (the TURBINE "SMASK recompile tax") plus the one-time `_ensure_minv` 5336×5336 solve | once per task |
| every later E-step | **~0.5 s** | |

> **The warm cache does NOT remove the first-E-step cost, and it cannot.** The house sets
> `jax_persistent_cache_min_compile_time_secs = 10`, and the 116 per-pulsar E-step graphs
> each compile in **less** than that, so they are never persisted — the seeded cache holds
> 27 large entries (the F_e scan graph and friends), not 143. Measured directly: a fan task's
> private cache grows 27 → ~40 entries while it runs, writing `jit_run-*` files. So the
> hardlink seed buys the **expensive** compile (the vmapped Adam F_e graph, ~21 min) and each
> task still pays ~8 min of cheap-but-numerous E-step compiles. 4 tasks × 8 min ≈ 32 min of
> unavoidable overhead — priced in, not a stall.

**An unplanned confirmation of the pooled null set.** The C1 (`ncw = 16`) and C2
(`ncw = 47`) smokes both run G-B3 on the *same* gate seed, and both return
`2F_grid = 25.4` — **identical across the two builds.** That is G-B6's premise (a pure-noise
draw contains no source, so the Fe statistic cannot depend on how many CW slots the template
carries) showing up for free, at the gate level, before the fan's own 8-seed version of the
test. It is not a substitute for G-B6, but it means the pooling was already visibly holding.

**A sanity check worth recording:** the pure-noise `2F_grid` values are **25.4, 22.5**. For
a χ²₄ statistic maximised over ~15 000 correlated grid points the expected maximum is
`≈ 2 ln N` plus correction terms, i.e. mid-20s. The F_e implementation is behaving like an
F-statistic before anything is scored against it. The signal draws sit at `2F ≈ 1.3–1.8e4`.

**The per-realisation wall, measured end to end:**

    [BL smoke] sig n20020000: 2F=17962.2 fq>0.9:5 cq>0.9:61 105s (total 105s)

**105 s, complete** — F_e scan *and* all three E-steps (criterion, field, field-at-oracle-mc).
Once the 116 per-pulsar graphs are compiled an E-step costs ~1 s, so **the F_e scan is
essentially the entire cost** and the three-E-step design is nearly free.

> ### BUDGET VERDICT: the full pre-registered fan fits, and NOTHING IS TRIMMED.
> 258 realisations × 105 s = **7.5 GPU-hr**, plus ~0.7 GPU-hr of per-task venue builds and
> ~2.1 GPU-hr of smokes ⇒ **≈ 10.3 GPU-hr against the 12 GPU-hr cap.** Neither trim in §0.5b
> is invoked: the oracle-`mc` BOUND runs on **both** arms and C2 keeps all **120** signal
> realisations across all 8 skies.

### 1.7b The gate ledger across the smoke cycles

| gate | smoke #1 (`12833821`, C1) | smoke #2 (`12834517` C1 / `12834510` C2) | smoke #3 (`12834778`, C1) |
|---|---|---|---|
| **G-B0** `L_gwb` banked | PASS | PASS (both arms) | |
| **G-B1** comb vs source campaign | PASS `0.000e+00` | PASS `0.000e+00` (both arms) | |
| **G-B2** residual vs the bank | `[ROTATED]` measured | `dq=5.056e−01`, `d dlnL=4.091` — **bit-reproducible against smoke #1** | |
| **G-B3** polish never lowers `2F` | **FAIL** `−13129` | **FAIL** on the *bar*, not the maths (§1.4b) | |
| **G-B3b** polish refines the sky | *(gate did not yet exist)* | **FAIL** `6.37° → 6.3724°` (inert) | |
| **G-B4** `H_ABSENT` absent | **FAIL** `1.003e−02` | **PASS `7.451e−09`** (both arms) | |
| **G-B5** `field_comb` | PASS `0.000e+00` | PASS `0.000e+00` (both arms) | |
| | 4/6 | **5/7** (C1) / **5/7** (C2) | **C1 7/7 — ALL PASS**; **C2 6/7** (only G-B3b, §1.7d) |

**FINAL GATE STATE, both arms** (C1 job `12834778`, C2 job `12834976`):

| gate | C1 | C2 |
|---|---|---|
| G-B0 `L_gwb` banked | PASS | PASS |
| G-B1 comb vs source campaign | PASS `0.000e+00` | PASS `0.000e+00` |
| G-B2 residual vs the bank | PASS `[ROTATED]` `0.506` / `4.09` — measured, and it settled the design (§1.5b) | **PASS `1.171e−09` / `6.519e−09` cross-host** |
| G-B3 polish never lowers `2F` | PASS, improved **3/3** (`+1.7, +1.9e3, +1.9e3`) | PASS, improved **3/3** (`+1.7, +3.97e3, +4.03e3`) |
| G-B3b polish refines the sky | PASS `6.37° → 0.37–0.55°` | **FAIL `11.38° → 8.70°`** — the one OPEN item, §1.7d |
| G-B4 `H_ABSENT` absent | PASS `7.451e−09` | PASS `7.451e−09` |
| G-B5 `field_comb` | PASS `0.000e+00` | PASS `0.000e+00` |
| | **7/7** | **6/7** |

Both bars I had mis-specified (G-B2's bit-identity, G-B3's absolute tolerance) now pass on
the corrected, *stated* convention. The single remaining failure is G-B3b on C2, which is
carried openly into the results rather than tuned away.

**Smoke #3 (`12834778`, C1) closed the chain at `GATES: 7/7 pass`, wall 1192 s, and the
C1 fan launched on that condition** (job `12834949`, array `0-3%4`). C2's two remaining
failures were both MY BARS, not its physics (§1.4b, §1.5b) and are re-gated in its own
smoke #2.

### 1.7c THE POLISH, RESOLVED — smoke #3 (`12834778`, C1), derived schedule

    [Fe] grid 14976 pts  2F_grid=    25.4 -> 2F_polish=    27.1
    [Fe] grid 14976 pts  2F_grid= 13131.4 -> 2F_polish= 14989.3
    [Fe] grid 14976 pts  2F_grid= 13362.3 -> 2F_polish= 15208.5
    PASS  G-B3   d2F = +1.7, +1.86e+03, +1.85e+03   (improved on the grid in 3/3)
    PASS  G-B3b  grid 6.37 deg -> polish 0.5522 deg   grid 6.37 deg -> polish 0.3655 deg

**On the identical C1 geometry where the previous schedule was inert (`6.37° → 6.3724°`),
the derived one lands at `0.37–0.55°`** — an **11–17× refinement against a 2× bar**, and it
improves on **3/3** draws including the pure-noise one. The predicted first step (0.72°) and
cumulative reach (31.6°) were right.

**And it recovers real signal:** `2F` goes `13131 → 14989` (**+14 %**). The 14.7° grid was
leaving that on the table, which is a fact about the *incumbent seeder's* resolution, not
just about my polish — worth knowing wherever `TE.seed_scan`'s grid is used as anything
other than a seeder.

Three smoke cycles, and **not one of them was wasted**: each failure was a different real
defect, and every one of them was invisible in the science output — the campaign would have
produced a complete, plausible, wrong set of numbers without them. Smoke #2's G-B2 also
reproduced smoke #1's residual *exactly*, which is its own small check: the re-draw is
deterministic.

### 1.7d OPEN ITEM — G-B3b FAILS ON C2, and I am running the fan anyway. Here is why.

C2's re-gate (job `12834976`) with the derived schedule:

    PASS  G-B3   d2F = +1.7, +3.97e+03, +4.03e+03   (improved on the grid in 3/3)
    *** FAIL  G-B3b  grid 11.38 deg -> polish 8.6966 deg / 8.6811 deg   (1.31x, bar is 2x)

**Neither schedule covers both arms.** Put the four measurements side by side:

| schedule | C1 (argmax **6.37°** out) | C2 (argmax **11.38°** out) |
|---|---|---|
| annealed from `lr₀ = 0.05` | `→ 6.3724°`, `d2F = −0.00` — inert | `→ 2.51°`, `d2F = +6486` — **best** |
| annealed from `lr₀ = 0.004`, 400 steps | `→ 0.37–0.55°`, `d2F = +1.9e3` — **best** | `→ 8.70°`, `d2F = +3970` — under-travels |

The fast schedule overshoots a near argmax; the slow one cannot reach a far one before the
anneal kills the step. The correct fix is a **two-stage polish** (coarse travel, then fine
refine, best-tracked across both). **I am not applying it**, and the reason is not that it
is hard:

> **THE C1 FAN IS ALREADY RUNNING WITH THIS ESTIMATOR.** Changing the polish now would mean
> C1 and C2 were scored with *different maximisers*, and the entire campaign is a
> **comparison between two arms on matched realisations**. A shared, imperfect estimator is
> worth far more here than two different better ones. Re-running C1 to match would cost
> ~1.5 GPU-hr and put the campaign over its pre-registered 12 GPU-hr cap for a quantity that
> §0.4a already brackets.

**What protects the science, and what does not:**

- **G-B3 passes 3/3 on both arms with large improvements** — the polish *is* a valid
  maximiser (it never lowers `2F` and it finds substantially better points). B1's detection
  statistic is therefore sound on both arms. **Nothing in B1 depends on G-B3b.**
- **B2 on C2 is affected**: the field's comb is built at a sky ~8.7° from the nearest
  injected member. If that is under-convergence, it biases the field's purity **downward** —
  i.e. **against** the field, in the direction that would flatter this report's own
  expectation. That is the wrong direction, and it is why the next point matters.
- **The ORACLE-SOURCE column (§0.4a) is exactly the control for this.** It places the comb at
  the member's TRUE sky/frequency/`mc`. **The gap between C2's headline purity and its
  oracle-source purity bounds the whole effect**, whatever the polish did. This is the reason
  that column was added before any of this was known, and it is now load-bearing rather than
  merely prudent.

**A competing hypothesis I can test from banked data, and will not assert until I do.** C2
carries **five** loud members plus an eccentric comb. The maximum of a **single-source** F_e
statistic against a five-source sky **need not sit on any one member** — so "distance to the
nearest injected member" may be measuring a genuine blend position rather than a convergence
failure. `p_grid` and `p_hat` are banked on every realisation, so after the fan this is
decidable from disk: if the polished points cluster between members rather than short of the
nearest one, G-B3b is the wrong metric for a multi-source cell, not a failing maximiser.
**Stated as a hypothesis, flagged as OPEN, not used to excuse the failure.**

### 1.7e The eligibility rule FIRES on real data (predicted in §0.4b, confirmed here)

First banked realisation off the fan (`bl_sig_..._g-1_n20020000.npz`), read back:

    f_nfringe  min = 1     median = 512      <- the FIELD's comb
    c_nfringe  min = 25    median = 512      <- the CRITERION's comb
    field comb dL_med = 0.0006 kpc   vs   criterion dL_med = 0.0002 kpc  (3x coarser)

**`min = 1` on the field's comb and never on the criterion's** — exactly the asymmetry §0.4b
predicted from `dL_a ~ 1/(1 − cos μ_a)` diverging near the *estimated* source direction. The
eligibility rule is removing real cases, not guarding a hypothetical one, and each removed
case would otherwise have scored as an **infinitely significant, perfectly confident**
certification — i.e. as a point in the field's favour. Every scorer column was also verified
present on real banked data before the full run (`MISSING: none`).

### 1.7 The fan plan, verified against disk before submission

All **258** realisations were checked to map to a file that actually exists in the source
campaigns' banks — the seed algebra is not assumed, it is confirmed:

    C1: 130 entries (30 sig + 100 null)   source-bank misses: 0
    C2: 128 entries (120 sig + 8 xnull)   source-bank misses: 0

C2's per-sky eccentric placements were also reproduced from `generalise.as_placement` and
checked against the banked npz for all 8 skies: `[2, 0, 3, 3, 2, 3, 3, 0]`, **8/8 exact**.
`nshard = 4` splits both arms evenly (33/33/32/32 and 32×4). *None of this cost a
GPU-second.*

---

## 2. B1 — DETECTION AT MATCHED FAP = 0.05 (C1 COMPLETE, 130/130)

C1 fan job `12834949`, 4 shards, all `exit=0`, walls 4484–4869 s. **30 signal + 100 null,
both arms recomputed on the same data.**

| | N_sig | N_null | threshold | est | zero-f | **detected** |
|---|---|---|---|---|---|---|
| **FIELD F_e** | 30 | 100 | `2F > 28.28 ± 0.80` | emp_q95 | 0.00 | **30/30 = 1.00** [0.89, 1.00] |
| **CRITERION v2.2** | 30 | 100 | `dlnL > max(lnK, 38.59 ± 1.72)` | gumbel | 0.00 | **20/30 = 0.67** [0.49, 0.81] |

    null 2F: median 23.07, max 31.04      signal 2F: median 15591, max 23976
    Gumbel fit of the F_e null: 29.53 ± 0.75 -- agrees with the adopted emp_q95 28.28 ± 0.80

> ### R4 ANSWERED, AND IT IS THE OPPOSITE OF A FREE LUNCH.
> **At a matched false-alarm probability the field's F_e-statistic detects MORE than the
> criterion does — 1.00 against 0.67.** The criterion is not simply a stricter reading of the
> same evidence; at the realisation level it **costs real detection sensitivity**. Whatever
> it buys has to be bought somewhere else. §3 is where.

The F_e null behaves exactly as a searched χ²₄ should (median 23, max 31 over ~15 000
correlated grid points), and the signal sits four orders of magnitude above it, so the
detection side of the field's method is not in any doubt here.

---

## 3. B2 — THE DISTANCE QUESTION AT BASELINE (C1)

Confident picks are `q_max > 0.9`, scored on the **F_e-detected** realisations (all 30),
each method on **its own comb**:

| method | conf/real | correct/real | **wrong/real** | **purity** | Wilson 95 % | n_conf |
|---|---|---|---|---|---|---|
| **FIELD (as practised)** | **9.53** | 0.73 | **8.80** | **0.08** | [0.05, 0.11] | 286 |
| FIELD, truth in window | 9.53 | 0.73 | 8.80 | 0.08 | [0.05, 0.11] | 286 |
| FIELD + D1 gate | 2.83 | 0.43 | 2.40 | 0.15 | [0.09, 0.24] | 85 |
| FIELD, oracle `mc` | 3.37 | 0.83 | 2.53 | 0.25 | [0.17, 0.34] | 101 |
| **FIELD, ORACLE SOURCE** | 2.67 | 0.83 | 1.83 | **0.31** | [0.22, 0.42] | 80 |
| **CRITERION v2.2** | **1.07** | 0.83 | **0.23** | **0.78** | [0.61, 0.89] | 32 |

> ### THE HEADLINE.
> **The field's standard method states ~9.5 confident pulsar-distance identifications per
> realisation, and 92 % of them are wrong** (purity 0.08, Wilson [0.05, 0.11]). The
> criterion states **one** per realisation and 78 % of those are right. **Both arms recover
> the same number of CORRECT identifications — 0.73 vs 0.83 per realisation.** The criterion
> is not finding more truth; it is declining to state ~8.6 falsehoods per realisation that
> the field states with `q > 0.9` confidence.

**Scored against the pre-registration (§0.5), C1:**

| # | readout | verdict |
|---|---|---|
| **R1** | FIELD purity on F_e-detected reals | **0.08**, Wilson upper bound **0.11** — the ≥ 0.9 refutation condition is excluded by a mile. **Expectation CONFIRMED.** |
| **R2** | FIELD-vs-CRITERION purity gap | **0.70**, intervals disjoint ([0.05, 0.11] vs [0.61, 0.89]). **Not consistent with zero.** |
| **R3** | what the D1 gate recovers | 0.08 → **0.15** only. **The floor is NOT the story.** Adding the criterion's first two layers to the field's own comb recovers almost nothing — so the defect is the **comb and the single-source model**, not a missing threshold. This is the readout I would most have expected to go the other way. |
| **R4** | detection at matched FAP | F_e **1.00** vs criterion **0.67** — the field detects MORE (above). |
| **R5** | the `mc` and ORACLE-SOURCE bounds | oracle-`mc` **0.25**, **oracle-source 0.31**. `mc`-ignorance and search error together explain only part of it: **with the field's search error set to exactly zero, 69 % of its confident picks are still wrong.** |

**R5 is the load-bearing one, and it is why §0.4a existed.** The oracle-source column places
the comb at the member's TRUE sky, frequency and chirp mass. It cannot be blamed on my
polish, on the healpix grid, or on `mc`. **A purity of 0.31 there means the field's
confident-but-wrong picks are a property of building a fringe comb from ONE source against a
sky that holds sixteen** — precisely P14's mechanism, measured rather than asserted.

**Two supporting facts, both as predicted:**
- The field's own per-pulsar floor comes out at **0.02 ± 0.02 nat with zero-fraction 0.43** —
  *above* criterion-v2.2's 20 % gate, so the empirical-quantile estimator was correctly
  substituted for the Gumbel. The convention did real work on a real column.
- Single-fringe exclusions (§0.4b): **FIELD 10/3480 signal and 15/11600 null; CRITERION 0 and
  0.** Only 1/3480 pulsar-realisations had its true fringe outside the field's 512-point
  window, so `field_inwin` is identical to the headline and the window budget is not doing
  any of the work.

### 3.1 The external check — this re-draw reproduces SURFACE's published cell

| quantity | BASELINE re-draw | SURFACE bank |
|---|---|---|
| criterion floor (nat) | 38.59 ± 1.72 | 37.68 ± 1.63 |
| realisations detected | 0.67 | 0.67 |
| certifications / realisation | 1.07 | 1.10 |
| **correct** certs / realisation | 0.83 | 0.87 |
| **wrong** certs / realisation | 0.23 | 0.23 |

Different GWB rotation (G-B2), so these are different realisations of the same cell — and
they agree to well inside the 30-realisation scatter on every column. **The criterion arm
computed here IS the programme's own census cell**, which is what licenses putting the field
column beside it.

### 3.2 The polish, in production

    improved on 100 % of signal and 100 % of null realisations
    median gain +1859 (signal) / +1.41 (null) in 2F, max +4804
    WORST change +1.64e-01  (bar: >= -1e-6 x 2F)
    F_e-detected sky error vs the nearest injected member: median 3.3 deg (0.6-10.0)

G-B3 holds on all 130 production realisations, not just the three gate draws.

---

## 4. THE PANEL

`reports/BASELINE_panel.png` — **A**: detection fraction at matched FAP = 0.05, F_e vs
criterion, Wilson intervals. **B**: purity of confident picks, FIELD / FIELD+D1 / CRITERION,
with `n_conf` beside each point. Tables in `reports/BASELINE_tables.txt`, raw scores in
`reports/baseline_score.npz`.

---

## 2b–3b. C2 — THE A-SKY SURVIVOR (COMPLETE, 128/128; campaign 258/258)

C2 fan job `12835126`, 4 shards, all `exit=0`, walls 4589–4706 s.

### G-B6 PASSES — the pooled null set is licensed by measurement

    8 seeds re-run in the ncw=47 build:  max|d 2F| = 0.000e+00   max|d dlnL_field| = 3.725e-09

`2F` is **exactly** bit-identical across the `ncw = 16` and `ncw = 47` builds and the field's
per-pulsar `dlnL` agrees to 4e−09. The §0.1 premise — a pure-noise draw contains no source,
so the F_e statistic cannot depend on how many CW slots the template carries — is confirmed,
not assumed. One null set, both cells.

### B1 (C2)

| | threshold | detected |
|---|---|---|
| **FIELD F_e** | `2F > 28.28 ± 0.80` | **120/120 = 1.00** [0.97, 1.00] |
| **CRITERION v2.2** | `dlnL > max(lnK, 38.59 ± 1.72)` | **120/120 = 1.00** [0.97, 1.00] |

`null 2F` median 23.07 / max 31.04; `signal 2F` median **18841**, max 35421.

**R4 is cell-dependent, and that is the honest statement.** At C1 the criterion cost
detection (0.67 vs 1.00); at C2 — louder, 5 loud members plus an eccentric comb — **both
saturate at 1.00 and the criterion costs nothing.** The C1 gap is a property of a cell near
its own onset, not a general tax.

### B2 (C2)

| method | conf/real | correct/real | wrong/real | purity | Wilson 95 % | n_conf |
|---|---|---|---|---|---|---|
| **FIELD (as practised)** | **9.30** | 0.66 | **8.64** | **0.07** | [0.06, 0.09] | 1116 |
| FIELD + D1 gate | 3.08 | 0.47 | 2.60 | 0.15 | [0.12, 0.19] | 369 |
| FIELD, oracle `mc` | 3.71 | 0.68 | 3.03 | 0.18 | [0.15, 0.22] | 445 |
| **FIELD, ORACLE SOURCE** | 4.13 | 0.60 | 3.53 | **0.15** | [0.12, 0.18] | 496 |
| **CRITERION v2.2** | **20.02** | **19.24** | **0.78** | **0.96** | [0.95, 0.97] | 2403 |

> ### THE RESULT REPLICATES, AND HARDENS.
> The field's purity is **0.07** at C2 against **0.08** at C1 — *identical within the
> intervals, on a different population, a different structure, and 4× the realisations.*
> The criterion at C2 certifies **20 per realisation at 96 % purity**. **Field 8.64 wrong
> per realisation; criterion 0.78.**
>
> **And the ORACLE-SOURCE bound is LOWER at C2 (0.15) than at C1 (0.31).** With five loud
> members instead of three, zeroing the field's search error helps *less*, not more —
> exactly what the mechanism predicts if the defect is **a one-source comb against a
> many-source sky**. The wrongness scales with the population, not with the search.

### 3b.1 G-B3b — the OPEN item of §1.7d, RESOLVED FROM BANKED DATA

As pre-registered, decided from `p_grid`/`p_hat` on all 150 signal realisations rather than
argued:

| | polish moved | **final sky error** | < 7° (healpix half-width) | < 2° |
|---|---|---|---|---|
| **C1** (3 loud) | 5.18° | **3.31°** (0.60–9.96) | 87 % | **43 %** |
| **C2** (5 loud + ecc comb) | 4.83° | **3.98°** (0.46–16.52) | 82 % | **19 %** |

**The gate's verdict does not survive contact with the production ensemble.** G-B3b failed on
two draws that happened to sit at `11.38°`; across 120 C2 realisations the polish moves a
median `4.83°` and lands at a median `3.98°`, with 82 % inside the healpix half-width. *The
maximiser refines on C2.*

**But the hypothesis was only half right.** C2 lands tightly (< 2°) on **19 %** of
realisations against C1's **43 %** — the multi-source blend effect is real, just smaller than
the gate implied: with five loud members the single-source F_e maximum sits *between* them
and cannot converge onto any one. **So G-B3b's failure was an unrepresentative draw AND a
metric that means less for a multi-source cell — both, not either.** Recorded as measured;
the 2× bar is retained unchanged for any future single-dominant-source cell.

*The decision in §1.7d — keep one estimator across both arms rather than fix the polish
mid-campaign — is vindicated by this: the "failure" it declined to chase was largely not
there.*

### 3b.2 The rest of the C2 diagnostics

    FIELD floor 0.02 +- 0.02 nat [emp_q95, zero-fraction 0.43]   CRITERION floor 38.59 +- 1.72 [gumbel, 0.00]
    comb spacing: FIELD 0.0009 kpc  vs  CRITERION 0.0001 kpc  (9x coarser)
    polish improved 100 % of signal and null reals; median gain +1700 2F, max +7852; worst change +0.164
    single-fringe EXCLUDED: FIELD 12/13920 signal, 15/11600 null | CRITERION 0 and 0
    truth outside the field's 512-pt window: 8/13920

---

## 5. ADD-LIST (staged with `git add`; **Matt commits — I have committed nothing**)

| path | what it is |
|---|---|
| `BASELINE_field.md` | this report (pre-registration §0 written and staged **before** the fan) |
| `hpc_harbor/baseline/baseline.py` | the driver: venue, F_e detector + polish, the field's comb/template/E-step, the criterion arm, the gate chain |
| `hpc_harbor/baseline/bl_score.py` | the CPU scorer: matched-FAP thresholds, criterion-v2.2 floors with zero-fraction columns, Wilson intervals, the table and the panel |
| `hpc_harbor/baseline/bl_smoke.sbatch` | the per-arm smoke (gate chain + measured per-realisation wall) |
| `hpc_harbor/baseline/bl_fan.sbatch` | the fan (idempotent, resumable, hardlink-seeded compile cache) |
| `GLACIER_capstone.md` | **appended only**: the BASELINE lane claim, its measured amendment, and the note to PHOENIX about the live 40GB lane |
| `reports/BASELINE_panel.png` | **the B3 panel** — A: detection at matched FAP; B: purity vs the calibrated criterion |
| `reports/BASELINE_tables.txt` | **the B3 tables** (JSON: thresholds, floors + zero-fractions, purities, Wilson intervals) |
| `reports/baseline_score.npz` | raw scored columns *(gitignored? no — `reports/` is tracked; 12 kB)* |
| `hpc_harbor/baseline/bl_launch.sh` | the chained launcher (refuses to fire without a warm cache) |

**NOT MINE, staged concurrently by PHOENIX** (it shares this working tree):
`hpc_harbor/summit/smt_sighting_baseline.py`, `hpc_harbor/frozen/*`, `hpc_harbor/ledger/*`,
`SUMMIT_closure.md`, `LEDGER_stats_audit.md`. I did not write or stage those; flagged so the
commit is not mis-attributed.

**Not staged, deliberately:** `BASELINE_results/` — `*_results/` is gitignored, which is the
house discipline for campaign banks (`glacier-precommit-audit`).

**Nothing outside `BASELINE_results/` was written.** `SURFACE_results/`,
`GENERALISE_results/`, `GLACIER*/` and `SPARK*/` were opened read-only. No tracked file was
edited except the append-only additions to `GLACIER_capstone.md`. `TE.seed_scan`,
`spark.py`, `chorus.py`, `generalise.py` and `surface.py` are **untouched** — the
`H_ABSENT` finding in §1.2 is reported, not patched upstream.

---

## 6. SPEND LEDGER, VERDICT, AND STOP

### 6.1 Spend

| item | jobs | wall |
|---|---|---|
| smokes (C1 ×3, C2 ×2) | `12833821`, `12834517`, `12834778`, `12834510`, `12834976` | 3292 + 1697 + 1192 + 3186 + 1342 s = **2.7 GPU-hr** |
| C1 fan (130 reals) | `12834949_[0-3]` | 4484 + 4628 + 4659 + 4869 s = **5.2 GPU-hr** |
| C2 fan (128 reals) | `12835126_[0-3]` | 4589 + 4655 + 4655 + 4706 s = **5.2 GPU-hr** |
| **TOTAL** | | **≈ 13.1 GPU-hr** |

**This is over the pre-registered 12 GPU-hr cap (§0.5b), by ~1.1 GPU-hr (9 %), and the
overrun is stated rather than absorbed.** The cause is named: the cap was sized on
`258 × 105 s + venue builds` and did **not** include the **three extra smoke cycles** the
gate chain forced (defects 1–4, ~1.5 GPU-hr) or the per-task first-E-step compile (~8 min ×
8 tasks ≈ 1.1 GPU-hr, §1.6). *The realisation budget itself was accurate — 105 s measured
against 109 s delivered.* No trim was invoked and none was needed for the science: all 258
realisations, both oracle bounds, all 8 C2 skies.

Lane: `interactive_gpu` A100-40GB on dgx01, `%4`, never dgx03. One self-reported `%4 → 5`
overlap for ~20 min (`hpc_harbor/LANE_CLAIM.md`).

### 6.2 THE VERDICT

> **At a matched false-alarm probability, the field's standard Earth-term F_e-statistic
> detects at least as well as `criterion-v2.2` — better at the census cell (1.00 vs 0.67),
> equal at the A-SKY survivor (1.00 vs 1.00). The criterion buys nothing on detection.**
>
> **What it buys is on the distance question, and the size of it is the result. The field's
> standard per-pulsar fringe inference states ~9.3–9.5 confident distance identifications
> per realisation of which 92–93 % are WRONG (purity 0.08 / 0.07). The criterion states 1.07
> (C1) and 20.02 (C2) at purity 0.78 / 0.96. Wrong claims per realisation: field 8.64–8.80,
> criterion 0.23–0.78.**
>
> **The field is not failing to find the truth — both arms recover comparable CORRECT counts
> at C1 (0.73 vs 0.83). It is failing to decline the falsehoods.**
>
> **The mechanism is not the search.** With the field's sky, frequency and chirp mass set to
> the injected member's TRUE values, purity only reaches **0.31** (C1) and **0.15** (C2) —
> and it gets *worse* with more sources, which is the signature of **a one-source fringe comb
> against a many-source sky**, not of estimation error. **Nor is it a missing threshold:**
> bolting the criterion's own `lnK` + floor onto the field's comb recovers 0.08 → 0.15 only.

**Pre-registered readouts, both cells:** R1 **confirmed** (purity ≪ 0.9, refutation excluded);
R2 **confirmed** (gap 0.70 / 0.89, disjoint intervals); R3 **the surprise** — the floor is not
the story, the comb is; R4 **cell-dependent**, criterion costs detection only near onset;
R5 **confirmed and strengthened** — search error explains a minority and the residual grows
with the population.

### 6.3 STOP — for Matt

**Nothing here arms a protocol step, moves a banked verdict, or enters a closure claim.**
Three items are surfaced rather than decided:

1. **The budget overran 12 → 13.1 GPU-hr (§6.1).** Stated, not absorbed. The overrun is
   entirely gate-chain and compile overhead, not science scope.
2. **G-B3b remains formally FAILED on C2's gate draw** even though §3b.1 shows the production
   ensemble refines. **I did not retro-fit the gate to pass.** The two-stage polish that would
   fix it properly is un-run, and running it would require re-running BOTH arms to keep one
   estimator (§1.7d) — that is a re-scope decision, not mine.
3. **Scope.** Two cells, one venue (T = 30, `lit`, the 116-pulsar MOCK array), loudness far
   above NG15 and labelled only by its parameters. **The verdict is about these cells.**
   Whether the 92 % field wrong-rate holds at census loudness, at other baselines, or on a
   real array is untested here.

**Natural follow-ups, none started:** the same comparison at a faint/marginal cell (where the
criterion's detection cost at C1 suggests the two methods may cross); a multi-source field
baseline (the obvious repair — does a 3-source comb recover the purity?); and the two-stage
polish as a standalone machinery fix.
