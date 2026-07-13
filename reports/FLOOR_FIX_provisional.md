# THE FLOOR FIX — PROVISIONAL, PENDING RE-SCORE

> ## ⛔ SUPERSEDED — 2026-07-13, BY `reports/RECUT_floors.md`.
>
> **The re-score landed. Every bound in this file is now resolved to an exact number, and this file
> is kept only for its trail.** Read `RECUT_floors.md` instead. Where the two disagree, RECUT wins.
>
> **How the provisional analysis scored, judged against the re-cut it demanded:**
>
> | this file said | the re-cut found | verdict |
> |---|---|---|
> | 57 ≤ N_onset ≤ 67; **do not quote 59** | **N_onset = 59** (2 died, 2 born) | **correct as a bound, loose as a forecast** — and the caution was right for the right reason: the *map* moved even though the *total* did not |
> | 2 onset cells AT RISK, bounds ≤ 1.20 / ≤ 1.37, **INDETERMINATE** | both **RETRACTED** — 0.77 and 0.90 | **correct, and conservative** — neither bound was reached |
> | 8 cells MAY IGNITE | **2 ignited**; the other 6 rose and stayed under | correct as a bound |
> | `h*` unbounded below (7/18): **SUSPENDED** | **REINSTATED at 7/18 — different seven** | suspension was right; the number survived, the names did not |
> | e = 0.3 switch-on: **INDETERMINATE, do not quote** | **REFUTED** (lit 0.70) / MARGINAL (vlbi 1.03) | **the embargo was correct and load-bearing** — the claim did not survive |
> | e = 0.5 is the defensible binding | **CONFIRMED, both tiers** | correct |
> | §5: `anchor_ladder.npz` offenders are transposed and **"nothing in the file says so"** | **half wrong — the file DID say so** (`offender_index`, a 48-element key, already encoded the exact permutation §5 reverse-engineered). The trap is real but narrower: it is sprung only by a re-cut that *ignores* the key | **the diagnosis was half wrong, and RECUT found it by checking rather than assuming.** The array is now permuted into metadata-row order and the naive re-cut is the correct one |
>
> **The one thing it got wrong, kept visible on purpose:** §5's recommended fix — *"or add an
> explicit index key"* — **asked for a key the file already had.**

**Session:** cronus doc/machinery session, 2026-07-13. **Status: SUPERSEDED (was PROVISIONAL).**
The floors below are final. **No count below is corrected, and none may be quoted as one** — the
corrected counts are in `RECUT_floors.md`.
Script: `CW_transition/floorfix_provisional.py` (CPU, seconds, no GPU).
Bank: `reports/floorfix_provisional.npz`. Every number here is read back off that bank.

---

## 0. THE CONVENTION, ADOPTED

Adopted verbatim from **ANCHOR §4** (`reports/ANCHOR_realdata_null.md`):

> **The D2 Gumbel floor is valid only where the nullN zero-fraction is ≲ 20 %. Above that, quote
> the empirical (1−α) quantile with a bootstrap error, and bank the zero-fraction beside it. The
> zero-fraction is a REQUIRED column, not a caveat.**

The offender statistic is `0.0` whenever a realisation has no cell passing layer 1 ⊕ layer 3. At
faint `h` that is *most* realisations. A Gumbel fitted to a point mass at zero is dragged **down**
toward it, understating the α = 0.05 bar by up to **2.8×** — and it errs **permissive**: a floor
that is too low lets pure-noise offenders through. This is the dangerous direction.

## 1. WHAT THIS SESSION COULD AND COULD NOT DO — the scope of inference

The floor is re-derivable from the banked summaries. **The count is not.** The per-realisation
signal `dlnL` columns live in `SURFACE_results/` and `CHORUS_results/` **on ACCRE** and are not in
this repo. So:

- every **FLOOR** below is final;
- every **COUNT** below is the one banked under the **old** floor, quoted **only to BOUND** the
  corrected count using the one thing that is certain — **the count is monotone non-increasing in
  the floor.** Floor up ⇒ count can only fall. Floor down ⇒ count can only rise.

**Nothing here interpolates a count to a new floor.** That would be a fabricated number, and it is
what the artifact-readback convention (§5) exists to forbid.

---

## 2. SURFACE — the blast radius is SMALL, and that is the finding

`surface_analysis.npz` banks `fN_zerofrac` and `fN_emp` per cell, so the corrected floor **value**
is recoverable for all 108 cells without the raw banks.

**The convention only replaces the floor where the zero-fraction exceeds 20 %. That is 15 of 108
cells. Of the 59 ONSET cells, it touches exactly 2.**

| | cells |
|---|---|
| banked verdicts (pre-fix) | **ONSET 59** · MARGINAL 4 · below 45 |
| zero-fraction > 20 % → Gumbel INVALID, floor := empirical q95 | **15 / 108** |
| zero-fraction ≤ 20 % → Gumbel valid, **floor and count stand untouched** | **93 / 108** |
| **ONSET cells UNTOUCHED → stand, definitively** | **57** |
| ONSET cells whose floor RISES → count can only fall → **AT RISK** | **2** |
| non-onset cells whose floor FALLS → count can only rise → **MAY IGNITE** | **8** |

> ### **The corrected onset count is BOUNDED: 57 ≤ N_onset ≤ 67.** (Banked, pre-fix: 59.)
> It cannot be pinned without the re-score. **Do not quote 59.**

**The pre-registered expectation is confirmed on the loud cells, definitively:** loud-cell onsets
have low zero-fractions, so their Gumbel is valid and neither floor nor count moves. **57 of the 59
onsets are untouched by the floor fix.** The exposure is entirely at faint `h`.

### The 15 touched cells — the entire blast radius

| h | T | tier | struct | zero-f | Gumbel | EMP q95 | ratio | count | @fl+sd | dir | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| −13.25 | 30 | lit | 2+14 | 0.78 | 3.08 | 6.36 | 2.07× | 0.40 | 0.40 | **UP** | below |
| −13.25 | 50 | lit | 2+14 | 0.73 | 5.33 | 8.89 | 1.67× | 1.43 | 1.37 | **UP** | **ONSET** |
| −13.25 | 40 | lit | 2+14 | 0.72 | 5.07 | 7.26 | 1.43× | 0.63 | 0.60 | **UP** | below |
| −13.25 | 30 | lit | 3+13 | 0.67 | 4.76 | 6.48 | 1.36× | 0.67 | 0.67 | **UP** | below |
| −13.25 | 40 | lit | 3+13 | 0.60 | 6.81 | 9.39 | 1.38× | 1.37 | 1.20 | **UP** | **ONSET** |
| −13.25 | 30 | vlbi | 2+14 | 0.60 | 5.33 | 5.72 | 1.07× | 0.40 | 0.27 | **UP** | below |
| −13.00 | 30 | lit | 2+14 | 0.46 | 12.57 | 12.73 | 1.01× | 0.13 | 0.10 | **UP** | below |
| −13.25 | 40 | vlbi | 2+14 | 0.42 | 8.25 | 7.11 | 0.86× | 0.37 | 0.27 | down | below |
| −13.25 | 30 | vlbi | 3+13 | 0.41 | 8.16 | 7.98 | 0.98× | 0.47 | 0.43 | down | below |
| −13.25 | 50 | lit | 3+13 | 0.41 | 12.06 | 11.46 | 0.95× | 1.10 | 0.87 | down | MARGINAL |
| −13.25 | 50 | vlbi | 2+14 | 0.33 | 9.30 | 8.05 | 0.87× | 0.90 | 0.63 | down | below |
| −13.00 | 40 | lit | 2+14 | 0.28 | 17.18 | 15.41 | 0.90× | 0.20 | 0.20 | down | below |
| −13.25 | 40 | vlbi | 3+13 | 0.27 | 10.24 | 8.92 | 0.87× | 0.63 | 0.60 | down | below |
| −13.00 | 30 | lit | 3+13 | 0.27 | 19.46 | 16.60 | 0.85× | 0.37 | 0.23 | down | below |
| −13.25 | 30 | lit | 5+11 | 0.21 | 20.17 | 16.69 | 0.83× | 0.60 | 0.47 | down | below |

### The two ONSET cells at risk

| cell | floor | banked count | corrected count | status |
|---|---|---|---|---|
| h = −13.25, T = 40, **lit**, 3+13 | 6.81 → **9.39** (1.38×) | 1.37 | **≤ 1.20** | **INDETERMINATE** |
| h = −13.25, T = 50, **lit**, 2+14 | 5.33 → **8.89** (1.67×) | 1.43 | **≤ 1.37** | **INDETERMINATE** |

Both are faint-`h`; both have rising floors, so their counts can only fall; neither bound reaches
≤ 1, so **neither is retracted and neither is confirmed.** They need the exact re-cut.

**These two are frontier cells** — they carry part of SURFACE §4's *"`h*` is not bounded below in
7 of 18 columns"* (‡). **That claim is therefore SUSPENDED pending the re-score**, not retracted.

**SURFACE banked no raw offender vectors**, so its empirical floors above carry **no bootstrap
error** — which the convention requires. ACCRE must bank the `nullN` offenders.

---

## 3. CHORUS — every floor in the campaign is invalid

`ch_floors.npz` **did** bank its raw offenders, so CHORUS's floors are **fully corrected here, with
bootstrap errors**.

> ### **All 26 of 26 CHORUS cells have zero-fraction > 20 %. Every CHORUS floor is invalid under the adopted convention. 23 of 26 rise, by up to 1.96×.**

Selected rows (full table in the bank):

| cell | zero-f | Gumbel ± fit | **ADOPTED** ± bs | ratio | count (old floor) |
|---|---|---|---|---|---|
| m0 e=0 lit | 0.73 | 4.31 ± 0.37 | **7.00 ± 1.10** | 1.63× | 0.70 |
| m0 e=0 vlbi | 0.45 | 7.59 ± 0.48 | **7.06 ± 0.40** | 0.93× | 0.43 |
| **m1 e=0.3 lit** | **0.73** | **7.39 ± 0.63** | **11.30 ± 1.02** | **1.53×** | **1.57** |
| **m1 e=0.3 vlbi** | **0.48** | **10.58 ± 0.83** | **10.78 ± 1.54** | **1.02×** | **1.13** |
| m1 e=0.5 lit | 0.81 | 4.61 ± 0.41 | **8.58 ± 0.87** | 1.86× | 3.60 |
| m1 e=0.5 vlbi | 0.51 | 9.55 ± 0.76 | **9.87 ± 0.91** | 1.03× | 2.43 |
| m1 e=0.7 lit | 0.69 | 8.51 ± 0.72 | **11.65 ± 0.76** | 1.37× | 7.83 |
| m2 e=0.7 vlbi | 0.64 | 11.88 ± 0.99 | **16.82 ± 3.17** | 1.42× | 6.73 |
| m3 eU lit | 0.80 | 5.77 ± 0.51 | **11.33 ± 1.10** | 1.96× | 5.23 |

### 3.1 THE e = 0.3 SWITCH-ON IS INDETERMINATE — and this governs what is said outside

The **e = 0.3 switch-on is the single cell the floor fix hits hardest relative to its margin.**

- **lit:** floor **7.39 → 11.30 nat (+53 %, 6.2σ of its own Gumbel fit error)**, against a banked
  count of **1.57** — a margin of just **+0.57** over the >1 bar.
- **vlbi:** floor 10.58 → 10.78, against a banked count of **1.13** — a margin of **+0.13**.

The floor rises, so **the count can only fall**, and it is falling toward a bar it sits barely
above. **Whether "one e = 0.3 member switches the count on" survives cannot be decided from this
repo.** It is not refuted. It is **not established**.

> ## EXTERNAL-QUOTE GUIDANCE (binding until the re-score lands)
>
> **Do NOT quote the e = 0.3 switch-on.** *"A single e = 0.3 member suffices"* rests on counts of
> 1.57 / 1.13 against a floor that has just risen 53 % / 2 %, and it may not survive.
>
> **The currently defensible switch number is e = 0.5.** Its counts (3.60 lit / 2.43 vlbi) clear
> the bar by **+2.60 / +1.43** against floor rises of 1.86× / 1.03×. Even a floor rise of that size
> cannot plausibly carry a count of 3.60 below 1.0 — but note this too is an argument from margin,
> **not a re-cut**, and it is stated as such. The e = 0.7 result (7.83 / 6.13) is safer still.
>
> CHORUS's *qualitative* headline is untouched: **the clock is NOT shared** (0 of ~120 lifted
> certifications are circular-attributed) and **eccentric structure transforms the count.** Neither
> depends on the floor. It is the *onset threshold in e* that is in play.

---

## 4. IGNITE-2 — one of the programme's two calibrated cells needs restating

ANCHOR §4 states *"IGNITE-2's two calibrated cells are safe."* **One word needs changing.**

| cell (STORY S6.2.1) | zero-frac | Gumbel | EMP q95 | status |
|---|---|---|---|---|
| (−12.75, 30, lit) | **0.00** | 38.86 ± 1.46 | 37.18 ± 3.32 | **GUMBEL VALID** — the estimator holds; floor and count (0.92) stand as banked |
| (−13.25, 30, vlbi) | **0.45** | 7.59 ± 0.48 | **7.06 ± 0.40** | **GUMBEL INVALID** — floor must be **RESTATED as the empirical quantile** |

The vlbi cell's count (0.54) is safe **only because its corrected floor FALLS** (count can only
rise) — **not because the estimator is sound there.** Its zero-fraction of 0.45 is over twice the
validity gate. **The floor VALUE must be restated as 7.06 ± 0.40 nat (empirical q95, bootstrap),
not 7.59 ± 0.48 (Gumbel).** Both cells' *verdicts* (both BELOW onset) are unchanged.

---

## 5. A BANK DEFECT, FOUND BY THE ARTIFACT-READBACK CONVENTION ON ITS FIRST USE

Re-deriving ANCHOR's ladder from its own raw offenders **did not reproduce ANCHOR's own published
numbers** — max deviation **79.9 nat**. The cause:

> **`anchor_ladder.npz` stores `offenders` in (rung-major, cell-minor) order, while every sibling
> column — `rung`, `h`, `tier`, `zero_frac`, `emp_q95` — is (cell-major, rung-minor). The array is
> the TRANSPOSE of its own metadata, and nothing in the file says so.**

Applying the permutation `perm[j] = (j % 6)·8 + (j // 6)` reproduces the banked `emp_q95` **and**
`zero_frac` to **exactly 0.000e+00**. So:

- **ANCHOR's published §3 ladder table is CORRECT.** The verdicts stand. Nothing is retracted.
- **The npz layout is defective.** ANCHOR §9 advertises *"every floor is re-cuttable from the bank
  without rerunning a GPU"* — it is, but a naive re-cut marries `offenders[j]` to the **wrong
  cell's label**, and the error is **silent**, because both are plausible floors.

**FIX ON ACCRE:** re-bank `offenders` in metadata-row order, or add an explicit index key. Until
then any re-cut from this file must apply the permutation, which
`CW_transition/floorfix_provisional.py` asserts rather than assumes.

*This is what the convention is for. The defect was invisible from the report, invisible from the
column names, and fatal to exactly the workflow the report recommends.*

---

## 6. WHAT ACCRE MUST SUPPLY (CPU-only; no GPU; the banks already exist)

1. **Re-score the counts** under adopted floors (empirical q95 + bootstrap where zero-frac > 20 %)
   for the **15 touched SURFACE cells** and **all 26 CHORUS cells**, from the existing
   per-realisation banks. This settles: the two at-risk onset cells, the eight that may ignite, the
   corrected onset count, the `h*`-unbounded-below claim, and the **e = 0.3 verdict**.
2. **Bank SURFACE's raw `nullN` offender vectors** — it is the only campaign that did not, and
   without them no bootstrap error exists on its empirical floors.
3. **Re-bank `anchor_ladder.npz`'s `offenders`** in metadata-row order (§5).
4. **Generate the canonical `b1_L_gwb.npz`** at `--cpus-per-task=8` and validate it against
   ANCHOR's g1 replay (80 banked `ig_nullN_*`, bit-identical). See `trackB_lgwb_gate.py`.

**Held pending that re-score:** the corrected onset-surface table; the STORY supersessions that
depend on it; every tag flip; `criterion-v2.2`.
