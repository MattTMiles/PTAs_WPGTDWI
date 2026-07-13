# THE FLOOR FIX — RE-CUT AND SETTLED

**Session:** RECUT, ACCRE, 2026-07-13. **Supersedes the PROVISIONAL status of
`reports/FLOOR_FIX_provisional.md`.** Every count below is **re-cut from the per-realisation
signal banks** against the adopted floor — none is bounded, interpolated, or inherited.

Scripts (CPU, no GPU, no new realisations): `CW_transition/recut_surface.py`,
`recut_chorus.py`, `bank_surface_offenders.py`, `rebank_anchor_ladder.py`, `make_lgwb_bank.py`.
Banks: `reports/recut_surface.npz`, `recut_chorus.npz`, `surface_nullN_offenders.npz`,
`anchor_ladder.npz` (re-banked), `b1_L_gwb_manifest.npz`.

---

## 0. THE CONVENTION, AND THE GATES THAT LICENSE THESE NUMBERS

Adopted verbatim from **ANCHOR §4**:

> **The D2 Gumbel floor is valid only where the nullN zero-fraction is ≲ 20 %. Above that, quote
> the empirical (1−α) quantile with a bootstrap error, and bank the zero-fraction beside it. The
> zero-fraction is a REQUIRED column, not a caveat.**

The onset test is **unchanged** — `ONSET` iff the count at **floor + its own error** exceeds 1.
Only the floor and the error change: above the gate the error is now the **bootstrap sd of the
empirical quantile**, so the onset test is made against a bootstrap error exactly as required.

**No corrected number was emitted until an uncorrected one was reproduced from the same raw
columns.** Both re-scores pass two readback gates:

| gate | what it proves | SURFACE | CHORUS |
|---|---|---|---|
| **A — floors** | recomputed offenders reproduce the *banked* Gumbel floor, sd, emp-q95, zero-fraction | **0.000e+00** | **0.000e+00** |
| **B — counts** | this scorer at the *old* floors reproduces the *banked* counts and verdicts | **0.000e+00** (108/108 verdicts) | **0.000e+00** |

Gate B is the load-bearing one: it proves the scorer used here **is** the one that produced the
published surface. Without it, a corrected count is just a different number.

---

## 1. SURFACE — THE CORRECTED ONSET COUNT IS **59**

> ### **N_onset = 59.** The provisional bound (57 ≤ N ≤ 67) is resolved.
> **The pre-fix banked count was also 59 — and that is a coincidence, not a confirmation.
> Two onsets died and two were born. The number is stable; the map is not.**

| | cells |
|---|---|
| zero-fraction > 20 % → Gumbel invalid, floor := empirical q95 ± bootstrap | **15 / 108** |
| zero-fraction ≤ 20 % → Gumbel valid, floor and count stand untouched | **93 / 108** |
| ONSET, untouched → stand definitively | **57** |
| ONSET **lost** (floor rose, count fell through the bar) | **2** |
| ONSET **gained** (floor fell, count rose through the bar) | **2** |
| **N_onset, re-cut** | **59** (MARGINAL 3, below 46) |

Of the 8 below-onset cells whose floors fall, **only 2 actually ignited**; the other 6 rose but
stayed under the bar. The provisional upper bound of 67 was correct as a bound and loose as a
forecast.

### 1.1 The four cells that flipped — named

| cell | zero-f | floor (Gumbel → adopted) | count (old → **re-cut**) | @floor+bs | was → **now** |
|---|---|---|---|---|---|
| h=−13.25, T=40, **lit**, 3+13 | 0.60 | 6.81 → **9.39 ± 0.79** (1.38×) | 1.37 → **0.77** | 0.67 | ONSET → **below** |
| h=−13.25, T=50, **lit**, 2+14 | 0.73 | 5.33 → **8.89 ± 1.17** (1.67×) | 1.43 → **0.90** | 0.67 | ONSET → **below** |
| h=−13.25, T=50, **lit**, 3+13 | 0.41 | 12.06 → **11.46 ± 0.51** (0.95×) | 1.10 → **1.37** | 1.10 | MARGINAL → **ONSET** |
| h=−13.25, T=50, **vlbi**, 2+14 | 0.33 | 9.30 → **8.05 ± 0.56** (0.87×) | 0.90 → **1.17** | 1.10 | below → **ONSET** |

**The two at-risk cells named in FLOOR_FIX §2 are both RETRACTED.** Neither bound (≤1.20, ≤1.37)
was reached: the true re-cut counts are **0.77** and **0.90**, comfortably below the bar. The
floor fix does bite, and it bites exactly where the provisional analysis said the exposure was —
faint `h`, high zero-fraction.

### 1.2 The 15 touched cells — the entire blast radius

| h | T | tier | struct | zero-f | Gumbel | **EMP q95 ± bs** | ratio | old | **re-cut** | @fl+bs | was → now |
|---|---|---|---|---|---|---|---|---|---|---|---|
| −13.25 | 30 | lit | 2+14 | 0.78 | 3.08 | **6.36 ± 1.23** | 2.07× | 0.40 | 0.23 | 0.13 | below → below |
| −13.25 | 50 | lit | 2+14 | 0.73 | 5.33 | **8.89 ± 1.17** | 1.67× | 1.43 | **0.90** | 0.67 | **ONSET → below** |
| −13.25 | 40 | lit | 2+14 | 0.72 | 5.07 | **7.26 ± 1.13** | 1.43× | 0.63 | 0.37 | 0.37 | below → below |
| −13.25 | 30 | lit | 3+13 | 0.67 | 4.76 | **6.48 ± 0.67** | 1.36× | 0.67 | 0.57 | 0.50 | below → below |
| −13.25 | 30 | vlbi | 2+14 | 0.60 | 5.33 | **5.72 ± 0.80** | 1.07× | 0.40 | 0.30 | 0.13 | below → below |
| −13.25 | 40 | lit | 3+13 | 0.60 | 6.81 | **9.39 ± 0.79** | 1.38× | 1.37 | **0.77** | 0.67 | **ONSET → below** |
| −13.00 | 30 | lit | 2+14 | 0.46 | 12.57 | **12.73 ± 1.06** | 1.01× | 0.13 | 0.13 | 0.10 | below → below |
| −13.25 | 40 | vlbi | 2+14 | 0.42 | 8.25 | **7.11 ± 0.42** | 0.86× | 0.37 | 0.53 | 0.50 | below → below |
| −13.25 | 30 | vlbi | 3+13 | 0.41 | 8.16 | **7.98 ± 0.54** | 0.98× | 0.47 | 0.47 | 0.43 | below → below |
| −13.25 | 50 | lit | 3+13 | 0.41 | 12.06 | **11.46 ± 0.51** | 0.95× | 1.10 | **1.37** | 1.10 | **MARGINAL → ONSET** |
| −13.25 | 50 | vlbi | 2+14 | 0.33 | 9.30 | **8.05 ± 0.56** | 0.87× | 0.90 | **1.17** | 1.10 | **below → ONSET** |
| −13.00 | 40 | lit | 2+14 | 0.28 | 17.18 | **15.41 ± 1.60** | 0.90× | 0.20 | 0.30 | 0.20 | below → below |
| −13.00 | 30 | lit | 3+13 | 0.27 | 19.46 | **16.60 ± 1.60** | 0.85× | 0.37 | 0.60 | 0.47 | below → below |
| −13.25 | 40 | vlbi | 3+13 | 0.27 | 10.24 | **8.92 ± 0.64** | 0.87× | 0.63 | 1.00 | 0.87 | below → below |
| −13.25 | 30 | lit | 5+11 | 0.21 | 20.17 | **16.69 ± 1.98** | 0.83× | 0.60 | 0.87 | 0.67 | below → below |

The 93 untouched cells — **including 57 of the 59 onsets** — keep their banked floor and count
exactly. SURFACE's loud-cell result is unmoved, as pre-registered.

### 1.3 The `h*`-unbounded-below suspension: **REINSTATED, with its membership changed**

SURFACE §4 ‡ / §11.5 claimed `h*` is not bounded below in **7 of 18** frontier columns (a column
= tier × structure × T; `h*` is unbounded when the faintest grid strain, h = −13.25, is ONSET).

> ### **RE-CUT: 7 of 18. The claim stands as published — but not on the same 7 columns.**

| | columns |
|---|---|
| lost the faint edge | lit 3+13 T=40 · lit 2+14 T=50 |
| gained the faint edge | lit 3+13 T=50 · vlbi 2+14 T=50 |
| unchanged (5) | lit 5+11 T=40, T=50 · vlbi 5+11 T=30, T=40, T=50 |

**The suspension is lifted and the claim is reinstated at its published strength (7/18).** The
count is coincidentally identical; the identity of the columns is not. Quote the number, but do
not quote the *same seven columns* — the figure and any per-column text must be regenerated from
`recut_surface.npz`.

**One caveat, stated.** The newly-gained `vlbi 2+14 T=50` column is an **isolated** faint-edge
onset: h = −13.25 certifies (1.17) while every louder h in that column sits below the bar
(0.33–0.90). It is a frontier column by the letter of the definition, but there is no contiguous
frontier running down to it. This is not a re-cut artifact — the floor grows as ≈ h^1.5–2, so a
faint cell can out-certify a loud one, and the *published* surface already showed this pattern in
`lit 2+14 T=50`. It is worth a sentence in SURFACE, not a retraction.

---

## 2. CHORUS — THE e = 0.3 SWITCH-ON **DOES NOT SURVIVE**

All 26 of 26 CHORUS cells sit above the validity gate (zero-fraction 0.33–0.81), so **every
CHORUS floor is restated** as the empirical q95 with a bootstrap error. 23 of 26 floors rise.

> ### **THE HEADLINE: "a single e = 0.3 member switches the count on" is REFUTED.**
> **lit collapses from 1.57 to 0.70 — below the bar. vlbi survives at 1.03, but only marginally:
> it does not clear its own floor's bootstrap error (0.60).**

| the switch-on cells | zero-f | floor (Gumbel → **adopted**) | ratio | count (old → **re-cut**) | @floor+bs | margin vs >1 | verdict |
|---|---|---|---|---|---|---|---|
| **m1 e=0.3 lit** | 0.73 | 7.39 ± 0.63 → **11.30 ± 1.02** | 1.53× (6.2σ) | 1.57 → **0.70** | 0.43 | **−0.30** | **below — REFUTED** |
| **m1 e=0.3 vlbi** | 0.48 | 10.58 ± 0.83 → **10.78 ± 1.54** | 1.02× (0.2σ) | 1.13 → **1.03** | 0.60 | **+0.03** | **MARGINAL — not confirmed** |
| m1 e=0.5 lit | 0.81 | 4.61 ± 0.41 → **8.58 ± 0.87** | 1.86× | 3.60 → **3.13** | 2.70 | +2.13 | **CONFIRMED** |
| m1 e=0.5 vlbi | 0.51 | 9.55 ± 0.76 → **9.87 ± 0.91** | 1.03× | 2.43 → **2.27** | 1.73 | +1.27 | **CONFIRMED** |

The provisional analysis called this **indeterminate** and forbade quoting it. The re-cut settles
it: **the count falls through the bar in `lit`, and in `vlbi` it clears the bar by +0.03 while
failing at floor + bootstrap error.** One eccentric member at e = 0.3 does **not** switch the
count on.

### 2.1 THE CORRECTED SWITCH-ON THRESHOLD IN e — the externally quotable number

> ### **With ONE eccentric member: the switch-on is at e = 0.5.**
> ### **With TWO OR MORE: it is at e = 0.3 (confirmed, both tiers).**
> **The docs' interim binding to e = 0.5 STANDS — and it was the right call.**

The threshold is **not a property of eccentricity alone**; it depends on how many eccentric
members carry it. That is the corrected, quotable statement.

| n_ecc | tier | e = 0.3 | e = 0.5 | e = 0.7 | switch-on |
|---|---|---|---|---|---|
| **1** | lit | 0.70 [0.43] **below** | 3.13 [2.70] CONFIRMED | 5.43 [4.90] CONFIRMED | **e = 0.5** |
| **1** | vlbi | 1.03 [0.60] MARGINAL | 2.27 [1.73] CONFIRMED | 5.77 [5.13] CONFIRMED | **e = 0.5** (0.3 not confirmed) |
| 2 | lit | 2.77 [2.07] CONFIRMED | 4.90 [4.63] CONFIRMED | 5.47 [4.27] CONFIRMED | **e = 0.3** |
| 2 | vlbi | 1.77 [1.43] CONFIRMED | 3.97 [3.20] CONFIRMED | 4.10 [2.73] CONFIRMED | **e = 0.3** |
| 3 | lit | 2.50 [2.17] CONFIRMED | 5.83 [4.97] CONFIRMED | 4.07 [3.60] CONFIRMED | **e = 0.3** |
| 3 | vlbi | 2.20 [1.93] CONFIRMED | 4.50 [3.97] CONFIRMED | 5.07 [3.93] CONFIRMED | **e = 0.3** |

*(count at the adopted floor, [count at floor + bootstrap error]. CONFIRMED = clears the bar at
floor + bs; MARGINAL = clears the bar only at the floor.)*

### 2.2 A published CHORUS sentence that is now FALSE

CHORUS §5: *"**Every** eccentric mix clears the >1 bar."* **This is no longer true.** `m1 e=0.3
lit` re-cuts to **0.70**. The sentence must be struck or amended to *"every eccentric mix clears
the bar except a single e = 0.3 member in the lit tier."*

### 2.3 What CHORUS keeps — and one more thing that does NOT survive

**Keeps:**
- **the clock is NOT shared** (0 of ~120 lifted certifications circular-attributed) — a
  *structural* result, independent of the floor;
- **eccentric structure transforms the count** (m3 eU lit: 3.40; m2 eU vlbi: 7.03 — CONFIRMED);
- **the 10× lever** survives comfortably: m0 e=0 → m1 e=0.7 is **14.8× (lit)** and **12.4×
  (vlbi)** on the re-cut counts (was 11.2× / 14.2×).

**Does NOT survive — a second casualty, found by checking rather than assuming:**

> **The "trade inverts at n_ecc = 3" claim is NOT stable under the re-cut.** Whether the count
> falls from n_ecc = 2 → 3 flips status in **3 of the 8 (e, tier) combinations**:

| e, tier | old m1→m2→m3 | inverts? | re-cut m1→m2→m3 | inverts? |
|---|---|---|---|---|
| e=0.3 lit | 1.57 · 3.27 · 4.57 | no | 0.70 · 2.77 · **2.50** | **YES** ← changed |
| e=0.5 lit | 3.60 · 7.33 · 7.00 | YES | 3.13 · 4.90 · **5.83** | **no** ← changed |
| e=0.7 vlbi | 6.13 · 6.73 · 6.67 | YES | 5.77 · 4.10 · **5.07** | **no** ← changed |
| e=0.7 lit · eU lit · eU vlbi | — | YES | — | YES (unchanged) |
| e=0.3 vlbi · e=0.5 vlbi | — | no | — | no (unchanged) |

Old: 5 of 8 inverted. Re-cut: 4 of 8 — **but not the same four.** The inversion is a difference
between two counts that both moved, and it is not robust to the floor correction. **Do not quote
"the trade inverts at n_ecc = 3" without re-deriving it from `recut_chorus.npz`.** It is not
refuted; it is no longer a clean claim.

The *onset threshold in e* and the *trade inversion* are what moved. Everything structural held.

### 2.4 All 26 cells, re-cut

| cell | zero-f | Gumbel ± fit | **ADOPTED ± bs** | ratio | old | **re-cut** | @fl+bs | wrong | grade |
|---|---|---|---|---|---|---|---|---|---|
| 0,e00,lit | 0.73 | 4.31 ± 0.37 | **7.00 ± 1.10** | 1.63× | 0.70 | 0.37 | 0.30 | 0.00 | below |
| 0,e00,vlbi | 0.45 | 7.59 ± 0.48 | **7.06 ± 0.40** | 0.93× | 0.43 | 0.47 | 0.43 | 0.10 | below |
| **1,e03,lit** | **0.73** | 7.39 ± 0.63 | **11.30 ± 1.02** | 1.53× | 1.57 | **0.70** | 0.43 | 0.00 | **below** |
| **1,e03,vlbi** | **0.48** | 10.58 ± 0.83 | **10.78 ± 1.54** | 1.02× | 1.13 | **1.03** | 0.60 | 0.00 | **MARGINAL** |
| 1,e05,lit | 0.81 | 4.61 ± 0.41 | **8.58 ± 0.87** | 1.86× | 3.60 | 3.13 | 2.70 | 0.03 | CONFIRMED |
| 1,e05,vlbi | 0.51 | 9.55 ± 0.76 | **9.87 ± 0.91** | 1.03× | 2.43 | 2.27 | 1.73 | 0.03 | CONFIRMED |
| 1,e07,lit | 0.69 | 8.51 ± 0.72 | **11.65 ± 0.76** | 1.37× | 7.83 | 5.43 | 4.90 | 0.07 | CONFIRMED |
| 1,e07,vlbi | 0.44 | 10.64 ± 0.82 | **10.95 ± 1.01** | 1.03× | 6.13 | 5.77 | 5.13 | 0.07 | CONFIRMED |
| 1,eU,lit | 0.76 | 7.15 ± 0.62 | **11.35 ± 0.60** | 1.59× | 4.13 | 2.60 | 2.30 | 0.03 | CONFIRMED |
| 1,eU,vlbi | 0.53 | 10.08 ± 0.80 | **11.44 ± 1.69** | 1.14× | 2.80 | 2.13 | 1.43 | 0.07 | CONFIRMED |
| 2,e03,lit | 0.60 | 9.24 ± 0.76 | **10.05 ± 1.41** | 1.09× | 3.27 | 2.77 | 2.07 | 0.03 | CONFIRMED |
| 2,e03,vlbi | 0.36 | 12.13 ± 0.91 | **11.36 ± 1.02** | 0.94× | 1.43 | 1.77 | 1.43 | 0.03 | CONFIRMED |
| 2,e05,lit | 0.73 | 7.83 ± 0.67 | **11.67 ± 0.41** | 1.49× | 7.33 | 4.90 | 4.63 | 0.10 | CONFIRMED |
| 2,e05,vlbi | 0.57 | 10.73 ± 0.87 | **12.55 ± 1.82** | 1.17× | 5.30 | 3.97 | 3.20 | 0.03 | CONFIRMED |
| 2,e07,lit | 0.71 | 9.91 ± 0.84 | **13.73 ± 1.42** | 1.39× | 8.53 | 5.47 | 4.27 | 0.00 | CONFIRMED |
| 2,e07,vlbi | 0.64 | 11.88 ± 0.99 | **16.82 ± 3.17** | 1.42× | 6.73 | 4.10 | 2.73 | 0.03 | CONFIRMED |
| 2,eU,lit | 0.72 | 9.50 ± 0.81 | **13.60 ± 1.25** | 1.43× | 8.70 | 5.30 | 4.33 | 0.03 | CONFIRMED |
| 2,eU,vlbi | 0.53 | 11.67 ± 0.93 | **12.12 ± 0.76** | 1.04× | 7.40 | 7.03 | 6.33 | 0.00 | CONFIRMED |
| 3,e03,lit | 0.65 | 9.06 ± 0.76 | **12.16 ± 0.99** | 1.34× | 4.57 | 2.50 | 2.17 | 0.00 | CONFIRMED |
| 3,e03,vlbi | 0.33 | 13.52 ± 1.00 | **12.93 ± 0.87** | 0.96× | 1.97 | 2.20 | 1.93 | 0.00 | CONFIRMED |
| 3,e05,lit | 0.78 | 5.71 ± 0.50 | **10.22 ± 0.72** | 1.79× | 7.00 | 5.83 | 4.97 | 0.03 | CONFIRMED |
| 3,e05,vlbi | 0.66 | 7.68 ± 0.64 | **10.69 ± 0.97** | 1.39× | 6.40 | 4.50 | 3.97 | 0.07 | CONFIRMED |
| 3,e07,lit | 0.77 | 7.65 ± 0.66 | **11.85 ± 1.06** | 1.55× | 5.97 | 4.07 | 3.60 | 0.03 | CONFIRMED |
| 3,e07,vlbi | 0.70 | 7.72 ± 0.65 | **10.96 ± 1.70** | 1.42× | 6.67 | 5.07 | 3.93 | 0.03 | CONFIRMED |
| 3,eU,lit | 0.80 | 5.77 ± 0.51 | **11.33 ± 1.10** | 1.96× | 5.23 | 3.40 | 3.00 | 0.10 | CONFIRMED |
| 3,eU,vlbi | 0.65 | 7.58 ± 0.63 | **8.43 ± 1.27** | 1.11× | 5.07 | 4.77 | 4.13 | 0.10 | CONFIRMED |

---

## 3. THE ANCHOR BANK DEFECT — NARROWER THAN REPORTED, AND NOW CLOSED

FLOOR_FIX §5 reported that `anchor_ladder.npz` stores `offenders` rung-major while its metadata
is cell-major, and concluded: *"The array is the TRANSPOSE of its own metadata, and **nothing in
the file says so**."*

> ### **The first half is right. The second is not — the file DID say so.**
> `anchor_ladder.npz` already banked **`offender_index`**, a 48-element string column whose entry
> *j* labels `offenders[j]` (`'R0|-13.25|lit'`, `'R0|-13.25|vlbi'`, …). Re-cutting by that key
> reproduces the banked `emp_q95` **and** `zero_frac` to **exactly 0.000e+00** — and the
> permutation §5 reverse-engineered, `perm[j] = (j%6)·8 + (j//6)`, is **identically** the mapping
> `offender_index` already encodes (verified: `PERM == index map` → `True`).

So the bank was **self-describing**, not undeclared. §5's recommended fix — *"or add an explicit
index key"* — asked for a key the file already had.

**The trap is still real, but it is narrower than stated:** it is sprung only by a re-cut that
*ignores* `offender_index` and assumes `offenders[j]` pairs with metadata row *j*. That re-cut is
wrong by up to **79.9 nat**, and silently so. FLOOR_FIX's own script fell into it and then
correctly climbed out.

**FIXED — the array, not the metadata.** `offenders` is now permuted into **metadata-row order**,
so the naive `offenders[j]` ↔ row *j* re-cut is the **correct** one. `offender_index` is reordered
to match (it now agrees with its siblings instead of cross-cutting them) and
`offender_index_orientation` states the convention in words. The metadata order was left alone
because it is the order of the published §3 table and of every downstream reader.

**ANCHOR's published §3 table reproduces from the re-banked file, exactly:**

| rung | cell | zero % | EMP q95 (re-cut from re-banked offenders) | published §3 | dev |
|---|---|---|---|---|---|
| R0 | −13.25 lit | 93.3 | 2.395 | 2.395 | 0.0e+00 |
| R0 | −13.25 vlbi | 80.0 | 4.982 | 4.982 | 0.0e+00 |
| R0 | −13.00 lit | 57.3 | 10.416 | 10.416 | 0.0e+00 |
| R0 | −13.00 vlbi | 19.3 | 9.297 | 9.297 | 0.0e+00 |
| R0 | −12.75 lit | 5.3 | 30.046 | 30.046 | 0.0e+00 |
| R0 | −12.75 vlbi | 0.0 | 30.687 | 30.687 | 0.0e+00 |
| R0 | −12.50 lit | 0.0 | 78.796 | 78.796 | 0.0e+00 |
| R0 | −12.50 vlbi | 0.0 | 78.625 | 78.625 | 0.0e+00 |

All **48** rows (not just the 8 R0 controls shown) reproduce to **0.0e+00**. **ANCHOR's verdicts
stand; nothing is retracted.** Original preserved as `anchor_ladder.npz.preRECUT.bak`.

---

## 4. SURFACE'S RAW OFFENDERS — BANKED, ORIENTATION DECLARED

`reports/surface_nullN_offenders.npz`: **108 cells × 100 nullN offenders**, float64, ~132 kB.
Readback-gated — these vectors reproduce the banked `fN_zerofrac` and `fN_emp` to **0.000e+00**,
so they *are* the vectors the published floors were fitted to.

Learning from §3, the orientation is **declared, not implied**: an `index` column labels row *i*,
in the same order as every metadata column (`h`, `T`, `k`, `tier`, `struct`, `n`, `zero_frac`,
`emp_q95`), and an `orientation` string states in words that `off_i` ↔ `index[i]` ↔ meta row *i*,
with **no transpose and no implied permutation**. The convention's bootstrap error now exists for
SURFACE.

---

## 5. THE CANONICAL `b1_L_gwb` BANK — GENERATED AND GATED

The frozen GWB square root that removes `--cpus-per-task` from the seed.

| | |
|---|---|
| generated | `CW_transition/make_lgwb_bank.py`, **ACCRE**, BLAS threads pinned to **8** *before* numpy import |
| shape | **3248 × 3248** (116 pulsars × 28 GWB components) |
| fingerprint | **`71677a810cbc7187`** · cpus = 8 · 81 MB |
| Phi spectrum | cond 2.5e+06, **1624 near-degenerate adjacent eigenvalue pairs** — exactly half the spectrum. *This is the hazard's cause, made visible.* |
| reconstruction | max‖L Lᵀ − Φ‖ = 6.5e−27 (rel 4.0e−15) |

### THE GATE — **PASS**

`anchor.py gate` replayed **ANCHOR's 80 banked `ig_nullN_*` T=15 realisations THROUGH the bank**:

| column | max &#124;diff&#124; vs banked |
|---|---|
| dlnL_det · lnK · qmax · mapk · ptrue · on_true | **0.000e+00** (all six) |
| offender statistic | **0.000e+00** |
| **realisations bit-identical** | **80 / 80** |

Provenance confirmed as `BANKED b1_L_gwb.npz fp=71677a810cbc7187 cpus=8` — **zero fallback
warnings in the gate log**, so the replay genuinely ran through the bank rather than silently
recomputing. **The bank is CANONICAL: it reproduces the repo's existing banked realisations
bit-for-bit.**

### Where the bank lives — and why it is NOT in git

**Decision (Matt, this session): the 84.4 MB bank is NOT committed.** It sits on ACCRE scratch at
`CW_transition/b1_L_gwb.npz` — exactly where `trackB_b1_core.L_GWB_BANK` looks for it — and is
now **gitignored** so it cannot be added by accident.

What *is* committed is `reports/b1_L_gwb_manifest.npz` (36 kB): the fingerprint
(`71677a810cbc7187`), the file SHA-256, the shape, `cpus=8`, and the full g1 gate evidence. Any
machine can therefore **verify** it holds the canonical bank by comparing fingerprints.

**It cannot regenerate one.** The basis is not reproducible off-ACCRE — that is the entire point
of banking it. To rebuild on ACCRE:

```
python CW_transition/make_lgwb_bank.py     # BLAS pinned to 8 threads before numpy import
python hpc_harbor/anchor/anchor.py gate    # 80-realisation g1 replay: bit-identical or STOP
```

The generator **refuses to overwrite** an existing bank, and it STOPs if two eigh calls disagree
at the same thread count. **A regenerated bank is not canonical until the g1 gate passes.**

> **Consequence, stated:** a machine without this file (e.g. cronus) **cannot reproduce any banked
> noisy number in the repo.** It will hit the `RECOMPUTED-UNSAFE` warning path and draw a rotated
> GWB realisation. That is the accepted cost of keeping 84 MB out of git history.

---

## 6. WHAT THIS SETTLES, AND WHAT IT COSTS

**Settled:**

1. **N_onset = 59** (was bounded 57–67). Two onsets retracted, two ignited.
2. **Both at-risk onset cells RETRACTED** — the corrected counts are 0.77 and 0.90.
3. **`h*` unbounded below: REINSTATED at 7/18** — same count, **different columns**.
4. **The e = 0.3 switch-on is REFUTED.** Corrected threshold: **e = 0.5 for one eccentric member,
   e = 0.3 for two or more.** The e = 0.5 binding stands.
5. **ANCHOR's bank was self-describing** (§5's diagnosis was half wrong); the array is re-banked
   and the published table reproduces exactly.
6. **`b1_L_gwb` is canonical and gated** — `--cpus-per-task` is no longer an undeclared input.

7. **CHORUS's "trade inverts at n_ecc = 3" is no longer a clean claim** — it flips in 3 of 8
   (e, tier) combinations. Re-derive, do not quote (§2.3).

**The cost, stated plainly:** CHORUS's loudest headline — *a single e = 0.3 member suffices* — is
gone. It was an artifact of a Gumbel fitted to a 73 %-zero point mass, which understated the bar
by 53 %. The floor fix was worth doing, and this is what it bought.

**Two claims were checked rather than assumed, and both moved:** the `h*` column membership (§1.3)
and the trade inversion (§2.3). Anything else in SURFACE/CHORUS that is a *difference or ratio
between two cells* should be re-derived from the re-cut banks before it is quoted — the counts
moved, so the differences did too.

**Held for Matt:** the corrected onset-surface table and figures, the STORY supersessions that
depend on N_onset and on the e-threshold, CHORUS §5's now-false "every eccentric mix" sentence,
CHORUS's trade-curve section, every tag flip, and `criterion-v2.2`.
