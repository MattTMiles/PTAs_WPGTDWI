# CHORUS — the MIXED-ECCENTRICITY POPULATION (does one eccentric source clock its circular neighbours?)

> ## ⚠ AMENDED 2026-07-13 BY THE FLOOR FIX (RECUT). READ THIS FIRST.
>
> **All 26 of 26 CHORUS cells sit above the floor-validity gate** (nullN zero-fraction 0.33–0.81),
> so **every floor in this campaign is restated** as the empirical q95 with a bootstrap error
> (ANCHOR §4). **23 of 26 floors rise.** Counts re-cut from the per-realisation `ch_sig` banks;
> gates A and B both **0.000e+00**. Primary source `reports/RECUT_floors.md` §2; banks
> `reports/recut_chorus.npz`, `ch_analysis_recut.npz`, `ch_floors_recut.npz`.
>
> ### THE HEADLINE THAT DID NOT SURVIVE
>
> > **"A single e = 0.3 member switches the count on" is REFUTED.**
> > **lit collapses 1.57 → 0.70 — below the bar. vlbi survives at 1.03, but MARGINAL: it clears the
> > bar at the floor and fails at floor + bootstrap error (0.60).**
> >
> > **The corrected, externally quotable statement:
> > ONE eccentric member → the switch-on is at e = 0.5.
> > TWO OR MORE → it is at e = 0.3 (CONFIRMED, both tiers).**
> > The threshold is **not a property of eccentricity alone** — it depends on how many members carry
> > it. The docs' interim binding to e = 0.5 STANDS, and it was the right call.
>
> **What survives, structurally:** the clock is **NOT shared** (0 of ~120 lifted certifications
> circular-attributed — independent of the floor); **eccentric structure transforms the count**;
> and the **~10× lever** holds at **14.8× (lit) / 12.4× (vlbi)** on the re-cut counts (was
> 11.2× / 14.2×). **What does NOT survive:** the "every eccentric mix clears the bar" sentence
> (§1 reading 1, §3 deliverable 1), and the **"the trade inverts at n_ecc = 3"** claim, which is
> **demoted to not-clean** — it flips status in 3 of 8 (e, tier) combinations (§1 reading 3, §3
> deliverable 3).
>
> **Figures `ch_surface.png`, `ch_trade.png`, `ch_pairs.png` are STALE** — regenerate from
> `recut_chorus.npz` before use. **`reports/F4_the_switch.png` and `F9_the_growth_path.png`, and
> any external text quoting the single-member e = 0.3 switch-on, are FLAGGED for correction.**

**Agent:** CHORUS · ACCRE · tag `criterion-v2.1` (`git rev-parse HEAD` → `6bec3d6` ✓;
discovery `136b270f`, ee `f73b8e0`) · **Date:** 2026-07-13 · **PURE EXECUTION** (no commits,
no tracked-file edits). **Amended 2026-07-13 (RECUT + cronus doc session); adopted at
`criterion-v2.2`.**

**Scratch paths (host):** code `hpc_harbor/chorus/` (`chorus.py`, `chorus_analysis.py`,
`ch_job.sbatch`, `ch_gate.sbatch`, `ch_README.md`), results `CHORUS_results/` (lean keyed npz;
raw statistics always banked: per-cell dlnL/lnK/qmax/mapk — **mapk is PRIOR-REFERENCED (k=0 at
the EM-prior mean), the IMPLEMENTABLE frame; `n_true_grid` is carried beside it and every
truth-referenced column is labelled ORACLE** — the D4 oracle-indexing lesson applied at bank
time), logs `hpc_harbor/logs/ch_*`, `chorus_*`. Lane
`-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1`,
**--cpus-per-task=8 (load-bearing — see §0.1)**.

---

## THE ANSWER IN ONE PARAGRAPH

**Eccentricity in a population transforms the certified count — but the clock is NOT
shared.** The all-circular census population sits below onset (**0.37 / 0.47** correct
certifications per realisation, lit / vlbi, under its own re-cut floors), exactly as
IGNITE-2 and SURFACE found. Replace ONE of the three loud members with an eccentric source
and **the count switches on at e = 0.5** (3.13 / 2.27) — **not at e = 0.3, which is REFUTED in lit
(0.70) and MARGINAL in vlbi (1.03)** — and reaches **5.4 / 5.8** per realisation at
e = 0.7: **14.8× (lit) / 12.4× (vlbi)** over the circular baseline, at fixed census loudness, while
the mix's own floor rises by up to 1.96× and its trials term grows ~8×. **With two or more
eccentric members the switch-on falls back to e = 0.3, confirmed in both tiers** — the threshold is
a property of the *mixture*, not of *e*. The surface peaks at
n_ecc = 2 in vlbi (7.03) but at n_ecc = 3 in lit (5.83), and **the published "the trade inverts by
n_ecc = 3 at high e" is DEMOTED to not-clean** — it flips status in 3 of the 8 (e, tier)
combinations under the corrected floors and must be re-derived, not quoted (§1 reading 3).
But the campaign's namesake question resolves NEGATIVE: in 20 exact seed-paired
(e = 0 vs e = 0.7) realisations, **every one of the ~120 certifications the eccentric arm
gains is attributable to the eccentric member's own comb template — zero are lifted
circular-member registrations** (pre-registered verdict: **MARGINAL** lit 1/10, **ABSENT**
vlbi 0/10) — and the soft loop adds nothing (all 30/30 trajectories flat: the
eccentric-seeded loop holds its 18-cert seeds exactly as the circular-seeded loop holds its
0–1). Certification in a mixture is a property of the SOURCE, not a shared array resource;
the single-source no-gos are rescoped anyway, because one moderately eccentric loud member
certifies a dozen pulsars single-handedly where the whole circular population certifies
none. The pre-registered scrambled-loop STOP fires (2/10 certify at some iteration, 1
keeps: a Δk_oracle = −266 comb noise-lock, seed-static — D1's hole, mixed flavour; the loop
self-cleans the other) — so the verdict is **STOP**, with the IGNITE-2 anatomy intact:
every failure is inherited from the criterion's seeds, none is generated by the loop. Two
findings travel beyond the campaign: the g1 gate caught a repo-wide reproducibility hazard
(**banked noisy realisations reproduce only at cpus-per-task = 8** — the NoiseDrawer's GWB
square root rotates with the BLAS thread count, §0.1), and the comb-null floor family
INVERTS the fN/fALL ordering (§1).

---

## 0. STAGE 0 — GATES (all three PASS; pre-registered, stop-at-first-failure)

**The machinery** (implementation completions C1–C9, fixed before any scoring, recorded in
`chorus.py`'s docstring): the 16-member census population (3 loud + 13 faint, pop seed 3000,
census loudness −13.25/−14.25) keeps its banked slot order; an eccentric member's own slot
carries its **n = 2 harmonic** and its other 31 harmonics (n = 1, 3..32, WEAVE-3.3 chirp-tied,
ATLAS's map + merger guard verbatim, N_HARM = 32 per ATLAS's validity notes) are appended
slots — graph shapes 16/47/78/109 for n_ecc = 0–3. A harmonic is ACTIVE iff
g(n,e) ≥ 1e-3·max g; only active slots enter the fringe grid (dL = min mode spacing), the
K_counted trials term, and the band-occupancy accounting (capacity 0.4·N·T·ΔF, P2-C band
convention). At e = 0 the active set is {n = 2} and the collapse to the banked machinery is
structural. Floors are refit per mix (spec C3 note), Gumbel-MLE α = 0.05, never max-of-N.

| gate | pre-registered condition | result |
|---|---|---|
| **g1** | e→0 limit reproduces a banked GEO draw's certified set bit-comparably | **PASS — bit-identical**, and stronger than required: all six raw columns (dlnL/lnK/qmax/on_true/mapk/n_true_grid) of banked `ig_sig_h1325_T30_lit_g0_n1200100` reproduce with max\|diff\| = 0.0 through the n_ecc=0 path; certified/detected sets identical at both the v1 (9.01) and fresh-lit (38.86) floors; member draw ncw-invariant (theta16 == theta47[:16 members], 0.0); e=0 grid collapse exact (dL, K_counted bit-equal, active = {n=2}); 47-slot ↔ 16-slot drift 3.7e-9 nat with identical certified set + mapk |
| **g2** | one single-eccentric-source cell reproduces ATLAS's κ at (e=0.65, Mc=1e9, f_orb=1e-8.5) to <5% | **PASS — exact**: κ_measured = **127.4660** vs ATLAS banked 127.4660 (rel 0.0000); both raw Fisher diagonals match ATLAS's banked `Fmc_cond` to all printed digits (7.051134e+01 / 4.339807e-03). ATLAS's setup replicated exactly (fiducial span, TMAX 7e8 s guard, mask=0 Asimov, loud0 geometry, other members inert) THROUGH the 47-slot population embedding |
| **g3** | floors refit per mix: counterpart-matched nulls, N≥100/cell, Gumbel α=0.05, fALL beside | **PASS**: (a) CHORUS's offender+Gumbel code reproduces IGNITE-2's banked fresh floors from the raw 540-null banks (fN 38.86/7.59, fALL 65.47/11.68, to the banked digits); (b) mixed-cell pilot null (m1e07_vlbi, 47 slots) runs the full path finite with zero false detections, occupancy accounting live (47 active slots / capacity 925) |

Jobs: gates 12507980 (1h18m, ALL PASS); diagnostics 12507594/12507755.

### 0.1 WHAT THE g1 GATE CAUGHT FIRST — a repo-wide reproducibility hazard (report to cronus)

The first gate attempt (job 12507521) FAILED g1a with a distinctive anatomy: dlnL shifted
O(0.1–1) nat on essentially every pulsar, qmax by up to 0.22, mapk flips only among
noise-degenerate fringes — while lnK, n_true_grid, and the certified/detected sets were
untouched. Diagnosis (jobs 12507594/12507755, then confirmed by configuration):
**`trackB_b1_core.NoiseDrawer` builds the GWB noise square root `L_gwb = V·sqrt(w)` by CPU
`np.linalg.eigh` of the 5336×5336 HD-correlated Φ_gwb, whose spectrum is near-degenerate;
the LAPACK eigenvector basis inside degenerate subspaces depends on the BLAS thread count.**
A job with `--cpus-per-task=16` therefore draws a *rotated* — different but
equal-distribution — GWB realisation at the same seed. Banked noisy realisations reproduce
bit-identically **only at the banked convention `--cpus-per-task=8`** (verified: 8 CPUs →
bit-identical on 3 cards across 2 nodes; 16 CPUs → reproducibly shifted, full anatomy in the
job-12507521/12507900 logs). Every CHORUS job is pinned to 8. **Durable fix is cronus's
call:** a deterministic single-threaded factorisation in NoiseDrawer, or banking `L_gwb`
itself. Until then, *cpus-per-task is part of the seed.*

### 0.2 Pre-registered cell constants (C8, drawn once, before compute)

e-distributions: e03/e05/e07 → all eccentric members at 0.3/0.5/0.7; eU → per-member
U(0.2, 0.8) with `default_rng(42000 + n_ecc)`: n_ecc=1 → [0.5801]; n_ecc=2 →
[0.5631, 0.6753]; n_ecc=3 → [0.3522, 0.2957, 0.6792]. Active-harmonic counts at census
parameters: e=0.3 → 8, e=0.5 → 17, e=0.7 → 32, e=0.8 → 32 (magnitude-truncated per ATLAS's
validity note; κ contour unaffected). Merger-guard clips at census parameters: 0.

Verdict language, fixed before any Stage-2 scoring: **clock-sharing REAL** if the paired
circular-attributed lift exceeds the floor-refit noise (quadrature of the two cells' fN fit
errors) for ≥ half the pairs; **MARGINAL** if for ≥ 1 but < half; **ABSENT** otherwise.
Attribution rule (fixed with it): a pulsar is ecc-attributed when its evidence needs
member 0's template (dlnL − dlnL_ct > 1 nat, where dlnL_ct is the member0-inert-template
rescore banked beside every paired realisation); circular-attributed otherwise.

---

## 1. STAGE 1 — THE MIX GRID (certified count vs mixture; 26 cells × 30 realisations)

Per-mix floors first (the mixture changes the null — spec C3 note, enforced): 100 nullN +
30 nullA + 30 nullL per cell (m0_vlbi reuses IGNITE-2's banked 150/270-null cell — identical
model, provenance in `ch_floors.npz`), **α = 0.05, nullN zero-fraction quoted — and it is the
zero-fraction that governs, because every cell in this campaign fails the validity gate.**

> ### **RE-CUT: ALL 26 OF 26 CELLS HAVE ZERO-FRACTION > 20 % (0.33–0.81).**
> **Every CHORUS floor is the empirical q95 with a bootstrap error. Not one Gumbel survives.
> 23 of 26 floors RISE, by up to 1.96×.** The comb's null is dominated by its point mass at zero —
> a Gumbel fitted to it understates the bar and errs **permissive**, and CHORUS's cells are the
> worst case in the whole repo. *This is the campaign the floor fix cost the most.*

**Two floor findings before any signal is scored** (both re-derived on the adopted floors):
1. **The comb raises its own bar.** The operative floor grows with the mix (lit tier: m0
   **7.00 ± 1.10** → m1e07 **11.65 ± 0.76** → m2e07 **13.73 ± 1.42** nat; vlbi: 7.06 → 10.95 →
   16.82) — the confusion cost of harmonic occupancy, measured in the null family itself. **The
   finding survives the re-cut; its slope is gentler** (the fix lifts the *circular* baseline
   floor hardest, because m0's null is the most silent of all). (Exception, and it also survives:
   the n_ecc = 3 high-e cells' floors *fall* again (m3eU vlbi 8.43) — at 109 active slots the
   per-pulsar fringe grid is so fine that K_counted explodes instead: the trials term takes
   over from the floor, K_sum ~11M vs m0's ~1M.)
2. **fALL inverts under the comb.** In every eccentric lit cell the blind-robust fALL sits
   BELOW fN (e.g. m1e07_lit 6.20 vs 8.51) — scrambled-comb templates noise-lock *less* than
   true-position ones, the opposite of the circular cells (m0_lit 8.68 vs 4.31). The
   fALL column is carried beside every count regardless (D1 convention). **NOT re-cut** — the
   fALL zero-fractions were never banked, so this finding stands on the pre-fix estimator. It is
   a *sign* statement, and the fix would have to be enormous to flip a sign; but it is not re-cut,
   and it is not quoted as if it were.

**THE COUNT SURFACE, RE-CUT** (correct = certified AND on the true fringe — ORACLE-verified;
mapk banked prior-referenced). **Bold** = the count at the adopted floor; `[·]` = the count at
floor + its bootstrap error, **which is the number the grade is issued on**. Grades:
**CONFIRMED** = clears the >1 bar at floor + error; **MARGINAL** = clears it only at the floor;
**below** = fails at the floor. Published counts in the `old` column. Full table
`recut_chorus.npz` / `ch_analysis_recut.npz`; figure `ch_surface.png` **is stale**.

| mix | z | floor lit (adopted) | **lit** [+bs] | old | grade | floor vlbi | **vlbi** [+bs] | old | grade | occ/925 |
|---|---|---|---|---|---|---|---|---|---|---|
| m0 (all circular) | .73/.45 | 7.00 ± 1.10 | **0.37** [0.30] | 0.70 | below | 7.06 ± 0.40 | **0.47** [0.43] | 0.43 | below | 16 |
| **m1 e=0.3** | .73/.48 | **11.30 ± 1.02** | **0.70** [0.43] | 1.57 | **below — REFUTED** | **10.78 ± 1.54** | **1.03** [0.60] | 1.13 | **MARGINAL** | 23 |
| **m1 e=0.5** | .81/.51 | 8.58 ± 0.87 | **3.13** [2.70] | 3.60 | **CONFIRMED** | 9.87 ± 0.91 | **2.27** [1.73] | 2.43 | **CONFIRMED** | 32 |
| m1 e=0.7 | .69/.44 | 11.65 ± 0.76 | **5.43** [4.90] | 7.83 | CONFIRMED | 10.95 ± 1.01 | **5.77** [5.13] | 6.13 | CONFIRMED | 47 |
| m1 eU(0.58) | .76/.53 | 11.35 ± 0.60 | **2.60** [2.30] | 4.13 | CONFIRMED | 11.44 ± 1.69 | **2.13** [1.43] | 2.80 | CONFIRMED | 38 |
| m2 e=0.3 | .60/.36 | 10.05 ± 1.41 | **2.77** [2.07] | 3.27 | CONFIRMED | 11.36 ± 1.02 | **1.77** [1.43] | 1.43 | CONFIRMED | 30 |
| m2 e=0.5 | .73/.57 | 11.67 ± 0.41 | **4.90** [4.63] | 7.33 | CONFIRMED | 12.55 ± 1.82 | **3.97** [3.20] | 5.30 | CONFIRMED | 48 |
| m2 e=0.7 | .71/.64 | 13.73 ± 1.42 | **5.47** [4.27] | 8.53 | CONFIRMED | 16.82 ± 3.17 | **4.10** [2.73] | 6.73 | CONFIRMED | 78 |
| m2 eU(0.56,0.68) | .72/.53 | 13.60 ± 1.25 | **5.30** [4.33] | 8.70 | CONFIRMED | 12.12 ± 0.76 | **7.03** [6.33] | 7.40 | CONFIRMED | 67 |
| m3 e=0.3 | .65/.33 | 12.16 ± 0.99 | **2.50** [2.17] | 4.57 | CONFIRMED | 12.93 ± 0.87 | **2.20** [1.93] | 1.97 | CONFIRMED | 37 |
| m3 e=0.5 | .78/.66 | 10.22 ± 0.72 | **5.83** [4.97] | 7.00 | CONFIRMED | 10.69 ± 0.97 | **4.50** [3.97] | 6.40 | CONFIRMED | 64 |
| m3 e=0.7 | .77/.70 | 11.85 ± 1.06 | **4.07** [3.60] | 5.97 | CONFIRMED | 10.96 ± 1.70 | **5.07** [3.93] | 6.67 | CONFIRMED | **109** |
| m3 eU(0.35,0.30,0.68) | .80/.65 | 11.33 ± 1.10 | **3.40** [3.00] | 5.23 | CONFIRMED | 8.43 ± 1.27 | **4.77** [4.13] | 5.07 | CONFIRMED | 63 |

**Readings, re-cut.**

**(1) THE SWITCH-ON THRESHOLD IN e — and it is NOT what was published.** The all-circular census
population sits below onset (0.37 / 0.47). The published sentence *"**Every** eccentric mix clears
the >1 bar — a single e = 0.3 member suffices at either tier"* is **FALSE and is struck.**
[SUPERSEDED → **every eccentric mix clears the bar EXCEPT a single e = 0.3 member, which is REFUTED
in the lit tier (0.70) and MARGINAL — not confirmed — in vlbi (1.03, dying at floor + bootstrap
error).**]

> ### THE CORRECTED, EXTERNALLY QUOTABLE STATEMENT
> ### **With ONE eccentric member the switch-on is at e = 0.5. With TWO OR MORE it is at e = 0.3 (CONFIRMED, both tiers).**
> **The threshold is not a property of eccentricity alone; it depends on how many members carry
> it.** The interim binding to e = 0.5 that the provisional floor fix imposed **STANDS — and it was
> the right call.**

| n_ecc | tier | e = 0.3 | e = 0.5 | e = 0.7 | **switch-on** |
|---|---|---|---|---|---|
| **1** | lit | 0.70 [0.43] **below** | 3.13 [2.70] CONFIRMED | 5.43 [4.90] CONFIRMED | **e = 0.5** |
| **1** | vlbi | 1.03 [0.60] MARGINAL | 2.27 [1.73] CONFIRMED | 5.77 [5.13] CONFIRMED | **e = 0.5** (0.3 not confirmed) |
| 2 | lit | 2.77 [2.07] CONFIRMED | 4.90 [4.63] CONFIRMED | 5.47 [4.27] CONFIRMED | **e = 0.3** |
| 2 | vlbi | 1.77 [1.43] CONFIRMED | 3.97 [3.20] CONFIRMED | 4.10 [2.73] CONFIRMED | **e = 0.3** |
| 3 | lit | 2.50 [2.17] CONFIRMED | 5.83 [4.97] CONFIRMED | 4.07 [3.60] CONFIRMED | **e = 0.3** |
| 3 | vlbi | 2.20 [1.93] CONFIRMED | 4.50 [3.97] CONFIRMED | 5.07 [3.93] CONFIRMED | **e = 0.3** |

**The cost, stated plainly:** the loudest headline this campaign produced — *a single e = 0.3 member
suffices* — was an artifact of a Gumbel fitted to a 73 %-zero point mass, which understated the bar
by **53 %** (7.39 → 11.30 nat, a 6.2σ move against its own quoted fit error). **The floor fix was
worth doing, and this is what it bought.**

**(2) At n_ecc = 1 the count is monotone in e** (lit: 0.70 → 3.13 → 5.43) — **14.8× (lit) /
12.4× (vlbi)** over the circular baseline at e = 0.7, *while the mix's own floor rises by up to
1.96×*: **the certification gain beats the confusion cost across the whole measured box.** This is
the campaign's surviving quantitative headline, and it is **re-derived, not inherited** (published:
11.2× / 14.2× — the ratio moved because *both* counts moved, which is exactly why RECUT §6 forbids
quoting a ratio without re-deriving it).

**(3) THE TRADE INVERSION IS DEMOTED — it is no longer a clean claim.**
[SUPERSEDED → *"the trade inverts between n_ecc = 2 and 3 at high e"* was true on the pre-fix
floors. Under the corrected floors **it flips status in 3 of the 8 (e, tier) combinations**, and
**not in the same direction**:]

| e, tier | published m1→m2→m3 | inverts? | **re-cut m1→m2→m3** | **inverts?** |
|---|---|---|---|---|
| e=0.3 lit | 1.57 · 3.27 · 4.57 | no | 0.70 · 2.77 · **2.50** | **YES** ← changed |
| e=0.5 lit | 3.60 · 7.33 · 7.00 | YES | 3.13 · 4.90 · **5.83** | **no** ← changed |
| e=0.7 vlbi | 6.13 · 6.73 · 6.67 | YES | 5.77 · 4.10 · **5.07** | **no** ← changed |
| e=0.7 lit · eU lit · eU vlbi | — | YES | — | YES (unchanged) |
| e=0.3 vlbi · e=0.5 vlbi | — | no | — | no (unchanged) |

Published: 5 of 8 inverted. Re-cut: **4 of 8 — but not the same four.** The inversion is a
*difference between two counts that both moved*, and it is not robust to the floor correction.
**The surface no longer even peaks in the same place: vlbi peaks at n_ecc = 2 (m2eU, 7.03) but lit
now peaks at n_ecc = 3 (m3e05, 5.83).** **DO NOT QUOTE "the trade inverts at n_ecc = 3" without
re-deriving it from `recut_chorus.npz`. It is not refuted; it is no longer clean.** The *mechanism*
(occupancy 109/925 and ~11× K_counted growth at m3e07) is untouched and is still the right place to
look — but the claim that it *wins* at n_ecc = 3 is now a coin-flip across the (e, tier) grid.

**(4) Wrong-cert rates stay ≤ 0.23/real everywhere** (cf. 23/50 realisations at IGNITE's
above-onset cells) — **the counts are not purity-bought**, and the re-cut moves them *down*, not
up (the wrong column falls with the count, as it must). fALL-scored counts track fN counts to
~10 % — **pre-fix estimator, not re-cut** (finding 2 above).

## 2. STAGE 2 — THE CLOCK-SHARING TEST (the campaign's reason to exist)

**(a) Paired conditional test.** 10 pairs/tier of (m0 vs m1e07) realisations sharing
noise seed AND physical truth distances (drawn on the circular-arm grid in both arms — the
pairing is exact by construction, C6), each banking the full raw columns plus the
member0-inert-template rescore `dlnL_ct` (recovery with the comb removed, data unchanged).

The certified count jumps in every pair (lit: certC 0–3 → certE 1–15; vlbi: 0–1 → 1–15).
**But the lift is the comb's own, not shared:** under the pre-registered attribution rule
(ecc-attributed iff dlnL − dlnL_ct > 1 nat), **zero of the ~120 lifted certifications
across all 20 pairs are circular-attributed.** Every pulsar the eccentric arm certifies is
certified through member 0's own comb template — its harmonic pulsar terms — not through a
lifted circular-member registration. The largest circular-attributed paired margin gains
sit at the floor-noise level: lit max +2.84 nat with 1/10 pairs above the 0.81-nat
floor-refit noise; vlbi max +0.21 with 0/10 above 0.95 nat.

> **Pre-registered verdict (conditional level): lit MARGINAL (1/10), vlbi ABSENT (0/10).**
> Certification in a mixed population is, at the fixed-source (E-step) level, a property of
> the SOURCE — the eccentric member certifies its own pulsar terms and drags the array's
> count up single-handedly — not a shared array resource.

**(b) The soft loop (spec §3 reference implementation), eccentric- vs circular-seeded.**
5 pairs/tier through the loop (B1Marg objective, nothing hard-locked, damped Newton +
backtracking, E-step re-score under the cell's fresh fN floor; eccentric members fit in the
PHYSICAL frame — 7 dims including e — through the tie; every per-iteration raw column
banked). The joint-fit channel is where SIREN's lag-diversity mechanism would operate (the
conditional test holds sources at truth and cannot see it). **It does not operate either:**
all 20 signal trajectories are FLAT — circular-seeded loops hold their sparse seeds (lit
0/0/0/0/1, vlbi 0/1/1/0/1 certs, all true), eccentric-seeded loops hold their large ones
(lit 18/6/2/2/6, vlbi 13/8/2/0/7) — no growth, no loss, no wrong-cert change, 2 iterations
each with every proposed step rejected (the soft M-step at truth has ≈ zero gradient, the
IGNITE-2 dynamics reproduced on the mixed problem). **The eccentric-seeded loop does not
consolidate further than the circular-seeded one; its entire advantage is in its seeds.**
One vlbi eccentric-arm realisation carries 3 seed-static wrong certs beside 4 true
(J1017-7156/J1101-6424/J1946-5403, |Δk_oracle| = 99–157, qmax = 1.0, present from iteration
0, untouched by the loop) — IGNITE's §10.8.3 noise-lock class alive under the comb model.

**(c) Scrambled-source null through the mixed pipeline (5/tier through the loop).**
Pre-registered condition: *nothing detects, or STOP + anatomy.* **STOP fires: 2 of 10
scrambled realisations certify at some iteration; 1 keeps its certification to the fixed
point.** The anatomy, in full:
- **The keeper** (lit, g0): J1640+2224, dlnL = 12.14 vs floor 8.51 ± 0.72, qmax = 1.000,
  mapk = −111 (prior-referenced), Δk_oracle = −266 — a confident noise-lock under a
  scrambled comb, present at iteration 0 and never touched by the loop. The D1
  wrong-counterpart hole, mixed-model flavour; seed-inherited, not loop-generated.
- **The self-clean** (vlbi, g3): J1658-5324 certifies at iteration 0 (dlnL = 21.3!,
  Δk_oracle = −269) and the loop's M-step + re-score DROPS it — IGNITE-2's self-cleaning
  behaviour reproduced in the mixed pipeline.
No scrambled trajectory grows; wrong-cert counts are flat in all 30/30 loop trajectories.
**Every STOP event is inherited from the criterion's seeds; none is generated by the loop**
— the IGNITE-2 anatomy, unchanged by the mixture.

> **PRE-REGISTERED CLOCK-SHARING VERDICT: MARGINAL (lit, 1/10 pairs above floor-refit
> noise) / ABSENT (vlbi, 0/10) — and the loop channel adds nothing.** The eccentric
> source's self-clock does NOT propagate to its circular neighbours at either the
> conditional or the joint-fit level. Certification in a mixed population is a property of
> the SOURCE. **And the count surface transforms anyway** — because one eccentric member's
> own comb certifies up to ~18 pulsars single-handedly. The single-source no-gos are
> rescoped by the mixture, but through the comb's own pulsar terms, not through a shared
> clock.

## 3. STAGE 3 — DELIVERABLES

**(1) The certified-count-vs-mix surface** — §1 table, **re-cut**; every number beside its floor's
zero-fraction, adopted floor and bootstrap error; `recut_chorus.npz` / `ch_analysis_recut.npz`;
figure `ch_surface.png` **STALE**. Headline: the circular census population is below onset
(**0.37 / 0.47**); **the switch-on is at e = 0.5 for ONE eccentric member and at e = 0.3 for TWO OR
MORE**; the surface peaks at n_ecc = 2 in vlbi (7.03) and at n_ecc = 3 in lit (5.83).
[SUPERSEDED → *"every eccentric mix clears onset; the surface peaks at n_ecc = 2 (8.7 lit / 7.4
vlbi) and inverts by n_ecc = 3 at high e"* — **the first clause is FALSE** (single-member e = 0.3
is REFUTED in lit, MARGINAL in vlbi) and **the third is DEMOTED to not-clean** (§1 reading 3).]

**(2) The clock-sharing verdict** — §2: **MARGINAL (lit) / ABSENT (vlbi)**, paired deltas
and the member0-inert attribution banked per pair (`ch_analysis.npz`: `pair_delta`,
`pair_delta_ct`, `pair_ecc_attr`), figure `ch_pairs.png`. Zero circular-attributed lifted
certifications in ~120 lifted certs across 20 exact pairs.

**(3) The capacity-vs-clock trade curve** — `ch_trade.png` **STALE**; raw material (occupancy,
capacity, K_sum, active counts) banked per cell and **unaffected** by the floor fix (occupancy and
K_sum are floor-independent).

[SUPERSEDED → *"the count rises steeply to ~50–70 active slots (5–8 % of the 0.4·N·T·ΔF capacity),
saturates, and **inverts by 109 slots (12 %)**; m3e07 posts 6.0/6.7 against m2e07's 8.5/6.7 and
m1e07's 7.8/6.1"* — **DEMOTED TO NOT-CLEAN.** Re-cut, m3e07 posts **4.07 / 5.07** against m2e07's
**5.47 / 4.10** and m1e07's **5.43 / 5.77**: the inversion **holds in lit and REVERSES in vlbi**,
and across the full (e, tier) grid it flips status in 3 of 8 combinations (§1 reading 3). **The
trade curve must be re-derived from `recut_chorus.npz` before any version of it is drawn or
quoted.** This is a difference between two counts that both moved — the class of claim RECUT §6
names explicitly.]

**What still stands here, and it is the mechanism rather than the crossing:** the binding cost at
high occupancy is the **K_counted trials term** (K_sum grows ~11× from m0 to m3e07) and the finer
joint fringe grid, **not the floor** — which *falls* at the n_ecc = 3 high-e cells while the count
falls with it. That statement is floor-independent and survives intact. **What is no longer
supportable is the claim that the trade CROSSES at n_ecc = 3, or that the crossing sits at
~8–12 % occupancy.** The capacity ceiling is real; where it bites is now an open number.

**(4) The updated join panel.**
- **vs the ATLAS corner:** ATLAS's single-source requirement for self-clocking was
  (e ≳ 0.6, Mc ≳ 1e9, f_orb ≳ 1e-8). In a POPULATION at census loudness — a full grid step
  below SURFACE's census-structure onset h\* = −12.50 — **one e = 0.5 loud member turns the
  certified count on (3.13 / 2.27 per realisation), and one e = 0.7 member posts 5.4 / 5.8**.
  The certification road's requirement relaxes from "louder or longer than anything
  in the modelled box" (IGNITE-2 §4) to **"the population contains one MODERATELY-TO-HIGHLY
  eccentric loud member (e ≳ 0.5), or two or more mildly eccentric ones (e ≳ 0.3)"**. The mixed
  case does not move ATLAS's absolute-σ_mc corner (that is a
  single-source Fisher statement, untouched here); it moves the COUNT criterion's onset.
  [SUPERSEDED → *"one e = 0.3 loud member already turns the certified count on (1.6/1.1)"* —
  **REFUTED by the re-cut**; the single-member threshold is **e = 0.5**. Note this pulls the
  requirement **back toward ATLAS's own corner** (e ≳ 0.6): the two thresholds, measured on
  entirely different statistics, now nearly agree — which is a *harder* requirement than the
  published CHORUS claimed, and a more coherent one.]
- **vs the SCOUT clock:** SCOUT priced certification as a loudness lottery
  (N̄ ≲ 0.01–0.1). CHORUS converts it to an ECCENTRICITY lottery at fixed census loudness:
  the question is now the fraction of census-class loud sources with **e ≳ 0.5 (single) or the
  probability of TWO OR MORE at e ≳ 0.3** — an astrophysical prior (Taylor/Farr handoff class),
  not a rarity of loudness draws. **The corrected threshold makes this lottery harder, not
  easier**, and *the two-or-more branch is now a genuinely different question from the
  one-member branch* — a population-multiplicity prior, not just an eccentricity prior.
- **vs SURFACE's structure lever:** SURFACE measured 3+13 → 5+11 (loudness structure) as a
  **2.5× median / 6.1× best** count lever at fixed (h, T, tier) (re-cut). CHORUS measures
  1-of-16-eccentric (waveform structure) as a **14.8× (lit) / 12.4× (vlbi)** lever at the same
  fixed loudness — **eccentric structure is the strongest single lever yet measured in the box**,
  and the re-cut *strengthens* the comparison rather than weakening it — with the caveat that the
  two levers were measured on different floors (each mix refits its own null, as it must).
- **The D1 hole travels, unchanged:** the scrambled keeper (|Δk_oracle| = 266) and the
  3-wrong-cert seed realisation (|Δk| = 99–157) are the same wrong-counterpart class as
  ever; the mixture neither opens nor closes it. Wrong-cert rates beside every count
  (≤ 0.23/real everywhere measured).

## 3.1 Caveats that travel
- **30 realisations/cell (5 skies × 6 weathers):** GEO's sky-dominance shows directly — the
  per-realisation cert count at m1e07 spans 2–18 across skies; per-cell rates carry the
  ±0.2-class-or-larger sky-sampling error. The dynamical loop statements are per-trajectory.
- **The loops start at the true source** (conditional convention, as IGNITE-2's did);
  mis-registered starts remain untested (§10.12 item 3), now also for the mixed model.
- **N_HARM = 32 truncation:** exact through e = 0.7 (ATLAS validity); the eU draws top out
  at e = 0.68. No merger-guard clips anywhere (n_clip = 0 banked per cell).
- **Floors at N = 100** (fresh cells; 150/270 for the reused m0_vlbi bank), nullN zero-fractions
  **0.33–0.81** quoted beside every count. **The Gumbel point-mass-at-zero caveat (IGNITE-2 §2)
  did not merely "apply at the high-zero-fraction lit cells" — it applies at ALL 26 CELLS, and it
  cost this campaign its headline.** Every floor here is now the empirical q95 with a bootstrap
  error. *A caveat that is true of every cell in the campaign is not a caveat; it is the result,
  and CHORUS filed it as the former.*
- **The fALL inversion under the comb** (§1 finding 2) is measured, not explained — a
  candidate anatomy (scrambled combs spread power off every fringe) is recorded as a
  question, not a claim.
- **Frames:** every banked mapk is prior-referenced (implementable); "correct" labels used
  n_true_grid and are ORACLE-labelled wherever they appear. No co-registration statistic
  was computed here (none is adopted in criterion-v2.1).
- **cpus-per-task = 8 is part of the seed** (§0.1) — any re-scoring of these banks from
  regenerated data must pin it.

## 4. EXECUTION RECORD

| item | value |
|---|---|
| lane | `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1`, cpus 8 |
| gates | **12507980** ALL PASS (g1 890 s, g2 3298 s, g3 485 s); first attempt 12507521 failed g1a → thread-count hazard found (§0.1); diagnostics 12507594 / 12507755 (bit-identical at 8 cpus, 3 cards / 2 nodes) |
| warm (4 shapes) | **12507981**: builds 214/190/380/436 s; first-in-process realisation 150/225/360/409 s; **warm realisation 1–2 s at every shape** |
| resume drill | **12508960_0** scancel'd at 7 m 32 s (3 npz banked) → resubmission **12509045** logged `already banked: 1; to run: 249` and completed (production kill→resume proven) |
| nulls array | **12509045 + 12509081** (16 shards, %8): 4 000/4 000 banked, ~40 min/shard |
| floors (CPU) | login node → `ch_floors.npz` (26 cells, N ≥ 100, m0_vlbi from the IGNITE-2 bank) |
| stage-1 array | **12509391** (8 shards, %6): 780/780 signal realisations |
| pairs | **12509392**: 40/40 paired realisations (member0-inert rescore columns banked) |
| sloop array | **12510901** (4 shards): 30/30 loops. First attempt 12509837 OOM'd (two B1Marg evaluator sets in one process) → one-shape-per-process + T_CHUNK 8; second attempt caught the tied-frame decode bug (§2b note in `chorus.py`) — the 11 affected m1e07 loops were **deleted and rerun**, never read |
| loop cost | 16-slot ~2 min; 47-slot 7–28 min (T_CHUNK = 8) |
| disk | ~4 900 lean npz ≈ 160 MB in `CHORUS_results/` + 3 figures + 2 summary npz — no budget concern |
| device log | every job logs GPU UUID + memory + foreign-process line; squat lottery drew clean on every task |

**Nothing was committed. Nothing was pushed. No tracked file was edited.** `CHORUS_mixed_pop.md`,
`CHORUS_results/`, `hpc_harbor/chorus/` are untracked. **STOP** (fourth mode, per the 2c
pre-registration, with the anatomy in §2c: both STOP events are criterion-seed events; the
loop generated neither).
