# LOTTERY — tier 2: validating the channel-budget switch-on proxy on the GPU

**Agent LOTTERY, ACCRE, 2026-07-17.** General-queue H200 (`batch_gpu`, `taylor_group_acc` /
`taylor_group_account_batch_gpu`), fair-share, **cpus-per-task=8**. Reserved `dgx_iacc` share and
every SPARK-3 job untouched. Staged, not committed.

> **Budget (measured).** 8 cells × (100 nullN + 20 sig) + on-demand basis controls =
> **960 realisations, 2 659 GPU-s ≈ 0.74 GPU-hr** on one H200 — **~1 %** of the ~60 GPU-hr cap.
> No trim needed. Per-realisation cost matches CHORUS (~0.5–0.8 s warm; ~450 s build+compile per
> `n_ecc` shape). Checkpoint-per-draw (keyed `lot_*.npz`); **resume-by-skip verified** (re-run
> after the 720-realisation batch: "already banked 720; to run 0").

---

## 0. WHAT WAS RUN, AND THE ONE HAZARD THAT SHAPED IT

Tier 1 flagged a disagreement band — mixtures with `n_active ≥ 30` (channel budget: **switch ON**)
whose eccentricities never trip the quoted threshold rule (single ≥ 0.5 / pair ≥ 0.3: **OFF**). The
proxy under test is **`n_active ≥ 30` ⟺ the certified count clears the >1 bar**. Each cell has a
deterministic `n_active`, so the analytic channel verdict is ON/OFF and the GPU measures the grade
(CONFIRMED / MARGINAL / below) under **criterion-v2.2, floors refit per cell** — CHORUS's exact
machinery (`chorus.run_realisation`) + the RECUT floor/score/grade convention verbatim
(empirical q95 + bootstrap where zero-fraction > 20 %; **all cells qualify**).

**The GWB-basis hazard, and how it was handled.** At T = 30 the extended GWB square root is 5336²,
but the frozen `b1_L_gwb` bank is the T = 15 baseline (3248²) — it shape-mismatches and NoiseDrawer
would raise. CHORUS's own T = 30 banks were drawn (2026-07-13) *before* that bank existed, via the
recompute-at-8-threads fallback. **LOTTERY reproduces that exact convention** (forces the recompute
path); all 8 cells were drawn on **one host (hgx03), one basis**, so every cell's floor and count
share it. This makes the 8 cells mutually comparable; comparison to CHORUS's dgx-banked numbers is
what the **basis control** (§2) quantifies.

---

## 1. THE REFIT MACHINERY IS VALIDATED — and the proxy is NOT a monotonic switch

The full same-basis (hgx03) calibration, sorted by `n_active` (figure `LOTTERY_tier2_calibration.png`,
panel A):

| cell | `n_active` | chan pred | thr pred | floor ± bs | count | @fl+bs | **grade** |
|---|---|---|---|---|---|---|---|
| m1 e=0.25 | 22 | OFF | OFF | 8.36 ± 0.67 | 1.40 | 1.20 | **CONFIRMED** |
| m1 e=0.30 (control) | 23 | OFF | OFF | 9.91 ± 0.86 | 1.20 | 1.05 | **CONFIRMED** |
| m2 e=0.25 | 28 | OFF | OFF | 11.60 ± 1.15 | 1.55 | 1.00 | MARGINAL |
| m1 e=0.45 | 29 | OFF | OFF | 9.20 ± 0.76 | 2.75 | 2.35 | **CONFIRMED** |
| m2 e=0.28 | 30 | **ON** | OFF | 9.79 ± 0.50 | 2.60 | 2.50 | **CONFIRMED** |
| m2 e=0.30 (control) | 30 | ON | ON | 9.93 ± 0.49 | 2.80 | 2.40 | **CONFIRMED** |
| m3 e=0.20 | 31 | **ON** | OFF | 10.36 ± 0.94 | 2.00 | 1.50 | **CONFIRMED** |
| m3 e=0.25 | 34 | **ON** | OFF | 13.65 ± 1.49 | 1.50 | 0.90 | MARGINAL |

**Two things fall out immediately:**

1. **The channel-budget `n_active ≥ 30` switch is not a switch.** On one common basis, *every* cell
   from `n_active = 22` to `34` clears the >1 bar at its floor (counts 1.20–2.80); 6 of 8 are
   CONFIRMED. The sub-30 cells the proxy calls **OFF** — 22, 28, 29 — **certify**. And the
   **highest**-channel cell, m3 e=0.25 (`n_active = 34`), is only **MARGINAL** — because its refit
   floor (13.65) is the highest in the set. **The count does not track `n_active` monotonically.**

2. **The controlling variable is the refit FLOOR, which the proxy ignores entirely.** The floor rises
   with eccentric structure (more active harmonics → a bigger fringe-template bank → a higher
   spurious null offender), and the grade is the race between count and floor. m3 e=0.25 has the most
   channels *and* the most signal, but its floor outran the count. **The channel budget models the
   numerator (registration channels) and omits the denominator (the floor). It cannot predict the
   grade.**

**Direction of the proxy's failure:** too **conservative** at low `n_active` (predicts OFF at 22–29,
measures ON) and too **optimistic** at the top (predicts CONFIRMED at 34, measures MARGINAL).

**The threshold rule is worse still:** it calls all six new cells OFF; all six clear the bar.

---

## 2. THE BASIS CONTROL — why the near-bar grade is fragile, and a headline that flips

CHORUS's own m1e03 / m2e03 cells, **re-derived on this host** (fresh seeds, hgx03 basis) vs their
banked dgx values (panel B):

| cell | floor dgx → hgx | count dgx → hgx | grade dgx → hgx |
|---|---|---|---|
| **m2 e=0.30** | 10.05 → **9.93** (Δ −0.12) | 2.77 → **2.80** (Δ +0.03) | CONFIRMED → CONFIRMED ✓ |
| **m1 e=0.30** | 11.30 → **9.91** (Δ **−1.40**) | 0.70 → **1.20** (Δ +0.50) | **below → CONFIRMED** ✗ |

- **m2 e=0.30 reproduces the banked numbers to 0.12 nat / 0.03 count.** The refit path here **is**
  the one that produced the published CHORUS counts — the scorer is correct, and cells well above the
  bar are host-robust.
- **m1 e=0.30 — the RECUT headline cell ("a single e = 0.3 member does NOT switch the count on")
  — FLIPS.** Its refit floor drops 1.40 nat across hosts and the count rises through the bar:
  **below → CONFIRMED.** The 1.40-nat move is ≈ 1 σ of the floor's own combined bootstrap error
  (0.86 ⊕ 1.02 ≈ 1.33) — i.e. **re-estimating this floor from a fresh 100-null draw, on a different
  host, moves it by about its own sampling error, and that is enough to reverse the verdict.**

> ### **THE CONSEQUENTIAL FINDING (for the CHORUS / RECUT owners).**
> **RECUT's "single e = 0.3, lit, REFUTED" — called "the single most expensive consequence of the
> floor defect in the whole repo" — sits on the floor's own error knife-edge. It is a below/CONFIRMED
> coin-flip at the ~1-nat level, NOT a settled refutation.** m2 e=0.30 (well above the bar) is solid;
> the e = 0.3 *single-member* result is not. This does not touch CHORUS's structural results (clock
> not shared, the 10× lever, the two-regime channel statement) — only the one near-bar grade.

---

## 3. THE VERDICT ON THE PROXY — the surface is NOT licensed at analytic cost

The mission's decision rule: *if the proxy tracks the measured grade within errors, the whole
P(switch-on) surface is licensed analytically; if it diverges, report where and which direction.*

> ### **It diverges. The analytic channel budget does NOT track the certified-count grade.**
> The `n_active ≥ 30` switch is not reproduced — sub-30 mixtures certify, the top-channel cell is
> suppressed — because certification is gated by the **refit floor**, which (a) rises with eccentric
> structure and (b) carries a ~1-nat sampling/host error the proxy omits. **A P(switch-on) forecast
> cannot be read off `n_active` alone; it must carry the per-cell floor and its error.**

**What IS licensed:** the *structural* reading — census-loudness mixtures with modest eccentricity
(e ≳ 0.2, one or more members) do reach the certification regime, and the refit **machinery** is
sound (m2 e=0.30 control). **What is NOT:** the sharp threshold, the monotonicity in `n_active`, and
any near-bar grade to better than ~1 nat of floor error. The tier-1 P(switch-on) **surface** stands
as the *channel-budget* statistic it is defined to be; it is **not** a calibrated predictor of the
GPU grade, and the two must be quoted as different objects — exactly the tier-1 lesson (the proxy is
one statistic, the graded count another) now measured, not asserted.

### The tier-1 validation column, filled

| tier-1 cell | analytic channel | measured grade (hgx03) | tracks? |
|---|---|---|---|
| V1 m3 e=0.25 (`n_active`=34) | ON → predict CONFIRMED | **MARGINAL** | ✗ (floor too high) |
| V2 m3 e=0.20 (31) | ON → CONFIRMED | **CONFIRMED** | ✓ |
| V3 m2 e=0.28 (30) | ON → CONFIRMED | **CONFIRMED** | ✓ |
| V4 m2 e=0.25 (28) | OFF → below | **MARGINAL** | ✗ (certifies at floor) |
| V5 m1 e=0.45 (29) | OFF → below | **CONFIRMED** | ✗ (certifies) |
| V7 m1 e=0.25 (22) | OFF → below | **CONFIRMED** | ✗ (certifies) |
| V6 m3 e=0.70 (109) — carried | ON → CONFIRMED | CONFIRMED (banked) | ✓ |
| control m1 e=0.30 (23) | OFF → below | **CONFIRMED** / below (dgx) | host-dependent |
| control m2 e=0.30 (30) | ON → CONFIRMED | CONFIRMED (both hosts) | ✓ |

**3 of 6 new cells contradict the channel verdict**, all in the sub-threshold band — the proxy is
systematically too conservative there, and the one place it predicts a clean ON at high `n_active`
(34) it overshoots. The concordant-high anchor and the above-bar control hold; the near-bar cells do
not.

---

## 4. WHAT THIS SETTLES, AND WHAT IT HANDS BACK

1. **The channel-budget proxy is not a calibrated switch-on predictor.** The P(switch-on) surface
   (tier 1) is the channel-budget statistic and must not be quoted as the certified-count switch;
   the two diverge, one-sidedly at low `n_active` and in the highest-channel corner.
2. **Certification is floor-gated.** The floor rises with eccentric structure; the highest-channel
   mixture can grade *below* a lower-channel one. Any forecast must carry the floor.
3. **The near-bar certified-count grade is fragile at ~1 nat** — floor sampling/host error flips
   below↔CONFIRMED. **RECUT's single-e=0.3-REFUTED headline flips on a second host** and should be
   re-examined against floor error (m2 e=0.30 and all structural CHORUS results are unaffected).
4. **The refit machinery is validated** (m2 e=0.30 reproduces the banked count to 0.03), so findings
   1–3 are properties of the physics + estimator, not a code discrepancy.

**Held for Matt:** whether the docs' P(switch-on) language should distinguish the channel-budget
surface from the (uncalibrated) certified-count switch.

> **RESOLVED — the m1 e=0.3 re-cut was decided and executed: see `LOTTERY_recut_band.md`.** A
> 5-state GWB-basis ensemble (dgx banked + hgx03 at 2/4/8/16 BLAS threads, four distinct `L_gwb`
> fingerprints) grades the cell **MARGINAL (band): floor 9.91–11.30, count 0.70–1.20** under a
> pre-registered two-sided rule — it neither confirms nor refutes, and the external quote
> (e = 0.5 single / e = 0.3 pair) is unchanged. Two refinements to §2 above: the flip is **not** a
> two-sided-error effect (at fixed basis the cell stays `below` at 0.97), and **basis rotation is
> benign** — across four bases on one host the count is identical (1.20 ×4) and the floor moves
> 0.18 σ; the systematic is the **host**, and it is **not** a common-mode offset (−11.4 % on this
> cell, +0.3 % on m2 e=0.3). Converges with SPARK-3 §4.4.
```
add-list (stage, do not commit):
  LOTTERY_tier2.md
  LOTTERY_tier2_calibration.png
  reports/lottery_tier2.npz
  hpc_harbor/lottery/lottery_tier2.py
  hpc_harbor/lottery/lottery_tier2_figure.py
  hpc_harbor/lottery/lot_tier2.sbatch
  (raw realisation banks: LOTTERY_results/ — 960 lean npz, ~80 MB — on ACCRE disk, GITIGNORED by
   the repo's `*_results/` rule, exactly as CHORUS_results; the committed bank is reports/lottery_tier2.npz)
```
**STOP.** Tier 1 + tier 2 complete; validation set executed, calibration measured, proxy verdict
delivered. No further submits without a fresh brief.
