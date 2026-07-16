# SPARK — the cascade's launch criterion: the pulsar-term-coherent detector, built and priced

**VERDICT: CASCADE-ALIVE.** The arithmetic is **not** closed. Coherence buys **2.96×** in the
reservoir's detection statistic and simultaneously **lowers** its null floor **2.6×**; both
push the same way, and recruitment is measured: **4/13 → 7/13** reservoir sources cross at an
*achievable* certified set, **4/13 → 12/13** at the coherence ceiling. The cost side (fringe
trials growth for the grown list) is **+0.578 nat** — an order below the gain.

**This is a statement about the ARITHMETIC, not a green light for the loop.** SPARK's
statistic is oracle-anchored and carries no trials factor; B0.2's search gap
(`project_progress.md:1450-1462`) is untouched and still blocks a self-found cascade. §5 is
the load-bearing caveat and §6 is the design that the ALIVE verdict actually licenses.

**Tree.** `MM_playground` @ `d87db93` (= `criterion-v2.2` + 1, the KINDLE stage-0 commit;
see AVALANCHE §0 for why the tag is not checked out detached). AVALANCHE's STOP is accepted
and stands; `reports/AVALANCHE_cascade.md` is the pre-flight record. SPARK executes its
proposed successor (AVALANCHE §7).

**REALISM.** Realistic spine — real pulsar positions, real TOA uncertainties, real published
distance priors, drawn white + per-pulsar red + GWB noise, honest per-cell null-calibrated
floors, Arm-B truth. **The spine is still a MOCK** (`data-spine-is-simulated`; AXIS, 1440 MHz
single-frequency, §10.15(a)): no real TOA is touched and the residuals *are* the injected
CW+CURN. A verdict here is a statement about an approximately-real IPTA's geometry and noise
budget.

**SCOPE OF INFERENCE.** One cell (h = −13.25, T = 30, lit, fiducial POP sky), **one signal
realisation × 13 reservoir sources**, floors on **100 nulls × 13 positions = 1300 draws** per
state. The floors are well-sampled; **the recruitment counts are a single noise draw** and are
not an ensemble — §5(d). L_gwb fingerprint **`9fd547b39b02c705`** on every job = the banked
CHORUS T=30 value; `cpus-per-task=8` throughout.

**COST.** 3 GPU jobs, **≈75 GPU-min total** (g0 1006 s; s2c 2820 s; plus the two superseded
runs, §4). Against AVALANCHE's ~80 GPU-hr STOP and the ~400 GPU-hr the full campaign would
have cost. The launch criterion was, as forecast, a few-GPU-hr experiment.

---

## 1. S1 — THE BUILD: the missing msd path, wired into the detector

AVALANCHE's ground (2) was that the cascade's `[E]→[D]` edge does not exist in code. Confirmed
and now repaired.

**The gap, precisely.** The incumbent seeder is the *fully-decohered limit by construction*:
`TE.seed_scan` (`:228`) → `TE.build_earth_single` (`:182`) → `build_fisher_amortised(...,
msd_factory=TE.EarthDelay)` (`:188`), and `TE.EarthDelay` (`:56`) is
`make_phase_connected_binary(pulsarterm=False)` — docstring: *"Earth-term-only … the
fully-decohered coherence limit"*. Certified pulsar terms **cannot enter it**.

**What was missing was the wiring, not the physics.** `trackB_b1_core.MaskedDelay` (`:92`)
already implements

```
delay_p = d_earth + m_p * (d_full − d_earth),     m_p = params[PMASK][ipsr]   # RUNTIME
```

with both limits gated (`m_p=1` ≡ `MultiSourceDelay(pterm=True)`, 0.00; `m_p=0` ≡
`TE.EarthDelay`, 0.00 — `trackB_b1_core.py:754-758`). It had never been given to a *detector*.
SPARK builds the ncw=1 coherent twin of `build_earth_single` via `C.build_b1_amortised`, so the
certified set enters through **`(pmask, Ldist)` as runtime args — one compiled graph for every
certified set** (`hpc_harbor/spark/spark.py`, `build_detectors` / `make_fstat_coh`).

**Soft / q-weighted per spec §3, never hard-locked:** `m_p = q_p ∈ [0,1]`, the pulsar's own
certification confidence — not a 0/1 lock. Uncertified pulsars sit at `m_p = 0`: **decohered,
not pinned to a wrong MAP fringe** (the hard-lock failure IGNITE-2 retired).

**A latent defect in the incumbent, reported not fixed:** `TE.build_earth_single` hard-codes
the DEFAULT GP counts (`N_COMPONENTS=14`, `rn=30`) regardless of T, and `TE.seed_scan` (`:233`)
hard-codes `T = 6.992e8` (the T=15 span) *"to avoid recompute"*. At T=30 the problem's own
counts are span-scaled (23/50) and the span is 37.14 yr, so the incumbent's noise model and
frequency grid do **not** match a T=30 problem. The seeder was only ever run at T=15, where
both are no-ops. SPARK builds both detectors at the **problem's** counts and uses the **real**
span.

### THE GATE CHAIN — all three pass

| gate | what it proves | result |
|---|---|---|
| **g0a** | SPARK's re-wired EarthDelay F-stat (data as runtime arg) ≡ `TE.make_fstat` verbatim — the re-wiring did not move the incumbent | **0.000e+00** (512 pts) |
| **g0b** | **THE GATE**: the coherent detector at **zero certified pulsars** reproduces the EarthDelay F-stat | **0.000e+00** on the statistic **and** 0.000e+00 on the profiled free params, full 14 976-pt grid |
| **g0c** | the data-driven loud list from the coherent detector's own `pmask=0` map, selection TRUTH-BLIND (F-stat + sky-exclusion NMS) | **3/3 loud captured** (6.37°, 6.39°, 16.67°; 2F 1524.3 / 1200.7 / 329.9; 40 seeds) |

g0b is **bit-identical**, better than the adopted bar (EMBER §2.2b: discrete exact, continuous
< 1e-6 — bit-identity was *not* demanded, since the two paths are different builders and
different graphs). `ll0` agrees to 0.000e+00 (939877.753465 both). Banked
`SPARK_results/spark_g0.npz`.

> **The detector is built and gated. AVALANCHE's ground (2) is CLOSED: `[D]` can now see the
> certified set.** Grounds (1) (the B0.2 search gap) and (3) (cost) are untouched — and (1)
> is why ALIVE is not a green light.

---

## 2. S2 — THE ARITHMETIC: the oracle ledger

**The statistic.** For each of the 13 reservoir sources k (log10_h = −14.25, exactly 1 dex
below the loud; `trackB_b1_core.py:85` `POP=(3, −13.25, −14.25)`):

```
2F_k(pmask) = 2·[ max_{cos_inc, log10_h, phase0, psi} logL(θ_k(free), data, pmask)
                  − logL(θ with k REMOVED, data, pmask) ]
```

with **sky, freq, and mc ORACLE-ANCHORED at truth**, **all other 15 sources at truth in the
model** (the 3 loud *and* the other 12 faint subtracted), the **same Adam profile as
`TE.make_fstat`**, and the floor **null-calibrated at matched state** (reservoir absent from
the data; RECUT `adopt()`; zero-fraction a required column).

**States.** `s0` = pmask 0 (the incumbent's fully-decohered limit). `sC_g` = q-weighted on
igniter g's banked certified set, certified distances at `L_true` (Arm B: certification is
exactly what pins `L_true = L0 + (n+u)·dL`). `sMAX` = **all 116 coherent at `L_true`** — the
ceiling; unreachable by construction.

### THE LEDGER (every column re-derived from the raw bank, `spark_readback.py`)

| state | N_cert | 2F med | 2F max | floor | ± | est | zf | bar | **clears** | gain (med) |
|---|---|---|---|---|---|---|---|---|---|---|
| **s0** (decohered) | 0 | 13.709 | 36.875 | 19.498 | 0.391 | gumbel | 0.000 | 19.889 | **4/13** | — |
| **sC_m1e07** | 14 | 17.226 | 46.705 | 16.300 | 0.322 | gumbel | 0.002 | 16.622 | **7/13** | **+1.936** |
| **sC_m1e05** | 10 | 12.990 | 49.335 | 15.749 | 0.308 | gumbel | 0.000 | 16.057 | **6/13** | −0.204 |
| **sC_m2e03** | 10 | 16.778 | 45.244 | 15.886 | 0.308 | gumbel | 0.002 | 16.194 | **7/13** | +0.739 |
| **sMAX** (ceiling) | 116 | **40.575** | 76.552 | **7.491** | 0.133 | gumbel | 0.008 | 7.624 | **12/13** | **+14.291** |

**RECRUITMENT vs s0: m1e07 +3, m2e03 +3, m1e05 +2, ceiling +8.**

All five floors re-cut from the banked raw nulls to **|dev| = 0.00e+00**. All five zero-fractions
are ≤ 0.008 → Gumbel valid at the ZF_GATE = 0.20; the emp-q95 branch is banked beside them
anyway.

### THE COST SIDE

**lnK growth per added source** (exact, CPU — `dL` is the MIN over the source list of the mode
spacing (`forge_b1.apply_geometry:91-94`), so a longer list can only shrink `dL` → more fringes
in the prior window → `K_counted` grows):

| n_src | K_sum | lnK_med | dL_med |
|---|---|---|---|
| 1 | 348 330 | 7.164 | 7.795e-04 |
| 3 (the loud) | 581 900 | 8.183 | 3.814e-04 |
| 4 (+1 recruited) | 622 702 | 8.291 | 3.642e-04 |
| 16 (the full census) | 993 132 | 8.761 | 2.167e-04 |

> **Growing the list 3 → 16 costs +0.578 nat** on the median pulsar's trials bar (K_sum ×1.707;
> the first recruit alone costs +0.108 nat). **Against a coherence gain of +14.29 (2F) = +7.15
> nat at the ceiling, and +1.94 (2F) = +0.97 nat at an achievable certified set.**

**THE ARITHMETIC IS NOT CLOSED. The gain exceeds the cost — at the ceiling by ~12×, and even
at an achievable certified set by ~1.7× on the first recruit. → CASCADE-ALIVE.**

---

## 3. THE MECHANISM — and a finding nobody predicted

Coherence helps through **two channels, both favourable, and the second is the surprise**:

**(a) The signal rises.** `sMAX` median 13.709 → 40.575 = **2.96×**. This is the expected
physics: a certified pulsar contributes a *second coherent copy* of the source (the pulsar
term at `L_true`), so the matched-filter power roughly multiplies.

**(b) The FLOOR FALLS — 19.498 → 7.491 (2.6×).** *Not* predicted, and it is half the effect.
With pulsar terms on at known distances the template is **rigid**: it must match a specific
two-component (Earth + pterm) structure at a fixed lever arm, and noise is far less able to
mimic it under the profile. **Coherent detection is more selective, not merely louder.**

**(c) At achievable certified sets the FLOOR channel dominates, not the signal.** The decisive
case is `sC_m1e05`: its signal median is **lower** than s0's (12.990 vs 13.709, gain **−0.204**)
and it still recruits **+2**, purely because its floor fell 19.498 → 15.749. With only 10 of
116 pulsars coherent the signal gain is in the noise, but the template-rigidity gain is not.
**Anyone pricing this cascade on signal gain alone would have called it dead.**

**The lever is per-pulsar and brutally sublinear in what is currently reachable.** The ceiling
(116/116) buys +14.29; the *most-certified realisation ever banked* (14/116) buys +1.94 — about
14% of the ceiling for 12% of the array, i.e. roughly linear in N_cert, with no super-linear
kick. **The channel-budget lever (J1, `n_active` ≈ 27) sets how many pulsars certify; it does
not amplify what each certified pulsar is worth.** The two levers multiply, they do not compound.

**Certified-set re-derivation (readback gate, unplanned and it passed).** Re-derived from raw
`(dlnL_det, lnK, qmax)` at the adopted criterion-v2.2 floor, my `N_cert` mean reproduces the
banked `corr + wrong` **exactly** for all three igniters: m1e07 5.500 = 5.433+0.067; m1e05
3.167 = 3.133+0.033; m2e03 2.800 = 2.767+0.033 (|dev| = 0.00e+00). Wrong-cert **rate** 1.05–1.21%
— the certified sets SPARK coheres on are ~99% correct counterparts, so the coherence gain
measured here is **not** contaminated by D1.

---

## 4. TWO SUPERSEDED RUNS — and how each was caught

Both ran, both are banked as trails, and **both were caught by reading their own output, not by
a gate I had pre-registered.** Recording them because the second one nearly produced a false
headline.

**(a) `spark_s2_fstat_trail.npz` — the ncw=1 F-stat ledger. Its `CASCADE-ALIVE` is RETRACTED.**
Two defects, either fatal:
- **Loud-sidelobe leakage.** The ncw=1 detector models ONE source against data holding 3 loud
  + 13 faint, so at a faint position the statistic is dominated by the **loud sources'
  sidelobes** — precisely the pathology TE's sky-exclusion NMS exists to defeat
  (`trackB_estimator.py:412-415`: *"loud2's true peak sank to rank-13"*). The tell was the
  **null floor: 145.3**, where a χ²₄-like statistic belongs near 9.5 — because the null
  *retains* the loud sources. Its single "clearing" faint source lay near a loud one: the
  detector was seeing a loud sidelobe, and my verdict logic scored that as recruitment.
- **mc mis-specification.** `TE.make_fstat` hard-codes `SEED_MC = 9.0` (`TE:163`). Harmless for
  the Earth term (the frequency barely evolves over T); **fatal for the pulsar term**, which
  looks back ~L/c (kyr) into the binary's past where the frequency was materially different.
  **Measured: coherence at SEED_MC LOWERED the statistic** (sMAX median 15.9 vs s0's 21.1) —
  the coherent model actively mismatches the data. That is a real and reportable fact about a
  coherent detector that does not know mc, but it is **anti-cascade**, and it broke SPARK's
  a-fortiori structure, which requires every convention to *favour* the cascade.

**(b) `s2b` — the oracle dlnL ledger. Killed by its own null.** It used `dlnL_k = logL(k in the
model AT TRUTH) − logL(k removed)`. Under the null (k absent from the data) inserting k at its
**true amplitude** can only *lower* the likelihood → the null collapsed to an all-≤0 sample: a
**100% zero point mass**, no floor, and a `scipy` `OverflowError` in the Gumbel fit. That is the
correct behaviour of an incorrect object: **dlnL-at-truth is not a detection statistic — it
already knows the true amplitude, which under the null does not exist.** A detector must
profile the unknown amplitude/phase; only then does the null acquire the χ²-like spread a floor
can be cut from. (The guard added to `gumbel_floor()` now returns NaN on a degenerate sample
rather than crashing — the zero-fraction gate routes any such cell to emp_q95 regardless.)

> The lesson SPARK would hand forward: **the null is the diagnostic.** Both defects announced
> themselves as an absurd floor (145.3) or an impossible one (a 100% zero mass), and neither
> would have been visible from the signal arm alone. A campaign that null-calibrates only at
> the end learns this too late.

---

## 5. WHAT "ALIVE" DOES **NOT** MEAN — the caveats that travel

**(a) ORACLE-ANCHORED, NO TRIALS FACTOR. The absolute `clears` counts are ceilings, not
forecasts.** 2F is read at each reservoir source's **true** sky/freq/mc, and the floor carries
**no search penalty**. A real detector must scan the grid (14 976 cells here, and a real search
is far finer); the trials factor would raise the bar substantially — plausibly into the low
30s for a χ²₄-like statistic at this grid size, which is above s0's best (36.9) only marginally
and would erase most of the `s0` 4/13. **The robust readout is the RELATIVE one — the gain and
the recruitment DELTA — not "4/13 reservoir sources are detectable today."** They are not.

**(b) THE SEARCH GAP IS UNTOUCHED, AND IT IS STILL THE BINDING CONSTRAINT.** B0.2
(`project_progress.md:1450-1462`) measured cold-start source recovery at **~0.5 scaled** against
a certification tolerance of **~1e-4 scaled** — 3–4 orders, at **zero noise**, failing
**confident-wrong** (q_max 0.5–0.99 on a shifted fringe). The registration needle is a cusp of
width **~1e-5°**; the seeder delivers **6–12°**. SPARK proves the *array* gets better when
pulsars certify. It proves **nothing** about whether a blind search can find what to certify
against. **A self-found cascade remains blocked**, and EMBER's boundary (wrong counterpart ×
motion, p = 0.002) still says a loop that moves under self-found counterparts manufactures.

**(c) THE ACHIEVABLE STATES ARE DELIBERATELY OPTIMISTIC.** `sC_g` uses the **most-certified
realisation ever banked** (N_cert = 14 / 10 / 10), not the mean (5.50 / 3.17 / 2.80). At the
*mean* certified set the gain is roughly a third of what is tabled, and `sC_m1e05`'s signal gain
is already negative. **The honest per-realisation expectation is recruitment of order +1, not
+3.**

**(d) N = 1 ON THE SIGNAL SIDE.** The floors rest on 1300 null draws each; **the recruitment
counts rest on a single noise realisation** × 13 sources. 4/13 → 7/13 is one draw. The *sign*
of the effect is secure (the ceiling's +8 and the 2.6× floor drop are far outside the floors'
errors, ±0.13–0.39); **the integer counts are not.** An ensemble is the first thing EPOCH
should buy.

**(e) ONE CELL.** h = −13.25, T = 30, lit, fiducial sky. No T=40 rung, no VLBI rung. The VLBI
two-sided reading (ANCHOR §7: tighter priors → smaller K → lower lnK bar → **higher** null
floors, ≈ +2.9 ± 1.0 nat at h = −13.25, and nothing measurable elsewhere) is **not** tested here
and remains the successor's question.

---

## 6. S3 — THE PRE-REGISTRATION (the ALIVE branch): the EM-MEDIATED loop

> **EPOCH — pre-registration.** SPARK licenses a cascade campaign **only in its EM-mediated
> form**, because that form closes by construction the one hole SPARK does not touch. The
> measured arithmetic is favourable — coherence buys 2.96× in the reservoir's detection
> statistic and drops its floor 2.6×, against a +0.578-nat trials cost for the grown list — but
> the recruitment step must never be **self-found**: a self-found counterpart is wrong by 3–4
> orders relative to the certification tolerance (B0.2), fails **confident-wrong** so the q>0.9
> layer cannot catch it, and a loop that moves under a wrong counterpart is exactly EMBER's
> manufacturing regime (motion sensitivity 1.00, p = 0.002; D1 open, D3/D4 rejected). The
> EM-mediated design replaces the self-found step with an **externally verified** one:
> **[D] the certified-coherent detector (SPARK's, gated g0a/g0b/g0c) flags a candidate over its
> per-iteration floor → [L] the loop's own sharpened localisation defines a sky error box →
> [A] DIRECTION-A COUNTERPART IDENTIFICATION: an external EM catalogue (host galaxy / AGN) is
> queried over that box; a candidate is admitted to the source list ONLY on a counterpart match,
> and the counterpart's catalogue position — not the loop's estimate — becomes the source's sky
> prior → [M] fit → [E] re-score and certify under v2.2 → repeat.** The wrong-counterpart mode
> is closed **by construction**, not by a purity layer (none exists): the admitted position is
> externally measured, so motion under it is *repair*, not manufacture — EMBER's boundary is
> respected rather than tested. This makes the campaign's governing quantity neither the channel
> budget nor the coherence gain but the **LOCALISATION–CONFUSION EXCHANGE RATE**: the area the
> loop must reach before the EM box contains ≲1 candidate host, versus the area coherence
> actually buys. **EPOCH's launch criterion is therefore a single cheap measurement, not a
> campaign: does the certified-coherent loop's localisation area for a *recruited* reservoir
> source shrink below the counterpart-confusion area of a real catalogue at that redshift?** If
> yes, the full ladder (T=40 × VLBI × the igniter set, with per-iteration floors, an ensemble
> ≥10 per cell to replace SPARK's N=1 signal arm, and the trials factor restored to the
> detection bar) is licensed and should be costed against AVALANCHE §0. If no, the cascade is
> **alive in the array and dead in the sky**, and that sentence — not a null loop result — is
> the paper's closure.

**Why this and not a purity layer.** D3 and D4 were both rejected; no purity layer exists, and
SPARK does not change that. EM mediation does not *detect* wrong counterparts — it **prevents**
them from entering, which is the only move available when the failure mode is confident-wrong.

---

## 7. WHAT SPARK CHANGES IN THE RECORD

1. **AVALANCHE's ground (2) is CLOSED.** The `[E]→[D]` feedback edge now exists, gated
   bit-identically at its decohered limit. Grounds (1) and (3) stand.
2. **"Coherence can at most double the power" is REFUTED as a bound.** Measured 2.96× in 2F,
   because the floor falls as well as the signal rising. My own pre-flight estimate (§AVALANCHE)
   that coherence buys ≤2× was **wrong**, and wrong in the cascade's favour — it priced only the
   signal channel and missed template rigidity entirely.
3. **The cascade is not arithmetically closed.** AVALANCHE's proposed kill-shot — "if the gain
   is smaller than the cost, nothing else need be built" — **did not fire**. The gain exceeds
   the cost at every certified set tested. The successor is licensed, in the EM-mediated form.
4. **The binding constraint is relocated.** It is not the arithmetic and not the channel budget.
   It is the **search/localisation gap** (B0.2), and EM mediation is the only measured way past
   it that does not require closing D1 first.
5. **`TE.build_earth_single` / `TE.seed_scan` are T=15-only** (hard-coded GP counts and a
   hard-coded 6.992e8 span). Reported, not fixed (HARD RULE). Any future T≠15 seeding must
   build at the problem's counts, as SPARK does.

---

## 8. FILES

```
hpc_harbor/spark/spark.py            S1 build + gate chain; s2/s2b/s2c ledgers (modes)
hpc_harbor/spark/spark_readback.py   CPU-only readback: every number re-cut from raw
hpc_harbor/spark/sp_g0.sbatch        gate      (dgx04, cpus=8, --exclude=dgx03, per-task cache)
hpc_harbor/spark/sp_s2.sbatch        superseded F-stat ledger (trail)
hpc_harbor/spark/sp_s2b.sbatch       superseded oracle-dlnL ledger (trail)
hpc_harbor/spark/sp_s2c.sbatch       THE ledger
SPARK_results/spark_g0.npz           gate evidence + both full-grid statistics
SPARK_results/spark_g0_grids.npz     grid checkpoint (resume)
SPARK_results/spark_s2c.npz          THE bank: raw 2F per (state, source) + raw nulls + lnK
SPARK_results/spark_s2_fstat_trail.npz   superseded, RETRACTED (§4a)
SPARK_results/spark_s2b.npz          not written (died in its own null, §4b)
```

Lean-npz throughout: **raw** 2F per (state, source) and raw null vectors are banked, never
verdicts — every floor in this report is re-cuttable on a CPU with no GPU, and
`spark_readback.py` does exactly that at **|dev| = 0.00e+00** on all five.
