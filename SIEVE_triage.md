# SIEVE — the triage battery

Two agents write into this file. **Sections merge by name; each agent appends only to its
own sections and never edits the other's.**

- **SIEVE-A** (banked-data / noisy half, ACCRE): **T2, T3, T6b, T7**
- **SIEVE-C** (deterministic half, cronus): its own sections

**REPORT-ONLY.** Nothing in the SIEVE-A sections arms a protocol step, moves a banked
verdict, or enters a closure claim. Bars-class findings are **posted and parked** per the
standing charter — they are stated with their number and left for Matt.

---

## SIEVE-A — STEP 0 RECORD

`git pull --ff-only` **did not run**: branch `glacier_lite` has no upstream
(`origin` has `master`, `MM_playground`, `add-sampling-notebook` only), so there is
nothing to fast-forward from. `git fetch --all` returned clean; HEAD is `ed78a1c`,
unchanged from session start. Read: `GLACIER_capstone.md` §S4.24–S4.24.2 (+ the
GLACIER-LITE addendum and the BASELINE lane claim), the criterion spec as *wired*
(`glacier_loop.py:533` for v2.2, `gl_v2_frontier.py` for frontier-v2,
`spark3.score_from_LNL` for the columns), and `LEDGER_stats_audit.md`, which landed at
02:13 mid-run and is cited where it bears.

**Lanes.** All four SIEVE-A tests ran on the **general CPU lane** (`batch`). The 0.5 GPU-hr
allocated to T6b was **not spent** — see §T6b. `SUMMIT`, `PHOENIX` (dgx03 A100-80GB),
`BASELINE` (dgx01 A100-40GB, `%4`) and the frozen-arm lane were not entered.
`GLACIER_results/` was read-only throughout; every write went to `SIEVE_results/` or
`reports/sieve_*.npz`.

---

## T7 — E-PROCESS SCOREBOARD

### VERDICT: **NO-GAIN** as a criterion component — KILL.
### But the run turned up a separate, strictly-dominating finding: **E0 — GLACIER's floor is cut on a different statistic from every other campaign's**. PROMOTE, posted-and-parked.

**Code** `hpc_harbor/ledger/sieve_t7_eprocess.py` · **bank** `reports/sieve_t7_eprocess.npz`
**Gate** re-derived criterion-v2.2 mask == banked `n_cert` on **706/706** banks (121 cells).

Three constructions were scored on every banked per-iteration history, all from the raw
columns (`dlnL_det`, `lnK`, `qmax`, `on_true`, `floor`, `zero_fraction`) with no floor
re-cut and no likelihood re-evaluated. **E1** is the brief's "vs the null floor
distribution": the banked pair (`zero_fraction`, `floor` = q95) two-point-calibrates an
exponential upper tail for the null max-offender, giving a conservative per-pulsar
p-value, converted by the standard κ=½ calibrator `e = 0.5·p^(−1/2)`. **E2** is
`exp(dlnL − lnK)`, the criterion's own trials bar on the e-scale. Bar throughout:
`e ≥ 1/α = 20`.

**The e-process does not discriminate manufacturing.** Of the 18 banked v2.2 certification
events (16 wrong / 2 on_true), E1 at α=0.05 leaves **2 survivors, both wrong**; E2 leaves
**13, all wrong**. **Both** of the campaign's true certifications die under both
constructions (r13p5/none/s3 psr8: E1 1.63, E2 5.23 at i0 and 2.47 at i4 — against a bar
of 20). On the pre-registered D2 population, E1 kills 11 of 13 and E2 kills 2 of 13, where
the **D2 rigidity gate kills 13/13 on R2**. The two E1 survivors are
`r13p25/e07/s0 psr62 @ i0` (E1 63.6) and `@ i5` (E1 33.3) — **the same pulsar that survives
LEDGER-B2's persistence rule** (`ledger_b2_persistence.npz`: 18 → 2 events, both wrong,
psr62 runs `[2,2]`). That convergence is the result: psr62 is manufactured with a fringe
peak-gap of 84.75 nat against `lnK` = 6.05 (E2 ≈ 1.5e34), and no statistic calibrated
against *noise* can touch it,
because the failure is **template misspecification**, not a fluctuation. R2's matched-power
test is a misspecification check and does something an e-process structurally cannot.
Adopting e-values for v2.3 would cost both true certs and buy 11 kills that D2 already
makes — **no gain**.

**A resolution limit worth recording.** For a cell whose nulls are silent
(`zero_fraction ≥ 0.95`, `floor = 0` — which is exactly the true-cert cell) the rule of
three caps the readable p at `3/n_null`, so with `n_null = 32` **no e-value above 1.63 is
obtainable from the bank at all**, regardless of construction. The true certs are not
"rejected" by E1 so much as **unresolvable** at the banked null resolution.

### E0 — GLACIER'S FLOOR IS NOT THE CAMPAIGN'S FLOOR

**Bars-class finding. Posted and parked; not acted on.** It is not in the capstone, not in
`LEDGER_stats_audit.md`, and not in any report in the tree.

The canonical offender statistic — `recut_surface.offender:75`, the function the floor
machinery was built around — is the **max over pulsars of `dlnL`**:

```python
def offender(dlnL, lnK, qmax):
    """Largest dlnL among pulsars passing layers 1+3. MAX OVER PULSARS -> Gumbel domain."""
    m = (dlnL > lnK) & (qmax > QBAR)
    return float(dlnL[m].max()) if m.any() else 0.0
```

**Every** other campaign in the tree uses that definition, verbatim or re-derived:
`anchor.py:361,469`, `chorus.py:605`, `ignite2.py:166`, `surface.py:370,448`,
`spark3.py:1063` ("RECUT recut_surface.py:75-78 verbatim"), `ember.py:378`,
`kindle_d7_fall.py:143`, `bank_surface_offenders.py:69`.

`glacier_loop._null_offenders:594` does **not**. It hand-rolls a different statistic:

```python
o = np.where(q_of > QBAR, np.maximum(dlnl - np.maximum(lnK, 0.0), 0.0), 0.0)
off.append(float(np.max(o)))
```

— the max over pulsars of **`(dlnL − lnK)₊`**. The gate is equivalent (`(dlnL−lnK)₊ > 0`
iff `dlnL > lnK`); the **returned value is smaller by `lnK`**. And the resulting floor is
then compared against `dlnL` (`glacier_loop.py:533`):

```python
cert[a] = (dlnL[a] > max(lnK[a] + 0.578, floor)) & (qmax[a] > QBAR)
```

So GLACIER's floors are **systematically low relative to the campaign convention by the
`lnK` of the argmax pulsar**, and the certification bar is permissive by the same amount.
Measured on the 18 banked certification events: median `lnK` = **6.08 nat** (6.32 nat over
the 16 events where the floor is actually the binding term, `floor > 0`); range 0.69–8.67.

**Measured directly on real null draws** (§T3, since GLACIER banks no null offender vectors
but GENERALISE does): cutting both statistics on the *same* 100 no-CW realisations of the
A-SKY survivor unit gives q95 = **59.47** (canonical, max `dlnL`) against **50.40**
(glacier_loop's `(dlnL−lnK)₊`) — a gap of **9.07 nat**. A bar cut with the `glacier_loop`
statistic and then applied to `dlnL` sits 9.07 nat low at that unit. This is a measurement
at the A-SKY venue, not at GLACIER's; for GLACIER the shift is bounded by its own `lnK`
column (median 6.08 nat at the certifying pulsars), and the 18→13 recut below is the
size of the effect actually visible in the GLACIER record.

This reads as an unintended divergence rather than a declared choice: `glacier_loop._stack()`
**imports the canonical `offender`** into its `FL` dict and then never calls it — the module
uses `FL["emp_quantile"]`, `FL["gumbel_floor"]` and `FL["adopt"]` (lines 517–520) and
nothing else. The right function is in the returned dict, unused.

Self-consistency can be restored from either end — cut the floor on `dlnL` (the canonical
function), or apply the existing floor to `(dlnL − lnK)₊`. **E0 does the latter**, because
it needs no null draws and no re-cut: `o[a] > max(0.578, floor)`, same banked floor, same
estimator, zero new compute. That gives:

| | events | wrong | on_true |
|---|---|---|---|
| criterion-v2.2 as wired | 18 | 16 | 2 |
| E0, floor on its own scale | **13** | **11** | **2** |

**Five wrong certifications are removed and both true certifications survive.** The five
are `r13p25/e07/s0` psr41@i0 (dlnL 30.48 vs matched bar 32.99), psr111@i0 (30.03 vs 30.26),
psr101@i3 (16.19 vs 23.00), psr1@i4 (17.69 vs 23.27), and `r13p5/e07/s0` psr43@i2
(9.44 vs 12.55). Unlike persistence (LEDGER-B2, which costs both true certs) and unlike
the e-process, **E0 is strictly dominating on the banked record**: it is the only tightening
tested here that removes wrong certs at zero cost in true ones.

**DESTINATION:** criterion-**v2.3** candidate term, alongside D2 rigidity. Matt's call —
it changes the meaning of a bar, so it is parked, and the 18→13 recut is quoted as a
measurement of the divergence, **not** as a proposed re-scoring of any banked verdict.

**Scope of the divergence, stated plainly.** It is confined to `glacier_loop`, so it
touches GLACIER and anything that inherits its floors (the Stage-1/Stage-2a cells, the
r13p9/r13p25/r13p5 ladders, and PHOENIX's frozen arm, which runs against those live cells).
Campaigns that call `recut_surface.offender` — ANCHOR, CHORUS, IGNITE-2, SURFACE, SPARK-3,
EMBER, KINDLE — are **not** affected. It does not touch the D2 rigidity verdicts (R1/R2 are
F-statistic tests that never read the floor), so **20/20 manufactured-set kills stand
unchanged**. What it does touch is every GLACIER statement of the form "cleared the floor".

---

## T3 — PARTIAL AMBIGUITY RESOLUTION vs CRITERION-v2.2

### VERDICT: **PAR-LOOSER**, decisively — and the reason is that its input is not calibrated. KILL for PAR; the calibration curve is the real product.

**Code** `hpc_harbor/sieve/sv_t3_par.py` · **bank** `reports/sieve_t3_par.npz`

The GNSS transplant is exact in form: the fringe integer is the ambiguity, the mode
spacing `dL` is the lattice, bootstrapped joint success over a subset is
`P_boot(S) = Π_{a∈S} P_a`, and PAR fixes the largest prefix (sorted by descending `P_a`)
clearing a level `P0`. **The full float covariance `Q` needed for a true LAMBDA
Z-transform is measured-unaffordable, not merely unbanked** — it is the same joint Hessian
SPARK-3 priced and refused (`spark3.faint_fisher_bounds`: two builds failed to return
inside a 1 h A100-80GB walltime). Declared substitution: `P_a := qmax_a`, the banked fringe
posterior, which `estep_per_target` already computes with every uncertified pulsar
**decohered** — the PTA analogue of decorrelation, and one that *removes* information, so
every joint rate quoted here is a **lower bound** on what a real LAMBDA fix would report.

**The calibration gate fails, and that is the finding.** On the A-SKY survivor unit
`AS_e03_h1275_k5_s4` (15 realisations × 116 pulsars = 1,740 samples). This is the
**cleanest possible case**, checked in the code rather than assumed: `generalise.py:373`
scores at `theta_base = theta_src` with `theta_src` the **drawn truth** (`gen_theta`), only
the distances moved to `L0` as the fringe machinery requires — **no M-step, no loop, no
wander, and the whole 16-source population present in the template**. Nothing is
mis-specified here except the E-step's own conditioning:

| qmax bin | n | mean qmax | empirical P(on_true) |
|---|---|---|---|
| [0.90, 0.95) | 64 | 0.9268 | 0.500 ± 0.063 |
| [0.99, 0.999) | 129 | 0.9963 | 0.620 ± 0.043 |
| [0.999, 1.000) | 1145 | 1.0000 | 0.853 ± 0.011 |
| **qmax > 0.9 (what v2.2 gates on)** | **1464** | **0.9942** | **0.786 ± 0.011** |

**`qmax` overstates the fringe-assignment success rate by 0.208 absolute in exactly the
region criterion-v2.2 uses — a ~19σ gap on 1,464 banked samples.** Every bin is
over-confident; none is calibrated.

This is the **empirical, banked-data confirmation of LEDGER's A1**, reached from the
opposite direction: LEDGER identified the mechanism at the lines (`estep_per_target`
evaluates the fringe posterior at ONE source point, so `q_max` is `P(fringe|θ̂)` charged as
`P(fringe|data)`) and measured it on a fringe toy (q 0.909 → 0.544 at half a fringe of
belief-induced peak motion). T3 measures the consequence on the actual banks and puts a
number on it. **The two are independent and they agree.**

The same read across **all 706 GLACIER per-iteration banks** (81,896 samples) is worse —
`qmax > 0.9` → on_true **0.305 ± 0.021**, and the most confident bin [0.999, 1.000) is
**0/90 correct**. Split by loop iteration it is **flat** (i0 0.233, i5 0.343, no trend), so
this is *not* accumulating M-step wander: it is present at the truth-anchored feed state.

**The two numbers are not the same estimator, and the gap should not be over-read.**
A-SKY scores through `sp.estep` (global pmask, whole population modelled); GLACIER scores
through `spark3.estep_per_target` with the carried census `H_ABSENT`. So the 0.786 → 0.305
difference bundles two things: the **census omission** (GLACIER carries 256 sources, ~250
of them unmodelled real signal sitting in the residual) and the **E-step variant** (the
per-target/global-pmask distinction already on record as the SPARK-2 `estep` confound).
What is common to both, and is the load-bearing claim, is that **`qmax > 0.9` does not
mean 90 % in either.**

**Consequences for PAR, and the head-to-head.** Because `P_boot` inherits the
miscalibration multiplicatively, PAR's predicted joint success is meaningless at these
subset sizes: at `P0 = 0.99` the safe subset averages **79.3 pulsars** with predicted
`P_boot = 0.992`, and the realised all-members-correct rate is **0/15**. Against
criterion-v2.2 on the same realisations:

- **census `r13p25/e07/s0` i0** (cascade): v2.2 certifies `{41, 62, 111}`; PAR at
  `P0=0.99` takes **14** pulsars, **PAR-LOOSER** by 11 — and **0 of the 14 are actually
  on_true**. Named PAR-only admissions: psr 24, 30, 39, 48, 51, 53, 54, 81 (all `q ≥ 0.997`).
- **census `r13p5/none/s3` i0** (the true certs): both select exactly `{psr 8}` —
  **AGREE**, and PAR's 1/1 is genuinely on_true.
- **A-SKY survivor** (pooled over 15 realisations, v2.2 scored at the **canonical floor
  59.47** cut from this unit's own 100 banked no-CW nulls — not at floor = 0): PAR takes
  **79.3 / 85.7 / 88.5** members at `P0 =` 0.99 / 0.95 / 0.90 against v2.2's **8.3**; mean
  PAR-only **71.0 / 77.4 / 80.2**, mean v2.2-only **0.0** at every level. **PAR is a strict
  superset in all 15 of 15 realisations.**

**Mechanism, named:** v2.2 carries an **amplitude** requirement (`dlnL` over the floor and
the trials term) that PAR has no analogue for; PAR carries only the multiplicity charge in
`Π P_a`. With `P_a` over-confident by 0.21, that charge under-bites, and the amplitude term
is the only thing standing between the criterion and PAR's 79-pulsar sets. **PAR does not
tighten anything here; it would loosen the criterion by ~71–80 pulsars per realisation at
the A-SKY unit and by 11–17 at the cascade cell.**

**By-product, feeding T7/E0.** Cutting the floor on the *same* 100 A-SKY nulls under both
offender definitions measures the `glacier_loop` divergence directly: canonical q95
**59.468**, `glacier_loop` statistic q95 **50.395**, gap **+9.072 nat**. Banked as
`asky_floor_canonical` / `asky_floor_glacier_stat` / `asky_floor_gap`.

**DESTINATION:** PAR **KILL** — not adopted, not pursued. The calibration curve
(`reports/sieve_t3_par.npz`, `reliability` / `glacier_reliability`) **PROMOTE** to the
criterion record as the measured cost of LEDGER-A1, and as the empirical case for the
`n_belief` sigma-point upgrade LEDGER specified. Bars-class (it is a statement about what
`QBAR = 0.9` actually buys): **posted and parked**.

---

## T2 — LISA-STYLE EVEN/ODD CROSS-VALIDATION

### VERDICT: the cross-validation **DOES NOT DISCRIMINATE** — KILL the split. The pre-registration is **REFUTED on its own terms**: the held-out refit fails the TRUTH-ANCHORED template 11–12 times in 16.
### But the quantity it wrapped does discriminate, and needs no split: **re-running frontier-v2's own data-support term on already-fed members**. PROMOTE that, at 2 likelihood calls per member.

**Code** `hpc_harbor/sieve/sv_t2_xval.py` · **banks**
`SIEVE_results/sieve_t2_xval_{r13p25_e07_s0_i5, r13p25_e07_s0_i0, r13p5_none_s3_i0}__*.npz`
· jobs `12835323/24/25`, wall 2663–3173 s each, **0 GPU-hr**.
**Gate G-T2a: bit-exact `0.000e+00` on every half built (even, odd, even-rebuilt) in all
three jobs** — the halves are the same venue, sparser not shorter.

**(B) NO-REFIT — the state scored directly on each half.** Counts are `PASS` = dlnL > 0 on
*both* halves, per fed slot (8 slots scored per cascade cell of 20–23 fed; 1 of 1 at
r13p5).

| cell | anchored | wandered |
|---|---|---|
| r13p25/e07/s0 **i5** (5 M-steps of motion, all certs wrong) | **8/8 PASS** | **2/8** |
| r13p25/e07/s0 **i0** | **8/8 PASS** | **1/8** |
| r13p5/none/s3 i0 (the true certs) | 1/1 PASS | **1/1 PASS** |

That separates cleanly, and in the right direction on both sides: the cascade cell's
wandered members go sharply **negative** (slot 261: −592.7 / −624.3; slot 260: −531.8 /
−534.7), while the same slots at the truth-anchored template are all strongly positive
(+426.9 / +413.5, +325.0 / +343.4). So the data really do contain those sources — the
negatives are template damage, not absent signal. And at the true-cert cell the wandered
state **passes**, which is the correct non-alarm: that certification is `on_true`.

**(A) HELD-OUT WITH REFIT — the actual cross-validation — fails.**

| cell | anchored | wandered |
|---|---|---|
| r13p25/e07/s0 i5 | **5/16** | 3/16 |
| r13p25/e07/s0 i0 | **4/16** | 2/16 |
| r13p5/none/s3 i0 | 2/2 | 2/2 |

The truth-anchored template — the best template that exists — fails its own held-out half
**11 of 16** times at i5 and **12 of 16** at i0. A test that rejects the truth cannot be
used to reject anything else. The mechanism is visible in the parameter-consistency column:
the two half-fits disagree by up to 6.4 M-step widths in `fgw` and 12.6 in `mc` (wandered
slot 256), and several widths come back non-finite — `mstep_quadratic` is a 2-axis
quadratic sweep on a half-sized, comb-multimodal surface with ~250 carried sources absent
from the model, so it overfits its own half and lands badly on the other. This is the same
multimodality LEDGER files as **C1**, met from a different direction.

**THE SPLIT IS DOING NO WORK — measured, not argued.** Across all 32 scored (cell, state,
slot) rows the even and odd halves **agree in sign on 32/32**. Every discrimination in
table (B) is already present in a single full-data evaluation; halving the data only halves
the SNR. So the LISA-style split is not the useful part of this test, and it costs three
extra venue builds (~30 min each) to obtain.

**What survives is cheaper than what was proposed.** The discriminating quantity is
`dlnL(k present) − dlnL(k absent)` at the **current** template with the other carried
members absent — which is *exactly* `run_cell` step (b)'s frontier-v2 data-support term
(`glacier_loop.py`, validated S4.20.1), except evaluated on **already-fed** members at the
current state instead of on candidates at feed time. The loop already computes this
quantity; it just never re-asks it after a member is fed. Cost to add: **2 likelihood calls
per fed member per iteration**, full data, no half-venue, no rebuild — genuinely free at
the scale of an iteration that already spends ~1400 s.

**DESTINATION:** **KILL** the even/odd split. **PROMOTE** the re-asked data-support term as
a **FORGE-B readout and a criterion-v2.3 candidate** — a fed member whose own-term support
has gone negative is a member the M-step has walked off its source. Pre-registration for
any adoption must be written before it is cut on the banked record, since the three cells
here were used to *find* it. Note the honest limit: it flags 13/16 wandered cascade members
but the cascade's certifications are on *pulsars*, not on the fed sources, so this is a
template-health monitor, **not** a replacement for D2 rigidity.

**SCOPE, declared:** 8 of 20–23 fed slots scored per cascade cell (slot order,
brightest-first), 1 of 1 at r13p5; `--max-fed 8`, banked as `n_fed_total` / `n_fed_scored`.
The `dlnL_full` reference column is absent by design (see the hazard note below).

### T2 method record — two bugs found before the numbers were believed

Stated because both would have produced a *confirmation* of the pre-registration:

1. **The contrast was degenerate** (submission 3, job `12834899`). `run_cell`'s frontier-v2
   writes the present-vs-absent pair as `setdiff1d(carried, [k])`, which is correct there
   because `k` is a *candidate* still in `carried`. T2 scores **already-fed** slots, for
   which `k ∉ carried`, so `th_on == th_off` and the statistic returned identically
   `0.000` — printing a clean "wandered FAIL, 0/2 PASS" that matched the pre-registration
   exactly and was pure artifact. The exact zeros are what gave it away. Fixed by building
   the absent set as `carried ∪ {k}`; the function now **asserts** the two templates differ
   rather than trusting it.
2. `bank_npz` was passed `stem` both positionally and as a keyword, so the science ran and
   the write failed.

**Operational finding, already banked (first submission `12834458`/`59`/`60`, all FAILED):
the CPU-lane map-count hazard is not confined to the E-step.** Three live
`build_b1_amortised` venues at ncw = 287 × 116 pulsars exhaust `vm.max_map_count`
(65530 on these nodes) and the third build dies inside XLA-CPU with
`INTERNAL: Failed to materialize symbols: {(<xla_jit_dylib_9>, {negate_power_fusion})}`
— raised from `MultiSourceDelay.__call__`, not from `estep_per_target`. **Two live venues
is fine**: full + even built and scored cleanly in every one of the three jobs; only the
odd build failed. A second failure (`12834779`) showed freeing the full venue was still not
enough: **compiling** its `logL` on top of a live `NoiseDrawer` (a 227 MB banked `L_gwb`
plus the refactored GWB block) exhausted the maps by itself — LLVM
`allocateMappedMemory failed: Cannot allocate memory`. Each JIT section is an mmap, so
compiles cost maps, not just RAM. The recorded hazard
(`cpu-lane map-count`, fixed in `fsky_stage0._install_evicting_ab` by evicting one
executable after use) therefore understates its own scope — it is a property of the
per-pulsar CW delay compilation generally, and it binds on **venue count**, not only on
the scoreboard. T2's fix is coarser than eviction and worked at the third attempt: **exactly one venue
live at a time** across four passes (full → even → odd → even rebuilt), and the full
venue's `logL` is **never compiled** — it is only asked for `inject_delay` and the noise
draw. Declared cost of that: the `dlnL_full` reference column is not produced. It was a
nice-to-have; the pre-registered test is on the halves and is unaffected.

Method, fixed before the runs: each half is a **stock venue build** with the discovery
pulsars subset to one TOA parity *after* the T-extension, so the halves are **sparser, not
shorter** — `ds.getspan` is unchanged, hence so are the span-scaled `rn_comp`/`gwb_comp`
GP counts, and the enterprise pulsars (which carry the injection and `theta_truth`) are
never subset. The realisation is built **once** on the full venue exactly as `run_cell`
does and then restricted, so both halves see the same banked noise draw rather than two
different ones. **Gate G-T2a**: the half `inject_delay(θ)` must equal the full
`inject_delay(θ)` restricted to that parity, **bit-exact**, on all 116 pulsars — no
downstream number is quotable otherwise. Cost measured by probe `12834315`:
`build_b1_amortised` = **640.5 s** per venue at ncw = 287, T = 30, so three builds ≈ 32 min.

---

## T6b — CRN: PAIRED vs INDEPENDENT

### VERDICT: **PAIRING WINS, LARGELY.** Variance ratio **7.7×** on the drain — the metric frozen-vs-live actually reads. **FLAG TO SUMMIT §2 CONVENTIONS** (bar was 1.5×).

**Code** `hpc_harbor/sieve/sv_t6b_crn.py` · **bank**
`SIEVE_results/sieve_t6b_crn__cn1486_noGPU.npz`, copied to the tracked
`reports/sieve_t6b_crn.npz` (the `*_results/` gitignore rule keeps campaign banks
ACCRE-local, and the per-run `yf`/`yl` arrays are what make the numbers below
re-derivable) · job `12834427`, wall 3839 s, **0 GPU-hr**.

6 seed pairs × 2 arms from one venue build (n_src = 16, T = 30, 3 iterations, circular
arm), each run through the frozen dial (per-slot freeze-after-first-fit) and the live
M-step on the **same** noise seed.

| metric | pairs | ρ(frozen, live) | Var_paired | Var_indep | **ratio** |
|---|---|---|---|---|---|
| `a_bg` — the drain | 6 | **0.957** | 1.699e−04 | 1.308e−03 | **7.70×** |
| `a_bg_sig` | 6 | 0.823 | 4.508e−07 | 1.394e−06 | 3.09× |
| `logL_end` | 6 | 0.986 | 5.707e+02 | 3.461e+04 | 60.6× |
| `fgw_hat` | 5 | 0.951 | 7.000e−04 | 1.314e−02 | 18.8× |
| `mc_hat` | 5 | 0.992 | 1.840e−03 | 1.353e−01 | 73.6× |
| `n_res` | 6 | 1.000 | **0** | 2.133 | ∞ |

Ratios are `Var_indep/Var_paired` with `Var_indep = Var(Y_f) + Var(Y_l)`, the
estimator-theory value; the crossed-pairs check on all 30 mismatched (i, j) agrees
closely where both are finite (`a_bg` 7.79× vs 7.70×, `logL_end` 62.6× vs 60.6×).
**Median over the five finite metrics: 18.8×. Minimum: 3.09×.** Every metric clears
1.5× by a wide margin, so **the SUMMIT §2 flag is raised**.

**What this buys.** Sharing the seed across arms cuts the variance of the arm difference
by 7.7× on the drain, i.e. **a paired design needs ~7.7× fewer realisations for the same
precision on `Δa_bg`** — directly, ~87 % of the GPU cost of any future frozen-vs-live
style comparison at equal precision. The paired `a_bg` differences are
`[−0.0260, +0.0005, −0.0134, −0.0151, 0.0000, −0.0312]` dex: a consistent, resolvable
frozen-minus-live offset that unpaired sampling at n = 6 would have buried, since the
across-seed spread of `a_bg` itself (~0.03 dex) is the same size as the effect.

**Two caveats, declared.**
1. `n_res` has **exactly zero** paired variance: the frozen and live arms resolved the
   same members on all 6 seeds (`[3,1,2,2,0,2]` both arms). At this toy scale the dial
   moves M-step *parameters* only and never changes the feed set, so the ratio is
   reported as ∞ but measures nothing. The five parameter/drain metrics carry the result.
2. `fgw_hat`/`mc_hat` are over **5** pairs, not 6: seed 4 fed nothing (`n_res = 0`), so
   those metrics are undefined for it in *both* arms. The quoted numbers come from a
   finite-pair-masked re-analysis of the banked per-run `yf`/`yl` arrays (they are in the
   bank, so it is reproducible); the driver's own analysis block has been corrected to
   mask non-finite pairs and print the surviving `n`, but was **not** re-run — that would
   have cost another 64 min of lane for numbers already banked.

**DESTINATION:** **PROMOTE** → SUMMIT §2 conventions: *pair the seeds across arms in every
two-arm comparison*. SIEVE-C reports whether the current drivers already do; PHOENIX's
frozen arm does pair by construction (it inherits the live arm's banked seeds), which this
number retrospectively justifies.

**Method note. The 0.5 GPU-hr was not spent, and the reason is a result in itself.** The toy's only
GPU-bound stage was `CertScoreboard.columns` → `spark3.estep_per_target`, which is the
documented XLA-CPU `vm.max_map_count` exhaustion hazard. But the certification columns do
**not** enter the arm difference under measurement — **the drain (`a_bg`) is what
frozen-vs-live actually reads** (capstone S4.15.1 item 2: "0.03–0.30 dex below baseline").
Dropping the E-step puts the whole test on the general CPU lane, which is why no claimed
GPU lane was entered. Metrics: `a_bg`, `a_bg_sig`, `logL_end`, and the M-step's own two
axes `fgw_hat` / `mc_hat` for the first-fed slot, plus `n_res`.

The freeze dial is **reimplemented** (behaviour-verbatim per-slot freeze-after-first-fit);
`hpc_harbor/frozen/frozen_mstep.py` is never imported, `FROZEN_results/` is never written,
and SIEVE uses its own seed base (`6_100_000`) so a toy can never be mistaken for a
PHOENIX cell. 6 seed pairs, `n_src = 16`, T = 30, 3 iterations; `Var_indep/Var_paired`
reported per metric with both the analytic (`Var(Y_f)+Var(Y_l)`) and crossed-pairs
estimates.

---
---

# SIEVE-A ADDENDUM — TRUE-CERT VALIDATION (2026-07-29)

*Matt's fairness requirement: the promoted tools must be tested where the loop WORKS, not
only where it fails. The GLACIER banks hold n = 2 true certifications — insufficient to
bound anything. V1/V2/V4 move the tests onto true-cert-rich substrate; V3 executes the
authorised E0 fix.*

## V1 — DOES THE E0 CANONICAL RECUT KEEP TRUE CERTIFICATIONS?

### VERDICT: **YES, but not for free — and my T7 claim of a strictly-dominating fix does not survive n = 1059.** True-cert survival **1044/1059 = 98.58%**; kill rate **1.42% [0.86%, 2.32%]** (Wilson 95%).

`hpc_harbor/sieve/sv_v1_truecert.py`, bank `reports/sieve_v1_truecert.npz`. CPU, banked
columns only — no likelihood re-evaluated, no campaign verdict moved.

**Substrate: 38 cells, three campaigns, 1059 true certifications** — GENERALISE arm A-SKY
(32 units = 4 cells × 8 skies, of which the survivor cell `e0.3 h-12.75 5+11` is 8),
CHORUS soft loops (4 cells), IGNITE-2 soft loops (2 cells). That is 530× the true-cert
count T7 had to work with.

### The correction to T7

T7 reported the E0 recut as "the only strictly-dominating tightening tested — 2→2 true
kept, 5 wrong removed." At n = 2 that was the only observable outcome; it was not
evidence of dominance. At n = 1059 the recut **does** cost true certifications, at
1.42% [0.86%, 2.32%]. The GLACIER result is fully consistent with that rate (expected
loss on 2 certs = 0.03), so nothing in T7's arithmetic was wrong — the claim of
dominance was over-read from a sample that could not support it. The fix is still worth
making; it is a trade, not a free lunch, and the trade is quantified below.

### The three cuts, scored on every realisation

| cut | what it is | certs | true | wrong | wrong rate |
|---|---|---:|---:|---:|---:|
| **A** canonical | what CHORUS/IGNITE-2/GENERALISE actually do, and what the patched GLACIER now does | 1141 | 1059 | 82 | 7.19% |
| **B** the defect | GLACIER's floor-scale mismatch, reconstructed on these banks | 1809 | 1664 | 145 | 8.02% |
| **C** E0 scale-matched | T7's recut: the banked floor applied on the scale it was measured on | 1124 | 1044 | 80 | 7.12% |

Per family, true-cert survival under C:

| family | cells | survival | net | kill rate [95%] |
|---|---:|---:|---:|---|
| A-SKY survivor (8 skies) | 8 | 369/382 = 96.6% | −13 | 3.40% [2.00%, 5.73%] |
| A-SKY all | 32 | 960/979 = 98.1% | −19 | 1.94% [1.25%, 3.01%] |
| CHORUS sloops | 4 | 57/62 = 91.9% | −5 | 8.06% [3.49%, 17.5%] |
| IGNITE-2 sloops | 2 | 27/18 = 150% | **+9** | 0% [0%, 17.6%] |

**C is not uniformly stricter than A.** It gains certifications in cells with a small
glacier-scale floor (IGNITE-2, +9) and loses them elsewhere. A and C agree on
**61329/61480 = 99.75%** of pulsar-rows and disagree in **both** directions. That
matters for V3: the scale-matched recut is a close stand-in for the canonical cut, not
an identity.

### Δ — the size of the defect, now measured properly

Δ = floor_canonical − floor_glacier over 38 cells: **median +6.97 nat**, mean +7.42,
range **[+1.96, +15.37]**, and **positive in all 38**. This supersedes the single
+9.07 nat figure T3 measured on one A-SKY unit's nulls.

The sign is not an empirical regularity — it is forced. `lnK = log(max(K,1)) ≥ 0`, so
the pulsar attaining the glacier maximum also satisfies `dlnL > lnK` and therefore lies
in the canonical mask, giving `off_glacier ≤ off_canonical` **on every draw**. For the
same reason both statistics vanish on exactly the same draws: the measured
`max|zf_canonical − zf_glacier|` across all 38 cells is **0.000e+00**, exactly. So
`zero_fraction`, and the degenerate-estimator branch that reads it, are untouched by
the correction. That identity is load-bearing in V3.

### Gates

- **G-V1a** — the canonical offender array recomputed from each cell's own null banks
  reproduces the campaign's *banked* offender array element-wise: **32/32 PASS** where a
  banked array exists (A-SKY `gen_ledger.offenders`), 6 n/a (CHORUS/IGNITE-2 bank a
  pooled array, not per-cell).
- **G-V1b** — the floor re-derived from that array with the campaign's own estimator
  reproduces the banked adopted floor: **38/38 PASS**.
- **G-V1d** — `off_glacier ≤ off_canonical` pointwise: no violations. This gate is what
  caught the first pass, which globbed `fnullA`/`fnullN` wrongly (the floor's own null
  set is `fnullN`; `fnullA`/`fnullL` are the *extra* draws behind the wider `fALL` floor)
  and printed an impossible Δ = −51 nat for IGNITE-2 as if it were a measurement.

**KILL / PROMOTE:** PROMOTE the canonical statistic (executed in V3). The 1.42% true-cert
cost is now on the record and is Matt's to accept or refuse.

## V3 — E0 FIX EXECUTED, AND EVERY GLACIER-LINEAGE VERDICT RE-CUT

### VERDICT: **PATCHED AND RE-CUT.** `glacier_loop._null_offenders` now computes the canonical offender. Banked verdicts: **18 → 13 certifications (E0 scale-matched, exact) or 18 → 11 (canonical, Δ-calibrated); wrong 16 → 11 or 9; both true certifications survive every variant.**

Code: `hpc_harbor/glacier/glacier_loop.py` — the defect-and-correction is written into the
function's own docstring, not just here, so the next reader of that code meets it there.
Re-cut: `hpc_harbor/sieve/sv_v3_recut.py`, bank `reports/sieve_v3_recut.npz`.

**Coverage: 730 banks** — 504 signal cells + 226 null/scrambled cells, spanning Stage-1
(`gl1_*`), the sky/array ladders (`gl2_r13p25`, `r13p5`, `r13p9`, `r13p9_w0p25`,
`r13p9_w0p5`), and the frozen arm (24 banks in `FROZEN_results/`, read-only — it has
landed). No campaign npz was rewritten; both columns ride in the SIEVE bank so any later
analysis must state which it used.

| cut | certs | true | wrong |
|---|---:|---:|---:|
| v2.2 **as banked** (the defect) | 18 | 2 | 16 |
| **E0 scale-matched** — exact from the banks | 13 | 2 | 11 |
| canonical, Δ = +6.97 (central) | 11 | 2 | 9 |
| canonical, Δ = +1.96 (most permissive) | 17 | 2 | 15 |
| canonical, Δ = +15.37 (most strict) | 5 | 2 | 3 |

Null/scrambled cells certify **0** under every cut, before and after.

### The context number that reframes the whole finding

**432 of 504 signal cells (86%) carry a DEGENERATE floor** (`floor = 0`,
`floor_est = emp_q95_degenerate`). There the criterion has already collapsed to the
trials bar `dlnL > lnK + 0.578`, and E0 is **inert**. The defect bites only on the 72
non-degenerate cells — which is exactly where all 18 certification events live. So the
correction is narrow in scope and total in effect on the events that exist.

### One thing this script got wrong first, and why it matters

The Δ-calibrated column initially added the pooled Δ to *every* floor, including the
degenerate ones — and that killed **both** banked true certifications (their floor is
0.000 and their `dlnL` is 2.35 and 1.60 nat, so a manufactured +6.97 bar erases them).
It is wrong because the zero-atom identity above forces Δ = 0 wherever the offender
sample is degenerate: an all-zero glacier sample **is** an all-zero canonical sample.
Fixed; both true certs now survive all five variants. Worth recording because the
plausible-looking "conservative" choice was the destructive one.

### The honest limit, and what is pre-registered for SUMMIT

`_null_offenders` banks only the scalar floor, never the per-pulsar null columns, so
**GLACIER's null draws are not recoverable from the banks and a bit-exact canonical
re-cut of banked verdicts requires re-running the nulls on a GPU.** That is
pre-registered, not performed. Until it is paid, the scale-matched column is the
defensible one and V1's 99.75% row-agreement is its error bar.

**G-V3a** criterion-v2.2 reproduced on **730/730** banks. **G-V3b** `off_gl ≤ off_can` on
**226/226** GLACIER null banks; on the 3 that are non-degenerate the GLACIER-native gap is
median **+5.79 nat** [+3.91, +6.53] — independent of V1's substrate and consistent with it.

Pre-registered for SUMMIT from iteration zero: (1) every GLACIER-lineage run started after
this commit uses the canonical statistic from iteration 0, and no cell may switch estimator
mid-trajectory; (2) banked verdicts are re-cut in report only; (3) the bit-exact re-cut
needs fresh nulls, costed above; (4) **D2 is untouched** — R1 (2F_coh ≥ 15.132) and R2
(Δ2F > 0) never read the floor, so the manufactured-set kills stand exactly as banked.

### G-V3c — POST-COMMIT VERIFICATION OF THE PATCH ITSELF (added after commit `657db13`)

The patched `_null_offenders` had been compile-checked but **never executed** — it only
runs inside a GLACIER cell, on GPU. A patch to a certification bar that has never run is
not a verified patch, so the new statistic was exercised directly against the canonical
function on 20 000 random 116-pulsar draws (including 2 858 draws carrying a deliberate
`inf` from a K = 1 pulsar):

| check | result |
|---|---|
| patched inline body vs `recut_surface.offender` | **0 mismatches / 20 000**, the only difference being the declared non-finite guard |
| zero-atom identity `(off_can == 0) ⇔ (off_gl == 0)` | **0 disagreements / 20 000** |
| `off_gl ≤ off_can` pointwise | **20 000 / 20 000** |

So the two identities V1 and V3 lean on — the zero atom, and the sign of Δ — are confirmed
on synthetic draws as well as on the 38 banked cells, and the patch computes the canonical
statistic exactly.

## V2 — FALSE-ALARM RATE OF THE T2 TEMPLATE-HEALTH MONITOR ON TRUE-ANCHORED MEMBERS

### VERDICT: **PROMOTION CRITERION MET.** False-alarm rate on truth-anchored **loud** members is **0/600 = 0.0000, 95% upper bound 0.64%** — under both forms of the contrast. PROMOTE to FORGE-B readout.

`hpc_harbor/sieve/sv_v2_monitor.py` + `sv_v2.sbatch`, job **12843873** (COMPLETED, 17m33s,
general CPU lane, 0 GPU-hr). Bank `reports/sieve_v2_monitor_c1.npz`.

**Substrate: GENERALISE arm A-SKY survivor cell, all 8 skies × 15 realisations × 16
members = 1920 member-scorings.** A-SKY is truth-anchored *by construction*, not by
inference: `generalise.py:376` builds the search template as
`theta_base = theta_src.copy(); theta_base[AI] = L0` — every source parameter at truth,
only the distances at the fiducial L0. That is the pure fringe problem with zero source
wander, which is exactly the condition under which a health monitor must not fire.

**CHORUS and IGNITE-2 sloops were deliberately excluded here** even though they hold the
true certifications V1 used. Their templates are *fitted*, so a member there is anchored
only to the extent the soft loop converged; scoring them would mix genuine wander into
the false-alarm count and inflate it. They can supply true certs; they cannot bound a
false-alarm rate.

### Two contrasts, because T2's promoted quantity sits between them

| contrast | definition | stratum | FAR | 95% CI | median dlnL |
|---|---|---|---:|---|---:|
| **incremental** (leave-one-out) | `logL(all 16) − logL(all but i)` | **LOUD (k=5)** | **0/600 = 0.0000** | [0.0000, 0.0064] | +6166 |
| | | FAINT (11) | 9/1320 = 0.0068 | [0.0036, 0.0129] | +19.3 |
| | | ALL | 9/1920 = 0.0047 | [0.0025, 0.0089] | +31.2 |
| **isolated** (member alone) | `logL(only i) − logL(none)` | **LOUD (k=5)** | **0/600 = 0.0000** | [0.0000, 0.0064] | +6321 |
| | | FAINT (11) | 528/1320 = 0.4000 | [0.3739, 0.4267] | +14.0 |
| | | ALL | 528/1920 = 0.2750 | [0.2555, 0.2954] | +45.3 |

Wilson intervals (Wald collapses to [0,0] at k = 0, which is the outcome the test exists
to detect). Per-sky spread is quoted in the bank because the 15 realisations inside a sky
share one geometry and are not independent: loud FAR is 0.00 in **every** sky.

**The stratification is not a convenience.** A member injected 1.5 dex below the loud set
carries no support for the data to find, and a monitor that says so is correct, not
false-alarming — calling that a false alarm would score the monitor on members the loop
would never feed. The `iso` FAINT rate of 40% is that statement, not a defect: a single
faint member against an empty template is genuinely unsupported about half the time.
**The promotion-relevant number is the loud one, and it is zero.**

### Gates

- **G-V2a** — refused to run unless the shared GENERALISE `L_gwb` bank already exists,
  so the job could not write into another campaign's bank directory. PASS.
- **G-V2c** — `logL(H_ABSENT = −30) − logL(−45) = +0.000e+00` exactly. −30 *is* the
  absence plateau on this venue. (spark3's −18 is not, per BASELINE at T = 30; this
  re-measures rather than inheriting the claim.)
- **G-V2d** — every present/absent template pair differs. This is the assertion added
  after T2's degenerate-contrast bug, and it is now load-bearing in two scripts.
- **Provenance** — `lgwb` fingerprint matched the banked realisation in **120/120**
  cases, so these are bit-reproductions of A-SKY's own data, not fresh draws from the
  same law. (8-physical-core pin; `fp=f92c9e36b460d6f5`.)

**KILL / PROMOTE:** **PROMOTE** the incremental (leave-one-out) contrast to the FORGE-B
readout — 2 likelihood calls per member per iteration, no split, no refit, FAR 0/600 on
anchored loud members. It remains a template-health monitor, **not** a D2 substitute:
D2 tests misspecification, this tests data support, and V1/V3 showed those fail in
different places.

## V4 — SIGMA-POINT E-STEP AT WORKING-REGIME BELIEF WIDTHS (SIEVE-C T5's working-regime arm)

### STATUS: **IN FLIGHT** — job **12845398**, banks at ~20:25 CT. This section is the PRE-REGISTRATION, written before any number was seen. The verdict goes below it.

`hpc_harbor/sieve/sv_v4_sigma.py` + `sv_v4.sbatch`. Bank
`SIEVE_results/sieve_v4_sigma_c1__*.npz` → `reports/`. CPU lane, 0 GPU-hr.

### The pre-registered question, and the decision rule

T3 measured `q_max > 0.9` to be uncalibrated: empirical `P(on_true | q_max > 0.9)` is
**0.786** on the truth-anchored A-SKY banks and **0.305** across GLACIER. LEDGER-A1 named
the mechanism — the E-step evaluates the fringe posterior at ONE source point, so `q_max`
is `P(fringe | θ̂)` charged as `P(fringe | data)`.

**The question is whether belief-averaging moves `q_max` toward calibration — not whether
it changes the certification count.** A variant that raised the count while leaving the
0.786 gap open has *failed* this test. Decision rule, fixed in advance:

- **ADOPT-FOR-THE-BELIEF-ARM** if `P(on_true | q_max > 0.9)` moves materially toward 0.9
  at a working-regime rung, with the Wilson interval excluding the incumbent value.
- **NO-GAIN** if the gap is unmoved (intervals overlapping the incumbent), regardless of
  what happens to the count.
- **HARMFUL** if `q_max` moves *away* from calibration, or if the reliability curve
  degrades at rungs where it was previously acceptable.

### The ladder is in fringes, and the conversion is exact

The pulsar-term phase goes as `f·L`, so a fractional frequency error is indistinguishable
from a fractional distance error: `δf/f ≡ δL/L`. One fringe is `δL = dL_a`. Hence

    σ(log10 fgw) = n_fr · dL_a / (L0_a · ln 10),   median over pulsars.

Rungs: **0** (the incumbent — A1's `sigma_points` collapses *structurally* to the point
rule, so this leg is an identity, not an approximation), **0.5** (sub-fringe), **2.0**
(few-fringe). `w0 = 1/3`, all weights strictly positive, and the average is taken in
**likelihood** space (logsumexp) — averaging log-likelihoods gives a geometric mean and
cannot widen a posterior at all, which is the failure mode this test would otherwise walk
straight into. The rule is imported from `ledger_a1_sigma_estep`, not reimplemented, so
A1's gates G-L1/G-L2 cover the arithmetic.

### Declared limitations, before the numbers

- **Belief on `log10_fgw` only.** `MSTEP_AXES` is `(I_FGW, I_MC)`, but only fgw has the
  closed-form fringe conversion above, so only there does a rung label mean what it says.
  `log10_mc` also moves the pulsar-term phase (the pulsar term is evaluated at a retarded
  time of order kyr, where the chirp has run) but its fringe-equivalent width would have
  to be calibrated numerically per cell; adding it with a guessed width would make every
  rung label a guess. **This is the conservative direction** — one believed axis widens
  the marginal less than two — so any improvement measured is a **lower bound**.
- **Reduced scope, and it is a cost wall, not a choice** (see below): ~7 realisations of
  120, `n_belief = 2`, 3 rungs. All 8 skies are still represented (realisation-major
  sweep) and the number actually scored is printed and banked, never silently truncated.

### The cost wall, measured across three submissions

| | |
|---|---|
| job 12844088 | **DIED** in the first sigma point — `Failed to materialize symbols` / LLVM `Cannot allocate memory` from `trackB_b1_core.estep:488`. 116 per-pulsar evaluators × ~4k mappings **cannot** fit under `vm.max_map_count = 65530` on any slate. |
| fix | `fsky_stage0._install_evicting_ab` ported in, plus `jax.clear_caches()` after the build: maps **30468 → 2036**. Eviction is mandatory here, not hygiene. |
| job 12844432 | Runs, but **368 s per E-step** (~3.1 s/pulsar), flat across calls. |
| hypothesis | `chorus.py:197` sets `jax_persistent_cache_min_compile_time_secs = 10`, above the ~3 s per-evaluator turnaround, so evicted evaluators are never cached. |
| job 12845138 | Threshold lowered to 0.2 s → **no change** (360.6 s vs 368.0 s; cache directory static at 3 entries / 376 MB). **Hypothesis refuted.** |
| conclusion | The eviction tax is **Python trace + lower** of the 116 per-pulsar terms, not XLA backend compilation. No executable cache can remove it. |

Consequence, declared: at ~330 s/E-step a 5-rung ladder buys ~3 realisations (~350 rows
per rung) inside the budget; a 3-rung ladder buys ~7 (~800 rows per rung). Since the
pre-registered readout is a **rate with a confidence interval**, rows are worth more than
rungs, so the ladder was trimmed to three and the declared sub-fringe-to-few-fringe band
is still spanned. The threshold knob is kept in the code and documented as *ineffective*,
so the next person does not re-derive it.

Measured cadence: venue 278 s; realisation 1 complete at 3494 s of the 25200 s budget.

### V4 VERDICT: **HARMFUL** — not merely no-gain. Sigma-point averaging moves `q_max` **away** from calibration, monotonically and by a large margin. Job 12845398 COMPLETED (7h42m), 8/120 realisations, all 8 skies, 928 rows per rung.

Bank `reports/sieve_v4_sigma_c1.npz`.

| belief (fringes) | n_eval | mean `q_max` | n(`q>0.9`) | **P(on_true \| q>0.9)** | 95% CI | gap vs 0.9 | on_true, all rows |
|---:|---:|---:|---:|---:|---|---:|---:|
| **0.00** (incumbent) | 1 | 0.9435 | 784 | **0.7997** | [0.7703, 0.8263] | −0.100 | 0.7381 |
| 0.50 (sub-fringe) | 5 | 0.9204 | 730 | **0.5041** | [0.4679, 0.5403] | −0.396 | 0.4461 |
| 2.00 (few-fringe) | 5 | 0.9190 | 710 | **0.2789** | [0.2471, 0.3130] | −0.621 | 0.2522 |

The intervals are disjoint and far apart. Against the decision rule fixed above, this is
the **HARMFUL** branch. The reliability curves say the same thing in more detail: at 2
fringes over-confidence grows in *every* bin (e.g. `q ∈ [0.90,0.95)`: mean q 0.932 vs
P(on_true) 0.095 — over-confident by 0.84, against 0.44 for the incumbent).

**Independent confirmation of T3 along the way.** The incumbent rung is the point E-step,
and it returns **0.7997 [0.7703, 0.8263]** — reproducing T3's 0.786 ± 0.011 on the same
substrate from a separate code path. The calibration gap is real and is now measured twice.

### The mechanism, and it is the interesting part

**`q_max` barely moves while accuracy collapses.** Mean `q_max` falls only 0.9435 → 0.9190
(−2.6%) while `on_true` falls 0.738 → 0.252 (−66%). So belief-averaging does **not widen**
the fringe posterior, which is what LEDGER-A1's diagnosis predicted it would do — it
**relocates** it. The rule therefore fails in the worst available way: it destroys the
fringe identification *without* lowering the confidence that is charged against the bar.
It is over-confident in a new way rather than a cured way.

Why: the fringe likelihood is **periodic with period one fringe**. The scaled UT places its
side points at ±√(m+λ) = **±1.73 SD** with **2/3 of the total weight off-centre** (m = 2,
w0 = 1/3, λ = 1; four points at 1/6 each). At the 0.5- and 2-fringe rungs those excursions
are ±0.87 and ±3.46 fringes — a full period or more — so two thirds of the likelihood
weight sits at essentially scrambled fringe phase. Averaging over a belief comparable to
the ambiguity spacing cannot resolve the ambiguity; it erases it.

### Scope — what this does and does not establish

- It does establish that on a **truth-anchored** template the rule is harmful at
  working-regime widths. Averaging around a *correct* centre can only dilute it, and
  A-SKY's centre is correct by construction (`generalise.py:376`).
- It does **not** establish what the rule does with the belief centred **off-truth**,
  where marginalisation could in principle cover the truth. That complementary arm is the
  natural follow-up and is **not** done here. The periodicity argument above suggests it
  will also scramble — a full-period average is phase-blind wherever it is centred — but
  that is an argument, not a measurement, and it is flagged as such.
- The finest non-zero rung measured is 0.5 fringes. The deep-sub-fringe regime (0.1) was
  in the 5-rung design that the cost wall removed, so the "nearly inert" end of the ladder
  is inferred from the structural collapse at 0 rather than measured at 0.1.
- One believed axis (fgw), 2 believed members, 8 of 120 realisations. Every one of those
  is a *conservative* limitation for a remedy — less averaging, not more — so none of them
  can explain away a harmful result.

**KILL / PROMOTE: KILL** the sigma-point E-step as a calibration remedy for `q_max`; do
not carry it into the belief arm on this evidence. **The `q_max` miscalibration itself
stays OPEN** — T3 found it, V4 confirms it, and the remedy LEDGER-A1 proposed for it does
not work at the widths that matter. That is the bars-class item with no fix in hand.
