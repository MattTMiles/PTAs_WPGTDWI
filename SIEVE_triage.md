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


---

<!-- ===== SIEVE-C sections (cronus), merged by union from the local copy ===== -->

# SIEVE_triage — the deterministic half (agent SIEVE-C, cronus)

**Session 2026-07-29. Venue: cronus, CPU only. `git pull --ff-only` clean (up-to-date at
`5461142`). Nothing committed, nothing pushed, no tracked file edited — new files staged
only; Matt commits.**

**HARD CONSTRAINT HONOURED.** cronus is RECOMPUTED-UNSAFE for noisy realisations.
**Every number in this report is zero-noise / Asimov.** No `rng` call touches a datum in
any of the three new compute modules; the two RNG calls that exist draw *nuisance
parameters* (T4's algebra gate) and *scrambled sky directions* (T4's manufacture ensemble),
neither of which is a noise realisation. The likelihood is the exact quadratic mismatch
`lnL = -0.5 Σ_p ||r_p(θ) − r_p(truth)||² / σ_p²`, so `lnL(truth) ≡ 0` and every ΔlnL quoted
is pure model-mismatch cost. The one test that genuinely needs fresh noise —
calibration-against-truth for T5 — is **deferred to SIEVE-A** and is not attempted here.

**BUDGET LAW — no cap approached.** T1 cap 2 GPU-hr, T5 cap 1 GPU-hr, T4 "CPU + light
4090": **0.00 GPU-hr spent on all four tests.** All three tests reduced to float64 numpy on
16×360 arrays (~18 min wall total) once the right structure was exploited — the pulsar-term
residual of pulsar *p* depends on `L_p` alone, so the distance block is **exactly separable**
and no 23-D optimiser is ever needed. This is reported as a *method* finding: the toy-scale
SIEVE questions are not GPU questions.

**THE TOY SPINE** (`hpc_harbor/sieve/sieve_toy.py`). 16 pulsars, top of the census by noise
rank (median TOA error ascending, restricted to pulsars holding both a canonical distance
prior and a tracked feather); **real** sky positions from `pos_t` in the tracked
`data_products/*.feather`; **real** distances and prior widths from
`CW_transition/best_distances.txt` (the canonical prior file). T = 18 yr, 20 epochs/yr.
*Toy choice, declared:* the measured per-pulsar spans are heterogeneous (6.5–20 yr); the toy
imposes a **uniform** T and cadence on every pulsar. That is a statement about the toy, not
about the real array's sampling.

Spine: J1909−3744, J1713+0747, J2234+0944, J0437−4715, J2241−5236, J1216−6410, J0125−2327,
J0613−0200, J1640+2224, J1614−2230, J2017+0603, J1811−2405, J1843−1113, J1853+1303,
J2010−1323, J0030+0451 (σ 7.4e−8 – 2.5e−7 s; L 0.156 – 2.499 kpc).

**A POINTER THAT DOES NOT RESOLVE (reported, not worked around).** The brief cites
`GLACIER_capstone.md §S4.24`. **No such section exists** — `GLACIER_capstone.md` (1154 lines)
runs to §S4.14 and the string `S4.2[0-9]` does not occur anywhere in the tracked `.md` set.
The "wander widths" the brief asks T5 to use are **`STORY.md §S4.2.4`** (+ §S4.2.3 for the
per-axis σ), which I have used and named explicitly. Likewise no file named as "the criterion
spec" exists; I took the criterion convention from `reports/GENERALISE.md`'s header
(criterion-v2.2: Gumbel valid only at zero-fraction ≤ 20 %, else empirical q95 ± bootstrap;
floors refit per cell). **None of T1/T4/T5/T6a is a floor-vs-count test, so criterion-v2.2 is
context here, not machinery.**

---

## T1 — WAVEFRONT CURVATURE

> ### VERDICT: **NEGLIGIBLE-EVERYWHERE** for the pipeline as it actually stands. The BEDROCK curvature rung is **KILLED**. Curvature is not a small effect that happens to be small — it is *exactly degenerate*, to Fresnel order, with the per-pulsar distance nuisance the programme already carries. It re-enters only if pulsar distances are ever pinned tighter than **7×10⁻⁴ of their canonical prior width**, which no pulsar in the census is.

**The mechanism, and it is an identity, not a coincidence.** Curvature and a pulsar-distance
error act on the observables through the *same single scalar* — the pulsar-term retardation
`dt_p`. Expanding the exact source-to-pulsar path
`dt_p = (|D ŝ − L_p p̂| + L_p − D)/c` to leading order in `L/D` gives the brief's Fresnel
phase `π L_⊥²/(λD)`, and that phase is reproduced **identically** by shifting the pulsar
distance by

```
    dL_eff,p  =  L_p² (1 + cos μ_p) / (2 D)
```

Because *any* observable that depends on the curvature only through `dt_p` cannot separate
it from a distance error, the degeneracy is exact — and it survives chirp, because a
frequency-evolving source still enters the pulsar term only via the retarded epoch.

**G0, the identity gate (measured).** Injecting the exact spherical front and scoring it
against a flat template with distances shifted by `dL_eff` costs, at SNR ≈ 97:

| D | SNR | max \|dL_eff\|/fringe | ΔlnL(sphere vs dL_eff-shifted flat) | ΔlnL(sphere vs Fresnel) |
|---|---|---|---|---|
| 16 Mpc | 96.7 | 0.0846 | **1.35e−07** | 1.35e−07 |
| 50 Mpc | 98.0 | 0.0271 | **1.42e−09** | 1.42e−09 |
| 200 Mpc | 98.4 | 0.0068 | **1.37e−11** | 1.37e−11 |

The two right-hand columns agree to every digit printed: the 1e−7 nat residual **is** the
Fresnel truncation error (the next order in `L/D`), not a failure of the identity. Against
`SNR²/2 ≈ 4.7e3` nat of signal, the identity holds to 3×10⁻¹¹ in relative terms.

**(a) + (b) — recovering curved data with FLAT templates.** ΔlnL cost of flatness and the
recovered source bias, in three distance regimes:

| D | regime | ΔlnL cost | Δsky | Δlog₁₀h | Δlog₁₀f | max\|ΔL\| |
|---|---|---|---|---|---|---|
| 16 Mpc | pinned | **1.9206** | 0.0010° | +6.1e−3 | +7.5e−6 | — |
| 16 Mpc | prior | **9.50e−3** | 0.0009° | +2.4e−3 | +8.6e−6 | 1.90 pc |
| 16 Mpc | free | **1.25e−3** | 0.0009° | +8.1e−4 | +8.6e−6 | 15.2 pc |
| 50 Mpc | pinned | 0.2020 | 0.0003° | +2.3e−3 | +2.4e−6 | — |
| 50 Mpc | prior | 7.65e−4 | 0.0003° | +6.7e−4 | +2.8e−6 | 0.013 pc |
| 50 Mpc | free | 8.18e−5 | 0.0003° | +1.9e−4 | +2.8e−6 | 15.2 pc |
| 200 Mpc | pinned | 0.01274 | 0.0001° | +6.0e−4 | +6.1e−7 | — |
| 200 Mpc | prior | 4.45e−5 | 0.0001° | +1.6e−4 | +6.9e−7 | 0.003 pc |
| 200 Mpc | free | 4.59e−6 | 0.0001° | +4.4e−5 | +6.9e−7 | 15.2 pc |

The pinned column scales as `D⁻²` exactly (1.9206/0.2020 = 9.51 vs (50/16)² = 9.77;
0.2020/0.01274 = 15.9 vs (200/50)² = 16.0), which is the Fresnel scaling and confirms the
fit is converged rather than optimiser-limited. **Even with distances pinned at truth,
curvature buys > 1 nat only inside D ≲ 22 Mpc at SNR ≈ 97** (cost = 1.92 × (16 Mpc/D)²,
∝ SNR²). Letting distances move under their *canonical* priors suppresses the 16 Mpc cost
**202×** to 9.5e−3 nat; unconstrained, **1540×** to 1.2e−3 nat. **Sky bias never exceeds
0.0010°** in any cell — three orders below the 0.188°/src break-even sky prior that the
R-posterior referendum identified, so curvature does not move that trade at all.

**Why it can never bite, per pulsar (16 Mpc, the worst case):**

| pulsar | L (kpc) | fringe (pc) | σ_L (pc) | dL_eff (pc) | cycles | **dL_eff/σ_L** | σ_L/fringe |
|---|---|---|---|---|---|---|---|
| J1909−3744 | 1.136 | 0.497 | 25.8 | 0.0018 | 0.0037 | 7.1e−05 | 51.9 |
| J1713+0747 | 1.053 | 0.554 | 60.9 | 0.0086 | 0.0154 | 1.4e−04 | 109.9 |
| J2234+0944 | 1.921 | 0.632 | 309.0 | 0.0535 | 0.0846 | 1.7e−04 | 488.6 |
| **J0437−4715** | 0.156 | 1.266 | **1.30** | 0.0009 | 0.0007 | **7.2e−04** | **1.03** |
| J0125−2327 | 1.348 | 0.880 | 91.2 | 0.0509 | 0.0579 | 5.6e−04 | 103.6 |
| J0030+0451 | 0.329 | 0.843 | 5.40 | 0.0029 | 0.0034 | 5.3e−04 | 6.41 |

**max dL_eff/σ_L over the whole spine = 7.2×10⁻⁴.** The curvature offset is at worst 0.07 %
of the tightest canonical distance prior in the census, and **every** prior spans ≥ 1 fringe
(min σ_L/fringe = 1.03, J0437 — the array's usual suspect, again the only pulsar anywhere
near a bar). Curvature is buried inside a nuisance the pipeline already marginalises.

**(c) — the GW-parallax check at 16 Mpc, curvature IN the template.** Profiling lnL over D:

| regime | peak | lnL(16 Mpc) − lnL(10⁵ Mpc) |
|---|---|---|
| pinned | **16 Mpc** | **+1.917 nat** |
| prior | 10⁵ Mpc | −1.05e−3 nat |
| free | 10³ Mpc | −5.62e−4 nat |

**GW parallax is real but unreachable.** With distances pinned at truth the profile peaks at
the true 16 Mpc and D is measured at 1.92 nat — the curvature phase *is* physically
recoverable, and this is the one place the effect exists. The moment distances are given
their canonical priors (let alone freedom) the profile is flat to ~10⁻³ nat and the truth is
not even favoured. The reason is structural: D enters through **one** parameter but its
signature is a per-pulsar pattern `∝ L_p²(1+cos μ_p)` living in a **16-dimensional** distance
nuisance space whose priors are 10³–10⁶× wider than the pattern. One parameter cannot be
separated from sixteen that already swamp it.

**KILL.** No promotion. The BEDROCK curvature rung should not be built.
*Scope lines that travel:* single non-evolving circular source, f = 10 nHz (λ = 0.972 pc),
SNR ≈ 97, 16-pulsar spine, local fits (no fringe hopping in the flat-recovery arm). Costs
scale as SNR² and as D⁻², so a source 10× louder at 16 Mpc reaches ~192 nat *pinned* — but
the pinned column is not the pipeline's regime, and the prior/free columns scale identically,
leaving the 202×/1540× suppression intact.

---

## T4 — CLOSURE OBSERVABLES

> ### VERDICT: **CLOSURE KILLS THE MANUFACTURE CLASS — but it is not the only thing that does.** The C3 double difference cancels the per-source phase and the per-pulsar achromatic clock exactly (to machine precision), and it drops the wrong-sky manufacture rate from **2/16 to 0/16**. However an explicitly modelled per-pulsar phase nuisance achieves the identical 0/16 at a statistically indistinguishable cost. Closure's value is **robustness and model-freedom, not new information** — and it has a real hole: it does **not** cancel a per-pulsar *chromatic* clock.

**Toy:** 6 pulsars (the spine's top 6) × 3 sources at f = {10.0, 15.1, 21.9} nHz, real
positions/distances, zero noise.

**PART 1 — the algebra, verified numerically.** Observable
`φ_pk = g_pk(sky, L_p, f_k) + Φ_k + c_p [+ 2π f_k τ_p]`, where `Φ_k` is the per-source phase
offset, `c_p` a per-pulsar **achromatic** phase offset, and `τ_p` a per-pulsar **time**
offset (a genuine clock error → chromatic phase `2π f_k τ_p`). Max deviation of each closure
class from its nuisance-free value:

| nuisance injected | C1 source-pair | C2 pulsar-pair | **C3 closure** | C4 chromatic |
|---|---|---|---|---|
| Φ_k only | 4.222e+00 | 3.64e−12 | **3.64e−12** | 3.66e−04 |
| c_p only | 3.64e−12 | 4.116e+00 | **7.28e−12** | 1.397e+08 |
| Φ_k + c_p | 4.222e+00 | 4.116e+00 | **3.64e−12** | 1.397e+08 |
| τ_p only | 3.117e−01 | 1.314e+00 | **4.048e−01** | 4.88e−04 |
| Φ_k + c_p + τ_p | 4.534e+00 | 4.204e+00 | **4.048e−01** | 1.397e+08 |

The geometric phases here are O(10⁴) rad and the C4 column is in units of seconds (O(10⁸)),
so **every "zero" above is machine precision in relative terms** (3.6e−12/1e4 ≈ 4e−16;
3.7e−4/1.4e8 ≈ 3e−12). Reading the kernels off the table:

- **C1** (source-pair, fixed p) kills `c_p`, keeps `Φ_k`.
- **C2** (pulsar-pair, fixed k) kills `Φ_k`, keeps `c_p`.
- **C3** (double difference) kills **both** — as advertised.
- **C3 does NOT kill a chromatic clock**: `τ_p` leaves a 0.405 rad residual. This is a
  *finding*, not a footnote — "per-pulsar clock-like term" is two different nuisances and
  closure only kills one of them. A real timing-clock error is chromatic.
- **C4** (the same double difference formed on `φ/f` rather than `φ`) kills `τ_p` **and**
  `Φ_k`, but **not** `c_p`. The two closure flavours have *complementary*, not nested,
  kernels; neither alone covers both clock classes.

**PART 2 — the money test (the STOP-#2 wrong-sky manufacture mode).** Data at truth; the
recovery template's source-0 sky scrambled isotropically (16 deterministic draws). Because
the observables are phases, everything runs in the **wrapped** (mod 2π) domain with pulsar
distances free over their **full canonical ±3σ prior**, which spans **8.8 to 2087 fringes**
per pulsar — i.e. **fringe hopping is allowed**, as it is in the real problem. This matters:
a first pass restricted to ±1 fringe rejected every scramble trivially and measured nothing.
Three statistics, identical freedom:

- `raw_blind` — wrapped residual with only `Φ_k` profiled out (**the pipeline's actual
  freedom — it carries no per-pulsar phase term**)
- `raw_nuis` — both `Φ_k` and `c_p` profiled out (the fair comparison)
- `closure` — C3, immune to both by construction

| arm | statistic | truth | scram min | scram median | **manufacture (scram beats truth)** |
|---|---|---|---|---|---|
| clean (no clock) | raw_blind | 0.0000 | 1.9474 | 18.288 | 0/16 |
| clean | raw_nuis | 0.0000 | 0.0366 | 0.1920 | 0/16 |
| clean | closure | 0.0000 | 0.0331 | 1.2458 | 0/16 |
| **clock injected** | **raw_blind** | **6.0990** | **4.9558** | 14.938 | **2/16** |
| clock injected | raw_nuis | 1.50e−13 | 0.0366 | 0.1920 | **0/16** |
| clock injected | **closure** | **0.0000** | 0.0331 | 1.2458 | **0/16** |

**The manufacture mode fires exactly when a per-pulsar systematic is present and unmodelled,
and only then.** In the clean arm all three statistics put truth on top and no scrambled sky
comes near — closure buys nothing. Inject an unmodelled achromatic per-pulsar clock and the
c_p-blind statistic is displaced from truth (χ² = 6.099 instead of 0) by enough that **2 of
16 scrambled skies fit the data better than the true sky** — a 12.5 % wrong-sky manufacture
rate, manufactured *not* by the sky freedom but by a systematic the model cannot represent.
Closure returns truth to exactly 0.0 and rejects all 16.

**But `raw_nuis` does the same job (0/16, truth 1.5e−13).** And the two agree numerically
where it counts: identical scram minima to four digits (0.0366 vs 0.0331 are different
statistics on different normalisations; the *ordering and the rejection* match). That is not
a coincidence — the closure is the *nuisance-profiled* raw statistic. Both remove 8 of the
18 phase dof and leave 4 after the 6 distances, so they are information-equivalent by
construction. **Closure is therefore a reparametrisation, not a new lever.**

**PROMOTE — narrowly, as an audit instrument, not as a detection statistic.**
→ *Destination:* a **C3+C4 closure-consistency diagnostic** run beside the existing
certification layer (never in place of it), whose job is to answer "is there an unmodelled
per-pulsar phase systematic?" — the one question that provably manufactures wrong-sky
certifications. Its value over `raw_nuis` is that it needs no nuisance model and costs 10
numbers instead of a 9-parameter profile; its cost is that C3 alone is blind to chromatic
clocks, so **both** flavours must be carried.
*Scope lines:* 6×3 zero-noise, phase-domain, one source's sky scrambled (not all three);
n = 16 scrambles, so 2/16 carries a ±1-event counting error and is a *demonstration that the
mode exists*, not a rate. **The rate itself, and its behaviour under noise, is SIEVE-A's.**

---

## T5 — SIGMA-POINT E-STEP, GATE HALF

> ### VERDICT: **GATE PASSES.** The sigma-point E-step reproduces the point E-step bit-comparably as the belief width → 0 (max TV = **1.11e−16** at exactly zero width, **7.92e−11** at 1e−8, converging as width² once the perturbation is sub-radian). At realistic-class source uncertainty the fringe posterior shifts by TV median **0.069–0.138** and the MAP fringe flips for **15 of 16** pulsars — but the underlying posterior is **near-uniform** (p(true fringe) 0.081 vs 0.059 uniform), so the flips are cheap and the *decision-relevant* shift is small. **PROMOTE the machinery, KILL the expectation that it rescues anything.**

**Why the toy had to chirp.** §S4.2.4 says 95–100 % of the pulsar-term phase wander is
`log10_mc`; a non-evolving toy has no `mc` at all. T5 therefore runs a GR-driven circular
binary with the exact closed form `f(t) = [f₀^(−8/3) − k t]^(−3/8)`, evaluated at the Earth
time and at the retarded epoch. Retardations are **171–16 010 yr**, so the pulsar term is a
genuinely earlier, lower-frequency epoch of the same binary. SNR = 128.5, 16 pulsars.

**Cross-check against §S4.2.4 (independent reproduction).** Per-dex gradients are
`|∂φ/∂log₁₀mc| = 0.94–6209 rad` and `|∂φ/∂log₁₀f| = 778–58 920 rad` — i.e. **per dex the
wander is dominated by f, not mc**. The banked "95–100 % is mc" is a *Fisher-weighted*
statement (gradient × that axis's posterior σ), and once weighted by σ(log₁₀Mc) = 0.501 dex
against the Earth-term-pinned σ(log₁₀f) = 0.144/389 dex, the mc share of the wander variance
is **72.8 – 100.0 %**, bracketing the banked 95–100 %. The low end (72.8 %) is a
different-pulsar-subset / fixed-Mc = 10⁹ M⊙ effect, not a contradiction. **The banked
statement is confirmed and its normalisation is now explicit: it is per-σ, not per-dex.**

**Sub-fringe artefact guard — the banked hazard is large here.** Scoring each fringe cell at
its **profiled** sub-fringe offset (ladder anchored on the *prior mean*, never on the
injected truth) versus the centre-pinned grid changes the fringe posterior by
**TV = 0.144 – 0.941, median 0.718.** That is a near-total rewrite of the posterior, and it
is a direct, independent confirmation of the banked sub-fringe-offset artefact at toy scale.
**Every T5 number here uses the guarded version.**

**GATE — width → 0:**

| belief-width scale | max TV vs point E-step | max \|Δw\| |
|---|---|---|
| 1e−0 | 1.667e−01 | 1.666e−01 |
| 1e−2 | 8.074e−02 | 4.039e−02 |
| 1e−4 | 3.784e−02 | 1.923e−02 |
| 1e−6 | 7.925e−07 | 3.982e−07 |
| 1e−8 | 7.925e−11 | 4.054e−11 |
| **0** | **1.110e−16** | **5.551e−17** |

**PASS.** At exactly zero width the sigma points collapse and the two E-steps agree to
machine epsilon. The 1e−6 → 1e−8 leg falls as width² — correct unscented behaviour. The
non-quadratic middle (1e−2, 1e−4) is *physical, not a gate failure*: at scale 1e−2 the
effective σ(log₁₀mc) = 5e−3 dex still drives up to 31 rad ≈ 5 fringes of pulsar-term wander,
so the perturbation is not small until scale ≲ 1e−6. Quoting this leg as a failure would be
a misreading; the gate is the width → 0 limit, and it passes.

**READOUT — the shift at §S4.2.3/§S4.2.4 widths:**

| σ(log₁₀Mc) | TV median (min–max) | MAP-fringe flips | mean p(true fringe) point → σ-point | wander at this width |
|---|---|---|---|---|
| 0.501 dex | 0.0689 (0.0021–0.1667) | 15/16 | 0.0814 → 0.0751 | 0.47–3111 rad = 0.075–495 fringes |
| 0.864 dex | 0.1352 (0.0086–0.2122) | 15/16 | 0.0814 → 0.0742 | 0.81–5360 rad = 0.13–854 fringes |
| 0.866 dex | 0.1376 (0.0066–0.2123) | 15/16 | 0.0814 → 0.0734 | 0.81–5380 rad = 0.13–856 fringes |

**Read this carefully, because the two headline numbers point opposite ways.** The MAP
fringe flips for 15 of 16 pulsars at every realistic width — which sounds decisive — but the
guarded fringe posterior is **already almost flat**: p(true fringe) = 0.0814 against a
uniform 1/17 = 0.0588. The argmax of a near-uniform distribution is not a meaningful
quantity, and the flips cost almost nothing. The honest metric is TV (0.069–0.138) and the
change in p(true fringe) (−8 % relative at the widest width), both **modest**. The
mechanism: at σ(log₁₀Mc) = 0.501 dex the pulsar-term wander is **up to 495 fringes**, so
averaging over the belief does not blur the posterior *between* neighbouring fringes — it
scrambles the ladder wholesale, and a wholesale scramble of an already-flat ladder returns
something flat. **Sigma-point averaging is the correct E-step, and it is correct for the
reason that makes it nearly pointless here: the fringe posterior carries so little
information at realistic source uncertainty that there is very little for the averaging to
change.** This is fully consistent with the banked D3 result that a single CW resolves
fringes for zero pulsars.

**PROMOTE (machinery), KILL (the hope).**
→ *Destination:* land the 5-point unscented E-step in the M-step/E-step driver as the
**default** — it is gate-clean, costs 5× one E-step, and removes an approximation that is
indefensible at σ(log₁₀Mc) ≈ 0.5–0.87 dex. Do **not** budget for it to move certified
counts; the readout says it will not.
**HANDOFF TO SIEVE-A (explicit):** calibration-against-truth — does the sigma-point E-step
produce *better-calibrated* fringe posteriors than the point E-step over noisy realisations?
— is **not attempted here** and cannot be, on banked data or otherwise, without noise draws
this box cannot make safely. SIEVE-A owns it. What SIEVE-A should inherit: the guarded
fringe ladder (prior-mean-anchored, sub-fringe-profiled), the 5-point unscented weights, and
the finding that TV is small *because the posterior is flat* — so the calibration test must
be powered against a near-uniform null, not against a peaked one.

---

## T6a — CRN AUDIT (code read, zero compute)

> ### VERDICT: **THE ARMS SHARE SEEDS. Both contrasts are fully PAIRED (common random numbers) by construction — census seed, noise seed, distance-truth seed, and scramble seed are all arm-independent.** This is deliberate and documented, not accidental. The paired-vs-independent variance demo is therefore a measurement of *how much* the pairing buys, not a check of *whether* it exists.

The design is stated in the module docstring itself —
`hpc_harbor/glacier/glacier_loop.py:10-12`: *"TWO IGNITERS (Matt's budget decision, Option
1), **PAIRED on census + noise seeds** (the CHORUS C6 discipline — the igniter contrast is
exact)"* — and the seed convention at `glacier_loop.py:69-70`: *"census 9_0xx_xxx (9_000_000
+ sky); noise 9_5xx_xxx (9_500_000 + 1000*sky + real); dist = noise + 10_000_000; scrambles
9_8xx_xxx."* **No term in any of those formulas is a function of the arm.**

**Line cites — every seed, and what it does not depend on:**

| seed | line | expression | arm-dependent? |
|---|---|---|---|
| census (Stage-1) | `glacier_loop.py:394` | `draw_population(SEED_POP_BASE + sky, ...)` | **no** |
| noise | `glacier_loop.py:421` | `noise_seed = SEED_NOISE_BASE + 1000*sky + real` | **no** |
| distance truth | `glacier_loop.py:422` | `dist_seed = noise_seed + 10_000_000` | **no** |
| scramble | `glacier_loop.py:443` | `default_rng(SEED_SCRAM_BASE + 1000*sky + real)` | **no** |
| bases | `glacier_loop.py:123-124` | `SEED_NOISE_BASE = 9_500_000`, `SEED_SCRAM_BASE = 9_800_000` | — |

**e07-vs-none.** The arm enters at `glacier_loop.py:404`,
`embed_igniter(pop, E_IGNITER[arm], VENUE_SPAN_S[T])` — it changes the *source slot list*
(e07 appends 31 harmonic slots; `none` is the same census all-circular,
`glacier_loop.py:120`), never a seed. The noise vector is therefore **bit-identical across
arms**, and this is provable rather than merely plausible: `NoiseDrawer.draw`
(`CW_transition/trackB_b1_core.py:592-613`) consumes its RNG stream as `c_gwb` (skipped
entirely when `"gwb"` is not in `components` — and GLACIER passes
`components=("white","rn")`, `glacier_loop.py:437`), then per-pulsar white, then per-pulsar
red. **The stream depends only on the pulsar set and TOA layout, never on the source list.**
Same pulsars + same T + same seed ⇒ same draw, whatever the arm did to the census.

**frozen-vs-live.** `mstep_policy` is a *venue axis of the output stem only* —
`glacier_loop.py:387-388` appends `_m<policy>` to the stem, with the explicit comment that
`"live"` appends nothing so *"incumbent stems"* stay bit-for-bit. The belief machinery is
constructed at `glacier_loop.py:466-471`, i.e. **after** the noise draw at line 437, so it
cannot perturb the RNG stream even incidentally. The seeds above are untouched by it.

> **CONCURRENCY, recorded (the §g2.4 hazard, live again).** `hpc_harbor/glacier/glacier_loop.py`
> **changed while this audit was running** — a concurrent session is staging FORGE-B belief
> work on this tree (`+177/−11` on `glacier_loop.py`, plus `CW_transition/trackB_b1_belief.py`,
> `reports/FORGE_B_belief.md`, `CW_transition/logs/forge_b_belief_gates.log`, and an `M` on
> `CW_transition/trackB_b1_softsource.py`, all of which appeared **in the index** and are
> **not mine**). I re-verified every T6a cite against the file as it now stands: lines
> 120/123/124/387/394/404/421/422/437/443 are **unchanged and still say what is quoted above**;
> the belief-book cite moved 456 → 466 and is corrected here. **The CRN verdict is unaffected** —
> no seed expression was touched. I did not stage, revert, or edit any of the concurrent
> session's files, and I did not race it.

**Two facts that travel beside the pairing, both from the same read:**

1. **Signal and scrambled-null cells at the same `(sky, real)` share the noise draw too**, since
   `noise_seed` does not depend on `scrambled` — and the scramble is applied to the *recovery
   template only*, never the data (`glacier_loop.py:440`, comment; `442-445`). The null is
   therefore *the same data with a scrambled recovery sky*, which is the strongest possible
   form of the control. Worth stating explicitly because it means the null and signal arms
   are paired as well.
2. **The pairing is exact within a `(T, wscale)` cell and NOT across T.** `white_scale`
   preserves the stream by construction (`trackB_b1_core.py:593-598`: the same normals are
   drawn at any scale, and `wscale = 1.0` is bit-identical). But across different `T` the
   same seed produces a *different-length* draw and hence a different realisation — so any
   T = 30 vs T = 40 comparison is **unpaired**, unlike every arm comparison. If SIEVE-A's
   variance demo touches the time ladder, that leg needs different bookkeeping.

**PROMOTE (as a measurement, already de-risked).**
→ *Destination:* SIEVE-A's paired-vs-independent variance demo can proceed directly on the
banked GLACIER arm pairs with **no plumbing change** — the CRN is already there and exact.
The demo's deliverable is the variance-reduction factor on ΔN and on the drain, per cell;
its scope must exclude the T axis per finding (2).

---

## SPEND

**0.00 GPU-hr** against caps of 2 (T1) + 1 (T5) + "light 4090" (T4). ~18 min of single-core
CPU wall across all three compute modules, plus one code read. No cap approached; no test
was trimmed, deferred, or reduced in scope for budget. The only deferral in this report is
the T5 noisy-calibration handoff, which is a *venue-safety* deferral (cronus is
RECOMPUTED-UNSAFE), not a budget one.

## ADD-LIST (staged with `git add`; **Matt commits — I have committed and pushed nothing**)

```
SIEVE_triage.md                              new   this report (T1/T4/T5/T6a)
hpc_harbor/sieve/sieve_toy.py                new   the zero-noise toy spine + CW model
hpc_harbor/sieve/sieve_t1_curvature.py       new   T1 wavefront curvature
hpc_harbor/sieve/sieve_t4_closure.py         new   T4 closure observables + money test
hpc_harbor/sieve/sieve_t5_sigmapoint.py      new   T5 sigma-point E-step, gate half
reports/sieve_t1_curvature.npz               new   T1 bank (gate, rows, D-profiles)
reports/sieve_t4_closure.npz                 new   T4 bank (algebra table, money ensemble)
reports/sieve_t5_sigmapoint.npz              new   T5 bank (gate, readout, wander shares)
```

No tracked file was edited. `hpc_harbor/glacier/glacier_loop.py` was **read only** for T6a.

**DO NOT DOUBLE-STAGE — the index is shared.** `git add` was run on the eight paths above and
**only** those. The index also currently holds, authored by a **concurrent session, not me**:
`hpc_harbor/glacier/glacier_loop.py` (M, +177/−11), `CW_transition/trackB_b1_belief.py`,
`CW_transition/trackB_b1_softsource.py` (M), `CW_transition/logs/forge_b_belief_gates.log`,
`reports/FORGE_B_belief.md`, `LOTTERY2.md`, `hpc_harbor/lottery/lottery2_realistic_pop.py`,
`reports/lottery2_realistic_pop.npz`, `reports/OBSERVERS_TABLE.md` (M). **Reconcile those with
their author before committing**; the remaining dirty working tree (CW_transition/*.npz,
SAMPLER_dev.md, notebooks, logs) is likewise not mine and not staged by me.

## STILL OWED / named (not dropped)

- **The brief's `§S4.24` pointer does not resolve** (§ above). I proceeded on `STORY.md
  §S4.2.4 + §S4.2.3` and named it everywhere. If that substitution is wrong, T5's readout
  widths change and the readout must be re-cut; the gate is unaffected.
- **T4's 2/16 is a demonstration, not a rate.** The wrong-sky manufacture rate under an
  unmodelled per-pulsar clock, and its noise behaviour, is SIEVE-A's.
- **T4's chromatic-clock hole (C3 leaves 0.405 rad on `τ_p`)** is unclosed by C3 alone. If
  the closure diagnostic is built, C4 must be built with it.
- **T5 calibration-against-truth** — SIEVE-A, on banked data. Power it against a
  near-uniform null.
- **T1's `free`-regime `max|ΔL| = 15.2 pc`** is a genuine walk along the exactly-degenerate
  direction inside the ±1-fringe local bracket, not a convergence failure — but it is
  bracket-limited, so the *unconstrained* suppression factor (1540×) is a lower bound on the
  true one.

---

# CLOSEOUT DELIVERABLES

## C1 — T1 paper text → *the draft's systematics discussion*

Drop-in paragraph. Numbers are the toy's (16 pulsars, 18 yr, SNR ≈ 97, `f = 10⁻⁸ Hz`);
re-cut against the real array before submission — the *scalings* and the identity are
array-independent, the absolute ΔlnL is not.

> **Wavefront curvature.** Continuous-wave searches model the source as a plane wave, but a
> source at finite distance `D` emits a spherical front, and the pulsar term samples it a
> transverse baseline `L_⊥ = L_p sin μ_p` off the Earth–source axis. Expanding the exact path
> `(|D ŝ − L_p p̂| + L_p − D)/c` to first order in `L_p/D` shows that the correction enters
> the likelihood *only* through the pulsar-term retardation `Δt_p`, in exactly the place the
> unknown pulsar distance already enters. Curvature is therefore not an independent
> systematic but an **exact reparametrisation of the per-pulsar distance nuisance**, with
> effective offset
>
> `dL_eff,p = L_p² (1 + cos μ_p) / (2 D)`,
>
> equivalently the Fresnel phase `π L_⊥²/(λ D)`. We verify the identity numerically: injecting
> a spherical front and recovering with a plane-wave template whose distances are shifted by
> `dL_eff,p` leaves a residual of 1.3 × 10⁻⁷ nat at `D` = 16 Mpc, falling as `D⁻²` to
> 1.4 × 10⁻¹¹ nat at 200 Mpc — identical, to all printed digits, to the truncation error of
> the first-order expansion itself, i.e. the shifted plane-wave template is the spherical
> template to the order we work at. Two consequences follow. First, the effect is
> **correctable at zero cost**: any pipeline that already fits pulsar distances absorbs it
> automatically. With distances pinned at their catalogue values the flatness penalty is
> 1.92 nat at 16 Mpc, but simply opening the distances under their standard Gaussian priors
> suppresses it by 202× to 9.5 × 10⁻³ nat, and the induced sky bias is ≤ 0.001° in every
> configuration tested. Second, the effect is **negligible for all but the nearest plausible
> source**: `dL_eff,p ≤ 0.053 pc` at 16 Mpc, which is 7 × 10⁻⁴ of the smallest distance
> uncertainty in the array and 0.08 of a fringe, and the uncorrected penalty falls below
> 1 nat beyond `D` ≈ 22 Mpc at this SNR (and further out still at lower SNR, since the penalty
> scales as `h² D⁻²`). The same identity forecloses the optimistic reading: because curvature
> acts *only* along the already-free distance direction, it carries no independent distance
> information. Profiling the likelihood over `D` with curvature in the template recovers the
> true 16 Mpc at +1.92 nat when distances are pinned, but that peak collapses to ≲ 3 × 10⁻³
> nat — flat, with the truth not even preferred — the moment the distances are allowed their
> priors. **There is no GW parallax to be had from a single continuous wave at these
> baselines.** We therefore adopt plane-wave templates throughout and quote the curvature
> systematic as bounded by the table below.

**One-line correction table** (16-pulsar toy, `f = 10⁻⁸ Hz`, SNR ≈ 97; `Δ` = flat-template
penalty, `resid` = what survives after applying `dL_eff,p`):

| `D` (Mpc) | max `dL_eff` (pc) | median | max `dL_eff/σ_L` | max cycles | ΔlnL pinned | ΔlnL prior-free | resid after correction |
|---|---|---|---|---|---|---|---|
| 16 | 0.0535 | 0.0077 | 7.2e−4 | 0.085 | 1.92 | 9.50e−3 | 1.35e−7 |
| 50 | 0.0171 | 0.0025 | 2.3e−4 | 0.027 | 2.02e−1 | 7.65e−4 | 1.42e−9 |
| 200 | 0.0043 | 0.0006 | 5.8e−5 | 0.0068 | 1.27e−2 | 4.45e−5 | 1.37e−11 |

Every column scales as `D⁻¹` (offsets) or `D⁻²` (ΔlnL) to 2–3 %, which is the convergence
check as well as the physics. Source: `reports/sieve_t1_curvature.npz`,
`hpc_harbor/sieve/sieve_t1_curvature.py`.

## C2 — T5 erratum note → **flagged for the ACCRE fold; capstone NOT edited here**

> **Erratum / clarification, §S4.2.4 — "95–100 % of the pulsar-term phase wander is
> `log10_mc`".** The statement is correct but **normalisation-dependent, and the normalisation
> is not stated in the text.** It is a *per-σ* (Fisher-weighted) decomposition: each axis's
> phase gradient multiplied by that axis's posterior width. Under that weighting the SIEVE-C
> toy independently reproduces an `mc` share of the wander variance of **72.8 – 100.0 %**,
> bracketing the banked 95–100 % (the low tail is a different-pulsar-subset and fixed
> `Mc = 10⁹ M⊙` effect, not a disagreement). Under the *per-dex* normalisation the ranking
> **inverts**: bare gradients are `|∂φ/∂log₁₀mc| = 0.94 – 6.21e3` rad against
> `|∂φ/∂log₁₀f| = 778 – 5.89e4` rad, giving an `mc` share of only **1.5e−6 – 1.1e−2**, i.e.
> per dex the wander is overwhelmingly **`f`-dominated**. Both statements are true of the same
> system; they differ only in whether the axes are compared per dex or per posterior sigma.
> The `mc` dominance in §S4.2.4 is bought entirely by `σ(log₁₀Mc) ≈ 0.5 dex` against an
> Earth-term-pinned `σ(log₁₀f) = 0.144/389 dex ≈ 3.7e−4 dex`. **Recommended fix:** add "per
> posterior σ" to the §S4.2.4 sentence, and note that the conclusion does not survive a
> reader who assumes per-dex. Anything downstream that inherits the "`mc` dominates" claim
> while using an *unpinned* or wider `f` posterior must re-derive the share — the crossover
> is at `σ(log₁₀f)/σ(log₁₀Mc) ≈ 10⁻²`, not far below the pinned value.
>
> Source: `reports/sieve_t5_sigmapoint.npz` keys `mc_share` (per-σ) and `mc_share_perdex`,
> `hpc_harbor/sieve/sieve_t5_sigmapoint.py`. Zero-noise toy; the *shares* are geometric and
> noise-free, so this does not defer to SIEVE-A.

## C3 — STAGE HYGIENE: files authored by SIEVE-C **only**

Per-author commits — the eight paths below are mine and complete. Everything else in the
index or the working tree belongs to FORGE-B or another concurrent session; reconcile with
their add-list separately.

```
SIEVE_triage.md
hpc_harbor/sieve/sieve_toy.py
hpc_harbor/sieve/sieve_t1_curvature.py
hpc_harbor/sieve/sieve_t4_closure.py
hpc_harbor/sieve/sieve_t5_sigmapoint.py
reports/sieve_t1_curvature.npz
reports/sieve_t4_closure.npz
reports/sieve_t5_sigmapoint.npz
```

I edited no tracked file, committed nothing, pushed nothing, and staged nothing outside this
list.
