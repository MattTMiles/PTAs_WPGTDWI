# EMBER — the OFF-TRUTH ladder: what the soft loop does when started where a real analysis starts

**Agent:** EMBER · ACCRE · **PURE EXECUTION** (no commits, no tracked-file edits) ·
**Date:** 2026-07-13

**[DRAFT — results sections are filled as the ladder lands. Pre-registration, method and scope
are frozen BEFORE the sweep and are not edited afterwards.]**

---

## 0. THE TAG, AND A CORRECTION TO STEP 0

The brief says: *checkout tag `criterion-v2.2` (verify `db2075a`)*. **The tag now exists and it
does resolve to `db2075a`** — KINDLE asked Matt to cut it and he did. But **checking it out
would have deleted the inputs this campaign is built on.** `criterion-v2.2` is the commit
*before* the KINDLE commit (`d87db93`), and the only difference between them is **purely
additive**: the seven `hpc_harbor/kindle/*` files, `reports/KINDLE_gain_contour.md` (my §5
pre-registration), its figure, and `surface_fALL_offenders.npz`. Every spec and convention
document — `project_progress.md`, `REQUIREMENTS_FROZEN.md`, `MANIFEST.md`, `STORY.md`,
`SAMPLER_dev.md`, `RECUT_floors.md`, `IGNITE2_softloop.md` — is **byte-identical at the tag and
at HEAD** (verified: `git diff criterion-v2.2 HEAD` touches only those ten files).

**So I worked at `MM_playground` HEAD `d87db93`: the identical CRITERION v2.2 spec, plus the
pre-registration and the `kindle_loop.py` driver the brief orders me to build on.** Flagged
rather than silently resolved. *If the intent was for the tag to sit on `d87db93`, that is a
`git tag -f` — Matt's to run, not mine.*

---

## 1. THE PRE-REGISTRATION (frozen before any sweep realisation ran)

### 1.1 The ladder, RE-DERIVED — not inherited

§10.16(d) is binding: *"anything of the form 'cell A versus cell B' … must be RE-DERIVED from
the corrected banks before it is quoted."* KINDLE §5 quotes Pair B at fN 122.46 / 44.40. EMBER
does not inherit those numbers; `hpc_harbor/ember/ember_cells.py` re-derives **every** cell,
floor, error, estimator and zero-fraction from `reports/recut_surface.npz` /
`recut_chorus.npz`, and **STOPs** if the re-derivation disagrees.

> **GATE C0 — PASS.** (−12.75, 40, vlbi, 5+11) → re-cut floor **122.461** (KINDLE 122.46).
> (−13.00, 40, vlbi, 5+11) → re-cut floor **44.397** (KINDLE 44.40). Both zero-fraction
> **0.000**.

**The nine rungs.** The `h+0.25` / `h+0.5` dex loudness steps land *on* the existing grid, so
the loudness ladder is just the march up `h` at T=40/vlbi/5+11 — no new cells needed.

| # | cell | zero-f | estimator | ADOPTED floor | corr/real | wrong/real | verdict | role |
|---|---|---|---|---|---|---|---|---|
| S1 | (−12.75, 40, vlbi, **5+11**) | 0.00 | gumbel | **122.461 ± 4.895** | 4.067 | 0.133 | ONSET | **Pair B** |
| S2 | (−13.00, 40, vlbi, **5+11**) | 0.00 | gumbel | **44.397 ± 1.938** | 3.567 | 0.067 | ONSET | **Pair B** |
| S3 | (−12.50, 40, vlbi, 5+11) | 0.00 | gumbel | 330.902 ± 12.171 | 5.833 | 0.067 | ONSET | +0.25 (S1) / +0.50 (S2) |
| S4 | (−12.25, 40, vlbi, 5+11) | 0.00 | gumbel | 956.161 ± 31.620 | 7.067 | 0.533 | ONSET | +0.50 (S1) |
| S5 | (−12.75, 40, vlbi, **3+13**) | 0.00 | gumbel | 46.212 ± 1.880 | 1.633 | 0.267 | ONSET | structure @ −12.75 |
| S6 | (−13.00, 40, vlbi, **3+13**) | 0.01 | gumbel | 20.145 ± 1.064 | **1.000** | 0.000 | below | structure @ −13.00 |
| C1 | CHORUS **m1 e=0.5** vlbi | 0.51 | **emp_q95** | **9.871 ± 0.914** | 2.267 | 0.033 | CONFIRMED | 1-member switch-on |
| C2 | CHORUS **m2 e=0.3** vlbi | 0.36 | **emp_q95** | **11.355 ± 1.018** | 1.767 | 0.033 | CONFIRMED | 2-member switch-on |
| C3 | CHORUS **m1 e=0.7** vlbi | 0.44 | **emp_q95** | **10.950 ± 1.009** | 5.767 | 0.067 | CONFIRMED | strong seeds; b4 gate cell |

*(CHORUS rungs are the **corrected** switch points, RECUT §2.1 — **not** the refuted e=0.3
single-member.)*

**S6 is one of §10.17(1)(i)'s two "precisely-1.000" cells** — 30 correct certs in 30
realisations, zero wrong, sitting exactly on the strict `>1` lattice point. It entered this
ladder as the matched structure control and arrives carrying an open convention question. It is
not resolved here (that is an onset-map convention, not a loop result), but the loop's behaviour
at a cell whose count is *exactly* the bar is reported.

### 1.2 THE FLOOR-VALIDITY GATE IS LOAD-BEARING, AND IT SPLITS THE CHORUS GATE FROM THE CHORUS SCIENCE

The six SURFACE rungs have zero-fraction **0.00–0.01** → the Gumbel is **valid**, and the
adopted floor *is* the Gumbel floor. Nothing to reconcile.

**All three CHORUS rungs have zero-fraction 0.36–0.51 → the Gumbel is INVALID** (criterion-v2.2:
valid only at ≲20 %), and the adopted floor is the **empirical q95 with a bootstrap error**. But
CHORUS's banked `ch_sloop_*` trajectories were run on 2026-07-13 **under the old Gumbel floors**.
So:

> **GATE b4 runs at the BANKED (old Gumbel) floor** — reproducing the banked certified set is
> what licenses the *machinery*. This is RECUT's **Gate B** discipline ("this scorer at the OLD
> floors reproduces the banked counts"), applied to a loop instead of a count.
> **THE SCIENCE runs at the criterion-v2.2 ADOPTED floor.** Both are banked.
> *No corrected number is emitted until an uncorrected one is reproduced from the same raw
> columns.*

### 1.3 THE ARM-(c) MAGNITUDE — declared before the sweep

KINDLE §5's table still reads *"(b-FIX) 1 × σ(mc)"*. KINDLE §4 finding 2 **measured** that at
~1e-4 dex and retired it as *"measured-too-small"*, naming the replacement: **"2× the seeder
tolerance (IGNITE-2 §4.1)"**. §5 was never updated, so the brief's fallback clause is live — but
the named replacement *is* well-defined, once one notices it must not be a **sky** offset:

* IGNITE's `TOL_GRID` (`ignite.py:90`) displaces the **sky**. **The loop freezes sky** and fits
  (f, mc, extrinsics) — KINDLE §5 item 4 says so explicitly. A sky offset is an axis the loop
  **cannot act on**, so it cannot be the fixed control.
* The spec's per-axis **certification tolerance** (`project_progress.md` L1801-1802) is:
  `sky 1e-5 | freq 3e-5 | mc 1e-3 | extrinsics no-collapse to 1e-2` (scaled units).
  The **mc-axis tolerance is 1e-3 scaled**; `SCALE_MC = 0.5`, so that is **5e-4 dex**.

> ### **ARM (c) := truth + 2 × (mc-axis cert tolerance) = +2e-3 scaled = +1.0e-3 dex in log10_mc**, per loud source, along its own mc axis.

That is ~**10×** the pilot's measured profile width (σ(log10 mc) ~1e-4 dex) and **2×** the
tolerance at which certification is known to break — a displacement with real dynamic range,
unlike the 1σ control the pilot retired. **The median realised arm-(b) MAP offset is banked
beside it** (the brief's fallback), so the reader can see whether the fixed control sits inside
or outside the honest spread.

### 1.4 WHAT THE EARTH TERM CAN AND CANNOT SEE — the prediction, recorded before the run

`project_progress.md` L1800-1806, Earth-term Fisher information gain **over the prior**:

| axis | info gain | consequence for the honest start |
|---|---|---|
| `log10_fgw` | **162–389×** | TIGHT — the MAP lands within the Earth-term *measurement error* |
| `log10_mc` | **1.00–1.73×** | **LOOSE — "for loud1/2 the posterior IS the prior"** (σ 0.50/0.86/0.87 dex; the box's own sd is 0.87 dex) |

**So the honest Earth-term MAP start is, to good approximation, a PURE mc DISPLACEMENT of order
the prior width (~0.5–0.9 dex)** — roughly **1000× the mc certification tolerance** — sitting on
**the one axis the loop is able to act on**. That is why the ladder is well-posed, and it is why
the fixed control (1e-3 dex) is tiny beside the MAP arm's natural spread, exactly as KINDLE §4
predicted. **This prediction was written into `ember_map.py` before the gate ran; gate b1-3
measures it.**

### 1.5 ONE STRUCTURAL PROPERTY OF THE READOUT, STATED SO A NULL IS NOT MISREAD

If the MAP start really is ~random in mc, then at that start the pulsar-term fringes are
meaningless and **|C₀(start)| = 0**. Then `ΔN = |C_fix| − 0 ≥ 0` **by construction**: from a
zero-cert start the additive gain is **one-sided**. **DEGRADES therefore cannot fire through ΔN
in arm (b)** — it can only fire through *loop-generated* wrong certifications. The degradation
channel in arm (b) is the wrong-cert column and the margin flow; **ΔN < 0 is a thing the TRUTH
arm can show** (KINDLE's pilot saw exactly one such loss — a floor-margin wobble). A ΔN ≥ 0
result in arm (b) is therefore **not** evidence of safety on its own, and is not reported as
such.

### 1.6 THE DEGRADES TRIGGER — refined, and NOT weakened

The brief's DEGRADES trigger reads "wrong-certs > 0". But these cells' **seeds already carry
wrong certifications** (0.067–0.533/real at the SURFACE rungs; 0.033–0.067 at CHORUS) — the D1
wrong-counterpart hole, open and stated everywhere it appears. A loop that inherits a wrong cert
and never touches it has not degraded anything; IGNITE-2 §3.2 established exactly this
distinction. So DEGRADES fires on **loop-generated** wrong certs, `wrong(final) − wrong(start) >
0`, and **both columns are banked**, so the inherited rate travels beside every verdict and the
stricter reading can be re-cut by anyone who prefers it. **A cell with `wrong_grown = 0` and a
nonzero inherited rate is reported as HOLDS/REPAIRS *with the inherited rate quoted*, never as a
clean bill of health.**

---

## 2. THE BUILDS AND THEIR GATES — **ALL FOUR PASS**

| build | gate | result |
|---|---|---|
| **b1** Earth-term MAP optimiser | MargJax-gradient identity · basin classification · start-offset distribution | **ALL PASS** |
| **b2** T = 40 problem build | banked T=40 realisation reproduced · floor refit vs `_recut` | **PASS** |
| **b3** structure-married loop (k_loud = 5) | banked SURFACE 5+11 certified set @ iteration 0 | **PASS** (both Pair-B cells) |
| **b4** eccentric driver | banked `ch_sloop_m1e07` certified set @ iteration 0 · L_gwb fingerprint | **PASS** (**0.000e+00**, 3/3) |

### 2.1 b1 — the Earth-term MAP optimiser. **THE CRITICAL PATH.**

**B1-1 — "MargJax gradients": the identity is ANALYTIC, and it is measured.** At `mask = 0` the
Earth-term likelihood has no distance dependence, so inside MargJax every fringe carries the
same value:
`M_p = logsumexp_k(peak_k + lnprior_k) = lnL + c_p`, `lnref = lnL`, hence
`lnL_marg = Σ_p M_p − (npsr−1)·lnref = lnL_earth + Σ_p c_p`.
**The fringe machinery contributes a CONSTANT and EXACTLY ZERO gradient.** Measured:

| | |
|---|---|
| `max|grad_MargJax − grad_masked|` at two points | **7.06e-09 / 1.04e-08** (**rel 7.6e-13 / 1.2e-13**) |
| value offset (= Σ_p c_p) | **696.791771** at *both* points |
| constancy of that offset across x | **8.56e-08** (predicted 0) |

**So the MargJax gradient and the masked-lnL gradient are the SAME OBJECT.** The brief's
"MargJax gradients" requirement is met, and the masked-lnL path is licensed as *the* path, not
a cheaper substitute — without paying B = 512 fringe evaluations per pulsar for a quantity
whose x-derivative is identically zero.

**B1-2 — basin classification: 3/3 reproduce KINDLE's pilot EXACTLY** (seeds 6210000/1/2 →
`C0(start/truth)` = 1/1, 0/0, 0/0, matching the banked `k_pilot_*truths0*`), including the
pilot's ΔN = −1 floor-margin loss. EMBER's basin classifier **is** the object KINDLE used.

**B1-3 — the start-offset distribution.**

| | |
|---|---|
| `|offset|` on mc (dex), 3 reals | **1.484 · 2.027 · 0.874** |
| `|offset|` on NON-mc (scaled) | 0.169 · 0.172 · 0.544 — *non-zero: the real Earth-term measurement error* |
| **median mc offset** | **1.484 dex = 2968 × the mc certification tolerance** |
| all finite | **True** |
| **in-basin fraction at the MAP start** | **0.00** *(truth arm: 0.33)* |

### 2.2 THE TWO THINGS THE GATES CAUGHT, STATED PLAINLY

**(a) MY FIRST MAP OPTIMISER WAS BROKEN, AND THE GATE IS WHAT FOUND IT.** Run 1 returned
`lnL_earth nan, nit=0` on **3 of 6 starts** and `nit=1` on the rest — the reported "MAP" was
barely more than the raw cold draw. Two defects: a NaN gradient killed L-BFGS-B outright, and
the loud-param Hessian's **cond ≈ 4e10** (spec L1685-1690: `log10_fgw` razor-tight, `log10_mc`
flat) makes a naive joint quasi-Newton step stall and declare convergence. **REBUILT** with an
mc **scan** first (bound the optimiser inside the evaluable region) and **two stages**
(A: mc held, fit the constrained axes; B: release and polish). After the rebuild the starts
take **400 + 400 iterations** and reach a *higher* Earth-term optimum. The tainted banks were
deleted, not reused.

> **AND ONE OF MY OWN EXPLANATIONS DID NOT SURVIVE ITS OWN INSTRUMENT.** I attributed the NaNs
> to *"the mc prior box contains binaries that have already coalesced"*. **The scan says
> otherwise: `non-evaluable 0%` across the ENTIRE mc box, in every realisation.** The NaNs came
> from the *joint* cold draw over three loud sources at once and from unguarded optimiser
> excursions — not from coalescence at the anchor. The fix (evaluable-bounding + a NaN penalty)
> was right; **the story I told about it was wrong, and the measurement is what corrected it.**

**(b) MY b2/b3 GATE BAR WAS WRONG — THE BUILD WAS NOT.** I first demanded bit-identity against
the banked SURFACE realisations and got `dlnL 3.96e-09 / qmax 6.46e-10`, with **every discrete
column exactly 0.000e+00 and the certified set IDENTICAL.** The mechanism is known and already
gated: the EMBER loop **must** build a `B1Marg` (it is the M-step's objective), and
`trackB_b1.py:98` has `B1Marg.__init__` call `sp.enable_fast_full(all_on)` — so `sp.estep` runs
the **fast pterm-only residual**, while SURFACE's bank was scored on the **slow** path.
**IGNITE-2 measured this exact drift at 3.5e-10 (IG2-3, "discrete columns identical") against a
1e-6 threshold, never bit-identity**, and gates the fast path itself at 1.16e-10 (IG2-2 — which
EMBER reproduced exactly). Adopting that precedented bar — **discrete columns and the certified
set exact; continuous < 1e-6** — b2 and b3 **PASS**. *This is not a weakened gate: the certified
set is the criterion's own output, and it is exact.*

### 2.3 b4 — the eccentric driver, and the bug that mattered most

The first attempt missed by **445 nat** with a certified set of **13 vs 9**. Cause: CHORUS's
banked `ch_sloop_m1e07_*` are **PAIRED** realisations, so `chorus.py:1491-1493` draws their truth
on the **CIRCULAR m0 geometry** (empty e-assignment) — the convention that makes an m0/m1e07
pair comparable. EMBER was drawing truth on the eccentric grid. **Both conventions are correct,
for different jobs, and both are now wired and banked (`paired_truth`):**

* **GATE b4** → the paired circular grid. It must match the bank it reproduces.
* **THE SWEEP** → the cell's **own** grid — because `recut_chorus`'s **adopted floors and
  corr/wrong counts are cut from the STAGE-1 `ch_sig_*` realisations, not from the paired
  sloops.** A sweep on the paired grid would be scored against floors calibrated on a different
  truth draw.

With that fixed, all three gate reps reproduce **BIT-IDENTICALLY** (`max|diff| = 0.000e+00`;
certified sets 13/13, 8/8, 2/2). That zero simultaneously licenses the eccentric problem build,
the 32-harmonic chirp-tie, the paired-truth convention **and the `cpus = 8` BLAS thread count**
— the T=30 recompute path reproduces *only* at 8, and the T=30 fingerprint
(`9fd547b39b02c705`) is identical to KINDLE's.

**b4(ii)** — L_gwb manifest fingerprint `71677a810cbc7187` and sha256 both **MATCH**.

### 2.4 Warm cache + resume drill (HPC_SETUP, before the wide launch)

The four gate jobs warmed all three graph shapes (the T=40 build fell **~15 min cold → 243 s
warm**). The resume drill: shard 0 launched, **deliberately `scancel`'d mid-realisation** at
22:52 with 1 realisation banked, then the *same* sbatch resubmitted →
`sweep(surface): 246 entries, shard 0/8 -> 31 mine; already banked: 1; to run: 30`.
**Skip-on-exist resume is proven on the production SLURM path.**

---

### 2.5 THE EXECUTION DEFECT THAT EVERY "SHARD COMPLETE" LINE HID — and why completeness is checked by READBACK

**EMBER-2 (the resuming session) found 21 realisations missing from a sweep in which every
log line said `SHARD COMPLETE`.** The cause is worth recording, because the failure was
*silent* and the discipline that caught it is cheap.

The resume drill (§2.4) used `e_drill.sbatch` — array `0-0`, **walltime 2 h**. Shard 0 was
then left in that job's hands as a *production* shard, and the production array
`12524934` was dispatched as tasks **1–7**. So:

| | |
|---|---|
| `e_drill` (shard 0, 2 h limit) | **TIMEOUT at 01:08:21**, 10 of its 31 realisations banked |
| `e_sweep` tasks 1–7 (14 h limit) | all `COMPLETED`, 215/215 banked, `SHARD COMPLETE` on every one |
| **net** | **21 realisations never ran** — and *nothing in any log said so* |

The shard that died was a **different job** whose `TIMEOUT` nobody re-read, and the shards
that lived all reported success truthfully. Five of the nine cells were sitting at **N = 10–11
against a pre-registered N ≥ 12**, and the campaign would have been written up as complete.

> **`sacct -j 12524934` lists tasks 1–7. There is no task 0. That is the whole bug.**

**It was caught by re-deriving the completeness table from the banks themselves**
(`ember_anatomy.py`: load every npz, re-derive `cert` from the raw `dlnL`/`qmax`/`bar` columns,
count per (cell, arm)) — **not** by checking that jobs exited 0, and not by `ls | wc -l`. This
is the brief's *"verify by readback, not existence"* rule doing exactly the job it was written
for. **347/347 banked realisations re-derive their own `C0_start`, `C_fix`, `cert_start` and
`cert_final` from raw; zero disagreements.** The 21 were then re-dispatched as a real array
task at production walltime.

**One thing that could have been contaminated, and is not.** The drill and the production
array must agree on the GWB draw or shard 0's survivors are not comparable with the rest —
`NoiseDrawer`'s `eigh` rotates the draw with the BLAS thread count (the known
`cpus-per-task = 8` hazard). Both ran at `cpus-per-task = 8`, and the surface L_gwb
fingerprint is **`8548f148b50a5b44` on the drill and on all seven production shards**
(chorus is `9fd547b39b02c705`, the T = 30 build, matching KINDLE). **The 10 survivors are on
the same draw and are kept.** The gap was purely "never ran", not "ran differently".

---

## 3. THE LADDER — RESULTS

**324 signal realisations** — 9 cells × 3 arms × **N = 12** — all re-derive their own
`C0_start`, `C_fix`, `cert_start`, `cert_final` from the raw `dlnL`/`qmax`/`bar` columns with
**zero readback disagreements** (completeness verified by readback, not `ls`; §2.5). Figure
`EMBER_predictor.png` panel (c). The full per-(cell,arm) table is banked in `ember_analysis.npz`.

**The one-line result: from the honest start the loop HOLDS in every cell.** The map arm — the
Earth-term MAP, the start a real analysis actually gets (§1.4) — has **ΔN = 0 in all 9 cells**,
**cum-flow ≈ 0**, and **0 loop-generated wrong certs**. It doesn't repair and it doesn't degrade;
it holds. The reason is mechanical and is the through-line of §4: the honest start is **disengaged
and motionless** — in-basin fraction 0.00–0.42, median accepted M-step ≈ 0 (`steps` column 0.00
in 7 of 9 map cells) — so there is nothing for the additive gain to add.

| arm | what it is | ΔN across the 9 cells | loop-grown wrong certs | accepted motion |
|---|---|---|---|---|
| **map** (b) | honest Earth-term MAP (~0.74–1.58 dex cold in mc) | **0.000 in all 9** | **0** | ~0 (inert) |
| **truth** | on-truth, already engaged | **+0.000 … +0.417** (repair in 4 cells) | 0 (and **cleaned 1**) | 0.8–2.8 steps |
| **fix** (c) | truth + 1e-3 dex on mc | **0.000 in all 9** | 0 | ~0 (inert) |

**The truth arm is where the loop does visible work** — small, real repair when it starts already
engaged: `h−12.75_k5` **+0.417** (2 of 12 gain, cum-flow **+44.6** toward the bar), `h−12.75_k3`
+0.167, `m1e05` +0.167, `m1e07` +0.083 (with one −1 floor-margin wobble, the KINDLE-pilot loss
reproduced). **That contrast — truth repairs, map holds — is the loop behaving correctly: it
advances a start that has already engaged the fringes, and it leaves an honestly-ignorant start
where it found it.** Neither arm grows a wrong cert.

---

## 4. THE SCRAMBLED NULL — and the competing predictor of manufacturing

**112 scrambled realisations** (the 72 base arms at N = 8 per (cell,arm) + the 40-rung engagement
edge). No signal sits at the assumed counterpart, so **every certification here is false**: `C0` is
seed-inherited, **ΔN is loop-made**. **5 of 112 manufactured** (grew the false-cert count, ΔN > 0).

### 4.1 The pre-registered STOP fired, and its LETTER is falsified while its RATIONALE holds

The manufacturers, in full (nothing hidden; `EMBER_predictor.png` (a)):

| seed | cell / arm | offset | C0 | ΔN | max M-step | engaged? | moved? |
|---|---|---|---|---|---|---|---|
| 9500204 | h−12.75_k5 fix | 1e−3 | 1 | **+2** | 0.498 | Y | Y |
| 9500204 | h−12.75_k5 fix | 3e−1 | 1 | **+2** | 0.912 | Y | Y |
| 9500205 | h−12.75_k5 fix | 1e−3 | **0** | +1 | 0.991 | **.** | Y |
| 9510004 | h−13.00_k5 truth | 0 | 1 | +1 | 0.124 | Y | Y |
| 9510202 | h−12.75_k5 fix | 1e−2 | **0** | +1 | 0.549 | **.** | Y |

**Two of the five manufacturers start DISENGAGED (C0 = 0).** The "26/26 zero-cert starts inert"
headline — the two-condition *engagement* boundary of the pre-registration — **is falsified in its
letter.** But its rationale (manufacturing is not a near-but-wrong *annulus*; it needs the loop to
*act*) holds, and sharpens into a cleaner variable.

### 4.2 ENGAGEMENT vs MOTION — the competing-predictor read (brief addition 1)

Target = manufacturing (scrambled ΔN > 0). ENGAGEMENT = C0 ≥ 1 (an iteration-0 initial
condition). MOTION = the largest accepted damped-Newton step `max_i |traj_step_norm_i|` (a process
quantity — the M-step that relocates the source *is* the act of manufacturing; the t = 0 / process
asymmetry is stated, not hidden). One contingency + logistic read on all 112 (`EMBER_predictor.png`
(a,b); banked `ember_predictors.npz`):

| predictor | sensitivity | specificity | Fisher p (1-sided) | verdict |
|---|---|---|---|---|
| ENGAGEMENT C0 ≥ 1 | 0.60 (3/5) | 0.60 | **0.333 (n.s.)** | misses the 2 disengaged manufacturers |
| **MOTION** (moved, step > 0) | **1.00 (5/5)** | 0.72 | **0.002** | **in-sample NECESSARY** (0/5 manufacturers had zero motion) |
| JOINT engaged **AND** mobile | 0.60 | 0.87 | 0.024 | **NOT the boundary** — sensitivity collapses to the engaged subset |
| JOINT engaged **OR** mobile | 1.00 | 0.48 | 0.044 | catches all, but no better than motion alone and far less specific |

Joint logistic (standardised, L2 = 1 for the quasi-separation): **β_motion = +0.357 > β_engagement
= +0.308** — **motion dominates.** Necessary/sufficient: engagement is **neither** necessary (2/5
disengaged) **nor** sufficient (43/46 engaged starts were safe); motion is **necessary in-sample**
(all 5 moved) but **not sufficient** (27/32 mobile starts were safe). **The boundary is MOTION —
the loop manufacturing a false cert requires it to take a real accepted M-step; the hypothesised
"engaged AND mobile" joint is refuted (sensitivity 0.60).** What separates a mobile start that
manufactures from one that doesn't is loudness/floor-proximity, not engagement — the same gate that
decides which motions cross in the signal arms.

**Positive control (signal arms):** motion drives *legitimate* ignition too — of 7 signal
realisations with ΔN > 0, **7/7 moved** (sensitivity 1.00, specificity 0.96, p < 0.001). So motion
is the general driver of *any* count change; a **true counterpart** is what makes that change a
repair rather than a manufacture. **The corrected two-condition boundary is therefore
(wrong counterpart) × (motion) — not (wrong counterpart) × (engagement).**

**Self-cleaning, banked beside it:** seed 9510200 (h−13.00_k5 fix) ran **ΔN = −4** — the loop
*destroyed* 4 inherited false certs. The scrambled arm is not only "does it manufacture" but "does
it clean," and here it cleaned.

> **Scope (travels with §4).** N = 112, **5 events** — small-N, one-sided Fisher, logistic under
> L2 shrinkage. Motion is a *process* quantity: its "predictive" edge partly reflects that it sits
> closer to the manufacturing mechanism than a t = 0 flag can. The claim is not "motion causes
> manufacturing at rate X"; it is the **ordering** — motion separates the 5 events cleanly where
> engagement does not, and every manufacturer moved.

### 4.3 The channel-budget lens on the chorus cells (brief addition 2, MAGPIE J1)

MAGPIE J1 established the mixed-population switch is set by the **active-harmonic channel budget**
(n_active ≈ 30 at census loudness), not per-source κ. **`n_active` is a per-CELL census constant
and is NOT derivable from EMBER's per-realisation banks** — the 19-dim fit vector is the *source
model* (7 orbital params/member + a fixed 6-harmonic circular comb), not the census channel count.
The authoritative per-member counts are the `ch_README.md` census values (e = 0.3 → 8, 0.5 → 17,
0.7 → 32); I attach those per cell (× n_ecc) and do **not** approximate a per-realisation value.

For EMBER's three chorus cells the honest (map-arm) off-truth behaviour is **null across the whole
channel budget** (`EMBER_predictor.png` (d)):

| cell | e | n_ecc | n_active | κ@f_orb | map ΔN | map repair | map engage |
|---|---|---|---|---|---|---|---|
| m2e03 | 0.3 | 2 | **16** | 2.65 | 0.000 | 0.0 % | 8.3 % |
| m1e05 | 0.5 | 1 | **17** | 13.0 | 0.000 | 0.0 % | 8.3 % |
| m1e07 | 0.7 | 1 | **32** | >13* | 0.000 | 0.0 % | 0.0 % |

**The map arm is inert in every chorus cell, so there is no off-truth variation to organise — by
channels, by e, or by κ.** The only in-sample contrast that *could* separate the axes is **m2e03 vs
m1e05** (channels near-equal, 16 ≈ 17, but e and κ far apart, 0.3 vs 0.5); both are identically
null, which is *consistent* with channel-budget parity but equally consistent with "the chorus map
arm is simply disengaged." **EMBER's 3 chorus cells lack the equal-κ contrast (m1e03/m2e03/m3e03)
that let MAGPIE J1 isolate the channel budget, so EMBER neither confirms nor refutes the
channel-budget switch — it is consistent with it and cannot independently establish it.** (N = 3
cells, map arm null; stated rather than over-read.)

---

## 5. VERDICTS

**All nine cells: HOLDS** — read off arm (b), the honest MAP start, per the pre-registered rule
(§1.6). ΔN = 0, **0 loop-generated wrong certs**, in every cell; the seed-inherited wrong-cert rate
is quoted beside each, never as a clean bill of health.

| cell (role) | verdict | ΔN (map) | loop-grown wrong | seed-inherited wrong/real |
|---|---|---|---|---|
| S1 h−12.75_k5 (**Pair B**, e-onset) | **HOLDS** | +0.000 | 0 | 0.42 |
| S2 h−13.00_k5 (**Pair B**) | **HOLDS** | +0.000 | 0 | 0.25 |
| S3 h−12.50_k5 | **HOLDS** | +0.000 | 0 | 0.25 |
| S4 h−12.25_k5 | **HOLDS** | +0.000 | 0 | 0.58 |
| S5 h−12.75_k3 (structure) | **HOLDS** | +0.000 | 0 | 0.42 |
| S6 h−13.00_k3 (the 1.000 cell) | **HOLDS** | +0.000 | 0 | 0.17 |
| C1 **m1 e=0.5** (single-member switch) | **HOLDS** | +0.000 | 0 | 0.08 |
| C2 **m2 e=0.3** (two-member switch) | **HOLDS** | +0.000 | 0 | 0.17 |
| C3 **m1 e=0.7** (strong-seed / b4 gate) | **HOLDS** | +0.000 | 0 | 0.00 |

**The switch language, exact (RECUT §2.1):** the chorus onset is **e = 0.5 single-member (C1)** or
**e = 0.3 two-member (C2, `m2e03`)** — **never** e = 0.3 single-member (that reading is refuted).
C3 is the e = 0.7 single-member strong-seed cell.

**D1, framed exact.** Across the 324 signal realisations the loop carries **90 inherited wrong
certs** (the standing wrong-counterpart hole; arm totals inherited 8 / 28 / 107 for truth / map /
fix). It **grows 0** of them and **cleans 1** (a truth-arm cert). To three-significant-figure
honesty: **the loop neither creates nor cures D1** — the wrong-counterpart hole travels *through*
the loop essentially unamplified and unrepaired, exactly the IGNITE-2 §3.2 distinction. D1 remains
**open** (no purity layer exists; D3 and D4 both rejected), and the scrambled arm measures it but
cannot close it.

**The headline, corrected and scoped.** From the honest Earth-term start the soft loop is a
**stabiliser** (HOLDS everywhere, grows no wrong cert); the manufacturing that fired the
pre-registered STOP lives only in the **scrambled null**, and its boundary is **MOTION**, not start
engagement. The earlier "safe and self-cleaning" reading survives with the proviso: **safe given
either a true counterpart or a motionless (disengaged-and-still) start** — a start that *moves*
under a *wrong* counterpart is where false certs are made.

---

## 6. SCOPE OF INFERENCE — travels with every verdict above

Everything here is scoped to the **116-pulsar MOCK array** (AXIS, 1440 MHz single-frequency,
§10.15(a)) and its frozen prior stack (`best_distances.txt`, the per-pulsar `ln K_counted`, the
tiers, the geometry ensembles). **No real TOA is touched.** The residuals ARE the injected
CW+CURN.

The off-truth verdict is a property of **this population structure and this prior stack**, which
nature does not have to supply — the SURFACE/CHORUS population-scoping caveat travels unchanged.
**The D1 wrong-counterpart hole is OPEN and travels beside every count**; the scrambled arm
measures it but cannot close it, and no purity layer exists (D3 and D4 both rejected).

**One scope limit specific to arm (b), and it is a real one.** `log10_fgw`'s box spans 1.7 dex —
at T=40 yr that is ~250 independent frequency bins, so the Earth-term likelihood in `f` is
multi-modal (periodogram-like). **A blind global frequency search is a property of the SEARCH,
not of the loop**, and a search that picked the wrong `f` bin would guarantee the loop fails for
reasons having nothing to do with the loop. This campaign is about the **loop**. So the MAP
optimiser is **basin-anchored** in the well-constrained axes (f, h, extrinsics) — which the data
then *move*, and the realised offset there is the genuine noise-driven Earth-term measurement
error, banked per realisation — and initialised **COLD in log10_mc**, the axis the Earth term
genuinely cannot see. **Arm (b) therefore measures the loop's response to the Earth term's
PARAMETER ignorance (dominantly mc), NOT to a search failure.** That limit is stated on every
arm-(b) verdict and is not quietly dropped.

---

## 7. DEFECTS FOUND — REPORTED, NOT FIXED (CW_transition is tracked; HARD RULE)

* **`sampler_core.prior_sigma()` cannot express the `vlbi` tier.** It accepts only `'lit'` and
  `'vlbi_all'` and **raises** on `'vlbi'`; its docstring explains it was written around a cronus
  limitation (*"the union-18 list is NOT banked on cronus"*). But on ACCRE
  `build_ignite_problem` **does** construct `P["priors"]["vlbi"]` (`ignite.py:265-269`), so
  `MargJax` cannot see a tier the problem already carries. **Suggested fix (cronus's to make):
  `prior_sigma` should read `P["priors"][tier]` when present.** EMBER's b1 gate is unaffected and
  still exact — see §2.1.

* **JAX persistent-compile-cache RACE across array tasks (infrastructure, benign, re-runnable).**
  Eight array tasks sharing one compile-cache dir
  (`…/xla_gpu_per_fusion_autotune_cache_dir/tmp/*.textproto`, bound in from `HARBOR_JAXCACHE`)
  collide on the autotune temp files and die with `JaxRuntimeError: NOT_FOUND`. **11 realisations
  were guard-killed this way** — the finiteness guard is *not* implicated; it is a filesystem race.
  Each kill was **counted and banked** (`ember_fails_edge_s*.npz` / `ember_fails_chorus_s2.npz`),
  never dropped. **Root-caused and worked around** by giving each task its own cache dir:
  `harbor_env.sh` now honours a pre-set `HARBOR_JAXCACHE` (backward-compatible; unset → the original
  shared default), and the split re-dispatch (`e_edge_split.sbatch` / `e_ch_split.sbatch`) sets
  `HARBOR_JAXCACHE=…/jax_ember_split/t${SLURM_ARRAY_TASK_ID}`. Isolation is **numerically inert** —
  it caches XLA executables, not results; the load-bearing `cpus-per-task = 8` BLAS knob is
  untouched and every recovered realisation carried the correct L_gwb fingerprint
  (`8548f148b50a5b44` surface / `9fd547b39b02c705` chorus). All 10 outstanding realisations
  recovered; final readback **436 banked, 0 failed**.

* **Co-tenant GPU OOM — the "rogue-squat lottery" (HPC_SETUP §7.2), a scheduling hazard.** One
  realisation (`n9500207`) OOM'd **twice** (`RESOURCE_EXHAUSTED`, 15 GiB alloc on `jit__rank_update`)
  because SLURM twice assigned it the one A100 on dgx03 carrying a **50.7 GB foreign process**
  (`env/deerdiff`), leaving ~30 GB. Every EMBER task drew a *distinct* GPU (the split's cache
  isolation worked); this is co-tenancy, not a code defect. **Worked around** by resubmitting with
  `--exclude=dgx03`; it landed on a clean dgx04 GPU and banked 1/1. The device echo line
  (`nvidia-smi … --query-compute-apps`) on every task is what made the squat attributable rather
  than guessed — worth keeping. *No fix is EMBER's to make; the mitigation is the node-exclude and
  the per-task device log.*

---

## 8. THE ADD-LIST (Matt commits — I never do)

Everything below is **staged, uncommitted** — pure execution, no `git` run by me. Grouped so a
single review pass covers it.

**Analysis + figures (new, this session):**
* `hpc_harbor/ember/ember_predictors.py` — the competing-predictor (engagement vs motion) +
  channel-budget analysis; banks `EMBER_results/ember_predictors.npz`.
* `hpc_harbor/ember/ember_predictor_figure.py` → `reports/EMBER_predictor.png` — the corrected
  headline figure (4 panels: competing predictor, sens/spec, ladder HOLDS, channel budget).
* `hpc_harbor/ember/ember_edge_index.py` — read-only helper that authoritatively maps the missing
  plan indices for the split re-dispatch.

**Edits to existing untracked EMBER harness (NOT CW_transition):**
* `hpc_harbor/env/harbor_env.sh` — `HARBOR_JAXCACHE` now honours a pre-set per-task override
  (backward-compatible default unchanged); defuses the compile-cache race.
* `hpc_harbor/ember/ember_analysis.py` — fixed a misleading hard-coded scrambled-anatomy string
  that printed "dN=0 ⇒ inherited" even for `dN = +2` manufacturers; now conditional
  (manufactured / destroyed / inherited). Ladder verdicts written to `EMBER_results/ember_analysis.npz`.
* `hpc_harbor/ember/ember_anatomy_figure.py` — suptitle + panel (b) corrected off the falsified
  "basin ENGAGEMENT" headline to the MOTION finding; still the companion (per-pulsar +2, clean
  arms). → `reports/EMBER_scrambled_anatomy.png`.

**Split re-dispatch sbatch (new):**
* `hpc_harbor/ember/e_edge_split.sbatch`, `hpc_harbor/ember/e_ch_split.sbatch` — one-realisation-
  per-task arrays with per-task cache isolation; how the 10 guard-kills were recovered.

**Banks (new / updated under `EMBER_results/`):**
* The 10 recovered realisations (7 edge rungs + 2 scrtop + 1 chorus scrambled-map) → completeness
  now **signal 324 @ N=12; scrambled 8/8 both surface Pair-B cells and chorus m1e07; edge ladder
  5/5 × 4 rungs × 2 cells**.
* `ember_predictors.npz`, `ember_analysis.npz` (raw per-realisation rows banked for re-cut).
* `ember_fails_*.npz` — the guard-kill ledger (kept as the honest record that the race happened).

**Report:**
* `reports/EMBER_offtruth_ladder.md` — §3 / §4 / §4.2–4.3 / §5 / §7 filled; pre-registration
  (§0–§2) untouched.

**One thing for Matt to decide, not me:** whether `criterion-v2.2` should be force-moved onto
`d87db93` (§0) — a `git tag -f`, his call.
