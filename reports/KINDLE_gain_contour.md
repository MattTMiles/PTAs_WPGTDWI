# KINDLE — the gain contour, and why the question had to be reasked

**Agent:** KINDLE · ACCRE · **PURE EXECUTION** (no commits, no tracked-file edits) ·
**Date:** 2026-07-13

**Tag note (read first).** The brief instructs checkout of tag `criterion-v2.2`. **That tag
does not exist** — locally or on `origin` (newest is `criterion-v2.1`). But the *content* of
criterion-v2.2 is present at branch HEAD `db2075a` "floor adopted, onset recut to 59":
`project_progress.md` §10.16 documents the convention and even names the not-yet-cut tag.
I worked at `db2075a` and flag the missing tag for Matt to cut. All floor-validity-gate,
zero-fraction, and re-cut conventions used here are read from §10.16 / `reports/RECUT_floors.md`.

**Scratch (host):** code `hpc_harbor/kindle/` (`kindle_loop.py`, `kindle.py`,
`kindle_d7_fall.py`, sbatch), results `KINDLE_results/`, logs `hpc_harbor/logs/k_*`. Lane
`-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1`,
**`--cpus-per-task=8`** (the NoiseDrawer BLAS-thread hazard, §10.14(g) — banked noisy
realisations reproduce bit-identically only at 8).

**Staged for Matt (explicit add-list at the end).** This report; `reports/surface_fALL_offenders.npz`
(the D-7 closure bank); `KINDLE_results/` (g1 diagnosis npz, D-7 re-cut, g0 + pilot banks);
the figure. Raw statistics always banked (per-iteration `iter_dlnL/iter_qmax/iter_mapk`,
prior-referenced frame; additive-gain and margin-flow columns per the redesign).

---

## 0. THE ANSWER IN ONE PARAGRAPH

**The brief's premise — "gain is EXACTLY 1.000 in all 60 finite honest runs, perfectly
marginal" — is an artifact of a degenerate statistic, not a measured property, and the
diagnosis (question B) had to be answered before the sweep (question A) could mean anything.**
The banked "gain" (`ignite2.py:501`) is a RATIO of certified-set cardinalities
`|C_it| / |C_{it-1}|` with a divide-by-zero guard that returns NaN. Measured across all **70**
banked loop trajectories (IGNITE-2's 40 + CHORUS's 30): every *ignition* (|C|: 0→1) divides by
zero and is **filtered out as NaN**; every *finite* value is an `n/n` **tautology** equal to
1.000 by arithmetic; and the certified count has **never once gone n→n+1 with n≥1** in 50 honest
runs, so ratio-gain > 1 is **unreachable**, not merely unobserved. Underneath, **25 of 30 honest
IGNITE-2 loops accept zero M-steps** — `dlnL` is bit-identical across iterations — because the
loop is **initialised at the true source, i.e. at its own fixed point**, where the gradient is
~0 and the line search correctly refuses to move. So the "1.000 marginality" is §10.17(1)(i)'s
own suspicion confirmed and then some: *not* a strict-inequality quantisation of a positive
sub-threshold gain, and *not* genuine inertness at a live ignition front — it is a **degenerate
readout of a loop run at rest at truth.** The gain contour, as posed, is the contour of a
quantity that cannot exceed 1 by construction. **The redesign (Matt, this session): retire the
ratio; read the loop with additive gain ΔN = |C_fix| − |C_0| and continuous margin flow; and
run the sweep from OFF-TRUTH starts — the only initialisation under which "does the array
self-amplify" is a well-posed question.** The two NaN'd 0→1 events are re-labelled as what they
are: **genuine ignitions.** [g0 license + off-truth pilot results: §3–4, pending GPU.]

**Second, independent deliverable — D-7, the last unfinished piece of the floor fix, is CLOSED,
and it closes RECOVERED rather than unrecoverable-from-disk.** §2.

---

## 1. THE MARGINALITY DIAGNOSIS (g1) — the ratio gain is degenerate by construction

### 1.1 What the statistic is

`ignite2.py:501`, verbatim:

```python
gain = (len(cert_sets[-1]) / len(cert_sets[-2])
        if len(cert_sets) > 1 and len(cert_sets[-2]) else np.nan)
```

Gain is `|C_it| / |C_{it-1}|` — a ratio of the sizes of the certified set at successive
iterations — with a guard that returns **NaN when the previous set is empty**. That guard fires
on exactly the event the whole campaign is looking for: an ignition from zero certifications.

### 1.2 The measurement (all 70 banked trajectories)

Banks: `reports/ig_sloop*.npz` (IGNITE-2, 30 signal + 10 scrambled) and
`CHORUS_results/ch_sloop*.npz` (CHORUS, 20 signal + 10 scrambled). Reading every per-iteration
`traj_gain` entry against the `traj_n_cert` transition that produced it:

| |C\| transition | what it is | ratio gain | count (honest runs) |
|---|---|---|---|---|
| 0 → 0 | stays cold | NaN (guard) | 38 |
| **0 → 1** | **IGNITION** | **NaN (guard)** | **2** |
| 1 → 0 | loses its one cert | 0.000 | 1 |
| n → n (n ≥ 1) | **holds** | **1.000 (tautology)** | 49 |

- **50 honest runs, 50 finite gain entries: 49 × 1.000 + 1 × 0.000.** Not "all 1.000" — the
  brief's "60 finite honest runs, all exactly 1.000" does not match the banks (there are 50
  honest runs, not 60; one reads 0.000; and the 1.000s are tautologies).
- **Every finite 1.000 is an `n/n` self-ratio.** The count did not change that iteration, so the
  ratio is 1 by arithmetic — it carries no information about gain.
- **Every ignition is NaN**, filtered out of the "finite" set. *The "finite honest runs" filter
  selects against the phenomenon it is meant to detect.*
- **The count has never gone n→n+1 with n ≥ 1 in 50 honest runs.** Ratio-gain can only exceed 1
  if an already-certifying realisation adds a second cert in one step; that has zero instances.
  **Ratio-gain > 1 is unreachable in the banked box, not marginally unmet.**

### 1.3 The continuous substrate — the loop is at rest, not poised

g1's second question: does the sub-threshold `dlnL` of the near-bar pulsars move *toward* the
bar (a positive continuous gain, quantised away by the discrete count) or not at all (genuine
inertness)? Measured on the 30 honest IGNITE-2 trajectories (margin `m_p = dlnL_p − max(lnK_p,
floor)`; ignition ⇔ `m` crosses 0):

- **25 of 30 loops accept ZERO M-steps** — every proposed step is rejected by the line search,
  so `dlnL` is **bit-identical across all iterations** and the continuous gain is **exactly
  zero**, not sub-threshold-positive.
- In the 5 loops that do move, margins shift **both ways**: of 3480 (pulsar × run) margins,
  83.3 % are exactly unchanged; of the rest, **278 move toward the bar, 302 away** (median
  moved-margin −0.019 nat). The nearest-to-bar uncertified pulsars (3/run, median 2.65 nat below
  the bar) are 75/90 exactly frozen, 8 toward, 7 away. **A reshuffle, not a climb.**

**Mechanism, and it is not subtle.** The loop **starts at the true source** (`ignite2.py:426`,
`theta_cur = G["theta"].copy()`). At truth the M-step gradient of `lnL_marg` is ~0 and the
3-candidate backtracking line search accepts a step only if the objective improves — so it
correctly **stays**. The fixed point is the seed set because the start is already the fixed
point. IGNITE-2 §4.1 flagged exactly this as untested: *"behaviour from a mis-registered but
unscrambled start … is untested here."*

### 1.4 The verdict on question (B), and the disambiguation §10.17 needs

**gain = 1.000 is a MEASURED PROPERTY OF A DEGENERATE STATISTIC, not a physical marginality.**
It is neither of §10.17(1)(i)'s two options — not a strict-`>`-inequality quantisation of a
positive continuous gain (the continuous gain is zero, not positive), and not inertness at a
live ignition line (there is no line here — the loop is at truth). It is a third thing: **the
ratio statistic is ill-posed at |C|=0, and the loop is initialised at its own fixed point.**

*One disambiguation this forces, for §10.17.* §10.17(1)(i) attaches the "exactly 1.000" to the
SURFACE **count** of 1.000 (30 correct certs in 30 realisations) sitting on the strict `>1`
onset bar — a genuine 1/30-lattice quantisation question about the ONSET test. That is a
**different statistic** from the loop **gain** the KINDLE brief names. Both are real; they are
not the same "1.000". The onset-bar quantisation (is the bar `>` or `≥`?) stands as an open
convention question for the onset map; the loop-gain degeneracy is settled here.

### 1.5 The redesign that makes question (A) well-posed (Matt, this session)

Recorded as Matt's call, per the IGNITE-2 §1.3 precedent that a criterion/convention change is
not a fenced agent's to make:

1. **Ratio gain RETIRED.** Everywhere: **additive gain ΔN = |C_fix| − |C_0|** (well-defined at
   zero; +1 *is* an ignition) and **continuous margin flow** (per-pulsar `dlnL − bar` movement,
   toward/away counts). Ratio gain is kept in the banks as a labelled degenerate trail column.
2. **Arm (a) — truth start** = the control, gated bit-identical against the banked sloop
   trajectories (§3, g0).
3. **Arm (b-MAP) — the criterion's own MAP source** (per realisation: sky at the counterpart;
   f, mc, extrinsics from the **Earth-term-only fit**) = the PRIMARY arm. The honest operational
   question — *does the loop, started where a real analysis actually starts, walk toward truth
   and pick up certifications?* The realised per-realisation start offset is **banked as a
   column** and supplies the x-axis. [Phase 2 — needs the Earth-term-only optimiser; §5.]
4. **Arm (b-FIX) — truth + 1 × σ(log10_mc)** along the mc eigendirection = the pre-registered
   control, same seeds as arm (a): *does a KNOWN displacement of KNOWN size get repaired?* One
   magnitude only. (The mc axis, not sky: the loop **freezes sky** and fits f/mc/extrinsics, so
   only an mc-direction offset is one the loop can act on.)
5. **Scrambled null through every arm.** Off-truth start + no signal is precisely where the
   hard-lock cascade was born (B1.3/IGNITE); the soft loop's silence there is a **deliverable,
   not an assumption.** STOP + anatomy on any scrambled detection.
6. **Sanity column:** |C_0| at the start vs at truth — is the start already inside the
   certified-seed basin? If the MAP starts already sit inside, that is a finding about the
   Earth-term fit, not a wasted arm.

**Pre-registered reading:** arm (a) inert everywhere = expected, maps nothing new; **arm (b) is
the experiment** — ΔN and margin flow vs the capability ladder, the ΔN>0-sustained contour
located or bounded. Scope-of-inference line on the verdict.

---

## 2. D-7 — THE `fALL` COLUMN, RE-CUT AGAINST THE VALIDITY GATE — **CLOSED (RECOVERED)**

The last open item of the floor fix (§10.16(g), §10.17 open-2): the "21 cells ignite under
`fALL`, all 5+11" claim (§10.13(f)) stood on the pre-fix Gumbel estimator, because
`surface_floors.npz` banks **no `fALL_zerofrac`** — the validity gate was never applied to the
`fALL` column, and RECUT banked only the nullN offenders.

**The recoverability finding that closes it as RECOVERED, not UNRECOVERABLE-from-disk.**
`surface.py:23` defines `fALL` as the Gumbel over `nullN + nullA + nullL` (N = 200/cell), and
the offender statistic (`surface.py:76`) is a **pure function of three columns** — `dlnL_det`,
`lnK`, `qmax` — that **every** one of the 21,600 `sf_null{N,A,L}_*.npz` banks carries. The
`fALL` offender *vectors* were never banked, but they are **fully recomputable from disk.**
Script `hpc_harbor/kindle/kindle_d7_fall.py` (CPU, no GPU, no new realisations), mirroring
`recut_surface.py`'s gate discipline exactly.

**The gates — no corrected number emitted until an uncorrected one reproduced from the same raw
columns:**

| gate | what it proves | result |
|---|---|---|
| **A-ALL** | recomputed oALL reproduces banked `fALL`, `fALL_mu`, `fALL_beta`, `fALL_sd`, `fALL_n`, `fALL_emp` | **0.000e+00** (all six) |
| **B-ALL** | this scorer at the OLD `fALL` floors reproduces banked `corr_fALL` / `corr_fALL_lo` | **0.000e+00** |
| **readback** | the banked offender vectors reproduce the re-cut `zero_frac` | **0.000e+00** |

**THE ADJUDICATION: the claim SURVIVES the floor fix exactly — 21 → 21, zero verdict changes.**

| | ONSET | MARGINAL | below | by structure |
|---|---|---|---|---|
| PUBLISHED (pre-fix Gumbel `fALL`) | 21 | 2 | 85 | 21 × 5+11 |
| **RE-CUT (criterion-v2.2 `fALL`)** | **21** | **2** | **85** | **21 × 5+11** |

**Mechanism — why the correction cannot bite here.** All 21 igniting cells have `fALL`
zero-fraction **0.00–0.01**, far below the 20 % validity gate, so the Gumbel is **valid** at
every one of them and is left untouched. Only **10 of 108** cells fail the gate (max
zero-fraction 0.71), and **none** was an `fALL` onset candidate — the `fALL` null family
(nullN+nullA+nullL, which includes the *scrambled-source* nullA/nullL that put loud offenders in
the tail) has a far lighter zero point mass than the `fN` family, exactly where the loud 5+11
cells live. **The `fN`-proxy argument (§10.13(f): matched cells at z = 0.00–0.03) was right —
and it is now a re-cut, not a proxy.** The best cells stand: (−12.00, 40, lit, 5+11) at 2.57,
(−12.00, 40, vlbi, 5+11) at 2.33.

**D-7 closes.** The `fALL` offenders are banked with orientation **declared, not implied** (per
§10.16(d)'s learning): `reports/surface_fALL_offenders.npz`, 108 × 200, with an `index` column
and an `orientation` string stating `off_i ↔ index[i] ↔ meta row i`, no transpose, no implied
permutation, columns 0:100 nullN / 100:150 nullA / 150:200 nullL. The re-cut table is
`KINDLE_results/kindle_recut_fALL.npz`.

---

## 3. STAGE 0 — g0: THE LICENSE — **PASSED (0.000e+00)**

The KINDLE loop body (`kindle_loop.py`) reproduces a banked IGNITE-2 soft-loop trajectory
**bit-identically** — and the target was chosen to be the ignition run
`ig_sloop_h1275_T30_lit_g0_n6240011` (n_cert 0→1), so the license and the demonstration are one
object. Job 12519844, dgx03, cpus=8, 970 s, clean GPU.

| | result |
|---|---|
| shared per-iteration + trajectory columns compared | **18** (`iter_dlnL/qmax/mapk`, `dlnL_final`, `qmax_final`, `mapk_final`, `cert_final`, `on_true_final`, `lnK`, `n_true_grid`, `traj_n_cert/wrong/W/sig_mc_dex/lnmarg/step_norm/src_mc_off`, `x_final`) |
| **worst max\|diff\|** | **0.000e+00** (0 nonzero) |

**The noise path.** The T=30 IGNITE geometry has **no canonical `b1_L_gwb` bank** — the RECUT
bank is (3248,3248) for the T=15 ANCHOR array; the T=30 build needs (5336,5336). IGNITE-2 ran
2026-07-12, *before* the bank existed, so its T=30 loops recomputed `L_gwb` via `np.linalg.eigh`
at cpus=8. KINDLE reproduces that path (`force_recompute_lgwb`, provenance line
`RECOMPUTED-UNSAFE fp=9fd547b39b02c705` in the log); the **0.000e+00 bit-identity is itself the
proof the thread count matched** (§10.14(g): the recompute reproduces only at cpus=8). *The
canonical bank did not exist for this geometry to begin with; recompute-at-8 is the convention,
not a fallback.*

**The ignition, read three ways** (from the licensed trajectory):

```
n_cert traj : [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
RATIO gain  : [nan, nan, nan, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   <- NaN ACROSS the ignition
ADDITIVE dN : [0,   0,   1,   1,   1,   1,   1,   1,   1,   1]      <- reads the ignition as +1
```

**This is g1 made concrete on a single trajectory:** the ratio statistic is blind to the exact
event the campaign exists to find (it NaN's the 0→1 step and then reports 1.000 tautologies);
the additive ΔN reads it correctly. **The KINDLE loop body IS the banked loop body**, so every
off-truth number the ladder produces is on the licensed machinery.

## 4. THE OFF-TRUTH PILOT — instrumentation validated; one design finding

At the VLBI onset cell (−13.25, 30, vlbi; fN 7.59), 6 reps × [arm (a) truth, arm (b-FIX)
1×σ(mc)] at matched seeds (the banked sloop convention, so the truth arm reproduces banked
trajectories) + scrambled-truth null. Jobs 12520039 / 12520542, dgx03, cpus=8. **The additive-ΔN,
margin-flow, start-offset, and seed-basin columns all bank correctly — the redesigned
instrumentation works end-to-end on the licensed machinery.** Three findings:

| | arm (a) truth | arm (b-FIX) 1×σ(mc) |
|---|---|---|
| M-step accepted ≥1 step | **2 / 6** | 1 / 6 |
| additive ΔN per real | [**−1**, 0, 0, 0, 0, 0] | [0, 0, 0, 0, 0, 0] |
| \|C_0(start)\| | [1, 1, 0, 0, 0, 0] | [**1**, 0, 0, 0, 0, 0] |
| realised `src_mc_off` (start→final, dex) | — | 0.0000→0.0001 |

1. **The loop is inert at truth — confirmed.** Only 2/6 truth-start reals accept any M-step
   (g1's 25/30-frozen picture, reproduced independently). The one that moved *lost* its cert
   (ΔN = −1) — a floor-margin wobble of a near-bar cert (IGNITE-2 §3.2's one lost-true-cert
   class), not a wrong-fix. **This is the direct evidence for the g1 claim that the loop is at
   rest at its start, not poised at an ignition front.**
2. **THE DESIGN FINDING: 1×σ(mc) is a negligible displacement.** The profile-width σ(log10_mc)
   is ~1e-4 dex, so the fixed offset moves the source `src_mc_off` by only ~1e-4 dex — far below
   anything the loop responds to. It changed a certification in exactly one place (rep with
   `|C_0(start)|` = 0 where truth had 1: the 1σ offset dropped a *floor-margin* cert and the loop
   did **not** walk it back), and moved nothing elsewhere. **So the pre-registered fixed 1σ(mc)
   control is too small to be a meaningful "known displacement" — which is itself the result: the
   meaningful off-truth perturbation is arm (b-MAP)'s Earth-term start with its *natural* spread,
   not a fixed σ-multiple.** Matt's design anticipated exactly this ("the MAP arm's natural
   spread covers the ladder; the fixed offset is the controlled A/B"): the A/B confirms the
   controlled displacement is repairable-or-inert only at the fragile floor margin, and hands the
   dynamic range to the MAP arm. **For the launched ladder, arm (b-FIX) should use a larger
   pre-registered magnitude (2× the seeder tolerance, IGNITE-2 §4.1) — 1σ(mc) is recorded here as
   measured-too-small.**
3. **Scrambled-truth reproduces IGNITE-2's STOP anatomy — and does so bit-for-bit.** 3
   scrambled-truth reals: 2 silent, **1 keeps a wrong certification** (`|C_fix|` = 1, `n_cert`
   1→1). The keeper is **J0437−4715, dlnL = 11.83, qmax = 0.906, Δk = −4** — a small-|Δk|
   noise-lock, **present at iteration 0 and untouched by the loop** (the loop froze after 2
   iters). This is **identical** to IGNITE-2 §3.2's documented scrambled keeper (same pulsar,
   same Δk = −4, same dlnL = 11.83, same qmax = 0.906) — the seed convention matches, so this is
   IGNITE-2's own realisation reproduced through the KINDLE driver: **a second bit-level license,
   for the scrambled path.** The verdict is the same: **STOP, with every failure inherited from
   the seed, none generated by the loop** — the D1 wrong-counterpart hole travels unchanged, as
   stated everywhere it appears. **Design note, banked as a lesson:** scrambled **+ fix_mc**
   (scramble the sky *and* offset mc) is a double perturbation that leaves the covariance
   degenerate and produces an all-non-finite E-step; the `_finite` guard caught it and — with the
   per-realisation `try/except` now in the driver — skips it without killing the shard (the
   blast-radius discipline HPC_SETUP §2 requires). **The scrambled arm is single-perturbation
   (scrambled-truth), exactly as IGNITE-2/CHORUS ran it.**

**Pilot verdict.** The redesign is executable and the instrumentation is sound (all new columns
bank; g0 licenses the signal path and the scrambled path reproduces IGNITE-2 bit-for-bit); the
loop's truth-start inertia is confirmed on fresh reals; the STOP anatomy (seed-static D1
noise-lock) is reproduced exactly; and the pilot *sharpened* the ladder design — the primary
signal is the MAP arm's natural start spread, the fixed-offset control needs a larger magnitude
(1σ(mc) measured too small), and the scrambled arm is single-perturbation. **No ΔN > 0 ignition
was seen in the pilot, but the pilot cell is below onset (truth-start certs 0–1) and the fixed
offset was sub-threshold — the ignition question is the launched ladder's, at the above-onset
Pair-B cells and louder, from MAP starts.** Banked: `KINDLE_results/k_pilot*.npz`,
`kindle_pilot_summary.npz`.

## 5. THE LADDER — pre-registered plan for the launched campaign

**Why launched, not inline.** The g1 redesign is settled and the driver
(`hpc_harbor/kindle/kindle_loop.py`) is built and g0-licensed (§3). But the full capability
ladder is a multi-component GPU campaign that must not be committed until the pilot (§4) shows
the off-truth arm expresses ΔN and the scrambled arm stays silent. Four engineering surfaces
separate the pilot from the ladder, and each is a real build, not a parameter change:

| # | surface | status | note |
|---|---|---|---|
| E1 | **T = 40 problem build** | `IG.build_ignite_problem(40)` exists (takes `T_label`) | Pair B is at T=40; the T=50 rows carry SURFACE's stated 35-yr-extrapolation caveat |
| E2 | **5+11 / 2+14 structure** | lives in `surface.py`, NOT `ignite.py` (which is census 3+13) | the loudness-structure axis must be married to the soft loop — the ladder driver builds the promoted population, not the census |
| E3 | **arm (b-MAP) Earth-term optimiser** | **NOT built** | the primary arm needs a mask=0 (EarthDelay-only) source MAP fit; `B1Split`/`B1Marg` expose the mask (HPC_SETUP G1), so it is buildable, but it is a new optimiser and its own gate |
| E4 | **eccentric cells** | CHORUS `chorus.py`, a different problem build | e=0.5 single / e=0.3 two-member / e=0.7; these run through the CHORUS soft loop, not ignite2's |

**The pre-registered ladder** (ascending the capability axes from the corrected onset cells;
≥12 realisations/cell, ≥5 sky draws among them; additive ΔN + margin flow + σ(log10_mc) + W +
wrong-certs banked per (cell, realisation, iteration); ratio gain kept as the degenerate trail
column):

1. **the two reserved purity-preferred cells** — Pair B, SURFACE §10.13(h), floor-fix-untouched
   (zero-fraction 0.00): **(−12.75, 40, vlbi, 5+11)** fN 122.46, **(−13.00, 40, vlbi, 5+11)**
   fN 44.40.
2. **+loudness** — the same two cells at **h + 0.25** and **h + 0.5 dex**.
3. **+structure** — the **5+11** and matched **3+13** variants at matched loudness (the SURFACE
   promotion axis).
4. **+eccentricity** — CHORUS **e = 0.5 single-member** and **e = 0.3 two-member** (the
   *corrected* switch points, §10.14(b) / RECUT §2.1 — NOT the refuted e=0.3 single), and one
   **e = 0.7** cell (strong seeds — does a big seed set change the ΔN regime?).

**Each cell runs three arms** — (a) truth control [gated vs banked], (b-MAP) primary [E3],
(b-FIX) 1×σ(mc) control — to fixed point or 12 iterations. **Scrambled null through both arms
at ≥2 cells (≥5 reals each): must stay silent; STOP + anatomy on any detection.**

**Cost.** IGNITE-2's banked per-realisation wall: ~750 s first-in-process (evaluator
materialisation), then 118 s (2-iter, step-rejected) to 320 s (full 10-iter). Off-truth arms
**accept steps**, so expect the 320 s class, not 118 s. Order-of-magnitude: ~14 cells × 3 arms ×
≥12 reals × ~300 s ≈ **35–45 GPU-hours** signal, plus scrambled, plus E3's optimiser gate.
A throttled array (`%8`) on the empty `dgx_iacc` lane clears this in a few wall-days.
**Blast-radius discipline:** one realisation per task, skip-on-exist resume, `--cpus-per-task=8`
(the hazard), device-log every task.

**The contour reading, pre-registered.** Arm (a) is expected inert everywhere (maps nothing
new — it confirms the seed set is the fixed point). **Arm (b) is the experiment:** ΔN and margin
flow vs the capability axes; the **ΔN > 0-sustained** crossing located or bounded. The verdict
will read either *"the self-amplifying regime begins at [cell class]"* or *"ΔN remains ≤ 0
everywhere in the box — the array does not self-amplify from an honest start within it, bounded
below by [the margin-flow extrapolation]."* The F15/F11-KINDLE panel is ΔN (discrete) and mean
margin flow (continuous) over the (loudness × structure × eccentricity) ladder, in the explainer
style.

---

## 6. SCOPE OF INFERENCE + THE ADD-LIST

**Scope of inference (mandatory on the verdict).** Everything here is scoped to the **116-pulsar
mock array** (AXIS, 1440 MHz single-frequency; §10.15(a)) and its frozen prior stack
(`best_distances.txt`, per-pulsar `ln K_counted`, the tiers, the geometry ensembles). The g1
diagnosis is a statement about the **banked soft-loop corpus** (70 trajectories, two campaigns,
all started at the true source). The redesigned sweep, when run, will measure the loop's
behaviour from **off-truth starts on this mock** — it does **not** touch a real TOA, and the
"self-amplifying array" verdict will be a property of *this population structure and this prior
stack*, which nature does not have to supply (the SURFACE/CHORUS population-scoping caveat
travels unchanged). **The D1 wrong-counterpart hole is open and travels beside every count**;
the scrambled arm measures, but cannot close, it.

**What is settled irrespective of the launched campaign:**
- **g1** — the ratio-gain "1.000 marginality" is a degenerate statistic + a fixed-point start,
  not a physical bound. Banked: `KINDLE_results/kindle_g1_diagnosis.npz`.
- **D-7** — the `fALL` column is re-cut against the validity gate and the "21 cells, all 5+11"
  claim survives exactly (21→21). The floor fix's last open item is closed. Banked:
  `reports/surface_fALL_offenders.npz`, `KINDLE_results/kindle_recut_fALL.npz`.

### THE ADD-LIST (Matt commits — I never do)

```
# D-7 closure (the last floor-fix item) — CPU, gated 0.0e+00
reports/surface_fALL_offenders.npz          # 108x200 fALL offenders, orientation declared
KINDLE_results/kindle_recut_fALL.npz        # the fALL re-cut table (108 cells)
hpc_harbor/kindle/kindle_d7_fall.py         # the re-cut script (gates A-ALL/B-ALL/readback)

# g1 diagnosis — the marginality answer
KINDLE_results/kindle_g1_diagnosis.npz      # 70-trajectory ratio-gain degeneracy + additive dN

# the KINDLE driver + g0/pilot (GPU) — g0 PASSED 0.0e+00, pilot COMPLETE
hpc_harbor/kindle/kindle_loop.py            # the re-instrumented soft loop (start-injection hook)
hpc_harbor/kindle/kindle.py                 # g0 + pilot driver
hpc_harbor/kindle/kindle_figure.py          # the F15/F11-KINDLE panel generator
hpc_harbor/kindle/kindle_pilot_analysis.py  # the pilot §4 table
hpc_harbor/kindle/k_g0.sbatch
hpc_harbor/kindle/k_pilot.sbatch
KINDLE_results/k_g0_h1275_T30_lit_truths0_g0_n6240011.npz   # the g0-licensed ignition trajectory
KINDLE_results/k_pilot*.npz                  # the off-truth pilot banks (12 signal + 3 scrambled)
KINDLE_results/kindle_pilot_summary.npz      # the pilot §4 summary

# the report + figure
reports/KINDLE_gain_contour.md
reports/KINDLE_gain_contour.png              # the F15/F11-KINDLE panel (g1 degeneracy + pilot)
```

**Doc / deck notes to stage (§10.17 + the deck):**
- **§10.17(1)(i)** — the named "exactly-1.000" question is answered as the **third option**: the
  loop *gain* statistic is degenerate-by-construction (ratio ill-posed at |C|=0) and the loop is
  initialised at its own fixed point. (Distinct from the SURFACE *count* 1.000-on-the-`>1`-bar
  quantisation, which §10.17(1)(i) also raises and which remains an onset-map convention
  question — the two "1.000"s are different statistics.)
- **§10.17 open-2 (D-7)** — mark **CLOSED (RECOVERED)**: the recompute-from-columns route beat
  the "unrecoverable-from-disk" allowance; 21→21, zero verdict changes, gates 0.0e+00.
- **The F8 trajectories caption** — the panel stands; retire the *"measured gain 1.000"* box in
  favour of the count-level statement: *"the certified count holds; it has never multiplied from
  n ≥ 1 in 50 honest truth-start runs — off-truth behaviour is the KINDLE campaign."*
- **The reserved final panel (F15/F11-KINDLE)** — is the ΔN + margin-flow contour, drawn from
  the launched ladder (§5), not the ratio gain.
- **Missing tag** — `criterion-v2.2` was never cut; HEAD `db2075a` is its content. Cut the tag.

**STOP.** Stage 0 delivered (g0 license, g1 diagnosis, g2/D-7 closure); the pilot de-risks the
launched ladder; the ladder is pre-registered and costed. The self-amplifying-array contour is
the launched campaign, awaiting the pilot verdict and the E3 optimiser.
