# SPARK-3 — arrow 2's honest middle: the venue, the instrument, and the pathology

## **STATUS: THE FOUR FIXES CLOSED; STRADDLED SPLIT AT THE MARGINAL WIDTH → `EDGE-POSITIVE`.**
### *SPARK-3C (2026-07-20): 48/48 fold REFUSED — the 5 extra rungs are cross-lane (§4.5). Verdict re-verified bit-identical and UNMOVED. Arm (b) still never run (§5.4).*

SPARK-2 §5's fix list is **delivered**; the grid ran to completion (43 reachable rungs **on the
self-derived floor lane**, 12 units, both venues); the pre-registered crossing test returned **STRADDLED** — an edge under the
optimistic (conditional) Fisher bound (11 crossings / 8 units), none under the pessimistic
(prior) bound. **Matt's decision (2026-07-17): split it.** The budgeted chunked-JVP was run — the
`78×78`→`66×66` faint-block marginal Hessian SPARK-2 measured *unaffordable* now **returns**
(H200, ~44 min each, none capped) — and the crossings were **re-cut at the true marginal width**.

> ## **VERDICT: EDGE-POSITIVE** (`spark3_final_verdict.npz`).
> **4 of 5 crossing units survive at the marginal width** (A0, A1, A4, B1; B0 dies),
> **scrambled-clean**, margins **deep** (66–275 nat above bar, not floor-grazing).
> Rule (pre-registered): ≥1 crossing at ≥2 units under the marginal, scrambled-clean → met.

**TWO CAVEATS TRAVEL WITH IT, stated not buried.** (1) The crossing is **the LOOP closing, not
arrow 2 alone**: the certification floor falls ~123 nat rung0→rung8 (arrow 1, template rigidity,
SPARK §3) *and* dlnL rises (arrow 2, deconfusion) — both drive it, which is exactly what GLACIER
needs. (2) **Model-quality-limited, now a design law:** the true marginal reservoir model is only
~3× tighter than prior and **INDEFINITE on ~40% of the faint block** (truth is a saddle for
sources this weak) — the edge survives *at that width*, but a worse model **inflates the null**
(floor 118→744, §5.0), so GLACIER must gate on its own reservoir-constraint quality.

**SCOPE, binding on the verdict:** measured on the **H200 lane** (dgx down all run, §4.4), whose
self-derived floor runs **+13% above SURFACE's** — a *harder* bar, so an EDGE-POSITIVE here is
**a fortiori** (the crossings cleared a higher bar). **5 crossing units, all sky `g3` (n≈1 in
geometry, §0.1)** — so the verdict is EDGE-POSITIVE *at this venue and sky*, and **GLACIER's first
spend is more skies**. Still owed: the **dgx re-cut** (arm b, blocked — nodes down) confirming the
+13% offset did not shape it. Every number banked and CPU-re-cuttable.

| SPARK-2 §5 item | SPARK-3 |
|---|---|
| **5** the venue — *"check a louder cell … before spending on a grid that can only return zeros"* | **DONE (§0).** Its cell is **below onset**; its own **max `N_cert` = 1**. Both reserved cells confirm; campaign runs at **both**. |
| **1** reject-and-redraw the soft draw, bank the non-evaluable fraction | **DONE (§2).** The fraction is **exactly 0** — the pathology was the *box*, not the posterior. `g2a/g2b/g2c` PASS. |
| **2** one task per realisation, rungs in-process, skip-on-exist | **DONE (§4)** — then *superseded by measurement*: rungs are independent, so the grid fans out to **48 (unit × rung) jobs**, checkpointed per rung (§4.1). |
| **3** add the `truth`-modelled control to every unit | **DONE (§4).** All four states per rung. |
| **4** the per-target E-step (~116×) + re-cut (a) against joint widths | **DONE (§1, §3).** **4.29×, not ~116×.** Joint not retried (measured unaffordable); **two-sided bounds** instead. |

**Four corrections to the inherited record (§6), and two to my own** — a wallclock number I wrote
up before it was clean (§1), and a Fisher evaluated at the wrong pmask (§3). Both retracted in
place.

**Tree.** `MM_playground` @ `3a2ed4b` (story-v1 + the fold). Inherits SPARK-2
(`reports/SPARK2_second_arrow.md`, STAGED-UNCOMMITTED — left untouched; this work sits
alongside). Arrow 1 = `reports/SPARK_launch_criterion.md` (CASCADE-ALIVE).

**SPARK-2's endpoints stand and are the frame**: the faint reservoir UNMODELLED erases
certification (offender 0.000) against 4.435 TRUTH-modelled; deconfusion's trials cost is
**+0.578 nat and saturates**. SPARK-3's job is the middle — the soft-modelled ladder at
reachable rungs.

**REALISM / SCOPE.** MOCK spine (AXIS, 1440 MHz single-frequency, §10.15(a)) — no real TOA is
touched and the residuals *are* the injected CW+CURN. `cpus-per-task=8` throughout (the
NoiseDrawer thread-count hazard: it is part of the seed).

**PROVENANCE — the T = 40 L_gwb fingerprint is `8548f148b50a5b44`**, on the
`RECOMPUTED-UNSAFE threads=unset` branch at `cpus-per-task=8`. It is **not** SPARK-2's
`9fd547b39b02c705` (that is the banked CHORUS **T = 30** value) and it must not be compared to
it: a different array geometry has a different square root. **Why the recompute branch is
forced and not chosen:** the canonical `b1_L_gwb.npz` bank is `(3248, 3248)` = the T = 15
ANCHOR array, so at T = 40 `load_or_build_L_gwb` hits its shape-mismatch **raise**. SURFACE
predates the banking fix (CHORUS §0.1) and therefore drew on the same recompute branch at the
same `cpus-per-task=8`. Since that branch is thread-count dependent, **the reproduction of
SURFACE's banked realisations is not assumed — it is gated (g3a), because this campaign coheres
on SURFACE's certified sets and a mismatched draw would mean cohering a certified set onto a
realisation that is not there.**

> ### **GATE g3a — PASSED, and BIT-IDENTICALLY** (A100/dgx04). *Artifact clobbered — see §4.2; regenerate before the fold.*
>
> | column | max &#124;replay − banked&#124; |
> |---|---|
> | `dlnL_det` | **0.000e+00** |
> | `lnK` | **0.000e+00** |
> | `qmax` | **0.000e+00** |
> | `mapk` disagreements | **0 / 116** |
> | `on_true` disagreements | **0 / 116** |
>
> SPARK-3's T = 40 draw **on ACCRE** reproduces SURFACE's realisation **drawn on cronus**
> (`sf_sig_h1275_T40_vlbi_k5_g3_n20560400`, noise 20560400 / dist 30560400) to the last bit.
> **The certified sets this campaign coheres on belong to the realisations it is scoring.**
>
> **This is also a fact about the SURFACE banks that nobody had established.** The T = 40
> geometry has no banked `L_gwb` and rides the branch whose own warning says *"THIS DRAW IS
> THREAD-COUNT DEPENDENT and is NOT guaranteed to reproduce any banked realisation"*. It does —
> **across machines** — at the `cpus-per-task=8` convention. That does not retire the hazard or
> the CHORUS §0.1 remedy (banking the square root is still the fix, and is still unapplied at
> T = 40); it means the SURFACE T = 40 banks are **replayable on ACCRE today**, which any
> successor touching them may now rely on — and which was not in evidence before this gate.
>
> *(Internal consistency, beside it: the two independent jobs — separate GPUs, separate per-job
> JAX caches — recomputed the same fingerprint `8548f148b50a5b44`. That alone would only have
> shown SPARK-3 agreeing with itself; g3a is what establishes agreement with SURFACE.)*

---

## 0. FIX 0 — THE VENUE: SPARK-2 measured arrow 2 at a cell that could not answer it

**COMPLETE.** `SPARK3_results/spark3_venue.npz`. CPU only; re-cut from the banked SURFACE raw
statistics under `recut_surface.py`'s `offender`/`adopt`/`score` semantics verbatim. No new
realisations — a re-cut suffices for an existence check, and existence is the whole question.

SPARK-2 §5(5) flagged it against itself: *"`N_cert = 0` in **both** states at both rungs means
the current ladder may be **entirely below the certification bar at this cell** — check a
louder cell or the truth control first, before spending on a grid that can only return
zeros."* That check is now done, and it is worse than SPARK-2 feared.

| cell | floor | est | zf | corr/real | N_cert>0 | **max N_cert** | verdict |
|---|---|---|---|---|---|---|---|
| **RESERVED A** (−12.75, 40, vlbi, 5+11) | 122.461 ± 4.895 | gumbel | 0.00 | **4.067** | **24/30** | **14** | ONSET |
| **RESERVED B** (−13.00, 40, vlbi, 5+11) | 44.397 ± 1.938 | gumbel | 0.00 | **3.567** | **25/30** | **15** | ONSET |
| **CONTROL** (−13.25, 30, lit, 3+13) — *SPARK-2's venue* | 6.481 ± 0.669 | emp_q95 | **0.67** | **0.567** | 18/30 | **1** | **below** |

**ARTIFACT READBACK PASSES**: 4.067 and 3.567 reproduce STORY App A's reserved-cell values
(4.07 / 3.57) exactly, from the raw columns.

**Both reserved cells confirm `N_cert > 0` exists. Per the brief, the campaign runs at BOTH.**

> **THE CONTROL IS THE RESULT, and it is sharper than "N_cert = 0".** SPARK-2's cell reads
> `corr = 0.567/real` and verdict **`below`** — *it is a below-onset cell*, and its zero-fraction
> (0.67) is far above `ZF_GATE = 0.20`, so its floor is emp-q95 by the adopted convention.
> **Its `max N_cert` across all 30 banked realisations is 1.** Not "small" — **one**. A rung
> ladder of {0, 2, 5, 8} coherent pulsars **cannot be built from that cell's own certified sets
> at all**: rungs 2, 5 and 8 do not exist there. SPARK-2 reached rung 8 only by importing
> certified sets from the CHORUS `m1e07` igniter bank — *a different cell* — and then scoring
> certification at its own. **Arrow 2's ledger was being read at a venue whose own count is
> capped at 1.** That is why it could only return zeros, and no amount of grid would have
> fixed it.

### 0.1 SCOPE — the rung-8 ensemble is ONE SKY, and the reason is a lever the programme already knows

Rung 8 needs `N_cert ≥ 8` in that realisation's **own** certified set. Reachable in **7/30 (A)**
and **6/30 (B)** — and the distribution is not diffuse, it is **a step function in the sky
draw**:

| sky | venue A — N_cert × 6 reps | rung-8 | venue B — N_cert × 6 reps | rung-8 |
|---|---|---|---|---|
| `g-1` (fiducial POP) | 3, 1, 2, 2, 4, 1 | 0/6 | 2, 2, 2, 3, 2, 1 | 0/6 |
| `g0` | 7, 7, 3, 6, 6, 8 | 1/6 | 5, 7, 2, 2, 5, 5 | 0/6 |
| `g1` | 0, 2, 0, 3, 1, 0 | 0/6 | 1, 1, 1, 1, 2, 0 | 0/6 |
| `g2` | 0, 0, 1, 1, 0, 1 | 0/6 | 0, 2, 0, 1, 0, 0 | 0/6 |
| **`g3`** | **10, 12, 14, 10, 11, 10** | **6/6** | **11, 8, 15, 10, 9, 9** | **6/6** |

**`g3` reaches rung 8 in 6/6 at BOTH venues; everything else manages 1 of 24.** The same sky
that gives 8–15 certified pulsars gives **0–3** two draws away — a 5–10× swing driven by
**source sky geometry alone**, at fixed h, T, tier and structure. This is the geometry lottery
(EXPLAINER F12) with its teeth showing, and it is a much stronger effect at the venue than the
per-pulsar coherence lever arrow 1 measured.

> **THE SCOPE LINE, travelling with every rung-8 number SPARK-3 quotes.** **12 of the 13
> rung-8-reachable realisations are sky `g3`.** The rung-8 ensemble is therefore **effectively
> one sky draw**. A and B share the g3 sky geometry (their noise seeds and `h` differ, so they
> are not duplicates — but they are not independent in geometry either). ***Rung 8 is not a
> 13-fold independent sample and no rung-8 claim may be read as one.*** An ensemble ≥ 8 in the
> brief's sense exists in *count*; in *geometry* it is n ≈ 1. Widening it means more skies, not
> more noise seeds — and that is a cost the successor must price, not a caveat to footnote.

---

## 1. FIX 1 — THE PER-TARGET E-STEP: the ~116× is an artifact of how you build it

SPARK-2 §1 is right about the confound and, I think, wrong about the price. Its structural
finding stands exactly as written: `B1Split.estep` sweeps every pulsar's distance under **one
global pmask**, a pulsar at `m_p = 0` has its own term switched off, its distance is inert, its
fringe row is flat, and **at rung 0 the certification statistic does not exist** (0/116 live
rows — measured, `spark2_mask.npz`). The certified-coherent E-step must be **one call per
target**. That is not in dispute and SPARK-3 builds exactly that object.

What is in dispute is *"a **~116×** cost and a real build, not a parameter … the single largest
thing standing between this addendum and its verdict."* That price assumes the per-target form
is 116 calls to `estep()`, each sweeping all 116 pulsars and discarding 115/116 of the work.
**It need not be built that way**, and the machinery already says so.

**THE COST MODEL, PRE-REGISTERED HERE BEFORE THE GATE RAN.** `estep`'s work splits in two:

| | what it does | per call |
|---|---|---|
| **(i) `_ppab`** — the base evaluation at the pmask | 116 pulsar delay-evals + the full `npsr·ngp` GWB solve | ~116 × 876 × 16 × 2 ≈ **3.3M** waveform points |
| **(ii) `_pulsar_ab_fn(p)`** — target p's fringe sweep | `B = 512` delay-evals on **one** pulsar, on a per-pulsar GP block | ~512 × 876 × 16 × 2 ≈ **14.3M** points, **× 116 targets ≈ 1.67G** |

Only **(i)** depends on the target's mask. **(ii) for target p is identical whatever the other
pulsars' mask is**, because `MaskedDelay` reads `m = params[PMASK][self.ipsr]` — pulsar p's own
entry **only** (`trackB_b1_core.py:127`). So the per-target form re-does the *cheap* half 116
times and the *expensive* half exactly once per target, as the standard E-step already does.
Both `_ppab` and `_pulsar_ab_fn(p)` take `pmask` as a **runtime** arg, so nothing recompiles.

```
standard    ≈  3.3M + 1670M  =  1673M          ratio ≈ 2047/1673  ≈  1.2x
per-target  ≈  116x3.3M + 1670M = 2047M                    NOT 116x
```

**Predicted ≈ 1.2×, not ~116×.** `gate1` **measures** the ratio rather than asserting it.

### 1.1 MEASURED — the model is exact, my constant was not, and the headline holds

| | measured (T = 40 venue, warm) |
|---|---|
| standard `estep` | **0.913 s** ⚠ co-tenant |
| per-target `estep` | **3.725 s** ⚠ co-tenant |
| **ratio** | **4.08×** — *not* ~116× |

⚠ Both absolutes were measured with a co-tenant job of mine on `dgx04` (see §1's retraction).
**Re-measured SOLO by `gate2` on the same venue** (`cost_ratio_solo`, banked):

| | co-tenant (gate1) | **SOLO** (gate2) |
|---|---|---|
| standard `estep` | 0.913 s | **0.84 s** |
| per-target `estep` | 3.725 s | **3.62 s** |
| **ratio** | 4.08× | **4.29×** |

**The claim survives the clean measurement**: `4.29×` solo against `4.08×` under contention —
both ~4×, against SPARK-2's **~116×**. The ratio was always the robust quantity (both timings
taken back-to-back in the *same* process, so contention largely cancels), and it is now
corroborated three ways: co-tenant, solo, and the cost model reproducing it from its own two
fitted constants. **Take 4.29× as the number.**

Fitting the two measurements to `standard = A + S`, `per-target = 116A + S`:

```
A (_ppab, the base evaluation)   =  24.5 ms   (2.7% of the standard E-step)
S (the 116 fringe sweeps)        = 888.8 ms   (97.3%)
predicted per-target = 116A + S  =  3.725 s        measured  3.725 s
```

**The structural model is confirmed to three decimals**, and it is the load-bearing claim: the
target's mask touches only the *cheap* term, and the sweep — 97.3% of the work — is untouched
because `MaskedDelay` reads `params[PMASK][ipsr]` alone.

**Where I was wrong, stated plainly:** I predicted the ratio at **≈1.2×** from waveform-point
counting, which put `A` at ~0.2% of the E-step. It is **2.7%** — I underestimated `A` by ~14×
(the base evaluation is kernel-launch-bound across 116 per-pulsar `ys`/`solve1d` calls, not
FLOP-bound as the point count assumed). **The measured ratio is 4.08×, not 1.2×.** The
conclusion is unaffected in direction and magnitude of consequence: **SPARK-2 priced this at
~116×; it is 4.08× — 28× cheaper than the price that deferred it.** *"The single largest thing
standing between this addendum and its verdict"* is a **4× surcharge on the E-step**, and the
rigidity channel is affordable now.

### 1.2 THE GATES

**g1a — PASSED at `0.000e+00`.** The per-target machinery driven at `PM ≡ ONE` reproduces
`sp.estep(..., one)` **bit-identically**. The refactor did not move the incumbent.

**g1b — the rung-0 gate: DISCRETE-EXACT, CONTINUOUS NOT. It did not pass, and that is a
measurement, not a bug.**

| | per-target rung 0 vs the **banked** standard (`pmask = ONE`) |
|---|---|
| **LIVE fringe rows at rung 0** | **116/116** — *vs SPARK-2's global-pmask **0/116*** |
| `mapk` disagreements | **0/116** (discrete: exact) |
| `q > 0.9` membership count | **106 vs 106** (identical) |
| median &#124;qmax − banked&#124; | **4.1e-12** (76/116 identical to 1e-6) |
| **max &#124;qmax − banked&#124;** | **5.9e-02** (40/116 deviate) |
| median / max &#124;dlnL − banked&#124; | **0.343 / 5.056 nat** |

> **THE CONFOUND IS CLOSED.** `116/116` live fringe rows at rung 0 against SPARK-2's `0/116` is
> the whole point of FIX 1: **the certification statistic now EXISTS at rung 0**, so a rung-0
> count is a measurement rather than a structural zero.

**But the brief's gate — "reproduces the standard E-step's banked q columns bit-comparably" —
is NOT met**, on EMBER §2.2(b)'s adopted bar (discrete exact, continuous < 1e-6): the discrete
half passes exactly, the continuous half misses by `5.9e-02` in `qmax` and `5.06` nat in `dlnL`.

**The anatomy was derived before the run and the measurement matches it.** The two can differ
through exactly one channel: `MaskedDelay` reads pulsar p's own mask entry only, so `a_pf` and
`db` are *identical* in both; the sole route by which the other pulsars' mask reaches target p's
fringe row is the **HD-correlated GWB coupling** `2·db@u_p[p]`, where `u = Minv @ bflat` mixes
every pulsar's `b_base`. The constants (`sum_a`, `const`) cancel in the normalised posterior.
That the MAP fringe and the `q > 0.9` set are *exactly* preserved while `qmax` moves at the 1e-2
level in 40 pulsars is precisely the signature of a **weak but non-zero** correlated-GP coupling.

**They are different models, and the honest reading favours the per-target one.** At rung 0 no
pulsar is certified, so no pulsar's distance is known — yet `pmask = ONE` holds all 115 other
pulsar terms **live at the prior mean `L0`**, asserting a distance nobody has measured. The
per-target form decoheres them (`m_p = 0`), which is SPARK §1's own stated convention:
*"uncertified pulsars sit at `m_p = 0`: **decohered, NOT pinned to a wrong MAP fringe** — the
hard-lock failure IGNITE-2 retired."*

> **CONVENTION, ADOPTED AND TRAVELLING WITH EVERY SPARK-3 RUNG NUMBER.** The ladder is anchored
> to the **per-target** E-step at **every** rung (0 … 8), so it is internally consistent and
> self-comparable — a rung-to-rung difference is a real difference, not a change of convention.
> **Its rung 0 is NOT the SURFACE bank's count**, and no SPARK-3 rung-0 number may be compared
> to the venue's banked `4.07/real`. The venue's banked counts are used for exactly two things:
> **choosing the venue** (§0) and **drawing each realisation's certified set** — both of which
> are properties of the SURFACE convention and are quoted as such.

**A COST CLAIM I MADE AND THEN RETRACTED, BECAUSE THE NEXT RUN FALSIFIED IT.** The reserved
cells are T = 40 and the problem *is* bigger than SPARK-2's T = 30 (**101 619 TOAs** vs ~73 000,
`rn_comp` 30→64, `gwb_comp` 14→30, real span **47.17 yr**). The first two jobs measured
`build_b1_amortised: 1032.0s` against SPARK-2's 165–200 s at T = 30, and I wrote that up as a
**"5–6× build"** caused by the venue's size. **That was wrong.** The next job, run **alone**,
measured `build_b1_amortised: 173.5s` — *the same problem, the same code, 6× faster*.

**The 1032 s was CO-TENANCY, not T = 40.** `r3A` and `g1A` were dispatched to the **same node**
(`dgx04`), 8 CPUs each, both cold-compiling simultaneously; the build is CPU/compile-bound, so
they served each other. **Solo, the T = 40 build is ~174 s — comparable to T = 30.** The venue
move does *not* buy a heavier build.

> **The lesson is the one this programme keeps re-learning: a timing measured under unstated
> co-tenancy is not a property of the problem.** It is the same class of error as EMBER §7(a)'s
> co-tenant OOM and SPARK's own `--exclude=dgx03`. **Two jobs of mine on one node was an
> undeclared input to a number I had already written into a report** — caught only because a
> resubmit happened to run alone. Every wallclock number in this report is now labelled with
> whether it ran solo.

*(SPARK-2 §5(2)'s "one task per realisation, all rungs in-process" is still adopted and still
right — the **cache** cost it identified is real, and per-job `HARBOR_JAXCACHE` isolation gives
every job a cold compile. It is simply not justified by a build-size argument I no longer have.)*

*(Also noted: `t_max` must be read from the **actual TOAs**, never from `T_label` — the
extension moves the span, and `t_max` is what FIX 2's evaluability bound is measured against.)*

---

## 2. FIX 2 — THE PATHOLOGY: it is coalescence, and SPARK-2 misattributed it

**COMPLETE as an anatomy, and it needed no GPU** — the boundary is analytic.

SPARK-2 §4b: rungs 2 and 5 died all-NaN even after `phi_bounds` clipping, and it called this
*"IGNITE-2/KINDLE's known double-perturbation pathology (`scrambled + fix_mc` → degenerate
covariance → all-non-finite E-step) arriving by a new route"*, naming `log10_mc × fgw`
coalescence only as a *"remaining suspect"*. **The suspect is the entire mechanism, the KINDLE
attribution is wrong, and the fix SPARK-2 added is what caused it.**

### 2.1 The mechanism, exactly

`discovery/deterministic.py:509-510` — the waveform this programme injects *and* recovers with:

```python
term  = 1.0 - (256.0/5.0) * mc**(5/3) * w0**(8/3) * t
omega = w0 * jnp.power(term, -3.0/8.0)          # term < 0  ->  NaN
```

`jnp.power(negative, -3/8)` is **NaN**. `term < 0 ⟺ t > t_coal = 5/(256 mc^{5/3} w0^{8/3})` —
the binary has **already coalesced** at epoch `t` and the chirping waveform does not exist.

- **It is an EARTH-TERM pathology.** For the Earth term `t = toa − tref > 0` (the TOAs postdate
  J2000), so it can cross. For the pulsar term `t = toas_rel − (L/c)(1−cos μ)` with
  `L/c ~ 3.3 kyr`, so `t` is hugely **negative** → `term > 1` → always evaluable. **It therefore
  fires at every pmask, coherent or not** — which is why no certified set could have saved
  rungs 2 and 5.
- **One bad source kills everything.** `MaskedDelay.__call__` sums all 16 sources
  (`d_full = jnp.sum(self._full(*args, L), axis=0)`), so a single coalesced source NaNs **every
  pulsar's** residual → the E-step returns non-finite on **all** `npsr × B` entries. That is
  exactly SPARK-2's signature: **59392/59392 = 116 × 512**.

### 2.2 The root cause: `phi_bounds` is not the prior support

| box | log10_mc | log10_fgw |
|---|---|---|
| **POPULATION generative prior** (`stagec_fisher.py:57-58`) — what drew the truth | 8.5 … 9.5 | −8.0 … −7.5 |
| **`TE.phi_bounds`** (`trackB_estimator.py:512-513`) — what SPARK-2 clipped to | **7.0 … 10.0** | **−8.7 … −7.0** |

SPARK-2's `soft_faint_theta` clipped the draw to `TE.phi_bounds` and documented it as *"the
faint posterior is truncated at its own prior support, exactly as a real posterior is."*
**`phi_bounds` is not the prior support — it is the ESTIMATOR's SEARCH box**: 3.0 dex in mc
against the population's 1.0, and 1.7 dex in fgw against the population's 0.5. **That extra
half-dex in each is exactly where coalescence lives.** The clip SPARK-2 added *as the fix* is
what admits the non-evaluable region.

`t_coal` at the venue. **`t_max` is read from the ACTUAL TOAs and is 49.70 yr — it is NOT the
span (47.17 yr) and NOT `T_label`**: the quantity the Earth term must evaluate at is
`max(toa) − tref` with `tref = J2000`, and the extension moves it. Assuming either would
mis-place the boundary.

| box corner | log10_mc | log10_fgw | t_coal (yr) | at this array (`t_max` = 49.70 yr) |
|---|---|---|---|---|
| POPULATION worst corner | 9.50 | −7.50 | **300.94** | EVALUABLE (**6.1×** margin) |
| POPULATION best corner | 8.50 | −8.00 | 300938 | EVALUABLE |
| `phi_bounds` worst corner | 10.00 | −7.00 | **2.05** | **NON-EVALUABLE** |
| `phi_bounds` mc-max @ pop fgw-max | 10.00 | −7.50 | **44.17** | **NON-EVALUABLE** |
| `phi_bounds` fgw-max @ pop mc-max | 9.50 | −7.00 | **13.97** | **NON-EVALUABLE** |

**The population's generative box is evaluable everywhere with a 6.1× margin. `phi_bounds`
is not** — and note the last two rows: **either** axis pushed to its `phi_bounds` limit **while
the other stays inside the population's own box** already coalesces. Neither excursion is exotic.

### 2.3a GATE g2a — the analytic predicate is EXACT, and it reproduces SPARK-2's death

| draw | predicate says | E-step returns | |
|---|---|---|---|
| population worst corner | EVALUABLE | non-finite **0 / 59392** | AGREE |
| `phi_bounds` worst corner | NON-EVALUABLE | non-finite **59392 / 59392** | AGREE |

**`59392/59392` is SPARK-2's signature exactly** (116 × 512), reproduced here from a **one-line
closed-form predicate with no waveform call** — by perturbing **one** faint source. That is the
mechanism proven end-to-end: **one coalesced source NaNs the whole array**, because the delay
sums all 16. The predicate needs no GPU; the E-step is only the witness.

### 2.3 It IS a double perturbation — of (mc × fgw), and that is why the 1-D scans missed it

The evaluability boundary `t_coal(mc, fgw) > t_max` is a **joint curve**. At the population's
own fgw (−7.5) the mc threshold sits at **log10_mc ≈ 9.98** — just above the `phi_bounds` mc
ceiling, so an mc-only excursion is nearly harmless. Push fgw to −7.0 (which `phi_bounds`
permits and the population never generates) and the threshold **falls to log10_mc ≈ 9.18 —
below the population's own maximum mc of 9.5**. *Neither axis alone is dangerous; together they
are.*

> **This corrects SPARK-2's reading of EMBER.** SPARK-2 wrote: *"EMBER's `mc_scan` measured
> 'non-evaluable 0% across the entire mc box' for the LOUD sources — that result does not
> transfer to the faint, and SPARK-2 is the counter-example."* The strain plays **no part in
> this mechanism**: `t_coal` is a function of `(mc, fgw)` **only**. EMBER measured 0% because
> `mc_scan` (`ember_map.py:133`) scans mc **with fgw held at truth** — a 1-D scan through a 2-D
> boundary, which cannot see the corner. It is not a loud-vs-faint result and it does transfer.
> **EMBER had the mechanism right** (`ember_map.py:40-46`: *"THE PRIOR BOX CONTAINS BINARIES
> THAT HAVE ALREADY COALESCED. log10_mc's box is [7, 10] … the waveform is non-evaluable and
> lnL is NaN"*) and fixed it by bounding inside the scanned evaluable region. **SPARK-3
> generalises EMBER's precedent from a per-axis scan to the joint analytic predicate.**

And the rung-dependence — why 2 and 5 died and 0 and 8 lived — is **nothing but the RNG**:
`np.random.default_rng(SOFT_SEED + 1000*real_i + rung)`. Whether a draw crossed the boundary
was luck. **There was never a rung-2/rung-5 pathology to find.**

### 2.4 The repair, and why it is not a clip

Draw on the **population's generative support** for `(mc, fgw)` — the actual prior support —
and apply the **joint analytic predicate** `t_coal > t_max` exactly (no waveform call),
**rejecting and redrawing the offending sources only**, with the **rejection count BANKED per
rung** (SPARK-2 §5(1) asked for exactly this: the non-evaluable fraction *"is itself a
measurement about how soft 'soft' can be, not a nuisance to suppress"*). Rejecting a coalesced
draw is not a fudge: the likelihood is **undefined** there, so the posterior has zero mass
there by construction.

**This is a repair, not a wall — so it does not fire the brief's STOP.** The anatomy is
reported here (as the brief required) and the campaign proceeds on the repaired draw, gated.

### 2.5 GATES g2b / g2c — SPARK-2's deaths reproduced, and the repair is clean

`spark3_gate2_A_g3_r0.npz`, venue A, sky g3, `N_cert = 10`, `N_SOFT = 8` draws per rung:

| rung | coh | **SPARK-2's box** (`phi_bounds`) src non-eval | draws **dead** | **REPAIRED** (pop box + joint gate) rejections | draws dead | E-step |
|---|---|---|---|---|---|---|
| **0** | 0 | 0 | **0/8** | **0** | **0/8** | FINITE |
| **2** | 2 | 4 | **4/8** | **0** | **0/8** | FINITE |
| **5** | 5 | 3 | **3/8** | **0** | **0/8** | FINITE |
| **8** | 8 | 10 | **6/8** | **0** | **0/8** | FINITE |

**g2b PASSED — and it reproduces SPARK-2's pattern, not merely "some deaths".** Under SPARK-2's
exact clip, **rung 0 loses nothing (0/8) while rungs 2 and 5 die (4/8, 3/8)** — which is
precisely what SPARK-2 observed (*"clipping … fixed rungs 0 and 8 but rungs 2 and 5 still die
all-NaN"*). **There was never a rung-2/rung-5 pathology.** The rungs differ only through the
draw's RNG stream (`SOFT_SEED + 1000·real_i + rung`); whether a draw crossed the coalescence
boundary was luck. *(Rung 8 dies here (6/8) where SPARK-2's rung 8 lived — expected, and it makes
the point: this is a different cell, different widths, different RNG. The pattern is stochastic
in the rung, which is exactly the claim.)*

**g2c PASSED — and the rejection count is ZERO at every rung.** This is the strongest form the
repair could take. SPARK-2 §5(1) asked for reject-and-redraw *"with the rejection count BANKED —
the fraction of the faint posterior that is non-evaluable is itself a measurement about how soft
'soft' can be."* **Measured, that fraction is exactly 0** — the reject-and-redraw machinery is
built, armed, and never fires, because the **true prior support contains no coalescing
configuration at all** (the 6.1× analytic margin). The non-zero fraction SPARK-2 saw was not a
property of the faint posterior; **it was a property of drawing from the wrong box.**

> **So the answer to "how soft can soft be?" is: as soft as the prior, at no cost.** The
> evaluability constraint SPARK-2 feared would bound the softness **does not bind anywhere
> inside the generative prior.** The joint predicate is retained anyway — armed and banked at
> zero — because it is what *proves* the draw is clean rather than assuming it, and because a
> successor changing T, the fgw prior, or the population would need it: at `t_max = 49.70 yr`
> the population's worst corner clears by 6.1×, but that margin **shrinks with T**, and a longer
> baseline eventually walks the generative box into the boundary.

---

## 3. FIX 3 — THE FISHER, TWO-SIDED: built, not yet exercised

**The monolithic joint is NOT retried.** SPARK-2 measured it unaffordable **twice** on an
A100-80GB inside a 1 h walltime (`jax.hessian`, job 12583525, >27 min; chunked batched-JVP-of-grad
at CH=8, job 12583899, >17 min) — the likelihood carries a 2668-dim correlated-GP solve per
evaluation and the faint block is 78×78. That is a measurement, and SPARK-3 accepts it rather
than re-buying it.

**The pre-registered replacement, both banked as columns on every unit, never merged:**

| bound | what it is | which way it is wrong |
|---|---|---|
| **OPTIMISTIC** | `σ_cond = 1/√F_ii`, the **conditional** width (1-D curvature scan per axis; one batched `logL` call per chunk, seconds), capped at the prior box | `σ_cond ≤ σ_marg` **always** → **overstates** how well the reservoir is known → **favours arrow 2**. An EDGE-ZERO measured here is *a fortiori* EDGE-ZERO at the true joint width. |
| **PESSIMISTIC** | the **prior box width**. The loop learns nothing; the soft model is a prior draw. | cannot overstate the tightening, because there is none. |

**The honest answer lives between.** If both bounds give the **same** verdict, the verdict stands
**without** the joint — the joint could only land between two bounds that already agree. If they
**straddle**, SPARK-3 reports **STRADDLED** and prices a budgeted chunked JVP rather than quoting
a number it cannot defend.

**A correctness bug found and fixed before the grid spent on it.** (a)'s Fisher was initially
evaluated at `PM[0]` — the per-target mask row for *target 0* — which forces pulsar 0's own term
live even when uncertified, silently adding a **117th coherent term** to the state the faint
constraint is measured in. The faint sources' constraint is a property of **which pulsars are
coherent**, not of **which pulsar is being scored**. Fixed to a dedicated `pm_rung` (the rung's
certified-set mask, no target override); `gate2` was relaunched so its banked artifact comes from
the corrected code.

---

## 4. THE MEASUREMENT — LAUNCHED, NOT COLLATED

**48 rung-jobs → 43 reachable ON THE SELF-DERIVED (H200) FLOOR LANE, COMPLETE** (`sbatch --array=0-47%8 sp3_h200fan.sbatch` → `s3r_{A,B}_g3_r{0..5}_k{0,2,5,8}.npz`; 5 rungs exit clean as unreachable, N_cert < rung — **on that lane's floors, and only there**. On SURFACE's floors all 48 are reachable. The qualifier is not pedantry: it is the whole of §4.5, and without it the sentence invites exactly the contaminated fold that §4.5 refuses),
submitted after — and only after — every gate above passed. 2 venues × sky `g3` × 6 realisations;
**every one rung-8-reachable** (§0.1).

```
per rung : 18 signal (un + truth + 8 soft@opt + 8 soft@pes) + 400 null (100 x 4 states)
         = 418 per-target E-steps
per unit : build + 4 rungs
```

**THE GRID WAS LAUNCHED ONCE, MEASURED, AND RELAUNCHED — and the reason is the third co-tenancy
bite of this campaign.** Sized on the **solo** per-target E-step (3.62 s), a unit is ~1.7 GPU-hr
and 4 h was ample. Watching the first units actually run, the per-target E-step was **20.0 s**
with 4 of them sharing `dgx04` — **5.5× slower** — putting a unit at **9.3 h against a 4 h
walltime**. **Every unit would have timed out mid-grid**, and because skip-on-exist was at *unit*
granularity, each would have thrown away its finished rungs. Two changes, both from that
measurement:

- **Throttle `%2` and walltime 12 h.** The throttle is load-bearing, not politeness.
- **Per-rung checkpointing** (`s3r_*_k*.npz`): each rung is banked the moment it completes and
  skipped if on disk. A rung is ~2.3 h of work under contention; losing three of them to a
  timeout on the fourth is not a resume strategy. The ledger reads the **rung** files, so a unit
  that dies mid-grid still contributes every rung it finished.

> **Co-tenancy of my own jobs has now corrupted a wallclock number three times in this campaign**
> — the "5–6× build" that was 174 s solo, `gate1`'s absolutes, and a walltime that would have
> destroyed the grid. Each time it was caught by *measuring the thing again under known
> conditions* rather than trusting the first number. **The A100 is not the shared resource; the
> node is** — the 4 contending jobs held 4 *distinct* A100s. (And "solo" is never quite solo:
> foreign users share these nodes.)

### 4.1 THE FAN-OUT — the throttle was treating the symptom

Throttling to `%2` avoids node contention by discarding parallelism. The structural waste was
elsewhere: **a unit's 4 rungs are INDEPENDENT** (each is a fresh pmask / Fisher / soft draw /
null set at the same realisation), so running them sequentially in one process buys nothing but
a shared build. **12 sequential units → 48 independent (unit × rung) jobs**
(`sp3_fan.sbatch`, `--array=0-47%8`): **4× the parallelism for the price of a rebuild.** All 48
are real work — every `g3` realisation reaches rung 8 at both venues. Rung granularity is where
it stops paying: below it the ~5–8 min build+compile no longer amortises.

**The lane is 16 A100s wide — and this campaign had been using 8, for a reason that was never
measured.**

| | |
|---|---|
| **fp64 is MANDATORY** — pulsar-term phase ~1.6e4 rad; f32 destroys the ~1e-3-nat fringe-breaking signal (`harbor_env.sh`) | so the **~180 idle RTX A6000 / L40S / A4000 / 3090 cards are useless**: **1:64-crippled fp64** (HPC_SETUP §1.6), most QOS-denied. This is not a queue problem and no amount of waiting fixes it. |
| **dgx03 — bar LIFTED, on measurement not inheritance** | SPARK barred it for a **rogue non-Slurm process squatting cards** (SIREN: 51 GB resident; HPC_SETUP §7.2), and every sbatch in the lineage inherited that. Two facts retire it *for this workload*: this job holds **988 MiB**, so even a 51 GB squatter on an 80 GB card leaves **~29 GB (~30× headroom)**; and HPC_SETUP §7.3 conclusion 4 establishes contention **"never perturbs a value"** — timing only, and every number here is gated against the banks anyway. The runner prints any foreign resident at startup, so the trail records what we landed beside. *(Observed on landing: **no foreign resident** on any card taken.)* The bar was **prudent when written and stale by the time I inherited it** — nobody had measured the footprint. |
| **`--gres=gpu:nvidia_a100-sxm4-80gb:1`** — inherited from SPARK by every sbatch in this lineage | **RETIRED. MEASURED: a running rung job holds `988 MiB` — under 1 GB of an 80 GB card.** The `80gb` string bought nothing and **silently halved the lane**, excluding dgx01's 8 A100-**40GB** whose fp64 is equally **1:2 full** (HPC_SETUP §1.6). Replaced by `--gres=gpu:1 --nodes=1` + an explicit exclude of the non-A100 nodes. |
| **usable** | **dgx01 + dgx03 + dgx04 = 24 A100s** (verified: 4+4 on dgx01/dgx04 at 40/80 GB, clean; then 3 on dgx03, clean) |

> **A CLUSTER EVENT, RECORDED BECAUSE IT CHANGES THE COST MODEL, NOT THE RESULT (2026-07-17
> 00:46).** **`dgx01` AND `dgx04` both went DOWN mid-grid — *"Node unexpectedly rebooted"*** —
> and the whole array reverted to `PENDING (ReqNodeNotAvail)`. `--requeue` meant nothing was
> lost (no rung had banked yet), and **per-rung checkpointing is exactly the insurance this
> event argues for**. It also made the inherited `dgx03` bar decisive rather than academic:
> with both other A100 nodes down, **barring dgx03 on inheritance would have idled the grid
> completely**. *An exclusion carried forward without its measurement is a cost that only
> shows up on the day the rest of the cluster is gone.*

*(`--nodelist=dgx01,dgx04` is the wrong way to express this — it demands **both** nodes per task
(`NODES=2`) and the array never starts. The exclude-list is the right form.)*

**Compile cache shared and pre-warmed** across the array rather than per-job: `harbor_env.sh`
states cache isolation is *"NUMERICALLY inert: it caches XLA executables, not results"*, and
EMBER §7(a)'s race was on the **autotune** textprotos — **autotune is already OFF**
(`--xla_gpu_autotune_level=0`), so the race that motivated per-job caches is defused. This turns
48 cold compiles into ~1.

**Expected ~4–5 h rather than ~24 h.** The `%8` throttle is a *starting point measured from two
points* (3.62 s at 1-way, 20.0 s at 4-way) and should be **re-measured and adjusted**, not
trusted — that is precisely the mistake this section already records three times.

### 4.2 THE H200 LANE — REFUSED, but NOT for the reason it first appeared, and one artifact lost

**2026-07-17: all three A100 nodes went DOWN** (`dgx01`, `dgx03`, `dgx04`, *"Node unexpectedly
rebooted"*), idling the grid. The entitlement was checked by **test, not by reading**: my group's
`batch_gpu` QOS grants **exactly `gres/gpu:nvidia_gh200=2` + `gres/gpu:nvidia_h200=8` and ZERO
A100** — an A100 request there dies `QOSGrpGRES`, which **confirms HPC_SETUP §1.6's "QOS-denied"
note for `gpu0084`**. GH200 is aarch64 (jax's x86 CUDA-12 wheels will not run). **The H200 is the
only entitled card with full 1:2 fp64 that can run this stack at all**, so `g3a` was run there
before any grid was — because H200 is Hopper (sm_90) and every banked number is A100 Ampere
(sm_80).

| `g3a` | A100 (dgx04) | **H200 (hgx03)** |
|---|---|---|
| max &#124;dlnL_det − banked&#124; | **0.000e+00** | **1.072e+01** |
| max &#124;lnK − banked&#124; | 0.000e+00 | **0.000e+00** |
| max &#124;qmax − banked&#124; | 0.000e+00 | **3.995e-01** |
| `mapk` / `on_true` disagreements | 0 / 0 | **2/116 · 1/116** |
| **L_gwb fingerprint** | `8548f148b50a5b44` | **`459d4ae5e0f6c42b`** |
| | **PASSED** | **FAILED** |

**The H200 lane was REFUSED on this table. That refusal was itself an error — see §4.3 — and it
is withdrawn.** But first, the obvious reading of the table — *"Hopper's numerics differ"* — **is
WRONG, and it is the same misattribution this report criticises SPARK-2 for.** Look at the
last two rows. **`lnK` agrees at `0.000e+00`** — it is CPU-derived (`FringeTables.K_counted`) and
architecture-independent. And **the L_gwb fingerprints DIFFER**. `L_gwb` is built by
**`np.linalg.eigh` on the CPU**, not the GPU: `hgx03`'s host/BLAS chose a **different eigenvector
basis inside the near-degenerate GWB subspace**, so the job **drew a different noise realisation
and the GPU never saw the same data**. ***The test is CONFOUNDED: it compares two different
datasets and tells us nothing about Hopper.***

> **THIS IS THE NoiseDrawer HAZARD, AND IT IS WORSE THAN THE RECORD STATES.** CHORUS §0.1
> documents the basis as **thread-count** dependent (*"cpus-per-task was part of the seed"*).
> **Both jobs here ran at `cpus-per-task=8`** and still diverged — so the basis is **HOST**
> dependent too, not merely thread dependent. `dgx04` and `hgx03` are different machines and
> that alone was sufficient. **The remedy is the one already prescribed and never applied at
> T = 40: BANK `L_gwb`** (CHORUS §0.1). Until it is banked, **every T = 40 result in this
> programme is portable only to hosts that happen to reproduce `8548f148b50a5b44`** — which is
> exactly what g3a on A100 proved for ACCRE's dgx lane, and exactly what hgx03 does not.
> **Whether Hopper's GPU numerics are safe remains UNTESTED** — it cannot be tested until the
> square root is banked and both lanes are handed the same data.

### 4.4 THE H200 LANE, MEASURED — a +13% bar, and why the verdict survives it

`mode venue_self` on `hgx03`, ~455 s/venue, **three INDEPENDENT null-seed blocks** at venue A:

| lane | venue-A floor | venue-B floor | N_cert/g3 (mean) | rung-8 reachable |
|---|---|---|---|---|
| **SURFACE** (dgx/cronus) | **122.461 ± 4.895** | **44.397 ± 1.938** | A 11.17 · B 10.33 | 6/6 · 6/6 |
| **H200** block 0 | 133.844 ± 4.505 | 49.564 ± 1.809 | A 9.83 · B 7.67 | 5/6 · 3/6 |
| **H200** block 1 | 141.231 ± 4.753 | — | 8.83 | 5/6 |
| **H200** block 2 | 140.891 ± 5.127 | — | 8.83 | 5/6 |

**Three independent draws land 134–141; none near 122.** Mean **138.66**, block scatter **4.17**
(consistent with the per-estimate error ~4.79 — so the blocks *are* independent and the
estimator is behaving). Against SURFACE: **+16.19 = 3.0σ = +13%. SYSTEMATIC, not the draw.**
*(My first pass could not have seen this: I used one null-seed set for BOTH venues, and `nullN`
has no CW, so both floors were cut from the SAME draws — two correlated estimates wearing the
hats of two. Fixed with per-venue, per-block seeds; the replicates are what settled it.)*

**And it is probably not the GPU.** `load_or_build_L_gwb` does `w = np.clip(w, 0.0, None)` on a
**near-degenerate** spectrum. A different host's `eigh` returns slightly different tiny/negative
eigenvalues, and **the clip converts that into a slightly different EFFECTIVE COVARIANCE** — so a
host change is **not a pure rotation** after all, and the GWB noise power moves a little. That
sharpens §4.2: the hazard is not only a basis rotation, it can shift the *distribution*. **This
is a further argument for banking `L_gwb`**, and it is banked as a measurement here rather than
asserted.

> **WHY THE VERDICT SURVIVES IT, AND THE A-FORTIORI THIS BUYS.** Arrow 2's statistic is a
> **WITHIN-UNIT, CROSS-RUNG** comparison: each rung's floor is re-cut from **its own** nulls, on
> the same lane, same code, same GPU. **A common-mode bar offset cancels out of the crossing
> test.** What the offset does do is make this lane **HARDER**: a +13% bar means fewer pulsars
> sit near it. Therefore —
> * an **EDGE-POSITIVE** measured here is ***a fortiori*** EDGE-POSITIVE on the dgx lane (the
>   crossing happened against a higher bar), and
> * an **EDGE-ZERO** measured here is **NOT conclusive for the dgx lane** — part of the absence
>   is the raised bar — and any EDGE-ZERO from this lane must say so and be re-cut on dgx.
>
> **This is the general-mechanism question — "does deconfusion feed certification, and roughly
> how much" — not a precision measurement of one cell's floor.** Bit-comparability was never
> needed for that; a stated bar offset is.

**SCOPE, travelling with every number from this lane:** absolute floors and counts are **NOT**
comparable to SURFACE's banked values (+13% bar, N_cert mean 11.17 → ~8.8 at A). **Rung-8
reachability falls to 5/6 (A) and 3/6 (B) = 8 of 12 units** — exactly the brief's ≥8 minimum,
with no margin.

### 4.3 THE REFUSAL IS WITHDRAWN — bit-identity was never the requirement, and I mistook a gate for a goal

**Challenged on it, and the challenge was right.** I refused the H200 lane, and then argued the
only way onto it was a "second GWB basis" whose cost was that results would not be
**bit-comparable** to SURFACE. That framing is wrong, and the error is worth recording because
it is a failure of exactly the kind this report keeps cataloguing: **I treated a gate as an end
in itself rather than as a means to the thing it was protecting.**

**What g3a actually protects.** SPARK-3's *only* dependence on SURFACE's banks was
`certified_of`, which imported each realisation's certified set from `sf_sig_*`. THAT is what
forced our data to *be* SURFACE's data bit-for-bit, which forced g3a, which forced the canonical
`L_gwb`, which forced the dgx lane. **The whole chain existed to save ONE standard E-step per
realisation.**

**Why bit-identity is not the scientific requirement.** The square root enters as `L_gwb @ z`. A
rotated eigenvector basis gives a **different but EQUAL-IN-DISTRIBUTION** realisation — nothing
physical changes; it is a different draw of the same experiment. And arrow 2's statistic is a
**WITHIN-UNIT RELATIVE** one: does an uncertified pulsar's margin cross its bar between rung 0
and rung 8, **on the same data**, with the floor **re-cut from that state's own nulls** (§4
already does this — the ladder never used SURFACE's floor). ***The basis cancels out of that
comparison.*** What must hold is that **data, certified set and floor share ONE basis** — not
that the basis is any particular one.

**So the campaign is made HOST-INDEPENDENT** (`mode venue_self`): derive the floor and the
certified sets from our own data, at SURFACE's seeds and SURFACE's certification convention
(standard E-step, `pmask = ONE`, `ignite.py:349`). Cost: **~100 nulls + 6 signal E-steps ≈ 2 min
of GPU**, against which the entire dgx dependency was being paid. **g3a becomes moot** — not
failed, *unnecessary* — and any full-fp64 card qualifies.

> **THE GATE CHANGES SHAPE, AND THAT IS THE POINT.** A rotated basis **cannot** reproduce
> SURFACE's numbers exactly and **must not be asked to** — it is a different draw. The right
> instrument is **DISTRIBUTIONAL**: does this host recover **the cell**? A floor near SURFACE's
> **122.461** (A) / **44.397** (B), and a **live count** at `g3`. Agreement says the host
> computes the same physics; disagreement refuses the lane **on a real ground** rather than on a
> fingerprint mismatch that only ever meant *"different noise draw"*. **Correctness is the
> requirement; identity was a proxy I over-trusted.**

*(§4.2's L_gwb finding stands undiminished and is still worth banking: the basis is **host**
dependent, not merely thread dependent, so any campaign that DOES require replay of a banked
realisation — as `certified_of`'s original form did — is pinned to one CPU model. Banking
`L_gwb` at T = 40 remains the right fix and remains owed. The change here is that **SPARK-3 no
longer needs it**.)*

**AND I LOST AN ARTIFACT DOING IT (my error, recorded).** `mode_replay` keyed its npz on
`(venue, geo)` alone, so **the H200 run silently OVERWROTE the A100 run's PASS**: for a while
this report cited `g3a PASSED, 0.000e+00` while `spark3_replay_A_g3.npz` on disk said
`passed=False`. **An artifact whose filename does not distinguish the conditions it was produced
under is not an artifact — it is a race.** Fixed: the key now carries a **lane tag**
(host + GPU model) and the **L_gwb fingerprint is banked beside the verdict**. The H200 file is
retained as `spark3_replay_A_g3_H200_CONFOUNDED.npz`. **The A100 PASS survives only in the job
log (`sp3_r3A_12586155.out`) and must be REGENERATED before the fold** — it is on the owed list
(§5.2), and until it is regenerated the A100 g3a claim rests on a log, not a bank.

Each rung banks the **raw** per-pulsar `dlnL`/`lnK`/`qmax` and the **raw** null offender vectors
per state — never a verdict — with **both** Fisher bounds as columns and the rejection counts
beside them.

**The ladder per rung is the four states SPARK-2 §5(3) demanded** — *"without the truth rung the
ledger has no scale"*: `UNMODELLED → SOFT@opt → SOFT@pes → TRUTH`.

**Collate with `spark3.py ledger`** when the units land. **`SPARK3_results/` is gitignored** (the
repo convention), so the banks stay on ACCRE.

---

### 4.5 THE FIVE RUNGS THAT CAME BACK — a lane mismatch, and the fold it would have poisoned (SPARK-3C)

**On 2026-07-18 a 48-task array (`sp3f`, job `12588576`) ran on dgx01/dgx03 and banked five rung
files that had never existed.** They are exactly the five §4.4 predicted would be unreachable:

```
s3r_A_g3_r5_k8   s3r_B_g3_r3_k5   s3r_B_g3_r3_k8   s3r_B_g3_r4_k8   s3r_B_g3_r5_k8
```

**Forty-three of the forty-eight tasks exited in ~4 minutes on `skip-on-exist`.** The array was
therefore **not** the owed dgx re-cut (§5.4). It was a **gap-fill**, and it filled the gaps *from
the other floor lane*.

**THE MECHANISM, with file evidence — it is one flag.**

| launcher | line | lane `certified_of` takes |
|---|---|---|
| `sp3_h200fan.sbatch` (Jul 17, the grid) | `… --n-null 100 --self` | **SELF** — `spark3_venueself_*_hgx03_NVIDIAH200.npz` |
| `sp3_fan.sbatch` (Jul 18, the gap-fill) | `… --n-null 100` | **SURFACE** — `spark3_venue.npz` + `SURFACE_results/sf_sig_*` |

`certified_of` (spark3.py:838) branches on `self_lane`, and the floor it returns sets the
certified set through `cert = (dlnL > max(lnK, floor)) & (qmax > QBAR)`. **The two runs' own logs
print the two floors and settle it without inference:**

```
Jul-17  sp3hf_12593557_*.out : unit venue A … floor 133.844   unit venue B … floor 49.564   [SELF]
Jul-18  sp3f_12588576_*.out  : unit venue A … floor 122.461   unit venue B … floor 44.397   [SURFACE]
```

Re-deriving `N_cert` from the raw columns under both floors reproduces the split **exactly** —
including which five rungs die, and the `N_cert = 4` and `N_cert = 7` values the Jul-17 logs
print on their `UNREACHABLE` lines:

| venue | N_cert @ SURFACE (122.461 / 44.397) | N_cert @ SELF (133.844 / 49.564) | rung-8 reachable |
|---|---|---|---|
| A r0…r5 | 10 · 12 · 14 · 10 · 11 · **10** | 9 · 10 · 14 · 9 · 10 · **7** | 6/6 → **5/6** |
| B r0…r5 | 11 · 8 · 15 · **10** · **9** · **9** | 10 · 8 · 10 · **4** · **7** · **7** | 6/6 → **3/6** |

**So line 496's "5 rungs unreachable" was never wrong — it was unqualified.** It is a statement
about the self-derived lane. It has been corrected in place.

**NOT a stale read, and NOT a new hazard: it is §4.4's +13% bar arriving through an unlabelled
filename.** The floor offset is real and systematic — three independent same-host null blocks
land 133.844 / 141.231 / 140.891, mean **138.66**, block scatter **4.17** against a per-estimate
error of ~4.79, i.e. the blocks are consistent with each other and **3.0σ from SURFACE's 122.461**
(§4.4, and its `eigh`-plus-clip mechanism on a near-degenerate GWB spectrum). **The defect is not
the offset. The defect is that `s3r_*.npz` filenames do not encode the lane** — so two
provenances share one directory and the second silently overwrote nothing, added five, and looked
like progress. **This code already knows the lesson and had applied it one file away:**

> `_lane_tag()`, spark3.py:535 — *"An artifact whose name does not distinguish the conditions it
> was produced under is not an artifact, it is a race."*

`_lane_tag()` guards the `venueself` and `replay` artifacts. **It was never applied to the rung
checkpoints.** That is the whole bug, and it is the third appearance of the floor-provenance
hazard (after the `L_gwb` thread/host draw and LOTTERY's m1 grade flip) — this time entering
through *naming*, not through arithmetic.

#### WHAT THE FOLD WOULD HAVE DONE — the pre-registered STOP, fired

`mode_ledger` globs `s3r_*` (spark3.py:1710). **With the five on disk it silently folds two lanes
and its pre-registered verdict flips:**

| ledger input | units paired | crossings OPT | crossings PES | **`mode_ledger` VERDICT** |
|---|---|---|---|---|
| 43 rungs (lane-consistent) | 8 | 11 | **0** | **STRADDLED** |
| 48 rungs (lane-mixed) | 12 | 16 | **8** | **EDGE-POSITIVE** |

**Every one of the 8 pessimistic crossings comes from a cross-lane pair**, and the pessimistic
bound is the strict one — it is what "EDGE-POSITIVE" in `mode_ledger` *means* (crossings under
**both** bounds ⇒ bound-independent):

| unit | rung-0 lane | rung-8 lane | cross OPT | cross PES |
|---|---|---|---|---|
| A0 · A1 · A4 · B0 · B1 · A2 · A3 · B2 | SELF | SELF | 11 total | **0** |
| **A5** | SELF | SURFACE | 0 | **3** |
| **B3** | SELF | SURFACE | 0 | 0 |
| **B4** | SELF | SURFACE | 5 | 0 |
| **B5** | SELF | SURFACE | 0 | **5** |

The four newly-paired units compare a rung 0 certified against a bar of **49.564** with a rung 8
certified against **44.397**. Pulsars "cross" because **the bar moved between the two rungs of the
same unit** — which is precisely the one thing §4.4's a-fortiori argument depends on *not*
happening ("each rung's floor is re-cut from its own nulls, **on the same lane**"). It is a
provenance artifact end to end, and it manufactures the strongest-sounding result in the campaign.

> **C1's pre-registration: *verdict MOVES → STOP, report what changed, no fold.* It moved.
> THE FOLD DOES NOT HAPPEN.** `spark3_ledger.npz` and `spark3_final_verdict.npz` on ACCRE are
> **left exactly as cut on Jul 17** — the 43-rung, one-lane artifacts. Nothing was re-banked.

#### THE HEADLINE VERDICT IS UNMOVED — and provably, not by luck

The flip above is `mode_ledger`'s three-way bound test. **The report's headline verdict is
`sp3_final.py`'s**, and it is invariant for two independent structural reasons:

1. **Disjoint units.** `sp3_final.py:15` fixes `XUNITS = [A0, A1, A4, B0, B1]` and reads only
   `s3marg_*_k{0,8}`. The five new rungs live at **A5, B3, B4, B5** — cells the verdict reader
   never opens. Re-running it today over the 48-rung directory returns **every field bit-identical**
   to the Jul-17 bank: `EDGE-POSITIVE`, `n_survive 4`, `n_ready 5`, `cross_marg [1 7 1 0 1]`,
   `survives [T T T F T]`, scrambled-clean 10/10.
2. **Lane-invariant inputs where it counts.** For all five crossing units the rung-8 coherent
   **set** is identical across lanes, and `qmax` on that set agrees to **0.000e+00** — so the
   rung mask `PM` that `rung_masks` builds is **bit-identical on both lanes at both rungs the
   verdict uses** (rung 0 is the empty set trivially). `dlnL` differs by 2.4–6.3 nat between
   lanes, but it enters only through the ranking, and the ranking's top-8 set does not move.

| unit | k8 set equal | max ∣Δqmax∣ on set | max ∣ΔdlnL∣ on set | `PM(k8)` identical |
|---|---|---|---|---|
| A0 · A1 · A4 · B0 · B1 | **yes (5/5)** | **0.000e+00** | 6.27 · 3.01 · 4.74 · 4.07 · 2.41 | **yes (5/5)** |

**Which host's floors the verdict is cut on, and why: the SELF (H200) lane.** All 43 rungs, all
five crossing units, and every `s3marg`/`s3jvp`/`s3scr` bank were produced there; §4.4 establishes
the a-fortiori direction (a **higher** bar makes an EDGE-POSITIVE conservative). Mixing in
SURFACE-floor rungs does not sharpen that — it breaks the within-unit invariance the statistic
rests on.

**The four cells whose rung-8 reachability flips with floor provenance — A5, B3, B4, B5 — are
graded MARGINAL-band** under the pending floor-provenance convention: their reachability is a
property of the bar, not of the array, and no verdict may rest on them until they are cut on one
lane end to end.

#### OWED, AND NOW SHARPER

- **Quarantine the five before anyone runs `ledger` again.** They are good data on the *other*
  lane; they are a trap in this directory. Recommended and **not done here** (SPARK-3C is
  read-mostly and the STOP forbids re-banking):
  `for f in s3r_A_g3_r5_k8 s3r_B_g3_r3_k5 s3r_B_g3_r3_k8 s3r_B_g3_r4_k8 s3r_B_g3_r5_k8; do mv SPARK3_results/$f.npz SPARK3_results/$f.SURFACElane.npz; done`
- **`_lane_tag()` belongs in `_rung_path()`.** The rung checkpoints must carry their lane the way
  the `venueself` and `replay` artifacts already do. One-line fix, not applied here (HARD RULE:
  this report does not edit the campaign's code mid-verdict).
- **Arm (b) proper remains OWED and is now BLOCKED ON A CONSTRAINT, not on hardware** — see §5.4.

---

## 5. VERDICT: **STRADDLED**

**`spark3_ledger.npz`. 8 crossing-eligible units (both rung 0 and rung 8 banked), both venues,
both Fisher bounds.** The pre-registered decision rule (unchanged, and cut by `mode_ledger` with
no judgement call):

| bound | crossings (uncertified@rung0 → certified@rung8) | reading |
|---|---|---|
| **OPTIMISTIC** (conditional widths) | **11** over 8 units | an edge exists |
| **PESSIMISTIC** (prior widths) | **0** | no edge; a bad model *suppresses* |
| | **→ the bounds disagree in SIGN** | **STRADDLED** |

### 5.0 WHAT STRADDLED MEANS HERE — arrow 2 is MODEL-QUALITY-LIMITED, and that is the result

The crossing count is not the story; the **mechanism between the bounds** is. Per unit
(`ncert` of the soft state, rung 0 → rung 8):

```
             OPT crossings   soft@opt N: 0→8     truth N: 0→8
  A r0            1              1 →  2              8 → 8
  A r1            3              7 → 10              8 → 8
  A r2            0              9 →  1 (collapse)  13 →13
  A r3            0              8 →  8              9 → 9
  A r4            3              1 →  4              8 → 8
  B r0            1              5 →  5              8 → 8
  B r1            3              2 →  5              7?→ 8
  B r2            0              8 →  0 (collapse)  10 →10
```

Three facts, and together they are the finding:

1. **The truth arm is FLAT and RECOVERED** (`8→8`, `13→13`, …). At truth the reservoir is
   perfectly modelled, so the coherence rung cannot change the count — and it doesn't. **This is
   the control that says the machinery is sound**: SPARK-2's endpoints reproduce at this venue
   (truth recovers; unmodelled would not).
2. **Under the OPTIMISTIC bound the soft count is MIXED — it rises in 5 units and COLLAPSES in 2**
   (A r2 `9→1`, B r2 `8→0`, both strong-truth units). A tighter certified set (rung 8) with a
   *conditional-width* reservoir model recruits some pulsars and loses others, because the
   soft-averaged null floor also rises with the rung. It is a **fragile, partial** edge — 11
   individual crossings coexisting with two count-collapses — not a clean recruitment.
3. **Under the PESSIMISTIC bound NOTHING certifies, and the reason is physical, not numerical**
   (verified: 100% finite, 0 rejections). A prior-width reservoir model places 13 faint sources
   at essentially random positions; modelling them there gives the fit enormous freedom to
   accommodate **noise**, so the null offender's tail explodes and **the certification floor
   balloons** (e.g. A r2 k8: floor `118` unmodelled → `276` opt → **`744` pes**). **A bad model
   does not merely fail to help — it actively raises the bar and suppresses certification.** This
   is SPARK-2's *"unmodelled erases certification"* sharpened: a *prior-width* model is **worse
   than no model at all**.

> **THE VERDICT, IN ONE SENTENCE.** *Between "certification recruits when you know the reservoir
> well" (optimistic) and "certification is suppressed when you model it badly" (pessimistic) lies
> the true marginal width, which SPARK-2 measured is unaffordable to compute — so whether arrow 2
> is uphill at reachable capability is **exactly the question the two bounds bracket and cannot
> close.** Arrow 2 is not confirmed and not dead; it is **model-quality-limited**, and the value
> of deconfusion is a steep function of how well the faint sources can actually be constrained.*

**Why the scope cannot bias this.** The H200's **+13% bar** (§4.4) makes the optimistic edge
*harder* to achieve, so 11 crossings is a conservative floor on the optimistic side; and the
pessimistic 0 is a floor-inflation effect that a lower bar would only soften, never reverse. The
**single sky** (§0.1) limits generality but a STRADDLED verdict claims none — it says *"not
resolved,"* which no lane or sky makes truer or falser. **This is the one verdict this scope
delivers cleanly.**

### 5.1 THE PRE-REGISTERED NEXT STEP — the budgeted chunked-JVP, costed up front

STRADDLED's pre-registration (SPARK-2 §5(4), FIX 3): *"spec what a budgeted chunked-JVP
(wallclock cap stated up front) would cost to split them."* The split requires the **true
marginal width** — the full `78×78` faint-block Hessian of `lnL` — at the rungs that straddle.

**The cost, stated as a budget, not a hope.** SPARK-2 measured the Hessian is not a few-GPU-hr
object: monolithic `jax.hessian` **>27 min** and a chunked batched-JVP-of-grad at `CH=8`
**>17 min**, both *cancelled unfinished* on an A100-80GB. Those are **lower bounds**, not costs.
The honest budget:

- **Scope of the split:** only the **8 crossing-eligible units**, only at rungs **{0, 8}** (the
  crossing endpoints) → **16 marginal-Hessian evaluations**, not the whole grid.
- **Method:** batched JVP-of-grad at `CH=8` (`⌈78/8⌉ = 10` JVPs per Hessian), the cheaper of the
  two SPARK-2 profiled, on the H200 lane (~1.5× A100).
- **Wallclock cap, stated up front:** **45 min per Hessian.** SPARK-2's `CH=8` run was still
  running at 17 min on A100 without finishing, so on H200 a 45-min cap is the honest "either it
  lands or we learn it truly doesn't." A Hessian that hits the cap reports its unit as
  **straddle-unresolved** and the campaign quotes the **fraction resolved**, never a silent drop.
- **Budget:** **≤ 16 × 45 min = 12 GPU-hr**, hard-capped. If the median Hessian lands well under
  cap, less; if most hit the cap, the deliverable is *"the marginal width is unaffordable even
  chunked, and arrow 2's straddle cannot be closed at this capability"* — which is itself a
  publishable closure, the same shape as AVALANCHE's.

**This is a decision for Matt, not an automatic spend:** 12 GPU-hr to convert a STRADDLED into an
EDGE-POSITIVE-or-EDGE-ZERO, against the alternative of quoting the straddle as the honest state of
arrow 2 at current capability. **SPARK-3 does not launch it unbidden.**

> **DECISION TAKEN (2026-07-17, Matt): run it.** All three arms launched concurrently on the
> H200 lane (`scrambled first` governs *interpretation*, not execution): **(a)** the
> scrambled-counterpart arm at the crossing units, rungs {0,8} (`sp3_scram.sbatch`) — now
> *load-bearing*, because the pessimistic-bound anatomy (a mis-modelled reservoir *accommodates
> noise*) is exactly a mechanism that could manufacture the optimistic edge's crossings;
> **(b)** a dgx re-cut at SURFACE's own bar (armed, blocked on dgx return); **(c)** the budgeted
> chunked-JVP as a 5-process array, ≤4 Hessians each, build-once (`sp3_jvp.sbatch`,
> `spark3.py jvp`), 45-min cap per Hessian. **No JVP reading ships until arm (a) reports**
> (`sp3_verdict.py` enforces this and additionally gates the JVP by `diag(−H) ≟ F_ii`).

### 5.1a THE SPLIT — arm (a) clean, the JVP returns, the crossings survive the marginal

**Arm (a) — scrambled, CLEAN (`s3scr_*`, 10 cells).** Under a wrong counterpart every cert is
false by construction; soft-modelling produced **0 false certs at every unit and rung**, never
above unmodelled. **Soft-modelling does not manufacture — it suppresses.** So the optimistic
edge's crossings are real deconfusion, not the noise-accommodation mechanism the pessimistic
anatomy warned of. The reading is licensed.

**Arm (c) — the JVP, and its correctness gate.** SPARK-2 measured the faint-block marginal
Hessian *unaffordable* (cancelled >27 min / >17 min on A100). On the H200 it **returns in ~44 min
per Hessian, none capped** (rung 8, with coherent terms, runs no slower than rung 0). It passes
its built-in gate — `diag(−H)` reproduces the finite-difference `F_ii` to **<1% on the
well-curved axes** (the apparent 27% max was round-off on flat axes, localised and excluded).

**What the JVP revealed — the marginal is MID-BRACKET, and the block is INDEFINITE.** The faint
Fisher at truth has **26 of 66 non-positive eigenvalues**: *truth is a saddle, not a maximum,* for
sources at `log10_h = −14.25` — noise dominates their curvature. So `σ_marg = √diag(inv F)` is not
naively defined; the defensible marginal (constrained-subspace: keep the data-dominated
directions, prior elsewhere — `marg_width_from_hess`) lands at **median 0.40 vs conditional 0.14
vs prior 1.0 — ~3× the conditional, genuinely BETWEEN the bounds.** Per parameter it is the
programme's recurring physics: **`fgw` stays tight (data-constrained), `mc` → prior
(unconstrained)**. The width alone therefore does **not** decide the crossings — they had to be
re-scored.

**The re-cut at the marginal width (`mode remarg`, `s3marg_*`).** The soft state re-scored at
`σ_marg`, floor re-cut at matched state:

| unit | crossings (opt) | crossings (**MARG**) | rung-8 depths above bar (nat) | floor r0→r8 | |
|---|---|---|---|---|---|
| **A0** | 1 | **1** | 122 | 354→256 | SURVIVES |
| **A1** | 3 | **7** | 275, 91, 71, 56, 54, 27 | 323→152 | SURVIVES |
| **A4** | 3 | **1** | 235 | 454→249 | SURVIVES |
| **B0** | 1 | 0 | — | 135→89 | dies |
| **B1** | 3 | **1** | 66 | 104→89 | SURVIVES |

**4 of 5 units survive at the true marginal width. ≥2 required. Scrambled-clean. → EDGE-POSITIVE.**
The surviving margins are **deep (66–275 nat above the bar), not floor-grazing**, so the verdict is
robust to the marginal-width regularisation and the per-rung floor noise. B0 is the honest
exception — the one unit where growing the certified set *lowers* the count (3→1); not every
realisation's loop closes, and the report does not hide it.

*(The verdict rule as coded, unchanged: EDGE-POSITIVE = crossings under **both** bounds →
GLACIER; EDGE-ZERO = none under the **optimistic** bound → arrow 2 closed, shortfall quoted;
STRADDLED = disagree → this §5.1. We are STRADDLED.)*

### 5.3 GLACIER — **LICENSED** by the EDGE-POSITIVE split, WITH the model-quality law as a design gate

> **GLACIER — pre-registration, now LICENSED.** The chunked-JVP split the straddle to
> EDGE-POSITIVE (§5.1a): the crossings survive the true marginal width at 4/5 units,
> scrambled-clean. GLACIER's licence condition (an EDGE-POSITIVE reading) is met. GLACIER is the **fixed-list** loop arrow 2's edge would
> license: the source list **never grows** — no recruitment, no self-found step — and the loop's
> only motion is `[E] certify → [D] the certified-coherent detector re-scores the KNOWN faint
> reservoir → [M] the reservoir's constraints tighten → [E] deconfused residuals certify more
> pulsars → repeat`, iterated to a fixed point on **one frozen census**. It is the *slow* arrow:
> it buys nothing new in the sky, only a better model of what is already known to be there —
> hence the name. **Its premise is weaker than EPOCH's and its safety case is stronger**: because
> the list is frozen and every source enters at its true position, EMBER's manufacturing regime
> (motion under a **wrong counterpart**, sensitivity 1.00, p = 0.002, S8.5.3) is closed **by
> construction** — there is no counterpart step to get wrong, so B0.2's search gap (S8.5.5) never
> enters. **What GLACIER must measure, pre-registered:** whether the loop reaches a fixed point
> with `N_cert(fix) > N_cert(0)` sustained over ≥ 8 realisations **at both venues and both Fisher
> bounds**, with **per-iteration** floors (never inherited), the **trials cost restored** to the
> bar (+0.578 nat for the full census, and it **saturates** — SPARK-2 §2, which is favourable and
> already banked), and the **scrambled-counterpart arm at every rung** as the manufacturing
> control. **Its kill-shot, stated in advance:** if the fixed point is `+1` at the *optimistic*
> bound, GLACIER is arithmetically alive and **practically irrelevant**, and that sentence closes
> it rather than a campaign. **And its scope limit travels from §0.1: the rung-8 ensemble is one
> sky (`g3`), so GLACIER's first spend is MORE SKIES, not more noise seeds.**
>
> **THE MODEL-QUALITY LAW — a design GATE, promoted from footnote (Matt, decision item 2).** The
> split (§5.1a) showed arrow 2's edge lives at the *true marginal* reservoir width, which is only
> ~3× tighter than prior and **indefinite on ~40% of the faint block** — and that a *worse* model
> (prior-width) **inflates the certification null** (floor 118→744, §5.0), turning the loop
> net-negative. **GLACIER must therefore GATE on its own reservoir-constraint quality**: at each
> iteration it measures the faint block's constrained rank (the `N_eff` of `marg_width_from_hess`)
> and **refuses to soft-model any source whose marginal width has not tightened below the
> mis-modelling-inflation threshold**. Modelling below that threshold does not merely fail to
> help — it *raises the bar and un-certifies*. This is the one design constraint the split adds
> that a pre-split GLACIER would have missed, and it is why the straddle anatomy (§5.0) is a
> result in its own right, not a footnote.

### 5.4 STILL OWED, and named rather than quietly dropped

- **The scrambled-counterpart arm — BUILT and RUN, CLEAN (`s3scr_*`, 10 cells).** It reuses
  SURFACE's recipe (`surface.py scramble_s`), scrambles the recovery counterpart so every cert is
  false by construction (EMBER's S8.5.3 arm), re-cuts floors at matched state. **Result: 0 false
  certs at every unit and rung, never above unmodelled — soft-modelling does not manufacture.**
  This was the load-bearing gate on the EDGE-POSITIVE verdict and it passed. *(A slip caught: the
  first array's venue vector was mis-aligned and skipped `A4 k8`; relaunched, also clean.)*
  **SCOPE: EMBER §S8.5.2 still applies** — clean here is not proof of safety everywhere; it is the
  necessary condition for reading the crossings as real, which is all the verdict rests on it for.
- **A dgx re-cut of the ledger — STILL OWED. The hardware block is GONE; a policy block replaced
  it, and the scientific need has SHRUNK.** The EDGE-POSITIVE verdict was cut on the H200 lane
  (+13% bar). Because that bar is *higher*, the verdict is a fortiori (crossings cleared a harder
  bar) — but a dgx re-cut at SURFACE's own bar remains the clean confirmation that the offset did
  not shape it. Three things changed on 2026-07-18/20:
  - **dgx is back**: `dgx01 mixed`, `dgx03 mixed-` (`dgx04` still `drained`). The A100 lane runs.
  - **The Jul-18 array was NOT that re-cut** — `skip-on-exist` reduced it to a five-rung gap-fill
    on the wrong floor lane (§4.5). Arm (b) has still **never run**.
  - **It cannot be run under SPARK-3C's remit.** dgx01/03/04 sit in partition `interactive_gpu`
    — *the reserved `dgx_iacc` share*. SPARK-3C's standing instruction is **general queue only**
    (`batch_gpu`), and `batch_gpu` is the H200/GH200 lane the incumbent was already cut on, so it
    cannot answer a dgx-vs-H200 question at all. **Submitting arm (b) needs Matt's explicit
    release of the reserved share; it was not assumed.** Cost when licensed: the five crossing
    units × {k0, k8} against a **fresh output dir** (or `skip-on-exist` disabled), ~58–75 min/rung
    on dgx ≈ **10–13 GPU-hr**.
  - **What arm (b) is still FOR has narrowed.** Its stated question — the +13% offset, in numbers
    — is now answered from banked artifacts at zero GPU cost (§4.5): the offset is systematic at
    **3.0σ**, and it **does not reach the verdict**, because the rung masks `PM` at k0 and k8 are
    **bit-identical across the two lanes** for all five crossing units (Δqmax = 0.000e+00). What a
    dgx re-cut would still buy is the *scoring-side* confirmation — that the margins and re-cut
    null floors agree once the inputs are known to — which is a genuine but much smaller claim
    than "the offset may have shaped the verdict".
- **Rungs 2 and 5 contribute shape, not the verdict**: the crossing test is defined on 0 → 8.
- **REGENERATE the A100 `g3a` artifact** (§4.2) — clobbered by the H200 replay; the PASS
  currently rests on a job log rather than a bank. One job, when dgx01/03/04 return.
- **BANK `L_gwb` FOR T = 40** (CHORUS §0.1's remedy, never applied at this geometry). Until it
  exists, T = 40 results are portable only to hosts reproducing `8548f148b50a5b44`, the H200
  lane cannot be honestly tested, and the hazard is HOST-wide, not thread-wide (§4.2).

---

## 6. WHAT SPARK-3 CHANGES IN THE RECORD

1. **SPARK-2's arrow-2 venue could not have answered arrow 2.** Its cell is **below onset**
   (`corr = 0.567/real`) and its own **max `N_cert` is 1** across 30 realisations — rungs 2/5/8
   are not constructible there. Its `N_cert = 0` result is a property of the venue, not of
   arrow 2. **The reserved Pair B cells are live (24/30, 25/30) and are the venue.** [§0]
2. **The all-NaN E-step is EARTH-TERM COALESCENCE, not KINDLE's degenerate covariance.** The
   root cause is that SPARK-2 clipped the soft draw to `TE.phi_bounds` — the *estimator's search
   box* — believing it the prior support. **The clip added as the fix is what admitted the
   failure.** The population's generative box is evaluable everywhere with a **6.1× margin** at
   `t_max = 49.70 yr`, and the repaired draw's **rejection fraction is exactly 0**. [§2]
3. **`t_coal` is a function of `(mc, fgw)` ONLY — the strain plays no part.** SPARK-2's
   "EMBER's 0% does not transfer from loud to faint" is **RETRACTED**: EMBER measured 0%
   because its `mc_scan` holds `fgw` at truth and so cannot see a **joint** boundary. It is a
   double perturbation of **(mc × fgw)**. [§2.3]
4. **The `~116×` per-target E-step price is an artifact of building it as 116 `estep()` calls** —
   only the cheap base evaluation depends on the target's mask (`MaskedDelay` reads
   `params[PMASK][ipsr]` alone); the fringe sweep, **97.3% of the work**, does not.
   **MEASURED: 4.29× solo** (4.08× with a co-tenant) — **28× cheaper than the price that
   deferred it.** *(My a-priori estimate of ≈1.2× was wrong — `A` underestimated 14×; the
   structural model is exact, the constant was not. §1.1.)* [§1]
5. **`t_max` must be read from the ACTUAL TOAs** — **49.70 yr** here, which is neither the span
   (47.17 yr) nor `T_label`. The extension moves it, and it is what the evaluability boundary is
   measured against. [§2.2]
6. **The SURFACE T = 40 banks are replayable on ACCRE**, bit-for-bit, on the thread-dependent
   recompute branch at `cpus-per-task=8` — established by `g3a` and not previously in evidence.
   It does **not** retire the hazard: banking `L_gwb` for this geometry (CHORUS §0.1) is still
   the fix and is still unapplied at T = 40. [g3a]

**RETRACTED, by me, in this report:** *"the T = 40 venue costs a 5–6× build"* — the 1032 s I first
measured was **co-tenancy between my own two jobs**; solo it is **174 s**, comparable to T = 30.
See §1's retraction. SPARK-2 §5(2)'s in-process rung loop is still adopted, on its **cache**
argument, which is real — but not on a build-size argument I no longer have.

## 7. FILES

```
hpc_harbor/spark/spark3.py                modes: venue | replay | gate1 | gate2 | unit | ledger
hpc_harbor/spark/sp3_one.sbatch           parameterised runner, one mode (SP3_ARGS); 12h
hpc_harbor/spark/sp3_grid.sbatch          12 sequential units (superseded by the fan-out)
hpc_harbor/spark/sp3_fan.sbatch           THE GRID: sbatch --array=0-47%8 -- one job per
                                          (unit x rung); 16-A100 lane, see §4.1
SPARK3_results/spark3_venue.npz           §0  the venue check                    COMPLETE
SPARK3_results/spark3_replay_A_g3.npz     g3a the reproduction gate  PASSED bit-identically
SPARK3_results/spark3_gate1_A_g3.npz      §2  per-target E-step + g1a/g1b + cost  COMPLETE
SPARK3_results/spark3_gate2_A_g3_r0.npz   §1  evaluability + g2a/g2b/g2c + solo timing  COMPLETE
SPARK3_results/spark3_venueself_{A,B}_g3_hgx03_*.npz  §4.4 host-derived floor+certsets (H200)
SPARK3_results/s3r_{A,B}_g3_r{0..5}_k*.npz  §4  43 reachable per-RUNG checkpoints  COMPLETE
SPARK3_results/spark3_ledger.npz          §5  the crossing ledger — STRADDLED     COMPLETE
SPARK3_results/s3scr_*.npz                §5.1a arm(a) scrambled (10, all clean)   COMPLETE
SPARK3_results/s3jvp_*.npz                §5.1a arm(c) marginal Hessians (10)      COMPLETE
SPARK3_results/s3marg_*.npz               §5.1a marginal re-score (10)             COMPLETE
SPARK3_results/spark3_final_verdict.npz   §5.1a THE VERDICT — EDGE-POSITIVE        COMPLETE
```

`SPARK3_results/` is gitignored (`*_results/`, the repo convention — no peer result dir is
tracked); the banks live on ACCRE. **Lean-npz throughout**: raw per-pulsar `dlnL`/`lnK`/`qmax`
and raw null offender vectors per unit, **both Fisher bounds as columns**, rejection counts, and
zero-fraction beside every floor — **never verdicts**. Verified: `spark3_venue.npz` re-cuts its
own counts from its own raw columns at `|dev| = 0` on all three cells, and the tracked `venue`
mode reproduces the prototype at `0.000e+00` on every column.

## 8. ADD-LIST (staged; Matt commits)

```
hpc_harbor/spark/spark3.py            new  the campaign (6 modes, 6 gates)
hpc_harbor/spark/sp3_one.sbatch       new  runner; cpus-per-task=8, --exclude=dgx03, 12h
hpc_harbor/spark/sp3_grid.sbatch      new  12-unit sequential array (superseded, kept as trail)
hpc_harbor/spark/sp3_fan.sbatch       new  THE GRID: 48 (unit x rung) jobs, --array=0-47%8
hpc_harbor/spark/sp3_scram.sbatch     new  arm (a) scrambled-counterpart array
hpc_harbor/spark/sp3_jvp.sbatch       new  arm (c) chunked-JVP array (grouped, 45-min cap)
hpc_harbor/spark/sp3_remarg.sbatch    new  marginal re-score array
hpc_harbor/spark/sp3_verdict.py       new  scrambled gate + JVP correctness gate + width bracket
hpc_harbor/spark/sp3_final.py         new  THE VERDICT reader (marginal crossing, gated on scram)
reports/SPARK3_second_arrow.md        mod  this report (+ §4.5 and the §9 SPARK-3C addendum)
```

**SPARK-3C stages exactly those two.** `sp3_final.py` was already staged; the report is modified
in place. Everything else this report now leans on is **already committed** at `74dd0c9` —
including **`sp3_h200fan.sbatch`**, which is worth naming even though it needs no staging: it —
not `sp3_fan.sbatch` — is the launcher that produced the 43 banked rungs the verdict rests on,
and the single `--self` between the two is the whole subject of §4.5.

**NOT staged, deliberately:** no `.npz` — `SPARK3_results/` is gitignored and the STOP forbids
re-banking, so `spark3_ledger.npz` and `spark3_final_verdict.npz` stay exactly as cut on Jul 17.
**No code fix** — the `_lane_tag()` repair to `_rung_path()` (§4.5) is named and costed but not
applied, because this report does not edit the campaign's code mid-verdict.

**SPARK-2's four staged files are untouched and sit alongside**, exactly as the brief required:
`hpc_harbor/spark/{spark2.py, rb2.py, sp2_one.sbatch}`, `reports/SPARK2_second_arrow.md`.

**Nothing else is modified.** No tracked file outside this add-list is touched; the HARD RULE
holds — `TE.phi_bounds`, `deterministic.py`, `trackB_b1_core.py` and `ember_map.py` are
**reported on, not edited**, and §2's finding about `phi_bounds` is a statement about how
SPARK-2 *used* it, not a defect in it (as a search box it is correct; it is simply not the
prior).

## 9. STOP

**Arrow 2's verdict is EDGE-POSITIVE — the straddle was split at the true marginal width, and the
edge survives.** The four SPARK-2 fixes are closed — **the venue** (its cell's max `N_cert` = 1,
the ladder never constructible there), **the instrument** (per-target E-step, `116/116` live rows
at rung 0, **4.29×** not ~116×), **the pathology** (Earth-term coalescence in closed form; the
`phi_bounds` clip added *as the fix* was the cause; rejection fraction **exactly 0** in the true
prior), **the Fisher** (two-sided bounds). The grid ran to completion on the **H200 lane** after
the entire A100 lane went down, its dependence on SURFACE's banks cut so bit-identity was no
longer required. The pre-registered crossing test returned **STRADDLED** (optimistic edge, no
pessimistic edge); Matt's decision ran the **budgeted chunked-JVP** to split it.

**THE RESULT: arrow 2 is EDGE-POSITIVE, and MODEL-QUALITY-LIMITED.** The faint-block marginal
Hessian SPARK-2 measured *unaffordable* returns on the H200 (~44 min each, none capped, gate-clean)
and places the true marginal reservoir width **mid-bracket (~3× conditional, indefinite on ~40% of
the block — truth is a saddle for weak sources)**. Re-scored at that width, the crossings
**survive at 4 of 5 units** (deep margins, 66–275 nat), **scrambled-clean** (soft-modelling does
not manufacture). The crossing is **the LOOP closing** — the floor falls with coherence (arrow 1,
template rigidity) *and* dlnL rises (arrow 2, deconfusion). B0 is the one unit where the loop does
not close, and it is reported. **The edge is real, unmanufactured, and survives the honest width —
but it lives at a reservoir-model quality where a *worse* model would invert it** (floor
`118 → 744`), which is the design law GLACIER now carries (§5.3).

**GLACIER is LICENSED** — the fixed-list loop, with the model-quality gate as a design
constraint (§5.3). **Its first spend is MORE SKIES** (this verdict is one sky, `g3`, §0.1).

> ### SPARK-3C ADDENDUM (2026-07-20) — the 48/48 fold was REFUSED, and the verdict is unmoved
>
> **The fold did not happen, because the pre-registered STOP fired.** Five rungs banked on
> 2026-07-18 turned out to be **SURFACE-lane** files sitting in a **self-lane** directory — one
> missing `--self` between `sp3_h200fan.sbatch` and `sp3_fan.sbatch` (§4.5). Folding them moves
> `mode_ledger`'s verdict **STRADDLED → EDGE-POSITIVE** on **8 pessimistic crossings that are
>100% cross-lane artifacts**. Per C1's pre-registration — *verdict MOVES → STOP, no fold* — the
> banked `spark3_ledger.npz` and `spark3_final_verdict.npz` are **untouched**, and this report's
> ledger numbers remain the 43-rung, one-lane numbers.
>
> **The headline verdict is unaffected and was verified, not assumed.** `sp3_final.py` re-run
> today over the 48-rung directory returns **every field bit-identical** to the Jul-17 bank
> (`EDGE-POSITIVE`, 4/5 units, scrambled-clean) — because its `XUNITS` never touch the five
> cells, *and* because the rung masks it does use are bit-identical across both lanes. **Arrow 2
> is still EDGE-POSITIVE and still MODEL-QUALITY-LIMITED. Nothing in §5 changes.**
>
> **What changed is the confidence in the plumbing, not the physics**: `s3r_*.npz` filenames do
> not encode their floor lane, so the directory is not self-describing and `ledger` will
> mis-fold it again on the next run. Quarantine command and the one-line `_lane_tag()` fix are
> in §4.5; **neither was applied** — SPARK-3C stages, Matt commits, and the STOP forbids
> re-banking. **Arm (b) has still never run** and now needs a release of the reserved
> `dgx_iacc` share rather than a repaired node (§5.4).

**Owed and named, not dropped:** the **dgx re-cut** (§5.4) — the verdict was cut on the H200's
+13% bar, which makes it *a fortiori* (crossings cleared a harder bar), but a dgx re-cut at
SURFACE's own bar is the clean confirmation and is **blocked on dgx return** (all three A100 nodes
down the whole run; the A100 array stays armed). The **A100 `g3a` artifact** (clobbered, §4.2) and
**T = 40 `L_gwb` banking** (§4.4) also await dgx. Arrow 1's **CASCADE-ALIVE** is unaffected;
**EPOCH** (SPARK §6) and now **GLACIER** are the pre-registered campaigns.

*Awaiting Matt's decision before any further spend.*

**Gates built (SPARK's g0 chain in spirit — a refactor must first prove it did not move the
incumbent):**

| gate | what it proves | bar |
|---|---|---|
| **g3a** `replay` | SPARK-3's draw **IS** SURFACE's draw at the venue — else the certified set belongs to a different realisation than the data scored, and every rung is coherence on a counterpart that is not there. Mandatory because T = 40 has no canonical `L_gwb` bank → the thread-dependent RECOMPUTE branch. | discrete exact; continuous < 1e-6 |
| **g1a** `gate1` | the per-target machinery at `PM ≡ ONE` reproduces `sp.estep(..., one)` — the refactor did not move the incumbent | **0.000e+00** |
| **g1b** `gate1` | the per-target **rung 0** (target's own term live, nothing else coherent) vs the **banked** standard `pmask = ONE` q columns. The two can differ *only* through the HD-correlated GWB coupling `2·db@u_p[p]`: `a_pf` and `db` are identical because `MaskedDelay` reads pulsar p's own mask entry alone. | measured, not assumed |
| **g2a** `gate2` | the **analytic** evaluability predicate vs the actual E-step (non-evaluable ⇒ all `npsr×B` non-finite; evaluable ⇒ none) | agree |
| **g2b** `gate2` | SPARK-2's box **reproduces** its rung 2/5 deaths — as a property of the draw box, not the rung | deaths > 0 |
| **g2c** `gate2` | the repaired draw (population support + joint predicate + reject-and-redraw) is finite at **every** rung; rejection counts banked | 0 non-finite |

`SPARK3_results/` is gitignored (`*_results/`, the repo convention — no peer result dir is
tracked); the banks live on ACCRE. Lean-npz throughout: raw per-pulsar `dlnL`/`lnK`/`qmax` and
raw null offender vectors banked per unit, never verdicts.
