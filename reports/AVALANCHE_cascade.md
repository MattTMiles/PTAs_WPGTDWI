# AVALANCHE — the multi-source cascade: PRE-FLIGHT STOP at gate g0

**Status: STOPPED BEFORE LAUNCH. No GPU time spent. No realisations banked.**
**Verdict: the campaign as briefed is NOT EXECUTABLE — and its central premise is
already REFUTED by a banked measurement in this repo (`project_progress.md:1450-1462`).**

**Tree.** `criterion-v2.2` EXISTS (contra the `kindle-gain-diagnosis` memory note, now
corrected below) and points at `db2075a` "floor adopted, onset recut to 59". HEAD of
`MM_playground` is `d87db93`, exactly one commit later — the KINDLE stage-0 commit that
*contains* `reports/KINDLE_gain_contour.md`, a required reading. HEAD therefore satisfies
"criterion-v2.2 or later" and strictly dominates the tag. **I did not check out the tag
detached:** it would drop KINDLE and put 41 dirty paths of prior-campaign work at risk.
Work performed on `MM_playground` @ `d87db93`.

**REALISM NOTE (the header the brief asked for).** Had this run, it would have run on the
realistic spine — real pulsar positions, real TOA uncertainties, real published distance
priors, drawn white + per-pulsar red + GWB noise, honest per-cell floors. **But the spine
is still a MOCK** (`data-spine-is-simulated`; AXIS, 1440 MHz single-frequency, §10.15(a)):
no real TOA is touched, and the residuals *are* the injected CW+CURN. A verdict here would
have been a statement about an approximately-real IPTA's *geometry and noise budget*, not
about real data. This distinction is load-bearing for §5 below.

**SCOPE OF INFERENCE.** Everything below is a statement about *the machinery and the
banked record*, established by reading source and the progress log. It is not a physics
measurement and is not offered as one. The one physics claim I rely on (§1) is quoted from
a banked prior measurement, not re-derived here.

---

## 0. BUDGET ESTIMATE (before spending, per the brief) — the gate the brief pre-registered

The brief's own STOP: *"if the per-iteration null refits push past ~80 GPU-hr, STOP and
propose the cheaper pre-registered fallback."* That gate fires, by a wide margin, on
arithmetic alone.

Cost basis, banked (`reports/KINDLE_gain_contour.md:339-342`): IGNITE-2's per-realisation
wall is ~750 s first-in-process (evaluator materialisation), then 118 s (2-iter,
step-rejected) to **320 s** (full 10-iter). Off-truth arms accept steps → the 320 s class.

| component | count | @300 s | GPU-hr |
|---|---|---|---|
| signal (4 capability cells × 4 igniters × ≥10 reals) | 160 | | **13** |
| scrambled through the full extended loop (2 cells × ≥5) | 10 | | 1 |
| **per-iteration floors — the brief's own novelty** | **≥100 nulls × 12 iters × 3 cells** | | **300** |
| iteration-0 floors elsewhere (13 cells × 100) | 1 300 | | 108 |

**The per-iteration null refits alone are ≈300 GPU-hr — ~3.75× over the brief's own
~80 GPU-hr STOP.** And this costing is *optimistic by a large factor*, because it prices
the loop **without** the [D] step. The detect step is not free and not built (§2): the
Earth-term F-stat scan is 47 freq × 192 healpix = **9 024 grid points**, each a 40-step
Adam profile over 4 free parameters with a gradient through all 116 pulsars
(`trackB_estimator.py:197-262`). A *pulsar-term-coherent* scan — the thing the brief
actually requires — multiplies that by the fringe marginalisation, **B = 512 evaluations
per pulsar** (`trackB_estimator.py:85`), per grid point, per iteration, per realisation,
**and per null**.

**The pre-registered fallback (full per-iteration floors at 2 cells + iteration-0 floors
elsewhere) does not rescue this**, because the fallback prices only the floors. It leaves
§1–§4 untouched. I am therefore **not** proposing the fallback: it is a cheaper way to run
a campaign that should not run.

---

## 1. THE CENTRAL FINDING — the cascade's premise is already refuted, and the project already pivoted because of it

This is the finding that matters. It is not a code defect and no amount of engineering
retires it.

AVALANCHE's cascade requires that a **self-found** source (one the [D] step detects, rather
than one placed at truth) can support **certification**. The brief is explicit and correct
that this is the crux: *"the detect step's new sources are SELF-FOUND counterparts, so
their wrongness rate is a first-class readout, not an assumption."*

**That wrongness rate is already measured. It is ~100%, and it is measured at zero noise.**

`project_progress.md:1450-1462`, the **B0.2 Asimov gate**, recorded as *"the make-or-break
finding"*:

> The F-stat seed scan **recovers 3/3 loud sources** (6–12°, <1 freq-bin, 2F 184–595 vs
> floor 15) → **seeding gate 3a PASS**. **B0.2 Asimov gate FAIL:** seed+EM converges but
> lands **median source err ~0.5 (scaled)** and **the ceilings collapse to 0**. The
> `perturb` scan measures WHY: the certification tolerance on the SOURCE parameters is
> **~1e-4 scaled**, set by the pulsar-term phase **LEVER ARM 2πL/dL ~ 1.6e4 rad** — any
> source error shifting that phase by ~1 rad re-registers the fringe pattern and the true
> fringe stops being the posterior mode. Achievable source precision ~0.5 scaled is
> **3–4 ORDERS coarser** → **the A2 census "3 certified pulsars" is a CONDITIONAL-on-truth
> ceiling, NOT achievable by any cold-start estimator, even at zero noise**. P0 confirms
> the posterior is **CONFIDENT-WRONG** at imperfect sources (q_max 0.5–0.99 on a SHIFTED
> fringe, k≠0), not diffuse — **a wrong-certification machine**.

Read that against **EMBER's safety boundary**, which the brief correctly identifies as the
governing hazard: **manufacturing = (wrong counterpart) × (MOTION)**, motion sensitivity
1.00 (5/5), Fisher p = 0.002; engagement refuted as the boundary (2/5 manufacturers started
disengaged).

The two findings compose, and the composition is the whole answer:

- A self-found source is **wrong by 3–4 orders** relative to the certification tolerance
  (0.5 scaled achieved vs ~1e-4 required). **It is a wrong counterpart, essentially always.**
- A loop that adds sources **moves by construction** — every accepted add is an M-step.
- Wrong counterpart × motion = **the manufacturing regime**, exactly.
- And the failure is **not** benign-diffuse: it is **confident-wrong** — q_max 0.5–0.99 on
  a shifted fringe. It therefore **passes the q>0.9 certification layer while being wrong**.

**AVALANCHE as briefed would reliably produce a CRITICAL verdict — a textbook snowball —
composed entirely of manufactured certifications, and its own q-based purity layer could
not detect this.** D1 (the wrong-counterpart hole) is **OPEN**; D3 and D4 were both
rejected; no purity layer exists. EMBER measured D1 travelling through the loop
*unamplified* — but, per EMBER §1.5 and §5, only because that loop was **inert**
(steps = 0.00 in 7 of 9 map cells, ΔN = 0.000 in all 9). **AVALANCHE's design deletes
precisely the property that made D1-unamplified true.** The growing source list is not a
new arm on a safe loop; it is D1's amplification channel, opened deliberately.

**The project already knows this and already pivoted.** `project_progress.md:1462` —
*"PIVOTED (with Matt) to the joint-registration reformulation."* The pivot's result
(`:1476+`) is that the registration needle **EXISTS and sits exactly at truth** — argmax at
truth at every scale, next-best cell −102 nat, an interior 2-D maximum — but it is an
**unresolved CUSP of width ~1e-5°** (array lever-arm prediction; the brief's 0.006° is the
J1713-class value, the array max is far larger). **The seeder delivers 6–12°.** The gap
between what a blind search delivers and what certification demands is **six to seven orders
of magnitude.**

That gap is why every campaign since has anchored rather than searched. EMBER's
basin-selection convention (`ember_map.py:60-74`) states it outright: the constrained axes
are **basin-anchored at the Earth-term detection's own basin**, and only `log10_mc` — the
one loose axis (Fisher info gain 1.00–1.73×, *"the posterior IS the prior"*) — is drawn
cold, **because a blind global search is a property of the SEARCH, not the loop.**

> **AVALANCHE's [D] step asks the loop to do the one thing this repo has measured it cannot
> do, and has excluded by convention in every campaign since 2026-07-06.**

---

## 2. GATE g0 CANNOT BE CONSTRUCTED — the detect step's feedback channel does not exist

Independently of §1, the [D] step is unbuildable as briefed.

The brief: *"[D] DETECT: **with the current certified set's pulsar terms coherent**
(soft/q-weighted per spec §3 — never hard-locked), run source detection over the residual
band: the F-stat seeder (Track B machinery, data-driven thresholds, exclusion radii)."*

**The F-stat seeder is Earth-term-only, structurally, by construction:**

- `trackB_estimator.py:228` `seed_scan` → `:182` `build_earth_single` → `:188`
  `build_fisher_amortised(..., msd_factory=TE.EarthDelay)`.
- `trackB_estimator.py:56` `class EarthDelay` — docstring: *"Earth-term-only multi-source CW
  delay (**the fully-decohered coherence limit**)"*; `:69`
  `make_phase_connected_binary(pulsarterm=False)`. The pulsar-distance parameter is kept
  *"inert (delay ignores it)"* purely to keep the theta layout aligned.
- Gated identity: `trackB_b1_core.py:754-758` — *"G1: mask=0 == EarthDelay (pterm off)"*,
  reproduced 0.00. **`EarthDelay` IS the mask=0 case.**

So the seeder sits at **the fully-decohered pole — the exact opposite of the coherent limit
whose improvement is the cascade's entire premise.** Certified pulsar terms **cannot enter
it**. Wired as briefed, [D] returns **the same source list at every iteration, regardless of
the certified set**. The [E] → [D] feedback edge — the one that closes the loop and makes it
a cascade rather than a sequence — **does not exist in the code**.

This is not a gap g0 would *fail*; it is a gap that makes g0 **meaningless**. Iteration-0
detection would reproduce the loud-source list (that part works — seeding gate 3a PASS,
3/3), and would then reproduce it again at every subsequent iteration. **The campaign would
return INERT with 100% reliability — a verdict about a missing code path, reported as a
statement about physics.** That is the worst available outcome: it is the one failure mode
that looks like a clean negative result.

**`scan_current` is not the missing piece.** `trackB_estimator.py:263` — despite the name,
it returns `lnL_p(n)`, the per-pulsar distance-fringe scan over `EV` at the current theta.
It is the **E-step**, not a sky/frequency source search. Nothing in the repo detects sources
with pulsar terms coherent.

**What building it would mean.** The pulsar-term mask machinery *does* exist as a runtime
argument (`trackB_b1.py:91` `B1Marg(P, data, mask, ...)`), so a coherent detection statistic
is constructible in principle. But it is a **new statistic**, not a wiring change: at each
of 9 024 grid points, marginalise the fringe posterior over each certified pulsar at
**B = 512** evaluations per pulsar — and then null-calibrate *that* statistic, per iteration.
See §0.

---

## 3. THREE FURTHER BLOCKERS, each independently launch-stopping

Confirmed at source, not inferred:

**(a) The loop fits exactly 3 sources; `N_LOUD` is a module constant with ~20 consumers.**
`trackB_b1_core.py:86` `N_LOUD = 3` (`:85` `POP = dict(ncw=16, seed=3000,
population=(3, -13.25, -14.25))`). Consumed to build the fit selector and mc dimensions:
`kindle_loop.py:132` `sel = np.array([nd + 8*k + j for k in range(C.N_LOUD) ...])`, `:138`
`mc_dims`. **A source list that grows past 3 silently fits only the first 3.** Because it is
imported as `C.N_LOUD`, monkeypatching mid-process desyncs `sel` from any already-built
`E`/`FT`. The whole point of AVALANCHE is recruiting from the 13-faint reservoir — i.e.
operating at 4+ sources. This needs a threading refactor, not a parameter.

**(b) The bar is computed ONCE, outside the loop — and a growing list makes it stale.**
`kindle_loop.py:168` `bar = np.maximum(lnK, floor_cell)`, before the iteration, from an `FT`
built once and a scalar `floor_cell`. But CHORUS measures that **growing the source list
moves both terms**: the floor rises with the mix (m0 7.00 → m1e07 11.65 → m2e07 13.73 nat)
and `K_sum` grows ~11× (m0 → m3e07). **Add a source inside the iteration without refitting
the floor and recomputing `K_counted`, and every post-growth certification is scored against
a bar belonging to a smaller population — which reads exactly as self-amplification.** This
is a second independent manufacturing route, distinct from §1, and it points the same way.
The brief's per-iteration floors are the right instinct and address the `floor_cell` half;
`lnK` is the half nobody has flagged.

**(c) Growing the slot count changes the graph shape → the documented OOM.**
`chorus.py:1475-1478`: *"ONE problem + ONE B1Marg per process lifetime (**two evaluator sets
OOM'd an 80 GB card, job 12509837**)"*, plus `B1.B1Marg.T_CHUNK = 8` for wide shapes.
Growth-*within*-a-process is precisely what OOM'd. A list that grows inside an iteration
needs one-shape-per-process with checkpoint/resume handoff — an architecture, not an edit.

**Corollary — the brief's "pure execution" framing does not hold.** KINDLE already
pre-registered the pieces AVALANCHE assumes: the eccentric cells are **E4** (*"a different
problem build; these run through the CHORUS soft loop, not ignite2's"*), the 5+11 structure
axis **"lives in surface.py, NOT ignite.py"** (E2), and the MAP start is **E3 —
NOT BUILT** (`kindle_loop.py:95` `inject_start` **raises `ValueError`** on `start_mode="map"`;
EMBER built `earth_map()` separately). AVALANCHE marries structure × eccentricity × soft loop
× a growing list × a coherent detector × a moving bar. **Six engineering surfaces, of which
KINDLE costed four as "a real build, not a parameter change" — and the sixth does not exist.**

---

## 4. THREE PREMISES IN THE BRIEF THE SOURCES DO NOT SUPPORT

Flagged because two would have propagated into the report as fact:

1. **"the 13-faint sub-threshold reservoir, SNR 1-3"** — the SNR band does not exist.
   `grep -rn "SNR"` over CHORUS + `chorus.py` returns nothing. The census split is pure
   **loudness**: `trackB_b1_core.py:85`, 3 loud at log10_h = −13.25, 13 faint at −14.25,
   exactly 1 dex down. The reservoir is real; the SNR characterisation is invented.
2. **"the 10× lever"** — **stale**. Post-RECUT it is **14.8× (lit) / 12.4× (vlbi)**
   (`CHORUS_mixed_pop.md:224`). RECUT §6 explicitly forbids quoting any cross-cell ratio
   without re-deriving it from the re-cut banks (§10.16(d); the rule EMBER's GATE C0
   enforces). The 10×/11.2×/14.2× figures are all pre-fix.
3. **"e=0.3 ×2 ignites"** — conflates two statistics KINDLE §1.4 explicitly separates.
   CHORUS measures a **count switch-on** (count > 1 bar); **"ignition"** is IGNITE/KINDLE
   vocabulary for **|C|: 0→1 on one trajectory**. Both got called "1.000" once already, and
   that conflation is what produced the retired gain statistic.

Also, for the brief's ladder: **"one e=0.3 member turns the count on" is REFUTED** post-RECUT
(lit 1.57 → **0.70**; vlbi → MARGINAL, dies at floor+bootstrap). The correct pre-registration
is CHORUS §2.1: **one eccentric member → switch-on at e = 0.5; two or more → e = 0.3.**
The brief's igniter set {e=0.7 single, e=0.5 single, e=0.3 ×2} is, correctly, exactly the
switch-on set — but it should be cited from §2.1, not from the refuted headline.

**And the all-circular control is not a negative control here.** The brief expects *"never
ignites, loop inert"*. But per §1 and §2 **every** arm is inert-by-construction in [D], and
per §1 every arm that *does* move manufactures. The control cannot discriminate; it would
return the same INERT as the treatment, for a reason unrelated to eccentricity.

---

## 5. THE VLBI TWO-SIDED READING — mechanism correctly identified, wrong source cited

The brief's pre-registered VLBI reading is **mechanically correct** and worth preserving for
the successor. One correction: it is **not in `RECUT_floors.md`** (read in full; RECUT's only
VLBI content is per-cell floors, and it never names the mechanism). It lives in
**`ANCHOR_realdata_null.md` §7, "THE VLBI PRICE ROW"** (`:322-326`):

> VLBI shrinks σ_d → fewer fringes in the prior window → **smaller `K_counted`** → a lower
> trials bar `ln K` → a pure-noise fluctuation clears layer 1 more easily → more and louder
> offenders → **a higher absolute floor**.

**The price: ≈ +2.9 ± 1.0 nat of null floor at h = −13.25, and nothing measurable elsewhere
in the box** (R2 +2.88 ± 1.01; R0 control +2.59 ± 1.41; R3 +2.76 ± 1.02 — errors of 3–6 nat
at louder h cannot resolve a price this size, *"and I will not quote one"*). Two things
travel with it: it is a **TIER effect, not a realism effect** (the R0 control already carries
it); and it is the **J0437 double edge** — J0437−4715 has the array's smallest trials factor
(ln K = 1.39, the array minimum) and moves +31 under R3, *"robustness to source error and
vulnerability to noise are THE SAME PROPERTY."*

**Cite ANCHOR §7 for the mechanism and RECUT §1.2/§2.4 for the adopted per-cell floors —
not RECUT for both.** The successor's recruitment readout (does the VLBI rung reach deeper
into the reservoir?) remains the right question and is preserved verbatim below.

---

## 6. WHAT THE CHANNEL BUDGET GIVES FOR FREE — the one piece that survives intact

MAGPIE J1's mechanism is **verified and reusable**, and it is the cheapest gate available to
any successor. Definition, exact — `chorus.py:134` + `:327`:

```python
G_ACTIVE_REL = 1e-3                       # :134
active = g >= G_ACTIVE_REL * g.max()      # :327
```

Harmonic n ∈ [1,32] of an eccentric member is **ACTIVE iff g(n,e) ≥ 1e-3 · max_n g(n,e)**
(Peters-Mathews power via `atlas_harmonics.g_n`). Counted by `active_slots()`
(`chorus.py:360-370`): 16 member slots, plus each active appended harmonic slot
(`base = 16 + 31*rank`), minus the member's own slot if its n=2 harmonic went inactive
(at high e the fundamental itself dies).

Re-derived directly from the 780 banked `CHORUS_results/ch_sig_*.npz` — reproduces J1 exactly:

| mix | n_active | grade |
|---|---|---|
| m0e00 | 16 | below |
| m1e03 | 23 | below / MARGINAL |
| **m2e03** | **30** | **CONFIRMED** |
| **m1e05** | **32** | **CONFIRMED** |
| m3e03 | 37 | CONFIRMED |
| m1e07 | 47 | CONFIRMED |

**Every ON cell has n_active ≥ 30; every OFF/marginal cell has n_active ≤ 23. Clean
separation, nothing in between** → the switch is placed at **n_active ≈ 27**, the midpoint of
an empirical gap (bracketed by [23, 30]) — **not a fitted constant**. κ cannot be the
controlling variable: m1e03/m2e03/m3e03 share **identical κ = 2.65** while grades flip
below → CONFIRMED as channels go 23 → 30 → 37. Point-biserial r(n_active) = **+0.53** vs
r(κ) = **+0.26**.

**The gift: `n_active` is a deterministic per-cell function of (n_ecc, e) — pure numpy via
`H.g_n`, no GPU, no likelihood.** A successor can compute the channel budget of **any**
candidate source list **in advance, before spending a second of GPU**. The J1 threshold can
be used to *predict* which recruitments could matter and price the campaign before running it.

MAGPIE's own caveat travels: κ and n_active are collinear along the e-axis (ρ = 0.73 vs 0.81)
and separate **only on the member-count axis at fixed e** — which is exactly where a cascade
campaign lives, so the contrast is available. (EMBER §4.3 correctly declined to read the
budget from loop banks: `n_active` is a per-**cell** census constant, not recoverable from a
per-realisation 19-dim source-model fit vector.)

---

## 7. WHAT WOULD MAKE THE QUESTION WELL-POSED — for the successor brief

The cascade question is **real and unanswered**, and the ingredients the brief cites are
genuinely favourable (per-pulsar upgrades serve every source; a 13-faint reservoir exists;
the channel budget is super-linear with a measured threshold at ~27; known positions buy back
confusion capacity). **The question is worth asking. It cannot be asked this way.**

The honest reformulation follows EMBER's own convention — **anchor, don't search**:

> **ORACLE-ANCHORED RECRUITMENT.** Do not ask the loop to *find* sources (measured
> impossible: 6–12° delivered vs ~1e-5° needle). Instead, anchor each candidate faint source
> at its **true basin** — exactly as `ember_map.py:60-74` anchors the constrained axes — and
> ask whether the certified-coherent array's detection statistic **crosses the per-iteration
> floor**. This measures the cascade's actual physics — *does the array's sensitivity grow as
> pulsars certify?* — while explicitly **not** claiming a blind search could find it.
>
> **Mandatory scope line:** *"conditional on basin anchoring; a statement about the array's
> sensitivity, not about a blind search. The search gap is 6–7 orders (B0.2, `:1450-1462`)
> and is not closed by this result."*

That version is well-posed, safe (anchored counterparts are *true* counterparts → EMBER's
boundary is not crossed → motion is permitted), and answers the question actually being
asked. It still requires four real builds, and they should be **separately gated and
separately costed**, not bundled:

1. **A pulsar-term-coherent detection statistic** (§2) — the novel science. Build and
   null-calibrate it *alone* first; it is the campaign, not a step in one.
2. **A moving bar inside the loop** (§3b) — **both** `floor_cell` *and* `K_counted`. The
   brief's per-iteration floors cover half; `lnK` is the unflagged half.
3. **`N_LOUD` threaded, not constant** (§3a) — ~20 consumers.
4. **One-shape-per-process + resume handoff** (§3c) — the 80 GB OOM is documented.

Do (1) first and in isolation. If the coherent detector's sensitivity gain, measured against
its own honest null, is smaller than the floor rise that adding a source costs (CHORUS: floor
7.00 → 13.73 nat, `K_sum` ~11×), **the cascade is arithmetically dead and (2)–(4) need never
be built.** That is a cheap, decisive, ~few-GPU-hr experiment, and it is the successor's
launch criterion. **AVALANCHE's real deliverable is that it identified this experiment and
declined to spend 400+ GPU-hr on a campaign whose verdict was determined before launch.**

**The forecast read the brief asked for, preserved for the successor:** whether the T=40 and
VLBI rungs flip the verdict remains the right question, and §5's two-sided reading is the
right way to read the VLBI rung (fewer certs at loud cells is *expected*, not an anomaly;
**recruitment depth into the reservoir** is the readout that matters). It cannot be answered
until (1) exists.

---

## 8. INFRASTRUCTURE THE SUCCESSOR MUST INHERIT (verified, unchanged)

- **`--cpus-per-task=8` is part of the seed** — not a performance knob. `NoiseDrawer` builds
  `L_gwb` via CPU `np.linalg.eigh` of a near-degenerate Φ_gwb (cond 2.5e+06, **1 624
  near-degenerate adjacent eigenvalue pairs — exactly half the spectrum**); the LAPACK
  eigenvector basis inside degenerate subspaces **rotates with the BLAS thread count**. 8 →
  bit-identical; 16 → reproducibly shifted (dlnL O(0.1–1) nat everywhere). Verify by L_gwb
  fingerprint on **every** shard (`8548f148b50a5b44` surface / `9fd547b39b02c705` chorus T=30
  / `71677a810cbc7187` b1 bank). **This check is what saved EMBER's 10 drill survivors.**
- **Per-task JAX compile-cache isolation** — `harbor_env.sh:24` honours a pre-set
  `HARBOR_JAXCACHE`; set it **before** sourcing (`e_edge_split.sbatch:23-24`). Eight array
  tasks sharing one cache dir race on autotune temp files → `JaxRuntimeError: NOT_FOUND`;
  **11 realisations guard-killed** in EMBER. Isolation is numerically inert (caches XLA
  executables, not results).
- **Co-tenant GPU OOM** — `--exclude=dgx03` and keep the per-task
  `nvidia-smi --query-compute-apps` device echo. That echo is what made a 50.7 GB foreign
  squat **attributable rather than guessed**.
- **Completeness by READBACK, never by exit code or `ls | wc -l`** — EMBER §2.5: **21
  realisations missing from a sweep in which every log line said `SHARD COMPLETE`**. A
  leftover drill sbatch held shard 0 and TIMEOUT'd; production dispatched tasks 1–7 and
  truthfully reported success. *"`sacct -j 12524934` lists tasks 1–7. There is no task 0.
  That is the whole bug."* Re-derive the completeness table **from the banks**
  (`ember_anatomy.py:26` `readback()`, which re-derives cert from raw:
  `cert_s_raw = (dlnL_start > bar) & (qmax_start > 0.9)` — *"do not trust the banked bool"*).
- **Lean-npz discipline** — bank raw `dlnL`/`lnK`/`qmax`, never verdicts. This is the reason
  every floor in the repo is re-cuttable on CPU with zero GPU (RECUT re-cut all 108 cells
  this way; `kindle_d7_fall.py` re-cut 21 600 banks). Forfeit it and you forfeit that.
- **Two-floor discipline** — gates run at the **banked (old Gumbel)** floor; the science runs
  at the **criterion-v2.2 adopted** floor. Bank both. *"No corrected number was emitted until
  an uncorrected one was reproduced from the same raw columns."*
- **The floor decision function** — `recut_surface.py:101` `adopt(gu, gu_sd, off, zf)`:
  `zf <= 0.20` → Gumbel + `C_SD·beta/sqrt(N)`, tag `"gumbel"`; else → emp q95 + bootstrap sd
  (B=4000, `SEED=20260713`), tag `"emp_q95"`. **The zero-fraction is a REQUIRED column, not a
  caveat.** A Gumbel fitted to a 73%-zero point mass **understated the bar by 53%** and
  manufactured CHORUS's loudest headline.
- **Adopted gate bar** (EMBER §2.2b): **discrete columns and the certified set exact;
  continuous < 1e-6.** Do *not* demand bit-identity — the loop must build a `B1Marg`, and
  `trackB_b1.py:98` `enable_fast_full(all_on)` runs the fast pterm-only residual where banked
  scores used the slow path (drift measured at 3.5e-10).
- **`prior_sigma('vlbi')` still raises** (standing defect, cronus's to fix): `sampler_core`
  accepts only `'lit'`/`'vlbi_all'`, though `build_ignite_problem` **does** construct
  `P["priors"]["vlbi"]` (`ignite.py:265-269`). Any VLBI-rung work trips this.
- **No numpy on the login node** — all work runs under `harbor_py` (apptainer); `HSYMT` is a
  bind-mount synthesised by `harbor_env.sh`, not a real host path.

---

## 9. MEMORY CORRECTION

`kindle-gain-diagnosis` records *"tag v2.2 never cut (HEAD db2075a)"*. **`criterion-v2.2`
now exists and points at `db2075a`.** The memory is stale and should be amended.

---

## STAGE + ADD-LIST

Per convention — **I stage, Matt commits.** No `git commit`, `push`, or `tag` was run.

```
reports/AVALANCHE_cascade.md        (new — this report)
```

**Nothing else.** No `AVALANCHE_results/` directory was created: there are no results, and an
empty results dir would misrepresent a pre-flight stop as an executed campaign.

## STOP

**The campaign is stopped before launch, on three independent grounds, any one of which is
sufficient:**

1. **Its premise is refuted by a banked measurement** (§1). Self-found sources are wrong by
   3–4 orders relative to the certification tolerance, **at zero noise**, and fail
   **confident-wrong** (q_max 0.5–0.99 on a shifted fringe) — passing the q>0.9 layer while
   wrong. Composed with EMBER's boundary (wrong counterpart × motion, p = 0.002), AVALANCHE
   would manufacture a CRITICAL snowball out of noise and could not detect that it had.
2. **Its central feedback edge does not exist in the code** (§2). The F-stat seeder is the
   **fully-decohered limit** by construction (`msd_factory=EarthDelay`). [D] cannot see the
   certified set, so it would return **INERT** with 100% reliability — a missing code path
   reported as physics.
3. **It is ~4× over its own pre-registered cost gate** (§0) — ≈300 GPU-hr of per-iteration
   null refits against an ~80 GPU-hr STOP, and that figure **excludes** the unbuilt detect
   step, whose coherent form multiplies a 9 024-point scan by B = 512 fringe evaluations per
   pulsar per certified pulsar.

**The pre-registered fallback is not proposed**, because it prices only the floors and leaves
§1–§3 untouched: it is a cheaper way to run a campaign that should not run.

**Recommended next action — the successor's launch criterion (§7):** build the
pulsar-term-coherent detection statistic **alone**, null-calibrate it against its own honest
floor, and measure its sensitivity gain against the floor rise that adding a source costs
(CHORUS: 7.00 → 13.73 nat; `K_sum` ~11×). **If the gain is smaller than the cost, the cascade
is arithmetically dead and nothing else need be built.** That is a few-GPU-hr experiment.
Await Matt's decision before any launch.
