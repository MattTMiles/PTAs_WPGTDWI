# FORGE_G_build.md — the SOFT SOURCE SIDE, built and unit-gated (agent FORGE-G, cronus)

*Dev/doc authority. Code STAGED (Matt commits). No git commit/push/tag from here. This box is
RECOMPUTED-UNSAFE for noise by design (no canonical `L_gwb` bank), so everything below is
ZERO-NOISE / ASIMOV: no campaign, no noisy realisations. The pure-logic gates ran on cronus and
PASS; the zero-noise GPU gates are staged for the sanctioned lane (jug-gpu).*

---

## 0. WHAT THIS BUILDS AND WHY

The spec-§3 loop (`CW_transition/trackB_estimator_spec.md` §3, the fringe-posterior-weighted
M-step) is soft on the **pulsar** side and hard on the **source** side. Every pulsar carries a
fringe posterior `q_p(n)` and the M-step is the q-weighted
`Q(theta) = Σ_p Σ_n q_p(n) lnL(theta, L_p(n))` (`trackB_b1.B1Marg`); a hard-locked fringe is
exactly the violation that caused the S8.1.2 cascade. But the **source** side is a bare point
vector: `trackB_b1_referendum.TargetedMarg` moves the loud sources' `(f, mc)` and **pins every
faint reservoir source at truth**. That asymmetry is the one hard edge left in the loop.

SPARK-2 and SPARK-3 measured what removing it is worth and named the law it must obey:

- **SPARK-2** (`reports/SPARK2_second_arrow.md`): the faint reservoir **UNMODELLED erases
  certification** (offender 0.000 vs 4.435 truth-modelled); the deconfusion **trials cost is
  +0.578 nat and SATURATES** (sublinear, bounded).
- **SPARK-3** (`reports/SPARK3_second_arrow.md`): the verdict is **STRADDLED** — arrow 2 is
  **MODEL-QUALITY-LIMITED**. Under an **OPTIMISTIC** (conditional-width) reservoir model
  deconfusion recruits (11 crossings); under a **PESSIMISTIC** (**PRIOR**-width) model a
  badly-constrained reservoir **actively raises the certification floor and suppresses the count**
  (floor 118 → 744 at one unit). §5.0: *a prior-width model is WORSE than no model at all.*
  **THE LAW: a soft component may feed the joint template only when its posterior is demonstrably
  tighter than prior width by a declared factor; below that it is carried but not fed.**

FORGE-G builds that law as tracked machinery: every source carries a **posterior**, the joint
consumes **soft components + pulsar fringe posteriors symmetrically**, and a **model-quality
gate** decides feed-back. This is the build only — no verdict, and GLACIER stays UNLICENSED
(SPARK-3 §5.3: STRADDLED is not EDGE-POSITIVE).

---

## 1. THE ONE CORE EDIT — SMASK, the source-side twin of PMASK (`trackB_b1_core.py`)

The pulsar side already hands the joint a per-pulsar coherence weight `m_p` (`PMASK`,
`MaskedDelay`: `delay_p = d_earth + m_p·(d_full − d_earth)`). The source side had no analogue: the
template is `Σ_s delay_s` over all sources, unconditionally. FORGE-G adds the symmetric
**per-source feed weight `SMASK`**:

```
delay_p = d_earth + m_p·(d_full − d_earth),   d_full/earth = Σ_s w_s · delay_s(source s)
```

- `w_s = 1` → source `s` contributes its full delay; `w_s = 0` → source `s` is **absent from the
  template, exactly**, via `jnp.where(w_s>0, w_s·delay_s, 0)` — SELECT, not multiply, so a carried
  source's waveform never enters the sum. Multiplying (`0·NaN = NaN`) would let one carried source
  in the coalescence region (SPARK-3 §2.1) NaN the whole array; selection makes `w_s=0` TRUE
  absence. CPU-validated bit-exact (§3.2 [4]).
- **ABSENT `SMASK` key → `w_s = 1` for every source** (`where(True,1·x,0) == x`, bit-exact), so
  **every pre-FORGE_G caller — FORGE, SPARK, IGNITE, CHORUS, SURFACE, softcascade, referendum — is
  bit-identical** (CPU-validated, §3.2 [1][2]). SMASK is opt-in.

Threaded through exactly as `enable_fast_full` already threads state: `B1Split.smask` (default
`None`), injected in `B1Split._params`, and `B1Split.set_smask(w)` / `B1Marg.set_smask(w)` clear
the jitted per-pulsar evaluators (they close over `smask` through `_params` at trace time). A
non-`None` smask **forces the masked residual path** — the pterm-only fast-full residual
(`MultiSourceDelay`) does not honour SMASK, so `set_smask` disables fast-full while a smask is
live. `B1Split.a_const` is the frozen `ldN + rednoise logdet` and is smask-independent; it is
computed once at `None` and reused (never recomputed under a smask via the smask-blind
`amo["logL"]`).

**This is the only change to a tracked, gated file.** It is additive and default-inert; `TE.phi_bounds`,
`deterministic.py`, `ember_map.py` are untouched.

---

## 2. THE SOFT SOURCE MODULE — `trackB_b1_softsource.py` (new)

| brief | object | what it is |
|---|---|---|
| **B1** | `SoftComponent` | one source carrying a POSTERIOR (mean + per-axis σ; Gaussian/Laplace, grid hook for multimodal). `at_prior_width()` = a reservoir component started at PRIOR width. |
| **B1** | `SoftSourceState` | the SET of components. `assemble(gate) → (src, smask)` hands the joint the template means + the gated feed weights — the source-side twin of the pulsar q-weights. Consumes soft sources + pulsar fringe posteriors symmetrically. |
| **B2** | `ModelQualityGate` | SPARK-3's law. `feeds(comp)` iff the binding concentration ratio (posterior/prior width, worst free axis) is **strictly** `< ratio_thresh`. |
| **B3** | `concentration_trace`, `scoreboard` | per-component concentration tracked per iteration; the v2.2 detection criterion runs ONLY as the end-of-cycle scoreboard, never inside the loop. NO binary `detected` flag in the loop. |
| **B4** | `condition_component`, `condition_state` | EM hard pins: one flag `None | 'sky' | 'sky+period'`. |
| — | `SoftCycle` | the zero-noise E-step + gated source hand-off (GPU; wraps `B1Marg`). |

### 2.1 The model-quality gate (B2) — the parameterisation, stated exactly

Concentration ratio `r = σ_post / σ_prior` per axis (B3 statistic), reduced to the **worst free
axis** (a component is only as well-modelled as its worst free axis — SPARK-3 §5.0: prior width on
*any* key axis lets the fit accommodate noise). Pinned axes (B4) are excluded (known exactly,
`r = 0`). The gate feeds iff `r < ratio_thresh`, **strict**:

| `ratio_thresh` | behaviour | why |
|---|---|---|
| **0** | **UNMODELLED** — nothing feeds | `r < 0` never true (strict; a delta posterior `r=0` still needs a positive threshold) |
| **∞** | **ALWAYS-FED** — everything feeds | `r < ∞` always true |
| `1/F` (default) | feed iff posterior is `> F×` tighter than prior | SPARK-3's declared factor |

`ratio_thresh = 1/(declared factor)`. The **default `SPARK3_GATE_RATIO = 0.5` (factor 2×)** is a
DECLARED parameter, provisional: SPARK-3 §5.0 shows the PESSIMISTIC bound is prior width (`r = 1.0`,
suppresses) and the OPTIMISTIC bound is the conditional width — the exact cross-over is the
marginal width *between* the two Fisher bounds, which SPARK-2 measured is unaffordable to compute.
**The machinery's content is the two LIMITS reproducing UNMODELLED / ALWAYS-FED; the intermediate
default is the one knob a campaign re-derives from its own Fisher.** Prior widths are the
**generative population** box (`stagec_fisher`: `log10_mc U(8.5,9.5)`, `log10_fgw U(-8.0,-7.5)`,
…), σ = half/√3 — **NOT `TE.phi_bounds`**, which is the estimator SEARCH box (SPARK-3 §2.2's
correction).

### 2.2 "Fed / carried / pinned" — three distinct states

- **fed** (`w_s = 1`): in the template, deconfuses, and a free parameter of the M-step.
- **carried, not fed** (`w_s = 0`): **absent from the template** but **present in the trials
  budget** (SPARK-2 §2's `dL = min` over the modelled list; +0.578 nat, saturating). SPARK-3's
  "carried but not fed".
- **pinned** (B4): fed, but the pinned axes are known exactly and **drop from the joint free
  dimension** (`SoftSourceState.free_dim`).

---

## 3. THE GATES (B5)

Run: `python trackB_b1_softsource.py gate` (jug-gpu) or `--cpu-only` (logic only, anywhere).

### 3.1 Pure-logic gates — RAN ON CRONUS, **ALL PASS**

```
-- Bg2: model-quality gate limits --
   threshold-0   fed 0/3  (unmodelled)   PASS
   threshold-inf fed 3/3  (always-fed)   PASS
   default r<0.5 (factor 2x) fed 2/3     PASS      (feeds tight+mid, CARRIES the prior-width one)
   monotone fed-count vs ratio_thresh [0,0,1,2,2,2,3,3]  PASS
   delta posterior fed@0=0, fed@default=1            PASS      (strict inequality)
-- Bg4: EM conditioning drops the pinned block --
   free_dim base 24                                  PASS
   free_dim one comp sky-pinned      22 (−2)         PASS      (cos_gwtheta+gwphi drop)
   free_dim one comp sky+period      21 (−3)         PASS      (+log10_fgw)
   free_dim one comp carried-not-fed 16              PASS      (not-fed → 0 free params)
```

### 3.2 SMASK core threading — RAN ON CRONUS (CPU), **ALL PASS**

The full B5 GPU gates need jug-gpu (below), but the load-bearing SMASK claims run cheaply on the
single-eval `B1Split.lnL` path (no 512-fringe sweep) and were validated on cronus CPU:

```
[1] amo.logL == split.lnL(smask None)                          |d| = 0.00e+00   incumbent unaffected
[2] split.lnL(smask ONES) == split.lnL(smask None)             |d| = 0.00e+00   BIT-EXACT (Bg5a core)
[3] split.lnL(loud-only)  gap vs full = -95.71 nat                              reservoir absent
[4] reservoir params -> JUNK, loud-only smask, lnL unchanged   |d| = 0.00e+00   BIT-EXACT ABSENT
[5] restore smask None == original                             |d| = 0.00e+00   cache-clear correct
```

`[2]` is Bg5a's kernel (`1.0·x == x`): always-fed smask≡1 is bit-identical to the incumbent.
`[4]` caught and drove a **fix**: `w_s·delay_s` computes the not-fed source's waveform then scales
by 0, and `0·NaN = NaN` — one carried source in the coalescence region (SPARK-3 §2.1's exact
mechanism) would poison the whole array. Replaced by `jnp.where(w>0, w·delay, 0)` (SELECT, not
multiply), so `w_s == 0` is TRUE absence: a carried source cannot NaN the template. `[2]` stays
bit-exact under the fix (`where(True, 1·x, 0) == x`).

### 3.3 Zero-noise B5 gates — RAN ON CRONUS (CPU, full 116×512 fringe machinery), **ALL PASS**

```
-- Bg5a: always-fed + truth-pinned == incumbent spec-§3 marginal loop --
   |lnL_marg(always-fed) - lnL_marg(incumbent)| = 0.000e+00   BIT-EXACT
   max|LNL(always-fed)   - LNL(incumbent)|      = 0.000e+00   BIT-EXACT
   max|q_max(always-fed) - q_max(incumbent)|    = 0.000e+00   BIT-EXACT
   (fast-full vs masked = 0.000e+00 < 1e-8)                         [Bg5a] PASS
-- Bg5b: zero-noise fixed point (start at truth, one cycle moves nothing) --
   max|src(fp) - truth|      = 0.000e+00 (<1e-6)
   census ceiling max|diff|  = 0.0004    (<0.02)                    [Bg5b] PASS
-- Bg5c: gate limits on the joint + w_s==0 bit-exact absence --
   threshold-0 smask all-zero (fully unmodelled)          True
   default smask == loud-only (reservoir carried, not fed) True
   |lnL_marg(reservoir@truth) - (@junk)| with w_res==0 = 0.000e+00  [Bg5c] PASS
```

| gate | proves |
|---|---|
| **Bg5a** | **ALWAYS-FED + truth-pinned == the incumbent spec-§3 marginal loop, BIT-EXACT** — on the full fringe posterior (lnL_marg, per-pulsar LNL, q_max), smask≡1 is `where(True,1·x,0) == x`; fast-full certified `== masked` separately. This is the mission's "spec-§3 loop reproduced bit-comparably at threshold-∞ + truth-pinned (the old behaviour as a limit)". |
| **Bg5b** | **zero-noise FIXED POINT** — start at truth, the source hand-off is stationary and the census ceilings reproduce. |
| **Bg5c** | gate limits on the real joint: threshold-0 → fully unmodelled (all-zero smask), default → reservoir carried-not-fed (loud-only smask); **w_s==0 is bit-exact absence** (a not-fed source's params are irrelevant even at the full-fringe level). |

**These ran on cronus CPU** (`JAX_PLATFORMS=cpu`; the box's GPU DNN library fails to init — the
broken-box condition — so the GPU path auto-defers, and `--cpu-only` runs the logic gates alone).
Correctness is proven here; a **jug-gpu re-run is owed only for the sanctioned-lane numbers and
speed**, not for the verdict — every gate is exact (`0.000e+00`) and lane-independent.

---

## 4. WHAT THIS DOES NOT DO (scope, travelling)

- **No campaign, no noisy realisations, no verdict.** Zero-noise/Asimov machinery + gates only.
- **GLACIER stays UNLICENSED** — SPARK-3's verdict is STRADDLED, not EDGE-POSITIVE (SPARK-3 §5.3).
- **The M-step optimiser is not re-built here.** `SoftCycle` provides the gated hand-off and the
  fixed-point property; the Newton step on `lnL_marg` is `trackB_b1_referendum`'s machinery.
- **The grid (multimodal) posterior representation is a HOOK**, not exercised — the Gaussian mean
  is the fed value in this build (Laplace). A multimodal axis (e.g. the `log10_fgw` periodogram,
  S8.5.3c) is where the grid rep earns its place; flagged, not built.
- **The default gate ratio is provisional** (§2.1). It is not a measurement.

---

## 5. ADD-LIST (staged; Matt commits)

```
CW_transition/trackB_b1_softsource.py   new   B1–B5: soft source state, model-quality gate,
                                              scoreboard, EM conditioning, zero-noise gates
reports/FORGE_G_build.md                new   this report
```

**Modified (tracked, gated — additive + default-inert, §1):**
```
CW_transition/trackB_b1_core.py         SMASK: per-source feed weight in MaskedDelay + B1Split
                                        (set_smask, _params inject). ABSENT key -> bit-identical
                                        to every existing caller.
CW_transition/trackB_b1.py              B1Marg.set_smask passthrough (clears batched evaluators)
```

Nothing else is touched. The SMASK edit is gated inert by Bg5a (always-fed smask≡1 == incumbent,
bit-exact) and by the existing `trackB_b1_core.mode_gate` G0–G8 (which every non-smask caller
still passes, smask defaulting to `None`).

---

## 6. STOP

The soft source side is built and unit-gated as tracked machinery: **SMASK** (the source-side twin
of PMASK, default-inert), **SoftComponent/SoftSourceState** (every source a posterior, reservoir at
prior width, symmetric hand-off), the **ModelQualityGate** (SPARK-3's law; the two limits gated),
the **soft detection scoreboard** (no in-loop flag), and the **EM conditioning pins** (position
block drops from the joint dimension). **ALL gates PASS on cronus** — logic (B2/B4), the SMASK
core threading (bit-exact inert + true absence), and the full zero-noise B5 suite (Bg5a/b/c, all
`0.000e+00` on the real 116×512 fringe machinery, CPU). Code staged, add-list above. **Awaiting
Matt: commit; a jug-gpu re-run is owed only for the sanctioned-lane numbers/speed, not the
verdict.**
