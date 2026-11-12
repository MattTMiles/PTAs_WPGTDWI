# FORGE_B_belief.md — the RUNNING-BELIEF M-STEP, built and gated (agent FORGE-B, cronus)

**third-arm machinery for the frozen-vs-live rider — ACCRE pulls and fans per the
pre-registration; fan gated behind the frozen-arm readback + this landing.**

*Dev/doc authority. Code STAGED (Matt commits). No git commit/push/tag from here. This box is
RECOMPUTED-UNSAFE for noise by design (no canonical `L_gwb` bank), so everything below is
ZERO-NOISE / ASIMOV or analytic fixture: no campaign, no noisy realisations. **All gates
g1–g5 + B3 PASS on cronus** (`CW_transition/logs/forge_b_belief_gates.log`).*

**Matt's directive, the arm's charter, banked verbatim:** *"give sources the same
running-belief treatment the pulsar fringes already enjoy; tighten only on earned
information; wander bounded by belief width by construction."*

---

## 0. STATE FOUND — and one STEP-0 gap Matt should know about before reading further

`git pull --ff-only` → **already up to date** (HEAD `5461142`, 2026-07-26).

**The three capstone sections the tasking sent me to are NOT IN THIS TREE.**
`GLACIER_capstone.md` ends at **§S4.14** (the dead-frontier defect, 2026-07-24). There is no
§S4.15, §S4.20 or §S4.24, and a repo-wide search for the terms they carry —
`running-belief`, `frozen-vs-live`, `frontier-v2`, `data-support`, `dual-ladder`, `wander` —
returns **nothing** in any tracked or untracked file. Likewise `reports/TURBINE_speed.md`
never landed (FORGE-G2 §0 hit the same thing). So I could not read:

| owed | what I did instead |
|---|---|
| §S4.15 M-step policy decision 2(d); the frozen-vs-live **pre-registration** | built all THREE arms in the driver (`--mstep-policy live\|frozen\|belief`) rather than assume a `frozen` arm existed to extend; `live` is the default and is byte-identical to the incumbent |
| §S4.24 dual-ladder verdicts + wander numbers | used the tasking's own numbers verbatim — **0.01–0.44 dex at ~0.2-nat curvature** — as g3's fixture. They are the only quantitative input the arm needed |
| §S4.20 frontier-v2, the data-support term | frontier-v2's exact term is unavailable, so I **declared** the data-support delta-lnL as the member's own present-vs-absent lnL at the feed point (`glacier_loop.data_support`, own-term-live per S4.14) and flagged it as declared, not adopted |

**Nothing in the build is blocked by the gap** — the machinery, the guard and all five gates
are specified by the tasking itself. But **the frontier-v2 data-support term is the one place
where I substituted my own definition for a pre-registered one.** If §S4.20's term differs,
it is a one-function change (`data_support`) plus a re-run of the feed-time numbers in §5;
nothing else in the arm depends on which delta-lnL is used, only that one exists.

Code found: `trackB_b1_softsource.py` (FORGE-G's `SoftComponent`/`SoftSourceState`/
`ModelQualityGate`) and `trackB_b1_core/trackB_b1` on the FORGE-G2 runtime-SMASK path, both as
their reports describe. `glacier_loop.mstep_quadratic` is the live M-step, unchanged since.
No prior FORGE-B artifacts anywhere.

---

## 1. THE THIRD ARM — what the loop was missing

The pulsar side has carried a running belief since spec-§3: every pulsar holds a fringe
posterior `q_p(n)` that **persists and is updated**, never re-derived. FORGE-G made the source
side soft in the *set* sense — each source carries a width, the `ModelQualityGate` decides
feed-back — but the **M-step itself** stayed a fresh re-maximisation: `mstep_quadratic` walks
from the current template value to the local likelihood mode, and nothing remembers where the
belief was or how much of it was earned.

The dual ladder is what makes that a defect rather than a taste: at **~0.2-nat class
curvature the refit wanders 0.01–0.44 dex per iteration**. At that curvature the likelihood
*cannot distinguish the wander from the truth*, so the live arm's motion is unearned. The
frozen arm has the opposite failure: it cannot move even on information that *is* earned.

**B1 — the belief state** (`Belief`, `BeliefBook`). Each fed component carries `centre` (the
value the template is fed), `anchor` (the centre **of record**), and `width` (per parameter).
It is **born at feed time from the frontier's own evidence**:

```
width_j = min( box_j , F_ii_j^(-1/2) * sqrt(max(1, DLNL_SUPPORT_REF / dlnL_support)) )
```

the curvature at the feed point (the same conditional width the frontier gated on, own-term-live
per S4.14) and the **data-support delta-lnL** — present-vs-absent for that member at that state.
The support term is the guard against a member that crosses on sharp-but-shallow curvature: fed
on thin evidence, it starts proportionately **wider**, capped at the prior box. The belief
travels in the component state (`SoftComponent.belief`, `SoftSourceState.attach_beliefs`),
**checkpoints** with it (`BeliefBook.pack/unpack`) and **banks per iteration**.

**B2 — the update rule.** Refits become **MAP against the belief**. `BeliefMStep.step` runs the
incumbent's coordinate-quadratic sweep unchanged — the arm changes *what is done with* the local
quadratic, never how it is measured — then combines the likelihood centre with the belief by
precision weighting and clips to the belief band. See §2: the guard is the whole subtlety.

**B3 — graduation.** `GraduationLadder` reads belief width **live**: `carried → fed (r<0.5) →
refine-eligible (r<0.2) → cert-eligible (r<0.1)`, same strict inequality and same worst-axis
reduction as `ModelQualityGate`, and the **first rung IS `SPARK3_GATE_RATIO`** — the ladder
extends the incumbent frontier, it does not move it. `GraduationLedger` banks the
belief-width trajectories as a first-class readout (**the sources-resolving-confidently
curve**), with the M-step-axis ratio beside it as a labelled diagnostic (§4).
"Eligible" is load-bearing: the ladder never certifies. The v2.2 criterion still runs **only**
as the end-of-cycle scoreboard (FORGE-G's no-in-loop-flag rule).

---

## 2. THE IMPLEMENTATION NOTE ON THE GUARD (the one load-bearing subtlety)

Iterations re-read the **same data**. A refit that narrows the belief on that re-read counts the
same photons twice and the belief collapses to a confidence the data never supplied — the exact
failure the arm exists to stop. The guard has **two halves**, and only the first is obvious.

**Half one — the width may narrow only on a licensed event.** `InformationLedger` is the
explicit detector for the four kinds that qualify, and nothing else qualifies:

| event | why it is genuinely new information |
|---|---|
| `fed_set` | what the joint models changed, so every other member's conditional information changed |
| `cert_set` | the E-step's PM rows cohere a different set, so the fringe posteriors carry different information |
| `background` | the drain moved beyond its own band (`\|dA_bg\| > K·sigma`), so the residual the sources are fitted against moved |
| `data_epoch` | new epochs — the only unambiguous case |

Another sweep, another iteration, a re-read, a **resume** — none are news. Licensing is
token-based: `observe` queues events, one M-step call `spend`s **all** queued events and is
thereby licensed **at most one** tighten. Multiple events in an iteration license one tighten,
not several. The **first** `observe` seeds the record and fires nothing — that information was
already spent by `Belief.from_frontier` at feed time. Every tighten is logged with the event
that licensed it (`Belief.tighten_log`, banked as `belief_event_*`).

**Half two — the MAP prior is ANCHORED, and this is the half that is easy to get wrong.** The
prior in the MAP update is `N(anchor, w_b²)` where `anchor` is the centre **of record**, which
moves *only* when a tighten is earned — **not** last iteration's working centre. Had the prior
chased its own previous posterior, repeated iterations at unchanged data would ratchet the
centre to the likelihood mode: the width would be honestly frozen while the centre quietly
completed the same walk, and the arm would decay into the live arm with extra bookkeeping. That
is double-counting through the centre rather than the width, it passes any width-only assertion,
and it is invisible in the banked widths. Anchored, repeated iterations converge to a fixed point
inside the band and stay there.

**The bound is then structural, not statistical.** The returned centre is hard-clipped to
`anchor ± WANDER_K·width`, so `|centre − anchor| ≤ K·width` holds for **any** likelihood,
including a pathological one, at **any** iteration count — "bounded by belief width **by
construction**", not in expectation. `n_clipped` banks how often the bound actually bound.

**Resume is part of the guard.** A resumed cell unpacks the banked belief and **never
re-derives** one; re-deriving would be a free tighten at unchanged information and would make a
resumed cell differ from a straight one (the drill's whole property). On resume the ledger is
brought up to date and its licences are then **discarded** — the tightens those events licensed
are already baked into the resumed belief. Resuming the belief arm from a pre-belief checkpoint
raises `CampaignStop` rather than silently re-deriving.

**Why the two limits are exact rather than close.** They are explicit branches, not `1/inf`
arithmetic: `width → inf` returns the likelihood centre untouched (so the whole numerical
trajectory is the incumbent's, bit-for-bit), and `width → 0` **skips the slot entirely** — a
delta belief has nothing to learn — so frozen reproduces including `n_eval == 0`.

---

## 3. THE GATES — ALL PASS (cronus, `CW_transition/logs/forge_b_belief_gates.log`)

Run: `python trackB_b1_belief.py gate` (add `--cpu-only` to skip g5).

### g1 LIMITS — the two existing arms are literal limits of this one

```
[A] live limit  max|src - glacier_loop.mstep_quadratic| = 0.000e+00, max|widths| = 0.000e+00, n_eval 24 == 24
[B] live limit  max|src - glacier_loop.mstep_quadratic| = 0.000e+00, max|widths| = 0.000e+00, n_eval 36 == 36
[C] live limit  max|src - glacier_loop.mstep_quadratic| = 0.000e+00, max|widths| = 0.000e+00, n_eval 24 == 24
[A/B/C] frozen limit max|src - feed point| = 0.000e+00, n_eval 0 (== 0, frozen costs nothing)
```

The live reference is **`glacier_loop.mstep_quadratic` itself**, imported and called — not a
lookalike. (The module also carries a verbatim transcription for when the driver is not
importable, and the gate asserts transcription == real function at `0.000e+00` before using it.)

### g2 THE GUARD — zero tightening without news, asserted exactly

```
6 quiet iterations: width array EXACTLY unchanged: True; n_tighten [0, 0] (all 0): True
re-centring remained free: max|centre - feed| = 1.500e-02 (> 0): True; stayed inside the belief band: True
event fired: ['cert_set'] -> tightens [1, 1] (exactly 1 each): True
licensed widths strictly narrowed on the M-step axes: True (unmeasured axes untouched: True)
every tighten logged with its event: True
5 further quiet iterations: width EXACTLY unchanged, n_tighten still 1: True
all four event kinds detected, in order: ['background', 'fed_set', 'data_epoch'] -> True
```

Equality is `np.array_equal` on the width arrays — **exact**, not tolerance. Re-centring is
shown still free in the same iterations (1.5e-2 of motion), which is what separates this from
the frozen arm.

### g3 WANDER BOUND — the walk-suppression demonstration

Fixture at the measured hazard's curvature class: `F` from `dlnL = ½F·dx²` at
`(0.44 dex, 0.2 nat)` → σ = 0.696 dex. The wander driver is a **deterministic seeded** random
walk of the analytic surface's mode across 120 iterations — a fixture, not a data realisation
(no `L_gwb`, no noise draw, nothing recomputed).

```
belief width 0.2887 dex, bound (K=1) 0.1443 dex
BELIEF arm: max displacement over 120 iters = 0.0087 dex (<= bound: True); tightens 0 (no news)
LIVE control (no belief, same fixture): max displacement = 1.6969 dex -- exceeds the bound
walk suppressed by x194.2
```

The live control is a **negative control in the FORGE-G2 SG5 sense**: it proves the fixture
bites. Without it "the displacement was small" would be unfalsifiable.

Note the honest reading of the belief width here: at 0.2-nat curvature the Laplace width
(0.696 dex) is **worse than the prior box**, so the belief is box-limited (0.2887 dex). The arm
does not manufacture a tight belief at flat curvature — it bounds the wander with the width it
actually has.

### g4 DEFAULT-INERT

```
[A] vs glacier_loop.mstep_quadratic: max|src| 0.000e+00, max|widths| 0.000e+00, n_eval 24 == 24
[B] ... 0.000e+00 / 0.000e+00, 48 == 48      [C] ... 0.000e+00 / 0.000e+00, 12 == 12
[D-nonconcave] ... 0.000e+00 / 0.000e+00, 12 == 12      (the non-concave branch too)
[empty book] a BeliefBook with no belief for the fed slots is inert: True
checkpoint round-trip exact: True; pre-belief checkpoint -> None (inert): True
```

### g5 RUNTIME-SMASK COMPATIBILITY — with belief state live

```
warm evaluator with belief live: fixed-pattern 112.70 ms, FLIP 109.89 ms (median of 20)
FLIP OVERHEAD = -2.80 ms (gate < 10 ms, lane-independent) -> PASS
jit cache size across flips: 1 -> 1 (unchanged -> NO RECOMPILE)
belief state stayed live and banked across every flip (20 rows, 3 beliefs): True
```

**Read this one carefully.** The no-recompile claim — the load-bearing one — is **asserted on
the jit cache size, not inferred from timing**, and passes exactly. The flip **overhead** is
`−2.80 ms`, i.e. zero within noise: flipping the pattern costs nothing beyond evaluating.

**The absolute `< 10 ms` budget is NOT met on this box and is OWED ON jug-gpu.** FORGE-G2's
10 ms is a jug-gpu number (6.35 ms there); on cronus **CPU** the per-pulsar evaluator alone
costs 112.7 ms, so no implementation — including the untouched incumbent — could meet it here.
The gate reports both numbers and gates the lane-independent one; the absolute sub-check is
declared deferred, not quietly re-scoped. The structural reason it will pass on the sanctioned
lane: **the belief is entirely numpy-side and never enters a trace**, so it cannot recompile.

---

## 4. A DEFECT THE BUILD CAUGHT — the graduation ladder's gate axes

Exercised against the **real b1 machinery** (zero-noise Asimov, `data_support` +
`fisher_conditional` → `Belief.from_frontier`):

| slot | dlnL_support | belief ratio, ALL axes | ratio, M-step axes only | frontier ratio |
|---|---|---|---|---|
| 0 (loud) | +680.350 | 0.0828 → **cert-eligible** | 0.0000 | 0.0828 |
| 3 (reservoir) | +1.614 | 0.9482 → **carried** | 0.0015 | 0.9482 |

My first wiring keyed the ladder on the M-step axes (`log10_fgw`, `log10_mc`) — the axes the
loop actually moves. On real numbers that reads **cert-eligible for a member the frontier scores
at 0.948**: believed to 1e-5 in frequency with a **prior-box sky**. That is SPARK-3 §5.0
exactly — prior width on *any* key axis lets the fit accommodate noise — and it is the S4.14
failure shape: a gate reading green off an instrument that cannot see the thing it is gating.

**Fixed:** gate axes default to **all axes**, matching the frontier's worst-free-axis reduction.
With the boxes equal (asserted in the driver — a divergence raises `CampaignStop`), a belief
born at support ≥ `DLNL_SUPPORT_REF` reads **exactly the frontier's own ratio at feed time**:
0.0828 and 0.9482 above are identical to `max(sig_opt/sig_pes)`. A **standing control** in the
B3 gate asserts both sides — a lopsided belief (1e-4 on the M-step axes, box elsewhere) must
read `carried` on the all-axis ladder and `cert-eligible` on the narrowed one — so re-narrowing
`gate_axes` fails the gate instead of reading green.

**The consequence, stated rather than engineered around:** while the M-step surface moves only
2 of 8 axes, the **deep rungs are structurally unreachable** for any member whose other six axes
are box-limited. That is a real property of the current loop, not of this arm. The M-step-axis
ratio is banked as a labelled **diagnostic** column (`grad_width_ratio_diag`) beside the binding
one, so the two are visible side by side and neither is silently substituted for the other.

---

## 5. WHAT THIS DOES NOT DO (scope, travelling)

- **No campaign, no noisy realisations, no verdict.** The fan is Matt's call per the
  pre-registration; the `HOLDS_CLEARED` refusal still fires for the belief arm (verified).
- **The driver's belief path has not been run end-to-end**, because it cannot be here: `run_cell`
  needs the canonical `L_gwb` bank and **cronus has none by design** (`b1_L_gwb_T30.npz` absent;
  `NoiseDrawer(strict=True)` would crash, correctly). What IS exercised on cronus: every unit of
  the arm (g1–g5), the driver's import/CLI/dispatch, and the feed-time evidence path
  (`data_support` + `fisher_conditional` → `Belief`) on the real b1 machinery (§4). **A
  driver-level three-arm cell is owed on ACCRE** — the natural shape is the resume drill
  (`drillcell` ×2 + `drillcmp`) at `--mstep-policy belief`, which also certifies the belief's
  resume identity.
- **The declared parameters are declared, not measured** (`WANDER_K = 1.0`,
  `DLNL_SUPPORT_REF = 1.0` nat, ladder `(0.5, 0.2, 0.1)`, `BG_BAND_K = 1.0`) — the
  `SPARK3_GATE_RATIO` tradition. The machinery's content is the LIMITS; these are the knobs a
  campaign re-derives from its own Fisher.
- **The data-support term is my definition, not §S4.20's** (§0). One function if it differs.
- **The frozen arm is my reading of "FROZEN-after-feed"** — fed members never move from their
  feed point — because the pre-registration is not in this tree. If §S4.15's frozen arm differs,
  g1's frozen limit is the thing to re-point.
- **No multimodal belief.** The Gaussian summary is the belief; FORGE-G's grid hook is still a
  hook. A multimodal axis (the `log10_fgw` periodogram, S8.5.3c) is where it would earn its place.

---

## 6. ADD-LIST (staged; Matt commits)

```
CW_transition/trackB_b1_belief.py            new   B1-B3 + g1-g5: Belief/BeliefBook,
                                                   InformationLedger, BeliefMStep,
                                                   GraduationLadder/Ledger, the gate battery
CW_transition/logs/forge_b_belief_gates.log  new   the full gate log (ALL PASS)
reports/FORGE_B_belief.md                    new   this report
```

**Modified (tracked, gated — additive + default-inert):**

```
CW_transition/trackB_b1_softsource.py   +63/-3    SoftComponent.belief (default None) +
                                                  set_belief/sync_from_belief; SoftSourceState
                                                  attach_beliefs/sync_from_beliefs/
                                                  belief_widths/n_belief. FORGE-G's Bg2/Bg4
                                                  gates RE-RUN and PASS unchanged.
hpc_harbor/glacier/glacier_loop.py      +188/-10  MSTEP_POLICIES ("live" default) + --mstep-policy;
                                                  data_support(); the belief-arm branch of the
                                                  M-step; feed-time belief at promote; belief
                                                  resume + bank block; stem tag _m<policy>.
                                                  mstep_quadratic itself UNTOUCHED.
```

**Inertness of the driver edit:** `mstep_policy="live"` (the default) never constructs a
`BeliefBook`, never calls the belief module, takes the `if book is None:` branch into the
unmodified `mstep_quadratic`, adds **no bank key**, and leaves the stem unchanged — so every
existing cell, null, ladder rung and drill is byte-identical. Non-live arms tag the stem
`_m<policy>` and therefore can never collide with, or resume from, an incumbent bank.

## 7. SIZE GLANCE

```
CW_transition/trackB_b1_belief.py       1259 lines (new; ~430 of them the g1-g5 battery)
CW_transition/trackB_b1_softsource.py    638 lines (was 579; +63/-3)
hpc_harbor/glacier/glacier_loop.py      1038 lines (was 872; +188/-10)
CW_transition/logs/forge_b_belief_gates.log  50 lines
diffstat (tracked files): 2 files changed, 238 insertions(+), 13 deletions(-)
```

---

## 8. STOP

The **third arm is built and unit-gated**: a persistent per-component belief born at feed time
from the frontier's own evidence, an M-step that **updates** it instead of re-maximising, the
**double-counting guard** in both its halves (width licensed only by an explicit event detector;
the MAP prior **anchored** so the centre cannot complete the walk the width was forbidden), the
wander bounded by belief width **by construction**, and a graduation ladder reading belief width
live with the trajectories banked as a first-class curve. **The two existing arms are exact
limits of it** — live at width → ∞ and frozen at width → 0, both `0.000e+00` against
`glacier_loop.mstep_quadratic` itself. **All gates PASS on cronus**; the only sub-check not met
here is g5's absolute-millisecond budget, which is a jug-gpu number this CPU lane cannot carry
and is declared owed, with the no-recompile claim asserted exactly on the jit cache.

**Awaiting Matt:** commit; the STEP-0 gap in §0 (three capstone sections and the frontier-v2
data-support term absent from this tree); and the ACCRE-side driver exercise (§5) — the resume
drill at `--mstep-policy belief` is the cheapest thing that certifies the arm's resume identity.
**Fan stays gated behind the frozen-arm readback + this landing, per the pre-registration.**
