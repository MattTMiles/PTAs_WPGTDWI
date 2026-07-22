# FORGE-G2 — SMASK as a jitted RUNTIME argument (TURBINE's proposal, built and gated)

**Agent:** FORGE-G2 on cronus (jug-gpu, RTX 4090, jax 0.10.1, x64, autotune off).
**Date:** 2026-07-22. Matt commits; nothing here is committed.
**Status:** DONE — **ALL GATES PASS**, every bit-exactness diff **0.000e+00**. GLACIER
Stage-1 unblocked: a frontier-update pattern flip is now an argument change (SG3: warm
flip 6.35 ms median vs the incumbent's 14 s for `_ppab` alone / 126 s through B1Marg).

## 0. State found (inventory before building)

* `reports/FORGE_G2_runtime_smask.md`: **absent** (tracked/staged/untracked — none).
* `reports/TURBINE_speed.md`: **absent from the repo entirely** — no file, no git history,
  no `*turbine*` anywhere. The tasking said the report is the source (TURBINE_results/
  lives on ACCRE); it never landed here. Built from the task brief's numbers (117
  evaluators recompiled per frontier update, 83–800 s/update, 17–32 GPU-hr over GLACIER;
  prototype 3.9 ms flips at max|diff| = 0.0) + `reports/FORGE_G_build.md` (the FORGE-G
  SMASK design + Bg gates). If TURBINE_speed.md lands later, diff its §1 against §2 below.
* Code: `trackB_b1_core.py` / `trackB_b1.py` held the **BAKED-constant** SMASK path —
  `B1Split.set_smask` re-jitted `_ppab` and cleared `_ab_fns`, `B1Marg.set_smask` cleared
  `_abb` and re-jitted `_ppab_b`; smask reached the evaluators by closure through
  `_params` at trace time. No runtime-arg implementation existed. Verdict: **ABSENT → built.**
* `git status`: no staged files from a prior FORGE-G2 attempt (working tree carries
  unrelated referendum/SAMPLER artifacts, untouched).

Recompile tax reproduced on this box before editing (pre-edit capture,
`logs/forge_g2_baked_capture.log`): a baked pattern flip cost **5–15 s for `_ppab`
alone** (1 of the 117), and **126 s for a flip + lnmarg(4)** through the 116 batched
per-pulsar evaluators. Warm evaluation after the recompile: ~10 ms — the work was
never the physics, it was XLA recompiling a constant.

## 1. What changed

### `CW_transition/trackB_b1_core.py`
* **`SMASK_SAFE`** (new module constant): the known-evaluable waveform point carried
  sources are pinned to inside the trace (benign prior-box corner: lowest chirp mass ×
  lowest frequency, zero angles).
* **`MaskedDelay.__call__`**: DOUBLE-WHERE guard (TURBINE's flagged hazard). Carried
  (w_s = 0) sources' waveform **inputs** are `where(gate, actual, SMASK_SAFE)`-sanitised
  *before* the waveform is evaluated; the outer `where` still selects exact 0 for their
  contribution. Fed sources' inputs pass through `where` unchanged → primal bit-exact.
  Without this, autodiff pushes a zero cotangent through a waveform evaluated at the
  carried source's ACTUAL params; junk/coalescence-corner params (SPARK-3 §2) give
  `0 × NaN = NaN` and one carried source poisons the whole gradient. Nothing in the
  branch inspects `w`'s values, so `w` may be a traced runtime argument.
* **`B1Split`**: `smask` threaded as an explicit runtime argument through `_params`,
  `_per_pulsar_ab_impl`, `_pulsar_ab_fn`, and every call site (`a_const_from`, `lnL`,
  `estep`). **`set_smask` no longer clears or re-jits anything** — it canonicalises the
  array, forces fast-full off (MultiSourceDelay does not honour SMASK, unchanged), and
  stores it. `smask=None` builds params WITHOUT the SMASK key → the trace is
  bit-identical to the incumbent/pre-FORGE-G path (default-inert). Only None ↔ array
  retraces (once per direction — the params pytree gains/loses a key); pattern changes
  are plain argument changes.
* Callers passing 3 args (`sampler_core.py`) are untouched: `smask=None` default.

### `CW_transition/trackB_b1.py`
* `make_pulsar_ab_batch` and `B1Marg._ppab_b` take smask as a runtime argument
  (`in_axes=(0, None, None, None)`); `estep_batch` passes `self.sp.smask` through.
* **`B1Marg.set_smask` clears nothing on a pattern change.** The batched evaluators are
  invalidated only when the RESIDUAL path swaps (first non-None smask while fast-full is
  on), because they close over `sp.ys` at make time — that is a one-time event per
  campaign, not per frontier update.

### New files
* `CW_transition/trackB_b1_smask_gate.py` — the FORGE-G2 gate battery (SG1–SG6 below).
* `CW_transition/b1_smask_baked_ref.npz` — pre-edit capture of the TRUE baked incumbent
  (including its single-where masked delay): per-pattern `(a_dep, b)`, split lnL at
  truth/perturbed/junk thetas, full-stack B1Marg lnmarg/posterior at fast-full-default,
  masked-None, and loud-only. This is the only bit-level record of the retired path;
  the gate script consumes it.
* `CW_transition/run_forge_g2_gates.sh` — sequential gate driver (shared GPU).

### Compile budget after this change
Per `B1Marg` lifetime: one trace per evaluator per params-STRUCTURE (smask absent /
smask present), i.e. at most 2× the old single-pattern compile count, paid once.
Frontier updates (pattern changes) cost **zero** compiles. The GLACIER Stage-1 loop's
per-update cost drops from 83–800 s to argument-passing.

## 2. Gates — ALL PASS

Full logs: `CW_transition/logs/forge_g2_smask_gate.log`,
`CW_transition/logs/forge_g2_softsource_gate.log`,
`CW_transition/logs/forge_g2_baked_capture.log` (pre-edit incumbent capture).
Diffstat: `trackB_b1_core.py` +54/−30, `trackB_b1.py` +14/−15.

### Task gate 2 — TURBINE bit-exactness vs the baked path (SG1, SG2)
Runtime-arg path vs the TRUE pre-edit baked incumbent (single-where and all), on
identical inputs — 6 patterns (none / ones / loud-only / zeros / random-0-1 /
fractional weights) × 3 thetas (truth / perturbed / junk-reservoir), at the
`(a_dep, b)` split level and split `lnL`:

| check | max\|diff\| |
|---|---|
| SG1 every pattern × theta, `a_dep`/`b`/`lnL` vs pre-edit capture | **0.000e+00** |
| SG2 freshly-baked constant-capture jit == runtime-arg call (loud/frac/rand01) | **0.000e+00** |

SG1 doubles as the double-where primal-inertness certificate: the capture is the
single-where incumbent.

### Task gate 3 — warm flip < 10 ms, no recompile (SG3)
* Per-pulsar evaluator, median of 20 alternating pattern flips (set_smask + eval,
  block_until_ready): **6.35 ms** vs 6.48 ms fixed-pattern — flip overhead ≈ 0
  (TURBINE prototype: 3.9 ms). Gate < 10 ms: **PASS**.
* Full 116-pulsar `_ppab` under alternating patterns: 7 ms median warm.
* **Jit cache sizes constant across all flips** (per-pulsar 1→1, `_ppab` 2→2;
  SG6 batched evaluators 3→3) — no recompile, asserted, not inferred from timing.
* Incumbent contrast (same box, pre-edit): 14–15 s per flip for `_ppab` alone,
  122–126 s for the first lnmarg after a flip.

### Task gate 4 — default-inert (SG4, SG6)
smask **None** vs the incumbent, bit-identical at **0.000e+00**: split `(a_dep, b)` +
`lnL`, AND full-stack `B1Marg` on both incumbent paths (fast-full default and the
masked path). `sampler_core.py`'s 3-arg calls hit the same None trace (default arg).

### Task gate 1 — FORGE-G Bg gates re-run on the new path
`trackB_b1_softsource.py gate`, unmodified, on the runtime-smask path:

| gate | result |
|---|---|
| Bg2 model-quality-gate limits (0/inf/default, monotone, delta) | PASS |
| Bg4 EM pin drops from the joint free dimension | PASS |
| **Bg5a** always-fed + truth-pinned == incumbent spec-§3 loop | **0.000e+00** on lnL_marg, LNL, q_max; fast-full vs masked 0.000e+00 |
| **Bg5b** zero-noise fixed point | src 0.000e+00, lnL_marg 0.000e+00, census ceiling \|d\| 0.0004 (<0.02) |
| **Bg5c** gate limits + w_s==0 bit-exact absence | all-zero/loud-only patterns exact; junk-reservoir invariance **0.000e+00** |

### The TURBINE grad guard — gradient at a carried→fed boundary (SG5)
Fixture: loud-only smask; the first CARRIED source's params set NON-EVALUABLE
(cos_gwtheta = 1.5 → arccos domain, log10_mc = 30, log10_fgw = −2 → coalescence
inside the span). Scalar = Σ delay² through the module's `MaskedDelay`, one pulsar:

* primal: double-where == single-where control, equal (1.254012e-10 both) — guard inert on values;
* `grad(theta)`: **all finite**; carried 8-block grad **exactly 0**; fed-block grads
  **bit-invariant** to the carried junk (max|diff| = 0.0 across two junk sets);
* `grad(smask)` at the w_s = 0 boundary: finite, exactly 0 on carried entries;
* **negative control**: the pre-guard single-where replica's gradient has **9
  non-finite entries** on the same fixture — the guard is load-bearing, the fixture bites.

### Residual caveats
* One-time cost that remains: the first non-None smask after a residual-path swap
  retraces each evaluator once for the array-carrying params structure (SG6 measured
  122 s for all 116 batched evaluators, paid ONCE per B1Marg lifetime). Every
  subsequent pattern change is free. None ↔ array switching is also cached both ways.
* Fractional weights (w_s ∈ (0,1)) trace the same graph; gated bit-exact vs baked
  (SG1/SG2 `frac`), but as before FORGE-G, only {0,1} patterns carry certified meaning.
* `NoiseDrawer` L_gwb bank warning in the logs is pre-existing (no canonical bank on
  cronus) and irrelevant here — every gate is zero-noise/Asimov.

## 3. Add-list

```
CW_transition/trackB_b1_core.py
CW_transition/trackB_b1.py
CW_transition/trackB_b1_smask_gate.py
CW_transition/run_forge_g2_gates.sh
CW_transition/b1_smask_baked_ref.npz
reports/FORGE_G2_runtime_smask.md
```
