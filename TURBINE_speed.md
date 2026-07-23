# TURBINE_speed.md — the GLACIER loop iteration, priced (agent TURBINE, ACCRE)

*Profiling/optimisation scout for GLACIER. READ-ONLY on all GLACIER/FORGE state: no tracked or
staged file modified; all measurements ran from `TURBINE_results/` (own scratch, own isolated
jax cache `~/.cache/jax_turbine`, starting empty) on the general `batch_gpu` H200 lane —
GLACIER's declared floor lane, matched at `cpus-per-task=8`, autotune off, so timings transfer.
The reserved dgx share untouched. Zero-noise data throughout; nothing science-bearing banked.
GPU spend: T1a job 12641063 (hgx03/H200, 29 m 27 s) + T1b job 12641096 (13 m 56 s) =
**0.72 GPU-hr of the ~4 budgeted**.*

**SCOPE CAVEAT THAT TRAVELS WITH EVERY ABSOLUTE NUMBER.** T1a profiled the *base* b1 problem
(`build_b1_problem`: T≈15-class feathers, rn 30, gwb 14 comps). The campaign venue (T = 40:
101 619 TOAs, rn 64, gwb 30) is ~3.9× heavier per E-step (fan logs: 3.58 s/per-target E-step
vs 0.92 s here). **Ratios and scalings below transfer; absolutes at the campaign venue carry
that ×~3.9.** ncw-scaling used the same population convention (`POP` seed 3000, 3 loud) at
ncw = 16/64/192.

---

## 0. THE VERDICT

**The budgeted 0.43 GPU-hr/rung survives GLACIER's component count on the E-step side — the
per-target E-step is ncw-INVARIANT (0.92 s at ncw=16 → 0.93 s at ncw=192, measured). The two
real threats to the 131 GPU-hr are (1) the SMASK recompile tax, ~7–13 min of pure XLA
recompilation per frontier update per cell at GLACIER scale — removable by a small, measured
bit-exact plumbing change — and (2) the scoreboard floor policy, which is 95.7 % of every
rung and is a pre-registration decision, not a code change.**

| rank | item | class | measured | projected effect on ~131 GPU-hr | size |
|---|---|---|---|---|---|
| 1 | **runtime-SMASK** (smask as jitted arg, not baked constant) | code, **LOW** (gate measured **0.0**) | tax 83→358→~440–800 s/update at ncw 16/64/192; runtime variant: flip = 3.9 ms, bit-exact, no warm penalty | **−(0.12–0.32) × N_it/10 of budget** (N_it = loop iterations/cell); also stops per-pattern cache bloat | small |
| 2 | **floor policy** (cheap path f) | policy, zero code (touches a HELD pre-reg line — Matt's call) | nulls = 400/418 = 95.7 % of rung (fan logs) | budget × (0.043 + 0.957 f); f=0.25 → **−94 GPU-hr** | none |
| 3 | **GPU packing** (2 cells/H200) | launch config, **NONE** | GPU util **mean 0.4 %, median 0 %**; aggregate throughput **1.33× at 2-way** (8-cpu convention kept), 1.53× at 4-way (6 cpus, off-convention) | **−25 % queue occupancy (~33 GPU-hr class) at 2-way** | none |
| 4 | SMASK-aware fast-full residual | code, LOW (existing gate pattern) | masked/fast-full sweep ratio 1.24× at ncw=16; sweep = 32 % of per-target at ncw=192 | ≲ −8 % | small |
| 5 | fringe-grid top-k truncation | code, **HIGH** (tolerance gate) | B 512→128 saves only 40 % of sweep at ncw=16 (launch-latency floor); sweep ≤ 32 % of E-step | ≲ −10 %, not worth the risk now | medium |
| — | mixed precision | **REJECTED** | f32 ulp = 0.977× the signal; pipeline is latency-bound (0.4 % GPU util) | — | — |

If 1–3 land: budget multiplier ≈ (0.043 + 0.957 f) with the recompile term deleted and
1.33× more cells per GPU-hour of queue occupancy. **If only "launch as-is" is on the table:
the budget still holds at ncw≤~200 provided the driver pays the smask tax once per iteration
(not per E-step) and keeps one cell per process — but item 1 is cheap enough that launching
without it burns a measurable double-digit GPU-hr for no reason.**

---

## 1. WHERE 0.43 GPU-hr/RUNG GOES (mined from the 48 H200 fan logs; zero GPU)

```
rung wall : n=43  min 1306  med 1347  mean 1498  max 1927 s   -> mean 0.416 GPU-hr
build     : n=48  med 167 s (once per task in the rung-per-task fan)
```

A rung = **418 per-target E-steps** (`spark3.mode_unit`): 18 signal + **400 null** (100 × 4
states). Per-target E-step: 1498/418 = **3.58 s on H200 ≡ the A100's solo 3.62 s** — the
H200's fp64 buys nothing; the E-step is kernel-launch/latency-bound (116 sequential
per-pulsar jitted calls + host-side numpy reduces), not FLOP-bound.

| component | per rung | share |
|---|---|---|
| **null refits (the floors)** | ~1434 s | **95.7 %** |
| signal-state E-steps | ~64 s | 4.3 % |
| build (amortised) | ~171 s/task | — |
| M-step, SMC ladder (h4d logs, A100) | **159 s** (64 particles, 11 stages, 63 sweeps) | ~10 % of a rung-equivalent |
| M-step setup (z_needle Hessian + brackets + quadrature) | ~15 min, once per unit | amortises |

T1a cost model (base venue, H200): `A` (`_ppab` base eval) = **4.9 ms**, `S` (116 sweeps) =
0.042 s, standard E-step 0.17 s, **per-target = 116·A + S = 0.92 s measured**, `lnmarg`
19.5 ms/particle (the M-step's primitive; a 32-chunk = 0.62 s). Cold compiles on this lane:
one sweep graph ~1.4 s, the batched `lnmarg` set 267 s, full cold E-step 165 s.

## 2. #1 — RECOMPILATION (candidate a): the frontier is a compile bomb; the fix is measured

**Mechanism (source: `trackB_b1_core.py:372-390`).** `PMASK` is a runtime argument — rung
masks recompile nothing (SPARK-3's FIX 1 relies on it). **`SMASK` is baked at trace time**:
`set_smask` clears all 116 per-pulsar sweep evaluators + `_ppab` and re-jits; the smask
constant is embedded in the HLO, so every new feed-weight pattern is a new compile *and a new
persistent-cache key*. GLACIER's whole mechanism is frontier updates (`carried → fed`, L.2)
— i.e. a new smask per loop iteration per cell.

**Measured tax per frontier update** (T1a PS/PN, hgx03, 8 cpus, empty→warm own cache):

| ncw | full tax (extrapolated 116 graphs + `_ppab`) | notes |
|---|---|---|
| 16 | **83 s** (measured full: 94.6 s incl. first E-step) | `_ppab` compile alone 22.7 s |
| 64 | **358 s** | `_ppab` 27.3 s; 1.29 s/sweep graph |
| 192 | **~440–800 s** | `_ppab` 43.0 s; 3.4 s/graph median (one 15.8 s outlier drives the 797 s mean) |

The persistent cache does NOT absorb it: flipping *back* to an already-seen pattern still
costs **63 s** at ncw=16 (only `_ppab` cache-hits, 22.7→3.1 s; the 116 sweep re-traces don't).
Campaign projection at ncw≈192: **144 cells × N_it iterations × ~0.12–0.22 GPU-hr/update ÷ 10
≈ 17–32 GPU-hr at N_it = 10** — 13–25 % of the whole budget, pure compile, and it grows
linearly in N_it and roughly linearly in ncw. It also multiplies cache volume (§7).

**The fix — smask as an argument of the jitted evaluators (the banked padded-N/boolean-mask
pattern; census frozen, shapes never change, mask-only):** `MaskedDelay` already handles a
*traced* smask as written (`params.get(SMASK)` → `jnp.where` SELECT — FORGE-G's NaN-safety
select is trace-compatible). Only the `B1Split`/`B1Marg` plumbing bakes it. T1a PR ran the
prototype on the real machinery:

```
flip w1 -> w2       : 3.9 ms          (baked path: 82.6 s)   -- zero recompile
values vs baked     : max|a,b diff| = 0.0 at BOTH patterns   -- bit-exact
warm cost           : 3.8 ms/probe vs baked 4.4 ms           -- no penalty
compile-once        : same cost as ONE baked pattern
```

**Risk class LOW, with one honest guard.** Forward paths (E-step, scoreboard, SMC `lnmarg`)
are select-protected and bit-exact. The one place runtime differs from baked in *kind*: a
**gradient** through `jnp.where` can propagate NaN from the untaken branch (baked constants
let XLA dead-code the carried source; a runtime gate does not). Any grad/HVP path that
touches `MaskedDelay` with a carried source in the non-evaluable region (SPARK-3 §2.1) needs
either the double-where guard in `MaskedDelay` or a grad-finiteness gate at a
carried-source-in-coalescence fixture. Verification gates for the change: (i) evaluator
(a, b) runtime-arg vs baked at fixed smask `== 0.0` (measured); (ii) Bg5a/Bg5c limits
unchanged; (iii) smask=None default path byte-identical (all pre-FORGE-G callers); (iv) the
grad-finiteness fixture. **Proposal only — the edit lands in `trackB_b1_core.py`/`trackB_b1.py`
(tracked, gated), which is FORGE/Matt's to apply.**

## 3. #2 — THE SCOREBOARD FLOORS (candidate d): 95.7 %, priced as a number

Per-iteration floors are **HELD** in the pre-registration (capstone L.1) — the cheap path
(full refits at key cells + carried-with-caveat elsewhere) is a scoped deviation to *decide*,
not an optimisation to slip in. The arithmetic at the fan's measured anatomy:

- Budget multiplier at full-floor cell fraction f: **0.043 + 0.957·f**
  (f = 96/144 → 0.68×; f = 0.25 → 0.28×, i.e. **~94 GPU-hr of 131 saved**).
- Orthogonal: n_null 100→50 halves floor cost (0.52×) at √2 error inflation
  (±3.5 → ±5.0 nat against floors of 118–744; zf = 0.00 in all 43 fan rungs, so the ANCHOR
  D2 permissive-bias regime is not in play at these states).
- Nulls share the iteration's compiled graphs — with §2 unfixed they at least add no *extra*
  recompiles; with §2 fixed they ride the compile-once economy.

**Each 0.25 of f forgone ≈ 31 GPU-hr.** The decision is Matt's.

## 4. #3 — PACKING (candidate c): the GPU is idle; sell its time back

Occupancy, measured (1 Hz, whole T1a job on a clean H200): **mean GPU util 0.4 %, median
0 %, 99 % of samples < 10 %, mean power 117 W of 700**. Peak VRAM at ncw=192: **1.4 GB of
141**. The device is a latency machine here, not a throughput machine — consistent with the
fan's H200 ≡ A100 equality and HPC_SETUP §7.3's contention-immunity.

T1b, measured — 1 vs 2 vs 4 concurrent warm-E-step processes on ONE clean H200 (100 iters
each, file-barrier-synchronised, warm shared cache):

| packing | per-process E-step medians | aggregate | vs solo |
|---|---|---|---|
| 1 × (8 cpus) | 0.171 s | 5.86 /s | 1.00× |
| **2 × (8 cpus each)** | 0.253 / 0.261 s | 7.78 /s | **1.33×** |
| 4 × (6 cpus each) | 0.441–0.449 s | 8.98 /s | 1.53× |

Sub-linear despite the idle device: the serialisation is host-side dispatch + CUDA-context
time-slicing, not SM contention. **The convention-preserving option is 2-way** (16 cpus/task,
8 per drawing process — fits the ~28-CPU GPU-share): 131 GPU-hr of work in ~98 GPU-hr of
allocation, task walltime ×~1.5 (each process runs 1.48× slower; size walltime accordingly).
The 4-way number carries a confound (6 cpus/process, off the drawing convention) and needs
32 cpus/GPU which the node REJECTS — not available at convention. 3-way at 24 cpus is
untested (interpolates to ~1.45×). Caveat: measured at the base venue (0.17 s E-steps); the
campaign venue's 3.6 s E-steps are still latency-bound so the ratio should hold or improve,
but a 10-minute 2-way spot-check can ride the campaign's warm-cache job before the fan
commits to it.

Constraints that travel with any packing:
- **cpus-per-task=8 per drawing process is load-bearing** — GLACIER has no banked `L_gwb` at
  its venue (header: RECOMPUTED-UNSAFE), so `eigh`'s BLAS thread count is part of the noise
  seed. hgx03 grants ~28 CPUs per GPU-share (a 32-CPU/1-GPU request is REJECTED — measured),
  so **≤3 packed processes at the pinned convention**. The clean alternative: bank `L_gwb`
  once at the venue with the existing tracked tool (`CW_transition/make_lgwb_bank.py`) on the
  pinned host/thread count — zero new code, provenance posture change, GLACIER's call.
- Slurm fences GPUs correctly on this lane; the standing co-tenant on hgx03 (a `sys/dashboard`
  job, 1 GPU / 28 CPUs, 8-day walltime) plus the known squat-lottery class are timing-only
  hazards; the fan's own 8-way node sharing showed no dgx04-style degradation (224 cores).

## 5. #4/#5 — SMALLER, AND ONE REJECTION

- **SMASK-aware fast-full** (`MultiSourceDelay` honouring SMASK, so a live smask keeps the
  pterm-only residual): masked/fast-full sweep ratio measured 1.24× (ncw=16); the sweep is
  ≤32 % of the per-target E-step at ncw=192 → **≤8 % saving**. LOW risk (the existing
  `gate_fast_full` bit-exact pattern extends verbatim). Worth bundling with §2's edit, not
  worth a launch delay alone.
- **Fringe-grid top-k truncation** (candidate b): B 512→128 saves only 40 % of the sweep at
  ncw=16 — the sweep is launch-latency-floored, and the sweep itself is a minority of the
  E-step. ≲10 % end-to-end for a HIGH-risk numerics change (within-fringe peaks move; needs a
  q_max/dlnL_det tolerance gate vs full grid). **Defer unless the budget still hurts after
  items 1–3.** B_CERT = 512 as a *criterion* constant is untouched either way.
- **Mixed precision (candidate e): REJECTED.** (i) f32 ulp at the pulsar-term phase
  (2πL/dL ≈ 1.6e4 rad) is 9.766e-4 rad = 0.977× the ~1e-3-nat fringe-breaking signal
  (HPC_SETUP H4a) — f32 in waveform assembly *is* the signal's size; (ii) the pipeline runs
  the GPU at 0.4 % — halving arithmetic width attacks a bound this loop never touches. The
  f64 Cholesky mandate (Schur) stands. No gateable f32 section is worth its risk.

## 6. NCW SCALING — the budget's base survives (T1a PN)

| ncw | build | A (`_ppab`) | sweep probe (12 psr) | per-target E-step | smask tax | peak VRAM |
|---|---|---|---|---|---|---|
| 16 | 183 s | 4.9 ms | 4.4 ms | **0.92 s** (measured) | 83 s | 0.6 GB |
| 64 | 229 s | 5.2 ms | 11.7 ms | 0.71 s (est) | 358 s | 1.0 GB |
| 192 | 441 s | 5.5 ms | 30.7 ms | **0.93 s** (est) | ~440–800 s | 1.4 GB |

The per-target E-step is **flat in ncw**: it is 116·A-dominated and A does not see the source
count (the sweep grows sub-linearly, 7× for 12× sources, and stays a minority). At the T=40
venue (×~3.9) the ncw≈200-class rung therefore stays ≈ the budgeted 0.43 GPU-hr — **the
131 GPU-hr base is sound at GLACIER scale on the E-step side.** What scales is the *compile*
(§2) and the build (441 s/process at ncw=192 — one more reason for one-cell-per-process,
never one-iteration-per-task; per-process graph materialisation is ~2 min even warm,
HPC_SETUP §7.4).

## 7. CACHE + LAUNCH HYGIENE (candidate f)

- **Shared warm campaign cache, not per-task**: SPARK-3's fan ran 48 tasks on one shared
  cache, zero races (the EMBER race was autotune textprotos; autotune is off).
  `gl_g0.sbatch`'s per-run cache is right for a gate, wrong for 144 cells.
- **NEW HAZARD, measured: H200 cache entries are huge.** `sp3_h200_shared` = **52 GB / 41
  entries** on the 186 GB `/home` quota (72 GB already in use pre-SPARK). Without §2's fix,
  every smask pattern mints new entries — a GLACIER campaign can plausibly hit the quota
  mid-fan. Prune stale lanes before launch; put a `du` + quota check in the fan prologue;
  with §2 fixed the growth collapses to one graph set.
- **Compile-once-loop-many** stands (SPARK-2): one cell per task, all iterations in-process.
- Walltime: measured × 1.3, rounded up to the 30-min backfill bucket (`bf_resolution=1800`).
- The lane starts jobs in seconds when GPUs are free (submit→start 27 s for the fan; 14 s for
  T1a) — `--test-only` estimates on this lane are noise (it quoted 2026-10-11 for a 30-min
  job while 7/8 H200s idled). Do not size the launch strategy on them.

## 8. T3 — NO-RISK, ZERO-CODE, READY

`TURBINE_results/gl_fan_TEMPLATE.sbatch` — a campaign-fan template embodying: shared warm
cache + quota check, cpus-per-task=8 kept load-bearing, 30-min-bucket walltime, `%6` throttle
(8-GPU lane minus the standing dashboard co-tenant), `--requeue`, foreign-resident print,
one-cell-per-task-all-iterations-in-process. Everything touching computation (§2, §5) stays a
proposal; §3 is a pre-registration decision.

## 9. STATE ON DISK / PROVENANCE

| thing | where |
|---|---|
| this report (staged) | `TURBINE_speed.md` |
| T1a profiler + raw results | `TURBINE_results/turbine_t1a.py`, `turbine_t1a.json`, `logs/t1a_12641063.*`, `logs/t1a_util_12641063.csv` |
| T1b packing worker + raw results | `TURBINE_results/turbine_worker.py`, `turbine_t1b_*.json`, `logs/t1b_12641096.*` |
| T3 template | `TURBINE_results/gl_fan_TEMPLATE.sbatch` |
| GPU spend | T1a 0.49 + T1b 0.23 — total **0.72 GPU-hr of ~4 budgeted** (sacct-verified) |

Nothing committed. No tracked or staged file touched. GLACIER's session, results dirs,
caches and the reserved dgx share untouched. Add-list: **`TURBINE_speed.md` only.**
