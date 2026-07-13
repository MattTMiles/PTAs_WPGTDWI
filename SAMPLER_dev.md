# SAMPLER — prong 3, actually run: the honest posterior machinery (DEV PHASE)

**Agent:** SAMPLER · cronus / RTX 4090 · **PURE EXECUTION** (no commits, no tracked-file edits;
canonical docs untouched per the doc-fence). **Date:** 2026-07-12.

**Code (untracked, `CW_transition/`):** `sampler_core.py` (the differentiable target),
`sampler_gates.py` (G0/G1), `sampler_unit.py` (G1a/G1b), `sampler_graddiag.py` (the G1 anatomy),
`sampler_run.py` (NUTS driver + g1/g2), `sampler_s4.py` (g4 / S-4), `sampler_sbc.py` (SBC harness).
**Banks:** `SAMPLER_results/`. **Logs:** `CW_transition/logs/sampler_*.log`.

---

## 0. THE HEADLINE

**The mission's premise was wrong in a way that mattered, and fixing it is the session's main
deliverable.** `trackB_b1.B1Marg.lnmarg` — the object the brief names as "your likelihood
object", the one deliverable R measured and IGNITE-2's soft loop ran on — **has no JAX
gradient**. Its fringe reduction (`marg_from_LNL`) is numpy (`np.maximum.reduceat` inside a
Python loop over 116 pulsars), and IGNITE-2's soft M-step accordingly used **finite-difference**
gradients (`reports/IGNITE2_softloop.md` §3: *"damped Newton (damping 0.5) with FD gradient"*).
So the pre-registered gate — *"verify grad(lnL_marg) against FD at 8 points, 1e-8, before
anything"* — **could not be run against the existing object at all.**

The target was therefore **rebuilt as a differentiable JAX twin** (`sampler_core.MargJax`),
gated **value-identical to `B1Marg`** and **gradient-exact against a monolithic `jax.grad`**.
Along the way the build measured three hard facts that govern every posterior this programme
will ever run, and that are the real content of the dev phase:

1. **The naive autodiff of `lnL_marg` does not fit on a 24 GB card.** A single reverse-mode
   graph over 116 pulsars × 512 fringes requests a **16.5 GiB single allocation** and OOMs.
   The target had to be restructured (fringe-chunked `lax.scan` + per-pulsar VJP assembly +
   an XLA-opaque `custom_vjp`, so a sampler may jit *around* `lnL_marg` but never *through* it).
2. **The 1e-8 FD bar is unreachable on this target, and the reason is the physics, not the code.**
   `lnL_marg` is assembled as `Σ_p M_p − 115·lnL_ref`, a difference of ~4.7e7-sized quantities
   yielding ~4e5 — so the *value* carries ~1e-7 nat of cancellation noise. Combined with the
   ~1e10 conditioning, **no single FD step certifies all 18 dims**: the step that resolves the
   registration dims (`log10_fgw`, `log10_mc`) is roundoff-swamped in the flat dims, and vice
   versa. **The gate was replaced with a strictly stronger one** (AD vs AD, below).
3. **Per-gradient cost on the 4090 is 1.09 s** (116 pulsars, B = 512, 18 active dims). That
   number, not the algorithm, is what sends SBC to ACCRE.

---

## 1. THE TARGET DENSITY (what is sampled, what is marginalised)

    p(theta_act | d)  ∝  exp[ lnL_marg(theta) ] · prior(theta_act)

    lnL_marg(theta) = lnL_ref(theta) + Σ_p LSE_n[ (LNL_p(n) − lnL_ref) + lnprior_p(n) ]
                    ≡ Σ_p M_p(theta) − (npsr−1)·lnL_ref(theta)      [the regrouping used here]

| quantity | treatment |
|---|---|
| fringe index `n_p` (116 pulsars) | **MARGINALISED** (the LSE). Never sampled as a discrete parameter. |
| within-fringe offset | **PINNED AT THE FRINGE CENTRE (u = 0) for 115/116 pulsars** — neither profiled nor marginalised. `A2.eval_grid` returns *fringe centres* whenever the window holds > `DENSE_FRINGE_MAX = 64` fringes, so `segment_max` is the **identity**. **See §9 — this is the session's sharpest finding, and it retracts an earlier claim in this report.** |
| GWB + per-pulsar RN GPs | **MARGINALISED ANALYTICALLY AND EXACTLY** (GPs in the covariance; hyperparameters frozen). |
| sky | **FIXED** at truth (the EM-counterpart convention). Freed only in g4. |
| faint 13 sources | fixed at truth (R's ACTIVE-DIM convention). |
| **ACTIVE (sampled)** | 3 loud × {`cos_inc`, `log10_mc`, `log10_fgw`, `log10_h`, `phase0`, `psi`} = **18 dims** |

Priors: uniform on the spec's physical box (`trackB_estimator._PHI_LO/_PHI_HI`), sampled through
a logit transform so no proposal leaves the box (the waveform nan-clips outside it); the
log-Jacobian is carried exactly.

**Distances are recovered post-hoc**, never sampled: at posterior draws the per-pulsar fringe
posterior `q_p(n | theta_s)` is read off and averaged over draws (`sampler_run.fringe_posterior_at_draws`).

---

## 2. GATES — machinery

### G0 — VALUE: `MargJax.lnmarg` == `B1Marg.lnmarg`. **PASS**

8 points (truth + a geometric perturbation ladder to 1e-1 scaled), 116 pulsars, B = 512:

    max |MargJax − B1Marg| = 5.146e-08 nat        on |lnL_marg| ~ 4.06e5   (rel ~1e-13)

and the identity cross-checks the banked deliverable R exactly:

    lnL_marg(truth) = 405686.343447        R banked: 405686.34   ✓

The residual 5e-8 is the float64 cancellation floor of the assembly itself, not a discrepancy.
The `B1Marg` fast-full path gate also re-passed at exactly **0.0**.

### G1a — GRADIENT: hand-assembled per-pulsar VJP == monolithic `jax.grad`. **PASS (1.1e-15)**

The mission's FD gate cannot certify this target (§0.2), so the chain rule is gated by the
**strictly stronger AD-vs-AD test** on a small array (npsr = 8, ncw = 2), where the monolithic
graph still compiles and the two gradients must agree to machine precision:

    |value difference|            = 0.000e+00
    max relative gradient error   = 1.131e-15      (bar 1e-10)   -> PASS

**Same code path**, only the array size differs. A wrong chain rule (a missed or duplicated
term) would show O(1) error here, not 1e-15.

### G1b — GRADIENT vs FD (the pre-registered form): FD-noise-limited, **reported not enforced**

Same small array, Richardson-extrapolated central differences, best step per dim:

| param | AD | Richardson FD | rel err |
|---|---|---|---|
| src0_cos_inc | 3.42152729e+00 | 3.42152722e+00 | 2.0e-08 |
| src0_log10_mc | −1.92825919e+03 | −1.92825918e+03 | **4.4e-09** |
| src0_log10_fgw | 1.91365013e+03 | 1.91365031e+03 | 9.4e-08 |
| src0_log10_h | −7.71146069e+00 | −7.71146057e+00 | 1.6e-08 |
| src0_phase0 | −1.12136926e-01 | −1.12136929e-01 | 2.5e-08 |
| src0_psi | 2.20840840e-01 | 2.20840873e-01 | 1.5e-07 |

Worst 1.5e-07 against the 1e-8 bar. **This is the FD's error, not the gradient's** — G1a proves
the gradient exact to 1e-15 on the identical code. Recorded as measured.

### G1 anatomy — why the full-array FD gate reads "8 % error", and why that is the 4e10 hazard

On the full 116-pulsar array, `sampler_graddiag.py` steps h over 1e-4 → 1e-7 (scaled):

| h | max abs gap | median abs gap |
|---|---|---|
| 1e-4 | 8.3e+04 | 2.5e-04 |
| 1e-5 | 1.9e+04 | 1.1e-03 |
| 1e-6 | 8.8e+02 | 8.4e-03 |
| 1e-7 | 5.3e+01 | 9.9e-02 |

The two families move in **opposite directions**:

- **Registration dims** (|grad| ~ 5e4): FD **converges onto AD** as h shrinks —
  `src0_log10_mc`: FD = −6.1e3 → 5.05e4 → 6.95e4 → **6.9703e4** against AD **6.9705e4**.
  At h = 1e-4 the FD is *garbage* because the step is far larger than the needle's curvature
  scale (σ ~ 1e-5–1e-6 scaled; R measured σ_sky 2.4–9.7e-6).
- **Flat dims** (|grad| ~ 1e0–1e2): FD matches AD to 1e-5–1e-6 at h = 1e-4 and **degrades** to
  garbage by h = 1e-7 (the ~1e-7 nat value noise divided by 2h).

**There is no h that certifies both.** The first G1 run's "8 % error on `src1_psi`" was exactly
this: a dim whose gradient is 1e-1 against a max of 7e4, read at a step tuned for neither.
*The gate did not fail; it measured the conditioning.*

**Fringe-group occupancy** (which decides whether the `segment_max` can even produce kinks):
K = 512 for **111 of 116** pulsars — one grid point per fringe, so the reduction is the identity
and the target is smooth there. Only **5 tight-prior pulsars** have K < 512 (min 25; up to 374
grid points inside one fringe) and actually engage the max.

### G2 — COST on the 4090 (116 psr, B = 512, 18 active dims)

| quantity | value |
|---|---|
| `value_and_grad` (F_CHUNK = 64) | **1.09 s** |
| `value_and_grad` (F_CHUNK = 8) | 6.98 s |
| naive single-graph `jax.grad` | **OOM** (16.5 GiB single allocation) |
| GPU peak, MargJax process | 0.6 GiB (F_CHUNK = 8) |
| GPU peak, `B1Marg` reference process | 15.4 GiB |
| first-build compile (116 per-pulsar gradient graphs) | 350–780 s (persistent cache warm: less) |

**Memory conventions this forced, all load-bearing:**
1. `F_CHUNK` — fringes per `lax.scan` step, each chunk `jax.checkpoint`-ed, so reverse mode
   recomputes one chunk instead of taping all 512.
2. **Per-pulsar VJP assembly** — no graph ever holds all 116 blocks.
3. **`custom_vjp` + `pure_callback`** — NUTS jits its potential; without this the sampler's own
   jit re-inlines the 116 blocks into the monolith and OOMs again.
4. **The reference and the twin cannot co-reside** (15.4 + build > 24 GB). The G0 gate runs
   `B1Marg` in its own process and banks the values.

---

## 3. GATE LADDER — posteriors

### STATUS AT REPORT TIME: **the posterior gates did not run to completion on cronus, and the
reason is a measured cost, not a bug.**

| gate | status | what exists |
|---|---|---|
| **g1** posterior concentrates at truth, σ(log10 mc) vs IGNITE-2 | **NOT RUN to completion** | driver built + launched (`sampler_run.py g12`); killed by the host-RAM OOM of §3.0, not relaunched (see the cost, below) |
| **g2** two chains, R-hat < 1.01 | **NOT RUN** | same driver; R-hat/ESS machinery in `sampler_run.diagnose` |
| **g3** noisy realisation, fringe posteriors from posterior-averaged sources | **NOT RUN** | machinery built (`fringe_posterior_at_draws`, `cert_from_qbar`) |
| **g4** S-4 revisited | **IN FLIGHT at report time** | `sampler_s4.py`, running; conditioning at init banked (below) |
| **SBC** N = 10 smoke | **NOT RUN** | harness built (`sampler_sbc.py`), incl. the fringe PIT |

**The arithmetic that stopped it, and it is the session's second real finding.** One gradient of
the fringe-marginalised target costs **1.09 s** (§2, G2). A NUTS chain of 300 warmup + 500 samples
at tree depth 6 needs ~2.4e4 gradients ≈ **7 GPU-hours** — *per chain*. g1/g2 needs two; SBC at
N = 10 needs ten; the mission's N ≥ 200 production run needs ~1400 GPU-hours. **The dev-phase
gates g1–g3 are themselves an overnight-to-multi-day job on one 4090**, which is why the brief's
own split (dev on cronus, ensembles on ACCRE) is right but drawn one notch too optimistically:
**even the single-realisation gates belong on ACCRE.** The handoff (§5) is therefore the primary
deliverable of this session, and it is sized from measured numbers rather than guesses.

### g4 — S-4 REVISITED

**PRE-REGISTERED READING (fixed before the chains finish).** The deliverable is the coverage
DIRECTION, not a precision number; at 2 + 2 realisations the binomial error is large and is stated
rather than hidden. The verdict is one of:
- **RESTORED** — scenario-C coverage consistent with nominal (0.90) within binomial error at this
  resolution ⇒ RING's S-4 undercoverage was the ESTIMATOR (profile slice + Laplace), and the
  standing stop-point closes.
- **STILL-UNDER** — coverage below nominal with the credible level at truth piling up near 1
  ⇒ the undercoverage is PHYSICS (a correlated background shifts which fringe/peak the MAP selects),
  and it survives an honest marginal posterior with the GWB analytically marginalised. Quantified.
- **INDETERMINATE-AT-THIS-DEPTH** — chains not converged / trees truncated at the depth cap ⇒ say
  what it would take and hand to ACCRE with the SBC production.

#### THE SADDLE DATUM — banked now, independent of any chain

Read off the log (`sampler_s4b.log`, seed 770001, scenario C, exact distances, GWB in the data):

    Hessian at init: eig(-H) in [-1.710e+09, 8.916e+08];  cond = 3.20e+08;  neg eigs 3/8

**Three of eight eigenvalues are NEGATIVE at truth: on a noisy realisation the exact-distance
posterior has a SADDLE at the true parameters, not a maximum.** This stands whatever the chains do,
and it is already an S-4-relevant statement: a credible region built by *expanding around* truth —
which is what both a profile slice and a Laplace approximation do — is describing a surface whose
curvature at truth is **indefinite**. The estimator's premise (truth ≈ mode ≈ centre of a locally
Gaussian bowl) fails before any coverage is counted.

**And it is the same physics as R's mc micro-dip / thin-shell, in a second formalism.** R found the
posterior mass sitting on a high-likelihood *shell just outside* the needle rather than at the peak;
here the curvature at truth is indefinite in 3 of 8 directions, so the mode is displaced off truth
by the noise. Two independent formalisms — R's evidence integral and this Hessian — reporting that
**truth is not where the posterior piles up.** The micro-dip was pre-registered in the mission as a
hazard the sampler must traverse; it is now also a *geometric* statement about the target.

#### g4 RESULT — **VERDICT: INDETERMINATE-AT-THIS-DEPTH, with the direction against RESTORED**

4 realisations (2 scenario-C + 2 controls), exact distances, sky free, all 8 loud-source params
sampled, GWB drawn into the data and **marginalised analytically and exactly** in the likelihood.
**Zero divergences in all four.** Banked: `SAMPLER_results/s4_scen{B,C}.npz`.

| scenario | noise | c(truth) — the credible level at truth | inside90 |
|---|---|---|---|
| **C** | white + RN + **GWB** | **0.903**, **0.813** | 1 / 2 |
| **B** (control) | white + RN | **0.053**, 0.753 | 2 / 2 |

Under a calibrated posterior `c ~ U(0,1)`. **Both scenario-C realisations land above 0.81**
(P = 0.2² = 0.04 under calibration, one-sided, N = 2 — suggestive, not decisive); the controls are
unremarkable (0.053, 0.753). `inside90` alone cannot separate the arms at this N: the binomial 95 %
interval on 1/2 spans [0.01, 0.99] and contains 0.90, so **RESTORED cannot be excluded on the count
statistic, and STILL-UNDER cannot be established on it either.** That is why the verdict is
**INDETERMINATE-AT-THIS-DEPTH** rather than either pole — the honest reading of 4 chains.

**What it would take** (the ACCRE ask): the continuous `c` statistic is far more powerful than
`inside90`. A **KS test of {c_i} against U(0,1) needs ~12–15 realisations per arm** to separate
"C is skewed high" from noise at 95 %; matching RING's binary `inside90` resolution (0.5 vs 0.9)
would need ~25–30 per arm. At **~38 min/realisation** on the 4090 (300 draws, depth 5, 0 divergences)
that is ~10 GPU-hours per arm — **one A100 afternoon**, and it rides along with the SBC production
(§5) since it shares the build.

**What is NOT indeterminate — and it is the finding that survives:** truth is a **SADDLE** of the
exact-distance posterior on **every one of the 4 noisy realisations** (negative eigenvalues 3, 4, 3,
5 of 8; `cond(-H)` = 3.2e8 / 2.2e8 / 4.9e7 / 1.3e9), **with and without the GWB**. So the saddle is
noise-induced, not GWB-induced. RING's S-4 asked whether its undercoverage was physics or estimator;
this says **the estimator's premise is broken independently of the answer** — a profile slice and a
Laplace region both expand around a point that is not a mode, on every realisation tested. The
GWB-specific question (does a correlated background *additionally* skew the marginal posterior?) is
what the ~15-realisation run must settle, and the 2 measured C-values lean toward "yes".

### g4 — what else is banked

The exact-distance (RING tier3) target, GWB drawn into the data and **marginalised analytically
and exactly** in the likelihood, sky free, NUTS over all 8 params of the loud source. Conditioning
of the posterior at truth, **read off the log** (`sampler_s4b.log`, seed 770001, scenario C):

    Hessian at init: eig(-H) in [-1.710e+09, 8.916e+08];  cond = 3.20e+08;  neg eigs 3/8

**Three of the eight eigenvalues are NEGATIVE: on a noisy realisation, TRUTH IS A SADDLE of the
exact-distance posterior, not a maximum.** That is measured, and it is already an S-4-relevant
statement independent of the coverage number: a credible region built by *expanding around* a
point that is not the mode — which is what a profile slice and a Laplace approximation both do —
is describing a surface whose curvature at truth is indefinite. RING's S-4 asked whether the
undercoverage was physics or estimator; this says the estimator's premise (truth ≈ mode ≈ centre
of a locally Gaussian bowl) **fails at truth on noisy data**, before any coverage is counted.

**The coverage numbers themselves are NOT reported: the run had not banked `s4_scen*.npz` at
report time.** Per §3.0, no number enters this report until it is read back off its bank.

### 3.0 AN EXECUTION FAULT, RECORDED IN FULL (because it nearly became a fabricated result)

The first g4 and g1/g2 runs were launched **concurrently**. Each builds the full 116-pulsar
enterprise + discovery likelihood, and `sampler_s4` additionally called `jax.hessian` on the
exact-distance target — forward-over-reverse through the whole array likelihood. The box
(62 GB) ran out of **host** RAM and the kernel OOM-killed the g4 process:

    kernel: Out of memory: Killed process 1783720 (python) total-vm:66831216kB anon-rss:34441160kB

It died **silently** — no traceback, no non-zero exit visible in the log, which ends mid-run at
`[C0] seed 770001`. The g1/g2 process died with it.

**The failure mode that matters:** a progress monitor surfaced plausible-looking S-4 numbers
(`c(truth) = 0.985 / 0.999 / 1.000`, "cond 3.96e10", "1500 draws, 158 s, 0 divergences") that
**appear nowhere in the log file and correspond to no banked npz**. They were relayed in-session
before being checked. They are **discarded, not reported**: the log contains zero `c(truth)`
lines, zero `Hessian at init` lines, and `SAMPLER_results/` contains no `s4_scen*.npz`. Every
number in this report is one that survives re-reading its own bank from disk.

*Convention this earns, and it is the same one the programme already has in a different dress:*
**a number is not a result until it is read back off the artifact it claims to come from.**
The banked-npz convention exists precisely so that verdicts cannot outrun their evidence; a
streamed log line is not a bank.

**Fixes applied:** (i) one 116-pulsar process at a time — these jobs are host-RAM-bound, not
just GPU-bound; (ii) `jax.hessian` replaced by a finite-difference Hessian **of the gated
gradient** (2·ndim gradient evaluations, a few hundred MB, instead of a whole-array
forward-over-reverse graph).

---

## 4. WHAT THIS CHANGES FOR PRONG 3

1. **The reference implementation of spec §3 now has a gradient.** IGNITE-2's soft loop — the
   adopted reference M-step — has been running on **finite-difference** gradients of a numpy
   object. An exact analytic gradient of the *same* objective now exists and is gated to 1e-15
   (G1a). `sampler_core.MargJax` is a drop-in for any future soft-loop M-step, and per STORY
   S10.1.6 (*"replace Adam with a Newton/trust-region step"*) that is the direction of travel.
   **This is reusable beyond the sampler.**
2. **The 4e10 conditioning is real, and it now has a second measurement.** It shows up as the
   inability of *any single FD step* to certify the gradient across the 18 dims (§2), and as
   `cond(-H) = 3.2e8` with **3 negative eigenvalues at truth** on the exact-distance noisy
   posterior (§3). The mass matrix must be seeded from the target's own Hessian; window-adapted
   diagonal mass will not find a 1e8–1e11 dynamic range inside any affordable warmup.
3. **The honest expectation stands** (STORY S10.1, deliverable (v)): given `f = 6.9e-7`, a correct
   sampler should report a **diffuse fringe posterior**, and *a correct diffuse answer is a
   result, not a failure.* Nothing measured here contradicts that; nothing here confirms it
   either, because the posterior gates have not run.

---

## 5. ACCRE HANDOFF SPEC (the at-scale SBC ensembles)

Sized from the **measured** cronus costs, not from an estimate.

### 5.1 The unit of work

One SBC realisation = prior draw (18 active dims) → arm-B fringe truth draw (116 distances) →
simulate (CW + white + RN + HD-GWB) → NUTS on `lnL_marg` → rank statistics + the fringe PIT.

| quantity | measured on the 4090 | note |
|---|---|---|
| `value_and_grad(lnL_marg)` | **1.09 s** | F_CHUNK = 64, 116 psr, B = 512, 18 dims |
| build + 116 gradient compiles | 350–780 s | **once per process**; JAX persistent cache makes reruns cheaper — the cache dir must be on a shared filesystem or every task pays it |
| GPU memory (MargJax) | < 1 GiB steady | after the restructure; the naive graph needs 16.5 GiB and OOMs |
| **host** RAM per process | **~34 GB peak** | the OOM-killer datum. **One 116-pulsar process per node**, and request ≥ 48 GB |

**Per realisation:** N_grad ≈ (n_warm + n_samp) × 2^depth_eff. At depth 6 and 300+500 draws with
a typical tree of ~30 leapfrog steps: ~2.4e4 gradients ≈ **7 GPU-hours**. This is the number that
decides the campaign: **SBC at N ≥ 200 is a ~1400 GPU-hour job** and does not fit on cronus.

### 5.2 Task batching (per HPC_SETUP conventions)

- **One realisation per SLURM task**, `--array=0-199`, seeds `9.0xx M + 100·i` (disjoint from
  every banked range). Skip-on-exist resume (`sbc_{i:03d}.npz` present → skip), drilled the way
  IGNITE-2 drilled its nulls array.
- Lane: the A100-80GB `interactive_gpu` lane IGNITE-2 used. The A100's larger memory permits
  **F_CHUNK = 256+**, which should cut the per-gradient wall further — *re-time it there; do not
  inherit the 4090 number.*
- **Bank the statistic, not the verdict** (the standing convention): each npz carries the full
  chain (`X`), per-dim ranks, the fringe PIT per pulsar, `q_bar_p(n)`, `n_true`, divergences,
  tree depths, step size, and the Hessian eigenvalues at init. Every one of those is needed to
  re-cut calibration on a different bar without re-running the GPU job.
- Checkpoint **per chunk** inside a chain, not only per realisation (`sampler_run.run_chain`
  already writes `chain_{tag}.npz` every chunk), so a pre-emption loses minutes, not hours.

### 5.3 Shapes

    sbc_{i:03d}.npz
      X          (n_samp, 18)      posterior draws, physical units
      rank       (18,)             #{s: X[s,j] < x_true[j]}       -> Uniform{0..S} if calibrated
      pit        (116,)            randomised PIT of the TRUE fringe under q_bar_p  -> U(0,1)
      qmax,mapk  (116,)            posterior-averaged fringe posterior summaries
      ptrue      (116,)            q_bar_p(n_true)
      n_true     (116,)            the drawn fringe truth (arm B)
      num_steps  (n_samp,)         leapfrog per draw (the cost model, banked)
      cond, div, eig               conditioning + sampler health

### 5.4 What the production run answers

The fringe PIT is **the calibration test the programme has never run** (STORY S10.1 (ii)): it is
the direct check of whether `q_max` is a probability or merely a score — the assumption under
*every* certification count in this repo.

---

## 6. CAVEATS THAT TRAVEL

- **The within-fringe offset is PROFILED, not marginalised.** `lnL_marg` takes the max over the
  B_CERT grid points inside each fringe. This is inherited from `B1Marg`/R unchanged (so every
  number here is directly comparable to `f = 6.9e-7` and to IGNITE-2), but it means the object is
  a *profile* in the sub-fringe distance and a *marginal* only in the fringe index. A fully
  marginal treatment would integrate within the fringe. **Not fixed here; flagged.**
- **The VLBI *onset* cell (h = −13.25, T = 30 yr, vlbi) cannot be rebuilt on cronus.** IGNITE's
  cell machinery (the T-extension: TOA append past each pulsar's last TOA at its own cadence +
  timing-design extrapolation with the R² = 0.99 gate + RN/GWB component rescaling) lived in
  ACCRE scratch (`hpc_harbor/`) and **is not in this repo**. The loudness axis is a no-op at
  h = −13.25 (the pop draw's own loudness) and the VLBI tier is a prior-column swap
  (`σ_d → min(lit, 1 pc)`), but **the T = 30 baseline is a substantial rebuild**. What is run
  here is therefore the **fiducial-span array**, and the cell is named as such — not relabelled.
  Also note IGNITE-2 **retracted the onset itself**: under D2-sized floors neither pre-registered
  cell clears onset, so "the VLBI onset cell" is a nominal address, not a certified one.
- **The VLBI tier as implemented here (`vlbi_all`) applies `min(lit, 1 pc)` to EVERY pulsar**,
  because IGNITE's GEO union-18 membership list is not banked on cronus. That is *more*
  optimistic than IGNITE's tier. Stated wherever used.
- **g4 runs in the B1 harness, not RING's.** 116 real pulsars (not the 30-pulsar ring), the
  pop-3000 loud source at `log10_fgw = -7.559` (RING's cell was `fgw = -8`), `log10_h = -13.25`,
  `log10_mc = 9.145`. It answers RING's *question* — does an honest marginal posterior with a
  properly-marginalised GWB restore coverage? — in a harness where the estimator leg is removed
  by construction (no profile slice, no Laplace, no grid). It is **not** a re-run of RING's cell,
  and the `inside90` numbers are therefore not numerically comparable to RING's 0.40–0.50; only
  the *verdict* (does coverage return to nominal?) transfers.
- **The SBC init convention.** Chains start at the realisation's own truth. SBC tests the
  POSTERIOR's calibration, not the sampler's global search; a cold init would confound the two,
  and the referendum (`f = 6.9e-7`) already answers the global-search question in the negative.
  **This means SBC as specced cannot detect a multimodality the sampler never visits** — stated,
  because it bounds what a passing SBC would license.

---

## 8. ARCHAEOLOGY ADDENDUM — WHY SNAP-TO-PEAK WORKED, AND WHAT IT WAS ACTUALLY DOING

*Read-only. Sources: `CW_node_sampling/` (MK2–MK9 + `EXPERIMENT_SUMMARY.md`, `README.txt`),
`reports/ignite_bank.npz` (2 070 realisations × 116 pulsars, banked statistics). Ran alongside
the g4 job; no GPU used.*

### 8.1 THE HYPOTHESIS, AND THE VERDICT

> *At loud SNR the fringe posterior `q_p(n)` is near-delta, so argmax == marginal, and per-draw
> re-snapping approximately profiles the discrete fringes out — snapping was collapsed sampling
> avant la lettre, valid exactly where `q_p` is concentrated.*

**VERDICT: the mechanism is CONFIRMED, the inference from it INVERTS.** `q_p` does become a delta
as loudness rises — precisely as claimed. But **the delta sharpens onto the WRONG fringe more and
more often as it sharpens.** Concentration and correctness move in *opposite* directions in h. So
"snapping ≈ marginalising" becomes a *better and better approximation of a worse and worse answer*.
The approximation error vanishes; the answer does not become right.

Measured, from the IGNITE bank (`kind == 'sig'`: the 1 200 pure signal realisations, per-pulsar;
`qmax` = posterior mass on the MAP fringe, `purity` = fraction of confident locks on the TRUE
fringe). Banked to `SAMPLER_results/archaeology_snap.npz`:

| h | median `q_max` | frac `q>0.9` | **purity of `q>0.9` locks** | P(confident **and** true) | median \|Δk\| of confident-WRONG locks |
|---|---|---|---|---|---|
| −13.25 | 0.164 | 0.063 | **0.785** | 0.049 | 9 |
| −13.00 | 0.356 | 0.144 | **0.612** | 0.088 | 57 |
| −12.75 | 0.631 | 0.308 | **0.455** | 0.140 | 85 |
| −12.50 | **0.936** | 0.537 | **0.366** | 0.197 | 98 |

(VLBI tier, T = 30, the tight-prior cell, same story: median `q_max` 0.309 → **0.996** while purity
falls 0.775 → 0.371.)

**Read the first and third columns together.** By h = −12.5 the fringe posterior is essentially a
delta (median `q_max` = 0.936) — argmax *is* the marginal, exactly as the hypothesis says — and
**~2 of every 3 confident locks are on the wrong fringe**, a median **98 fringes** away. Meanwhile
at census loudness, where the posterior is *diffuse* (median `q_max` = 0.164), the few confident
locks that do occur are **78 % correct**. **Confidence and correctness are anti-correlated in
loudness.** The delta-function limit is real; it is a limit of *confidence*, not of *correctness*.
This is the same object IGNITE named ("tight local width + wrong global registration = confident
nonsense", §10.8.3) seen from the sampler's side, and it is why *"snap was collapsed sampling"* is
true and **does not license snapping**.

### 8.2 WHAT THE OLD SAMPLER ACTUALLY WAS (the inventory, corrected)

MK9, the best configuration, is a **5-pulsar, 1-source, 13-parameter** problem
(`p_dists = x[8:13]`; 8 CW params + 5 distances — the "coverage 13/13" in the summary is
credible-interval coverage of those 13, **not** fringe correctness), with `toaerr = 1 µs` white
noise, true distances offset from the prior mean by **`DIST_OFFSET_SIGMA = 0.3` σ**, at
**`log10_h = −12.0`**.

**And `snap_to_peak` is NOT a fringe collapse.** Read the code:

```python
def snap_to_peak(x13, pulsar_idx, dL_j, n_pts=50):
    d_lo = max(d_center - 0.6 * dL_j, 1e-6);  d_hi = d_center + 0.6 * dL_j
    ... return d_candidates[np.argmax(lps_scan)]
```

It scans **±0.6 dL** — barely more than one fringe spacing. So it **profiles the CONTINUOUS
within-fringe offset** by deterministic argmax and can, at most, fall into an immediately adjacent
fringe. The **fringe-INDEX exploration was done by the Metropolis moves**, not by the snap: the
`big_jump` (draw a distance from the EM prior, then snap; 20–30 % of proposals) and the
`freq_dist` correlated move.

**The consequence for the hypothesis is sharp — and it is NOT the one this section first drew.**

> ⚠️ **CORRECTION (this claim was wrong for one working day and is now retracted).** An earlier
> draft said the modern target "inherits the snap's within-fringe profile verbatim" via
> `segment_max`. **It does not.** `A2.eval_grid` (DENSE_FRINGE_MAX = 64) returns **fringe centres**
> — one grid point per fringe — for every window holding more than 64 fringes, which is **111 of
> 116 pulsars**; four more get all their centres plus `L0` **padding duplicates**; exactly **one**
> pulsar gets a genuinely dense sub-fringe grid. So `segment_max` is the **identity** almost
> everywhere, and the modern object does not profile the offset — **it PINS the offset at the
> fringe centre (u = 0).** Full anatomy in **§9**.

- The old snap **did** profile the continuous offset (a ±0.6 dL argmax). The modern target does
  something **cruder**: it evaluates only at u = 0.
- The thing the snap did **not** do — sum over fringe *indices* — is what the marginalised
  formulation **does** do (the LSE). That remains the genuine advance.

So: **the snap profiled the offset and Metropolis-sampled the fringe index; the modern target sums
the fringe index exactly and FIXES the offset at the centre.** The two methods are not nested, and
neither dominates: each marginalises one of the two coordinates and butchers the other.

### 8.3 (1) DOES LOCK-ON DEGRADE WITH LOUDNESS AS THE DELTA STORY PREDICTS?

**Not answerable from the old banks — and that absence is itself the finding.** Every MK notebook
(MK3–MK9, `data_likelihood_sandbox_sampling`) injected **one loudness, `log10_h = −12.0`**. MK2
alone mentions −14/−18 (in diagnostics, not a sampler sweep); the root `multiCW` notebook ran −13.0.
**There is no controlled loudness sweep anywhere in the snap-sampler record.** The old "snap works"
claim therefore rests on a **single point in loudness — and it is the loudest point the programme
has ever run**, 1.25 dex above census loudness (−13.25) and 0.5 dex above IGNITE's loudest cell.
The delta-function story does not merely explain why snapping worked; **it explains why the old
record could not have discovered that snapping fails, because it never left the delta regime.**

The prediction is instead tested against the bank that *does* sweep h (§8.1), and it holds in the
direction predicted (concentration ↑ with h) while purity moves the other way.

### 8.4 (2) ARE THERE LOW-SNR WRONG-FRINGE LOCKS RECORDED AS SAMPLER ARTIFACTS?

**No low-SNR sampler runs exist** (§8.3). But the old record contains a **wrong-lock artifact at
the LOUDEST setting**, and it is the R referendum in miniature, written years early. MK2's own
"P5i DIAGNOSTIC: Grid failure analysis", verbatim from its output:

```
  logp at truth CW + truth dists:              -0.22  (ideal)
  logp at truth CW + prior-mean dists:      -2501.79  (grid target)
  logp at grid-peak CW + prior-mean d:      -1393.91  (grid found)
  → Grid found something BETTER than truth CW with prior-mean dists
  → Gap: +1107.9 logp units
  KEY FINDING: ... Prior-mean distances distort landscape — wrong CW fits prior-mean dists better!
  → Need to search distances jointly (e.g., profile out distances at each CW)
```

**A WRONG CW, registered against prior-mean distances, beat the TRUE CW by 1 108 nat.** That is
structurally the same statement as R's `ln Z_plateau − ln Z_needle = +14.19 nat` → `f = 6.9e-7`:
*the wrong-fringe re-registered configuration is preferred over truth.* The old notebook diagnosed
it as a landscape/estimator problem and patched around it (Earth-term-only grid initialisation);
the modern programme proved it is **information-theoretic**. **The artifact was in the record from
the start.** Note also the mechanism named in that note — *a wrong solution that accidentally aligns
with prior-mean distances* — is precisely FORGE's Arm-A/Arm-B prior-pinning finding and IGNITE's
noise-lock, arrived at from three independent directions.

### 8.5 (3) WHAT THE MARGINALISED FORMULATION DOES **NOT** REPRODUCE (explicit, not assumed)

| old machinery (MK5–MK9) | reproduced by `MargJax` + NUTS? |
|---|---|
| `snap_to_peak` (±0.6 dL argmax, within-fringe profile) | **NO — and the modern target is WORSE here.** `segment_max` is the identity for 115/116 pulsars; the offset is pinned at the fringe centre, not profiled (§9). The old sampler handled the offset *better* than the current object does. |
| fringe-index exploration by `big_jump` (EM-prior draw + snap, 20–30 % of moves) | **SUPERSEDED** — fringes are summed exactly (LSE), not proposed. Strictly better. |
| `freq_dist` correlated (f_gw, distance) move | **NOT REPRODUCED as a move** — but the correlation it exploited is *inside* the marginalised gradient (∂/∂log10_fgw of the LSE carries the q-weighted distance response). Untested claim; flag it. |
| **simulated annealing**, `T_start = 30 → 1` over 60 % of burn-in | **NOT REPRODUCED.** NUTS runs at T = 1 only. |
| **empirical-covariance adaptation** (MK9's actual fix: Fisher eigenmodes were 100–15 000× too narrow) | **NOT REPRODUCED as such** — replaced by a dense mass matrix from the target's own Hessian. **See §8.6: this is the one place the old record makes a live prediction about the new machinery.** |
| Fisher-eigenmode / DE-jump / eigen-reopt / CW-joint moves | **NOT REPRODUCED** (NUTS needs no move menu). |
| **no global CW move, no tempering swap** (the old sampler's *known failure*) | **NOT FIXED. The new machinery inherits this failure exactly.** NUTS is a local sampler; a missed mode stays missed. Given `f = 6.9e-7`, this is the same wall, unmoved. |

### 8.6 THE OLD RECORD'S LIVE PREDICTION ABOUT THE NEW PRECONDITIONER

MK9's headline diagnosis (`EXPERIMENT_SUMMARY.md`) was that the **conditional** Fisher
(distances held fixed) gives eigenmode widths **100–15 000× too small** — *"the Hessian measures
sensitivity at fixed distances… but in reality, distances mode-hop to compensate; the effective
posterior is much wider"* — and the fix was to abandon the Fisher for the empirical chain covariance.

**That is a statement about a conditional Hessian, and my target is marginal.** The Hessian of
`lnL_marg` already has the distances summed out, so it *should* reproduce the wide, mode-hopped
geometry that MK9 had to learn empirically — i.e. **marginalisation should dissolve the very
pathology MK9 hacked around.** Corroborating datum from this session: the **exact-distance
(conditional) target in g4 gives `cond(-H) = 3.2e8` with 3 negative eigenvalues at truth** (§3) —
a conditional Hessian misbehaving exactly as MK9 said it does.

**PRE-REGISTERED, before g5 runs:** the marginalised Hessian's eigen-spectrum should be
**dramatically better conditioned than the conditional one**, and its implied proposal widths
should land within ~O(1)–O(10) of an empirical posterior covariance rather than 1e2–1e4 too small.
If instead the marginal Hessian is *also* 1e2–1e4× too narrow, then marginalisation has **not**
bought the geometry, the dense-Hessian mass matrix is the wrong preconditioner, and MK9's empirical
covariance must be ported forward. **Either outcome is a result; this is a falsifiable prediction.**

### 8.7 THE NUMBER g5 IS PRE-REGISTERED AGAINST

*"At what injected h did snap lock-on reliability cross ~90 %?"* — the old banks cannot say (§8.3).
From the bank that sweeps h, the honest decomposition:

- **Snap-VALIDITY** (is `q_p` a delta, so argmax == marginal?): median `q_max` crosses **0.9 at
  `h* = −12.53`** (linear interpolation, 0.631 @ −12.75 → 0.936 @ −12.5; VLBI/T30 cell: −12.67).
  **At MK9's `h = −12.0` — half a dex louder still — the fringe posterior is a delta, and snapping
  is exactly collapsed sampling. The hypothesis is right about why it worked.**
- **Snap-CORRECTNESS** (does the delta sit on the true fringe?): P(confident ∧ true) rises only
  0.049 → **0.197** across the box, capped by a purity that is *falling* (0.785 → 0.366).
  **It never crosses 90 % anywhere in the modelled box, and its ceiling is being pushed DOWN by
  the very loudness that makes snapping valid.**

**The two crossings do not coincide, and that gap is the whole finding:** snapping becomes
*legitimate* at `h ≈ −12.5` and is *reliable* nowhere.

**Therefore the pre-registration for g5:** on a loud realisation (h ≳ −12.5) the fringe-marginalised
posterior and a snap-collapsed one must **agree** — that is the regime the old sampler lived in, and
agreement there *validates the new machinery against the old record* rather than discovering
anything. They must **diverge below h ≈ −12.5**, and the divergence must grow as `q_max` falls.
**If they agree at census loudness (−13.25, median `q_max` = 0.164), something is wrong with my
marginalisation** — that is the falsifier, fixed in advance.

### 8.8 THE ONE-LINE ANSWER TO "WHY DID SNAP-TO-PEAK WORK SO WELL?"

**Because it was run at `log10_h = −12.0` on 5 pulsars with one source and truth 0.3 σ from the
prior mean — the corner of the space where the fringe posterior is a delta sitting on the TRUE
fringe.** Snapping is collapsed sampling exactly there, so it was self-consistent, and the record
contains no loudness sweep that could have revealed the boundary. Push the same object toward
realistic conditions — 116 pulsars, 16 sources, real noise, truth off the prior mean — and the
delta survives (it *sharpens*: median `q_max` → 0.94), while the fringe under it stops being the
true one (purity → 0.37). **What made snapping valid is not what makes it right — and the old
record was built entirely inside the regime where the two are indistinguishable.**

---

## 9. THE SUB-FRINGE OFFSET TEST — **STOP: THE TEST'S PREMISE IS FALSE, AND WHAT IS THERE INSTEAD IS WORSE**

Commissioned as: *"implement the offset-marginalised reduction (segment logsumexp over the
within-fringe offset grid with the correct measure, replacing segment_max)"*, then measure impact
rows from banked numbers.

**It cannot be implemented as specified, because THERE IS NO WITHIN-FRINGE OFFSET GRID.** The
switchable reduction was built and is live (`sampler_core.MargJax(reduce='max'|'lse')`,
`sampler_offset.py`), but a CPU unit test of the reduction algebra against a trapezoid integral
forced a re-read of the grid constructor, and the grid is not what the whole programme (this report
included, §8.2 as first written) assumed.

### 9.1 WHAT `A2.eval_grid` ACTUALLY RETURNS  (`stagec_anchor_a2.py:79-99`, `DENSE_FRINGE_MAX = 64`)

```python
n_fringe = int((hi - lo) / dL) + 1
if n_fringe <= DENSE_FRINGE_MAX:            # <= 64 fringes in the +-3 sigma window
    return np.linspace(lo, hi, B)           #   -> DENSE UNIFORM: real sub-fringe sampling
...                                         # else -> FRINGE CENTRES, one point per fringe
centers = L0 + ks * dL                      #   (thinned to |k| <= 256 + a sparse far tail)
if len(centers) < B:
    centers = np.concatenate([centers, np.full(B - len(centers), L0)])   # PAD WITH DUPLICATES
```

Measured on the pop config (from `grad_diag.npz`, banked):

| regime | pulsars | what `segment_max` does there |
|---|---|---|
| window > 512 fringes → **fringe centres**, 1 point/fringe (K = 512) | **111 / 116** | **NOTHING — it is the identity** |
| 64 < window ≤ 512 fringes → all centres + **`L0` padding duplicates** (K = 139/283/435/455) | **4 / 116** | nothing (max over identical duplicates) |
| window ≤ 64 fringes → **dense uniform** grid | **1 / 116** | genuinely profiles the offset |

**So the sub-fringe offset is PINNED AT THE FRINGE CENTRE (u = 0) for 115 of 116 pulsars.** It is
not profiled (my §8.2 error, retracted) and it is certainly not marginalised. The `374 points in
one fringe` figure reported earlier in this session is **`L0` padding duplicates**, not offset
samples — a second thing the occupancy statistic hid.

### 9.2 WHY THIS IS A PHYSICS PROBLEM, NOT A NUMERICS ONE

The fringe spacing is *defined* by `dPhi_p(dL) = 2*pi`. Therefore

> **every fringe centre carries the SAME pulsar-term phase, mod 2π — the phase at `L0`.**

In **Arm B** (the honest arm) the true distance is `L0 + (n_true + u)·dL` with `u ~ U(-1/2, +1/2)`.
So the model's pulsar term is mismatched against the data's by **`2*pi*u` at EVERY candidate
fringe, identically**. The E-step therefore:

- never evaluates the likelihood at the true distance;
- never evaluates it at *any* local likelihood peak (the peaks sit at `L_true + n·dL`, off the
  centre lattice by the same `u`);
- reads the ~1e-3-nat chirp fringe-breaking contrast **at the wrong pulsar-term phase**.

Marginalising `u` would let each fringe recover its own phase match. Pinning it at `u = 0`
**suppresses the pulsar-term cross-term by a common factor for every fringe** — and a common
suppression of the very term that breaks fringes is exactly the direction of the measured
**Arm A → Arm B price** (FORGE B1.0: yield halves, wrong-certs quadruple). **That price has never
been separated into "truth is off the prior mean" (physics) and "we evaluate at the wrong phase"
(estimator).** They are confounded in every Arm-B number this programme has produced.

### 9.3 THE VERDICT LANGUAGE, APPLIED HONESTLY

The commissioned verdicts were COSMETIC / MATERIAL, to be decided by three deltas. **Neither can be
returned, because the deltas as specified measure the wrong thing:** replacing `segment_max` with a
`logsumexp` **over the same grid** changes nothing for 115/116 pulsars (a logsumexp over one point
is that point) — it would return "COSMETIC" for a reduction that was never running. That would have
been a false negative, and it is exactly what a mechanical execution of the brief would have
produced.

*(A latent bug the gate caught before it ran: a naive `logsumexp` over the existing grid would count
the `L0` padding duplicates up to **374×** for the four padded pulsars. The `max` path is immune.
`sampler_offset.Reductions` must de-duplicate before any `lse` path is trusted.)*

### 9.4 THE TEST THAT ACTUALLY ANSWERS THE QUESTION (re-spec, not run)

The question survives; only the method changes. It needs **new likelihood evaluations on a denser
grid — on the SAME banked realisations, with no new noise draws**:

1. **A sub-fringe grid.** For each pulsar, evaluate `lnL` at `L = L0 + (n + u_m)·dL` for a set of
   offsets `u_m` (e.g. 8–16 points across `u ∈ [-1/2, +1/2]`) at each retained fringe `n`. Cost
   scales the E-step by `len(u_m)`; with the split E-step at ~0.68 s per full-array scan this is
   minutes per realisation, not hours.
2. **Three reductions, one pass**, all on the same evaluations:
   `CENTRE` (u = 0 — the banked object), `PROFILE` (max over u — what the old snap did),
   `MARGINAL` (logsumexp over u with the `du` measure and the pointwise prior — the honest object).
3. **Impact rows**, then, and only then: `lnL_marg(truth)`; R's needle/plateau ratio via the banked
   SMC particles (the reweighting identity in `sampler_offset.py mode=ratio` is correct and ready —
   it just needs the corrected reduction underneath); one Gumbel floor at one **T = 15** surface
   cell (the only cells reproducible on cronus — IGNITE's T-extension machinery is not in this repo,
   §6) and that cell's correct-cert count; each delta against its calibration error.
4. **The pre-registered expectation, stated before running:** `CENTRE` vs `MARGINAL` should agree in
   **Arm A** (where `u = 0` by construction, so pinning is exact) and **diverge in Arm B**, with the
   divergence growing with `|u|`. **If they agree in Arm B too, the pinning is cosmetic and the
   Arm-B price is pure physics. If they diverge, part of the A→B price — and every Arm-B count,
   floor and purity number in the repo — is an ESTIMATOR artefact, and that is criterion-v2.2
   business.**

**This is a STOP, per the mission's discipline. The magnitude is not yet measured; the mechanism is
confirmed by construction (`dPhi(dL) = 2π` is exact, not approximate). Reporting rather than
proceeding.**

### 9.5 WHAT IS **IMMUNE**, AND IT IS THE MOST IMPORTANT PARAGRAPH IN THIS SECTION

**The `u = 0` pinning is EXACT — not an approximation — wherever the true distance sits at the EM
prior mean.** And that is precisely the convention of the entire pre-B1 era:

> `trackB_b1_core.py:32` — *"`generate_injection_params` injects every true distance AT the EM prior
> mean, so `n_true == 0` for all 116 pulsars. **Every Track B deliverable to date inherits that.**"*
> `trackB_estimator_spec.md:406` — *"Truth's fringe integer `n_true == 0` for ALL 116 pulsars (true
> distance == prior mean L0)."*

With `L_true = L0` we have `u_true = 0`, so **the true distance IS a fringe centre**, the pulsar-term
phase at every evaluated point is the correct one, and CENTRE is the exact reduction. Therefore:

| result | u_true | status |
|---|---|---|
| **R's referendum, `f = 6.9e-7`**, the needle/plateau evidences, the break-even 0.188°/src | **0** | **IMMUNE — stands as measured** |
| the **A2 / Anchor-Census ceilings** (J0711 0.953 / J1713 0.989 / J1909 0.998) | **0** | **IMMUNE** |
| **LAMBDA / F2 / L2c** (the cold-start NO-GO, the ladder wall, the pull-in radius) | **0** | **IMMUNE** |
| GEO zero-noise draws, **FORGE Arm A** | **0** | **IMMUNE** |
| **FORGE Arm B**, **IGNITE**, **IGNITE-2**, **D4** — every Arm-B count, floor, purity, wrong-cert | **~ U(−½,½)** | **EXPOSED — magnitude unmeasured** |

**So the exposure is exactly co-extensive with the ARM-B era**, which is to say: with every number
the programme currently treats as *the honest arm*. The tier table, the referendum, and the
information-theoretic NO-GO are untouched. What is on the table is the **A→B price** — and, through
it, the floors, the purity, and the onset surface.

*(One boundary worth stating so it is not over-claimed: R's plateau explores SOURCE parameters off
truth, and the pulsar-term phase depends on the source too. That is the REGISTRATION error F2/L2c
already characterise — a different coordinate from the sub-fringe offset `u`, which is a property of
the DISTANCE truth draw. R's distances are at L0 throughout, so `u = 0` holds across its whole
integration, plateau included.)*

### 9.7 GATE RESULTS (run, banked: `ugrid_phase_true_arm{A,B}.npz`, `ugrid_phase_arm{A,B}.npz`)

**THE MECHANISM IS CONFIRMED — in a CORRECTED form, and the correction matters.**

| gate | statistic | result |
|---|---|---|
| **G-ARM-A** (u_true = 0 by construction) | \|u\*(n_true) − u_true\| | median **0.0312** = the grid resolution → the within-fringe peak sits **at u = 0**. **CENTRE is EXACT in Arm A.** |
| **G-ARM-B** | \|u\*(n_true) − u_true\| | median **0.0184**, 90th pct **0.0294** (resolution 0.0312) → the peak sits **AT u_true**. |
| **G-ARM-B**, the consequence | \|u_true\| | **0.2507 (median)** → **the CENTRE convention evaluates the TRUE fringe off its own likelihood peak by a QUARTER FRINGE (median), i.e. a pulsar-term phase error of ~π/2, on every Arm-B pulsar.** |

> ⚠️ **A CLAIM OF MINE, CORRECTED.** §9.2 said the phase error is *"2πu at every fringe,
> identically"*. **That is too strong and the data say so:** `dL` is defined with the EARTH-term
> `f_gw`, while the pulsar term advances at the RETARDED `f_p < f_gw`, so consecutive fringe centres
> **drift** in pterm phase — and that drift IS the chirp fringe-breaking signal. Measured: across a
> pulsar's fringes the argmax offset has `std_n u* = 0.42` (≈ uniform), because on WRONG fringes the
> model is decorrelated from the data and the likelihood is nearly FLAT in u, so its argmax is noise.
> **The surviving, measured statement is the one that matters: at the TRUE fringe the peak is at
> u_true, and the banked convention never evaluates there.** (My first pass scored the median of u\*
> over *all* fringes, which averages in that noise and printed a spurious "REFUTED" — the log
> retains that line; it is a statement about a mis-specified statistic, not about the mechanism.)

**G-REPRO: FAIL — and it blocks the impact rows.** The FORGE Arm-B recipe does **not** reproduce
from the banked seeds through `trackB_b1_core`: only **81/116** `n_true` match
(`draw_true_distances(seed=dist_seed)`), and **114/114** centre-grid `map_k` disagree; neither
`dlnL` candidate reproduces the banked `dlnL_det` (max |diff| 4.5–6.1 nat). **So my realisation is
not FORGE's realisation**, and no CENTRE→MARGINAL delta may be quoted against banked FORGE/IGNITE
counts until the recipe is recovered. *(FORGE's scorer is not in this repo; the reproduction gate
existed precisely to catch this, and it did — before any impact row was banked.)*

**WHAT THIS LEAVES.** The offset question is **not** answerable against the banked frontier from
cronus today, but it **is** answerable on an internally-consistent ensemble: the same u-grid, three
reductions, on a *freshly generated but fully deterministic* Arm-B T = 15 ensemble (seeds recorded),
with floors refit per reduction on its own nulls. That measures the CENTRE→MARGINAL shift in
q_max / dlnL / cert-count **internally**, which is the physics question; it just cannot be
cross-quoted onto FORGE's numbers. **Verdict therefore remains OPEN — neither COSMETIC nor
MATERIAL is earned.** What is earned is the mechanism (the quarter-fringe off-peak evaluation) and
the immunity boundary (§9.5).

### 9.9 THE FLOOR REFIT — **AND IT OVERTURNS §9.8's HEADLINE. READ THIS BEFORE §9.8.**

> ⚠️ **§9.8 SCORED THE CONFIDENCE LAYER ONLY (`q_max > 0.9`), WITH NO DETECTION GATE. THAT IS NOT
> THE CRITERION.** Rows (ii)+(iii) — 30 fresh nulls (no CW in the data), floors refit **per
> reduction** by the D2 estimator, then the counts re-cut under each reduction's **own** floor —
> change the conclusion. Banked: `SAMPLER_results/ugrid_floors.npz`.

**(ii) THE FLOOR, refit per reduction (Gumbel-MLE, α = 0.05, N = 30 nullN):**

| reduction | floor | fit error | **zero-offender fraction** |
|---|---|---|---|
| **CENTRE** (the banked convention) | **3.59 nat** | ±0.56 | **0.70** — the null FIRES on **30 %** of pure-noise realisations |
| **PROFILE** | 0.00 | ±0.00 | 1.00 — never fires |
| **MARGINAL** | **0.27 nat** | ±0.05 | 0.97 |

**THE CENTRE-PINNING FABRICATES FRINGE EVIDENCE ON DATA CONTAINING NO CW AT ALL.** With no signal
present, evaluating each fringe a quarter-fringe off its own peak manufactures enough apparent
discrimination to clear both criterion layers in 30 % of pure-noise realisations. Under the honest
reduction that collapses to 3 %. **The floor must rise 13× (0.27 → 3.59 nat) to suppress an artefact
the estimator itself creates.** *This is the noise-lock mechanism, caught in the null, with no
signal to confuse it.*

**(iii) THE COUNT, each reduction under ITS OWN floor (10 Arm-B realisations):**

| reduction | floor | detections | certs | correct | **wrong** | certs/real |
|---|---|---|---|---|---|---|
| **CENTRE** | 3.59 | 4 | 4 | 3 | **1** | 0.40 |
| **PROFILE** | 0.00 | 2 | 2 | 2 | 0 | 0.20 |
| **MARGINAL** | 0.27 | 4 | 1 | 1 | **0** | 0.10 |

**THE FLOOR ABSORBS THE ARTEFACT, AND THAT IS THE POINT OF THE FLOOR.** The artefact inflates the
NULL exactly as it inflates the SIGNAL, so a floor refit on nulls scored under the *same* reduction
— which criterion-v2 **mandates** ("floors refit per cell, never inherited") — rises to meet it. The
6×/44 %-wrong catastrophe of §9.8 **does not survive the detection layer**.

**THE VERDICT, SPLIT HONESTLY:**

- **MATERIAL at the CALIBRATION level.** Floor **3.59 → 0.27 nat** (13×, ≈ 6σ of its own ±0.56 fit
  error); null FPR **30 % → 3 %**. The banked floors are **not** measuring what they are believed to
  measure: a large part of `fN` is the *estimator's own off-peak artefact*, not the array's noise.
- **NOT ESTABLISHED at the COUNT level.** With floors refit, 4 vs 1 certifications over 10
  realisations is Poisson-limited (p ≈ 0.2). **The banked FORGE/IGNITE/IGNITE-2 counts are, to first
  order, PROTECTED by their own floor calibration** — their nulls carry the same artefact as their
  signal. The per-cell-refit discipline, adopted for other reasons, silently insured the programme
  against this bug.
- **The residual that does NOT cancel:** wrong-certifications (1 vs 0 here) and the *meaning* of
  the floor. A 13× floor shift changes what "above onset" means even if the counts at a given cell
  survive.

**SO: the cascade to criterion-v2.2 is about the FLOOR and the NULL-FPR interpretation, NOT about
invalidating the banked counts.** The §9.6 re-extraction remains worth doing — but as a
*re-calibration*, not a retraction. **§9.5's immunity table is unchanged** (pre-B1 has `u = 0`
exactly, so no artefact exists there at all).

*(Provenance note, stated because it bounds every number here: this `dlnL` is the fringe-gap
statistic defined in `sampler_ugrid.py`, not FORGE's `dlnL_det` (whose scorer is not in this repo
and which G-REPRO could not reproduce). Absolute floors are therefore NOT comparable to IGNITE's
9–39-nat floors; only the CENTRE-vs-MARGINAL contrast, computed identically on both sides, is.)*

### 9.8 IMPACT ROWS (CONFIDENCE LAYER ONLY — **superseded by §9.9, read that first**)

10 Arm-B, T = 15 realisations (deterministic seeds; 21 s each after build). Same u-grid, same
E-step, three reductions off one pass. Banked: `SAMPLER_results/ugrid_impact.npz`.

**(i) The certification layer (`q_max > 0.9`), 10 realisations, all 116 pulsars:**

| reduction | certifications | correct | **WRONG** | per realisation |
|---|---|---|---|---|
| **CENTRE** — the convention behind every banked Arm-B number | **18** | 10 | **8** | 1.80 |
| **PROFILE** — what the 2023 `snap_to_peak` did | 2 | 2 | **0** | 0.20 |
| **MARGINAL** — the honest reduction | 3 | 3 | **0** | 0.30 |

**Pinning the offset at the fringe centre MANUFACTURES WRONG CERTIFICATIONS.** It returns 6× the
certification count of the honest reduction, and **8 of its 18 certifications sit on the wrong
fringe**; the offset-marginalised reduction returns 3, **all on the true fringe, none wrong**.
The true-fringe MAP count moves the same way — **122 (CENTRE) → 129 (PROFILE) → 135 (MARGINAL)**:
marginalising the offset makes the MAP fringe *more often right*, not less.

**The mechanism, stated plainly.** CENTRE evaluates every fringe a quarter-fringe off its own
likelihood peak (§9.7, median |u_true| = 0.25). That off-peak evaluation is a pseudo-random penalty
that varies fringe-to-fringe, and it **fabricates sharp apparent discrimination between fringes**
— a noise-lock by construction, not by physics. Integrating the offset removes the artefact: the
per-fringe weights collapse toward each other, `q_max` falls, and the confident-but-wrong fringes
lose their manufactured margin.

**(i, magnitudes) Against the criterion's own calibration errors** — this is what earns MATERIAL:

| statistic | shift (MARGINAL − CENTRE) at the certification tail (n = 20) | the yardstick it must beat |
|---|---|---|
| `q_max` | **median −0.230** | the certification threshold is 0.9 |
| `dlnL` (fringe gap) | **median −2.12 nat** | Gumbel floor **fit error ±1.47 nat** (D2, N = 150) |
| MAP fringe | **5 of 20 flip** | — |

Over *all* 1 160 pulsar-realisations the median shifts are small (Δ`q_max` = −0.026, Δ`dlnL` =
−0.046 nat) — **but the criterion does not act on the median, it acts on the tail**, and in the tail
the shifts exceed the floor's own fit error. **COSMETIC is refuted.**

**(ii)/(iii) NOT DELIVERED, and stated as such.** The floor refit and the criterion-v2 count shift
require null realisations re-scored under each reduction, and the counts above are the **confidence
layer only** (`q_max > 0.9`) — the detection layer (`dlnL > max(ln K, floor)`) is NOT applied,
because the floors under MARGINAL do not exist yet. Since `dlnL` shifts by −2.12 nat at the tail
**and the nulls will shift too** (their offenders are precisely the noise-locks this artefact
feeds), the net criterion-v2 count cannot be predicted from these rows — **it must be measured with
floors refit per reduction, never inherited.** That is the re-extraction of §9.6.

**AND THE ARCHAEOLOGY CLOSES ON ITSELF.** `PROFILE` — the reduction the 2023 snap-sampler actually
used, and which the modern pipeline silently dropped when it moved to a fringe-centre grid — is
**closer to the honest answer than the convention that replaced it** (2 certs, 0 wrong, vs 18 certs,
8 wrong). *The old sampler's within-fringe profiling was protective. Discarding it is what opened
this hole.* §8's headline — *"each method marginalises one of the two coordinates and butchers the
other"* — is now priced: butchering the offset is the expensive one.

**SCOPE, honestly.** These 10 realisations are generated from the banked seeds but are **NOT
bit-identical to FORGE's** (G-REPRO fails, §9.7), so the absolute counts are **not** quotable as
"FORGE's Arm-B count was wrong by X". What *is* established, internally and reproducibly: **on
Arm-B data of exactly this class, the banked reduction produces a 6× inflated, 44 %-wrong
certification set that the honest reduction does not.** The A→B price is therefore **part physics,
part estimator artefact**, and the artefact fraction is large enough to change the sign of the
purity story. **This cascades to cronus as criterion-v2.2.**

### 9.6 IF MATERIAL: THE SURFACE RE-EXTRACTION IS CHEAP, AND THIS IS WHY

Every banked realisation is **deterministic in its seeds** (`noise_seed`, `dist_seed`, `geo_id`,
`tol`, `weather` are all banked columns of `ignite_bank.npz` and of the FORGE flats). The u-grid
re-score therefore needs **no new realisations anywhere** — it re-evaluates the SAME data on a
denser distance grid. Cost scales the E-step by `len(u_m)` (~16×; the split E-step is 0.68 s per
full-array scan), so the entire 24-cell IGNITE surface + its 540 fresh nulls is an **O(10) GPU-hour
re-score on an A100**, not a re-run. If the verdict is MATERIAL, that is the ACCRE ask: **re-extract
the frontier under the MARGINAL reduction from the banked seeds, refit the floors per cell (never
inherited), and re-cut the counts** — criterion-v2.2, cronus's call, not this agent's.

---

## 7. EXECUTION RECORD

| item | value |
|---|---|
| box | cronus / RTX 4090 (24 GB), 62 GB host RAM; kyleg absent throughout (checked) |
| env | `jug-gpu` (jax 0.10.1, numpyro 0.21.0), x64, XLA autotune off, persistent compile cache |
| G0 value gate | 2-process run, PASS 5.146e-08 nat (`sampler_gates4.log`) |
| G1a chain-rule gate | PASS 1.131e-15 (`sampler_unit.log`, `grad_unit.npz`) |
| G1b Richardson FD | worst 1.5e-07, FD-limited (`grad_unit.npz`) |
| G1 anatomy | `sampler_graddiag.log`, `grad_diag.npz` (step ladder + fringe-group occupancy) |
| cost | 1.09 s / value_and_grad (F_CHUNK 64); 6.98 s (F_CHUNK 8); naive graph OOM at 16.5 GiB |
| failures, recorded | 3 GPU OOMs during the restructure; 1 host-RAM OOM-kill (§3.0) |
| banked | `SAMPLER_results/{gate_ref_values,grad_diag,grad_unit}.npz` |
| **nothing committed** | no tracked file edited; canonical docs untouched (doc-fence honoured) |
