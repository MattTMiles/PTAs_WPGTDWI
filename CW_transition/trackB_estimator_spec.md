# Track B — B0.1 Estimator specification (spec BEFORE code)

**Deterministic-Annealing EM (DAEM) fringe-identification estimator.**
Reaches the Anchor-Census P(true) ceilings from a COLD start (no truth beyond the EM
priors). Ueda & Nakano (1998), *Deterministic annealing EM algorithm*, Neural Networks
11(2):271 — annealing here is in pulsar-term COHERENCE (Stage A.2), not a generic
temperature.

## 0. Model and notation

- Sources: `N_CW` phase-connected CW, 8 params each
  `theta_s = [cos_gwtheta, gwphi, cos_inc, log10_mc, log10_fgw, log10_h, phase0, psi]`.
  Stack `theta = {theta_s}` (8·N_CW params). GP hyperparams (GWB A/γ, per-pulsar RN)
  FROZEN at truth, as everywhere in Stage C.
- Distances: `L_p`, p = 1..116. Likelihood is `lnL(theta, {L_p})`, the amortised
  zero-noise fast-scan discovery logL (`stagec_fisher.build_fisher_amortised`).
- Fringe grid (per pulsar, per current theta): spacing
  `dL_p = c / [ f_gw (1 − cos μ_p) ]` = `min_s compute_mode_spacing(source s, psr p)`
  (the tightest source sets the grid). Fringe centres `L_p(n) = L0_p + n·dL_p`, n ∈ ℤ,
  `L0_p` = EM prior mean `pe.pdist[0]`.
- EM prior: Gaussian `N(L0_p, σ_EM,p²)`, σ from the canonical `best_distances.txt`
  (+ LIT_INJECT), exactly the A2 priors. Window `±3σ_EM`, `K_p ≈ 6σ_EM/dL_p` candidate
  fringes.

## 1. Latents and the star topology

Latent per pulsar = fringe index `n_p`. Within-fringe offset (the sub-fringe distance)
handled separately by snapping to the local conditional max (§5). **Star topology:
distances are conditionally independent given theta** (all inter-pulsar coupling flows
through the shared source params) ⇒ the E-step factorises exactly per pulsar; NO
cross-pulsar message passing.

## 2. E-step (per-pulsar fringe posterior) — the A2 Bayesian criterion, verbatim

For the current source estimate theta, for each pulsar p evaluate the likelihood at the
`B = 512` fringe-centre distances in the ±3σ_EM window (the existing vectorised scan,
`scan_all`; fringe-centre eval, NEVER the deprecated prominence peak finder):

```
lnL_p(n) = lnL(theta, L_p = L_p(n))          # only L_p varies; other L fixed at current est.
log w_p(n) = β · [ lnL_p(n) − lnL_p(n*) ]  −  (L_p(n) − L0_p)² / (2 σ_EM,p²)
q_p(n)  =  softmax_n log w_p(n)
```

`n*` = current MAP fringe (reference). The Gaussian term is the A2 EM prior tail
(`classify_full`: `logw = (lnL_u − ll) − offs²/(2σ²)`). `β = β(t_c)` is the annealing
coherence weight (§4); at β = 1 this is exactly the A2 `P_true` softmax, and
`P(true)_p = q_p(n_true)`. Star topology ⇒ this is the exact per-pulsar posterior given
theta.

## 3. M-step (soft fringe-weighted source update)

Update theta to maximise the expected complete-data log-likelihood

```
Q(theta) = Σ_p  E_{q_p}[ lnL(theta, L_p(n)) ]  =  Σ_p Σ_n q_p(n) · lnL(theta, L_p(n))
```

- **Truncation:** q_p concentrates on a few fringes; keep the top-K per pulsar carrying
  ≥ 0.999 mass (K small when identified, up to the full window when degenerate — in the
  degenerate/uniform limit the pulsar-term averages out and the pulsar contributes to
  theta mainly through its Earth term, which is the correct behaviour).
- **Extrinsic profiling (linear params):** a monochromatic CW's Earth-term residual is a
  linear combination of {F⁺·cos Φ_E, F⁺·sin Φ_E, F×·cos Φ_E, F×·sin Φ_E}; given the
  NON-linear shape params (sky, log10_fgw, log10_mc, cos_inc) the amplitude/phase/psi
  extrinsics (log10_h, phase0, psi, and the projection) enter linearly and are profiled
  by a whitened linear least-squares solve (QuickCW-style). The remaining NON-linear
  params get jax gradient steps (LBFGS/Adam) on Q. **Fallback if profiling is fiddly on
  the discovery internals:** full-theta gradient ascent on Q (all 8·N_CW), which is
  correct but slower; the profiling is an accelerator, not a correctness requirement.
- **Distances:** represented by the soft mixture over fringe CENTRES during the M-step;
  after convergence each `L_p` is snapped to its local conditional max (§5).

## 4. Annealing schedule (DAEM in coherence)

Per-baseline decoherence prior (Stage A.2): `σ_φ,p² = τ_p / t_c`, τ_p a per-pulsar
coherence time-scale. Map to the pulsar-term coherence weight

```
β(t_c) = exp( − σ_φ² / 2 ) = exp( − τ̄ / (2 t_c) )   ∈ (0, 1]
```

applied to the pulsar-term-induced likelihood variation in the E-step (§2) and,
consistently, to the fringe-varying part of Q in the M-step. **Small t_c ⇒ β → 0**:
pulsar term damped, fringe structure washed out, the objective is near-unimodal in the
Earth-term source params (sky/freq/strain found first). **t_c → ∞ ⇒ β → 1**: pulsar term
re-coheres, fringes condense around the found solution — the full A2 likelihood.

Schedule (geometric ladder, reported): `t_c ∈ {t_c^min · r^k}`, e.g.
`β = [0.05, 0.1, 0.2, 0.4, 0.7, 1.0]` (6 rungs); a few EM iterations per rung, then
polish at β = 1 to convergence. The exact β ladder + iters-per-rung reported with the
B0.2 results.

## 5. Within-fringe snap + convergence

- **Snap:** at each pulsar's MAP fringe, 1-D parabolic/Newton refine of `L_p` on the
  conditional lnL (3-point at centre ± δ, existing 512-eval machinery) → sub-fringe
  distance. Reported distance = snapped `L_p`.
- **Convergence:** stop when BOTH
  `‖Δtheta‖_∞ < tol_θ` (per-param scaled; default 1e-3 of the param’s prior width) AND
  `max_p TV(q_p^new, q_p^old) < tol_q` (default 1e-3). Report EM iteration count.

## 6. Precomputation split (QuickCW-style) — HARD REQUIREMENT

The E-step scan varies ONLY `L_p` (the pulsar-term phase Φ_p(t; L_p)). Everything not
depending on `L_p` — the Earth-term residual, all OTHER pulsars' blocks, and the frozen
outer Cholesky of `(Φ_gwb⁻¹ + blockdiag(FᵀN⁻¹F))` (already cached by
`build_fast_scan_likelihood`) — is computed ONCE per E-step (per theta). Only pulsar p's
pulsar-term projection `F_pᵀ N_p⁻¹ y_p(L_p)` is recomputed per fringe. **The chirp
(pulsar-term frequency evolution over L/c) is RETAINED** — it is the ~1e-3-nat
fringe-breaking signal; a pure-phase-shift precompute would destroy it, so the pterm
block is re-evaluated per fringe (cheap: one pulsar, not the array). Baseline vs split
E-step cost reported on one config (N_CW=16, 116 psr, 512 fringes).

### 6.1 STATUS — split DEFERRED-TO-B1 (2026-07-06, algebra verified against the code)

The low-rank split is NOT implemented for P1. P1's registration map is served by the
batch-only single-vmap E-step (`trackB_p1_map.batched_scan`): all trials share one
fixed `data_obs`, so `fisher_logL` batches over a flat (trial×pulsar×fringe) stack vs
one data leaf — ONE compiled `(4096, n_theta)` shape (avoids attempt-1's pathological
`(D,4096,244)` double-vmap). GATE: batched == D=1 `scan_current` reference **0.00e+00
(bit-exact)** on 8 grid trials. Measured throughput **0.544 ms per full-array eval** →
169-trial 13×13 stage ≈ 91 min. DECISION (2026-07-06): keep batch-only, size P1 grids so
each stage < 45 min (9×9 2-D coarse/zoom1 + high-res 1-D truth cuts + uniqueness patch),
DEFER the split to B1 (≥20 realisations, where the one-time precompute amortises).

**The split algebra, now traced through `discovery/matrix.py` and confirmed exact.**
`fisher_logL` (stagec_fisher.py) is a TWO-LEVEL Woodbury: `make_kernelterms_vary`
(matrix.py:2631, taken because our residuals `y_p(params)` are callable) returns
`(a, b, c)` where the per-pulsar common-rednoise GP is marginalised INNER (its own frozen
per-pulsar factor `cf_p = matrix_factor(FᵀN⁻¹F_p + Π_rn,p)`), and the OUTER global GWB
GP (HD-ORF, `hd_orf`) enters ONLY through the frozen outer factor
`cf_cached = matrix_factor(Φ_gwb⁻¹ + blockdiag(c))`. Distances enter `fisher_logL` at
exactly two per-pulsar-LOCAL places:
```
fisher_logL(θ) = Σ_p a_p(L_p)  +  ½[ b(θ)ᵀ cf_cached⁻¹ b(θ) − ldP − logdet_cf ]
  a_p(L_p) = −½( yᵀN⁻¹y_p − FtNmy_rn,p · cf_p⁻¹ FtNmy_rn,p )   [+ frozen ldN,ldP,logdet_inner]
  b_p(L_p) = TtNmy_gwb,p − TtNmF_p · (cf_p⁻¹ FtNmy_rn,p)         (T = GWB basis, F = RN basis)
```
`c` (terms[2]) is distance-INDEPENDENT (depends only on `cf_p`, `TtNmF`, both frozen) →
`cf_cached` is built ONCE (already is). So sweeping pulsar p's fringe changes ONLY `a_p`
(scalar) and `b_p` (n_gwb vector); every other pulsar's block and `cf_cached` are fixed.
With `M ≡ cf_cached`, `u ≡ M⁻¹ b_base` (one solve per trial) and the per-pulsar blocks
`(M⁻¹)_pp` (precompute ONCE per grid stage, M frozen), each fringe becomes a rank-block
update:
```
Q(L_p) = Q_base + 2 Δb_p·u_p + Δb_pᵀ (M⁻¹)_pp Δb_p ,   Δb_p = b_p(L_p) − b_p(base)
lnL(L_p) = lnL_base + [a_p(L_p) − a_p(base)] + ½[Q(L_p) − Q_base]
```
Per-fringe cost collapses from a full `(npsr·n_gwb)` triangular solve + 116-pulsar delays
to ONE single-pulsar kernelterm (`yprod_p`, matrix.py:2646) + two `n_gwb` dot products
(~1000× per fringe). Per trial: 1 array kernelterm (baseline) + 1 outer solve; per grid
stage: the `(M⁻¹)_pp` precompute. This is the B1 production accelerator; implement with a
bit-exact gate against `batched_scan` before use. **The chirp is retained** — `b_p(L_p)`
re-evaluates pulsar p's residual (Earth+pterm incl. frequency evolution over L/c), so the
~1e-3-nat fringe-breaking signal is preserved; only the array coupling is amortised.

### 6.2 STATUS — split BUILT + all gates PASS (2026-07-07). `CW_transition/trackB_split.py`

Implemented as `class Split` (frozen per-pulsar rednoise Woodbury reconstructed by mirroring
`discovery.matrix.make_kernelterms_vary`, reusing `dsm.*` primitives so the algebra is
bit-identical by construction; `a_const` fixed by difference vs `fisher_logL` at truth). It
exposes `.lnL(theta)`, `.estep(theta_base, EV, AI) -> LNL[npsr,B]`, and the candidate scorer
`.build_tables(source, EV, AI)` + `.score_candidates(assignments[K,npsr]) -> lnL[K]`. Enabling
edits: `stagec_fisher.build_fisher_amortised` now returns an additive `internals` dict (fml/vsm/
cf_cached/ys/frozen/...); `trackB_estimator.build_problem(build_earth=False)` skips the second
amortised compile. The batch-only `trackB_p1_map.batched_scan` stays the correctness ORACLE
(Amendment 6): no shared code.

GATES (pop config, N_CW=16, 116 psr; jug-gpu):
- **0 scalar recon:** `split.lnL == fisher_logL` **0.00e+00** (truth, two plateau, random).
- **1 E-step values:** per-pulsar max|ΔlnL| **1.16e-10** (0/116 > 1e-8); **soft posteriors**
  max|Δq_p(n)| **9e-12** (< 1e-10).
- **2 source gradient:** max rel err **0.00e+00 / 4.05e-13** (frozen M NOT in the differentiated
  graph — safe for the M-step).
- **3 small-signal:** adjacent-fringe ΔlnL agree to **5.8e-11** at a ~1.17e-3-nat step (chirp
  retained).
- **4 scorer values:** 16 random assignments **1.16e-10** (< 1e-8).

PROFILING (honest, Amendment 5): full E-step (116×512) = **0.68 s vs 32 s batch-only = 47×**
(per-pulsar yprod+solve ~0.31 s, rank-update+overhead ~0.37 s — bounded by the single-pulsar
yprod recompute; **NOT 1000× for the E-step**). Candidate scorer: **37 µs/candidate at K=1e5**
(overhead-bound ~1 ms/candidate at K=1e3) — 1e6 candidates in ~37 s, the number the LAMBDA
feasibility probe is budgeted off. `B_CERT=512` (certification never lower); `B_SEARCH` may be
smaller for search phases.

WIRED into `trackB_p2_pipeline.py` behind `--use-split` (default OFF); `soft_estep(..., split=)`
swaps `batched_scan -> split.estep`. Customers = the LAMBDA feasibility probe and B1 (both
separate deliverables — NOT started).

## 7. Multi-start (basins)

≥ 8 cold starts, sources drawn from the generative prior (`generate_injection_params`,
scenario="realistic"; NO truth info). Staged:
1. few coarse DAEM iterations (low-β rungs) for ALL starts;
2. kill basins with `Q` clearly below the running best (gap ≫ per-iter ΔQ);
3. polish survivors (full β ladder → β=1) to convergence.
Report basin statistics: fraction of starts reaching the best optimum (within ΔQ tol),
and whether the ceiling-reproducing basin is found from ≥ 1 cold start.

## 8. Gates

- **B0.2 Asimov (census pop draw, zero-noise, cold):** converged q_p must reproduce the
  A2 ceilings — J0711-6830, J1713+0747, J1909-3744 at P(true) within ~0.05 of census
  (0.953 / 0.989 / 0.998), and NO extra pulsar at P>0.9 the ceiling forbids. Then
  equal-strain N_CW=8, N_CW=16 rows vs the A2/D4 numbers.
- **B1 noisy:** ≥20 batched noise realisations; certified set, P(true) dists for the 3
  ceiling pulsars, and the WRONG-certification count; SBC reliability curve; BH-FDR@0.05.

## Reuse (no re-derivation)
`stagec_anchor_a2.py` (scan_all, classify_full softmax, eval_grid, prep_draw,
make_geometry_injection, LIT_INJECT priors), `stagec_fisher.build_fisher_amortised`
(amortised zero-noise logL + inject_data + fisher_logL_batched), `stagec_hessian` HVP,
`stagec_d4.assemble_hessian`, D4 draw seeds (2000+d, pop 3000). XLA autotune OFF, jax
persistent cache. NEVER the prominence peak finder.

---

## B0.1 BUILD NOTES / CORRECTIONS (empirical, from the debug harness)

1. **Infra:** cache PASS (fisher_logL 2.6s vs 465s cold). E-step machinery PASS:
   `validate-estep` reproduces the A2 population ceilings EXACTLY at truth
   (J0711 0.953 / J1713 0.989 / J1909 0.998, diff 0.000; certified set = the 3).

2. **M-step annealing correction (IMPORTANT).** The coherence-BLEND M-step originally
   specced (`(1-beta)*lnL_earth + beta*lnL_full`) was found to BIAS the optimum off
   truth: the Earth-only (decohered) likelihood's maximum is not at the true source
   params (its residual absorbs pulsar-term power), so ascending it at low beta drags
   the sources away and the certification collapses and does not recover within the
   schedule. Confirmed by the from-truth diagnostic: loudErr drifted 0.15->1.0+, f
   DECREASED, ceilings 0.95/0.99/1.0 -> 0.
   FIX: the M-step ascends the FULL likelihood only, at confidence-gated fringes
   (`q_max > 0.5` to reassign a distance, else hold at the prior-mean/last value). Its
   maximum IS truth (zero-noise Asimov), so truth is a stable fixed point: the
   from-truth diagnostic now holds EXACTLY across the whole beta ladder (f=405413.5
   constant, dphi=0.00, loudErr=0.00, ncert=3, ceilings stable). Annealing therefore
   lives in the E-step temperature beta (posterior flattening) + the confidence gate,
   NOT in a decohered-model blend. NB: temperature on the full likelihood does not move
   its argmax, so pure temperature gives no OPTIMIZER smoothing; the only correct
   smoothing is soft marginalisation of the fringe latents (the expensive soft M-step).

3. **Cold-start basin problem (the open item).** With the corrected full-only M-step,
   a single cold start (sources from the prior) does local ascent on the rugged
   multi-source CW likelihood: it climbs slowly (f 403509 -> ~403900 over 5 iters) and
   stays ~1500 nats below truth (405413), recovering ZERO certifications. This is the
   expected blind multi-source CW recovery difficulty: the annealing smoothing that a
   biased earth-blend provided is gone, and there is no informed initialisation.
   The E/M machinery is sound; what is missing is a GLOBAL SEARCH / SMOOTHING for the
   source parameters. Three candidate fixes (decision pending):
     (A) DAEM SOFT M-step: weight lnL_full over the (temperature-flattened) fringe
         posterior; correct + provides real annealing smoothing, but ~100x costlier
         (grad through sum_p topK fringe variants) -> likely infeasible for B1's >=20
         noisy realisations on the shared 4090.
     (B) F-STATISTIC / matched-filter SEEDING of the loud sources (QuickCW-standard):
         a coarse frequency x sky pre-scan (max over amplitude/phase) places the loud
         sources near truth; the validated M-step then polishes and the E-step
         certifies. Cheap, PTA-standard, aligns with the "QuickCW-style" framing.
     (C) Brute-force many more multi-starts (weakest; may not scale to 16 sources).

---

## B0.2 RESULT — Asimov gate FAILS, and WHY (the make-or-break finding)

**F-stat seeding gate (3a): PASS.** The truth-blind matched-filter scan (47 freq x 192
healpix sky = 9024 pts, amplitude/phase/inc profiled, 393 s) recovers ALL 3 loud sources:
J0711-host loud0 6.4 deg / 0.03 freq-bins / 2F=595; loud1 6.4 deg / 0.05 / 2F=419;
loud2 11.9 deg / 0.64 / 2F=184. All within the M-step capture range; 2F floor ~15.

**B0.2 Asimov gate (pop, seed + EM): FAIL.** The seeded EM converges (basin 1/1) but lands
at median source error ~0.53 (scaled) and the ceilings collapse to ZERO
(J0711 0.000 / J1713 0.002 / J1909 0.000 vs 0.953/0.989/0.998). Not an optimiser bug:

**CERTIFICATION-TOLERANCE SCAN (`perturb`) — the mechanism.** Perturb the truth loud-source
params by a scaled magnitude delta and run the E-step ONLY (P_true at truth+delta):

| delta (scaled) | J0711 | J1713 | J1909 | ncert | mode |
|------|------|------|------|--|--|
| 0     | 0.953 | 0.989 | 0.998 | 3 | all-loud |
| 1e-4  | 0.937 | 0.000 | 0.000 | 1 | all-loud |
| 3e-4  | 0.049 | 0.000 | 0.000 | 0 | all-loud |
| 1e-3+ | 0.000 | 0.000 | 0.000 | 0 | all-loud |
| 1e-4  | 0.954 | 0.957 | 0.999 | 3 | freq-only |
| 3e-4  | 0.941 | 0.000 | 0.000 | 1 | freq-only |
| 3e-3  | 0.646 | 0.000 | 0.000 | 0 | freq-only |

**The certification tolerance on the source parameters is ~1e-4 (scaled)** — e.g. sky
position to ~1e-4 rad (~0.006 deg) for J1713/J1909. It is dominated by SKY/geometry (the
all-loud column dies faster than freq-only); frequency alone tolerates ~3e-4. This is set
by the pulsar-term phase 2*pi*L/dL ~ 1.6e4 rad: any source-param error that shifts this
phase by O(1) rad re-registers the whole fringe pattern, so the TRUE fringe (at L0) is no
longer the posterior mode and P_true collapses.

**Achievable source precision is ~0.5 (scaled)** (seed+EM) — 3-4 ORDERS OF MAGNITUDE
coarser than the ~1e-4 certification tolerance. Even a PERFECT cold-start estimator cannot
close this: no PTA localises a CW source to ~0.006 deg / ~3e-5 in log10 f.

**INTERPRETATION (Track B headline candidate).** The Anchor-Census certification "ceilings"
are computed at EXACT truth source parameters and are a CONDITIONAL artifact of that
conditioning: they collapse under source-parameter errors ~1e-4, far below any achievable
source precision. They are a CEILING (necessary condition), NOT an achievable measurement.
A realistic estimator's source-parameter uncertainty smears the pulsar-term phase across
many fringes and destroys certification. This is the make-or-break unknown Track B was
built to test, answered at ZERO noise (Asimov, best case) — noise (B1) can only make it
worse, so the noisy arm as specified would confirm the same at higher cost.

**IMPLICATION for B1 and the census result.** The census "3 data-driven certified pulsars"
is a conditional ceiling, not a forecast of what a pipeline achieves. Before running the
B1 noisy ensembles, the science question changes: the right quantity is the JOINT
source+distance posterior (marginalising source-parameter uncertainty into the fringe
posterior), not the conditional-on-truth ceiling. DECISION PENDING with Matt.

---

## P1 STATUS AT SESSION END (registration-map, joint-registration pivot)

- **P1 attempt 1 ABANDONED at 64/169 coarse trials.** mode_p1 coarse grid (13x13 sky
  scan of loud0 about truth, +/-20 deg, D=4 padded batched E-step). Per-trial cost
  ~33 s under shared-GPU contention (kyleg co-tenant) -> ~90 min/grid, untenable for the
  staged zoom (coarse -> zoom1 -> zoom2). Timeout at 2280 s (64 trials).
- **Partial baseline saved: `trackB_p1_partial.npz`** — trial_sky_positions (169x2,
  deterministic, recomputed) + timing provenance. **profiled_lnl / registration_count
  were IN-RAM only and NOT persisted before the kill -> stored NaN/-1, UNRECOVERABLE.**
  The next fast-path rebuild must recompute these and reproduce the grid positions
  exactly (deterministic), then this npz is superseded. Log: `logs/p1_coarse_attempt.log`.
- **Why slow:** p1_eval_trials loops trials in Python, each calling A2.scan_all
  (116x512 double-vmap) + a per-trial Python fringe_posterior (116 pulsars). The D=16
  batching (attempt 0) hit a pathological (16,4096,244) double-vmap compile; D=4 padded
  (attempt 1) compiles but is still ~33 s/trial.
- **Next session plan (batched rebuild, run as nohup background w/ npz checkpoints):**
  1. Make the TRIAL the runtime axis of ONE compiled graph: batch trial source-params as
     data-tuple args (QuickCW split), so the per-grid-stage E-step is one vmapped call,
     not a Python trial loop. Target << 1 s/trial.
  2. Vectorise fringe_posterior (drop the 116-pulsar Python loop -> array ops over
     (trials, pulsars, fringes)).
  3. Staged zoom coarse->zoom1->zoom2 with an npz checkpoint after EACH stage
     (trial_sky_positions, profiled_lnl, registration_count, n_registered_true), so a
     kill never loses computed trials again.
  4. Gate: batched rebuild must reproduce any recomputed coarse trials to 1e-8 vs a
     reference scan_current(D=1) eval before trusting the map.
- **Decisive question P1 still must answer:** is there a unique JOINT-lnL / true-fringe
  registration NEEDLE at truth, width ~ the 1-rad prediction (2*pi*L/dL ~ 1.6e4 rad ->
  ~0.006 deg)? secondary partial-registration side modes? does the needle need all 116
  pulsars or do the certified few carry it? Then P2 (joint solve + B0.2 re-gate) / P3 doc.

---

## STEP 2 — the RTK / LAMBDA integer solve (P2 verdict: continuous methods insufficient)

**Why (P1 + probes + P2, 2026-07-07).** The census ceilings correspond to a UNIQUE global
optimum at truth (P1 needle), but the objective is a sharp CUSP on a noisy PLATEAU (integer
fringe ambiguity). Continuous optimisers fail: the HARD-snap joint is a gauge-conspiracy
plateau (probe 1); the SOFT marginal is a sharp cusp at truth (STEP 1a line scan) but a
gradient-ascent M-step (a) climbs the fixed-q Q that the ~98 degenerate pulsars displace and
(b) falls off the cusp tip onto the plateau (test-mstep: 0.05 scaled drift, census->0).
Fixing the fringe INTEGERS turns the cusp into a smooth quadratic in the source params ->
the source localises to mas by Newton. This is exactly GPS-RTK carrier-phase ambiguity
resolution; the estimator is an integer least-squares (LAMBDA), not gradient descent.

### The correspondence (paper architecture figure)

| GPS / RTK carrier-phase positioning | PTA phase-connected CW distance/certification |
|---|---|
| carrier phase (ambiguous, precise)  | pulsar-term phase Phi_p = 2*pi*L_p/dL_p |
| integer ambiguity N_i               | fringe index n_p (which cycle of the pulsar term) |
| code / pseudorange (coarse, unambiguous) | F-stat / Earth-term matched filter (degrees, no fringes) |
| float solution (real-valued N)      | soft-EM source+distance estimate (~deg; the plateau landing) |
| float covariance Q_N                | fringe-index covariance from the Fisher / soft posterior |
| LAMBDA Z-transform (decorrelate)    | decorrelate the fringe-integer covariance before search |
| integer least-squares (bounded search) | bounded enumeration of fringe integers in the float neighbourhood |
| conditional baseline given fixed N  | source sky/freq re-solved with carrier distances fixed (smooth quadratic) |
| ratio test (best vs 2nd-best)       | certification: dlnL(best vs runner-up fringe set) vs ln K, and A2 P_true |
| partial ambiguity resolution (fix the confident subset) | fix only the needle-CARRIERS (q_max>0.9, census-3 first); keep the ~98 degenerate pulsars soft |
| baseline vector (position)          | source localisation (sky/freq) |
| fixed-solution mm positioning       | mas-class source localisation == distance certification (same measurement; exchange rate = the lever arm 2*pi*L/dL) |

### STEP 2 plan (partial-ambiguity, feasibility-gated)

- **2a Partial set.** Resolve integers ONLY for the concentrated carriers (emergent q_max>0.9,
  bootstrapped census-3 first). Never search the full 116 x K lattice; the degenerate majority
  stays marginalised (they carry Earth-term info, not fringe info).
- **2b Search.** Float covariance from the existing Fisher/soft-posterior; Z-transform to
  decorrelate; bounded enumeration over the float neighbourhood; each integer candidate scored
  by the batched E-step (0.544 ms/eval) after a conditional source re-solve; certification =
  ratio test cross-calibrated against the A2 Bayesian P_true.
- **2c Feasibility probe (GO/NO-GO).** Seeded by STEP 1b's actual float landing zone: does the
  bounded integer search + conditional source re-solve recover TRUTH's carrier integers with a
  ratio-test margin? Asimov first; report search volume, wall-clock, margin. THEN decide the
  full build.

---

## LAMBDA feasibility probe — RESULT (2026-07-07). VERDICT: NO-GO (earned at F2; see amendment)

**READ THE F1/F2 AMENDMENT BELOW FIRST.** The initial NO-GO (this section) was seeded by a
truth-blind soft-EM float that was CONTAMINATED (an M-step compaction changed the gradient; the
loud2 slot was mis-seeded; the EM walked good seeds off immediately). Those contaminate the L0
float and everything the L2 BLIND search read from it. The verdict SURVIVES but is re-derived
cleanly in the amendment via the lever-arm ladder (F2), independent of the contaminated float.
BANKED and unaffected: the L2b/L2c fixed-integer pull-in result (< 1e-4 scaled, cond ~4e10,
truth a real sharp neg-def max). Read this section for the machinery + L2b/L2c; read the
amendment for the corrected anatomy and the earned verdict.

Truth-blind, zero-noise Asimov, pop config (N_CW=16, 116 psr), jug-gpu + `--use-split`.
Files: `trackB_L0_float.py`/`.npz`, `trackB_L1_spec.py`, `trackB_L2_probe.py`/`.npz`,
`trackB_L2b_basin.py`/`.npz`, `trackB_L2c_newton.py`/`.npz`; logs/L0_float.log, L2_probe.log,
L2b_basin.log, L2c_newton.log. Truth's fringe integer `n_true == 0` for ALL 116 pulsars
(true distance == prior mean L0), so "matches truth" = winner integer 0 (diagnostic only).

### RTK correspondence, instantiated with the measured numbers

| GPS / RTK stage | PTA instantiation | Probe measurement |
|---|---|---|
| code / pseudorange | F-stat / Earth matched filter | loud0/1 seeded 6.4 deg, loud2 17.7 deg (its 11.9-deg capture-seed fell outside the top-16, outranked by noise seeds 2F 285-350) |
| float solution (real N) | truth-blind soft-EM landing (L0) | soft-EM DIVERGES: loud sources -> 4.31 scaled off truth by it19; ZERO carriers (q_max spectrum tops out 0.52); census MAP fringes n=163/222/-2 (truth 0), q_max 0.29/0.17/0.31, ZERO posterior mass on n_true |
| float covariance / Z-transform | fringe-index covariance | not built (probe scale): census candidate sets (>1% mass) already span n in [-182,226] and DO NOT contain n=0 -> Z-transform moot, truth outside the neighbourhood |
| integer least-squares (bounded search) | staged `score_candidates` at float source (L2) | stage-1 census lattice 11x26x13 = 3718 (padded to 28175 for the overhead floor); winner integers -113/-143/-19 (NOT 0), margin 0.198 nat, matchesTruth=False |
| conditional baseline given fixed N | source re-solve, integers FIXED (L2/L2b/L2c) | LBFGS from float: stays 101.9 deg / 4.31 scaled off. **Pull-in radius test (integers FIXED at TRUTH): 0** — neither LBFGS nor Newton locks (sky<0.006 deg) from ANY offset >=1e-4 scaled |
| fixed-solution mm positioning | mas-class source == certification | NOT reached: P_true(census) = 0.000/0.000/0.000 (ceilings 0.953/0.989/0.998); certified set = {} |
| ratio test | dlnL(best vs runner-up) | 0.198 nat at WRONG integers; dlnL(n_true vs best-wrong) at the final wrong source = -13.6/-10.9/-8.4 (truth fringe is BELOW the wrong fringe once the source is off) |

### The mechanism (L2c, decisive): the needle EXISTS but its basin is curvature-microscopic

At EXACT truth (integers fixed, 24 loud source params): `|grad|_inf = 0.000`, Hessian
eigenvalues in **[-5.965e11, -14.43]** — strictly NEGATIVE-DEFINITE. So truth IS a genuine
local max (the P1 needle is real; the L2 auto-label "needle-absent" was an LBFGS artifact and
is WITHDRAWN). The objective is locally quadratic, but the condition number is **~4e10**: the
registration direction (sky/freq, driven by the pulsar-term phase 2*pi*L/dL ~ 1.6e4 rad) has
curvature ~6e11 while the faint-source / nuisance directions are ~O(10). The quadratic bowl is
valid only within `delta ~ sqrt(2/6e11) ~ 1e-6` scaled; beyond that the pulsar-term phase wraps
to an adjacent fringe and the local max Newton climbs is a DIFFERENT, fringe-shifted source
position. Measured (`trackB_L2c`, damped-Newton, integers fixed at truth, faint frozen):

| start offset (scaled) | Newton lands (sky deg) | dlnL vs truth | locked? |
|---|---|---|---|
| 1e-4 | 0.046 | -622 | no |
| 3e-4 | 0.114 | -923 | no |
| 1e-3 | 0.285 | -936 | no |
| 3e-3 | 0.834 | -965 | no |
| 1e-2 | 1.71 | -922 | no |
| 1e-1 | 17.7 | -1406 | no |

Even from 1e-4 scaled (~the certification tolerance) with the CORRECT integers, Newton lands
~0.05 deg out and ~600 nat below truth. **Conditional-re-solve pull-in radius < 1e-4 scaled**
(curvature/lever-arm-limited, ~1e-6 in principle).

### Verdict (pre-registered language: iii FAIL)

The LAMBDA / integer-least-squares cold-start is **INFEASIBLE on zero-noise Asimov**, and the
FAIL is structural, not an optimiser bug — it fails at BOTH gates the method requires:
1. **Integer neighbourhood.** The truth-blind continuous float lands O(1)-few scaled off truth
   (here 4.3, diverging; even a charitable ~0.5 landing shifts the pulsar-term phase by
   hundreds of fringes), so the bounded integer search around the float MAP integers (163/222/-2)
   never contains truth's n=0 for the sky-sensitive carriers. Blind search fixes WRONG integers.
2. **Source pull-in.** Even GIVEN truth's integers, the conditional source re-solve locks only
   within < 1e-4 scaled of truth (registration cusp curvature ~6e11, lever arm 2*pi*L/dL ~ 1.6e4);
   the float is 4+ orders of magnitude outside this basin. Fixing integers does NOT deliver a
   navigable smooth quadratic from any achievable float.

This is the same precision cliff B0.2 found (certification tolerance ~1e-4 vs achievable source
precision ~O(1)), now shown to defeat the DISCRETE integer solve as well: neither the integer
neighbourhood nor the source basin is reachable from the cold float. The needle is real but its
basin of attraction is set by the pulsar-term lever arm and is unreachably small. GPS-RTK works
because the code/pseudorange float lands within one carrier wavelength; the PTA "pseudorange"
(F-stat, ~deg) lands ~1e4 carrier cycles away, and no fixing step spans that. NO-GO for the full
LAMBDA build. B1 noisy ensembles would only confirm this at higher cost.

**Caveat / one avenue not closed.** The soft-EM float here DIVERGED (loud -> 4.3 scaled, zero
carriers) rather than settling near the ~0.5 scaled B0.2 reported; a float held closer to truth
(e.g. frozen at the F-stat seed ~6-18 deg, or a stronger source-localisation prior) would still
fail gate 2 (pull-in < 1e-4 << 0.1 rad ~ 6 deg) but would give the integer search a less
degenerate seed. This does not change the NO-GO (gate 2 is the hard wall), but it is the only
place a future attempt could look; it is not worth building.

---

## LAMBDA probe — F1/F2 AMENDMENT (2026-07-07, verdict CONTESTED then RE-EARNED)

The float used above DIVERGED (loud -> 4.3 scaled), contaminating the L0/L2-blind chain. Three
bounded follow-ups. Files: `trackB_F1a_gate.py/.npz`, `trackB_F1c_descent.py/.npz`,
`trackB_F2_ladder.py/.npz`; logs/F1a_gate, F1c_descent, F2_ladder.

### F1 — L0 divergence anatomy (the float was contaminated three ways)

- **F1a (compaction NOT exact — my "mathematically exact" claim WITHDRAWN).** One M-step, full
  348-config scan vs compact 72, identical state: **|dQ| = 4.20e4 (rel 4.4e-3), max|dsrc| = 0.60
  raw = 0.20 scaled** (phase0). F1a's first-guess root cause (nan_to_num-after-sum poisoning) was
  TESTED and REFUTED in the close-out (see REPAIRS LOG, Repair 1): the fix is a no-op, gradients
  agree to 2.2e-16 at one step -> the real cause is multi-step FP-order AMPLIFICATION through the
  ~4e10-conditioned Adam (numerical chaos), NOT NaN poisoning. Either way the compact M-step used
  for the L0 float is a DIFFERENT trajectory than the pipeline's full path -- the "exact" claim is
  withdrawn; the M-step is simply numerically unstable.
- **F1b (NOT a label switch; the "4.3" is amplitude).** L0-endpoint fitted-vs-truth loud sky 3x3:
  best permutation 130.5 deg vs naive-diagonal 138.0 deg (marginal -> no clean swap). The 4.3
  scaled max is **log10_h (amplitude, extrinsic, irrelevant to fringe registration)**: per-loud
  scaled offsets log10_h = 3.23/1.96/4.31. REGISTRATION-critical offsets (sky, freq) are
  src0 ~0.18, src1 ~0.10, **src2 ~1.2** scaled -> two loud sources held near the seed, src2 alone
  walked off. Co-source confusion: two fits crowd truth1 (8 deg, 12 deg), truth0/truth2 left with
  only distant fits (28 deg, 47 deg).
- **F1c (walk-off is IMMEDIATE + seeding gap).** The F-stat seed cleanly places only 2 of 3 loud
  (loud0/1 at 6.4 deg = 0.056/0.045 scaled sky); the THIRD loud slot is a SPURIOUS 2F=350 noise
  seed at 108.95 deg (loud2's real 17.7-deg seed is buried at rank-13 among the faint fillers ->
  `seeds_to_phi` keeps only the top-16 by 2F). The soft-EM then walks the two good seeds OFF
  IMMEDIATELY (it0: 6.4 -> 18-19 deg), and by it4 the source-to-truth assignment PERMUTES (sources
  swap nearest-truth: co-source confusion, the flagged P1-conditioning risk). The descent never
  approaches ~2 deg -- the coarse F-stat floor IS ~6.4 deg and EM makes it worse. So the L0 float
  is a product of (mis-seeding + immediate co-source walk-off + the compact-vs-full gradient
  change), NOT a trustworthy "landing zone."

### F2 — the lever-arm ladder spectrum (the criterion the probe skipped; EARNS the NO-GO)

Per (pulsar, loud-source) registration tolerance = scaled source offset that shifts the
pulsar-term phase `Phi_p = 2*pi f L_p (1-cos mu_ps)/c` by 1 rad (auto-diff of Phi wrt the source
sky/freq params, divided by the param scale). All 348 pairs finite. **Sky (the binding axis):**

| sky tol (scaled) | loosest | median | tightest |
|---|---|---|---|
| value | 1.85e-3 | 3.81e-5 | 1.34e-6 |

Pairs with sky tol >= 0.05 (F-stat float floor): **0**; >= 1e-2: **0**; >= 1e-3: 3; >= 1e-4: 44.
Frequency tolerances are looser (loosest 2.52e-2, median 1.04e-4) but IRRELEVANT: the sky float
error (0.05 scaled ~ 6.4 deg) alone shifts the phase by `0.05 / 1.85e-3 ~ 27 rad` at the LOOSEST
pulsar (~1300 rad at the median) -> total de-registration regardless of frequency precision. So
sky binds.

**The ladder does NOT span.** The registration ladder (1.85e-3 -> 1.34e-6, internally continuous,
max gap 0.30 dex) sits ENTIRELY below the achievable blind sky float. There are ZERO rungs between
the float floor (~0.05, or even a charitable 0.03) and the loosest rung (1.85e-3): a **~27x
(1.4 dex) wall at the TOP of the ladder**. The GPS-RTK wide-lane cascade fixes the loosest baseline
first, tightens, then descends the ladder -- but here NO baseline is loose enough to register from
a degree-class sky float, so the cascade CANNOT FIX ITS FIRST RUNG. It never bootstraps.

### F3 — NOT run (pre-registered gate: cascade attempted ONLY if F2's ladder spans; it does not).

### EARNED VERDICT: NO-GO (iii FAIL), re-derived independent of the contaminated float

The two ends bracket the impossibility, both measured directly (no reliance on the diverged float):
- **TOP (F2):** the loosest pulsar-term registration baseline tolerates only 1.85e-3 scaled sky;
  the best blind sky localisation (F-stat) is 0.05 scaled (6.4 deg). 27x wall, zero rungs above.
  No first rung can be fixed -> the wide-lane cascade cannot start.
- **BOTTOM (L2b/L2c, BANKED):** even GIVEN the exact integers, the conditional source re-solve
  pull-in is < 1e-4 scaled (cusp curvature ~6e11, cond ~4e10; truth a real neg-def max). 
The ENTIRE registration ladder lives at <= 1.85e-3 scaled while any achievable float lives at
>= ~0.03-0.05 scaled -- a clean >=1.4-dex separation, with the razor cusp at the bottom. GPS-RTK
works because its code/pseudorange float lands within one carrier cycle; the PTA F-stat float
lands ~27-1300 pulsar-term cycles from the loosest baseline and there is no looser lane to bridge
the gap. NO-GO for the full LAMBDA build; F1's M-step instability (see REPAIRS LOG: it is
numerical CHAOS in the ~4e10-conditioned Adam, NOT the NaN-poisoning F1a first guessed) and the
loud2 seeding gap are real but SEPARATE issues that do NOT change the ladder wall. Do NOT start B1.

---

## TRACK B CLOSE-OUT — REPAIRS LOG (2026-07-07)

### Repair 1 — the soft-EM M-step gradient (F1a follow-up): NOT what F1a claimed

F1a found the 7.7x-compaction path (compact-72) and the full-348 path disagree on one M-step by
|dQ|=4.20e4, max|dsrc|=0.20 scaled, and attributed it to `nan_to_num`-AFTER-the-config-sum
poisoning (one NaN config, incl. a w=0 one, zeroing the whole gradient). **That diagnosis was
WRONG.** Applied the fix (per-config `nan_to_num` BEFORE the weighted sum, trackB_p2_pipeline.py
`_cfg_grad`) and re-ran the gate: the result was **BYTE-IDENTICAL** to the pre-fix baseline
(Q_full=9456390.744704, |dQ|=4.20e4, max|dsrc|=0.605) -> the fix is a NO-OP -> there are NO NaN
config gradients to neutralise. The real anatomy (trackB_step1_diag): **[0]** the compact and full
NON-ZERO config sets are IDENTICAL (sorted weights match, |dsum(w)|=1.8e-15); **[1]** at ONE Adam
step full-vs-compact `max|dsrc| = 2.2e-16` -- the GRADIENTS agree to machine zero. So the F1a
divergence is **multi-step AMPLIFICATION** of a machine-epsilon FP summation-order difference
(full interleaves the nonzero configs with hundreds of exact zeros; compact sorts by weight)
through the ~4e10-conditioned Adam M-step -- i.e. the soft-EM M-step is **numerically CHAOTIC**,
not NaN-poisoned. This is consistent with F1c's chaotic label-permuting descent. The `nan_to_num`
relocation is RETAINED as hygiene (correct in principle; harmless here) but is **NOT credited as
the repair**. CONSEQUENCE: the soft-EM source trajectory (and any wall-height ceiling read off it)
is TRAJECTORY-SENSITIVE at the ~0.2-scaled level; the Step-3 number below is quoted with that
uncertainty. This does not touch F2 (geometry) or L2b/L2c (fixed-integer) -- both float-independent.

### Repair 2 — the F-stat seeder (STEP 2): sky-only sidelobe suppression

The single-source F-stat map throws cross-frequency sidelobes of the loud sources (2F 285-350)
that the old sky-AND-freq NMS kept (freq differed), burying loud2's true peak at rank-13. Fix
(trackB_estimator.py `pick_seeds`): SKY-ONLY exclusion at 25 deg around each accepted stronger
peak. GATE (pop draw, trackB_step2_seeder): 3/3 loud captured; loud2 promoted from rank-13/17.70deg
to rank-11/11.88deg (its true peak); loud0 6.37deg, loud1 6.39deg. Sweep: excl in [15,30]deg all
give 3/3 (chose 25); excl <=12deg over-suppresses (2/3). Data-driven; no truth in selection.

### Repair 3 / STEP 3 — the wall's height (repaired float, both repairs in)

One truth-blind float (repaired sky-only seeder 3/3: loud0/1/2 at 6.37/6.39/11.88 deg; M-step
nan_to_num hygiene; 30-iter budget). The soft-EM does NOT converge to truth -- it CHAOTICALLY
walks the loud sources OFF the F-stat seed (co-source confusion + the ~4e10 M-step instability):
registration-critical (label-matched to nearest of all 16 fitted sources) sky offset climbs from
the seed floor ~0.05 scaled (it0) to ~0.21 scaled (endpoint), transiently touching one census
certification (J1713 P_true=0.977 at it1) then dropping it. **Repaired-float REGISTRATION CEILING
= 0.05 (seed floor, best) to 0.21 (endpoint) scaled sky** (loud0 0.098 / loud1 0.208 / loud2
0.130 at the endpoint; median 0.11), TRAJECTORY-SENSITIVE at the +/-0.2 level flagged in Repair 1.
Extrinsics (log10_h ~2, phase0/psi ~O(1) scaled) are large and SEPARATE (they do not register
fringes). **WALL = ceiling / loosest-rung = 0.05..0.21 / 1.85e-3 = 27x..112x = 1.4..2.05 dex,
ZERO rungs in the gap.** Matches the pre-registered de-risk prediction (float floor 0.01-0.05 ->
0.7-1.4 dex; the chaotic endpoint pushes it to 2 dex). The repaired seeder is a clean improvement
(loud2 11.88 deg vs 17.70; 40 candidate seeds vs 501; census briefly reachable) but does NOT move
the wall -- the soft-EM cannot hold the seed, and even the seed floor is 27x above the loosest rung.

## TRACK B TERMINAL VERDICT (2026-07-07)

**NO-GO for COLD-START POINT ESTIMATION of CW distance certification on this array** -- earned,
bracketed FLOAT-INDEPENDENTLY:
- **TOP (F2, geometry only):** loosest per-(pulsar,loud) sky registration baseline 1.85e-3 scaled;
  ZERO baselines tolerate >= 1e-2. Repaired-float ceiling 0.05..0.21 scaled -> WALL 1.4..2.05 dex,
  ZERO rungs. The GPS-RTK wide-lane cascade cannot fix its first rung.
- **BOTTOM (L2b/L2c, BANKED, float-independent):** even GIVEN exact integers, conditional source
  re-solve pull-in < 1e-4 scaled; truth a genuine sharp neg-def max (eig [-5.96e11,-14.4]), cusp
  curvature ~6e11 (lever arm 2*pi*L/dL~1.6e4).
- Continuous methods floor 100-500x short (P2); the discrete solve cannot reach its first rung (F2).
  The census ceilings (J0711/J1713/J1909, P 0.953/0.989/0.998) remain REAL but CONDITIONAL --
  achievable only given source knowledge at the ~1e-3-scaled level THIS ARRAY CANNOT SELF-GENERATE.

**TERMINAL CHARACTERIZATION -- RESOLVED by deliverable R (posterior referendum, 2026-07-08):
the NO-GO deepens to INFORMATION-THEORETIC.** The honestly fringe-marginalised POSTERIOR does NOT
concentrate at the needle -- it smears across the plateau. Referendum number **f = Z_needle /
(Z_needle + Z_plateau) = 6.9e-7** (see the DELIVERABLE R section below). Point estimation was
already NO-GO; R shows INFERENCE does not rescue it either -- prong-3 (guided/anchored sampling)
cannot win from the data alone, because the marginal evidence itself prefers the wrong-fringe
re-registered plateau over truth. The transient J1713=0.977 hint is real per-pulsar but overwhelmed
in the aggregate evidence (the plateau MAP certifies census-3 at P_true = 0/0/0). Only MORE DATA (or
external sky knowledge; see the design theorem lever iii) helps.

## DELIVERABLE R -- THE POSTERIOR REFERENDUM (2026-07-08). VERDICT: (ii) INFORMATION-THEORETIC NO-GO.

Module `trackB_R_referendum.py` (jug-gpu, x64, autotune off; Matt commits, never committed here).

**The object -- honest count-once STAR-TOPOLOGY fringe-marginalised likelihood:**
`lnL_marg(theta) = lnL_ref(theta) + sum_p LSE_n[ (LNL_p(n) - lnL_ref) + lnprior_p(n) ]`,
theta = LOUD source params; faint (13) FIXED at truth; base distances = L0; LNL_p(n) via the split
E-step; lnprior_p(n) = -(n dL_p)^2/(2 sig_p^2) (lit priors, A2 Gaussian tail). This is
log(sum over the fringe lattice) under the additive/star approximation; lnL_ref counted ONCE (the
multi-counted `sum_p logsumexp[full lnL]` is NOT a valid evidence -- it counts lnL_ref 116x, and is
what the STEP-1a linescan's -8000/0.003 numbers were). GATES (all PASS): split==fisher 0.00;
batched E-step==oracle 1e-10; exact lnL_ref==fisher_logL(truth) 0.00; exp(-m_p)==A2 P_true census
(0.953/0.989/0.998) to 3.7e-8; additive-approx error 0.14 nat at the 116-pulsar MAP config.
Identity: `lnL_marg(truth) = lnL(truth) - sum_p ln P_true_p = 405413.51 + 272.83 = 405686.34`
(the 272.83 nat is the fringe ENTROPY the plateau competes with).

**ACTIVE-DIM reduction (measure-consistent, tractable).** f is computed over the SKY registration
dims (2/loud, 6 total); freq + extrinsics (amp/phase/psi/inc) FIXED at truth. Justification: their
Laplace-marginalisation factors are ~equal on needle and plateau (source detected across the whole
+-2 deg basin) and CANCEL in f; F2 established SKY is the binding registration axis (loosest sky tol
1.85e-3 vs freq 2.5e-2). The full 24-D broad-extrinsic box is intractable (a narrow good-fit ridge
collapses the SMC acc to 0.04 AND dilutes the plateau). Flagged favorable conditioning -- the real
(less-conditioned) case is only WORSE, so verdict (ii) is conservative-safe.

**R1 -- Z_NEEDLE (truth-placed local quadrature; diagnostic truth use LABELLED).** FD Hessian over
the 6 sky dims (batched adaptive phase-aware steps) -> Laplace + per-eigenvector 1-D quadrature
(batched, doubling-gated). **ln Z_needle = 405617.64 (quadrature) / 405617.84 (Laplace)** -- agree
0.2 nat, doubling-gate STABLE (0.000), 6/6 positive-curvature eigs, needle sigma_sky 2.4-9.7e-6
scaled (razor-sharp; marginalisation does NOT broaden the L2c cusp), r_gauss~0.75-1.11. Z_needle =
g0 - 68.7. (The +-6sigma-corner separability residual 11.5 nat is at weight e^-36 -> immaterial;
quad==Laplace confirms robustness.) `trackB_R_znaddle.npz`.

**R2 -- Z_PLATEAU (tempered SMC, vectorised through the split; cov-preconditioned RW).** Box:
truth-centred +-2 deg sky/loud (the 6-deg F-stat seed offset is an nside=4 healpix GRID artifact,
NOT a float bias -> truth-centred; this box doubles as the EM-targeted / known-sky scenario). Needle
sub-box excised (0 particles hit it). N=256, 2 seeds, ESS-adaptive tempering (beta->1 in 5 stages),
per-stage checkpoints. **ln Z_plateau = 405631.83 +- 0.053** (2 seeds agree to 0.05 nat -> evidence
ROBUST; the high-beta acc dip to 0.01 biases Z_plateau DOWN if anything -> f smaller -> verdict
firmer). Plateau ceiling L~405659 (28 nat below the needle peak 405686 but 40 nat ABOVE the needle's
integrated evidence 405618); posterior collapses to the high-L shell just outside the needle.
`trackB_R_zplateau_s{0,1}.npz`, `trackB_R_zplateau_summary.npz`.

**R3 -- THE REFERENDUM NUMBER.** ln Z_plateau - ln Z_needle = **+14.19 nat** ->
**f = Z_needle/(Z_needle+Z_plateau) = 6.9e-7**. The needle is 28 nat higher at its PEAK but its 6-D
sky volume fraction (~e^-73) loses to the plateau's +-2 deg volume by 14 nat. **BREAK-EVEN sky prior
theta* = 0.188 deg per loud source** (V*/V_2deg = 6.9e-7): the needle wins ONLY if the sky is
pre-localised to < 0.19 deg/source -- 32x tighter than the F-stat blind floor (~6 deg); an EM
counterpart (host galaxy) supplies exactly this (design-theorem lever iii). Census-3 P(true fringe)
at the posterior peak (plateau MAP) = **0/0/0**. `trackB_R_referendum_result.npz`.

**R4 -- VERDICT (ii), pre-registered.** f = 6.9e-7 <= 0.05: the honestly marginalised data do NOT
prefer truth; the NO-GO is INFORMATION-THEORETIC, not search-only. Prong-3 does not unpark on the
data alone. Design theorem + E-track (lever ii: do eccentric harmonics populate the empty lane band?)
+ lever (iii) EM-targeted sky are the whole forward story. Only more data / external sky helps.

**R5 T-crossing pencil (banked, `PENCIL_t_crossing.md`), quoted with the verdict.** The Earth-term
blind float ceiling reaches the F2 loosest registration rung (1.85e-3 scaled) only at **T ~ 11 kyr**
(baseline sigma_sky ~ T^-1/2) or, at fixed T, per-source **SNR ~ 289** (h ~ 10^-12.3, ~27x the
current loudest source); the L2c pull-in (1e-4) at T ~ 3.75 Myr. Consistent with the >=1.4-dex
F2/L2c wall, now expressed as a timescale -- the wall does not close on any observational horizon.
The referendum's break-even 0.188 deg is the INFERENCE-side counterpart of this pencil's
point-estimation wall: both say the current 15-yr array cannot self-localise the sky finely enough,
by optimisation OR by marginalisation.

**R POSTMORTEM (bounded, CPU-side on the banked plateau samples; `trackB_R_postmortem.npz`).** Q:
does any single pulsar's EM prior + lever arm near-decide its OWN fringe even while the joint sky is
NOT? Per-pulsar over the 512 plateau samples (both seeds): distinct best-fringe indices, marginal
P(true fringe), and effective wrong-solution count exp(H[fringe posterior]).

| pulsar | K_cens | sig_EM_pc | #distinctMAP | MAP=0 frac | P_true(plateau=full) | exp H | needle P_true |
|---|---|---|---|---|---|---|---|
| **J0437-4715** | **4** | **0.25** | 5 | **0.887** | **0.568** | **3.41** | 0.682 |
| J2222-0137 | 20 | 0.97 | 9 | 0.066 | 0.125 | 12.1 | 0.577 |
| J1909-3744 | 72 | 3.0 | 37 | 0.021 | 0.019 | 11.0 | 0.998 |
| J0711-6830 | 558 | 21.2 | 66 | 0.146 | 0.117 | 34.6 | 0.953 |
| J1713+0747 | 1264 | 40.0 | 68 | 0.025 | 0.014 | 33.9 | 0.989 |

**Answer: YES, exactly one -- J0437-4715.** On the plateau (joint f=7e-7), 89% of samples still MAP to
the TRUE fringe (k=0), P_true=0.57 (barely below its needle 0.68), exp H=3.4. Mechanism: tight
composite EM prior (0.25 pc -> only K=4 candidate fringes) x nearby-pulsar loose lever arm -> its
fringe is robust to the joint-sky wander. The CENSUS-3 (J0711/J1713/J1909), which certify AT truth
(needle P_true 0.95/0.99/0.998), SMEAR on the plateau (P_true 0.12/0.014/0.019, exp H 11-35, 37-68
distinct wrong fringes) -- their fringe info is DATA-driven (loud-source-broken, needs the joint sky
right), not prior-pinned, so it evaporates once registration fails. J2222 (K=20) is intermediate
(smears; prior tight in pc but K not small enough). Clean separation: PRIOR-PINNED fringe info
(J0437) survives joint failure; DATA-DRIVEN fringe info (census-3) does not. This is design-theorem
lever (iii) made concrete -- J0437 (the sole Anchor-Census K<=3 pulsar) is the seed a conditional /
anchored pipeline can fix from its prior alone, independent of the failed joint sky. It does NOT
rescue the referendum (one pulsar cannot phase the array; f stays 7e-7) but it names the anchor.

## THE DESIGN THEOREM (constructive -- a paper section, not a limitation)

Certification requires registration LANES below ~2e-3 scaled that a blind float can reach, OR a
blind float above the loosest existing lane. tol ~ 1 / (2*pi f L_p (1-cos mu_ps)/c) per sky axis.
Three constructive levers:

**(i) WIDE LANES FROM NEARBY PULSARS (measured).** tol ~ 1/L_p. The 3 loosest F2 rungs are the
NEAREST pulsars: J0711-6830 (0.106 kpc, tol 1.85e-3), J1630+3734 (0.089 kpc, 1.39e-3), J0437-4715
(0.155 kpc, 1.02e-3) -- vs the median array pulsar L=1.38 kpc (farthest 25 kpc). ~13x nearer ->
~13x looser lane. A dedicated ~100 pc, well-timed target sub-array raises the top lanes ~10x toward
the float ceiling; build the certification target list by crossing the readable-pulsar-term
sub-array WITH DISTANCE (J0711 is already a census pulsar; J0437 the anchor K<=3).

**(ii) WIDE LANES FROM ECCENTRIC HARMONICS (E-track's sharpened job).** For an eccentric binary the
beat lanes at the orbital FUNDAMENTAL sit at n-times the tolerance of the n-th emitting harmonic
(fundamental lever arm = 1/n of the harmonic's). NEW SHARPENED OBJECTIVE for the E-track: measure
whether eccentric-source lanes POPULATE the 2e-3 -> float-ceiling (~0.05) gap the circular ladder
leaves empty.

**(iii) LOWERING THE WALL FROM ABOVE.** Better Earth-term localisation shrinks the float floor
toward the lanes: EM counterparts (a known sky removes the sky axis entirely), the anchor-paper
external-distance regime, and the inverse-direction Q1 computation (source localisation FROM
certified distances). Each buys dex on the float side.

Certification is feasible where lever (i)/(ii) raises a lane above the lever-(iii)-lowered float.

## B1 STEP 1 (2026-07-09) -- THE TARGETED (f, mc) LADDER. CASCADE DOES NOT IGNITE AT ANY TIER.

B1's brief specified the targeted pipeline as "sky FIXED at truth, FREQUENCY free (1-D scan to
needle resolution), extrinsics fit". The B1 PILOT (`trackB_b1_pilot.py`) shows that presumes a
false premise. Files: `trackB_b1_core.py` (gated machinery), `trackB_b1_pilot.py`,
`trackB_b1_ladder.py`; `b1_pilot_m{1,2,3,4}.npz`, `b1_ladder.npz`, `b1_loopgain.npz`,
`b1_cascade_tiers.npz`; logs/b1_{core_gate,pilot,ladder,cascade_tiers}.log. Zero-noise Asimov,
truth-placed (LABELLED). Env jug-gpu, autotune off, x64.

### The premise that failed: log10_mc is a THIRD registration axis

`make_phase_connected_binary` builds the pulsar term at the RETARDED time with the full PN chirp
(`evolve_phase(tp)`, discovery/deterministic.py:522). Its phase therefore depends on log10_mc.
Measured on the pop draw (tau_p median 3.52 kyr, min 0.019, max 86.0): the exact pterm-minus-earth
GW phase has median 1.302e4 rad (max 1.735e5); the chirp changes it by a median factor 0.955 vs
the non-chirped 2 pi f tau. **F2's ladder (`trackB_F2_ladder.py:47`) used the NON-chirped phase
`2 pi f L (1-cos mu)/c`, whose gradient wrt log10_mc is IDENTICALLY ZERO -- F2 was structurally
blind to an mc lane and could not have found one.** Reproduction gate: the pilot's exact-phase
autodiff reproduces F2's chirp-free freq ladder to 3 digits (loosest 2.518e-2 / median 1.035e-4
vs F2's 2.52e-2 / 1.04e-4), so the disagreement is the chirp term, not the code.

**M1 -- 1-rad registration ladders (SCALED, loosest / median / tightest):**

| axis | loosest | median | tightest | # >= 1e-2 | # >= 1e-3 |
|---|---|---|---|---|---|
| sky | 1.852e-3 | 4.472e-5 | 3.847e-6 | 0 | 3 |
| log10_fgw (chirped) | 2.518e-2 | 1.435e-4 | 2.062e-5 | 2 | 27 |
| **log10_mc** | **6.047e+1** | **7.836e-4** | 1.345e-5 | 53 | 161 |
| freq, F2's chirp-free object | 2.520e-2 | 1.035e-4 | 5.578e-6 | -- | -- |

Including the chirp TIGHTENS the freq ladder by a median 0.72x (F2's freq numbers were mildly
optimistic; its SKY numbers are unaffected and stand).

**M2 -- per-axis CERTIFICATION tolerance** (E-step only at truth, one axis perturbed at a time;
the delta at which the census-3 P(true) collapses). Directly comparable to B0.2's ~1e-4 scaled,
which perturbed all 8 params at once and could not attribute:

| axis | cert tol (scaled) |
|---|---|
| sky | 1e-5 (J1713 0.989 -> 0.018) |
| log10_fgw | 3e-5 |
| log10_mc | 1e-3 |
| extrinsics (cos_inc, log10_h, phase0, psi) | NO collapse to 1e-2 (census 0.953/0.989/0.998 flat) |

So the 8 source params split CLEANLY: 4 extrinsics (harmless), 2 sky (supplied by the counterpart),
and **(log10_fgw, log10_mc) = the two registration axes a targeted pipeline must still find.**
This also MEASURES the assumption R1/R2 asserted (that extrinsic Laplace factors cancel in f).

**M3 -- what the Earth term actually supplies** (Earth-only Asimov, HVP-assembled Hessian, symmetry
residual 4.32e-15, inverted WITH the generative priors -- never pinv, which would report sigma=0 for
an unconstrained axis):

| axis | sigma_prior | posterior loud0/1/2 | info gain |
|---|---|---|---|
| log10_fgw | 1.636 | 1.008e-2 / 4.201e-3 / 6.092e-3 | 162x / 389x / 269x |
| **log10_mc** | 1.732 | 1.002 / 1.727 / 1.731 | **1.73x / 1.00x / 1.00x** |

sigma(log10_mc) = 0.501 / 0.864 / 0.866 dex. For loud1 and loud2 the posterior IS the prior: the
Earth-term chirp (Delta_phi_E = pi fdot T^2 ~ 0.05 rad for the loudest) is below the noise.
Per-axis targeted wall: freq 203x, **mc 1730x**.

### STEP 1A -- the joint ladder: dense, but NO FIRST RUNG

Pulsar a registers when its pulsar-term phase wander over the current (f,mc) box falls below 1 rad:
`R_a = 1 / sqrt( sum_k [(sig_f,k g_f,ak)^2 + (sig_mc,k g_mc,ak)^2] )` (R_a >= 1 = registers now).

- **95-100% of the wander is log10_mc** (J0711 95.0%, J1713 100.0%, J1909 99.9%, J0437 99.3%).
- max R = 2.143e-2 (J0711-6830, 0.106 kpc); median 2.036e-4; min 7.441e-6.
- **IGNITION: 0 pulsars at R >= 1.** Loosest rung needs a 47x box shrink; the binding census
  pulsar (J1713) needs **3305x**.
- The ladder SPANS internally (42 rungs between float and census target, max log10 gap 0.387 dex
  < F2's 0.7 criterion) -- so, unlike sky, the failure is NOT rung spacing. **It is the mc BOX.**

### STEP 1B -- loop gain (conditional (f,mc) Fisher at a FIXED fringe subset)

Fixing a subset S of fringes and re-solving (f,mc) IS a masked likelihood. `MaskedDelay` makes the
per-pulsar pulsar-term mask a RUNTIME arg, so F_S = -H[lnL(theta; data_S, mask_S)] over the 6
registration params is one call per subset -- no optimiser, hence no L2c pull-in problem. Asimov
data injected with the SAME mask => residual(truth)=0 => H = -Fisher exactly.

| subset fixed | sig_f | sig_mc | sig_f/tol | sig_mc/tol | next-rung R = GAIN |
|---|---|---|---|---|---|
| earth only (mask=0) | 1.592e-3 | 1.710e+0 | 53.1 | 1709.6 | **0.04** |
| top-1 (J0711) | 1.347e-3 | 2.204e-1 | 44.9 | 220.4 | 0.12 |
| top-3 | 1.013e-3 | 5.766e-2 | 33.8 | 57.7 | 0.15 |
| top-6 | 3.709e-4 | 7.720e-3 | 12.4 | 7.7 | 0.19 |
| top-12 | 2.720e-4 | 4.240e-3 | 9.1 | 4.2 | 0.36 |
| top-24 | 1.407e-4 | 1.196e-3 | 4.7 | 1.2 | 0.70 |
| **top-48** | 5.073e-5 | 2.663e-4 | 1.7 | 0.3 | **1.34** |
| top-116 | 1.363e-5 | 1.530e-5 | 0.5 | 0.0 | -- |
| **census-3 only** | 1.403e-4 | 1.053e-3 | 4.7 | 1.1 | **4.50** |

**KILLED HYPOTHESIS (recorded, not quietly dropped).** Pre-registered guess: since
`R_a ~ 1/(sig_mc g_mc,a)`, a pulsar is registrable exactly to the degree it is BLIND to mc, so no
fixing order can bootstrap ("the wide lane is wide because it is uninformative"). **REFUTED by the
table**: fixing J0711 alone cuts sigma_mc 7.8x; the census-3 alone cut it 1600x. Loose rungs carry
mc information through their AMPLITUDE/SNR, not only through g_mc.

**The real structure is BISTABILITY.** Loop gain is 0.04 at the float, crosses 1 between 24 and 48
fixed fringes, and is 4.50 given the census-3. A self-consistent registered state EXISTS and is a
strong attractor (this is the P1 needle, seen from the (f,mc) side); it is simply unreachable from
the Earth-term float. Certification is SELF-REFERENTIAL: ~30 registered fringes are needed to
measure the chirp mass that lets any fringe be registered.

### STEP 1C -- cascade seeded from each pre-registered conditioning tier

| tier | conditioning | sigma_f (scaled) | sigma_mc (scaled) | best free R | ignites? |
|---|---|---|---|---|---|
| A | sky only | 0.0101/0.0042/0.0061 | 1.002/1.727/1.731 | 0.0214 | NO |
| B | + EM period, sigma_P/P = 1e-3 | 0.00145 x3 | unchanged | 0.0220 | NO |
| C | + host D_L (sigma_mc 0.036/0.022/0.021 dex) | 0.00145 x3 | 0.0727/0.0435/0.0409 | **0.271** | NO |

Tier C's mc box is set by the CW amplitude relation
`log10 h = (5/3) log10 Mc + (2/3) log10 f - log10 D_L + const`, i.e. by the ARRAY's own
sigma(log10_h) (0.060/0.036/0.034 dex from M3) plus sigma(log10 D_L)=0.005 -- NOT by anything the
counterpart supplies beyond the redshift.

- **An EM period buys nothing** (B vs A: 0.0214 -> 0.0220). Frequency is not the binding axis;
  a 7x tighter f moves the best rung by 3%. Consistent with the 95-100% mc share.
- **Only the distance moves anything**, and it moves it 12.7x against the 47x needed -- Tier C
  misses ignition by **3.7x**. (A pre-run prediction of "misses by 1.4x" was WRONG: it scaled R by
  the MEDIAN mc shrink (34x); R is a quadrature over the 3 loud sources and the loosest rung J0711
  is dominated by loud0, whose sigma_mc shrinks only 13.8x because loud0 has the worst Earth-term
  sigma(log10_h). The binding shrink is the WORST source's, not the median's.)

### STEP 1 VERDICT + the three doc items

**The wide-lane cascade cannot start in the targeted scenario either.** The ladder has no first
rung at any physically-supplied conditioning tier. Per the pre-registration (Matt, 2026-07-09) the
cascade is therefore NOT B1's pipeline; STEP 2 (the (f,mc) referendum at tiers A/B/C) decides
whether any tier is legitimate conditioning for B1.0-1.4.

(a) **F2 CHIRP-BLINDNESS.** F2's ladder is VALID for sky (reproduced to 3 digits) and blind to mc
    BY CONSTRUCTION. Its freq ladder is mildly optimistic (median 0.72x). Its sky verdict and the
    Track B TERMINAL VERDICT are unaffected -- both rest on the sky axis and on L2b/L2c.

(b) **R RETROACTIVE NOTE.** R's ACTIVE-DIM reduction fixed freq AND log10_mc at truth and argued
    their Laplace factors cancel between needle and plateau. For log10_mc that argument is WRONG
    (needle sigma ~1e-3 scaled vs plateau ~1.7). The error is FAVORABLE: including mc as an active
    dim shrinks Z_needle by its extra ~e^-(tens) volume fraction, so **f = 6.9e-7 stands and is
    FIRMER than reported.** The blind verdict (ii) is unchanged. What does NOT follow, and what R
    never claimed but the design theorem's lever (iii) implies, is that removing the sky makes
    certification possible -- (f, mc) then carry their own volume contest. That is STEP 2.

(c) **PHYSICAL IDENTIFICATION: mc-registration IS the evolution/timestamp channel.** The pulsar
    term is a kyr-baseline timestamp of the source's phase. What registers it is not the source's
    frequency but its frequency EVOLUTION: the fringe index is set by the accumulated chirp over
    tau_p, so the array must know fdot (hence Mc) to place any fringe. The 22-yr Earth term cannot
    measure fdot (Delta_phi_E ~ 0.05 rad). This is exactly the E-track's eccentric-harmonic
    mechanism in CIRCULAR form -- the harmonics there and the chirp here are the same physical
    handle (frequency structure resolved over the pulsar-term lag). **The E-track and the targeted
    pipeline have merged**: lever (ii) is not an alternative to lever (iii), it is the missing
    ingredient lever (iii) needs.

## B1 STEP 2 + 1D + 2C (2026-07-09/10) -- THE TARGETED REFERENDUM, THE LAST DOOR, AND THE
## BREAK-EVEN CURVE. VERDICT: THE TARGETED SCENARIO HAS ITS OWN INFORMATION WALL.

### HEADLINE (the mechanism, not the number)

**Conditioning the (f, mc) PRIOR BOXES barely moves the evidence -- the plateau does not fill them.
What moves the evidence is LIKELIHOOD STRUCTURE, which the eccentric harmonic comb supplies and no
prior box can.** Measured across three tiers and a four-point response curve, below.

Files: `trackB_b1_referendum.py`, `trackB_b1_softcascade.py`, `b1_step2_table.py`;
`b1_referendum_tier{A,B,C}.npz`, `b1_softcascade.npz`; logs/b1_ref_tier*.log,
logs/b1_softcascade.log. Zero-noise Asimov; truth-placed needle quadrature LABELLED.

### The object and the tiers

R's count-once star-topology fringe marginal, verbatim, with the SKY of every loud source FIXED
at truth (the EM-counterpart baseline) and ACTIVE dims = (log10_fgw, log10_mc) x 3 loud = 6.
Extrinsics fixed at truth -- and this is now MEASURED, not asserted: the pilot's M2 shows census
P(true) is flat in (cos_inc, log10_h, phase0, psi) out to delta = 1e-2 scaled, so their Laplace
factors do cancel between needle and box (the justification R1/R2 assumed).

Identity reproduced: `lnL_marg(truth) = 405686.3434 = lnL(truth) 405413.51 + fringe entropy 272.83`.

### TIER C (sky + EM period + host D_L) -- the most generous physically-defensible conditioning

`sigma(log10_mc) = (3/5) sqrt(sigma_logh^2 + ((2/3) sigma_logf)^2 + sigma_logDL^2)`, i.e. set by
the ARRAY's own sigma(log10_h) plus the host redshift -- NOT by anything else the counterpart
supplies. FROZEN at 4 seeds (see the gate convention below):

| quantity | value |
|---|---|
| ln Z_needle (bracketed quadrature, `\|dlnJ\| = 0.055` STABLE, 6/6 brackets closed) | 405629.6337 |
| ln Z_box (4 seeds: 405633.906 / 405632.985 / 405631.878 / 405633.369) | 405633.035, std 0.859, s.e.m. 0.429 |
| d = ln Z_needle - ln Z_box | **-3.4008 nat** |
| **f(C)** | **0.0323 +- 0.0134**; 2-sigma band [0.0055, 0.0591] |
| pre-registered verdict (f - 2sig >= 0.95?) | **FAIL by 16.1x** |
| gap / sampler scatter | 3.96x |
| break-even lambda_mc = exp((d - logit 0.95)/3) | 0.1206 -> deficit 8.29x **[SUSPENDED -- see the saturation note]** |

**~97% of the targeted posterior's evidence sits on the wrong-fringe plateau even with the sky
exact, the frequency pinned by an EM period, and the chirp mass constrained by the host distance.**

A and B have strictly wider mc boxes (the generative prior; ln V_box = -7.72 vs C's -17.82), so
they cannot pass where C fails; they are measured for the table, not for the verdict.


### SATURATION NOTE (2026-07-09) -- THE DEFICIT FACTOR IS SUSPENDED, THE VERDICT IS NOT

The break-even above assumed `Z_box ∝ lambda_mc^3` (the uniform-plateau-density approximation R
used for theta* = 0.188 deg). **Tier A's SMC falsifies that assumption for these boxes.** Tier A's
6-D box is e^{10.1} ~ 24000x larger in volume than Tier C's, yet

    ln Z_box(A) = 405632.432  (density 405640.148 + lnV -7.716)
    ln Z_box(C) = 405633.035  (density ~405651.7   + lnV -17.816)

agree to within the 0.86-nat sampler scatter. Z_box is an INTEGRAL, so Z_box(A) >= Z_box(C)
necessarily; measuring them equal means **the plateau's evidence has SATURATED -- it is confined
well inside Tier C's box** (consistent with the log-scan: lnL_marg falls 919 nat by the box edge).
Enlarging the mc box adds volume carrying negligible likelihood, so Z_box stops responding.

CONSEQUENCES, separated:
- **f(C) = 0.0323 +- 0.0134, FAIL by 16.1x: UNAFFECTED.** A measured ratio of two measured
  evidences; no extrapolation enters it.
- **The soft-cascade FAIL: UNAFFECTED.** Independent of any box argument.
- **DEFICIT = 8.29x: SUSPENDED.** It extrapolates Z_box across exactly the range over which Z_box
  is now shown to be insensitive. The EXISTENCE of a deficit is safe; the NUMBER is not.
- The two-term price (sigma_h x11.3 AND sigma_logDL <= 1.0%; or kappa >= 8.29) inherits the
  suspension: those numbers are functions of the deficit and must be requoted after it is measured.

THE FIX IS A MEASUREMENT, NOT AN ARGUMENT: run Z_box on a deliberately SHRUNKEN mc box and read the
true break-even off the response curve Z_box(lambda_mc). Until then, quote only: no static tier
concentrates the targeted posterior (f(C) = 0.032, failing 0.95 by 16x), and the soft cascade
cannot bootstrap.

**TIER-GRADIENT (the emerging headline).** Conditioning the (f, mc) PRIOR BOXES barely moves the
evidence at all, because the plateau does not fill them. What moves it is LIKELIHOOD STRUCTURE --
which is what an eccentric harmonic lever buys and what a prior box cannot.


### THE THREE-TIER TABLE (frozen 4-seed protocol; `b1_step2_table.npz`)

| tier | conditioning | ln Z_needle | ln Z_box | d (nat) | f | +-2sig | gate |
|---|---|---|---|---|---|---|---|
| A | sky only | 405629.637 | 405632.017 | -2.380 | 0.0847 | +-0.131 | FAILED |
| B | + EM period | 405629.634 | 405632.619 | -2.986 | 0.0481 | +-0.0227 | FAILED |
| C | + host D_L | 405629.634 | 405632.734 | -3.101 | 0.0431 | +-0.0369 | FAILED |

(Tier C frozen at 4 seeds reads f = 0.0323 +- 0.0134; the table's auto-ingest used the completed
5-seed npz, f = 0.0431 +- 0.0185. Both FAIL identically; the discrepancy is recorded, not resolved
by preference.)

- **Z_needle is tier-independent to 0.003 nat** (405629.6367 / .6337 / .6337): three boxes, two
  independent bracket algorithms, one local integral. This also retroactively confirms that the
  first (broken-quadrature) f(A) = 0.769 was inflated by ~2.2 nat.
- **The tier gradient is FLAT, and mildly INVERTED**: f = 0.085 -> 0.048 -> 0.043 as conditioning
  TIGHTENS. Counterpart information does not help. Z_box(A) >= Z_box(C) is required (an integral
  over a 24000x larger volume cannot shrink); measured -1.02 +- 0.95 nat, i.e. consistent with
  EQUALITY at ~1 sigma. That equality IS saturation.
- All three tiers: gate FAILED on the range statistic, verdict FAIL on the +-2sigma band by
  4.4x-13x. Two instruments, two statements.

### STEP 2C -- THE BREAK-EVEN RESPONSE CURVE (`trackB_b1_breakeven.py`, `b1_breakeven_curve.npz`)

Z_box measured at lambda_mc in {1 (banked), 0.3, 0.12, 0.05}, 2 seeds each, f-box held at Tier C's.
Needle excision 0.0% at EVERY lambda -> all four points are clean plateau evidences.

| lambda_mc | mc box (dex, median) | ln Z_box | +-sem | f |
|---|---|---|---|---|
| 1.0 | 0.0652 | 405633.035 | 0.429 | 0.032 |
| 0.3 | 0.0196 | 405631.754 | 0.616 | 0.107 |
| 0.12 | 0.0078 | 405630.535 | 0.174 | 0.289 |
| 0.05 | 0.0033 | 405628.910 | 0.664 | 0.673 |

**(a) SATURATION SCALE.** Z_box begins responding at lambda ~ 0.3: **the plateau's own chirp-mass
extent is ~0.02 dex.** A newly measured quantity -- the width of the wrong-fringe plateau in Mc.
(The knee is marginal at 0.3: a 1.28-nat drop against a 1.23-nat threshold; the response is
unambiguous by lambda = 0.12.)

**(b) TRUE BREAK-EVEN.** ln Z_box never falls to the f = 0.95 target (405626.689) inside the
measured range -> **lambda* < 0.05: a BOUND, not a value.**

**(c) CORRECTED DEFICIT: > 20x**, replacing the suspended 8.29x. **The suspension was vindicated:
the real price is WORSE than the proportionality implied.** Even shrinking the mc box 20x reaches
only f = 0.673 -- it still fails. (A log-linear extrapolation of the last two points would give
lambda* ~ 0.015, deficit ~66x. NOT REPORTED AS A RESULT: extrapolating through an assumed scaling
is exactly the error this curve exists to correct.)

### THE PRICE, REQUOTED OFF THE CURVE

sigma(log10_mc) delivered by Tier C = 0.0364 / 0.0217 / 0.0205 dex. Break-even needs a box below
lambda* < 0.05 of that, i.e. **sigma(log10_mc) < ~0.003 dex, a > 20x improvement.**

Setting sigma_h -> 0 leaves a FLOOR of 0.00301 dex set entirely by the assumed sigma(log10 D_L) =
0.005 dex (host z to ~1%). At the corrected (>20x) deficit that floor is at or ABOVE the requirement
for every loud source: **strain alone cannot close it, and neither can any counterpart-supplied prior
box.** The two-term price (sigma_h x11.3 AND sigma_logDL <= 1.0%) was computed off the 8.29x deficit
and is SUPERSEDED; requote as: sigma(log10_mc) must improve by > 20x, which no combination of
(position, period, host distance) delivers.

**WEAVE's kappa >= 8.3 at e >~ 0.6 INHERITS THE SUSPENSION.** The MECHANISM stands -- an eccentric
harmonic comb supplies likelihood structure a prior box cannot. The THRESHOLD re-prices off the
curve: kappa must now exceed the >20x bound, not 8.29. The e-value that delivers it is the E-TRACK's
Fisher calculation, NOT made here. Do not quote an eccentricity until the E-track produces it.

### THE (SUPERSEDED) PRICE ARITHMETIC, kept for the audit trail

| | loud0 | loud1 | loud2 |
|---|---|---|---|
| sigma(log10_h) delivered by the array | 0.0604 | 0.0359 | 0.0337 |
| sigma(log10_mc) delivered (Tier C) | 0.0364 | 0.0217 | 0.0205 |
| sigma(log10_mc) needed for f = 0.95 | 0.0044 | 0.0026 | 0.0025 |

Setting sigma_h -> 0 leaves a FLOOR of **0.00301 dex**, set entirely by the assumed
sigma(log10 D_L) = 0.005 dex (host z to ~1%). That floor is BELOW loud0's requirement but ABOVE
loud1's and loud2's by 1.15x and 1.22x. **So the deficit is not closed by strain alone: even a
perfectly measured amplitude leaves the two quieter loud sources short, because the host distance
becomes the binding term.** Closing it requires BOTH:
- sigma(log10_h) improved ~11.3x. With sigma_h ~ 1/SNR ~ T^{-1/2} that is T x 128 ~ **2840 yr**;
  at fixed T it is ~1 dex louder in strain (log10_h ~ -12.2 vs the population's -13.25).
- sigma(log10 D_L) <= 0.0044 dex, i.e. **<= 1.0% in D_L** (vs the 1.16% assumed).

Or, substituting for the pair: an eccentric-harmonic lever **kappa >= 8.29**. (The e-value that
delivers kappa = 8.29 is the E-TRACK's measurement, NOT made here; do not quote an eccentricity
until the E-track's Fisher map produces it.)

### STEP 1D -- THE SOFT-CASCADE PROBE (the last door), pre-registered, one probe one verdict

STEP 1's R_a >= 1 is the HARD criterion (when can a fringe be LOCKED). It does not ask what a
SOFT, posterior-weighted mixture leaks: at R ~ 0.27 a soft fix spreads weight over ~4 fringes
rather than locking one, and a 4-fringe mixture still carries chirp-mass information. Tier-C
conditioning, Asimov, <= 5 iterations, NO hard fixes anywhere. Mc update = the width of the
mixture-marginalised posterior (`lnL_marg` + tier prior). Truth (n_true = 0) scores W only; it
never steers the loop.

| iter | sigma(log10_mc) dex, loud0/1/2 | S = sum_p max_n q_p - S_prior | W = sum_p [1 - q_p(n_true)] | # pulsars whose W GREW |
|---|---|---|---|---|
| 0 | 2.03e-3 / 4.00e-4 / 5.0e-5 | +12.525 | 114.821 | -- |
| 1 | 2.5e-4 / 2.75e-2 / 9.4e-4 | +12.968 | 114.873 | 54 |
| 2 | 8.3e-4 / 1.9e-4 / 8.88e-3 | +12.456 | 113.473 | 36 |
| 3 | 2.5e-4 / 1.7e-4 / 1.32e-3 | +12.447 | 113.479 | 61 |
| 4 | 5.0e-5 / 1.2e-4 / 1.23e-3 | +12.467 | 113.459 | 65 |
| 5 | 5.0e-5 / 1.1e-4 / 5.6e-4 | +12.778 | 113.125 | 70 |

per-iter shrink factors [0.057, 0.306, 3.335, 1.428, 1.153]; cumulative 3.773x; **monotone = False**;
**W grew**. **VERDICT: FAIL -- the door is checked and closed.**

Cumulative shrinkage exceeds 1.3x but is NOT a leak: the per-iteration factors swing two orders as
the point moves between local modes, which is exactly why the pre-registration demanded monotonicity.
Independently, the number of pulsars concentrating on FALSE fringes climbs 54 -> 70 -- the soft
analogue of the GPS wrong-fix, softened but not cured.

**WHY it cannot compound (the mechanism).** sigma_mc is ALREADY 1e-4..1e-3 dex at iteration 0, far
NARROWER than the 0.0026-0.0044 dex needed. The local mode was never the problem. Median q_max sits
flat at 0.067 -> 0.070 (~16 effective fringes per pulsar, unchanged) and S is flat at +12.5 -> +12.8.
The missing information is not local curvature; it is WHICH of ~16 fringes, and that choice never
sharpens. Local width != global concentration: the referendum already showed the evidence lives on
the plateau.

### VERDICT -- and the frontier statement

**Targeted certification of a CIRCULAR source has its own information wall.** The pulsar term is a
kyr-baseline TIMESTAMP of the source phase; it cannot be read without the CLOCK RATE. The clock rate
is fdot, i.e. Mc. A 22-yr Earth term cannot measure fdot (Delta_phi_E ~ 0.05 rad; information gain
over the prior 1.00-1.73x). The best physically-defensible EM conditioning -- position + period +
host D_L -- still leaves sigma(log10_mc) **8.29x too loose**, the soft cascade cannot bootstrap it
(probe FAIL, monotonicity and wrong-fix both), and no static tier concentrates the posterior
(f(C) = 0.032, failing 0.95 by 16x).

The wall is a FUNCTION with a PRICE, not an impossibility:
  (i)   sigma(log10_h) x 11.3 better  [T x 128 ~ 2840 yr, or ~1 dex louder], **AND**
  (ii)  sigma(log10 D_L) <= 1.0%  [the counterpart is CO-BINDING; strain alone floors at 0.00301 dex]
  (iii) OR an eccentric-harmonic lever kappa >= 8.29, which substitutes for both.

Lever (ii) of the design theorem is therefore not an alternative to lever (iii) -- **it is the
missing ingredient lever (iii) needs.** The E-track's eccentric map is now the PRICED next
experiment: it must deliver kappa >= 8.29.

**R's blind verdict is UNAFFECTED.** The micro-dip and the (f, mc) volume contest live in the dims
R held fixed at truth; R's sky-plane Z_needle had 6/6 positive curvature and quad/Laplace agreement
at 0.2 nat. **f = 6.9e-7 stands**, and including mc as an active dim would only shrink Z_needle
further, so it stands FIRMER.

### The needle is a thin SHELL, not a point (doc note)

The Hessian of lnL_marg at truth over the (f, mc) dims has **5/6 positive curvature eigenvalues and
one negative** (lambda = -1.32e9), stable under Richardson (`|H(h)-H(h/2)|/max|H| = 2.4e-2`) and at
steps 0.3x the sharpest eigen-sigma. It is NOT a saddle: lnL_marg falls 159 nat at 1x base and
919 nat at the box edge along that same eigenvector. It is a **MICRO-DIP at truth**: lnL_marg =
lnL_ref + sum_p m_p with m_p >= 0; lnL_ref is maximal at truth (zero-noise Asimov) but every m_p
GROWS as its pulsar de-registers, so a sub-fringe offset buys more fringe entropy than it costs.
The local max sits ~1e-5 scaled from truth, 0.12 nat higher. CONSEQUENCE: Z_needle is a well-defined
local integral but **NOT a Laplace integral** -- quadrature only, and negative eigenvalues must never
be clipped (clipping turns a non-curving direction into a wide Gaussian and INFLATES Z_needle, i.e.
biases toward the answer that lets B1 proceed).

## CONVENTION -- SMC EVIDENCE GATES (adopted 2026-07-09, after two self-defeating gates)

Two gates written during B1 STEP 2 were **structurally unsatisfiable**, in the same way: the
prescribed remedy moved the gated statistic AWAY from its target. Both are recorded because the
failure mode is not obvious and will recur.

1. **`acc >= 0.25` while the RW scale adapts toward 0.234.** The scale update
   `s <- s*exp(acc - 0.234)` drives acceptance to the RW optimum 0.234, i.e. BELOW the gate.
   Every high-beta stage exhausted `max_mcmc` (observed: sweeps=14, acc=0.24) and the seed loop
   then added seeds -- which cannot change a kernel property. FIX: adapt toward `adapt_acc`
   (0.35), strictly above `target_acc` (0.25), with a startup assert. After the fix the same
   stage ran at sweeps=3, acc=0.29: the strengthened sampler is both correct AND cheaper.
2. **Seed "spread" defined as `max - min`.** A RANGE grows monotonically with the number of
   seeds, so "add seeds until spread <= 0.3 nat" is defeated by its own remedy. Observed on
   Tier C: 0.920 nat (2 seeds) -> 2.028 nat (3 seeds) -> 2.028 (4 seeds).

**CONVENTION GOING FORWARD.**
- The SMC mixing gate is the sample **std** at a FIXED, pre-registered seed count -- never a
  range, never a statistic whose remedy inflates it.
- The **s.e.m.** (std/sqrt(nseed)) is what propagates into the verdict's error bar.
- A gate must never prescribe a remedy that moves the gated statistic the wrong way. If the
  gate cannot be met by the remedy, the gate is measuring the sampler, not the science.
- A FAILED mixing gate and a FIRM verdict are TWO INSTRUMENTS, TWO STATEMENTS. Report both.
  A mixing scatter of s nat cannot overturn an evidence gap of d nat when d >> s; say so with
  the numbers rather than suppressing the verdict or laundering the gate.

## CONVENTION -- LOGICALLY-REDUNDANT MEASUREMENTS RETAIN AUDIT VALUE (adopted 2026-07-09)

A measurement whose OUTCOME is implied by another is not thereby uninformative: it tests the
ASSUMPTIONS the primary measurement conditions on. Tier A and Tier B were argued twice (by me) to
be "logically foreclosed" by f(C) -- their mc boxes are strictly wider, so f(A), f(B) < f(C) must
hold. Run anyway for the table, Tier A's very first seed FALSIFIED the uniform-plateau-density
assumption on which the Tier-C break-even (and hence the whole price paragraph) rested: its box is
24000x larger in volume yet ln Z_box matched C's to within the sampler scatter, i.e. the plateau
had SATURATED. Saturation is invisible from a single box BY CONSTRUCTION -- it can only be seen by
comparing boxes. **Never drop a redundant arm to save time when the primary result is conditioned
on an untested proportionality.**

Corollary, adopted with it: a break-even / threshold is a RESPONSE CURVE, never a single point
extrapolated through an assumed scaling. Measure Z at >= 3 box scales, locate the saturation knee,
and report a BOUND when the crossing lies outside the measured range.

## CONVENTION -- SMC WALL-CLOCK ACCOUNTING (2026-07-09)
SMC seed wall time is ~98% `sweeps x (one G call)`: Tier C s0 37 sweeps -> 1924 s predicted vs
1978 s actual; Tier A s0 56 -> 2912 vs 2976. Delays are MAX-SWEEP GRINDING at high beta, not GPU
contention and not idle-waiting. Diagnose slowness by counting sweeps before blaming the device.

## CONVENTION -- FENCED AGENTS NEVER QUOTE IN-FLIGHT ARTIFACTS (2026-07-11)
A fenced (read-only, no-doc-write) fleet agent reports against BANKED state only -- committed docs
and frozen npz -- never against a run still in progress or an uncommitted working-tree artifact.
An in-flight number can still move (a seed loop can add seeds, an SMC stage can re-mix), so a fenced
quote of it is a snapshot that may be falsified by the time the writer folds it in. The writer (this
session) is the ONLY party that reconciles fleet reports against the live tree and edits the
canonical docs; the fleet writes standalone codenamed reports (SCOUT / SHOVEL / GALLERY / WEAVE /
...) and the writer folds their BANKED numbers in. This is why the three ACCRE census/Q1/siren
reports are consolidated here by the writer rather than self-committed by the agents that produced
them. (Companion conventions this week: SMC evidence gates = sample **std** at a fixed pre-registered
seed count + s.e.m., never a range or a remedy that moves the gated statistic the wrong way, above;
and logically-redundant measurements are RETAINED for audit value -- both in the sections above.)

## STATUS -- TRACK B CLOSED
- Verdict TERMINAL and COMPLETE (deliverable R done 2026-07-08: f=6.9e-7, INFORMATION-THEORETIC
  NO-GO). Machinery BANKED, all gated + reusable by B1 / E-track:
  low-rank split (47x E-step + us scorer), candidate scorer, batched E-step oracle, REPAIRED
  F-stat seeder (sky-only NMS 25 deg, 3/3), repaired-hygiene soft-EM float. NB the soft-EM M-step is
  numerically CHAOTIC (cond ~4e10) -- for any future use, replace Adam with a Newton/trust-region
  step or profile the extrinsics (spec sec 3).
- Queue head = E-TRACK (spec: `CW_transition/WEAVE_etrack_merged_spec.md`, added this commit) +
  the sharpened lever-(ii) lane-measurement objective.
- B1 REFRAMED: noisy ensembles test the CONDITIONAL pipeline (given the source-localised regime),
  NOT cold-start certification.

## CONTAMINATION ASTERISK (which pre-repair results carry it)
The FIRST NO-GO write-up's float-DEPENDENT chain carries a contamination asterisk: L0 (diverged
float), L1 (search-space spec off it), L2 blind search + re-solve (seeded by it). UNAFFECTED and
banked: L2b / L2c (fixed-integer pull-in; geometry + Hessian at truth) and F2 (pure geometry). The
Repair-3 float supersedes the L0/L1/L2 float but reaches the SAME wall.

---

# THE CERTIFICATION CRITERION (criterion-v1, adopted 2026-07-12) — CANONICAL

This section is **THE criterion**. Every certification number in this repo is stated under it.
It supersedes the raw Bayesian bar (`P_true > 0.9`) and the two-layer gate of `FORGE_b1_loop.md`
§9; both are preserved below as superseded-with-trail. Code: `CW_transition/trackB_criterion.py`
(fits the floor and emits the table from the banks; asserts the invariants; no new realisations).

## The three layers

For pulsar `a`, with `dlnL_a` = the likelihood-only fringe gap (best minus runner-up peak,
prior-free), `K_counted,a` = the counted candidate fringes, `q_max,a` = the E-step fringe
posterior's modal mass:

    DETECTION      dlnL_a > max( ln K_counted,a , DLNL_FLOOR )      DLNL_FLOOR = 9.01 nat
    CERTIFICATION  q_max,a > 0.9   (strict: > 0.99)   applied ONLY within detections

- **Layer 1, `ln K` — the trials factor.** The fringe-breaking evidence must beat the number of
  fringes it chose among. Relative, per-pulsar. This is GEO's `ident` / "flat count".
- **Layer 2, `DLNL_FLOOR` — the absolute floor.** NEW in criterion-v1. Layer 1 alone is a
  likelihood-ratio test whose bar collapses for a tightly-EM-prior'd pulsar: J0437-4715 has
  `ln K = 1.39 nat`, the array minimum, so a pure-noise fluctuation clears it. The floor is the
  absolute evidence a detection must carry regardless of how few fringes it had to beat.
- **Layer 3, `q_max` — the Bayesian bar.** Retained, but demoted: it is now a *quality* bar on
  something already established as a detection, not itself a detection statistic.

**Certification is defined on `q_max`, NOT on `P_true`** (binding invariant; `trackB_criterion.py`
asserts it). `P_true > 0.9` does not reproduce the banked certifications: the cells where
`q_max > 0.9` but `P_true <= 0.9` are exactly the **wrong-certifications** (2 in Arm A, 8 in Arm B).
Within detections `q_max == P_true` — i.e. once a cell passes the detection gate, the modal fringe
*is* the true fringe. Scoring on `P_true` would silently define wrong-certs out of existence.

## Derivation of DLNL_FLOOR = 9.01 nat

Fitted against the **89 banked null realisations** (`reports/flat_*.npz`), 27 of which are nulls:

| null bank | n | construction |
|---|---|---|
| `null` (910000) | 12 | 3 loud sources scrambled (the original B1.2 line) |
| `nullL` | 5 | 3 loud scrambled, 13 faint at truth |
| `nullA` | 5 | ALL 16 sources scrambled (honest Cornish-Sampson) |
| `nullN` | 5 | **NO CW in the data**, recovery at the true source |

At `floor = 0` (i.e. the FORGE §9 two-layer gate) **9 null cells still certify**. The criterion is
the *smallest floor that zeroes all of them*. Because detection uses a strict `>`, that floor is
exactly the largest `dlnL` among those 9 cells:

| kind | pulsar | dlnL | ln K | q_max |
|---|---|---|---|---|
| **nullN** | **J1713+0747** | **9.009** | 6.752 | 1.000 | ← the binding cell: pure noise, no CW |
| nullL | J1713+0747 | 7.211 | 6.752 | 0.999 |
| null | J1713+0747 | 6.778 | 6.752 | 0.999 |
| null | J1909-3744 | 5.831 | 4.190 | 1.000 |
| nullA | J1909-3744 | 4.679 | 4.220 | 0.974 |
| null | J0437-4715 | 3.599 | **1.386** | 0.998 |
| nullL | J0437-4715 | 3.304 | 1.792 | 0.984 |
| null | J0437-4715 | 1.960 | 1.792 | 0.997 |
| nullA | J0437-4715 | 1.679 | **1.386** | 0.994 |

`floor_min = 9.0094` nat → **adopted `DLNL_FLOOR = 9.01`** (rounded UP, so the zeroing property is
preserved). FORGE §9.3 predicted "~8 nat"; the measured value is **9.01**, and the cell that sets it
is a **`nullN` J1713 fluctuation on data containing no CW at all** — the floor is set by pure noise,
not by a mis-modelled source.

**SMALL-K ANATOMY — verified.** J0437-4715's residual false alarms (the §9.3 pathology: `ln K = 1.39`
beaten by noise) sit at `dlnL = 3.60 / 3.30 / 1.96 / 1.68`. **All four die at any floor ≥ 3.60 nat**,
i.e. with **5.4 nat to spare** below the adopted floor. The small-K false-alarm channel is closed
decisively, not marginally.

## THE MARGIN — and it is thin

| quantity | value |
|---|---|
| null-side margin | **0 by construction** — the floor IS the largest null `dlnL` |
| signal-side margin | **0.29 nat** — lowest surviving real detection (GEO, J1909, `dlnL = 9.30`) |
| Arm-A margin | 0.49 nat (lowest survivor `dlnL = 9.50`, J1909) |
| **Arm-B margin** | **NEGATIVE** — Arm B's largest `dlnL` is **8.0**, *below* the null max 9.01 |

**State this margin whenever the floor is quoted.** The floor is the **maximum of a 27-realisation
null sample** — the noisiest order statistic there is. A 28th null realisation could exceed it. The
criterion is honest about what it is: a floor calibrated to the nulls we have, with a 0.29-nat gap to
the nearest thing it is trying to keep. It is not a comfortable separation, and no amount of restating
makes it one.

The consequence, stated plainly: **noisy Arm B contains no cell whose likelihood gap exceeds the worst
pure-noise false alarm in the null banks.** That is not a thresholding artefact — it is the finding.

## REGISTRATION-TOLERANCE SENSITIVITY — UNMEASURED, and it is a real gap

The check flagged before adoption: the null banks pass through the seeded pipeline, and this spec
(§ "certification tolerance", L286) puts the source-parameter certification tolerance at **~1e-4
scaled**. Does the floor move with the registration/seeder tolerance?

**It cannot be answered from the banks.** Every one of the 27 null realisations carries
`tol_scale = 0.0` — a single point, no grid. The floor is therefore calibrated **at perfect
registration only**.

The one available off-tolerance datum is `wrongpos` (`tol_scale = 5.0`, i.e. 5× the certification
tolerance, n=2, B1.4): its sole certification is **J0437-4715, `dlnL = 4.41`, on the TRUE fringe** —
a *correct* certification, and the 9.01-nat floor **kills it**. So at 5× tolerance the floor is
already destroying true positives, and we have no null bank at that tolerance to say whether it is
also still suppressing false ones.

**Open, and it gates any above-onset claim: bank nulls on a `tol_scale` grid and re-fit.** Until then
the floor is stated as calibrated at `tol_scale = 0` and the tolerance-dependence is an acknowledged
unknown, not a solved one. IGNITE should carry this.

## Superseded, with trail

| criterion | GEO/draw | Arm A/real | Arm B/real | null/real | status |
|---|---|---|---|---|---|
| Bayesian `P_true > 0.9` | 4.50 ± 1.48 | 2.87 ± 1.48 | 1.43 ± 1.05 | **0.8 – 2.8** | **SUPERSEDED** — no detection statistic; carries a source-independent floor |
| two-layer `dlnL > ln K` ⊕ `q>0.9` (FORGE §9) | 1.35 ± 0.82 | 0.33 ± 0.54 | 0.13 ± 0.43 | **0.2 – 0.4** | **SUPERSEDED** — null still fires, from small-K |
| **three-layer (criterion-v1)** | **0.275** | **0.067** | **0.000** | **0.000** | **CANONICAL** |

Neither superseded number was wrong as *measured*; both were wrong as *interpreted*. The Bayesian
count was never a detection statistic, and the two-layer gate's `ln K` bar is not a floor for a
pulsar whose `K` is 4.

## CONVENTION — confidence without a detection statistic is prior-pinning in disguise

Adopted 2026-07-12, from the two-layer lesson.

`P_true > 0.9` / `q_max > 0.9` is a statement about **where posterior mass sits among candidate
fringes**, not about **whether there is anything there to find**. A pulsar with a tight EM prior
concentrates >0.9 of its mass on the MAP fringe **from the prior tail alone** — measured in `nullN`:
pure noise, no CW in the data, and it still certifies. The Bayesian bar cannot tell that apart from a
detection, because it never asked the question.

**Rule: every confidence bar must sit downstream of a detection statistic that can return zero.**
A criterion that cannot fire on a null is not a criterion. Before quoting any posterior-mass
threshold, state what had to be *detected* first, and show the null.

Corollary, and the reason this cost a re-score: the same prior-pinning that makes J0437-4715 the
*robust* certifier (survives sky redraws, noise, off-prior-mean truth, and a wrong counterpart) is
what makes it fire on pure noise. **Robustness to source error and vulnerability to noise are the
same property viewed from two sides.** Tiny `K` cuts both ways, always.

## CONVENTION — summary files carry raw statistics, not only verdicts

Adopted 2026-07-12. FORGE's B1.0 production npz stored booleans (`cert90`, `wrong90`) but not the
`dlnL` they were derived from. When §3's null anatomy demanded a re-cut on a *different* criterion,
the verdicts were useless and the run had to be re-executed on the cluster purely to extract an
array that already existed in memory (job 12496241, 420 s). **Bank the statistic, not the verdict:**
a verdict answers one question, a statistic answers the ones you have not thought of yet. Any
summary npz must carry the continuous quantities (`dlnL`, `lnK`, `q_max`, `P_true`, `on_true`) that
a re-scoring could need.

## CONVENTION — reports and their empirical basis travel together

Adopted 2026-07-12. `FORGE_b1_loop.md` §9 was written against `flat_*.npz`; the report reached
`reports/` and the banks did not, and the criterion could not be fitted until they followed.
**A report is not landed until the arrays it is scored from are landed with it.**

---

# THE CERTIFICATION CRITERION (criterion-v2, adopted 2026-07-12) — CANONICAL

**This section supersedes criterion-v1 above**, which is preserved in full as
superseded-with-trail. criterion-v1 was fitted to 27 nulls at one loudness, one baseline, one
tolerance; IGNITE (`reports/IGNITE_onset_map.md`, banks `reports/ignite_bank.npz`) swept the box
around it and found the floor is not a constant of the pipeline. Three decisions follow. They are
**executed here with rationale, not reopened.** Gates: `CW_transition/criterion_v2_gates.py`
(8/8 PASS + census triple bit-identical; banked npz only, no GPU, no new realisations).

    DETECTION      dlnL_a > max( ln K_counted,a , floor(h, T, tol) )
    CERTIFICATION  q_max,a > 0.9  (strict 0.99)   applied ONLY within detections
    PURITY         co-registration statistic R_a  — PRE-REGISTERED TEST, NOT YET ADOPTED (D3)

Layers 1 and 3 are unchanged from criterion-v1. Layer 2 becomes a **function** (D2), calibrated
against a **named null family** (D1). Layer 4 is **defined and pre-registered but not in force**
(D3): no number in this repo may be quoted under it until IGNITE-2 reports.

## D1 — THE NULL FAMILY: counterpart-matched nulls are the operative calibration

**DECISION. The operative calibration family for the targeted programme is `fN` —
counterpart-matched nulls: pure noise / no CW in the data, recovery at the TRUE source position.**
The all-null family `fALL` (adding wrong-counterpart scrambles: `nullA` all-16-scrambled, `nullL`
loud-scrambled) is retained permanently as the **blind-robust column** and travels beside every
`fN` number in every onset table, in the docs and in the npz.

**Rationale, recorded.** A *targeted* analysis faces exactly the counterpart-matched null. The
premise of the conditional pipeline is that a real counterpart exists and its sky position is
known — that is the scenario, by construction. The false alarm it can actually suffer is **noise
mimicking fringe-breaking under the correct source model**. A sky-scramble null asks whether the
pipeline can be fooled by a source that is *not there* — a **blind-search** question, and the
targeted analysis does not ask it. Calibrating a targeted criterion against a blind-search null
imports a bar the targeted scenario never has to clear.

**The consequence of the alternative, recorded permanently and never suppressed:** under `fALL`,
**there is no onset anywhere in the modelled grid.** Best cell = 0.24 certifications/realisation,
of which **0.22 correct** — against the >1 bar, at *every* one of the 24 cells (h ∈ {−13.25,
−13.0, −12.75, −12.5} × T ∈ {15, 20, 30} yr × {lit, vlbi}). The scrambled-source null's
noise-lock grows ∝ h², so the wrong-counterpart-robust floor rises **faster than the signal** and
closes the very window it was meant to guard. **This is not a footnote. It is the price of D1,
and it is a physics decision, not a technicality:** the targeted programme's onset exists *because
of* the null family it is calibrated against. Any onset number quoted outside this repo carries
its `fALL` column or it is not quoted.

**What D1 gives up, and where it is bought back.** `fN` presumes the counterpart is right. It
therefore has **no defence against a wrong counterpart** — and IGNITE measured exactly that hole:
the Stage-2 scrambled-source loop **detects in 2 of 5 realisations** under the `fN` floor
(`dlnL` up to ~15 nat > `fN` = 5.46). Defending that hole was the `fALL` floor's job, and D1 fires
the `fALL` floor. **D3 is the replacement defence** — and D3 buys it back *without* paying the
window-closing price, because co-registration rejects a wrong counterpart on **geometry**
(the pulsars' implied source solutions do not agree) rather than on **amplitude** (a floor tall
enough to outrun the noise-lock). D1 and D3 are one design, split across two decisions; D1 is
adopted now, D3 is on test. **Until D3 reports, the wrong-counterpart hole is OPEN and stated.**

## D2 — THE FLOOR IS A FUNCTION, AND ITS ESTIMATOR CHANGES

**DECISION. The constant `DLNL_FLOOR = 9.01 nat` is RETIRED.** It was never a property of the
pipeline; it was the census-loudness value of a function. The floor is now
`floor(h, T, tol)`, **refit per cell, never inherited**.

### D2.1 The loudness law and its mechanism

    floor_fN   ∝ h^1.66      (measured across the grid; per-(T,tier) fits span 1.5–1.7)
    floor_fALL ∝ h^1.88      (per-(T,tier) fits span 1.7–2.0)

**Mechanism, recorded.** The E-step evaluates a model whose pulsar-term amplitude ∝ h against the
data. The per-fringe log-likelihood carries a **matched-filter cross term** that is *linear in the
model amplitude* — so on data containing **no CW at all**, the null `dlnL` fluctuations still grow
with h. With a *scrambled* source meeting loud real data the noise-lock grows ∝ h², which is why
`fALL` scales more steeply than `fN`. Making the source louder raises the bar almost as fast as it
raises the signal, and **the certified count is non-monotone in h** (T = 20 vlbi: 0.72/real at
h = −13.25 falls to 0.38/real at h = −12.5 — a 10× louder source *lowers* the honest count).
Onset is therefore **baseline-driven, not loudness-driven**: T^{5/2} `fdot`/coherence leverage
beats the h^{1.66} floor race; louder alone does not.

Concretely, the h-law lives in the **Gumbel scale** of the null offender distribution:
`beta` = 2.1–2.4 nat at h = −13.0, 4.2–7.0 at h = −12.75, 13–24 at h = −12.5. The floor is
`mu + beta·z`, and both `mu` and `beta` carry h^1.66. **`9.01` is the census-loudness value of
`beta·z + mu`, nothing more.**

### D2.2 The estimator: fixed-FPR tail fit, NOT the sample maximum

**This is the part of D2 that the data forced, and it is not cosmetic.** criterion-v1's floor is
the **maximum of N nulls** — the smallest value zeroing all null certifications. That estimator
does not converge:

- The null offender statistic is itself a **max over pulsars**, so it sits in the **Gumbel domain**
  by construction, not by assumption. There, `sd(max_N) → 1.283·beta`, **INDEPENDENT of N**, while
  `E[max_N] = mu + beta·ln N` **creeps up without bound**.
- Measured (G7): `sd(max_N)` = 8.91 / 8.68 / 8.79 / 8.74 nat at N = 10 / 30 / 100 / 1000. **Flat.**
- **Therefore banking more nulls does not stabilise the criterion-v1 floor — it inflates it.**
  IGNITE §7.5's "more nulls per cell is the single cheapest credibility purchase" is **half right,
  and this is the correction**: more nulls buy nothing *with the max estimator*. The estimator had
  to change first.
- Worse, the max-of-N floor has **no fixed false-alarm rate**: it is implicitly the
  `1 − 1/(N+1)` quantile, so **its stringency is an accident of how many nulls happened to be
  banked** (27 → 30 → 750 are three different criteria wearing one name).

**ADOPTED ESTIMATOR.** `floor(h, T, tol)` is the **(1 − α) quantile of the per-cell null offender
distribution at a STATED per-realisation false-alarm rate α**, estimated by a **Gumbel
(block-maximum) tail fit** over the cell's nulls:

    floor(h, T, tol) = mu_hat + beta_hat · z(α),     z(α) = -ln(-ln(1-α)),     α = 0.05

`α = 0.05` is adopted; criterion-v1's max-of-27 was implicitly α ≈ 1/28 = 0.036, so this is a
mild, **explicit, N-independent** loosening of a bar that was previously set by accident. α is now
a stated dial and must be quoted with the floor.

**Sizing N (the answer to IGNITE's ±5-nat-per-30 scatter).** With the fitted estimator,
`sd(floor_hat) = c · beta / sqrt(N)` with **c = 2.80** measured (α = 0.05) — it *does* shrink:

| N | sd(floor_hat), onset-cell scale (beta ≈ 6.9) |
|---|---|
| 30 | 3.54 nat |
| 100 | 1.94 nat |
| 400 | **0.96 nat** |

Because `beta` itself scales as h^1.66, an **absolute** 1-nat target is loudness-dependent and
over-specifies the loud cells (a 1-nat target against a 45-nat floor is meaningless). Two rules,
both binding:

- **SCALE-FREE RULE (operative): `N ≥ 100` counterpart-matched nulls per cell**, giving
  `sd(floor_hat) < 10 %` of the floor **at any loudness** (G7: N ≥ 89 suffices; 100 adopted).
- **ABSOLUTE RULE (onset cells): `N ≥ 150`**, giving `sd(floor_hat) < 1 nat` at the onset-cell
  scale (beta ≈ 4.2 at (−12.75, 30, lit)). At the loud h = −12.5 cells a 1-nat absolute target
  would demand N ≈ 2 000–5 000 and is **explicitly not adopted** — the scale-free rule governs there.

Two-pass procedure: bank 30 nulls, estimate `beta_hat`, then size N per cell. Cost at IGNITE's
measured throughput (~1–2 s/realisation warm at T = 30) is ~2 GPU-hours for 150 nulls × 24 cells —
**the sizing is cheap, and the reason it was never done is that nobody had checked that the
estimator converges.**

### D2.3 The tolerance axis

`tol` = registration offset in units of the spec-L286 certification tolerance (1e-4 scaled).
Measured at the fiducial cell (IGNITE §3, `tol` ∈ {0, 1, 2, 5}):

- **The null floor is FLAT-to-mildly-RISING in tol, and small** — `fN` = 0.00 → 0.00 → 2.06 →
  4.37 nat across tol = 0 → 1 → 2 → 5. The `fALL` spread ({8.48, 14.03, 8.09, 4.37}) is
  **sampling noise, not tol dependence** — four independent 30-null redraws of one statistic.
- **It is the TRUE-POSITIVE channel that dies of mis-registration**: true two-layer
  certifications/realisation fall 0.14 → 0.05 → 0.10 → **0.00** by tol = 5. And **no per-tol refit
  floor kills a single surviving true positive** (0 own-floor kills at every tol).

**The K1 hole is closed, and it inverts.** criterion-v1's tolerance caveat feared an *inflating
null*. The measured failure is the opposite: the null barely moves, the signal evaporates. The
banked "9.01-nat floor kills the correct `wrongpos` J0437 certification (`dlnL` = 4.41)" pathology
was an **artifact of applying a tol = 0 floor to a tol = 5 realisation** — calibrated at its own
tolerance, the floor (4.37 nat) sits *below* that survivor. **Refit the floor at the analysis's own
tolerance and the pathology disappears.**

### D2.4 criterion-v1's 0.29-nat margin — superseded, with trail

criterion-v1 quoted a **0.29-nat signal-side margin** (floor 9.01 vs lowest surviving real
detection, GEO J1909 `dlnL` = 9.30) and correctly flagged it as thin. **It is now annotated as
WITHIN CALIBRATION NOISE and carries no evidential weight.** The floor it is measured against is a
max-of-27 order statistic whose sampling scatter is **±5 nat** (IGNITE §3, four independent 30-null
redraws: {8.48, 14.03, 8.09, 4.37}); a 0.29-nat gap against a ±5-nat ruler is not a margin, it is a
rounding error. The same annotation applies to **every** IGNITE onset-cell margin (0.01–2.0 nat).

**This does not overturn a single criterion-v1 verdict.** Arm B's largest `dlnL` is 8.0 against a
floor whose scatter band is ~4–14 nat: Arm B detects **0.000** under any floor in that band, and
the headline — *the noisy conditional pipeline, honestly gated, detects nothing at census
loudness* — is **robust to the calibration noise, because it never depended on the margin**. What
dies is the *precision* of the margin, not the *sign* of the result. **Nothing measured under
criterion-v1 is retracted; the criterion-v1 table stands and is gated bit-identically (G1–G2).**

## D3 — THE PURITY LAYER: pre-registered TEST, adoption CONDITIONAL

**criterion-v1's purity property is dead above onset.** It claimed "the gate perfectly purifies
what is left" (0 wrong-certs, both arms). IGNITE shows that was a **census-loudness artifact**: at
the (−12.75, 30 yr, lit) onset cell, **23 of 50 realisations carry a wrong certification** — the
same noise-lock that raises the floor gives *wrong* fringes floor-beating gaps. **Fringe
correctness — the one discriminator that survived real noise at census loudness — degrades exactly
where the count turns on.** A count without a purity statement is now meaningless.

### D3.1 The statistic — the needle's co-registration, made a statistic

The needle: **the true source solution is the unique point at which every pulsar's fringe comb
co-registers.** A pulsar on the wrong fringe is not merely "wrong" — it demands a *different
source* from the one its neighbours demand. Make that a statistic.

Each pulsar `a` has a pulsar-term phase `Phi_a(theta)` and a **lever**
`g_a = grad_theta Phi_a` in the source parameters `theta = (log10 f, log10 mc)`
(SIREN's levers, analytic in the lag `tau_a`: `dPhi/dlog10 f ∝ tau_a`,
`dPhi/dlog10 mc ∝ fdot·tau_a²`). Choosing fringe `k_a` instead of the true `k_a*` is a phase error
of exactly `2π(k_a − k_a*)`, which can only be absorbed by displacing the source. So `a`'s MAP
fringe **implies a source displacement**, minimum-norm:

    u_a = 2π · (k_a − k_a*) · g_a / |g_a|²          [the source shift pulsar `a` demands]

**Cheapest sufficient form — a 1-dof leave-one-out chi-square.** For candidate certification `a`,
build the inverse-variance-weighted source solution `u_R` implied by a **reference set R** of the
*other* pulsars, and test concordance:

    R_a  =  (u_a − u_R)^T (Sigma_a + Sigma_R)^-1 (u_a − u_R)       ~ chi2(2) under co-registration
    VETO: certify `a` only if  R_a < chi2_crit(2, p = 0.01)

**Computable from the bank with no new likelihood evaluations.** The common (unknown) true-source
displacement **cancels in the difference** `u_a − u_R`, so only the banked *relative* fringe
offsets are needed: `ignite_bank.npz` carries `mapk` and `n_true_grid` per pulsar per realisation,
and `g_a`, `Sigma_a` are analytic in `tau_a` (or come free from the amortised Fisher the E-step
already computes). **No re-run of the E-step, no GPU.** That is the cheapest sufficient form.

**Why it should have teeth here.** The wrong fringes IGNITE certifies are not `±1` neighbours —
measured `|Δk|` at the onset cell runs **25 to 395**. A wrong fringe therefore demands a source
displacement of `2π·Δk / |g_a|`: orders of magnitude away from what the correct pulsars demand.
**Watch the lever-dependence** (the mechanism to check, pre-registered as such): a long-lag pulsar
has a huge `|g_a|` and absorbs a given phase error with a *small* source shift, so its wrong
fringes are the *cheapest* to hide. Power is not uniform across the array, and the test must report
it per-pulsar, not pooled.

### D3.2 The reference set — and the ceiling it imposes (measured, pre-registered)

Two variants, both to be tested:

- **`R_det`** — reference set = the other **DETECTED** pulsars (layers 1+2) in the same
  realisation. This is the literal form.
- **`R_all`** — reference set = **all 115 other pulsars' fringe posteriors**, information-weighted.
  A pulsar that fails detection still has a `q(k)` posterior and a lever; it constrains the source
  weakly but non-trivially.

**`R_det` has a hard, measured ceiling, and it is recorded before the test is run** (G8, from the
bank, at the (−12.75, 30, lit) onset cell):

| | wrong certs | true certs |
|---|---|---|
| total | 23 | 77 |
| statistic **DEFINED** (≥1 other detection) | **20 (87 %)** | **69 (90 %)** |
| **UNDEFINED** (candidate is the only detection) | 3 | 8 |

**The statistic is undefined for a singleton detection — there is nothing to co-register against.**
Consequences, both binding on the pre-registration:

1. Undefined **must default to PASS**. Defaulting to FAIL would kill 8 of 77 true certs at this
   cell (and 15 of 58 at the vlbi onset cell, and 22 of 47 at h = −12.5), **violating criterion (b)
   outright before the statistic does any work.**
2. Therefore **`R_det` can kill at most 20 of the 23 wrong certs = 87 %**, and **cannot reach
   100 %**, at *any* threshold. At the loud h = −12.5 cell the ceiling is **43 %**. This is a
   property of the *detected-set restriction*, not of the statistic.

`R_all` has **no such ceiling** — its reference set is never empty — and it is the variant expected
to carry the test. **`R_det` is retained as the control precisely because its ceiling is known in
advance: if `R_det` lands at ~87 % and `R_all` does not beat it, the co-registration idea is not
what is doing the work.**

### D3.3 THE PRE-REGISTRATION (binding; adoption is conditional on it)

Applied to IGNITE's **banked Stage-1 cells** (`reports/ignite_bank.npz`), with the threshold fixed
at `p = 0.01` **before looking**, and reported for `R_det` and `R_all` **separately and both**:

- **(a) KILL THE WRONG CERTS.** At the (−12.75, 30 yr, lit) onset cell the veto must remove
  **≥ 90 % of the 23 wrong certifications.** *Recorded in advance: `R_det` cannot pass this — its
  ceiling is 87 % (20/23). Passing (a) therefore REQUIRES `R_all`.* This is not a moved goalpost;
  it is the ceiling being stated before the run, per convention.
- **(b) PRESERVE THE TRUE CERTS.** **≥ 90 % of true certifications survive** the veto at both
  onset cells ((−12.75, 30, lit) and (−13.25, 30, vlbi)). Undefined ⇒ PASS.
- **(c) FIRE ZERO ON THE NULLS.** **Zero surviving certifications on all three null banks**
  (`nullN`, `nullA`, `nullL`) at every cell. *Note, and it is the substantive part:* under `fN`
  this is already true by floor construction, so (c) alone is not evidence. The **binding** report
  is the veto's **rejection rate on the `nullA`/`nullL` DETECTIONS** (1–4 per realisation at the
  onset cells) — that number, and only that number, measures whether the purity layer buys back the
  wrong-counterpart robustness **D1 gave up**. It must be **high**, and it must be quoted.
- **(d) KILL THE SCRAMBLED LOOP.** The veto must reject the Stage-2 **scrambled-source loop's
  2-of-5 detections** (`ig_loop*` banks). This is the D1 hole in its most concrete form; if the
  purity layer cannot close it, D1's stated hole stays open.

**ADOPTION IS CONDITIONAL ON (a)–(d). If they are not all met, the purity layer is NOT adopted,
the wrong-certification rate above onset stands as measured (23/50), and every above-onset count in
this repo continues to carry its purity number beside it.** No partial adoption, no
threshold-tuning after the fact: the `p = 0.01` cut and the reference-set variants are fixed here,
in advance. **This is a test, not a plan to adopt.**

## Superseded, with trail (criterion-v1 → criterion-v2)

| element | criterion-v1 | criterion-v2 | status |
|---|---|---|---|
| Layer 2 floor | constant **9.01 nat** | **`floor(h, T, tol)`**, refit per cell | **SUPERSEDED** — 9.01 was the census-loudness value of an h^1.66 function |
| floor estimator | **max of N nulls** | **Gumbel-tail (1−α) quantile, α = 0.05** | **SUPERSEDED** — max-of-N has sd flat in N and no fixed FPR |
| nulls per cell | 27 (one cell) | **N ≥ 100** (scale-free); **≥ 150** at onset cells | **SUPERSEDED** |
| null family | `fALL` (implicit, unnamed) | **`fN` operative + `fALL` blind-robust column** | **DECIDED (D1)** — was never a stated choice |
| purity | "the gate perfectly purifies what is left" | **FALSE above onset** (23/50 wrong); purity layer on test | **REFUTED (D3)** — census-loudness artifact |
| 0.29-nat margin | quoted as thin-but-real | **within ±5-nat calibration noise** | **SUPERSEDED-WITH-TRAIL** — carries no evidential weight; no verdict changes |

**Neither criterion-v1 number was wrong as measured; the criterion-v1 table is gated
bit-identically here (G1–G3) and stands.** What was wrong was the *interpretation*: a floor fitted
at one loudness was read as a constant of the pipeline, and a purity property measured below onset
was read as a property of the gate.

## CONVENTION — a calibrated threshold must state its false-alarm rate and its sampling scatter

Adopted 2026-07-12, from D2. criterion-v1's floor was the max of the nulls that happened to be
banked. That is not a threshold, it is an order statistic wearing a threshold's clothes: its
stringency drifts with sample size (`E[max] = mu + beta·ln N`) and its scatter never shrinks
(`sd(max) = 1.283·beta`, independent of N). **Any threshold this project adopts must state
(i) its target false-alarm rate α, (ii) the estimator used to reach it, and (iii) the sampling
scatter of that estimator at the N actually banked.** A threshold whose sampling scatter exceeds
the margin it is defending is not defending anything — and until it is measured, you cannot know
which case you are in. **Corollary: "bank more nulls" is not automatically a credibility purchase.
It is one only if the estimator converges.**

---

# CRITERION-v2.1 (adopted 2026-07-12, tag `criterion-v2.1`) — CANONICAL

**The operative criterion is UNCHANGED from criterion-v2. What v2.1 adds is two rejections, each by
its own pre-registration, and the conventions they forced.**

    DETECTION      dlnL_a > max( ln K_counted,a , floor(h, T, tol) )     [D1 family, D2 estimator]
    CERTIFICATION  q_max,a > 0.9  (strict 0.99)   applied ONLY within detections
    PURITY         NONE. Both candidate layers were tested and BOTH ARE REJECTED.

## D3 — PURITY LAYER (per-pulsar): **TESTED, REJECTED** (IGNITE-2, `reports/IGNITE2_softloop.md` §1)

Scorecard against the pre-registration of D3.3: **(a) PASS** (23/23 = 100 % of wrong certs killed;
the `R_det` control landing at exactly its pre-recorded 87 % ceiling, which shows the co-registration
idea is what does the work); **(b) FAIL** — 2/77 = **3 %** (lit) and 39/58 = **67 %** (vlbi) of TRUE
certifications survive, against a **≥90 %** bar; **(c)** 42/42 = 100 % rejection of `nullA`/`nullL`
detections; **(d) FAIL** — 5 of 6 scrambled-loop certifications killed, one survives (J1909-3744,
Δk = −4, `R_all` = 4.65 against the 9.21 bar). **NOT ADOPTED. No partial adoption; no threshold was
tuned.**

**The anatomy (recorded).** Above onset the array-wide fringe field is itself poisoned: the
leave-one-out reference `u_R` is dragged by the same confident wrong fringes everywhere in the array,
so `u_R` is displaced from zero in essentially every realisation, and a TRUE certification
(`u_a` = 0) then fails concordance **with its own poisoned reference**. The veto measures *"this
realisation's fringe field is discordant"* — true above onset for every realisation — rather than
*"this pulsar disagrees with the others"*, which is what certification purity needs. **A per-pulsar
veto built on a reference this contaminated cannot be selective at any threshold: the failure is
structural, not a tuning artifact.** The layer fails hardest exactly where it was needed most (3 %
survival at the wide-prior lit cell vs 67 % at the tight-prior vlbi cell).

## D4 — THE REALISATION-LEVEL DISCORDANCE GATE: **DESIGNED, PRE-REGISTERED, TESTED, REJECTED**

Full report `reports/D4_discordance_gate.md`; code `CW_transition/criterion_v2_1_d4.py` (which
carries the pre-registration text machine-readably) + `CW_transition/run_d4_score.py`; scores banked
to `reports/d4_score.npz`. CPU-only, banked npz only, no new realisations.

D3 left exactly one live lead — its **(c) = 42/42**: the co-registration statistic rejects
wrong-counterpart *detections* perfectly at the **realisation** level even where it destroys true
certifications per pulsar. D4 promotes that lead to a gate.

**THE STATISTIC** (chosen on the banked distributions *before* any condition was scored — the
cheapest sufficient aggregate, one 2×2 solve per realisation, no leave-one-out loop):

    S_det = J^T I^-1 J  over the DETECTED set D,   I = sum_D g_b g_b^T / s2_b,  J = sum_D g_b (2 pi dk_b)/s2_b
          = chi2(u=0) - min_u chi2(u)   -- the GLRT for "the pulsars the data actually registered
            co-register at a source OTHER than the assumed counterpart".  FLAG => veto the realisation.

Selected because it is **the only aggregate whose true-signal distribution concentrates at the null
value** (median **0.0** at both onset cells): the `max/mean/min/frac` of the per-candidate `R_all`
inherit D3's poisoned reference in full (true-signal median `min_R` = 1.4e4 against a 9.21 bar), and
`S_ref` over all 116 pulsars puts **pure noise ABOVE the wrong-counterpart population**. **This
inverts D3's variant ranking** — there `R_all` carried the test and `R_det` was the control — because
the detected set is the *clean* subset, and building the reference from it is what removes the
poisoning. It is the formalisation of IGNITE-2 §1.4's sentence: *"the R_det control kills 6/6 — but
for the degenerate reason that every detection is discordant with every other under a scrambled
source, which is exactly the realisation-level (not per-pulsar) signal."*

**VERDICT: (i) FAILS IN ALL EIGHT PRE-REGISTERED COMBINATIONS** (2 dk-conventions × 2 thresholds ×
2 onset cells). Best catch at a ≤10 % false-flag rate: **90.3 %** against the **≥95 %** bar. The one
setting that catches 97.5 % (ORACLE / χ²-bar, lit cell) flags **44 %** of true-signal realisations
against the **≤10 %** bar. **Adoption required (i) AND (ii). NOT ADOPTED.**

**THE ANATOMY — both failures are ONE statement: `S_det` is a `|Δk|` detector, and `|Δk|` is not the
difference between a right and a wrong counterpart.**

- **The misses.** Every wrong-counterpart realisation the gate misses is a noise-lock that landed
  **within ±1 fringe of truth** (median max|Δk| among detections = **1**, against **137** (lit) /
  **13** (vlbi) for the caught ones). The limit case is decisive: one missed realisation has
  **Δk = 0** — a wrong counterpart whose surviving detection sits on the **true** fringe. The fringes
  co-register *because they are right*; the **source** is wrong. **A co-registration statistic tests
  the fringes, not the counterpart. NO co-registration statistic can close the D1 hole in general.**
- **The false flags.** At the (−12.75, 30, lit) onset cell **13 of 36 (36 %)** detecting TRUE-signal
  realisations have an **impure detected set** (≥1 detection on a wrong fringe). The gate's
  in-sample false-flag rate there is **36.1 %**. **These are the same number** — the gate faithfully
  measures the cell's own impurity and cannot beat it. At the vlbi cell, impurity 12 % → false-flag
  rate **0 %**.

> **THE SCISSORS. D3 failed because the REFERENCE was poisoned; D4 fails because the POPULATION IT
> MUST PROTECT is itself poisoned. Same disease, one level up.** Above onset a true-signal
> realisation and a wrong-counterpart realisation *contain the same kind of object* — a confident
> noise-locked fringe — and a geometry test cannot tell which of them is the counterpart.

**(iii) THE D1 HOLE'S CLOSURE TEST — the one genuine positive, and it is reported in full.** All
**three** scrambled-loop keepers are flagged by the realisation-level form, **including the
small-|Δk| J0437-4715 (Δk = −4) case that defeated the per-pulsar statistic** (`R_all` = 4.65 →
MISSED; `S_det` = 55.9 → FLAG); B1937+21 (Δk = +21) → `S_det` = 1 728; J0711-6830 (Δk = +231) →
3.2e5. **The hole is closable on every instance this campaign holds — and no gate that closes it
survives condition (ii).** That is the hole's status, stated exactly: not "no statistic sees these
events", but "the statistic that sees them cannot distinguish them from the impurity the true
population already carries at the only cells where the count turns on."

**STATUS OF THE D1 WRONG-COUNTERPART HOLE: OPEN, PERMANENTLY STATED, and now known to be structurally
un-closable by co-registration.** The wrong-certification rate travels beside every above-onset count:
**14/50 realisations at the lit onset cell under the fresh floors** (23/50 under the retired
max-of-10 floors).

## CONVENTION — a statistic evaluated against truth is an ORACLE until its implementable form is scored

Adopted 2026-07-12, from D4. The fringe grid is indexed about the **EM-prior mean** (`k = 0`), so
`dk = mapk − n_true_grid` is referenced to the **TRUE** fringe — which a real analysis does not know.
D4 therefore scored **both** forms: **ORACLE** (`dk = mapk − n_true_grid`, D3's convention) and
**IMPLEMENTABLE** (`dk = mapk`, referenced to the prior mean, with the `(1−q_max)` factor on the
prior-variance term **dropped — forced by the change of reference, not tuned**, because the prior's
distance error is present however confident the fringe posterior is).

**Measured: the implementable form is 2–4× weaker** (catch 25–52 % vs 43–97.5 %), because
`σ_EM/dL` is **O(150–800) fringes** in the lit tier — **the EM prior is wide enough to absorb almost
any source displacement, so a wrong counterpart does not look displaced relative to a prior that was
never going to localise the fringe anyway.**

> **THIS CAVEAT TRAVELS BACKWARD ONTO D3.** Every D3 number — (a) = 100 %, (c) = 42/42 included —
> was computed in the **oracle** convention. **No co-registration number in this repo may be quoted
> as an achievable power without its implementable-form value beside it.**

**The constructive corollary.** The oracle/implementable gap **closes with σ_d**: D4-OBS is 1.6×
stronger in the VLBI tier (51.6 %) than in the lit tier (32.9 %). This is the *same* lever RING
identified (only sub-3-pc σ_d matters) and the same lever the onset map rewards. **Sub-3-pc
distances are doubly load-bearing: they buy detections AND they buy wrong-counterpart robustness.**

## CONVENTION — every quoted onset carries its floor's N and its fit error

Adopted 2026-07-12, from IGNITE-2 §2. Fresh D2-sized floors (N = 150, Gumbel α = 0.05) at the two
pre-registered onset cells land **8 / 2 nat ABOVE** the max-of-10 floors IGNITE ran under, and under
them **neither cell clears onset** (0.92 and 0.54 correct certifications/realisation, against the
>1 bar). **IGNITE's h\* was partly an artifact of the retired floor estimator, and NO
properly-calibrated onset exists anywhere in the modelled box.** An onset number quoted without its
floor's `N` and fit error is not quoted. (The 10-null floors at the other 22 cells carry ±2–18 nat
fit errors and cannot support an onset claim.)

## CHORUS — PRE-REGISTRATION (the mixed-eccentricity population campaign; queued, no compute yet)

**The open question IGNITE-2 flagged and nothing in this campaign has touched: every result above is
for a SINGLE-POPULATION source model. Nature supplies a MIXTURE.** CHORUS measures whether an
eccentric minority changes the certification arithmetic for the circular majority.

**Axes:** (fraction eccentric) × (e-distribution) × (N_CW). **Deliverables:** (1) *certified count vs
mix* — the onset surface re-cut as a function of the eccentric fraction, under criterion-v2.1 with
the fresh-floor convention; (2) **THE CLOCK-SHARING TEST — the campaign's reason to exist:** does a
single e ≈ 0.7 source, whose comb self-clocks (ATLAS's corner), lift the *circular* sources' pulsars
over the floor? The mechanism under test is that the eccentric source's harmonic comb pins `f_gw`
and `mc` for the array (SIREN's lag-diversity argument says short-lag pulsars pin frequency and
thereby free the long lags to carry chirp mass), so the clock it supplies may be a **shared** array
resource rather than a private one — i.e. certification may be a property of the POPULATION, not of
the source. If it is, every single-source no-go in this repo is scoped to a premise nature does not
satisfy. (3) **The capacity-vs-clock trade curve:** more sources raise the trials/confusion floor
(the `ln K_counted` term and the noise-lock that sets `floor(h)`) while the eccentric member lowers
the registration cost — CHORUS measures where the two cross.

**Machinery (all banked, nothing new to build):** WEAVE chirp-tied harmonic stacks for the eccentric
members + criterion-v2.1 scoring (fN operative, fALL blind-robust column, Gumbel α = 0.05 floors at
N ≥ 150, refit per cell — the mixture changes the null, so **floors are refit per mix, never
inherited**) + the soft loop (spec §3, `Q = Σ_p Σ_n q_p(n)·lnL`) as the reference implementation.
**Gates to be specified in the launch prompt.** The pre-registered STOP conditions and the
adoption/rejection conditions are to be fixed **before** any compute, per standing convention.

**Standing caveats that CHORUS inherits and must restate:** the EOB-tier validity limit (ATLAS —
the e ≳ 0.85 corner that would clear the absolute 0.003-dex floor is exactly the toy-invalid one);
the D1 wrong-counterpart hole (OPEN, un-closable by co-registration); the oracle/implementable
caveat on any co-registration statistic; and the 15-reals/cell sky-sampling error (GEO: the sky draw
dominates yield variance).
