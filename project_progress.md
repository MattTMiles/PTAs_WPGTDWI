# Pulsar Distance Likelihood Project — Progress Tracker

*Living document. Update the running log at the bottom each working session.*
Last updated: 2026-06-25

---

## 1. The project in one paragraph

PTAs constrain pulsar distances through the **phase-connected** continuous-wave (CW)
signal model: the pulsar-term phase depends on distance via `L(1 - cos μ)/c`, so the
CW likelihood is periodic in each pulsar distance with fringe spacing
`dL = c / [f_gw (1 - cos μ)] = λ_gw / (1 - cos μ)`. A single CW leaves the distance
degenerate across fringes; multiple CWs (different `dL` per source) break the
degeneracy because only the true distance phases all sources simultaneously. We are
investigating whether — and when — this can be turned into a usable distance
measurement, with an eye to the transition between resolvable CWs and a stochastic
background, and to phasing up the array.

Repo: `github.com/MattTMiles/PTAs_WPGTDWI`
Branches needed: `MattTMiles/discovery`, `MattTMiles/enterprise_extensions`.

## 2. The three prongs

1. **Can pulsar distances be optimised out?** (search/estimation, not sampling)
2. **The GWB ↔ CW transition region** — barely explored anywhere, high value.
3. **Improving sampling of pulsar distances when a CW is present.**

## 3. Collaborators & their positions (do not lose track of these)

- **Matt Miles** (Vanderbilt) — lead. Building JUG (JAX pulsar timing) in parallel;
  wants this integrated with PTA analysis pipelines that already sample noise.
- **David Hogg** (Flatiron) — **strongly wants prong 1 pushed and believes it can be
  done without sampling.** Position: "MCMC is the method of last resort"; know the
  likelihood before sampling; don't sample well-determined parameters; start with
  brute-force 2-D distance–distance grid scans on the two best pulsars, build up.
  Working on sky models that don't resolve into sources (spherical harmonics).
  Key remark: PTA angular resolution is `min(λ/D, λ/[cT])` → only need sky structure
  to spherical-harmonic degree ~10. **ACTION OWNER on keeping prong-1 optimisation
  alive: us, on Hogg's behalf.**
- **Will Farr** (Flatiron) — framed the three sub-questions (info / search / Bayes)
  and that **"there is provably no point doing 2 before 1, nor 3 before 2."** Says
  the pulsar term is almost always recoverable given noise assumptions; the *distance*
  is formally unrecoverable only in the truly stochastic / isotropic / short-
  correlation-time limit. The transition between the two limits is the interesting,
  poorly-understood regime. Suggested the integer-cycle / offset-phase distance jump
  (fast/slow split), ref arXiv:2506.10846; Neil Cornish used it but never published.
- **Steve Taylor** (Vanderbilt) — paper co-author elsewhere; localisation paper
  arXiv:2603.10120 partially explores the phase-connected effect. Likes hybridising
  optimisation inside a sampler (illegal moves during burn-in, then freeze for
  detailed balance). Wants it embedded in standard PTA frameworks for uptake.
- **Mihir Shetty** (NYU) — built JaxPINT (PINT ported to JAX) to model nonlinear
  timing effects without marginalising the distances. Produced 1-D and 2-D distance
  sweeps and a **Fisher-information-vs-N_CW** scaling plot (the "independence line").
  Open question he flagged: how to quantify when the false-peak-elimination
  transition happens.
- **Konstantin Leyde** (Flatiron) — interested in folding the (periodic) distance
  space; to look at interpreting the likelihood surface.
- Also cc'd: Max Isi, Andy Casey, Yacine Ali-Haïmoud.

## 4. Shared infrastructure (status: working)

- `CW_lnL_check/cw_helpers.py` — enterprise simulation + discovery likelihood.
  - `build_enterprise_pta`, `generate_injection_params` (scenarios: single,
    close_freq, close_sky, well_separated, realistic), `simulate`.
  - `MultiSourceDelay` — N CW sources sharing one `p_dist` per pulsar, vmapped.
  - `build_disco_likelihood` / `build_noisefree_likelihood` — Woodbury-marginalised
    (timing model + red noise + GWB GPs).
  - **`build_fast_scan_likelihood`** — caches the outer `(N_psr·2·n_comp)³` Cholesky
    once (GP hyperparams frozen at truth), ~50–100× speedup; enables N_psr=116 runs.
  - `scan_pulsar_distance`, `analyze_peaks`, `compute_mode_spacing`,
    `find_best_wrong_mode_in_prior`, `compute_joint_best_wrong_in_prior`.
- Metric throughout: **joint Δln L = lnL(truth) − lnL(all pulsars at best wrong
  fringe within ±3σ EM prior).**

## 5. Status by prong

### Prong 1 — optimise distances out
- **Explored:** Joint optimiser (scan conditional likelihoods across EM prior, keep
  top-N peaks, fix easiest→hardest, then let pulsars hop modes). Works with loud CWs,
  no noise: 20/20 recovered. With stochastic noise + GWB: ~75–80% of distances.
- **Not explored / open:** Honest *marginal* competing-mode metric (current metric is
  conditional, see Issues); coordinated multi-pulsar decoy search; Hogg's pure
  2-D grid → build-up programme as a standalone estimator; coherence proxy that
  avoids solving the whole sky per grid point.
- **Hogg wants this to be THE focus.** Keep a standing item to push the estimator.

### Prong 2 — GWB ↔ CW transition  ← CURRENTLY ACTIVE
- **Explored:** nb 05 = single realistic operating point (h~−14…−13.5, GWB at NG15,
  1 µs WN, 116 psrs). Mihir's Fisher-info-vs-N_CW plot. Mock-population Ncorr heatmap
  suggesting single-CW phase-up essentially never achievable for realistic pops.
- **Active work (this session):** self-contained Fisher-information transition
  calculation — conditional vs **marginal** distance information as sources go from
  few/loud (resolvable) to many/faint (confusion → stochastic). See §6.
- **Not explored / open:** systematic sweep at fixed total power; sensitivity to
  GWB mis-specification (frozen-at-truth GP is optimistic here); heterogeneous noise;
  the practical CW-strength / distance-precision threshold (goal "y").

### Prong 3 — sampling
- **Explored:** 3-phase annealed snap-sampler (`CW_node_sampling`): cool T 5000→1
  over 15k steps; 5k adapt (learn empirical CW covariance); 5k production. Proposals
  50% eigenmode / 30% distance(prior draw + ±0.6 dL grid of 30) / 20% joint CW.
  Newton-snap distances to conditional MAP after CW moves. Single & multi-CW (N=3).
- **Known failure:** locks distances reasonably but **struggles to find CW
  parameters**; no global CW move, no tempering swap (single-chain annealing can't
  recover a missed mode).
- **Not explored / open:** Farr integer-cycle/offset-phase distance parametrisation;
  global CW big-jump (draw f/sky from prior); parallel tempering ladder.

## 6. Active workstream — Prong 2 transition (details)

**Question (Farr's framing, info-only):** as we move from a few resolvable CWs to many
sub-threshold sources approaching an isotropic stochastic background, how does the
information about a pulsar's distance behave?

**Hypothesis:**
- *Conditional* distance info (all CW phases known) ≈ linear in N_source — Mihir's
  independence line.
- *Marginal* distance info (per-source phase/amplitude marginalised) rises while
  sources are resolvable, then **collapses** once sources confuse (≳1 per 1/T
  frequency bin) — Farr's zero-information stochastic limit.
- The gap between the two curves *is* the transition physics, and the turnover sets
  the practical threshold (goal "y").

**Method:** faithful non-evolving phase-connected CW residual (Earth + pulsar term),
realistic TOA grid + white noise, autodiff Fisher matrix over {L, per-source phase,
per-source log-amplitude}; marginal L-information = Schur complement. Sweep N at fixed
total characteristic strain across the band.

**Files:** `prong2_transition.py` (model + joint-array Fisher), figure
`prong2_distance_information_transition.png`.

**First result (toy model, 15 psr, 15 yr, 100 ns WN, band 3–20 nHz):**
- Confirmed the key mechanism: for a *single* pulsar, distance is perfectly degenerate
  with each source's pulsar-term phase → distance is only an *array* measurement, where
  the global Earth term pins the source phases. (Single-pulsar marginal info ≈ 0.)
- Conditional distance info ∝ N (reproduces Mihir's independence line); marginal info
  rises then turns over.
- `marg/cond` falls monotonically: ~0.95 (N=1) → ~0.53 (N=50) → ~0.18 (N=96).
- Fixed-total-power: conditional info ~flat while marginal collapses — the literal
  "CW fragmented into a background loses distance info" (Farr's limit), info-only.
- **`marg/cond` is nearly identical in fixed-per-source vs fixed-total** → the
  *recoverable fraction* is set by source crowding/geometry, not power distribution.
- Caveat: toy amplitude normalisation + uniform 100 ns WN + idealised geometry. Next:
  real array (sky + EM σ from feathers), heterogeneous noise, then port onto the real
  discovery likelihood (marginal info via autodiff of `build_fast_scan_likelihood`,
  marginalising CW params) to cross-check against the conditional joint-Δln L metric.

**Scaled result (Stage A, 116 psr, 15 yr weekly = 782 TOA, N log-grid 1→1000, 30 seeds,
both power modes; GPU, jax 0.4.28 venv, ~21 min total):**
- `marg/cond` (array-median): **0.99 (N=1) → 0.87 (N≈100) → 0.53 (N≈400) → 0.08 (N=1000)**;
  0.5-knee at N≈410 for the fiducial (15 yr, 3–20 nHz). Conditional ∝ N holds out to
  N=1000 (independence line); marginal rises, turns over, collapses.
- **Mode-independence is EXACT.** `max |Δ(marg/cond)|` between fixed-per-source and
  fixed-total = **0.000** at every N (identical to machine precision); no break found up
  to N=1000. The recoverable *fraction* is purely geometric — independent of how CW power
  is distributed. (In fixed-total, conditional info is ~flat ≈2e11 while marginal
  collapses: the literal "CW fragments into a background, distance info lost" picture.)
- **Knee tracks N\* = T·Δf in scaling, offset by the array.** Varying T (10/15/20 yr,
  band 3–20 nHz): knee = 288/412/555, N\* = 5.4/8.0/10.7 → **knee/N\* ≈ 52 (constant)**.
  Varying band (3–12/3–20/3–40 nHz, T=15 yr): knee = 255/412/716 → knee/N\* = 60/51/41
  (drifts). So the **T-scaling** of the confusion onset is captured by T·Δf; the
  **bandwidth-scaling** only roughly. The ~50× offset is the array-resolution boost: 116
  pulsars pin the global source phases ~50× better than naive one-source-per-1/T-bin
  counting, pushing confusion onset to N≈50·(T·Δf). **Practical threshold (goal "y"):
  distance info survives marginalisation up to N ≈ 50·T·Δf for a 116-pulsar array.**

**Optimisation (this session, GPU):** padded source array to fixed N_max + boolean mask
→ one compiled XLA shape serves all N (no per-N recompile; padded sources contribute
exactly zero, verified padded-vs-unpadded to 1.9e-15). Schur source-block inverse moved
from ridge-regularised solve to **eigh pseudo-inverse** (rcond 1e-10·max-eig) — exact on
the rank-deficient padded/confusion block. Source-extras + distance gradients computed
**analytically** (closed-form, validated vs jacfwd to 6e-16) instead of a 2·N_max-tangent
jacfwd; per-call 16 s → 0.7 s (~22×), which is what lets N_max=1000 fit in 24 GB. Fisher
assembled per-pulsar (F_LL diagonal) to keep memory O(n_src²) not O(n_param·n_data).

**Stage A.1 probes (robustness + physics readouts; ~40 min run):**
- **rcond robustness.** Swept the eigh-pinv threshold rcond ∈ {1e-8, 1e-10, 1e-12}:
  0.5-knee = **422.6 for all three (0.0% spread)**. Marginal is insensitive to the cut →
  production rcond = **1e-10**. The collapse is a real loss of source-block rank, not a
  regularisation artifact.
- **Marginal is an ARRAY effect, not per-pulsar.** F_ss is built as (1/σ²) Σ_a H_aᵀH_a
  over ALL pulsars (the shared global Earth term pins the source phases). Verified: a
  pulsar marginalising on its OWN data block alone gives marginal = **8.2e-10 · cond ≈ 0**
  (its distance is degenerate with its own source phase/amp — all single sinusoids span
  the same 2-D {sin,cos} space as the 2 nuisances), while the full-array marginal is
  0.93·cond. So marg≠0 requires the cross-pulsar sum. (Confirms the single-pulsar
  degeneracy claim above, now demonstrated rather than asserted.)
- **σ_L in parsec (J0437-like, L=156.8 pc, idx 0).** fixed_total: σ_L = 1/√I_marg rises
  **0.0034 → 0.0076 pc** across N=1→1000 (≈3.3× degradation as the CW fragments);
  fixed_persource σ_L only *improves* (0.034 → 0.0014 pc, per-source power keeps adding
  info). **Neither crosses 1 pc nor the ~0.25 pc EM/VLBI prior in N≤1000** — the schematic
  ζ=h/2πf normalisation is far too optimistic in absolute terms. Only the *shape* (3.3×
  relative degradation) transfers; absolute crossings need the real discovery-likelihood
  amplitudes + heterogeneous noise (the porting step). **Caveat for the writeup: do not
  quote absolute σ_L from the toy.**
- **Per-source optimal SNR at the knee** (N=423): median **≈195** (√ of Σ snr² over TOAs
  ×116 pulsars) — again schematic-amplitude-inflated; relative use only.
- **Population amplitude mode** (k=3 loud + (N−k) faint, h_loud/h_faint=10): marg/cond
  departs the uniform curve only mildly — first |Δ|>0.02 at **N=158, max Δ=−0.025**. The
  exact amplitude-degeneracy IS broken, but the 3 loud sources dominate the recoverable
  *fraction*, so the big effect is on absolute conditional info: cond **plateaus**
  (~2.8e10 at N=1000 vs 2.1e12 uniform) because the faint sources carry h²/100 power.
  Lesson: mode-independence of marg/cond is an artifact of equal amplitudes; with a
  realistic luminosity function the curve is set by the handful of loud sources.
- marg/cond at key N (fiducial): N=1 → **0.991**, N=10 → 0.986, N\*=8 → 0.987,
  knee=423 → **0.500**, N=1000 → **0.084**.

### Two distinct transitions (do not conflate)

A. SOURCE CONFUSION / RESOLVABILITY  ← what Stage A measured
   - Controlled by N vs N* = T·Δf (number of resolvable 1/T frequency slots),
     with a ~50× array-resolution boost from 116 pulsars seeing each source
     through different (1-cos mu).
   - An ESTIMATION limit: the array can't separately determine each source's
     phase once sources crowd within a resolution element, so marginalising the
     phases destroys distance info. Conditional info (phases known) never
     collapses (∝N forever); all loss is in resolvability.
   - RECEDES with longer T. Knee scales as T·Δf (verified, knee/N* ≈ 52).

B. SIGNAL COHERENCE / DECOHERENCE  ← Farr's transition, PROBED in Stage A.2
   - Controlled by source coherence time t_c vs Earth-pulsar light-travel time
     tau_p = L(1-cos mu)/c (thousands of yr for a kpc pulsar).
   - INTRINSIC: a stochastic field with t_c << tau_p has no expected Earth–pulsar
     phase correlation, so distance info -> 0 regardless of how well sources are
     resolved.
   - ~~INDEPENDENT of T~~ — **CORRECTED by Stage A.2.** Was the pre-registered guess;
     the Bayesian-Fisher measurement shows decoherence loss DOES grow with T (via
     SNR accumulation, see Stage A.2). It exists at N=1 (the real distinction from
     confusion) but is not T-independent. Our monochromatic sources have infinite
     t_c, so the Stage-A model structurally cannot reach this limit without dphi.

DISCRIMINATOR (corrected): the clean separator is the **N-axis, not T**. At N=1 there
is zero confusion, yet a full decoherence transition exists (Stage A.2 Exp1). Both
transitions happen to scale with T (confusion via T·Δf, coherence via SNR²∝T), so
growing T does NOT isolate them — N=1 + finite t_c does.

CAVEAT: in a realistic SMBHB background the two roughly coincide (the background's
stochasticity IS many unresolved binaries), but they are set by different ratios.

METRIC NOTE: for spiky (realistic) populations, marg/cond is misleading — it stays
high while absolute information is small and carried by a few loud sources
(Stage A.1: pop conditional plateaus at 2.8e10 vs 2.1e12 uniform, ratio barely
moves). Switch the headline to absolute sigma_L once Stage C gives trustworthy units.

**Stage A.2 — the COHERENCE transition (model + measurement).** Gave each source a
finite coherence time via an INDEPENDENT per-(source,pulsar) pulsar-term phase offset
dphi_{s,p} (default 0 = coherent), with a Gaussian decoherence prior of variance
sigma_phi² = tau_{s,p}/t_c, tau_{s,p}=(L_p/c)(1−cos mu). Bayesian Fisher
F_total = F_data + Pi (Pi=diag(1/sigma_phi²) on the dphi block). cond = dphi pinned;
marg = Schur-complement {phase0, log10h, dphi}. Run: `python prong2_transition.py
coherence` (N=1, 116 psr, 15 yr weekly; ~9 min, dominated by the confusion-knee contrast).
- **Gate (i):** t_c→∞ (sigma_phi²→0 ⇒ dphi pinned) reproduces the coherent N=1 value
  **0.9931** to **3.4e-6** (asserted). NB this is the 116-psr number; the 15-psr toy was
  ~0.95.
- **Exp 1 (headline, clean isolation):** at N=1 there is ZERO confusion, so the
  marg/cond vs tau_p/t_c curve is Farr's decoherence uncontaminated. It falls
  **0.99 → 0** across tau_p/t_c = 1e-3 → 1e3; **0.5-crossing at tau_p/t_c = 0.0144**
  (linear form). Decoherence bites early (t_c ≈ 70·tau_p) because dphi_{s,p} is *exactly*
  anti-parallel to ∂r/∂L per baseline ⇒ marg/cond = 1/(1 + (tau/t_c)·SNR²).
- **Gate (ii) shape-robustness:** saturating form sigma_phi²=(π²/3)(1−e^{−tau/t_c})
  moves the 0.5-crossing to 0.0044 — a **0.30× shift** (factor ~3). The transition
  LOCATION (tau_p/t_c ≪ 1, order-of-magnitude) is robust; the exact coefficient and
  SHAPE are model-dependent. **Absolute t_c values are model-dependent — only the
  existence and the location-scaling (∝ tau_p) are the claim.**
- **Exp 2 (falsified hypothesis, reported honestly):** fixed t_c (mid-transition),
  vary T=10/15/20/30 yr → coherence marg/cond = **0.62/0.52/0.45/0.35**, i.e. it FALLS
  with T as **1/(1+T/16 yr)**. The pre-registered guess that decoherence is
  T-independent is **wrong**: the constant per-baseline offset is increasingly resolved
  as SNR²∝T grows, so the penalty is (tau/t_c)·SNR². Both transitions therefore scale
  with T; **the clean discriminator is the N-axis (Exp1), not T-scaling** (confusion knee
  over the same T: 288/412/555/756, ∝T·Δf).
- **Files:** coherence model in `prong2_transition.py` (`coherence` CLI),
  `prong2_coherence.npz`, `prong2_coherence_transition.png` (2-panel: Exp1, Exp2).

## 7. Cross-cutting issues / caveats flagged

1. **Conditional vs marginal.** `compute_joint_best_wrong_in_prior` finds each
   pulsar's best wrong mode independently (others + CW params at truth), then sets all
   wrong at once. Over-states separability: forbids decoy cooperation and conditions on
   perfect CW recovery. Same trap as σ_MLE < σ_post in the JUG binary work. Need a
   marginal / coordinate-ascent version.
2. **Frozen-at-truth GP hyperparams = upper bound on distance info.** Acceptable for
   prongs 1/3; wrong simplification for prong 2 (GWB is part of what we separate from).
3. **Sampler distance move can't migrate fringes** (only ±0.6 dL around a fresh prior
   draw; ~1/100s chance to hit the right fringe). Likely a primary cause of slow
   distance convergence, independent of the CW-finding problem.
4. **No global CW proposal / no tempering** → missed CW modes unrecoverable.
5. **Homogeneous 1 µs WN (nb 05)** overstates weak pulsars; real arrays are dominated
   by a few good pulsars — also the lever for the Deller/VLBI targeting science case.

## 8. Running log

- **2026-06-25** — Onboarded full repo + email chain. Wrote this tracker. Built
  prong2_transition.py (joint-array conditional/marginal distance Fisher) and produced
  the first transition figure; logged the single-pulsar degeneracy finding and the
  marg/cond collapse. Next: scale on GPU + real array + port onto discovery likelihood.
  Logged issues 1–5. (Claude + Matt)
  prong-2 as active workstream; started the conditional-vs-marginal distance-
  information transition calculation. Logged issues 1–5. (Claude + Matt)
  
- 2026-06-25 (Stage 0, cronus/4090) — Toy transition reproduced GPU-bit-identical to
  CPU (marg/cond 0.950 @N=1, 0.178 @N=96; cond/N flat = cond ∝ N). 5.3s total,
  compile-dominated, sub-second/N → Stage A is cheap. Built isolated jax 0.4.28 venv.
  Added plot_prong2.py (was missing). Flagged discotech-GPU risk for Stage C.

- 2026-06-25 (Stage A, cronus/4090) — Optimised + scaled prong2_transition.py for GPU.
  Padding+mask (one compiled shape ∀N), eigh pseudo-inverse Schur, analytic gradients
  (22× faster; jacfwd kept as validator). Correctness gates: padded==unpadded 1.9e-15,
  analytic==autodiff 6e-16. Full run 116 psr / 15 yr weekly / N=1→1000 / 30 seeds / both
  modes + T,band knee diagnostics in ~21 min. Findings: marg/cond 0.99→0.08, 0.5-knee
  ≈410; marg/cond EXACTLY mode-independent (Δ=0.000, no break to N=1000); knee ∝ N\*=T·Δf
  under T (knee/N\*≈52 const) but offset ~50× by array resolution (band-scaling only
  rough). Goal-"y" threshold: distance info survives to N≈50·T·Δf for 116 psr. Refreshed
  3-panel figure + npz. Files now under CW_transition/. (Claude + Matt)

- 2026-06-25 (Stage A.1, cronus/4090) — Robustness + physics probes on prong2. (1) rcond
  sweep {1e-8,1e-10,1e-12}: knee=422.6, 0% spread → kept production rcond=1e-10. (2)
  Asserted F_ss = Σ over ALL pulsars; single-pulsar-only marginal=8.2e-10·cond≈0 →
  non-zero marginal is a genuine cross-pulsar (shared-Earth-term) effect. (3) σ_L(pc) for
  J0437-like: fixed_total degrades 0.0034→0.0076 pc, persource improves; neither crosses
  1 pc / 0.25 pc EM in N≤1000 (toy amplitude too optimistic — DON'T quote absolute σ_L).
  (4) median per-source SNR at knee≈195 (schematic). (5) population mode (3 loud+faint,
  10×): marg/cond departs uniform only at N=158 (Δ=−0.025); real break is cond plateauing
  on the loud sources → mode-independence is an equal-amplitude artifact. Figure now 2×3
  (added σ_L + rcond panels, population overlays); npz extended. (Claude + Matt)

- 2026-06-25 (Stage A.2, cronus/4090) — Built the COHERENCE transition (Farr): per-
  (source,pulsar) pulsar-term phase offset dphi with Gaussian decoherence prior
  sigma_phi²=tau_p/t_c; Bayesian Fisher F_data+Pi; marg = Schur out {phase0,logh,dphi}.
  Gate (i) t_c→∞ reproduces coherent N=1 = 0.9931 to 3.4e-6 (asserted). Exp1 (N=1, zero
  confusion = clean isolation): marg/cond 0.99→0 vs tau_p/t_c, 0.5-crossing at
  tau_p/t_c=0.0144; gate (ii) saturating form shifts it 0.30× (location robust, shape
  model-dependent). Exp2 FALSIFIED the T-independence guess: coherence marg/cond falls
  1/(1+T/16yr) because the constant baseline offset is resolved as SNR²∝T → both
  transitions scale with T; clean discriminator is N (coherence full at N=1 where
  confusion absent), not T. Corrected §6 "two transitions" accordingly. New files
  prong2_coherence.{npz,png} via `prong2_transition.py coherence`. (Claude + Matt)

## 9. Environment (cronus)
- GPU box cronus = NVIDIA RTX 4090, driver 550.120.
- Toy/Fisher work (prong2_transition.py): isolated venv with jax 0.4.28
  (CUDA 12.4 / cuDNN 8.9), driver-compatible. Latest jax[cuda12] FAILS cuDNN init
  on this driver. Activate: source <scratch>/env/bin/activate.
- discovery lives in the shared `discotech` env (separate jax) — DO NOT modify it.
- OPEN RISK for Stage C: confirm discotech's jax can init the GPU; if it shares the
  broken jax, Stage C Hessian is CPU-bound or needs a pinned discovery venv.