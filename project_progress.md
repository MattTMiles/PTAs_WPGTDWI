# Pulsar Distance Likelihood Project — Progress Tracker

*Living document. Update the running log at the bottom each working session.*
Last updated: 2026-07-08

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
  few/loud (resolvable) to many/faint (confusion → stochastic). See §6. **TOY PHASE
  CLOSED** (Stages A/A.1/A.2/A.3): confusion and coherence transitions mapped and shown
  to COMPOUND (not factorise). Next = Stage C (real discovery likelihood, absolute σ_L).
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

### Two distinct transitions (do not conflate) — REVISED after Stage A.2

A. SOURCE CONFUSION / RESOLVABILITY  (Stage A)
   - Controlled by N vs N* = T·Δf, with ~50× array-resolution boost (116 psr,
     varied 1-cos mu). An ESTIMATION limit: can't separate source phases once
     crowded. Conditional info never collapses (∝N); loss is all in resolvability.
   - Requires N>1. Knee scales with T·Δf.

B. SIGNAL COHERENCE / DECOHERENCE  (Stage A.2) = Farr's transition
   - STOCHASTIC decoherence only (random per-baseline phase wander, var σ_φ²);
     NOT deterministic chirp (a chirp is recoverable info and HELPS — opposite sign).
   - Closed form (N=1): marg/cond = 1/(1 + SNR²·σ_φ²).
   - MODEL-INDEPENDENT LAW: distance info halves when rms Earth-pulsar phase wander
     reaches ~1/SNR rad (verified: linear & saturating σ_φ² forms cross at the same
     SNR²σ_φ²≈1; the t_c value is model-dependent, do NOT quote it as physical).
   - Two absolute regimes: σ_L → dL/(2π·SNR) (coherent, SNR-limited) ;
     σ_L → dL·σ_φ/(2π) (decoherent, SNR-INDEPENDENT floor).
   - Consequence: beyond SNR ~ 1/σ_φ, more strain does NOT improve distances.
     This is the rigorous form of "a strong CW can phase up the array" — it can't,
     past the wander floor.
   - Present at N=1 (zero confusion). 

BOTH transitions scale with T (confusion via T·Δf; coherence via SNR²∝T) — my earlier
"coherence is T-independent" was WRONG (physical wander is T-independent, but its
DETECTABILITY scales with SNR). T-scaling does NOT separate them.
DISCRIMINATOR: the N-axis. Orthogonal knobs = N (confusion) and t_c (coherence).
CAVEAT: in a realistic background the two ARRIVE TOGETHER (a background IS many
unresolved sources), so N=1+finite-t_c is a conceptual isolation, not a real scenario.

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

**Stage A.3 — the N x t_c MAP (toy capstone): confusion x coherence do NOT factorise.**
Combined model: N fixed-per-source sources, each with a per-(source,pulsar) phase offset
dphi_{s,p} ~ N(0, tau_{s,p}/t_c). conditional = {phase0,log10h,dphi} all pinned;
marginal = Schur out all three; R(N,Y)=marg/cond. Narrow band b3-12 (N*=4.26, knee~222
reachable at N<=256). Run: `python prong2_transition.py ntcmap` (~90 s).
- **Efficient marginalisation (required):** dphi_{s,p} enters only pulsar p, so F_data's
  dphi-block is BLOCK-DIAGONAL over pulsars (116 independent N x N blocks). Marginalise
  dphi by 116 per-pulsar N x N solves -> F_eff over {L,phase0,log10h}, then the dense 2N
  block. **Gate:** structured == full dense single-Schur Fisher at N=8 to **5.9e-11**
  (asserted).
- **Edge gates (anchor the 2-D map to the 1-D stages):** y->0 (t_c->inf) slice reproduces
  the Stage-A fixed-per-source confusion R_conf(N) to **4.2e-8**; x=1 (N=1) slice
  reproduces the Stage-A.2 coherence path to **4.7e-15**. Both asserted.
- **Invariant axes:** x = N/N* (confusion), y = Y (coherence). Y normalised to the N=1
  half-info point (Y := t_c_50/t_c, t_c_50~1.3e13 s) so Y=1 is the coherence midpoint --
  this absorbs the O(geometry) antenna 1/(1-cos mu) weighting that makes the bare
  SNR²·sigma_phi² mis-centred (consistent w/ Stage-A.2 gate ii: location is the invariant,
  absolute t_c model-dependent). R is median-over-realisations of marg/cond (the stable
  reduction; do NOT use ratio-of-medians -- they differ at N=1, a bug we hit and fixed).
- **HEADLINE — they INTERACT, sub-multiplicatively (COMPOUND).** The R=0.5 contour is
  **CURVED** (bows toward the origin), not a rectangle -> the thresholds are not
  independent, and R falls BELOW the factorised product in the confusion-dominated
  interior. Crowded sources AND decohered pulsar terms destroy distance information
  faster than either limit alone -- the realistic SMBHB background (many unresolved +
  finite-coherence binaries) sits in this compounding region, so the practical threshold
  is WORSE than the independent-threshold estimate.
- **ROBUSTNESS RE-CHECK (no re-grid).** The first pass measured the residual with
  median-of-ratios per cell, which does NOT commute with the product. Recomputed two
  COMMUTING ways from the saved per-cell median info: (1) Rbar = median(I_marg)/median
  (I_cond), test |Rbar - Rbar(N,0)Rbar(1,Y)|; (2) info-level (I_marg vs I_cond·Rbar(N,0)
  Rbar(1,Y), normalised by I_cond = same as (1)). Result: peak deviation drops
  **0.331 -> 0.202** (at N/N*=27, Y=0.14) but stays **strongly negative (< -0.1)** through
  the confusion-dominated interior -> **compounding is PHYSICAL, not a reduction
  artefact**; median-of-ratios merely inflated the magnitude ~50%. The commuting map also
  reveals SIGN STRUCTURE: a weak POSITIVE lobe at N/N*<1 (mild de-confusion -- random
  pulsar-term phases slightly aid separation before confusion sets in) and the dominant
  NEGATIVE (compounding) lobe at N/N*>~3. Verdict kept: COMPOUND, corrected peak ~0.20.
  (Used median info, not mean, because mean needs per-seed data = a re-grid, prohibited;
  median is an equally-valid commuting reduction and the sign/magnitude are decisive.)
  Figure: `prong2_Ntc_map_corrected.png`; arrays added to `prong2_Ntc_map.npz`
  (`Rbar`, `resid_commuting`, `max_factor_dev_commuting`).
- **Files:** `prong2_transition.py` (`ntcmap` CLI), `prong2_Ntc_map.npz`,
  `prong2_Ntc_map.png` (heatmap of R + factorisation residual).

### Stage C — real discovery likelihood (COMPLETE, 2026-07-02)

Ported marginal distance info onto `build_fast_scan_likelihood` (jug-gpu env, jax
0.10.1 on the 4090); absolute sigma_L in pc; fringe-identification split. Deliverables
D0-D6 all done.

**>>> STAGE C HEADLINE (goal "y" in real units, honestly).**
For the phase-connected CW model on the real NG15-like 116-pulsar array:
1. **Within-fringe distance precision is excellent and cheap:** marginal sigma_L
   ~0.07-0.27 pc (equal strain, improving with N_CW), below the (wide, median ~235 pc)
   feather EM priors for essentially all pulsars already at N_CW=1. **But this is NOT a
   distance measurement** -- it is the curvature within one fringe.
2. **Fringe IDENTIFICATION is the binding constraint, and it is hard.** A SINGLE CW
   identifies ZERO fringes (dlnL_a ~ 0 vs ln K ~ 8; D3). Equal-strain sources break the
   degeneracy only slowly: genuine measurements (class-i) first appear at N_CW=16
   (median 5/116, tail-driven; D4). A REALISTIC luminosity function (3 loud + 13 faint)
   gives class-i = 0 even at N_CW=16 -- only loud sources break fringes, so it behaves
   like ~3 loud CWs. **A handful of loud SMBHBs cannot phase up the array.**
3. **The conditional metric barely overstates here (marg/cond ~1.0-1.04 up to N_CW=16):**
   the real likelihood is deep in the RESOLVABLE regime, far below the toy confusion knee
   (~50 T Df ~ hundreds). The real likelihood does NOT reproduce the toy's
   resolvable->confusion turnover in the achievable range; what it ADDS is the discrete
   fringe-ID layer the toy Fisher never saw.
4. **Frozen-GWB optimism is small (<=9% for +/-0.5 dex; D5).** Dominant realism levers
   are the source luminosity function and real-noise-vs-Asimov, not the GWB freeze.

**WHAT CHANGED vs the project's prior (conditional joint-Delta lnL) metric.**
- The old headline metric (`compute_joint_best_wrong_in_prior`) is (i) CONDITIONAL
  (CW params at truth, no marginalisation) and (ii) computed with a find_peaks
  prominence floor of 0.5 that returns "no wrong mode" whenever the fringe modulation
  is < 0.5, leaving that pulsar at truth -> **systematically OPTIMISTIC** about
  separability. Blast-radius triaged across all 6 CW_lnL_check nbs; nb-05 rechecked with
  a direct fine-grid scanner: corrected joint dlnL is 7-45% MORE negative (less
  separable) at every N_CW=1-8 (sign unchanged, margin firmer).
- Stage C replaces this with (i) the honest MARGINAL sigma_L (Schur-complement of the
  full joint CW block out of the zero-noise Fisher, H=-F exactly, eigh-pinv discipline,
  rcond-stable) AND (ii) a SEPARATE, prominence-free fringe-identification test
  (dlnL_a vs ln K_a). The key conceptual correction: **within-fringe precision and
  fringe identification are DIFFERENT questions**; the old conditional metric conflated
  "the likelihood curves sharply here" with "the distance is measured", and its
  prominence bug additionally over-credited identification. In real units the distance
  is measured only when BOTH hold (class-i) -- which at achievable N_CW is a small,
  loud-source-dominated minority, and zero for a realistic population.

**Files created (all under `CW_transition/`):**
`stagec_d0_derivative_gate.py`, `stagec_hessian.py` (D1), `stagec_fisher.py` (D2, the
amortised zero-noise Fisher + sigma_L), `stagec_d3.py` (D3 fringe split),
`stagec_d4.py` (D4 N_CW sweep), `stagec_d5.py` (D5 GWB freeze),
`stagec_nb05_recheck.py` (blast-radius closure). Result archives: `stagec_d0_results.npz`,
`stagec_d2_results.npz`, `stagec_d3_results.npz`, `stagec_d4_results.npz`,
`stagec_d5_results.npz`. Env: jug-gpu (jax 0.10.1), XLA autotune disabled (shared GPU).

**Deliverable detail (gate numbers, tables) follows.**

**D0 — derivative validation gate (PASS, 2026-07-02).** 5 psr / 1 CW (scenario
"single", log10_h=-13.5, GWB at NG15 -14.6, equad -6, 14 components, seed 1234),
at injected truth; lnL(truth)=22495.785026. Script
`CW_transition/stagec_d0_derivative_gate.py`, results `stagec_d0_results.npz`.
- Gate (a) grad autodiff vs central FD over [5 distances + 8 CW params]:
  **max rel err 6.1e-8** (< 1e-5, PASS).
- Gate (b) Hessian autodiff vs FD-of-autodiff-grad: distance block **2.0e-8**,
  distance-CW cross **6.4e-7** (< 1e-4, PASS; CW block info-only 6.9e-8).
- Gate (c) cho_solve path differentiated: rebuilding with frozen GWB log10_A
  +0.1 changes grad by max 657 (7.6e-2 of max|g|) — the cached-Cholesky solve is
  in the differentiated graph. PASS.
- **Method note (matters for any future FD checks):** naive steps 1e-5·|x| FAIL
  (rel err up to 1.6 on log10_fgw) — NOT an autodiff problem but FD truncation:
  the pulsar-term phase Phi = 2·pi·L/dL ~ 1.6e4 rad makes the logL oscillatory, so
  the natural scale is the 1-rad phase scale (dL/2pi for distances; 1/Phi_max for
  sky/log10_fgw/log10_mc). With 1e-2 of that scale + Richardson extrapolation, all
  central-FD errors shrink 4.00x from h to h/2 (pure O(h^2) truncation) and land at
  1e-8–1e-7 vs autodiff.
- Warnings under jax 0.10.1: the known `cho_solve` FutureWarning (batched-1D
  deprecation; value+derivatives now both validated), JAXopt deprecation,
  pkg_resources deprecation. No behaviour changes observed.

**D1 — Hessian infrastructure at scale (PASS, 2026-07-02).** Module
`CW_transition/stagec_hessian.py`. `stagec_hessian(problem, n_cw, method=...)`
returns the Hessian of the fast-scan logL at truth over [116 distances + 8·N_CW
CW params]. Two paths: `dense` (jax.hessian, full matrix) and `hvp`
(jax.jvp-over-jax.grad, chunked unit tangents; assembles ONLY the distance block
H_dd 116×116 and dist-CW cross H_dc 116×8N_CW — the blocks D2 needs).
- **GATE (5 psr / 1 CW): dense vs HVP [H_dd|H_dc] = 2.5e-14** (< 1e-8, PASS);
  dense vs the D0 FD-validated Hessian = 7.8e-14. HVP is exact, not approximate.
- **116-psr benches** (nb-05 realistic config; per-process peak device memory via
  memory_stats, BFC allocator, preallocate off):

  | N_CW | method | n_sel | peak GiB | warm eval | cold (compile) |
  |------|--------|-------|----------|-----------|----------------|
  | 1  | hvp   | 124 | 0.44 | 1.09 s | 340 s |
  | 4  | hvp   | 148 | 0.62 | 1.50 s | 407 s |
  | 16 | hvp   | 244 | 1.12 | 2.67 s | 508 s |
  | 1  | dense | 124 | 2.35 | 0.32 s | 442 s |
  | 4  | dense | 148 | 3.42 | 0.38 s | 489 s |
  | 16 | dense | 244 | **OOM** (single 12.32 GiB alloc) | — | — |

  Dense warm eval is faster (one fused call) but its forward-over-reverse
  transient scales ~n_sel² and OOMs by N_CW=16 on the 24 GB card (as the D1
  carry-forward note predicted — the transient, not the final matrix, is the
  limit). HVP is the production path: peak grows ~linearly in N_CW (assembles
  only 116 distance rows), 1.12 GiB at N_CW=16 → headroom to N_CW≈100+.
  **Use dense only for the small gate; use hvp for D2+.**
- Distance-block diagonal is finite for all 116 pulsars; count with NEGATIVE
  curvature (well-constrained direction) rises 78→95→110 as N_CW 1→4→16 (more
  sources phase up more pulsars). The remaining positive-diagonal pulsars are
  near-degenerate distances — expected; D2's Fisher = −H handles sign per pulsar.
- **Co-tenancy note (cronus is a SHARED GPU):** runs coincided with an external
  job (`smoke_ensemble.py` ~5 GB) + another user (~1.4 GB). XLA autotuning
  profiling OOMs under a co-tenant spike and kills the compile, so the bench
  disables it (`XLA_FLAGS=--xla_gpu_autotune_level=0`) — default kernels, robust
  to co-tenancy, peak numbers unaffected (per-process BFC pool). Cold-compile
  times (340–510 s) are one-time per (N_CW, method) shape; warm eval is what D2+
  pays per realisation.

**D2 — conditional vs marginal sigma_L, single realistic CW (PASS, 2026-07-02).**
Module `CW_transition/stagec_fisher.py`; results `stagec_d2_results.npz`.
116 psr, N_CW=1, nb-05 realistic operating point (real sky/TOA errors, GWB GP at
NG15 A=-14.6 gamma=13/3, h in [-14,-13.5], f in 10-32 nHz), 10 CW-GEOMETRY draws.

- **Amortisation (carry-forward note 2).** The fast-scan residual is a callable
  `y(params)=data-delay(params)`; only `data` is baked. `build_fisher_amortised`
  makes `data` a RUNTIME arg (via discovery's callable-residual vary path), so ONE
  compiled Hessian graph serves every draw. **Compile once 465 s, then 0.10 s warm
  per draw.** GATE: amortised zero-noise Hessian == D1 dense Hessian (reordered) to
  **4.6e-16** on the 5-psr case.
- **Zero-noise Fisher (carry-forward note 1).** Data = injected CW at truth, NO
  noise draw, but the real heterogeneous noise covariance (white + frozen GWB GP +
  frozen RN GP, rn=30 comp) still weights the likelihood. residual(truth)=0 =>
  H = -Fisher exactly. **Neg-curvature count = 116/116 on EVERY draw** (NSD by
  construction; zero degenerate pulsars). NB this differs from D1's noisy OBSERVED
  Hessian (78/95/110 for N_CW=1/4/16) — the D1 counts were the observed Hessian on
  noisy data, not the Fisher; on zero-noise data the count is the full 116.
- **sigma_L definitions (per pulsar, both conditional on other distances):**
  conditional_i = 1/sqrt(F_dd[i,i]) (CW known); marginal_i = 1/sqrt(F_dd[i,i] -
  (F_dc F_cc^{-1} F_cd)_ii) = 1/sqrt(diag of the Schur complement) (CW
  marginalised). Only the 8x8 CW block F_cc is inverted, via the toy's eigh
  pseudo-inverse (Stage A.1). **rcond sensitivity (1e-8 vs 1e-12): 2.7e-2**
  (stable). [An earlier pass used the FULL 116x116 pinv(F_marg) — marginalising
  over the OTHER distances too — which is near-singular for a single CW and gave a
  254% rcond spread; that is NOT the metric. The per-pulsar Schur-diagonal is the
  toy-consistent, stable quantity.]

  10 best pulsars (median over draws, PARSEC):

  | pulsar | marg_pc | cond_pc | EM_pc | marg/cond | marg<EM |
  |--------|---------|---------|-------|-----------|---------|
  | J1713+0747 | 0.097 | 0.093 | 19   | 1.06 | YES |
  | J1045-4509 | 0.100 | 0.099 | 68   | 1.01 | YES |
  | J0711-6830 | 0.119 | 0.118 | 21.2 | 1.01 | YES |
  | J0030+0451 | 0.125 | 0.124 | 3.6  | 1.01 | YES |
  | J2317+1439 | 0.125 | 0.116 | 290  | 1.09 | YES |
  | J1640+2224 | 0.150 | 0.146 | 95   | 1.04 | YES |
  | J2043+1711 | 0.150 | 0.146 | 110  | 1.03 | YES |
  | B1937+21   | 0.151 | 0.136 | 200  | 1.10 | YES |
  | J2145-0750 | 0.164 | 0.161 | 22   | 1.00 | YES |
  | J1125-6014 | 0.167 | 0.165 | 197.6| 1.02 | YES |

- **HEADLINE: median 116 of 116 pulsars have marginal sigma_L below their EM prior**
  (range 115-116 across draws). Median marginal sigma_L per draw ~0.11-0.55 pc; best
  pulsars ~0.10 pc. EM priors span ~1-290 pc.
- **marg/cond distribution: median 1.008, 16/84 = 1.001/1.033, max 12.95.** At
  N_CW=1 the CONDITIONAL METRIC DOES NOT OVERSTATE — 116 pulsars pin the single CW's
  8 params so tightly that marginalising them barely inflates per-pulsar distance
  precision. (Consistent with the toy's N=1 array-marginal info ratio ~0.93 =>
  sigma ratio ~1.03.) The overstating is expected to appear as N_CW grows and
  sources confuse — that is the D4 sweep.
- **CRITICAL CAVEAT (stated for D3): sigma_L here is WITHIN-FRINGE only.** "marginal
  sigma_L < EM prior" is NOT yet a distance measurement — it presumes the correct
  distance fringe is identified. Fringe identification is the separate discrete
  problem quantified in D3. Do not read the 116/116 headline as "116 distances
  measured".

**D3 — the fringe split (THE HONEST HEADLINE) (2026-07-02).** Module
`CW_transition/stagec_d3.py`; results `stagec_d3_results.npz`. Same 116 psr, N_CW=1,
10 zero-noise (Asimov) geometry draws as D2 (seeds 2000+d).

- **Threshold (given).** Pulsar a fringe-identifiable iff joint dlnL(true fringe vs
  best wrong fringe in +/-3 sigma_EM) > ln(K_a), K_a = candidate fringes in window =
  2*3*sigma_EM/dL_a (likelihood gap must beat the trials factor). Also a stricter
  dlnL > ln(K_a)+3.
- **Method note (a bug caught + fixed).** `find_best_wrong_mode_in_prior`'s
  find_peaks has a prominence FLOOR of 0.5, but for a single CW the ENTIRE fringe
  modulation across the prior is < 0.5 (see below) -> it registers NO competing peak
  and returns nan, which naively reads as "perfectly identified". That is BACKWARDS.
  Replaced with `fringe_deltalnL`: evaluate lnL directly at the fringe CENTERS
  L0+k*dL (the likelihood is periodic in the pulsar-term phase 2*pi*L/dL), fixed
  512-eval budget per pulsar for one XLA shape; report the TRUE K.
- **Result: a single CW resolves fringes for essentially NO pulsar.** The
  true-vs-best-wrong-fringe gap is **dlnL_a: median 0.0000, 84th pct 0.0002, max
  0.0028** across all 116 psr x 10 draws, versus **ln K_a: median 7.94 (range
  2.48-11.05)**. The gap is ~1000x below the trials threshold. Direct scan confirms
  it: for J1713 the adjacent fringe sits only 0.06 below truth and the global best
  wrong fringe is 0.004 below; for J0437, 0.0005. The chirp (pulsar-term frequency
  evolution over L/c) breaks the fringe degeneracy only at the ~1e-3 lnL level here
  -- far too weak to identify the fringe.
- **CLASSIFICATION (median over draws): class (i) = 0/116, class (ii) = 0/116
  (range 0-1), class (iii) = 116/116 (range 115-116).** Class-(i) [fringe
  identifiable AND marginal sigma_L < EM prior] is EMPTY. The 2 draws with one
  "identifiable" pulsar are K_a-small edge cases (few fringes in a narrow prior,
  ln K_a ~ 0), not physical chirp-breaking. Under the stricter dlnL>lnK+3 the count
  is identical (0 median).

  Representative K_a table (median over draws):

  | pulsar | K_a | ln K_a | dlnL_a | dL_pc | marg_pc | EM_pc | ident |
  |--------|-----|--------|--------|-------|---------|-------|-------|
  | J1713+0747 | 254 | 5.54 | 0.001 | 0.448 | 0.097 | 19  | no |
  | J0437-4715 |  12 | 2.48 | 0.001 | 0.489 | 0.249 | 1   | no |
  | J1909-3744 | 183 | 5.21 | 0.000 | 0.424 | 0.181 | 13  | no |
  | J1744-1134 | 154 | 5.04 | 0.000 | 0.361 | 0.199 | 9.3 | no |
  | J0030+0451 |  44 | 3.78 | 0.003 | 0.485 | 0.125 | 3.6 | no |

- **SNR cross-check (amendment 3), draw seed 2002.** Fringe-count arithmetic
  confirmed: J1713 K = 2*3*19 pc / 0.448 pc = 254 (matches). sigma_L scaling: the
  coherent SNR-limited law sigma_L = dL/(2*pi*SNR_pterm) holds to ~10% when SNR is
  the TIMING-MARGINALISED pulsar-term SNR (ratios sigma_pred/sigma_Fisher: J1713
  0.87, J1909 0.81, J0437 0.95, J1744 0.90). Using the RAW white-noise pulsar-term
  SNR overpredicts precision ~1.6x (ratio ~0.6) because the timing model (and, at the
  array level, the frozen GWB+RN GPs) absorb the low-frequency pulsar-term power --
  the same frozen-GP optimism D5 will bound. Per-pulsar pulsar-term SNRs are sub-unity
  (~0.24-0.69) even for the loud draw.
- **J0437-4715 bookkeeping (amendment 4).** IN the array. D2: marginal sigma_L
  0.249 pc < EM prior 1 pc (within-fringe curvature below prior, rank 39/116 -- 38
  pulsars have tighter within-fringe sigma, so it misses the D2 top-10). D3: NOT
  fringe-identifiable (dlnL 0.001 vs ln K 2.48) -> class (iii). NB the feather EM
  prior for J0437 is 1 pc, already LOOSER than its true VLBI distance (~0.16 pc); the
  whole array's feather priors are wide (median 235 pc), so "marg < EM" is an easy
  bar and still not a measurement without the fringe.

- **HEADLINE (goal "y" in real units, honestly): a single realistic CW yields ZERO
  genuine pulsar-distance measurements across a 116-pulsar array**, despite D2's
  116/116 within-fringe curvature below the EM prior. Within-fringe sigma_L without
  fringe identification is NOT a distance measurement. Fringe identification requires
  breaking the single-CW degeneracy -- i.e. MULTIPLE CWs with different dL per source
  (the D4 sweep) -- or an external distance prior tighter than dL (~0.4 pc), which is
  below even VLBI for all but the nearest pulsars.

**LANGUAGE RULE (enforced henceforth): no "N distances measured/improved" without the
class-(i) qualifier. D2's 116/116 is "within-fringe curvature below EM prior", a
necessary-not-sufficient condition; the measurement count is the class-(i) count = 0
at N_CW=1.**

**D4 — N_CW sweep on the real likelihood (2026-07-02).** Module
`CW_transition/stagec_d4.py`; results `stagec_d4_results.npz`. N_CW=1,2,4,8,16 equal
strain (h=-13.75 all sources), 10 Asimov zero-noise geometry draws each; plus one
population draw (3 loud + 13 faint, 10x) at N_CW=16. Uses the amortised zero-noise
Fisher (HVP-assembled full Hessian so N_CW=16 fits) + the D3 direct 512-eval fringe
scan. Fisher-validity flag: pulsar with marginal sigma_L > dL/4 is "fringe-unresolved"
(within-fringe width non-Gaussian) and excluded from sigma_L quotes.

| N_CW | class-i (med/min/max) | med marg sigma_L (pc) | marg/cond med(84,max) | dlnL_a med(84,max) | ln K_med |
|------|----------------------|-----------------------|-----------------------|--------------------|----------|
| 1  | 0/0/0   | 0.268 | 1.01 (1.0,13) | 0.000 (0.00,0.1)  | 7.86 |
| 2  | 0/0/0   | 0.184 | 1.01 (1.0,10) | 0.002 (0.01,4.2)  | 8.28 |
| 4  | 0/0/1   | 0.128 | 1.02 (1.1,12) | 0.029 (0.13,5.1)  | 8.49 |
| 8  | 0/0/2   | 0.093 | 1.02 (1.1,6)  | 0.246 (0.87,7.9)  | 8.65 |
| 16 | 5/2/10  | 0.067 | 1.04 (1.1,11) | 0.946 (2.97,28.4) | 8.79 |

(Fisher-valid count reported ~116/116 every draw -- **CORRECTED 2026-07-03: this was a
units bug** (`marg` [kpc] compared to `dL/4` [pc], 1000x too loose). True Fisher-valid
(marg_pc < dL_pc/4) medians are **29/16/28/43 of 116** at N_CW=1/4/8/16 (pop 71/116): the
majority of pulsars have marg sigma_L > dL/4, i.e. the within-fringe width is a large
fraction of a fringe (non-Gaussian) -- consistent with the sub-unity pulsar-term SNRs (D3).
The bug affected ONLY the reporting flag and the "median marg over valid pulsars" quote;
**no classification / certification count changes** (class-(i) and the fringe-ID test never
used `valid`). All 116/116 NSD is unaffected. Patched in stagec_d4.py:115.)

- **Class-i (goal "y") switches on GRADUALLY and TAIL-DRIVEN, not sharply.** Median
  class-i is 0 through N_CW=8 (max 1-2) and first becomes non-zero at **N_CW=16
  (median 5, range 2-10 of 116)**. The MEDIAN pulsar is never identifiable even at
  N_CW=16 (dlnL_a median 0.95 << ln K 8.8); it is the tail (max dlnL_a grows
  0.1 -> 4.2 -> 5.1 -> 7.9 -> 28.4) that crosses the ~8.8 trials threshold. So a few
  best-geometry pulsars become genuine measurements at N_CW=16 while most stay
  fringe-degenerate. The fringe-ID gap dlnL_a scales steeply in N_CW (median
  0.000 -> 0.946, ~x1000 over 1->16) -- roughly one order of magnitude per doubling --
  consistent with each extra source adding an independent fringe-breaking constraint.
- **REALISTIC POPULATION KILLS IT: 3 loud + 13 faint (10x) at N_CW=16 gives class-i =
  0** (vs median 5 for 16 equal). Only the loud sources carry fringe-breaking power,
  so identification behaves like N_CW ~ 3 -> back to the degenerate regime. Within-
  fringe sigma_L is actually BETTER (0.051 pc, the loud sources tighten the curvature)
  but ZERO distances are identifiable. This is the equal-strain-artifact lesson from
  the toy (Stage A.1: "mode-independence is an equal-amplitude artifact; a realistic
  luminosity function is set by the handful of loud sources") now in real units and
  with the fringe-ID layer: **the real SMBHB sky (few loud + many faint) sits in the
  fringe-degenerate regime -- a single-digit number of loud sources cannot phase up a
  116-pulsar array.**
- **marg/cond stays ~1 (median 1.01 -> 1.04); the conditional optimism does NOT open
  up over N_CW<=16.** The toy predicted marg/cond collapses at the confusion knee
  (N* ~ 50*T*Deltaf ~ hundreds-thousands); at N_CW<=16 the real likelihood is deep in
  the RESOLVABLE regime, so there is plenty of within-fringe information and
  conditional ~ marginal. The real likelihood therefore does NOT reproduce the toy's
  resolvable->confusion turnover in this range (we never approach the knee). What it
  ADDS beyond the toy Fisher is the DISCRETE fringe-identification requirement, which
  the toy (within-fringe info only) never saw -- and THAT, not confusion, is the
  binding constraint on real distance measurements at achievable N_CW.
- **Within-fringe sigma_L improves monotonically** 0.27 -> 0.067 pc (equal strain,
  N_CW 1->16) as sources accumulate distance information -- but this is
  necessary-not-sufficient (LANGUAGE RULE: not a measurement without class-i).

**D4 headline: genuine pulsar-distance measurements require MANY equal-strain CWs
(onset ~N_CW=16, and only for a best-geometry tail of ~5/116); a realistic spiky SMBHB
population (few loud sources) yields ZERO even at N_CW=16. The binding constraint at
achievable N_CW is fringe identification, not within-fringe precision or source
confusion.**

**D5 — frozen-GWB sensitivity (caveat 7.2) (2026-07-02).** Module
`CW_transition/stagec_d5.py`; results `stagec_d5_results.npz`. At N_CW=8 equal strain
(pre-onset, most sensitive), recomputed the D2/D3 pipeline with the frozen GWB
log-amplitude at truth and truth +/- 0.5 dex (rebuild the amortised Fisher per shift;
Asimov CW data is GWB-independent, only the frozen GP prior amplitude moves). Same 10
geometry draws (seeds 2000+d).

| frozen GWB log10_A | med marg sigma_L (pc) | med dlnL_a | class-i med(min,max) |
|--------------------|-----------------------|------------|----------------------|
| -15.10 (truth-0.5) | 0.0930 | 0.238 | (1, 0, 2) |
| -14.60 (truth)     | 0.0933 | 0.235 | (0, 0, 2) |
| -14.10 (truth+0.5) | 0.0953 | 0.214 | (0, 0, 2) |

- **Movement: marginal sigma_L x1.00 (-0.5 dex) to x1.02 (+0.5 dex); median dlnL_a
  x1.01 to x0.91; class-i unchanged (median 0-1, range 0-2).** A +/-0.5 dex error in
  the frozen GWB amplitude moves within-fringe precision by <=2% and the fringe-ID gap
  by <=9%, and does not change the (already ~zero at N_CW=8) class-i count. Louder
  frozen GWB (+0.5) mildly degrades (absorbs more low-frequency power overlapping the
  pulsar term); quieter (-0.5) is indistinguishable from truth.
- **Where realism effort should go (the decision D5 was meant to inform).** The
  frozen-GWB-amplitude optimism is SMALL: ~2-9%. By contrast the D3 SNR cross-check
  showed the timing-model + GP absorption already costs ~1.6x in pulsar-term SNR
  (~2.6x in information / dlnL, ~1.6x in sigma_L) -- an order of magnitude larger than
  the GWB-freeze sensitivity. And D4 showed the population LUMINOSITY FUNCTION is
  decisive (equal-strain N_CW=16 class-i median 5 -> realistic 3-loud+13-faint class-i
  0). **Conclusion: un-freezing / marginalising the GWB amplitude is NOT where the
  optimism lives for this prong; the dominant realism levers are (1) the source
  luminosity function and (2) real-noise vs Asimov. The frozen-at-truth GP is an
  acceptable approximation for the distance question at the +/-0.5 dex level.**

**D5 caveat scope:** this varies only the GWB log-AMPLITUDE at fixed spectral index and
fixed (frozen) per-pulsar red noise; a full treatment would also marginalise gamma and
the RN, but the amplitude is the leading GWB term and its <=9% effect bounds the concern.

**D5 add-on — nb-05 joint-dlnL blast-radius closure.** Module
`CW_transition/stagec_nb05_recheck.py`. Reran nb-05's joint-DeltalnL (116 psr, seed
1234, realistic draws, NOISY simulate, GWB NG15, 1 us WN) with the prominence bug
FIXED: best-wrong distance from a DIRECT FINE-GRID scan (up to 8000 pts, max beyond
0.5 dL, no find_peaks / no prominence floor) instead of the peak finder. (First tried
D3's fringe-CENTRE eval but nb-05 is noisy -> the best-wrong peak shifts off centre, so
centres under-find it; the fine grid is the faithful fix.)

| N_CW | original joint dlnL | corrected joint dlnL | corrected more negative |
|------|---------------------|----------------------|-------------------------|
| 1 | -60.83  | -65.24  | 7%  |
| 2 | -76.73  | -99.47  | 30% |
| 3 | -122.84 | -132.67 | 8%  |
| 4 | -137.05 | -159.57 | 16% |
| 5 | -177.20 | -208.04 | 17% |
| 6 | -186.63 | -230.45 | 23% |
| 7 | -195.59 | -253.40 | 30% |
| 8 | -174.76 | -254.23 | 45% |

(N_CW=9,10 not reached before the wall-clock cap; the 1-8 trend is monotonic and
decisive.)
- **The bug was OPTIMISTIC, as triaged.** The corrected joint dlnL is MORE NEGATIVE at
  every N_CW (7-45%): the peak finder left fringe-degenerate pulsars at truth,
  understating how well a wrong-fringe combination fits, i.e. OVERSTATING separability.
  The fix makes the wrong config fit better -> dlnL more negative.
- **nb-05's SIGN verdict is unchanged, its margin firmer.** Joint dlnL is strongly
  NEGATIVE for all N_CW=1-8 both ways: in nb-05's noisy realisation the truth distances
  are NOT separable from wrong-fringe combinations at these strains -- a wrong config
  beats truth by 60-250 lnL. The correction does not flip any conclusion; it removes an
  optimistic bias of order 10-45%. The bug bit LESS than in the Asimov D3 case because
  nb-05's noise itself creates spurious peaks above the 0.5 floor, partially masking the
  degeneracy.
- **Closure scope:** nb-05 (the key science nb) is now rechecked. 01_single is also
  suspect and 00/02 mixed at N_CW=1 (see triage), NOT rerun -- flagged in the doc for a
  decision; the corrected fine-grid scanner is the tool to use if they are revisited.

(D6 close-out synthesis is at the top of this Stage C section; Stage C COMPLETE.)

### Prong-2 close-out (P2-A..D, 2026-07-02)

**P2-A — DETECTABILITY vs RANGEABILITY (the field-facing headline figure).** Module
`CW_transition/stagec_p2a.py`; `stagec_p2a_results.npz`, `p2a_detect_vs_range.png`.
Same D4 draws (116 psr, equal strain h=-13.75, 10 Asimov draws) + the 3-loud+13-faint
population draw. DETECTABLE = # sources with Earth-term-only matched-filter SNR > 5
against the real marginalised noise (SNR_s^2 = (s|s) = 2[logL(theta_off,0) -
logL(theta_off, r_earth_s)] from the amortised likelihood with the CW model amplitude
killed). RANGEABLE = D4 class-i.

| N_CW | detectable SNR>5 (med) | (>3, >8) | rangeable class-i (med) |
|------|-----------------------|----------|-------------------------|
| 1  | 1/1   | (1, 0)  | 0 |
| 2  | 1/2   | (2, 1)  | 0 |
| 4  | 3/4   | (4, 2)  | 0 |
| 8  | 7/8   | (8, 3)  | 0 |
| 16 | 13/16 | (16, 7) | 5 |
| pop (3 loud+13 faint) | 3 | (7, 3) | 0 |

- **THE GAP (headline sentence):** the sky becomes VISIBLE almost immediately -- most
  sources are individually detectable from N_CW>=2 (SNR>3 catches ~all; SNR>5 catches
  most: 13/16 at N_CW=16) -- but it becomes RANGEABLE only at N_CW=16 (5/116 pulsars,
  and only the best-geometry tail), and NOT AT ALL for a realistic population (3
  detectable loud sources, 0 rangeable). **Detectability and rangeability are different
  transitions separated by ~an order of magnitude in N_CW; for the real SMBHB sky
  (few loud sources) the array is a detector, not a ranger.**
- **GATE (cross-check, one loud source log10_h=-13, f=20 nHz, mc=9.5):** Earth-term
  SNR_full (marginalised white+timing+RN+GWB) = 23.49; white matched-filter
  sqrt(sum (r/sigma)^2) = 36.80; ratio 0.638. This is NOT within 20% of the NAIVE white
  formula -- but 0.638 = 1/1.567 REPRODUCES THE D3 timing+GP absorption factor (~1.6x)
  to <1%. Two independent SNR computations (D3 pulsar-term, P2-A Earth-term) give the
  same absorption -> the machinery is validated; the naive white formula is the one
  that is 1.6x optimistic (it omits timing-model marginalisation). Reported as
  gate-satisfied-in-substance (code correct; the >20% offset is the known, precisely
  reproduced physical absorption, not an error).

**P2-B — coherence-axis physical grounding (analytic; v2 after a first-pass rejection).**
Module `CW_transition/stagec_p2b.py`; `p2b_coherence_grounding.png`,
`stagec_p2b_results.npz`. Coherence Y-axis unit = SNR*sigma_phi (Stage A.2: distance
info halves at SNR*sigma_phi = 1). Pulsar-term lag tau_p = (L/c)(1-cos mu) ~ 3.3 kyr for
L=1 kpc.
- **FIRST-PASS ERROR (stated for the record):** v1 computed the pulsar-side wander as
  sigma_phi = 2 pi f sigma_res (treating red/DM noise as a time-base SHIFT of the pulsar
  term) = ~6e-15 rad -- WRONG projection. Additive noise perturbs the pulsar-term PHASE
  by the NOISE-TO-SIGNAL ratio sigma_phi = sigma_res / A_CW, A_CW = h/(2 pi f) the CW
  pulsar-term timing amplitude. The missing division by A_CW is a factor 1/h = 5.6e13
  (~13.8 orders). Corrected sigma_phi = sigma_res*2 pi f/h = 1/SNR_pterm ~ O(1) rad.
- **Item 2 CORRECTED (pulsar-side, real feather white noise sigma_TOA*sqrt(2/N) in the
  10 nHz bin; A_CW = per-pulsar pulsar-term RMS; h=-13.75):** chain per named pulsar --
  J1713+0747 sigma_res 40.5 ns / A_CW 32.1 ns -> sigma_phi 1.26 (SNR_pterm 0.79);
  J0437-4715 93.3/126.2 -> 0.74 (1.35); J1909-3744 51.2/6.7 -> 7.64 (0.13). **Array
  sigma_phi median 2.13 rad (16/84 = 0.66/8.14), SNR_pterm median 0.47.** So sigma_phi ~
  O(1) rad = 1/SNR_pterm: this IS the measurement/SNR floor, already inside the
  marginalised-noise Fisher (D2-D5), sitting at the knee BY CONSTRUCTION
  (SNR_pterm*sigma_phi = 1). Red/DM noise is an SNR-axis effect, NOT an independent extra
  coherence term. (Consistent with D3's sub-unity per-pulsar pulsar-term SNR.)
- **Item 1b CORRECTED (GW chirp accumulated phase Delta_phi = pi*fdot*tau_p^2, tau_p=3
  kyr, across the (Mc,f) plane):** 0.02 rad (Mc=1e8, 3 nHz) to 2.9e4 rad (Mc=1e9.5, 30
  nHz); e.g. Mc=1e9 @ 10 nHz -> 76 rad. LARGE but DETERMINISTIC and MODELLED
  (evolve=True) -> recoverable info, OFF the coherence axis (NOT plotted as a band).
- **Item 1a CORRECTED (required stochastic df/f to cross the knee, as an SNR family):**
  df/f = 1/(2 pi f tau_p SNR); crossings at df/f = 3.4e-5 (SNR=5), 8.4e-6 (SNR=20),
  1.7e-6 (SNR=100). Magnitude of real environmental df/f left explicitly OPEN.
- **Figure (2 panels, clipped):** (1) SNR*sigma_phi vs df/f family with knee crossings;
  (2) the 116 per-pulsar sigma_phi histogram (median 2.13 rad, knee at 1).
- **VERDICT:** the coherence axis is NOT a SEPARATE limit beyond SNR on the pulsar side
  -- red/DM noise gives sigma_phi = 1/SNR_pterm, the measurement floor already in D2-D5.
  The genuinely independent coherence question is SOURCE-side stochastic df/f, which
  crosses the knee at df/f ~ 1e-6..1e-5 over kyr lags; whether real SMBHB environments
  reach that (gas/stellar coupling over kyr) is an astrophysics question this project's
  methods cannot answer -> HANDOFF (Taylor/Farr).

**P2-C — array-boost scaling law (toy).** Module `CW_transition/stagec_p2c.py`
(+ 200-psr rerun); `p2c_array_boost.png`, `stagec_p2c_results.npz`. Reran the Stage A
confusion-knee measurement (fiducial T=15 yr, band 3-20 nHz -> N*=8.05,
fixed_persource, isotropic sky, 30 seeds, N grid 1-1500) at N_psr = 15, 30, 60, 116,
200 (200 = synthetic isotropic positions, same protocol).

| N_psr | knee/N* |
|-------|---------|
| 15  | 6.47  |
| 30  | 12.92 |
| 60  | 26.73 |
| 116 | 52.75 |
| 200 | 91.45 |

- **GATE: N_psr=116 -> knee/N* = 52.75, reproduces the Stage A ~52 (PASS).**
- **Power law: knee/N* = 0.40 * N_psr^1.03.** Low-end slope (15-60) = 1.02, high-end
  (60-200) = 1.02 -> **NO saturation; the array boost is LINEAR in N_psr through 200.**
- **Forecasting law (converts the Stage A constant into a prediction): the confusion
  knee sits at N_knee ~ 0.40 * N_psr * N* = 0.40 * N_psr * T * Delta_f.** More pulsars
  push the confusion wall out PROPORTIONALLY (each doubling of the array doubles the
  resolvable-source count before marg/cond collapses); no diminishing returns up to 200
  pulsars. The Stage A "~52x for 116" is thus the N_psr=116 value of a linear law, not a
  universal constant.

**P2-D — loose ends (two items).** Module `CW_transition/stagec_p2d.py`;
`stagec_p2d_item1.npz`, `stagec_p2d_item2.npz`.
- **Item 1 (strain reconciliation, closes the D6 "were the old runs louder"):**
  per-source OPTIMAL SNR (full Earth+pulsar CW, marginalised noise), median over 80
  sources: h=-13.75 (D4 equal strain) = **10.7**; h=-13.0 (nb-03 loud) = 59.9;
  **h=-12.0 (CW_node_sampling optimizer, MK5/MK7) = 599.5**. Ratios: optimizer/D4 =
  **56.2x** (= 10^1.75 exactly -> SNR scales linearly with strain), nb-03/D4 = 5.6x.
  **The old prong-1 optimizer's 20/20 distance recoveries ran at SNR ~ 600 per source
  (h=-12, frequently noise-free) -- ~56x louder than the realistic D4/nb-05 regime
  (SNR ~ 11). The earlier "distances recoverable" successes were a far easier regime;
  this reconciles them with Stage C's realistic-population zero.**
- **Item 2 (nb-01 recheck, corrected fine-grid scanner):** nb 01_single, single CW,
  scenario="single", NOISY simulate, seed 1234, three strains:

  | log10_h | original joint dlnL | corrected joint dlnL |
  |---------|---------------------|----------------------|
  | -14.0 | -56.50  | -68.87  |
  | -13.5 | -46.64  | -101.79 |
  | -13.0 | **+296.09** | **-108.35** |

  Optimistic at every strain, and **at h=-13 it FLIPS THE SIGN**: the buggy peak finder
  reports joint dlnL = +296 (reads as "truth preferred -> distances SEPARABLE"), the
  corrected scanner gives -108 (a wrong-fringe combo BEATS truth -> NOT separable). The
  bug inverted the verdict for the loud single-CW case; corrected confirms the Stage C
  physics (single CW fringe-degenerate, not separable even at h=-13; cf D3). 00/02
  remain flagged-not-rerun per the D6 decision.

**>>> PRONG-2 CLOSURE (2026-07-02).**
Prong 2's INFORMATION-THEORETIC content is COMPLETE under BOTH interpretations:
- *distance-information transition* (Stages A/A.1/A.2/A.3 toy + Stage C D0-D6 real
  likelihood): within-fringe marginal sigma_L, the fringe-identification split, the
  N_CW sweep, GWB-freeze sensitivity, and the array-boost scaling law
  (knee/N* = 0.40 N_psr^1.03, P2-C);
- *detectability vs rangeability* (P2-A): the sky is visible almost immediately but
  rangeable only at N_CW~16 (best-geometry tail), never for a realistic population.
Both converge: for the real SMBHB sky (few loud + many faint) a 116-pulsar array is a
DETECTOR, not a RANGER; genuine distance measurement is gated by fringe identification,
which a realistic population does not supply.

ONE OPEN QUESTION, formulated and assigned OUTWARD: whether real SMBHB environmental
evolution produces a STOCHASTIC df/f above the P2-B thresholds (~1e-6..1e-5 over kyr
lags) -- an astrophysics question (gas/stellar coupling) NOT answerable by this
project's methods -> HANDOFF (Taylor/Farr). The pulsar-side coherence component IS
resolved in-house (P2-B item 2: red/DM = the 1/SNR measurement floor, in D2-D5 already).

Section-7 caveats resolved/bounded/assigned: 7.1 conditional-vs-marginal RESOLVED
(honest marginal; the old metric shown optimistic & sign-flipping, corrected scanner is
the standing tool); 7.2 frozen GP BOUNDED (D5: <=9%); 7.5 homogeneous WN SUPERSEDED
(real feather noise used); old-loud-run strain RESOLVED (P2-D item 1: 56x). 7.3/7.4 are
prong-3 sampler items (out of scope here).

**PRONG 2 IS CLOSED AS A COMPUTATION AND OPEN AS ONE HANDOFF** (the environmental-df/f
question to Taylor/Farr).

**>>> TOY PHASE CLOSED (historical anchor).** Stages A (confusion), A.1 (robustness/readouts), A.2 (coherence),
A.3 (interaction map) establish the information-theoretic skeleton. **Next = Stage C:**
port marginal distance info onto the real discovery likelihood
(`build_fast_scan_likelihood`, autodiff Hessian marginalising CW params), real
heterogeneous noise + GWB GP (drop the frozen-at-truth optimism), and report **absolute
sigma_L in pc** (the toy's schematic ζ=h/2πf normalisation makes only shapes/ratios
trustworthy, never absolute precision). Cross-check the marginal-info threshold against
the conditional joint-Δln L metric. Open risk (from §9): discotech's GPU jax may need a
pinned discovery venv.

### Anchor Census — real EM priors vs the feather priors (COMPLETE, 2026-07-03)

**Question.** Stage C used the array's stored "feather" distance priors (median sigma_EM
~235 pc) and certified ZERO distances at realistic populations, because a distance counts
as measured only if its fringe is IDENTIFIED: joint dlnL(true vs best-wrong fringe in the
+/-3 sigma_EM window) must beat the trials factor ln K, K = 6 sigma_EM/dL. The feather
priors are far wider than the best current EM measurements. HYPOTHESIS: with REAL priors
(VLBI / timing parallax / orbital-Pbdot), a few "anchor" pulsars may have K ~ 1-10 -- cheap
to certify, able to phase-lock the array (cf. arXiv:2603.28897, which ASSUMES such anchors).
This census audits whether they exist on the real 116-pulsar array. It does NOT build the
sequential estimator (the next experiment). Env: jug-gpu. Priors from the canonical
`/home/mattm/soft/scripts/MISC/get_distance.py` (live Cornell parallax catalogue, cached).

**A0 -- real priors (canonical script).** Imported get_distance.py (fetched 156-pulsar
Cornell catalogue), built PulsarData per array pulsar from MPTA-DR3 + NG-CHIME par files
(3 B->J aliases: B1855+09=J1857+0943, B1937+21=J1939+2134, B1953+29=J1955+2908), wrote the
canonical `best_distances.txt` (116 rows, FROZEN in CW_transition/). Parallax-class =
{Catalog_PX, VLBI, Timing_PX}; everything else KEEPS the feather prior (flagged; NO DM
distance ever used as a Gaussian sigma_EM).
- **Coverage of 116:** Catalog_PX 30, VLBI(legacy) 1, Timing_PX 38 -> **69 parallax-class**;
  47 feather-kept (44 wide DM-defaults, 3 no-match). 80/116 have no Cornell entry (mostly
  southern MPTA); every unmatched pulsar reported, none silently dropped.
- **GATE J0437-4715 PASS:** real distance 0.1563 kpc (~0.157 expected); Catalog_PX
  sigma_EM 1.32 pc, Deller+2008 VLBI PX=6.396+/-0.054 mas.
- **Caveat:** get_distance prefers VLBI-technique over a more precise timing PX, so for
  several bright MSPs the "real" prior is WIDER than the feather (J1713+0747: real 60.9 pc
  vs feather 19 pc). The tighten-hypothesis holds strongly for a few (J2222-0137 R=0.018,
  J2241 0.15) but NOT universally.

**A1 -- dual-prior K table (fiducial = D3 single-CW config, dL median over 10 draws).**
Reported K under two prior choices side by side:
- K_canonical: sigma_EM exactly as best_distances.txt delivers (collaboration-standard).
- K_optimal: sigma_EM = MIN sigma among geometric-quality methods {all Cornell techniques,
  VLBI, Timing_PX, secure PBDOT_Shk}. Secure PBDOT = Pbdot SNR>10 AND PM measured AND masses
  known (GR subtracted) -- **NO pulsar qualifies** (J0437 fails on absent M1), so PBDOT never
  contributes. Inverse-variance combination reported informationally; K uses min-sigma.

  | pulsar | K_opt | K_canon | K_feath | sig_opt pc | dL pc | sig_opt method |
  |--------|-------|---------|---------|-----------|-------|----------------|
  | J0437-4715 | 11.88 | 16.21 | 12.28 | 0.97 | 0.489 | Catalog Timing (Reardon24) |
  | J2222-0137 | 13.56 | 13.56 | 746 | 0.97 | 0.431 | Catalog VLBI (Ding24) |
  | J1744-1134 | 62.0 | 288.6 | 154.6 | 3.73 | 0.361 | Timing_PX |
  | J0030+0451 | 66.9 | 66.9 | 44.5 | 5.41 | 0.485 | Catalog VLBI |

- **Counts (of 116): K<=3 = 0 ; K<=10 = 0 ; K<=30 = 2 (J0437, J2222)** -- under BOTH
  canonical and optimal columns. Smallest K anywhere = **J0437 K_opt = 11.88**.
  Named K<=10 list is EMPTY under every prior.
- **HAND-CHECKS.** J0437: K_opt = 6*0.9675/0.4887 = 11.88; K_canon = 6*1.320/0.4887 = 16.21;
  K_feath = 6*1.000/0.4887 = 12.28. J1713 REVERSAL (not averaged away): K_opt = 6*40.0/0.4479
  = 536, K_canon 816, K_feath 255 -- real sigma (40-61 pc) is 2-3x WIDER than feather (19 pc)
  because get_distance's VLBI-preference is looser than the timing PX baked into the feather.
- **>>> GATE OUTCOME (J0437 K_optimal ~3): the canonical script CANNOT reach it.** The
  script's tightest J0437 is Cornell Reardon24 timing 0.97 pc -> K = 11.88. The famous
  156.79 +/- 0.25 pc (Reardon+2016) is a JOINT-timing composite (annual+orbital parallax +
  Pbdot + kinematics), not in the catalogue and not reproducible by any single method (a
  perfect Pbdot gives ~2.8 pc, SNR-limited). At sigma=0.25 pc, K = 6*0.25/0.489 = 3.07 (the
  gate arithmetic is right; only the value is unavailable from the script). **Decision (with
  Matt): inject published composites for A2's best case, run A2 under BOTH literature and
  script-optimal priors.** Literature injections (WITH refs): J0437 156.79+/-0.25 pc
  (Reardon+2016), J1909-3744 1152+/-3 pc (Reardon+2021 Shklovskii). J2222 catalog Ding24
  0.97 pc is already the best. Net: literature injection yields exactly ONE K<=3 pulsar
  (J0437, K=3.07).

**A2 -- reclassification with real priors (controlled: same seeds 1000/2000+d/3000; truth at
feather distance, only sigma_EM width changes).** Reran D3 (N_CW=1) + D4 (N_CW=4,8,16 equal
strain + 3-loud/13-faint population), 10 draws, THREE prior columns (feather, script-optimal,
literature) classified from ONE scan/draw (sigma_L and per-fringe lnL are prior-independent).
Optimizations: vectorized fringe scan across all pulsars x all draws (GATE: vectorized ==
per-pulsar loop bit-identical, max|.| = 0.0); JAX persistent compilation cache; Hessian once
per draw shared across the 3 columns. Four classifications reported SEPARATELY.

  ANCHOR-CERTIFIED flat (dlnL > ln K_counted; K_counted = true wrong-fringe centres in
  window), med(min,max):

  | N_CW | feather | script-real | literature |
  |------|---------|-------------|------------|
  | 1  | 0(0,1) | 0(0,1) | 0(0,1) |
  | 4  | 0(0,1) | 0(0,1) | 0(0,1) |
  | 8  | 0(0,1) | 0(0,1) | **1(0,3)** |
  | 16 | 5(2,10) | 5(2,11) | **6(3,11)** |

  CERTIFIED-BAYES (P_true>0.9, softmax[lnL + ln Gaussian-prior] over candidate fringes --
  the honest small-K criterion; prior tails penalise wrong fringes, ~1.9 nats for J0437's
  adjacent fringe at 0.25 pc):

  | N_CW | feather | script-real | literature |
  |------|---------|-------------|------------|
  | 1  | 0(0,1) | 0(0,1) | **1(0,1)** |
  | 4  | 0(0,1) | 0(0,1) | 0(0,1) |
  | 8  | 2(0,5) | 2(0,5) | **3(1,5)** |
  | 16 | 18(11,21) | 18(12,21) | 18(12,21) |

- **PRIOR-CERTIFIED (K<1, prior pins the fringe with no CW): 0 at every N_CW, every prior.**
  No pulsar's EM prior is tighter than one fringe spacing; even J0437 literature (0.25 pc)
  has K_counted = 2-3. **There is NO prior-alone anchor on this array.**
- **IMPROVED (flat-ident AND marg sigma_L < sigma_EM real):** tracks flat-ident closely
  (feather/script 5, lit 6 at N_CW=16); harder for tight-prior anchors as expected.
- **POPULATION (3 loud+13 faint) N_CW=16 -- the headline cell, fully audited:** flat = **0
  (all priors, incl. literature)**; bayes = **3** (J0711-6830, J1713+0747, J1909-3744).

  | pulsar | P_true(lit) | dL pc | sigma_lit pc | K_counted | dlnL | cause tag |
  |--------|-------------|-------|--------------|-----------|------|-----------|
  | J0711-6830 | 0.953 | 0.228 | 21.2 | 558 | 3.52 | DATA-driven (dL_need 51 >> dL) |
  | J1713+0747 | 0.989 | 0.190 | 40.0 | 1264 | 5.25 | DATA-driven (dL_need 96 >> dL) |
  | J1909-3744 | 0.998 | 0.243 | 3.0 | 72 | 4.02 | DATA-driven (dL_need 7.2 >> dL) |

  All 3 are **DATA-DRIVEN** (median geometry, the CW dlnL from the 3 loud sources breaks the
  fringe), NOT geometry-wide, NOT prior-tight. **They certify under ALL THREE prior columns
  including feather** (P_true 0.95/0.995/0.97) -- so the population bayes-3 is a
  LOUD-SOURCE effect, prior-INDEPENDENT, NOT a composite-prior payoff. Sampling check: J1909
  (sigma_lit 3 pc -> K=72 < 512) is FULLY sampled and robust; J1713/J0711 have wide windows
  (real sigma 19-40 pc -> K 558-1264 > 512 sampled) so the softmax denominator is
  under-sampled -- but the near-fringe competitors (which dominate) are captured and the far
  fringes carry large dlnL (loud-source-suppressed), so P_true is a safe over/robust estimate,
  not inflated. NB: none of the 3 has a tight real prior; the population certification does
  NOT come from the narrow composite windows -- it comes from the loud sources.
- **CAUSE-TAG AUDIT (required; per-draw dL vs sigma).** Every Bayesian certification is
  tagged prior-tight/geometry-wide vs data-driven. A pulsar can reach P_true>0.9 with NO CW
  help if its draw has dL large enough that <2 wrong fringes fall in the +/-3 sigma window
  (prior-only bound: at sigma, dL_need = sigma*sqrt(2*2.89) = 2.40*sigma for P_true=0.9 at
  N_CW=1). This is the "geometry-wide" artifact (antipodal pulsar, 1-cos mu -> 0 -> dL huge,
  but the pulsar-term SNR is ~0 so it carries no real information).
- **J0437 onset triple, CORRECTED (bar, prior-column, per-draw fraction of 10):**

  | prior | flat-ident 1/4/8/16 | bayes 1/4/8/16 | cause of the low-N certs |
  |-------|---------------------|----------------|--------------------------|
  | literature (0.25pc) | 3/3/5/7 | 4/2/5/6 | N_CW=1,4 GEOMETRY-WIDE (all bayes-draws dL>=0.69pc > dL_need 0.60pc; incl. a dL=12.3pc antipodal draw) |
  | script (0.97pc) | 1/1/0/2 | 1/1/0/3 | N_CW=1 carried by the single dL=12.3pc antipodal draw |
  | feather (1.0pc) | 1/1/0/2 | 1/1/0/3 | same dL=12.3pc antipodal draw |

  The median-geometry dL=0.489pc with sigma=0.25pc gives prior-only P_true=**0.77 < 0.9**, so
  J0437 CANNOT certify on median geometry at N_CW=1 -- the apparent N_CW=1 lit certs are all
  favorable-geometry (large-dL) draws. **J0437's first DATA-DRIVEN (median-geometry,
  CW-carried) certification is at N_CW=8** (bayes 5/10, dL~0.27-0.38pc, dlnL~2.7-3.5 > ln K),
  reaching 6-7/10 at N_CW=16. So the honest onset is (bar=bayes-data-driven, prior=literature,
  N_CW=8); the feather/script onset is N_CW=16 (median-geometry). The composite prior advances
  the DATA-DRIVEN onset one octave (8 vs 16) and adds geometry-wide certs at N_CW=1-4 that are
  NOT anchor physics.
- **Flat vs Bayesian, same draws:** under the conservative FLAT lnK bar real priors barely
  move certification (literature median +1 pulsar at N_CW=8 and N_CW=16). The honest BAYESIAN
  criterion certifies ~3x more (~18 vs ~6 at N_CW=16) but at N_CW=16 the CW likelihood
  dominates and all three priors converge (18/18/18) -- real priors help only where the CW is
  weak (low N_CW), and there mostly via geometry-wide draws for J0437.

**>>> ANCHOR CENSUS VERDICT (2026-07-03, cause-tag-audited).**
(i) **The strict anchor set (K<=1, certifiable by the EM prior alone) is EMPTY** -- zero
pulsars, under every prior including the hand-injected published composites. The tightest,
J0437 with the Reardon+2016 composite (0.25 pc), still has K_counted = 2-3 at median geometry
(prior-only P_true = 0.77 < 0.9); only ONE pulsar reaches K<=3 at all (J0437, K=3.07 at the
fiducial), and none lands in the K~1-10 "cheap anchor" band. **The prior-alone-anchor
assumption of arXiv:2603.28897 FAILS on the real array** -- there is no pulsar the EM prior
can phase-lock without help from the CW itself.
(ii) **Under the conservative FLAT (ln K) bar the realistic population certifies ZERO, at
every prior including composites** -- this is the honest conservative floor and it stands
unchanged from Stage C. Equal-strain flat onset is N_CW=16 (feather median 5, literature 6).
(iii) **Under the honest BAYESIAN bar** (posterior over fringes with the Gaussian prior tail),
the realistic 3-loud+13-faint population certifies **3 pulsars -- J0711-6830, J1713+0747,
J1909-3744 -- all DATA-DRIVEN by the loud sources and prior-INDEPENDENT (they also certify
under the feather prior)**, and equal-strain N_CW=16 certifies ~18 (prior-independent at high
N_CW). J0437's first genuine (data-driven, median-geometry) certification is at N_CW=8 under
the literature prior (bayes 5/10) vs N_CW=16 under feather; its apparent N_CW=1-4 certs are
geometry-wide artifacts (favorable large-dL draws), not anchor physics.
**Bottom line for the sequential estimator.** The array has NO prior-alone anchors; the
composite prior buys at most J0437 as a marginal, geometry-dependent, N_CW>=8 anchor. The
real seed set is the DATA-DRIVEN Bayesian certifications -- the ~3 loud-source-broken pulsars
that survive a realistic population (J0711/J1713/J1909) and the ~18 at equal-strain N_CW=16 --
NOT the EM-prior anchors. A sequential scheme must treat fringe identification as the binding
constraint and bootstrap from loud-source-broken pulsars, subject to the cause-tag audit
(exclude geometry-wide draws) and the wide-window softmax-sampling caveat (K>512 pulsars
undersampled; tighten if any becomes load-bearing).

**STATUS UPDATE (2026-07-06, Track B P1).** The census "3 data-driven certified pulsars"
(J0711/J1713/J1909) and the ~18 at equal-strain N_CW=16 are a CONDITIONAL-on-truth ceiling
(computed at exact truth source params; B0.2 showed the certification tolerance on those params
is ~1e-4 scaled ~ 0.006 deg sky, set by the 2*pi*L/dL lever arm). ACHIEVABILITY from a cold
start is the Track B question. **P1 (2026-07-06) established EXISTENCE:** a unique joint-lnL /
registration needle sits exactly at truth (the combs co-register ONLY at truth -> the joint
source+fringe optimum is unique and correct), so the ceiling is achievable IN PRINCIPLE.
**REACHABILITY is the open sub-question (P2):** the profiled(HARD-snap) joint-lnL is a smooth
funnel only > ~2 deg, then a near-flat noisy plateau to ~0.05 deg, then DEAD FLAT + a 0.006 deg
delta cusp -- a FLAT GAP between the F-stat basin and the needle that a local search cannot
cross. P2 tests whether the SOFT (marginal, logsumexp) surface has gradient across that gap; if
not, the ceilings are conditional artifacts UNREACHABLE by cold-start local search (confirming +
geometrically explaining B0.2). Achievability verdict = P2's result (pending).

**Canonical priors file going forward: `CW_transition/best_distances.txt` (frozen, A0), with
the two published composite injections (J0437 Reardon+2016, J1909 Reardon+2021) recorded in
`stagec_anchor_a2.py::LIT_INJECT`. Per-pulsar arrays: `anchor_a0_priors.npz` (A0),
`anchor_a1_Ktable.npz` (A1 K table, both columns), `anchor_a2_results.npz` (A2 counts).**

**Files (all in CW_transition/):** `stagec_anchor_a0.py`, `best_distances.txt` (frozen),
`anchor_a0_priors.npz` (A0); `stagec_anchor_a1.py`, `anchor_a1_Ktable.npz` (A1);
`stagec_anchor_a2.py`, `anchor_a2_results.npz` + `anchor_a2_diag.npz` (A2, the latter holds
per-draw per-pulsar dL/marg/sigma/P_true/dlnL/ident/bayes for the cause-tag audit);
`stagec_d4.py` (valid units-bug patch). Gates: J0437 A0 0.1563 kpc; A1
hand-checks (J0437 K arithmetic, J1713 reversal); A2 vectorized==loop 0.0, feather N_CW=1
reproduces D3 (certified 0, dlnL ~1e-3). Env: jug-gpu; compilation cache at
`~/.cache/jax_stagec_cache`. Side-finding: D4's Fisher-`valid` flag had a kpc-vs-pc units bug
(reported 116/116; true N_CW=1 valid ~41/116) -- affects only the within-fringe-Gaussian flag,
not certification.

### Track A — anisotropic-confusion distance information in the cross-pulsar covariance (F0–F2, 2026-07-03)

**Question (Hogg's contention / the confusion-regime referendum).** In the confusion regime the
GWB is finite-population shot noise, anisotropic at low multipoles. Does that sky *structure*
carry pulsar-distance information in the cross-pulsar covariance — recoverable without resolving
sources — or does Farr's isotropic zero-information theorem effectively hold (channel empty or
N-fragile)? The null (Farr holds) is the default; verdict language is kept no softer than the
numbers support.

**Method.** Frequency-domain cross-pulsar covariance including pulsar terms with finite distances,
C_ab(f) = ∫ P(f,n̂)[Fp_aFp_b+Fx_aFx_b](1−e^{−2πifτ_a})(1−e^{+2πifτ_b}) dΩ/4π, τ_p=(L_p/c)(1−cos μ_p).
Per-pulsar distance Fisher under the **Whittle likelihood** (amendment 2): I_LpLp = ½ Σ_k
Tr[C_k^{−1} (dC_k/dL_p) C_k^{−1} (dC_k/dL_p)] summed over the discrete Fourier bins f_k=k/T only
(T=15 yr, f_k=2.11 nHz·k), C_k = S_gw(f_k)·C_geo(f_k) + diag(white). dC/dL_p is rank-2
(e_p r_pᵀ + r_p* e_pᵀ); a closed-form trace (validated vs dense Tr[BdCBdC] to 1e-12) gives the
per-pulsar Fisher in O(n²). **Real 116-pulsar array**: real sky positions + real heterogeneous
white noise (median TOA err 74 ns–5 µs, med 464 ns → PSD 2σ²Δt, weekly). Fiducial GWB
h_c(f)=A(f/f_yr)^{−2/3}, **A=2e-15**, S_gw=h_c²/(12π²f³). Sky drawn as N equal-strain isotropic
point sources → empirical real c_lm to l=10 → band-limited P, normalised to **unit monopole so
total power is N-independent** (ties anisotropy to shot noise, c_lm∝1/√N; verified 3.836→1.211
for N=100→1000, ratio 3.17≈√10). ≥10 sky realisations per N. Sky integral: Gauss–Legendre×uniform
quadrature, separable real-Ylm basis; **complex64 GEMM with fp64 outer accumulation** (the
between-chunk oscillatory cancellation is carried by the fp64 accumulator; the RTX 4090's 1/64
fp64 rate made pure-fp64 at the mandated resolution ~hours/run). Files:
`CW_transition/trackA_covariance.py` (F0), `trackA_fisher.py` (F1/F2, modes validate|gates|f1|f2),
`trackA_realarray.npz`, `trackA_f1_results.npz`, `trackA_f2_results.npz`, `trackA_f2_scaling.png`,
`trackA_gates.npz`.

**Amendments implemented.** (1) baseline subtraction: headline is the anisotropic ENHANCEMENT
I_LL_aniso(N)=I_full(shot-noise sky,N)−I_iso(matched total power); iso row reported explicitly,
never asserted zero. (2) Whittle bins as above. (3) aliasing n_mu≥10·ftau0 (see gates). (4)
auto-term sanity (below). (5) FRINGE CAVEAT: every σ_L is WITHIN-FRINGE covariance information —
C oscillates in L_p with the same fringe period as the resolved-source likelihood; the fringes
moved into the covariance, they did not disappear. Fringe identification/certification applies to
this channel too.

**Gates.**
- **Trace-formula vs dense** 1e-12; **U-factorisation C vs F0 covariance()** 3e-16 (fp64) / 3e-5
  (production complex64); C Hermitian 1e-9. Real-Ylm orthonormality 1.6e-14; iso p≡1 exactly.
- **F0 spec correction stands** (amendment 1 origin): the isotropic cross-pulsar distance
  derivative is NOT zero — it is the finite-distance HD correction, cross/auto = 4.5e-2 / 9.6e-3 /
  2.85e-3 at ftau0=10/30/100, suppressed ~ftau0^{−1.2}. This is the baseline that is subtracted.
- **G1 — enhancement exactly zero off-anisotropy (GATE):** a shot-noise sky with all l≥1 c_lm
  zeroed reproduces the iso baseline to **0.00e+00** (identical code path). The enhancement can
  appear ONLY when c_lm≠0. PASS.
- **Aliasing (amendment 3).** F0 already shows the pulsar-term oscillation needs n_mu≳4·ftau0
  (at ftau0=100, n_mu=2·ftau0 gives the wrong cross 0.925; 4·ftau0 converges, drift ~1e-14 auto /
  1e-10 cross). Production n_mu=16384: **bin1 RESOLVED at 10.4·ftau0** (ftau0max=1411);
  bins 2/3/4 CAPPED at 5.8/3.9/2.9·ftau0 **for the farthest pulsar** (L=6.49 kpc; would need
  n_mu=28k/42k/56k). Crucially the cap bites only the ~handful of distant pulsars: a **typical**
  L≈1.3 kpc pulsar has bin4 ftau0≈1105, i.e. n_mu=16384 is ~15·ftau0 — fully resolved. Doubling
  gate on bin4 for the worst pulsar (J1903+0327, L=6.49, 2.9→5.8·ftau0, n_mu 16384→32768):
  covariance entries converge (**rel Δ|C_auto|=1.7e-5, Δ|C_cross|=7.4e-5**) but its per-pulsar
  **Fisher row moves ~100%** — so the capped-bin Fisher is unreliable *for the few L≳4 kpc pulsars*
  and reliable for the majority. This bounds the absolute-I_LL / σ_L uncertainty to the
  far-pulsar tail; the α-scaling and verdict are unaffected (the enhancement's resolution artefact
  is common to iso and shot-noise skies and cancels in I_full−I_iso, and its N-dependence is set
  by c_lm∝1/√N, not resolution).
- **Whittle observability table (amendment 2), iso sky.** Per-bin auto dC_aa/dL_a tracks the F0
  (π/3)·ftau0 law with constant ratio **0.787** across k=1..6 (i.e. auto grows ∝ f), while
  cross_max/auto falls ~1/ftau0:

  | k | f [nHz] | ftau0(medL) | auto_med | (π/3)ftau0 | ratio | cross_max | cross/auto |
  |---|---|---|---|---|---|---|---|
  | 1 | 2.11 | 276 | 2.28e2 | 2.89e2 | 0.787 | 2.40 | 1.05e-2 |
  | 2 | 4.23 | 552 | 4.55e2 | 5.78e2 | 0.787 | 2.61 | 5.73e-3 |
  | 3 | 6.34 | 828 | 6.83e2 | 8.68e2 | 0.787 | 2.20 | 3.22e-3 |
  | 4 | 8.45 | 1105 | 9.11e2 | 1.16e3 | 0.787 | 2.13 | 2.34e-3 |

  Observability correction: the Fisher SUMS these discrete bins; a smooth-f integral over the band
  would misestimate the information because the pulsar-term factor oscillates hundreds of cycles
  between adjacent Whittle bins (period 1/τ ≈ 1e-11 Hz ≪ 2.1 nHz spacing).
- **Auto-term sanity (amendment 4).** dC_aa/dL_a = (π/3)·ftau0 grows unboundedly with ftau0 (=f·L),
  but in the Fisher it enters as amp²(π/3 ftau0)²/(amp·C_aa+N)²; where GW power dominates the
  numerator and denominator both scale with S_gw(f), and where white noise dominates the ratio is
  suppressed. Net: the per-bin auto contribution is bounded, and the growth is genuine
  WITHIN-FRINGE curvature (amendment 5), not free distance information. Consequence: iso Fisher is
  **NOT bin-1-dominated** — the auto ∝ f growth offsets the red GW spectrum, so info splits
  **47.6% bin1 / 52.4% bins2–4** (measured at full resolution). This is why the capped-bin
  resolution matters and is checked by G3.

**F1 — I_LL per pulsar (kpc⁻², WITHIN-FRINGE).**
- (i) **ISO baseline**: median 13.7, max 119.2 (the finite-distance HD + auto covariance info that
  exists even for a perfectly isotropic background).
- (iii) **Maximal-anisotropy controls** (single l=2 / l=4 mode at unit monopole, min(P)=0 envelope):
  enhancement median 7.2 / 6.1, max 232 / 154 — the upper envelope.
- (ii) **Shot-noise skies** (enhancement I_full−I_iso, psr-median [best pulsar]):
  N=100 → 14.95 [J0437 168.6, 16–84 = 90–307]; N=300 → 5.43 [63.8, 7.5–117];
  N=1000 → 0.84 [13.0, **−2.1 to 28.8**]; N=3000 → 0.33 [10.8, **−3.1 to 15.9**].
  By N≳1000 the 16–84% band includes zero — the enhancement is consistent with zero within
  shot-noise scatter.

**F2 — the verdict computation.**
- **Scaling:** I_LL_aniso(N) ∝ N^{−α}, **α = 1.17 [1.09, 1.27]** (psr-median; best pulsar 1.23),
  as expected from shot-noise anisotropy (c_lm∝1/√N ⇒ variance∝1/N).
- **Resolved channel (Stage A, fixed-total):** knee N≈423, post-knee marginal info ∝ N^{−2.00}.
- **Enhancement / iso baseline:** 1.09 (N=100) → 0.40 → 0.06 → 0.024 (N=3000).
- **σ_L (WITHIN-FRINGE, best pulsar J0437-4715):** N=300 → full-covariance 77.5 pc (iso baseline
  98.7 pc, anisotropy-only 125 pc); N=1000 → 96.1 pc (iso 98.7, aniso-only 423 pc). Compare the
  resolved channel's ~2 mpc (Stage A `j_sigma_pc`): 4–5 orders of magnitude worse.

**>>> TRACK A VERDICT (2026-07-03): OUTCOME (ii) — channel real but formal; the null effectively
holds at realistic N.**
The anisotropic covariance channel is *real*: at low source count (N≲300) the enhancement rivals
or exceeds the iso baseline and follows a clean α≈1.17 power law. Its decay rate (α=1.17) is
*slower* than the resolved channel's N^{−2} past the knee — a formal rate-crossover. **But it is
not a useful crossover:** by N≈1000–3000 the enhancement is only 2–6% of the iso baseline and its
16–84% band already includes zero; the implied σ_L is ~80–100 pc and WITHIN-FRINGE (amendment 5),
versus the resolved channel's ~2 mpc. Per amendment 5 the covariance channel inherits the
fringe-identification problem, which the Anchor Census showed **fails** on this array (no
prior-alone anchor; certification is loud-source-driven). So the confusion-regime optimisation
route is **not** rescued by anisotropy: information survives with a measurable exponent but is
astronomically useless at realistic confusion N and sensitivity. Graduating it would need σ_L to
improve by ~4–5 orders — i.e. resolving sources, which is the resolved regime. (What would move
it: much longer T and lower white noise raise the per-bin SNR of the higher Whittle bins, but the
1/N anisotropy decay is set by shot statistics, not instrument.)

**Caveats.** (a) bins 2–4 ran CAPPED (5.8/3.9/2.9·ftau0 vs the strict 10×); they carry 52% of the
iso Fisher, so the absolute I_LL / σ_L values carry a resolution uncertainty bounded by G3 — the
α-scaling and the verdict are robust (driven by c_lm∝1/√N, resolution-independent). (b) The F2
figure plots the resolved channel (toy amplitude normalisation) and the covariance channel
(real-GWB normalisation A=2e-15) on shared axes; the large absolute vertical gap is partly a
normalisation artefact — the *meaningful* comparisons are the decay slopes and the physical σ_L,
not the absolute Fisher heights.

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

- 2026-06-25 (Stage A.3, cronus/4090) — Toy CAPSTONE: the N x t_c map. Combined confusion
  + coherence; dphi-block marginalised by per-pulsar N x N solves (block-diagonal over
  pulsars), validated vs dense single-Schur at N=8 (5.9e-11). Edge gates: y->0 == Stage-A
  confusion (4.2e-8), x=1 == Stage-A.2 coherence (4.7e-15). Headline: the transitions do
  NOT factorise — max|R - R(N,0)R(1,Y)| = 0.331 at (N/N*=12, Y=0.27), R<product, so
  decoherence COMPOUNDS confusion (curved 0.5 contour, not a rectangle). Hit+fixed a
  ratio-of-medians vs median-of-ratios bug (differ at N=1). Marked the TOY PHASE CLOSED;
  Stage C (discovery likelihood, real noise, absolute sigma_L) is next. New files
  prong2_Ntc_map.{npz,png} via `prong2_transition.py ntcmap`. (Claude + Matt)

- 2026-06-25 (Stage A.3 robustness, cronus/4090) — Re-checked the factorisation residual
  with COMMUTING reductions (no re-grid; used saved per-cell median info). Peak deviation
  0.331 -> 0.202 but stays < -0.1 in the confusion-dominated interior -> compounding is
  PHYSICAL, median-of-ratios inflated it ~50%. Commuting map shows a weak positive
  (de-confusion) lobe at N/N*<1 and the dominant negative (compounding) lobe at N/N*>3.
  Verdict kept (COMPOUND). prong2_Ntc_map_corrected.png + extra npz keys. (Claude + Matt)

- 2026-07-02 (Stage C, D0, cronus/4090, jug-gpu) — Derivative validation gate on
  build_fast_scan_likelihood under jax 0.10.1: grad 6.1e-8, Hessian dist-block 2.0e-8 /
  cross-block 6.4e-7, cho_solve-path-in-graph confirmed (grad moves 7.6e-2 rel under
  +0.1 dex frozen-GWB perturbation). ALL PASS. Key lesson: naive FD steps (1e-5·|x|)
  falsely fail by O(1) — oscillatory logL needs phase-aware steps (1e-2 of the 1-rad
  phase scale) + Richardson; h→h/2 error ratio 4.00 everywhere proves truncation-
  dominance. cho_solve FutureWarning persists, harmless. New files
  CW_transition/stagec_d0_derivative_gate.py, stagec_d0_results.npz. (Claude + Matt)

- 2026-07-02 (Stage C, D1, cronus/4090, jug-gpu) — Hessian infrastructure at scale.
  stagec_hessian.py with dense + HVP paths. GATE dense-vs-HVP 2.5e-14 (<1e-8), dense-vs-D0
  7.8e-14. 116-psr benches: HVP peak 0.44/0.62/1.12 GiB for N_CW=1/4/16 (production path,
  ~linear in N_CW); dense 2.35/3.42 GiB then OOM (12.32 GiB transient) at N_CW=16 — matches
  the predicted forward-over-reverse blow-up. Warm eval 1–2.7 s (HVP) / 0.3–0.4 s (dense).
  Shared-GPU hazard found + handled: co-tenant spike (Matt's smoke_ensemble.py + kyleg) OOMs
  XLA autotuning; disabled it (autotune_level=0), robust + peak unaffected. New file
  CW_transition/stagec_hessian.py. (Claude + Matt)

- 2026-07-02 (Stage C, D2, cronus/4090, jug-gpu) — Conditional vs marginal sigma_L,
  single realistic CW, 116 psr, 10 geometry draws. Built amortised zero-noise Fisher
  (stagec_fisher.py): residual DATA made a runtime arg via discovery's callable-residual
  vary path -> ONE compile (465 s) then 0.10 s warm/draw. Gate: amortised H == D1 dense H
  to 4.6e-16. Zero-noise => H=-Fisher, neg-curv 116/116 every draw (vs D1's noisy observed
  78/95/110). Per-pulsar marginal = 1/sqrt(Schur-complement diagonal) (eigh-pinv on 8x8 F_cc
  only; rcond spread 2.7e-2 stable) after an initial full-116x116-pinv version gave 254%
  spread (wrong metric, fixed). HEADLINE: median 116/116 pulsars have marginal sigma_L < EM
  prior; median sigma_L ~0.1-0.55 pc. marg/cond median 1.008 (16/84 1.001/1.033) — at N_CW=1
  the conditional metric does NOT overstate (array pins the 1 CW), matches toy N=1 ~1.03.
  CAVEAT flagged: within-fringe only, not a distance measurement until D3 fringe ID. New files
  CW_transition/stagec_fisher.py, stagec_d2_results.npz. (Claude + Matt)

- 2026-07-02 (Stage C, D3, cronus/4090, jug-gpu) — The fringe split. Per-pulsar fringe
  identifiability: dlnL(true vs best-wrong fringe in +/-3sig_EM) vs ln K_a (K_a=2*3sig_EM/dL).
  Caught+fixed a metric bug: find_best_wrong_mode_in_prior's prominence floor (0.5) exceeds
  the whole single-CW fringe modulation (<0.5) -> returns "no wrong peak" which naively reads
  as "identified" (backwards). Replaced with direct lnL eval at fringe centers (fixed 512-eval
  shape). Result: dlnL_a median 0.0000 (max 0.0028) << ln K_a median 7.94 -> class(i)=0/116,
  class(iii)=116/116 (median); 2/10 draws had 1 K-small edge case. SNR cross-check: sigma_L =
  dL/(2pi SNR_pterm) holds to ~10% with the timing-marginalised pulsar-term SNR (raw white
  overpredicts 1.6x -> GP/timing absorbs low-freq power). J0437: in array, marg 0.249<EM 1 pc
  but NOT fringe-identifiable -> class iii. HEADLINE: a single realistic CW gives ZERO genuine
  distance measurements; D2's 116/116 is within-fringe only. Enforced the language rule. New
  files CW_transition/stagec_d3.py, stagec_d3_results.npz. (Claude + Matt)

- 2026-07-02 (Stage C, D4, cronus/4090, jug-gpu) — N_CW sweep (1,2,4,8,16 equal strain,
  10 draws) on the real amortised Fisher + direct fringe scan; HVP-assembled Hessian for
  N_CW=16. Also triaged the D3 prominence-floor bug's blast radius across the 6 nbs (05
  realistic + 01 single suspect/optimistic; 03/04 loud safe; 00/02 mixed at N_CW=1). Results:
  class-i (genuine measurements) 0/0/0/0(max1)/0(max2)/5(N_CW=16, range2-10) -- tail-driven
  onset at N_CW=16, median pulsar never identifiable (dlnL 0.95<<lnK 8.8). Population draw (3
  loud+13 faint) at N_CW=16 -> class-i=0 (faint sources don't break fringes; behaves like
  ~3 loud). marg/cond stays ~1 (1.01->1.04): real likelihood is in the RESOLVABLE regime
  (far below toy confusion knee ~50 T Df); does NOT reproduce toy turnover -- instead fringe
  identification is the binding constraint the toy never saw. Within-fringe sigma_L improves
  0.27->0.067 pc. Fisher-valid reported ~116/116 (UNITS BUG, corrected in the Anchor Census
  2026-07-03: true valid 29/16/28/43 of 116 at N_CW=1/4/8/16; certifications unaffected).
  All NSD. New files CW_transition/stagec_d4.py,
  stagec_d4_results.npz. (Claude + Matt)

- 2026-07-02 (Prong-2 close-out P2-D + CLOSURE, cronus/4090, jug-gpu) — Loose ends.
  Item 1 strain reconciliation: per-source optimal SNR median h=-13.75 -> 10.7, h=-13 -> 59.9,
  h=-12 (CW_node_sampling optimizer) -> 599.5; optimizer/D4 = 56.2x (=10^1.75). Old optimizer
  20/20 recoveries ran ~56x louder (SNR~600), reconciled with realistic zero. Item 2 nb-01
  recheck: bug optimistic at all strains, SIGN-FLIPS at h=-13 (+296 original -> -108 corrected
  = separable -> not-separable). Wrote the PRONG-2 CLOSURE: info-theoretic content complete
  under both interpretations (distance-info transition A-C; detectability-vs-rangeability P2-A);
  one open handoff (environmental df/f > P2-B thresholds -> Taylor/Farr); pulsar-side coherence
  resolved in-house; all Sec-7 caveats resolved/bounded/assigned. CLOSED AS COMPUTATION, OPEN
  AS ONE HANDOFF. Files stagec_p2d.py, stagec_p2d_item{1,2}.npz. (Claude + Matt)

- 2026-07-02 (Prong-2 close-out P2-C, cronus/4090, jug-gpu) — Array-boost scaling. Reran
  Stage A knee at N_psr=15/30/60/116/200 -> knee/N* = 6.47/12.92/26.73/52.75/91.45. GATE:
  116 = 52.75 reproduces Stage A ~52 (PASS). Power law knee/N* = 0.40 N_psr^1.03, slopes
  1.02 low / 1.02 high -> LINEAR, no saturation to 200. Forecasting law: N_knee ~ 0.40 N_psr
  T Delta_f (confusion wall scales linearly with array size). N_psr=200 first OOM'd under a
  co-tenant (jit_median 10.5 GiB); rerun clean on the freed GPU. Files stagec_p2c.py,
  p2c_array_boost.png, stagec_p2c_results.npz. (Claude + Matt)

- 2026-07-02 (Prong-2 close-out P2-B v2, analytic) — Coherence-axis physical grounding;
  REDONE after a first-pass rejection. First-pass error: pulsar-side sigma_phi computed as
  2 pi f sigma_res (time-shift) = 6e-15 rad; correct is noise/signal sigma_res/A_CW =
  2 pi f sigma_res/h = 1/SNR_pterm (factor 1/h = 5.6e13, ~13.8 orders). Corrected: pulsar-side
  sigma_phi median 2.13 rad = 1/SNR_pterm (SNR_pterm median 0.47) = the SNR floor already in
  D2-D5, at the knee by construction, NOT independent coherence. GW chirp accumulated phase
  Delta_phi=pi fdot tau_p^2 = 0.02..2.9e4 rad (Mc,f) -- modelled, off-axis. Required stochastic
  df/f to cross knee: 3.4e-5/8.4e-6/1.7e-6 (SNR 5/20/100), magnitude OPEN. 2-panel clipped
  figure (df/f family + 116 per-pulsar sigma_phi). Verdict: pulsar-side = SNR floor (in-house);
  source-side environmental df/f the open question -> HANDOFF (Taylor/Farr). Files stagec_p2b.py,
  p2b_coherence_grounding.png, stagec_p2b_results.npz. (Claude + Matt)

- 2026-07-02 (Prong-2 close-out P2-A, cronus/4090, jug-gpu) — Detectability vs
  rangeability figure. DETECTABLE = Earth-term-only matched-filter SNR>5 (marginalised
  noise, via (s|s)=2[logL(off,0)-logL(off,r_earth)] on the amortised likelihood); RANGEABLE
  = D4 class-i. Result: detectable 1/1/3/7/13 (N_CW=1/2/4/8/16), rangeable 0/0/0/0/5;
  population (3 loud+13 faint) detectable 3, rangeable 0. Sky visible ~immediately, rangeable
  only at N_CW=16 (tail) and never for a realistic population. GATE: Earth-term SNR vs white
  matched filter ratio 0.638 = 1/1.57 = the D3 timing+GP absorption reproduced <1% (machinery
  validated; naive white is 1.6x optimistic). Files stagec_p2a.py, stagec_p2a_results.npz,
  p2a_detect_vs_range.png. (Claude + Matt)

- 2026-07-02 (Stage C, D6 close-out, cronus/4090, jug-gpu) — Consolidated the Stage C
  section: headline (within-fringe precision cheap but fringe ID is the binding
  constraint; single CW = 0 measurements, equal-strain onset N_CW=16 tail-driven,
  realistic population = 0), the "what changed vs the conditional joint-dlnL metric"
  paragraph (marginal not conditional; separate prominence-free fringe-ID test; old
  metric conflated curvature with measurement AND was optimistically biased), and the
  files list. Also ran the nb-05 add-on (D5): corrected joint dlnL 7-45% more negative
  than original across N_CW=1-8 (bug optimistic; sign unchanged). Stage C COMPLETE.
  (Claude + Matt)

- 2026-07-02 (Stage C, D5, cronus/4090, jug-gpu) — Frozen-GWB sensitivity at N_CW=8:
  frozen GWB log10_A at truth +/- 0.5 dex. Movement tiny: marg sigma_L x1.00-1.02, med
  dlnL_a x0.91-1.01, class-i unchanged (0-1). Much smaller than the D3 timing+GP absorption
  (1.6x SNR). Decision: GWB-amplitude freezing is NOT the dominant optimism; realism effort
  -> luminosity function (D4) + real noise, not un-freezing the GWB. Verified D4 marginal
  Schur-complements the FULL joint 8N_CW block (F_cc=F[nd:,nd:], not per-source). New files
  CW_transition/stagec_d5.py, stagec_d5_results.npz; nb-05 joint-dlnL recheck (add-on) in
  stagec_nb05_recheck.py. (Claude + Matt)

- 2026-07-03 (Anchor Census A0-A3, cronus/4090, jug-gpu) — Audited whether the real array
  has K~1-10 anchor pulsars with corrected EM priors (vs the wide feather priors). A0: ran
  canonical get_distance.py (live Cornell fetch, cached), matched to the 116 array pulsars
  by J-name (3 B->J aliases), froze best_distances.txt; coverage 69 parallax-class (30
  Catalog_PX / 1 VLBI / 38 Timing_PX) + 47 feather-kept; GATE J0437 0.1563 kpc / 1.32 pc
  (Deller08). A1 dual-prior K table: K<=3 = 0, K<=10 = 0, K<=30 = 2 (J0437, J2222) under
  BOTH canonical and min-sigma-optimal; smallest K_opt = 11.88 (J0437). GATE FAILED as
  stated -- get_distance's tightest J0437 is 0.97 pc (K~12), the Reardon+2016 0.25 pc
  (K=3.07) is a joint-timing composite absent from the catalogue; decision (Matt) to inject
  published composites for the best case. J1713 reversal confirmed (real 40-61 pc > feather
  19 pc, both directions). A2 reclassification (same seeds; feather/script/literature
  columns from one scan/draw): PRIOR-CERTIFIED (K<1) = 0 everywhere (no prior-alone anchor);
  flat ANCHOR-CERTIFIED lit 0/0/1/6 (N_CW 1/4/8/16) vs feather 0/0/0/5; CERTIFIED-BAYES
  (P_true>0.9, honest small-K) lit 1/0/3/18; population (3loud+13faint) N_CW=16 flat 0 (all
  priors) / bayes 3 (J0711/J1713/J1909, all DATA-driven & prior-independent -- certify under
  feather too). CAUSE-TAG AUDIT (per-draw dL vs sigma): J0437's N_CW=1-4 bayes certs are all
  GEOMETRY-WIDE (favorable large-dL draws, dL>=0.69pc incl. a 12.3pc antipodal; median dL
  0.489/sigma0.25 gives prior-only P_true 0.77<0.9); its first DATA-DRIVEN median-geometry
  cert is N_CW=8 (lit) / N_CW=16 (feather). VERDICT: no strict (K<=1) anchor exists even with
  composites (2603.28897's assumption fails); J0437 alone reaches K<=3 (only via Reardon
  0.25pc, and only as a geometry-dependent N_CW>=8 anchor); flat bar => realistic population
  certifies ZERO (conservative floor); Bayesian bar => population certifies 3 (loud-source
  data-driven, prior-independent) -- the real seed set for the sequential estimator, NOT
  prior-alone anchors. Also FIXED a units bug in D4/D2 Fisher-valid (marg[kpc] vs dL/4[pc]);
  true valid 29/16/28/43 of 116 (N_CW 1/4/8/16), certifications unaffected. Patched
  stagec_d4.py:115. OPTIMIZATIONS (approved): vectorized
  scan across all pulsars x draws (GATE vectorized==loop 0.0, ~35x faster than the per-pulsar
  loop), JAX persistent compilation cache, Hessian once/draw. Canonical priors file =
  best_distances.txt + LIT_INJECT (J0437 Reardon16, J1909 Reardon21). Also found+noted a
  units bug in D4's Fisher-valid flag (kpc vs pc). New files CW_transition/stagec_anchor_a0.py,
  stagec_anchor_a1.py, stagec_anchor_a2.py, best_distances.txt (frozen), anchor_a0_priors.npz,
  anchor_a1_Ktable.npz, anchor_a2_results.npz. (Claude + Matt)
- 2026-07-03 (Track A F0–F2, cronus/4090, jug-gpu) — Tested whether an anisotropic confusion
  background carries pulsar-distance info in the cross-pulsar covariance (Hogg's contention /
  the confusion-regime referendum). Built the pulsar-term cross-covariance C_ab(f) and per-pulsar
  distance Fisher under the Whittle likelihood on the real 116-psr array with real heterogeneous
  white noise and shot-noise anisotropic skies (real Ylm to l=10, c_lm∝1/√N, unit monopole).
  Amendments 1–5 all implemented (baseline subtraction, Whittle bins, aliasing gate, auto-term
  sanity, within-fringe caveat). Result: the anisotropic ENHANCEMENT I_LL_aniso(N) ∝ N^−1.17
  [1.09,1.27] — real and measurable at low N (N≲300 it rivals the iso baseline) but decays to
  2–6% of baseline by N=1000–3000 (16–84% band includes zero), with σ_L ~ 80–100 pc WITHIN-FRINGE
  for the best pulsar vs the resolved channel's ~2 mpc. Decay is slower than the resolved N^−2
  (formal rate-crossover) but astronomically useless and fringe-limited → VERDICT (ii) channel
  real but formal; the null effectively holds, the confusion-regime optimisation route is not
  rescued by anisotropy. Gates: G1 monopole=iso enhancement 0.00 exactly; C00 matches F0; Whittle
  auto ∝ (π/3)ftau0 (ratio 0.787), cross/auto ~1/ftau0; iso info splits 47.6%/52.4% bin1/bins2–4
  (auto ∝ f offsets the red spectrum). Perf note: RTX 4090 fp64 is 1/64 fp32 → used complex64 GEMM
  with fp64 outer accumulation; a naive einsum path allocated a 30 GiB temp that OOM-retry-looped
  (fixed with a hard cp cap + explicit contraction). New files CW_transition/trackA_fisher.py,
  trackA_realarray.npz, trackA_f1_results.npz, trackA_f2_results.npz, trackA_f2_scaling.png,
  trackA_gates.npz. (Claude + Matt)
- 2026-07-06 (Track B B0 — cold-start estimator + the precision-cliff finding, cronus/4090,
  jug-gpu) — Built the DAEM fringe-identification estimator (CW_transition/trackB_estimator.py)
  to reach the Anchor-Census P(true) ceilings from a COLD start. **B0.0 PASS**: jax persistent
  compilation cache verified cross-process (fisher_logL 2.6 s vs 465 s cold; batched/HVP 38/52 s);
  A2 vectorized scan gate max|loop-vec|=0.0. **B0.1 built + validated**: E-step = A2 Bayesian
  softmax, reproduces the pop ceilings EXACTLY at truth (J0711 0.953 / J1713 0.989 / J1909 0.998,
  diff 0.000). M-step: the coherence-BLEND annealing ((1-b)lnL_earth+b lnL_full) was found to
  BIAS truth (earth-only optimum absorbs pterm power -> drifts away, f decreases, ceilings->0);
  CORRECTED to full-lnL ascent at confidence-gated fringes, which holds truth as an exact fixed
  point (from-truth dphi=0.00, ceilings stable across the whole beta ladder). Adam loop JITted
  (lax.fori_loop) + bounds/nan guards. F-stat SEED scan (truth-blind matched filter, 47 freq x
  192 healpix sky, amp/phase profiled) recovers 3/3 loud sources (6-12 deg, <1 freq-bin, 2F
  184-595 vs floor 15) -> seeding gate 3a PASS. **B0.2 Asimov gate FAIL, and this is the
  make-or-break finding**: seed+EM converges but lands median source err ~0.5 (scaled) and the
  ceilings collapse to 0. The `perturb` scan measures WHY: the certification tolerance on the
  SOURCE parameters is ~1e-4 scaled (~0.006 deg sky for J1713/J1909; sky/geometry-dominated,
  freq tolerates ~3e-4), set by the pulsar-term phase LEVER ARM 2*pi*L/dL ~ 1.6e4 rad -- any
  source error shifting that phase by ~1 rad re-registers the fringe pattern and the true fringe
  stops being the posterior mode. Achievable source precision ~0.5 scaled is 3-4 ORDERS coarser
  -> **the A2 census "3 certified pulsars" is a CONDITIONAL-on-truth ceiling, NOT achievable by
  any cold-start estimator, even at zero noise**. P0 (miscalibration anatomy) confirms the
  posterior is CONFIDENT-WRONG at imperfect sources (q_max 0.5-0.99 on a SHIFTED fringe, k!=0),
  not diffuse -- a wrong-certification machine; the source error is absorbed by per-pulsar integer
  shifts (GPS-RTK structure). **PIVOTED (with Matt) to the joint-registration reformulation** (the
  honest object = joint posterior over source params + fringe integers, LAMBDA/RTK; certification
  <-> localization exchange rate via the same lever arm): P0 done; P1 (registration-needle map)
  attempt 1 abandoned at 64/169 coarse trials (~33 s/trial, shared-GPU; per-trial Python loop
  untenable). Partial baseline CW_transition/trackB_p1_partial.npz (grid positions + timing;
  per-trial lnL/reg were in-RAM, NOT persisted -> NaN, unrecoverable); log
  CW_transition/logs/p1_coarse_attempt.log. Next session: batched trial evaluation (trials as
  runtime args to one compiled graph, vectorized fringe_posterior, staged zoom w/ per-stage npz
  checkpoints) to test whether a unique needle exists at truth (width ~0.006 deg), then P2 (joint
  solve + B0.2 re-gate) / P3 doc. If no needle -> that is the Track B verdict (ceilings are
  conditional artifacts). Edited stagec_fisher.build_fisher_amortised (+optional include_pterm /
  msd_factory, non-breaking) for the Earth-only decohered likelihood. Full estimator spec +
  findings in CW_transition/trackB_estimator_spec.md. NOT committed. (Claude + Matt)

- **2026-07-06 (Track B P1 — the registration NEEDLE, cronus/4090, jug-gpu) — THE NEEDLE
  EXISTS.** Rebuilt the P1 registration map batched (`CW_transition/trackB_p1_map.py`).
  **Blocker first:** the "salvaged 44 trials" gate baseline did not exist --
  `trackB_p1_partial.npz` had profiled_lnl all-NaN / registration_count all -1 (in-RAM
  only, never persisted; killed at 64 not 44). Redefined the gate (with Matt) vs a live
  D=1 `scan_current` reference on the deterministic grid. **Batched rebuild = single vmap
  of fisher_logL over a flat (trial×pulsar×fringe) stack vs the ONE fixed data_obs**
  (all trials share it), avoiding attempt-1's pathological (D,4096,244) double-vmap.
  GATE: batched == D=1 reference **0.00e+00 (bit-exact)**. Measured **0.544 ms/eval →
  91 min per 13×13 stage** (no occupancy win; the full-array eval IS the wall). DECISION
  (Matt): keep batch-only, size grids < 45 min, DEFER the low-rank QuickCW split to B1
  (algebra traced+verified vs discovery/matrix.py, recorded spec §6.1: distances enter
  only per-pulsar-local a_p/b_p, GWB coupling confined to the frozen outer cf_cached ->
  rank-block update, ~1000×/fringe, pays across B1's ≥20 realisations).
  **RESULT (pop draw = 3 loud + 13 faint, N_CW=16, lit priors; loud0 sky scanned):**
  a UNIQUE joint-lnL registration needle sits EXACTLY at truth at every scale --
  coarse ±20° (9×9, data-driven): argmax at truth, next-best cell −102 nat, only the
  centre within 50 nat; zoom1 ±2°: argmax 0.00° from truth; uniqueness patch ±0.05°
  (centre = zoom1 argmax, DATA-DRIVEN): argmax at truth, **interior 2-D maximum** (not a
  ridge/saddle), floor −60…−79 nat. **Width:** the profiled joint-lnL needle is an
  unresolved CUSP -- −26 nat at one 0.0019° sky step, −126 nat at one 0.25 freq-bin step
  -- FWHM below the finest cut in BOTH sky axes and frequency, consistent with the
  ARRAY lever-arm prediction max_p 2π L/dL → needle width **~1e-5°** (the brief's 0.006°
  is the J1713-class value; the array max is far larger). **Registration needle (the
  physical certification tolerance):** true-reg (q>0.5 & fringe k=0) = **18/116 at truth**
  vs floor 2–3 (max-off 6–7); the census-3 (J0711/J1713/J1909) register **3/3 at truth,
  0–1 off** -- and the registration HALF-width ≈ **0.006° sky (18→9)**, reproducing B0.2's
  ~1e-4-scaled certification tolerance as a needle. **4c (needs-all-116?):** the needle is
  carried by ~18 loud-source-broken pulsars, census-3 cleanly among them; not just 3, not
  all 116. All stages 100% finite (per-stage npz + finiteness-abort checkpoints). Labels:
  the 1-D width cuts are TRUTH-PLACED; coarse/zoom1/patch-centre are DATA-DRIVEN (P1 is
  the physics referendum; the truth-blind end-to-end zoom is P2's gate, not yet run).
  **VERDICT: the census ceilings are ACHIEVABLE IN PRINCIPLE -- the combs co-register
  ONLY at truth, so a joint source+fringe solve HAS a unique global optimum at truth
  (the GPS-RTK / LAMBDA structure is real).** The open question is now purely
  SEARCH/REACHABILITY. **TERRAIN (corrected from the 2-D radial means, not the central
  axes):** the profiled(HARD-snap) joint-lnL has THREE regimes -- (1) smooth radial funnel
  > ~2 deg (ring means -186/-301/-418/-549 nat at 5/10/15/20 deg; the F-stat Earth-term
  envelope), (2) a near-FLAT noisy plateau ~2 deg -> ~0.05 deg (ring means -72/-76/-78/-87
  at 0.25/0.5..2 deg, noise +/-11, weak inward tilt ~8 nat/deg), (3) DEAD FLAT floor +
  a 0.006 deg delta cusp < ~0.05 deg (patch ring means -68/-69/-70/-69, noise +/-5). So
  there is a FLAT GAP (~2 deg -> 0.006 deg) between the F-stat basin and the needle where
  the hard-snap has no gradient -> a local funnel-descent search STALLS at ~2 deg and never
  reaches the cusp. This is the geometric form of B0.2 (achievable ~degrees, cert tolerance
  ~0.006 deg, flat in between). **COROLLARY (lever arm 2*pi*L/dL as the governing quantity):**
  the lnL cusp width ~1e-5 deg ~ few x 10 mas -- a LOCKED registration localises the source
  at mas class, so certification and localisation are the SAME measurement, exchange-rate =
  the lever arm; but the same lever arm makes the target needle mas-narrow, which is why the
  gap is the problem.
  **PROBE 1 RESULT (2026-07-06, corrected interpretation) -- the plateau is a GAUGE
  CONSPIRACY, not absent physics.** The flat plateau is NOT a per-source flatness: it is
  SIMULTANEOUS RE-REGISTRATION FREEDOM. In P1's all-snapped HARD surface every one of the
  116 combs is free to shift together, so for any displaced source there EXISTS a globally
  consistent WRONG-fringe distance set that refits nearly as well (joint profiled lnL back
  to the ~-70 nat floor) -- a gauge conspiracy of the integer fringe latents. The moment
  the array is ANCHORED (other distances held near truth) the conspiracy is broken and the
  source-sky gradient is STEEP and smooth across the WHOLE gap: with all other distances at
  truth, moving loud0's sky costs **-82750 nat at 2 deg** (monotone to the cusp), and the
  true-fringe posterior weight QTRUE = sum_p q_p(k=0) rises smoothly **4.9 (2 deg) -> 9.6
  (0.5) -> 20 (0.1) -> 23.3 (0.05) -> 24.45 (truth)** -- a coherent gradient spanning the
  entire 0.05-2 deg gap that the hard-snap joint hides. **So the funnel is CONDITIONAL on
  anchoring; the information is present (soft/posterior) but gauge-hidden in the
  simultaneous fringe freedom.** (Caveat: PROBE 1 held others at truth = full anchoring;
  the loop-gain question is how MANY anchors restore a followable gradient -> PROBE 2.)
  **PROBE 2 RESULT (2026-07-07) -- LOOP CLOSES; the gap is SOFT-NAVIGABLE.**
  Anchor-fraction sweep (hold k pulsars at truth fringes, self-consistent 3-pass E-step,
  JOINT lnL objective = the real algorithm's surface; anchors = strongest true-fringe
  registrars, which turn out to BE the census set: idx 88/62/19 = J1909/J1713/J0711, q0
  1.00/0.99/0.95). HARD gap-drop(0.05->2 deg) grows monotonically with anchors:
  **k=0->50, k=1->50, k=3->132, k=6->172, k=12->244, k=24->324 nat**; STRICT monotone-across-
  the-full-gap closes at **k=6**. Two decisive refinements of the P1 picture:
  (a) **the SELF-CONSISTENT (3-pass) E-step already recovers gradient inside ~1 deg even at
  k=0** -- HARD is monotone cusp->1 deg (0->-66 nat); P1's single-pass snap had flattened this
  to ~-70. Only the 1-2 deg far edge stays a noisy plateau (+/-7 nat), and that zone is
  covered by the F-stat + coarse funnel (>2 deg, strong gradient).
  (b) **QTRUE (sum_free q_p(k=0), the soft-EM Q-function's driver) is MONOTONE across the
  WHOLE gap from k=0** (6.25 at 2 deg -> 24.45 at truth), i.e. the SOFT/posterior objective is
  followable end-to-end with NO anchors; anchoring only sharpens the HARD-MAP surface.
  census-3 anchors (the data-driven set the pipeline actually bootstraps, prior-INDEPENDENT)
  give a 132-nat gradient monotone to 1 deg. **VERDICT = pre-registered FAVORABLE branch:
  loop closes at census-class anchoring (k=3-6, NOT k>>3; anchors are data-driven, not
  prior-limited) -> P2 = SOFT-EM ALTERNATION is viable; the census ceilings are ACHIEVABLE
  by a cold-start soft-EM pipeline (modulo the real B0.2 truth-blind gate).** The flat P1
  plateau was a HARD-snap + single-pass + no-anchor artifact on three counts; none survive
  the soft, self-consistent, lightly-anchored treatment.
  Files: trackB_p1_map.py, trackB_p2_probe.py (probe1: line/freqpolished),
  trackB_p2_probe2.py (probe2: anchor sweep), trackB_p2probe2_anchorsweep.npz,
  trackB_p2_anchorsweep.png,
  trackB_p1_{coarse,zoom1,patch,cut_cos_{wide,fine},cut_gwphi_{wide,fine},freqzoom}.npz,
  trackB_p1_registration_needle.png, logs/p1_{gate,production}.log,
  logs/p2_probe{,2}.log. NOT committed. (Claude + Matt)

- **2026-07-07 (Track B P2 pipeline -- soft-EM M-step de-risk FAILS; needle needs the
  INTEGER solve) -- the LOCK leg.** Built the truth-blind soft-EM cold-start pipeline
  (`trackB_p2_pipeline.py`: F-stat seed -> soft E-step w/ emergent anchors -> soft M-step
  = Adam ascent of the EM Q-function Q(src)=sum_p sum_n q_p(n) lnL(src,L_p(n)) over the
  8*N_CW source params -> certify). Memory: the vmap-all Q-gradient OOMs at 32 GiB; fixed
  with a scan over mixture configs (grad of a sum = sum of grads, one array-lnL backward
  at a time). **DE-RISK GATE (test-mstep, from-truth fixed-point + certification):
  DECISIVE NEGATIVE.** From TRUTH the soft M-step drifts the REGISTRATION params
  (sky ~0.05-0.10, freq ~0.04-0.05 scaled) and census-3 P_true collapses 0.953/0.989/0.998
  -> ~0 (ncert 3->0), across ALL THREE M-step weightings tried: equal, q_max-confidence,
  and hard-gate-to-emergent-anchors (q_max>0.9). Mechanism, now concrete: with all pulsars
  the ~98 fringe-degenerate ones outvote the ~18 confident anchors and pull sky/freq off
  truth; hard-gated to the 3 census anchors the 24 loud-source params are UNDER-determined
  and wander. **ROOT CAUSE (fundamental, not tuning): the fringe-marginalised (soft)
  objective's optimum is BROAD and DISPLACED from truth; a continuous gradient-ascent
  M-step converges only to ~0.01-0.05 scaled -- a precision floor 100-500x COARSER than
  the ~1e-4-scaled (mas) certification needle.** Continuous optimisation of a smoothed
  objective cannot reach a cusp. This is the LOCK-leg failure and it CONFIRMS B0.2
  geometrically: soft-EM (probe 2) NARROWS the reachable precision (~0.5 -> ~0.05 scaled)
  but does NOT close the gap to the needle. **P2 VERDICT: the census ceilings are NOT
  achievable by continuous soft-EM; the needle is reachable ONLY by FIXING the fringe
  integers first (RTK/LAMBDA) -- then the source localises to mas via the sharp registration
  surface (probe 1: -82000 nat / 2 deg anchored).** So the honest cold-start bridge is the
  discrete integer solve, NOT a smoothed continuous optimiser. Reconciles the whole arc:
  P1 needle exists (unique optimum) + probe1 gauge-conspiracy + probe2 loop-gain
  (anchored gradient) + P2 soft-EM precision floor => certification == localisation ==
  fixing the integers; that is an INTEGER least-squares (LAMBDA) problem, not gradient
  descent. Files: trackB_p2_pipeline.py, logs/p2_testmstep.log. NOT committed. (Claude + Matt)

- **2026-07-07 (P2 STEP 1a line scan -- MECHANISM CORRECTED; not a displaced optimum).**
  `trackB_p2_linescan.py`, `trackB_p2_linescan.{npz,png}`. Static scan of the SELF-CONSISTENT
  soft marginal (others->MAP per point) along loud0 gwphi & cos_gwtheta through truth. RESULT:
  HARD (max_n) AND SOFT (logsumexp_n) marginals BOTH peak SHARPLY and EXACTLY AT TRUTH
  (argmax 0.0000 scaled both axes; -8000 nat already at +/-0.003 scaled), and HARD~=SOFT to
  ~45 nat. **The prior entry's "soft optimum BROAD and DISPLACED ~0.05, precision floor
  100-500x" was WRONG mechanistically.** The soft marginal is NOT broad and NOT displaced --
  it is a sharp CUSP at truth on a noisy plateau (the same integer-ambiguity structure as P1;
  marginalisation does NOT smooth the peak because the true-fringe term dominates). The
  test-mstep M-step drift (0.05, census->0) is an OPTIMISER artifact, two causes: (1) the
  M-step holds q FIXED through 40 Adam steps and the FIXED-q Q (unlike the self-consistent
  marginal) IS displaced by the ~98 degenerate pulsars' wrong-fringe terms; (2) Adam falls
  off the sharp cusp TIP onto the plateau and wanders. A fully self-consistent soft-EM
  (E-step per micro-step) would track the sharp marginal -- expensive, not impossible.
  **CORRECTED VERDICT (direction UNCHANGED, mechanism fixed): continuous optimisation is
  defeated by the CUSP-ON-PLATEAU (integer ambiguity), not by a smooth precision floor;
  fixing the fringe integers (LAMBDA/RTK) turns the cusp into a smooth quadratic -> Newton
  converges. The census ceilings remain reachable only via the integer solve.** STEP 1b (the
  end-to-end soft-EM run) will document the plateau STALL + the float landing zone (integer
  support per pulsar = LAMBDA's search space). Next decision: LAMBDA scoping + feasibility
  probe (STEP 2). NOT committed. (Claude + Matt)

- **2026-07-07 (low-rank QuickCW split BUILT + all gates PASS -- the ~47x/ E-step + the
  microsecond candidate scorer).** `CW_transition/trackB_split.py` (`class Split`). Motivated
  by the honest cost review: every Track B test is bound by the full-array fisher_logL
  (0.544 ms/eval) x 116x512 = ~32 s per source position; deferring the split (spec 6.1) was
  the false economy. Built per the SPEEDUP_BRIEF + 7 amendments. Reconstructs the frozen
  per-pulsar rednoise Woodbury by mirroring discovery.matrix.make_kernelterms_vary (reuses
  dsm.* -> bit-identical by construction; a_const via difference vs fisher_logL at truth);
  distances enter only per-pulsar a_p/b_p, GWB coupling in the frozen outer cf_cached (verified
  6.1 algebra). Enabling edits (additive/safe): stagec_fisher.build_fisher_amortised returns an
  `internals` dict; trackB_estimator.build_problem gains build_earth flag (default True; split
  uses False). GATES (pop, N_CW=16, 116 psr): (0) scalar recon split.lnL==fisher_logL 0.00e+00;
  (1) E-step per-pulsar max|dlnL| 1.16e-10 (0/116>1e-8), soft-posterior max|dq| 9e-12; (2)
  source-gradient 0.00e+00/4e-13 (frozen M not differentiated); (3) small-signal adjacent-fringe
  diff 5.8e-11 at a 1.17e-3-nat step (chirp retained); (4) scorer 16 assignments 1.16e-10. All
  PASS. PROFILING (honest): E-step 0.68 s = **47x** vs 32 s (bounded by single-pulsar yprod
  recompute, NOT 1000x); candidate scorer **37 us/candidate at K=1e5** (~1e6 in 37 s -- LAMBDA's
  budget number). B_CERT=512 enforced. WIRED into trackB_p2_pipeline.py behind --use-split
  (default OFF). Batch-only batched_scan stays the correctness ORACLE (no shared code). Customers
  = LAMBDA feasibility probe + B1 (NOT started -- separate deliverables). Files: trackB_split.py,
  SPEEDUP_BRIEF.md, logs/split_{gate0,gate1,post}.log; edited stagec_fisher.py, trackB_estimator.py,
  trackB_p2_pipeline.py, trackB_estimator_spec.md sec 6.2. NOT committed. (Claude + Matt)

- **2026-07-07 (Track B LAMBDA feasibility probe -- GO/NO-GO -> NO-GO; the Track B verdict,
  cronus/4090, jug-gpu).** Seeded by the truth-blind soft-EM float, does a bounded integer
  (fringe) search + conditional source re-solve recover truth's carrier integers on zero-noise
  Asimov? **VERDICT: FAIL (iii).** L0 (`trackB_L0_float.py`): the truth-blind soft-EM (F-stat
  seed -> soft M-step, `--use-split`, per-iter checkpoints) DIVERGES -- loud sources drift to
  4.31 scaled off truth by it19, ZERO carriers emerge (q_max spectrum tops 0.52), census MAP
  fringes land at n=163/222/-2 (truth n=0 for all 116 psr, dist==L0) with ZERO posterior mass on
  the true fringe. (Soft M-step sped 7.7x -- 652->85 s/iter -- by compacting the config scan to
  its nonzero-weight rows, exact.) L1: no carriers; census candidate sets (>1% mass) span
  n in [-182,226], none contain 0. L2 (`trackB_L2_probe.py`): staged `score_candidates` at the
  float fixes WRONG integers (census -113/-143/-19, margin 0.198 nat), P_true=0/0/0 (ceilings
  .953/.989/.998), certified={}. L2b/L2c -- the decisive mechanism: at EXACT truth (integers
  fixed) |grad|=0 and the 24-loud-param Hessian is neg-def with eig [-5.96e11,-14.4] (cond ~4e10)
  -> truth IS a sharp local max (P1 needle real; L2's "needle-absent" auto-label was an LBFGS
  artifact, withdrawn), but the registration-cusp curvature ~6e11 (lever arm 2*pi*L/dL~1.6e4)
  makes the quadratic bowl ~1e-6 scaled: **conditional re-solve pull-in radius < 1e-4** -- neither
  LBFGS nor damped-Newton locks (sky<0.006deg) from ANY offset >=1e-4 even WITH truth's integers
  (from 1e-4: lands 0.046deg / -622 nat). FAIL is structural at BOTH gates: (1) the O(1)-scaled
  cold float puts truth's n=0 outside the bounded integer neighbourhood; (2) even given truth's
  integers, the source basin (<1e-4 scaled) is 4+ orders inside any achievable float. Same
  precision cliff as B0.2, now shown to defeat the DISCRETE solve too. GPS-RTK works because its
  float lands within one carrier cycle; the PTA F-stat float lands ~1e4 cycles away and no fixing
  step spans that. NO-GO for the full LAMBDA build; B1 would only reconfirm at higher cost. Spec
  sec "LAMBDA feasibility probe" has the correspondence table + numbers. Files: trackB_L0_float,
  trackB_L1_spec, trackB_L2_probe, trackB_L2b_basin, trackB_L2c_newton (+ .npz + logs). NOT
  committed. B1 NOT started. (Claude + Matt)

- **2026-07-07 (LAMBDA probe -- verdict CONTESTED then RE-EARNED via the lever-arm ladder;
  F1/F2/F3, cronus/4090).** The prior NO-GO was seeded by a soft-EM float that DIVERGED, so
  three bounded follow-ups. **F1 (float contaminated 3 ways):** (a) the 7.7x M-step compaction is
  NOT exact -- one M-step full-348 vs compact-72 differs |dQ|=4.2e4, max|dsrc|=0.20 scaled;
  root cause = `make_soft_mstep` does nan_to_num AFTER the config sum, so ONE NaN per-config grad
  (even a w=0 config, 0*NaN=NaN) zeros the WHOLE M-step gradient -> the full path is often
  NaN-zeroed, compaction drops the poisoning rows -> a different (real) gradient. My L0 float used
  the compact path, not the pipeline's; "exact" claim WITHDRAWN. Correct fix = per-config
  nan_to_num before summing (not applied; F1 report-only). (b) L0 endpoint is NOT a label swap
  (best perm 130.5 vs diag 138.0 deg); the "4.3 scaled" is log10_h (amplitude, extrinsic) --
  registration-critical offsets are src0 0.18 / src1 0.10 / src2 1.2 scaled. (c) descent trace:
  F-stat seeds only 2/3 loud (6.4 deg); the 3rd loud slot is a spurious 2F=350 noise seed at 109
  deg (loud2's real seed buried at rank-13 > top-16), and the soft-EM walks the two good seeds off
  IMMEDIATELY (it0 6.4->18 deg), labels permute by it4 (co-source confusion). **F2 (the ladder --
  EARNS the NO-GO):** per-(pulsar,loud) sky registration tolerance (phase 2*pi f L (1-cos mu)/c
  shifts 1 rad): loosest 1.85e-3, median 3.8e-5, tightest 1.3e-6 scaled; ZERO pairs tolerate >=
  0.05 (or >= 1e-2). Freq looser (2.5e-2) but sky binds (0.05 sky err -> 27 rad phase at the
  loosest pulsar). Ladder does NOT span: 27x (1.4 dex) wall between the blind float floor (~0.05)
  and the loosest rung (1.85e-3), zero rungs between -> the wide-lane cascade cannot fix its first
  rung. **F3 NOT run** (pre-registered: only if the ladder spans). **EARNED VERDICT: NO-GO**, now
  independent of the contaminated float: TOP (F2) loosest registration baseline 1.85e-3 scaled vs
  best blind sky float 0.05 scaled; BOTTOM (L2b/L2c BANKED) fixed-integer pull-in < 1e-4. The whole
  ladder < 1.85e-3, any float >= ~0.03-0.05 -- >=1.4-dex separation with the razor cusp beneath.
  F1's M-step NaN-poisoning bug + loud2 seeding gap are real but SEPARATE repairs that don't move
  the wall. Files: trackB_F1a_gate, trackB_F1c_descent, trackB_F2_ladder (+ .npz + logs). Spec
  "LAMBDA probe -- F1/F2 AMENDMENT". NOT committed. B1 NOT started. (Claude + Matt)

- **2026-07-07 (TRACK B CLOSE-OUT -- repairs, wall height, TERMINAL VERDICT; cronus/4090).**
  Repairs + the quoted number. **Repair 1 (M-step gradient):** F1a's NaN-poisoning diagnosis was
  TESTED and REFUTED -- the per-config nan_to_num fix (trackB_p2_pipeline `_cfg_grad`) is a NO-OP
  (byte-identical gate), and trackB_step1_diag shows full-vs-compact gradients agree to 2.2e-16 at
  ONE step [1] with identical config sets [0]; the F1a 0.2-scaled divergence is multi-step FP-order
  AMPLIFICATION through the ~4e10-conditioned Adam (the M-step is numerically CHAOTIC, matching
  F1c). nan_to_num relocation KEPT as hygiene, NOT credited. **Repair 2 (seeder):** sky-only NMS
  exclusion 25 deg (trackB_estimator `pick_seeds`) kills loud0/1 cross-frequency sidelobes; GATE
  3/3 loud, loud2 promoted rank-13/17.70deg -> rank-11/11.88deg, 40 cand seeds vs 501. **Step 3
  (wall height):** one repaired truth-blind float -- soft-EM does NOT converge, chaotically walks
  the seed off (co-source confusion + M-step chaos); registration-critical ceiling (all-16 matched)
  = 0.05 (F-stat seed floor) -> 0.21 (endpoint) scaled sky, transiently touching J1713 P_true=0.977
  at it1. WALL = ceiling/loosest-rung(1.85e-3) = 27x..112x = **1.4..2.05 dex, ZERO rungs**
  (trajectory-sensitive +/-0.2 per Repair 1). **Lever (i) tol~1/L CONFIRMED:** the 3 loosest F2
  rungs are the nearest pulsars -- J0711 (0.106kpc,1.85e-3), J1630 (0.089kpc,1.39e-3), J0437
  (0.155kpc,1.02e-3) vs median L=1.38kpc. **TERMINAL VERDICT: NO-GO for COLD-START POINT
  ESTIMATION** (earned, float-independent: TOP F2 loosest 1.85e-3 vs float ceiling 0.05-0.21;
  BOTTOM L2b/c pull-in <1e-4, banked). Census ceilings REAL but CONDITIONAL (need ~1e-3-scaled
  source knowledge the array can't self-generate). TERMINAL CHARACTERISATION RESOLVED by deliverable
  R (2026-07-08, entry below): the fringe-marginalised posterior SMEARS -> f=6.9e-7 -> NO-GO deepens
  to INFORMATION-THEORETIC (prong-3 does not win on the data alone; the J1713=0.977 hint is
  overwhelmed in the aggregate evidence). DESIGN THEOREM recorded (levers i nearby-pulsar
  lanes / ii eccentric harmonics=E-track's sharpened job / iii lower the float via EM+anchor+Q1).
  Track B CLOSED; machinery banked (split, scorer, E-step oracle, repaired seeder, float); queue
  head = E-TRACK; B1 reframed to the CONDITIONAL (source-localised) pipeline. Contamination
  asterisk on the pre-repair L0/L1/L2 float chain; L2b/c + F2 unaffected. Files: trackB_F1a_gate,
  trackB_step1_diag, trackB_step2_seeder, trackB_step3_ceiling, trackB_L0_float(repaired) (+ .npz +
  logs). Spec REPAIRS LOG + TERMINAL VERDICT + DESIGN THEOREM. NOT committed. (Claude + Matt)

- **2026-07-08 (Track B deliverable R -- THE POSTERIOR REFERENDUM; VERDICT (ii)
  INFORMATION-THEORETIC NO-GO, cronus/4090, jug-gpu).** Point estimation was settled NO-GO; R
  answers the terminal question -- does the fringe-marginalised POSTERIOR concentrate at the needle
  (truth) or smear across the plateau? Built `trackB_R_referendum.py` on the low-rank split. Object =
  honest COUNT-ONCE star-topology fringe-marginal `lnL_marg(theta) = lnL_ref + sum_p LSE_n[(LNL_p(n)
  -lnL_ref)+lnprior_p(n)]` (loud source params vary; faint 13 FIXED at truth; base dist L0). GATES
  all PASS: split==fisher 0.00; batched E-step==oracle 1e-10; exact lnL_ref==fisher_logL(truth) 0.00;
  exp(-m_p)==A2 P_true census (0.953/0.989/0.998) 3.7e-8; additive/star approx error 0.14 nat at the
  116-pulsar MAP config. `lnL_marg(truth)=405686.34 = lnL(truth) 405413.51 + fringe-entropy 272.83`.
  - **Measure-consistent reduction:** f computed over the SKY registration dims (6; freq+extrinsics
    FIXED at truth). Their Laplace factors CANCEL in f (source detected across the basin); F2 shows
    SKY binds. The 24-D broad-extrinsic box is intractable (narrow good-fit ridge -> SMC acc->0.04 +
    dilutes the plateau). Favorable conditioning -> verdict conservative-safe.
  - **R1 Z_needle (truth-placed quadrature, labelled):** ln Z_needle = **405617.64** (quad) /
    405617.84 (Laplace), agree 0.2 nat, doubling STABLE, 6/6 pos-curv, needle sigma_sky 2.4-9.7e-6
    scaled (razor-sharp; marginalisation does NOT broaden the L2c cusp). = g0 - 68.7.
  - **R2 Z_plateau (tempered SMC, truth-centred +-2deg sky, cov-preconditioned, needle excised):**
    ln Z_plateau = **405631.83 +- 0.053** (2 seeds agree 0.05 nat -> robust; high-beta acc-dip biases
    it DOWN -> f firmer). Plateau ceiling ~405659 (28 below the needle peak, 40 ABOVE Z_needle);
    posterior collapses to the high-L shell just outside the needle. (6-deg F-stat seed offset =
    nside=4 healpix GRID artifact, not a bias; truth-centred box = the EM-targeted scenario.)
  - **R3 the REFERENDUM NUMBER:** ln Z_plateau - ln Z_needle = +14.19 -> **f = 6.9e-7**. The needle is
    28 nat higher at its PEAK but its 6-D sky volume (~e^-73) loses to the +-2deg plateau by 14 nat.
    **Break-even sky prior theta* = 0.188 deg/loud source** (V*/V=6.9e-7): needle wins only if the sky
    is pre-localised < 0.19 deg -- 32x tighter than the F-stat floor (~6 deg); an EM host supplies it
    (lever iii). Census-3 P(true) at the posterior peak (plateau MAP) = 0/0/0.
  - **R4 VERDICT (ii):** the honestly marginalised data do NOT prefer truth; NO-GO is
    INFORMATION-THEORETIC, not search-only -- prong-3 does not unpark on the data alone. The transient
    J1713=0.977 is real per-pulsar but overwhelmed in the aggregate evidence. Only more data / external
    sky helps. Design theorem (lever ii E-track eccentric harmonics; lever iii EM sky) is the forward
    story.
  - **R5 T-crossing pencil (banked, PENCIL_t_crossing.md), quoted:** Earth-term float ceiling reaches
    the F2 loosest rung (1.85e-3) at T~11 kyr (sigma_sky~T^-1/2) or SNR~289 (h~10^-12.3, 27x current)
    at fixed T; L2c pull-in at T~3.75 Myr. The break-even 0.188 deg is the INFERENCE-side counterpart
    of the pencil's point-estimation wall -- both say the 15-yr array cannot self-localise the sky
    finely enough, by optimisation OR marginalisation.
  - **POSTMORTEM (bounded, banked plateau samples; `trackB_R_postmortem.npz`):** does any single
    pulsar's EM prior+lever arm near-decide its OWN fringe while the joint sky is not? YES, exactly
    ONE -- **J0437-4715** (tight 0.25pc prior -> K=4; MAP=true-fringe in 89% of plateau samples,
    P_true 0.57, exp H 3.4). The CENSUS-3 (J0711/J1713/J1909, certify at truth) SMEAR on the plateau
    (P_true 0.12/0.014/0.019, exp H 11-35) -- their fringe info is DATA-driven (needs joint sky), not
    prior-pinned. Clean split: prior-pinned (J0437, survives) vs data-driven (census-3, evaporates).
    Names J0437 as the concrete lever-(iii) anchor seed; does NOT change f (one pulsar can't phase
    the array).
  Files (CW_transition/): `trackB_R_referendum.py`, `trackB_R_znaddle.npz`,
  `trackB_R_zplateau_s{0,1}.npz`, `trackB_R_zplateau_summary.npz`, `trackB_R_referendum_result.npz`,
  `trackB_R_postmortem.npz`; logs/R_{gate,r1,r2,r3,postmortem}.log. Spec: "DELIVERABLE R" section
  (+ R POSTMORTEM) + terminal-verdict RESOLVED. NOT committed. (Claude + Matt)

- **2026-07-09 (B1 machinery + PILOT + STEP 1: the targeted (f,mc) ladder; cronus/4090, jug-gpu).**
  B1 was to test the CONDITIONAL (EM-counterpart-targeted) pipeline under real noise. Built the
  machinery, then a pilot re-scoped the deliverable. **MACHINERY (`trackB_b1_core.py`, 8/8 gates
  PASS):** (i) `MaskedDelay` puts a PER-PULSAR PULSAR-TERM MASK in the likelihood as a RUNTIME arg
  (delay = earth + m_p*(full-earth)), so ONE compiled graph serves the Earth-only fit, the
  phase-up fit with a certified subset's pterms on, and the full model -- this is what B1.3's
  loop-gain experiment needs; (ii) DATA is a runtime arg all the way into the low-rank split
  (`B1Split`; the banked `trackB_split.Split` BAKES data_obs into its jitted per-pulsar evaluators,
  which would have meant >=20 compiles); (iii) real noise draws -- heterogeneous white (the exact
  N diagonal, dug out of `vsm.Ns[p]` = WoodburyKernel_novar -> `.N.N`) + per-pulsar RN GP + the
  HD-correlated GWB GP; the linear timing-model GP is NOT drawn (marginalisation device, prior rms
  100 ns). GATES: masked lnL(mask=1)==stagec_fisher 0.00e+00; mask=0==EarthDelay 0.00e+00; split
  a_const invariant to (theta,data,mask) 2.9e-10; B1Split.lnL==masked logL 3.5e-10; B1Split.estep
  ==banked Split.estep **0.00e+00** (0/116 >1e-8); zero-noise E-step at truth == A2 census ceilings
  (0.9534/0.9887/0.9984); white draw var/N_diag 0.985; GWB draw corr matches the Phi-implied HD.
  lnL(truth|zero-noise)=405413.51 (== R). Noise is NOT a perturbation: median per-pulsar CW rms
  458 ns vs drawn noise 2005 ns (white 1414, RN 805). **DECISION with Matt: ARM A (L_true=L0, the
  convention every prior Track B deliverable inherited from `cw_helpers.py:228`) is a CONTINUITY
  row only; ARM B (n_true drawn prior-weighted over the SAMPLED fringe centres, within-fringe
  offset ~ U(-1/2,1/2)) is the HEADLINE -- the first Track B computation with truth off the prior
  mean.** Arm-B draw gated: 116/116 sampled in-window, |n_true| median 84 max 251.
  **PILOT (`trackB_b1_pilot.py`) -- the re-scope.** The brief's "sky fixed, FREQUENCY free (1-D
  scan), extrinsics fit" presumes one remaining registration axis. There are TWO. The model chirps
  the pulsar term (`evolve_phase(tp)`), so **log10_mc is a registration axis**; F2's ladder used
  the chirp-free phase, whose d/dlog10_mc is identically zero, so F2 was blind to it BY
  CONSTRUCTION (the pilot's exact-phase autodiff reproduces F2's chirp-free freq ladder to 3 digits
  -- 2.518e-2/1.035e-4 vs 2.52e-2/1.04e-4 -- so the disagreement is the chirp, not the code).
  Ladders (scaled, loosest/median): sky 1.85e-3/4.47e-5; log10_fgw 2.52e-2/1.44e-4; **log10_mc
  6.05e+1/7.84e-4**. Per-axis CERT tolerance (E-step, one axis at a time): sky 1e-5, freq 3e-5,
  **mc 1e-3**, extrinsics NO collapse to 1e-2 -> the 8 source params split cleanly into 4
  harmless extrinsics + 2 sky (counterpart) + 2 registration axes. Earth-term Fisher (HVP; NEVER
  pinv -- it reports sigma=0 for an unconstrained axis): info gain over the prior is 162-389x for
  log10_fgw and **1.00-1.73x for log10_mc** (sigma 0.50/0.86/0.87 dex; for loud1/2 the posterior IS
  the prior). Targeted wall: freq 203x, **mc 1730x**.
  **STEP 1 (`trackB_b1_ladder.py`) -- the cascade F2 pre-registered but never ran, now on (f,mc).**
  **95-100% of the pulsar-term phase wander is mc.** Registration radius R_a (box shrink to reach
  1 rad): max 2.14e-2 (J0711, 0.106 kpc), median 2.04e-4. **IGNITION: 0 pulsars** (loosest rung
  needs 47x, binding census pulsar J1713 needs 3305x). The ladder SPANS internally (42 rungs, max
  gap 0.387 dex < F2's 0.7) -- so unlike sky the obstruction is the mc BOX, not rung spacing.
  **LOOP GAIN** (conditional (f,mc) Fisher at a fixed fringe subset = a masked likelihood, so no
  optimiser and no L2c pull-in problem; Asimov injected with the SAME mask -> H=-Fisher exactly):
  0.04 at the float, 0.12 (top-1), 0.70 (top-24), **1.34 (top-48)**, **4.50 (census-3 only)**.
  A pre-registered structural hypothesis -- "loose rungs are loose BECAUSE they are blind to mc, so
  nothing bootstraps" -- was **REFUTED** (fixing J0711 alone cuts sigma_mc 7.8x; loose rungs carry
  mc info via SNR, not only via g_mc). The real structure is **BISTABILITY**: the registered state
  exists and strongly attracts (the P1 needle from the (f,mc) side) but is unreachable from the
  float. Certification is SELF-REFERENTIAL -- ~30 registered fringes are needed to measure the
  chirp mass that lets any fringe be registered. **CASCADE from each conditioning tier: NONE
  ignites.** A (sky only) R=0.0214; B (+EM period sigma_P/P=1e-3) R=0.0220 -- **an EM period buys
  NOTHING**, a 7x tighter f moves the best rung 3%; C (+host D_L -> sigma_mc 0.036/0.022/0.021 dex
  via log10 h = (5/3)log10 Mc + (2/3)log10 f - log10 D_L, i.e. set by the ARRAY's own
  sigma(log10_h)) R=**0.271**, missing ignition by 3.7x. (A pre-run prediction of 1.4x was WRONG --
  it used the MEDIAN mc shrink; R is a quadrature over the 3 loud sources and the loosest rung is
  dominated by loud0, the source with the worst sigma(log10_h).) DOC ITEMS recorded in the spec:
  (a) F2 chirp-blindness (its SKY ladder + the Track B terminal verdict stand; its freq ladder is
  0.72x optimistic); (b) **R retroactive note -- R fixed mc at truth and its cancellation argument
  is wrong for mc, but FAVORABLY: including mc shrinks Z_needle, so f=6.9e-7 STANDS and is FIRMER**;
  (c) **mc-registration IS the evolution/timestamp channel** -- the pulsar term is a kyr-baseline
  timestamp and what registers it is fdot, not f; the 22-yr Earth term cannot measure fdot
  (Delta_phi_E ~ 0.05 rad). That is the E-track's eccentric-harmonic mechanism in CIRCULAR form:
  **the E-track and the targeted pipeline have merged** -- lever (ii) is the missing ingredient
  lever (iii) needs, not an alternative to it. NEXT: STEP 2, the (f,mc) referendum at tiers A/B/C
  (`trackB_b1_referendum.py`); pre-registered, the loosest tier with f>=0.95 defines B1's legitimate
  conditioning. B1.0-1.5 NOT started. Files: trackB_b1_core.py, trackB_b1.py, trackB_b1_pilot.py,
  trackB_b1_ladder.py, trackB_b1_referendum.py; b1_pilot_m{1,2,3,4}.npz, b1_ladder.npz,
  b1_loopgain.npz, b1_cascade_tiers.npz; logs/b1_*.log. NOT committed. (Claude + Matt)

- **2026-07-09 (B1 STEP 2 + STEP 1D -- THE TARGETED REFERENDUM AND THE LAST DOOR; the targeted
  scenario has its OWN information wall, and it is PRICED. cronus/4090, jug-gpu).** STEP 1 closed
  with no cascade ignition at any tier. STEP 2 asked the evidence question R asked blindly, now with
  the sky supplied: does the fringe-marginalised posterior concentrate at truth over the two
  REMAINING registration axes (log10_fgw, log10_mc)? Object = R's count-once star-topology marginal
  verbatim (`lnL_marg(truth) = 405686.3434 = lnL(truth) + fringe entropy 272.83`, reproduced).
  Extrinsics fixed at truth is now MEASURED not assumed (pilot M2: census P(true) flat in
  inc/h/phase0/psi to 1e-2 scaled, so their Laplace factors cancel).
  **TIER C** (sky + EM period sigma_P/P=1e-3 + host D_L; mc box set by the ARRAY's own
  sigma(log10_h) via log10 h = (5/3)log10 Mc + (2/3)log10 f - log10 D_L), frozen at 4 seeds:
  ln Z_needle = 405629.6337 (bracketed quadrature, |dlnJ|=0.055 STABLE, 6/6 brackets closed);
  ln Z_box = 405633.035 (std 0.859, s.e.m. 0.429); d = **-3.4008 nat**; **f(C) = 0.0323 +- 0.0134**,
  2-sigma band [0.0055, 0.0591] -> **FAIL by 16.1x** (gap/scatter = 3.96). Break-even
  lambda_mc = 0.1206 -> deficit 8.29x **[SUSPENDED: Tier A's Z_box shows the plateau has SATURATED -- its box is 24000x larger in volume yet ln Z_box agrees with C's to within the 0.86-nat scatter, so Z_box is NOT proportional to the box volume and the break-even extrapolation is invalid. f(C) and the soft-cascade FAIL are unaffected; the deficit NUMBER must be re-measured on a shrunken mc box]**. ~97% of the targeted posterior's evidence sits on the
  wrong-fringe plateau even with the sky EXACT. **THE PRICE HAS TWO TERMS**: sigma(log10_mc)
  delivered 0.0364/0.0217/0.0205 dex vs needed 0.0044/0.0026/0.0025; but setting sigma_h -> 0 leaves
  a FLOOR of 0.00301 dex set by the assumed 1% host distance -- BELOW loud0's need, ABOVE loud1's and
  loud2's by 1.15x/1.22x. So strain alone cannot close it: need sigma(log10_h) x11.3 (T x128 ~
  **2840 yr**, or ~1 dex louder) **AND** sigma(log10 D_L) <= 1.0% -- the counterpart is CO-BINDING.
  Or an eccentric-harmonic **kappa >= 8.29** substituting for both (the e-value is the E-TRACK's to
  measure; not quoted here).
  **STEP 1D SOFT-CASCADE PROBE (the last door; pre-registered, one probe one verdict).** R_a>=1 is
  the HARD lock criterion; a soft posterior-weighted mixture at R~0.27 spreads over ~4 fringes and
  might still leak Mc. Tier-C conditioning, Asimov, 5 iters, NO hard fixes; Mc update = width of the
  mixture-marginalised posterior; truth scores the wrong-fix column only. Result: per-iter sigma_mc
  shrink [0.057, 0.306, 3.335, 1.428, 1.153] -- cumulative 3.77x but **NON-MONOTONE**; false-fringe
  mass W = 114.8 -> 113.1 while the number of pulsars whose W GREW climbs **54 -> 70** (the soft
  analogue of the GPS wrong-fix). **VERDICT FAIL: THE DOOR IS CHECKED AND CLOSED.** Mechanism:
  sigma_mc is ALREADY 1e-4..1e-3 dex at iter 0, far narrower than the 0.0026-0.0044 needed -- the
  local mode was never the problem; median q_max is flat 0.067->0.070 (~16 effective fringes/pulsar,
  unchanged) and S flat +12.5->+12.8. The missing information is not local curvature but WHICH of ~16
  fringes, and that choice never sharpens. Local width != global concentration.
  **FRONTIER STATEMENT:** the pulsar term is a kyr-baseline TIMESTAMP; it cannot be read without the
  CLOCK RATE; the clock rate is fdot, i.e. Mc; a 22-yr Earth term cannot measure fdot (Delta_phi_E ~
  0.05 rad, info gain over prior 1.00-1.73x). **Design-theorem lever (ii) is not an alternative to
  lever (iii) -- it is the ingredient lever (iii) is missing.** The E-track's eccentric map is the
  PRICED next experiment: it must deliver kappa >= 8.29.
  **R UNAFFECTED (one line):** the micro-dip and the (f,mc) contest live in dims R held at truth;
  R's sky-plane Z_needle had 6/6 positive curvature and quad/Laplace agreement at 0.2 nat --
  **f = 6.9e-7 stands, and FIRMER** (including mc shrinks Z_needle).
  **NEEDLE = a thin SHELL, not a point.** Hessian of lnL_marg at truth over (f,mc): 5/6 positive,
  one NEGATIVE (-1.32e9), Richardson-stable (2.4e-2) at steps 0.3x the sharpest eigen-sigma. Not a
  saddle (lnL_marg falls 159 nat at 1x base, 919 at the box edge along that eigenvector) but a
  MICRO-DIP at truth: lnL_marg = lnL_ref + sum_p m_p, m_p >= 0 grows as a pulsar de-registers, so a
  sub-fringe offset buys more entropy than it costs; the local max sits ~1e-5 scaled away, 0.12 nat
  higher. Z_needle is therefore quadrature-only, NEVER Laplace, and negative eigenvalues must never
  be clipped (clipping inflates Z_needle -- biased toward letting B1 proceed).
  **METHOD FAILURES RECORDED (both self-defeating gates, spec CONVENTION section).** (1) `acc>=0.25`
  gated while the RW scale adapted toward 0.234 -> every high-beta stage exhausted max_mcmc
  (sweeps=14, acc=0.24) and the remedy (add seeds) cannot change a kernel property; fixed by adapting
  toward 0.35 with a startup assert -- the fixed sampler is both correct AND cheaper (same stage:
  sweeps=3, acc=0.29). (2) seed "spread" defined as max-min is a RANGE, which GROWS with seeds, so
  "add seeds until spread<=0.3" is defeated by its own remedy (0.920 -> 2.028 at seed 2). Tier C's
  SMC gate is therefore recorded **FAILED-AS-SPECIFIED**, while the verdict is **SEPARATELY FIRM**
  (0.86-nat scatter vs a 3.40-nat gap). New convention: mixing gates = std at a FIXED pre-registered
  seed count; s.e.m. for the verdict; never a remedy that moves the gated statistic the wrong way.
  Also caught before use: an unresolved-peak quadrature returning its own grid spacing (sum lnJ fell
  by exactly -ln2 per doubling) -- fixed by bracketing the -30-nat contour from the FUNCTION, not the
  curvature; and a first soft-cascade run whose sigma_mc came out 0/nan because inv(Hpos+Pi) is not
  PD and because 3e-6 steps measure the razor, not the marginal width (replaced by profile
  half-widths). Files: trackB_b1_referendum.py, trackB_b1_softcascade.py, b1_step2_table.py;
  b1_referendum_tierC.npz, b1_softcascade.npz; logs/b1_ref_tier*.log, logs/b1_softcascade.log.
  B1.0-1.5 and the E-track map NOT started. NOT committed. (Claude + Matt)

- **2026-07-10 (B1 STEP 2 CLOSE-OUT: three-tier table + break-even RESPONSE CURVE; the suspended
  deficit was too OPTIMISTIC. cronus/4090, jug-gpu).** Completed Tiers A and B on the frozen 4-seed
  protocol and measured the break-even as a CURVE rather than extrapolating a point.
  **THREE-TIER TABLE:** f(A)=0.0847+-0.131 (sky only), f(B)=0.0481+-0.0227 (+EM period),
  f(C)=0.0431+-0.0369 (+host D_L); ln Z_needle tier-independent to **0.003 nat** across three boxes
  and two independent bracket algorithms. All three gates FAILED on the range statistic; all three
  verdicts FAIL on the +-2sigma band by 4.4x-13x. **The tier gradient is FLAT and mildly INVERTED**
  (f falls as conditioning tightens): Z_box(A) >= Z_box(C) is required by integration over a 24000x
  larger volume, and is measured at -1.02 +- 0.95 nat, i.e. EQUAL -- which IS saturation.
  **BREAK-EVEN CURVE** (lambda_mc = 1/0.3/0.12/0.05, 2 seeds each, needle excision 0.0% everywhere):
  ln Z_box = 405633.035 / 405631.754 / 405630.535 / 405628.910; f = 0.032/0.107/0.289/0.673.
  **(a) SATURATION SCALE: the plateau's own chirp-mass extent is ~0.02 dex** -- a newly measured
  quantity. **(b) lambda* < 0.05 -> a BOUND.** **(c) CORRECTED DEFICIT > 20x**, replacing the
  suspended 8.29x -- **the suspension was vindicated and the real price is WORSE than the
  proportionality implied**; even a 20x-shrunken mc box only reaches f=0.673. (A log-linear
  extrapolation gives ~66x; NOT reported -- that is the very error the curve exists to correct.)
  **CONSEQUENCE:** sigma(log10_mc) must improve >20x, i.e. below ~0.003 dex, while the sigma_h -> 0
  floor set by a 1% host distance is 0.00301 dex: **no combination of (position, period, host
  distance) delivers it.** The two-term price (sigma_h x11.3 AND D_L <= 1.0%) is SUPERSEDED.
  WEAVE's kappa>=8.3 / e>~0.6 threshold INHERITS the suspension -- mechanism stands, threshold
  re-prices above the >20x bound; the e-value is the E-track's to compute.
  **HEADLINE (mechanism, not number): conditioning the (f,mc) PRIOR BOXES barely moves the evidence
  because the plateau does not fill them; what moves it is LIKELIHOOD STRUCTURE, which the eccentric
  harmonic comb supplies and no prior box can.** Method conventions added to the spec: (i) mixing
  gates = std at a FIXED seed count, never a range (max-min GROWS with seeds: 0.920 -> 2.028 -> 4.087
  observed); (ii) a break-even is a RESPONSE CURVE, never a point through an assumed scaling;
  (iii) logically-redundant arms retain AUDIT value -- Tier A, twice argued "foreclosed", falsified
  the uniform-density assumption the whole price rested on, and saturation is invisible from one box
  BY CONSTRUCTION; (iv) SMC wall time is ~98% sweeps x (one G call) -- count sweeps before blaming
  the device. Files: trackB_b1_breakeven.py, b1_breakeven_curve.npz, b1_step2_table.npz,
  b1_referendum_tier{A,B,C}.npz; logs/b1_ref_tier{A,B}.log, logs/b1_breakeven.log,
  logs/b1_step2_table.log. Finiteness verified on VALUES (the only NaNs are BY DESIGN: lnZn_lap is
  refused when curvature is non-PD; loopgain's next-rung R is undefined when no free pulsars remain).
  B1.0-1.5 and the E-track map NOT started. NOT committed. (Claude + Matt)

## 9. Environment (cronus)
- GPU box cronus = NVIDIA RTX 4090, driver 550.120.
- Toy/Fisher work (prong2_transition.py): isolated venv with jax 0.4.28
  (CUDA 12.4 / cuDNN 8.9), driver-compatible. Latest jax[cuda12] FAILS cuDNN init
  on this driver. Activate: source <scratch>/env/bin/activate.
- discovery lives in the shared `discotech` env (separate jax) — DO NOT modify it.
- OPEN RISK for Stage C: confirm discotech's jax can init the GPU; if it shares the
  broken jax, Stage C Hessian is CPU-bound or needs a pinned discovery venv.