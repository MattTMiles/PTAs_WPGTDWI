# STORY — the programme's account, as a numbered-claim skeleton

**Status:** tracked master document. Single-writer rule: **cronus doc sessions only.** Every future
campaign writes into this file; fenced agents never do.

**What this is.** The comprehensive architecture of the pulsar-distance programme, stated as numbered
claims. Nothing here is newly derived — every number is transcribed from a banked source. This is
scaffolding for the paper(s), not analysis.

**Tag discipline (binding).**

- `[MEASURED: source §section — numbers inline]` — the claim is banked. Its source is named. Its
  numbers travel with it.
- `[PENDING: CAMPAIGN — what it will supply]` — the claim is a hole with a named owner.
- `[DISPUTED: …]` — two banked sources disagree, or a claim's status is genuinely ambiguous. **Flagged,
  not resolved.** A future session resolves these with a measurement, never with a preference.
- `[SUPERSEDED → …]` — the claim was true as measured and is now scoped, retracted, or replaced. **The
  trail is kept.** Nothing is silently deleted; the sequence of wrong answers is part of the result.

**Standing conventions that govern every number below** (full text in
`CW_transition/trackB_estimator_spec.md`):

- **Every quoted onset carries its floor's `N` and fit error, or it is not quoted.**
- **A calibrated threshold states its false-alarm rate α and its sampling scatter.** An order statistic
  is not a threshold.
- **The D2 Gumbel floor is valid only where the null's zero-fraction is ≲ 20 %.** Above that, quote the
  empirical (1−α) quantile with a **bootstrap** error. **The zero-fraction is a REQUIRED column, not a
  caveat.** *(criterion-v2.2, adopted from ANCHOR §4 — a Gumbel fitted to a point mass at zero errs
  PERMISSIVE.)*
- **Confidence without a detection statistic is prior-pinning in disguise.**
- **Count and correctness are quoted together.** An above-onset count without its purity number is
  meaningless.
- **Bank the statistic, not the verdict.**
- **Bank the ORIENTATION, not just the array.** An index column and a sentence stating what pairs with
  what — the alternative is a silent transpose. *(RECUT §3/§4.)*
- **A difference or ratio between two cells is RE-DERIVED, never carried forward**, when either cell's
  floor has moved. *(RECUT §6 — two such claims were checked after the floor fix; both moved.)*
- **A break-even is a response curve, never a point extrapolated through an assumed scaling.**
- **A statistic evaluated against truth is an ORACLE until its implementable form is scored.**

**Primary sources.** `project_progress.md` (canonical tracker, §10 = the ACCRE consolidation);
`CW_transition/trackB_estimator_spec.md` (the criterion + Track B); `reports/` (GEO, RING, SIREN,
ATLAS, FORGE, IGNITE, IGNITE2, D4 + their banks); `CW_transition/{ABACUS,SCOUT,LANES,WEAVE,PENCIL,
GALLERY}*.md`.

---

## S1 — THE PROBLEM AND THE REFRAME

### S1.1 The fringe problem

**S1.1.1** The phase-connected CW likelihood is **periodic in each pulsar distance**, with fringe
spacing `dL = c / [f_gw (1 − cos μ)] = λ_gw / (1 − cos μ)`. A single CW leaves the distance degenerate
across fringes; multiple CWs with different `dL` per source break the degeneracy in principle, because
only the true distance phases all sources simultaneously.
[MEASURED: project_progress §1 — fiducial `dL ≈ 0.36–0.49 pc` for the best-timed pulsars]

**S1.1.2** The pulsar-term phase is a **lever arm**: `Φ_p = 2π L_p / dL_p ≈ 1.6 × 10⁴ rad`. Every
result in this programme — the certification tolerance, the needle width, the mas-class localisation
equivalence, the unreachable basin — is this one number wearing different clothes.
[MEASURED: spec "B0.2", trackB L2c — cusp curvature ~6e11, Hessian condition number ~4e10]

**S1.1.3** For a single pulsar, distance is **perfectly degenerate** with each source's pulsar-term
phase: distance is only ever an *array* measurement, in which the shared global Earth term pins the
source phases. A pulsar marginalising on its own data block alone retains `marginal = 8.2e-10 × cond ≈ 0`;
the full-array marginal is `0.93 × cond`.
[MEASURED: project_progress §6 Stage A.1 — asserted, not assumed]

### S1.2 The reframe: located ≠ certified

**S1.2.1** **Within-fringe precision and fringe identification are different questions, and the old
metric conflated them.** Marginal within-fringe `σ_L` is excellent and cheap: **116/116 pulsars** have
`σ_L` below their EM prior at `N_CW = 1` (median 0.11–0.55 pc; best ~0.10 pc; EM priors span 1–290 pc).
This is **not** a distance measurement. It presumes the correct fringe.
[MEASURED: project_progress Stage C D2 — `marg/cond` median 1.008 (16/84: 1.001/1.033), rcond spread 2.7e-2]

**S1.2.2** **A single realistic CW identifies ZERO fringes across a 116-pulsar array.** The true-vs-best-wrong
fringe gap is `dlnL_a`: median **0.0000**, 84th pct 0.0002, max 0.0028 — against `ln K_a` median **7.94**
(range 2.48–11.05). The gap is ~1000× below the trials threshold. Class-(i) [fringe identified AND
`σ_L` < EM prior] = **0/116**.
[MEASURED: project_progress Stage C D3, 10 Asimov draws]

**S1.2.3** **LANGUAGE RULE, enforced everywhere since.** No "N distances measured/improved" without the
class-(i) qualifier. D2's 116/116 is "within-fringe curvature below EM prior" — necessary, not sufficient.
[MEASURED: project_progress Stage C, D3 close]

**S1.2.4** The metric that produced the old optimism was **buggy in a way that inverted a verdict**.
`find_best_wrong_mode_in_prior`'s `find_peaks` prominence floor (0.5) exceeds the *entire* single-CW
fringe modulation, so it returned "no competing peak" — which naively reads as "perfectly identified",
exactly backwards. Corrected with a prominence-free direct scan at fringe centres.
[MEASURED: project_progress Stage C D3 method note]
[SUPERSEDED → the conditional joint-`ΔlnL` metric is retired; the corrected fine-grid scanner is the
standing tool. Blast radius: nb-05 corrected 7–45 % more negative at every `N_CW` = 1–8 (sign unchanged);
**nb-01 at h = −13 FLIPS SIGN: +296.09 → −108.35** — "separable" → "not separable".]

**S1.2.5** The old prong-1 successes (20/20 distances recovered) ran at **SNR ≈ 600 per source**
(`h = −12`), **56.2× louder** (= 10^1.75 exactly) than the realistic D4/nb-05 regime (SNR ≈ 10.7 at
`h = −13.75`). The successes were real; the regime was not.
[MEASURED: project_progress P2-D item 1]

### S1.3 The three prongs: as posed → as answered

**S1.3.1 Prong 1 (can distances be optimised out?) — as posed by Hogg: yes, without sampling.
As answered: NO from a cold start, and the NO is information-theoretic, not merely a search failure.**
[MEASURED: spec Track B TERMINAL VERDICT + DELIVERABLE R — `f = Z_needle/(Z_needle+Z_plateau) = 6.9e-7`]

**S1.3.2 Prong 2 (the GWB↔CW transition) — as posed by Farr: where does distance information die?
As answered: CLOSED as a computation.** Two distinct transitions were mapped (confusion, coherence),
shown to compound, and then shown to be *the wrong binding constraint* for real data — the real
likelihood never reaches the confusion knee at achievable `N_CW`; **fringe identification binds first**.
One handoff remains open (S7.6).
[MEASURED: project_progress "PRONG-2 CLOSURE" 2026-07-02]

**S1.3.3 Prong 3 (sampling distances when a CW is present) — as posed by Taylor: better proposals.
As answered: the honest posterior does not concentrate at truth, so prong 3 cannot win on the data
alone.** The referendum is not a statement about samplers; it is a statement about the evidence a
sampler would be asked to find.
[MEASURED: spec DELIVERABLE R — R4 verdict (ii)]
[PENDING: SAMPLER — the actual posterior machinery + SBC; see S10]

**S1.3.4** Farr's ordering — *"there is provably no point doing 2 before 1, nor 3 before 2"* — was
followed, and it **held**: the information answer (2) is what closed (1) and pre-empted (3).
[MEASURED: project_progress §3, collaborator positions]

---

## S2 — THE ANATOMY

### S2.1 The confusion knee, and its calibration against Boyle & Pen

**S2.1.1** Marginal/conditional distance information collapses as sources crowd. Array-median
`marg/cond`: **0.99 (N=1) → 0.87 (N≈100) → 0.53 (N≈400) → 0.08 (N=1000)**; 0.5-knee at **N ≈ 410**
for the fiducial (116 psr, 15 yr, 3–20 nHz). Conditional information stays ∝ N to N = 1000 (Mihir's
independence line); the marginal rises, turns over, collapses.
[MEASURED: project_progress §6 Stage A, 30 seeds]

**S2.1.2** The knee tracks `N* = T·Δf` with a large **array-resolution boost**. The boost is **linear
in `N_psr` with no saturation through 200 pulsars**: `knee/N* = 0.40 · N_psr^1.03` (fit A = 0.3990 ±
0.008, b = 1.0263 ± 0.005, RMS 0.007 dex). Forecasting law: **`N_knee ≈ 0.40 · N_psr · T · Δf`**.
[MEASURED: project_progress P2-C — knee/N* = 6.47/12.92/26.73/52.75/91.45 at N_psr = 15/30/60/116/200;
gate: 116 → 52.75 reproduces Stage A's ~52]

**S2.1.3** **The measured law CALIBRATES Boyle & Pen (2012).** Identifying `T·Δf ≡ F` (their bin count),
the two laws share the form `∝ N_psr · F`. Our coefficient **0.40 = 2/5** is their **2/7** with their two
source-sky-position parameters removed (7 → 5 parameters per source), because our metric *supplies* the
source sky and asks only whether pulsar distances survive. **`2/5 ÷ 2/7 = 7/5 = 1.400`, matched to <1 %.**
The residual ~14 % in the raw count at 116 psr is the soft (half-information) vs hard (full-deconfusion)
threshold definition plus finite-F `(1 − 1/2F)`.
[MEASURED: ABACUS §3–§5 — RECONCILED; ledger: BP known-dist 266.8, BP unknown-dist 250.2, mapped
known-sky 373.5, measured knee 424.6]

**S2.1.4** **Mode-independence of `marg/cond` is an equal-amplitude artifact.** Fixed-per-source and
fixed-total power give identical `marg/cond` to machine precision (max |Δ| = 0.000 at every N). But under
a realistic luminosity function (3 loud + 13 faint, h ratio 10) the *conditional* information plateaus
(2.8e10 vs 2.1e12 at N = 1000) because the faint sources carry h²/100 — **the curve is set by the handful
of loud sources.** This lesson recurs, in real units, as S3.2.4.
[MEASURED: project_progress Stage A + A.1 — first |Δ| > 0.02 at N = 158, max Δ = −0.025]

### S2.2 The coherence law and the wander floor

**S2.2.1** **Model-independent law: distance information halves when the rms Earth–pulsar phase wander
reaches ~1/SNR rad.** Closed form at N = 1: `marg/cond = 1/(1 + SNR²·σ_φ²)`. Verified across two
functional forms (linear and saturating), which cross at the same `SNR²σ_φ² ≈ 1`; the *location* is the
invariant, the coefficient and shape are model-dependent (0.5-crossing moves 0.0144 → 0.0044, a 0.30×
shift). **Absolute `t_c` values must never be quoted as physical.**
[MEASURED: project_progress Stage A.2, gates (i) 3.4e-6 and (ii)]

**S2.2.2** **The wander floor: beyond SNR ~ 1/σ_φ, more strain does NOT improve distances.** Two regimes:
`σ_L → dL/(2π·SNR)` (coherent, SNR-limited) and `σ_L → dL·σ_φ/(2π)` (decoherent, **SNR-independent
floor**). This is the rigorous form of "a strong CW can phase up the array" — it cannot, past the floor.
[MEASURED: project_progress §6, "two distinct transitions" B]

**S2.2.3** A pre-registered guess was **falsified and is recorded as such**: decoherence is NOT
T-independent. Coherence `marg/cond` falls as `1/(1 + T/16 yr)` (0.62/0.52/0.45/0.35 at T = 10/15/20/30 yr)
because the constant per-baseline offset is increasingly *resolved* as SNR² ∝ T. Both transitions scale
with T; **the clean discriminator is the N-axis, not T-scaling.**
[MEASURED: project_progress Stage A.2 Exp 2]

**S2.2.4** **Confusion × coherence COMPOUND — they do not factorise.** The R = 0.5 contour is *curved*
(bows toward the origin), and R falls below the factorised product through the confusion-dominated
interior. A realistic SMBHB background sits in the compounding region, so the practical threshold is
**worse** than the independent-threshold estimate.
[MEASURED: project_progress Stage A.3 — peak factorisation residual 0.202 under commuting reductions;
sign structure: weak positive (de-confusion) lobe at N/N* < 1, dominant negative lobe at N/N* ≳ 3]
[SUPERSEDED → the first-pass residual 0.331 used median-of-ratios, which does not commute with the
product; the commuting recomputation gives 0.202. **The verdict (COMPOUND) is unchanged; the magnitude
was inflated ~50 %.**]

**S2.2.5** **The pulsar-side coherence axis is NOT an independent limit — it IS the 1/SNR measurement
floor.** Additive red/DM noise perturbs the pulsar-term phase by the noise-to-signal ratio,
`σ_φ = σ_res/A_CW = 2πf σ_res/h = 1/SNR_pterm`. Array median `σ_φ = 2.13 rad` (16/84: 0.66/8.14), median
`SNR_pterm = 0.47` — i.e. `SNR_pterm · σ_φ ≈ 1` **by construction**: the array already sits at the knee.
[MEASURED: project_progress P2-B item 2 (corrected)]
[SUPERSEDED → P2-B v1 computed `σ_φ = 2πf σ_res` (a time-base *shift*) = ~6e-15 rad. **Wrong projection,
by a factor 1/h = 5.6e13 (~13.8 orders).** Stated for the record; the corrected chain is the claim.]

**S2.2.6** The genuinely independent coherence knob is **source-side stochastic `df/f`**, which crosses
the knee at `df/f = 3.4e-5 (SNR 5) / 8.4e-6 (SNR 20) / 1.7e-6 (SNR 100)` over kyr lags. Whether real
SMBHB environments (gas/stellar coupling) reach that is **an astrophysics question this project's methods
cannot answer.**
[MEASURED: project_progress P2-B item 1a — HANDOFF (Taylor/Farr), open]

### S2.3 The readable sub-array

**S2.3.1** **Only ~a quarter of the array can hear its own pulsar term.** Of 116 pulsars, **30** have
`SNR_pterm > 1` at realistic strain (σ_φ below the 1-rad knee). Leaders: J0613−0200 (3.04), J0751+1807
(2.53), J1045−4509 (2.42), J1012+5307 (2.35). Named chain: J1713+0747 (σ_φ 1.26 / SNR 0.79),
J0437−4715 (0.74 / 1.35), J1909−3744 (7.64 / 0.13).
[MEASURED: GALLERY §1.5 / P2-B — array median SNR_pterm 0.47]

**S2.3.2** The **timing model + frozen GPs absorb low-frequency pulsar-term power**, costing a factor
**~1.6× in SNR** (~2.6× in information). Two independent computations agree to <1 %: D3's pulsar-term
SNR cross-check (ratios 0.81–0.95 against the coherent law) and P2-A's Earth-term gate (23.49 vs 36.80
white matched-filter = 0.638 = 1/1.567). **The naive white formula is the optimistic one.**
[MEASURED: project_progress Stage C D3 + P2-A gate]

### S2.4 Detectable ≫ rangeable

**S2.4.1** **The sky becomes VISIBLE almost immediately and RANGEABLE almost never.** Detectable
(Earth-term SNR > 5) vs rangeable (class-i), median over 10 Asimov draws:

| N_CW | detectable (SNR>5) | rangeable (class-i) |
|---|---|---|
| 1 | 1/1 | 0 |
| 2 | 1/2 | 0 |
| 4 | 3/4 | 0 |
| 8 | 7/8 | 0 |
| 16 | 13/16 | **5** (range 2–10 of 116) |
| **pop (3 loud + 13 faint)** | **3** | **0** |

**For the real SMBHB sky, a 116-pulsar array is a DETECTOR, not a RANGER.**
[MEASURED: project_progress P2-A]

**S2.4.2** Class-i switches on **gradually and tail-driven**, not sharply. The *median* pulsar is never
identifiable even at `N_CW = 16` (`dlnL_a` median 0.95 ≪ `ln K` 8.8); it is the **tail** that crosses
(max `dlnL_a` grows 0.1 → 4.2 → 5.1 → 7.9 → **28.4** across N_CW = 1→16, roughly one order of magnitude
per doubling).
[MEASURED: project_progress Stage C D4]

**S2.4.3** **`marg/cond` stays ≈ 1 (median 1.01 → 1.04) through N_CW ≤ 16.** The real likelihood is deep
in the *resolvable* regime, hundreds of sources below the toy confusion knee. **The real likelihood does
not reproduce the toy's turnover in the achievable range** — what it adds is the discrete fringe-ID layer
the toy Fisher never saw, and **that**, not confusion, is the binding constraint.
[MEASURED: project_progress Stage C D4 — the single most important structural correction of the toy phase]

**S2.4.4** Frozen-GWB optimism is **small and bounded**: ±0.5 dex in the frozen GWB log-amplitude moves
`σ_L` by ≤ 2 % and the fringe-ID gap by ≤ 9 %, and does not change the class-i count. **The dominant
realism levers are the source luminosity function and real-noise-vs-Asimov — not un-freezing the GWB.**
[MEASURED: project_progress Stage C D5, N_CW = 8, 10 draws]

### S2.5 Track A closure — the anisotropic covariance channel

**S2.5.1** **Hogg's contention tested and answered: outcome (ii) — the channel is real but formal; Farr's
null effectively holds at realistic N.** The anisotropic enhancement in the cross-pulsar covariance obeys
a clean power law `I_LL_aniso ∝ N^{−α}`, **α = 1.17 [1.09, 1.27]** (from shot-noise anisotropy,
`c_lm ∝ 1/√N`) — *slower* than the resolved channel's `N^{−2}` past the knee, a formal rate-crossover.
[MEASURED: project_progress Track A F2]

**S2.5.2** **It is not a useful crossover.** By N ≈ 1000–3000 the enhancement is **2–6 %** of the isotropic
baseline and its 16–84 % band already includes zero; the implied `σ_L` is **77–96 pc for the best pulsar
(J0437)** against the resolved channel's ~2 mpc — **4–5 orders worse** — and it is **WITHIN-FRINGE**, so
it inherits the same fringe-identification problem that S3.1 shows fails on this array. **Anisotropy does
not rescue the confusion-regime route.**
[MEASURED: project_progress Track A F1/F2 + amendment 5]

**S2.5.3** Gates that make the verdict load-bearing: the enhancement is **exactly zero** off-anisotropy
(monopole-only shot-noise sky reproduces the iso baseline to **0.00e+00**, identical code path); the
isotropic cross-pulsar distance derivative is **not** zero (it is the finite-distance HD correction,
suppressed ~`ftau0^{−1.2}`) and is subtracted as a baseline, never asserted away.
[MEASURED: project_progress Track A gates G1, F0]

**S2.5.4** Caveat that travels: Whittle bins 2–4 ran **CAPPED** at 5.8/3.9/2.9·ftau0 (vs the strict 10×)
and carry **52 %** of the iso Fisher, so absolute `I_LL`/`σ_L` carry a resolution uncertainty for the
few L ≳ 4 kpc pulsars (worst pulsar's Fisher row moves ~100 % under a doubling). **The α-scaling and the
verdict are unaffected** (they are set by `c_lm ∝ 1/√N`, not by resolution).
[MEASURED: project_progress Track A caveat (a)]

---

## S3 — THE CENSUS AND THE GEOMETRY

### S3.1 Anchors do not exist

**S3.1.1** **The strict anchor set (K ≤ 1 — certifiable by the EM prior alone) is EMPTY.** Zero pulsars,
under *every* prior column including hand-injected published composites. **The prior-alone-anchor
assumption of arXiv:2603.28897 FAILS on the real array.**
[MEASURED: project_progress Anchor Census A2 — PRIOR-CERTIFIED (K<1) = 0 at every N_CW, every prior]

**S3.1.2** Counts of 116, under both canonical and min-σ-optimal columns: **K ≤ 3 = 0; K ≤ 10 = 0;
K ≤ 30 = 2** (J0437−4715, J2222−0137). Smallest K anywhere = **J0437, K_opt = 11.88.** The named K ≤ 10
list is **empty**.
[MEASURED: Anchor Census A1 — canonical prior file `CW_transition/best_distances.txt` (frozen, 116 rows,
69 parallax-class + 47 feather-kept)]

**S3.1.3** **Only ONE pulsar reaches K ≤ 3 at all, and only via a literature composite the canonical
script cannot produce.** J0437 at Reardon+2016's `156.79 ± 0.25 pc` gives `K = 6·0.25/0.489 = 3.07`. The
catalogue's tightest J0437 is 0.97 pc (K = 11.88); the 0.25 pc value is a *joint-timing composite*
(annual + orbital parallax + Pbdot + kinematics), not reproducible by any single method.
[MEASURED: Anchor Census A1 gate outcome + `stagec_anchor_a2.py::LIT_INJECT`]

**S3.1.4** **Real priors are not uniformly tighter than the feather priors — some are WIDER.** J1713+0747:
real σ = 40–61 pc vs feather 19 pc (`get_distance` prefers a VLBI technique over a more precise timing
parallax), so `K_lit = 536–816` vs `K_feather = 255`. **The tighten-hypothesis holds for a few
(J2222 R = 0.018, J2241 R = 0.15) and NOT universally.**
[MEASURED: Anchor Census A0/A1 — the J1713 reversal, hand-checked in both directions]

**S3.1.5** **The composite prior buys, at most, one marginal geometry-dependent anchor.** J0437's first
**data-driven** (median-geometry, CW-carried) certification is at `N_CW = 8` under the literature prior
(bayes 5/10) vs `N_CW = 16` under feather. Its apparent `N_CW = 1–4` certifications are **geometry-wide
artifacts** — favourable large-`dL` draws (including a `dL = 12.3 pc` antipodal draw), which the cause-tag
audit exists to catch. At median geometry (`dL = 0.489 pc`, σ = 0.25 pc) the prior-only `P_true = 0.77 < 0.9`:
**J0437 cannot certify on median geometry at N_CW = 1.**
[MEASURED: Anchor Census A2 — the J0437 onset triple, cause-tag-audited]

**S3.1.6** **The real seed set is the DATA-DRIVEN certifications, not the prior anchors.** Under the
Bayesian bar, the realistic population certifies **3** — J0711−6830, J1713+0747, J1909−3744 — all
loud-source-broken, all **prior-INDEPENDENT** (they certify under the feather prior too). Under the
conservative flat (`ln K`) bar it certifies **ZERO**, at every prior including composites.
[MEASURED: Anchor Census A2 verdict (ii)+(iii)]
[SUPERSEDED → both of these numbers are pre-criterion. See S3.2 (the count is a *distribution*, not a
number), S5.2 (the Bayesian bar is not a detection statistic at all), and the trail in S3.2.5.]

**S3.1.7** **Exactly one pulsar's fringe information survives the collapse of the joint sky, and it is the
prior-pinned one.** On the R plateau (joint `f = 7e-7`, i.e. registration has failed), **J0437−4715**
still MAPs to the true fringe in **89 %** of samples (`P_true = 0.57` against its needle value 0.68,
`exp H = 3.4`). The census-3 **evaporate** (`P_true` 0.117/0.014/0.019, `exp H` = 11–35, 37–68 distinct
wrong fringes). **Clean split: prior-pinned fringe information survives joint failure; data-driven fringe
information does not.**
[MEASURED: spec R POSTMORTEM — the sharpest single statement of the anchor/census dichotomy]

### S3.2 The certified-count distribution

**S3.2.1** **The census's count is a draw from a distribution, not a property of the array.** Across
**40 independent isotropic sky redraws** at fixed population loudness, the Bayesian count is
**4.50 ± 1.48** (median 5, range **1–9**); strict (P > 0.99) **1.57 ± 0.98** (range 0–4, with 5/40 draws
certifying nothing); flat (`dlnL > ln K`) **1.38** (range 0–4). **The census's single draw (3) sits at the
25th percentile.** *"Is it 3 ± 1 or 0–6?" is answered: neither.*
[MEASURED: GEO §2 — 40 draws, zero-noise Asimov, conditional-at-truth ceilings, literature priors]

**S3.2.2** **The trail of the count, and it is the spine of the whole criterion story:**

| stage | count | what changed |
|---|---|---|
| census (one sky draw) | **3** | the original number |
| GEO Bayesian (40 skies) | **4.50 ± 1.48** (1–9) | the sky is a random variable |
| two-layer `dlnL > ln K` ⊕ `q > 0.9` (FORGE §9) | **1.35 ± 0.82** | a *detection* statistic added; ~70 % of the Bayesian count was prior/trials-driven |
| **criterion-v1 three-layer** (+ absolute floor 9.01 nat) | **0.275/draw** | the null could still fire at the two-layer gate; the floor zeroes it. The Bayesian 4.5 was **~94 %** prior/trials-driven |
| criterion-v2 (floor is a *function*) | *(not re-cut at zero noise)* | 9.01 is now known to be the census-loudness value of `floor ∝ h^1.66`; GEO's 0.275 stands **as a census-loudness number** with its floor's ±5-nat scatter attached |

| **criterion-v2.1 ZERO-NOISE CEILING** (flat gate, layers 1 + 3, **no floor term**) | **1.350 ± 0.82/draw** | **the floor is RETIRED at zero noise, as a category error** — see below |

[MEASURED: GEO §2; FORGE §9.2; project_progress §10.0 THE FINAL TABLE; IGNITE §4.1; SURFACE §8]
[~~DISPUTED (D-3)~~ → **CLOSED by SURFACE §8.** At zero noise there are **no fluctuations for a
noise-floor to gate**: GEO's data are Asimov-at-truth, and the 9.01-nat floor was fitted to *noisy*
nulls. Applying it to a noiseless statistic is **not a mis-sized number — it is a CATEGORY ERROR.**
The zero-noise ceiling is therefore quoted under the **flat gate alone** (layers 1 + 3) and is
**1.350 ± 0.82/draw**, reproducing FORGE §9.2's independent two-layer number (1.35 ± 0.82) to three
digits. ***`0.275` is RETIRED as a floor-concept category error, not as a wrong measurement. Its
sign was always safe; its value was never meaningful.*** And the caveat is **closed, not
inherited**: **0 of 4640** (draw, pulsar) cells have `dlnL < 0`, so at zero noise the MAP fringe
**is** the true fringe and `q_max ≡ P_true` identically — **GEO's count is IMPLEMENTABLE, not an
oracle.** The D4 oracle/implementable gap does **not** bite at zero noise.]

**S3.2.3** **Surviving detectors under criterion-v1, per pulsar:** GEO zero-noise — J1909−3744 (0.225),
J0437−4715 (0.025), J1713+0747 (0.025). FORGE Arm A — J1909−3744 (0.067). **FORGE Arm B — none.**
[MEASURED: project_progress §10.0]

**S3.2.4** **THE HEADLINE, and it is a hard one.** Under a detection gate calibrated so the null cannot
fire, the noisy conditional pipeline with truth off the prior mean (**Arm B — the honest arm**) detects
**NOTHING**: its largest fringe gap (`dlnL = 8.0`) is **below the worst pure-noise false alarm in the null
banks (9.01)**. **The conditional ceiling under real noise, honestly gated, is zero. Not small — zero.**
[MEASURED: project_progress §10.0; FORGE §9.4]
*(Robust to the floor's calibration noise: Arm B detects 0.000 under **any** floor in the ±5-nat scatter
band 4–14 nat. The headline never depended on the margin.)*

**S3.2.5** **The A→B price, and it survives every re-score because it is RELATIVE.** Arm A (truth at the
prior mean — the convention every prior Track B deliverable silently inherited from `cw_helpers.py:228`)
reproduces the census count under real noise (2.87 ≈ 3). **Arm B (truth drawn off the prior mean — the
honest case) HALVES the count (2.87 → 1.43) and QUADRUPLES the wrong-certification rate (0.067 → 0.267
per realisation; 3 confident P > 0.99 wrong certs vs 0).** Registration-from-the-prior-mean was worth
~2× in yield and was suppressing essentially all confident wrong certs.
[MEASURED: FORGE §1 — 30 realisations/arm = 10 skies × 3 noise weathers]

### S3.3 The selection function

**S3.3.1** **THE NAMES ARE A MEASURE-ZERO OUTCOME.** The census triple {J0711, J1713, J1909} is reproduced
in **0 of 40** skies. All three certify together in only 6/40; **34/40 draws have ≥1 census name fail**;
Jaccard(draw, census) = **0.384 ± 0.132**. *"The count is plausibly robust; the names are not" was half
right — and the wrong half was the reassuring one.*
[MEASURED: GEO §3]

**S3.3.2** **Per-pulsar certification frequency (40 skies):** J1909−3744 **38/40 (0.95)**; J0437−4715
**32/40 (0.80)**; J1713+0747 27/40 (0.68); J1744−1134 20/40 (0.50); J0711−6830 **16/40 (0.40)**. Union
(certifies in ≥1 draw) = **18 pulsars**. **J1909 is an 18/20 pulsar and J0711 is a 6/20 pulsar — they live
at opposite ends of the same census sentence.**
[MEASURED: GEO §3]

**S3.3.3** **The census omitted the array's best-measured pulsar by bad luck.** J0437−4715 certifies in
32/40 redraws — more often than J1713, twice as often as J0711 — and is absent from the census list only
because seed-3000 happened to be one of its 8/40 failures. **It has the array's smallest K
(K_lit = 3.07).**
[MEASURED: GEO §3]

**S3.3.4** **SKY-CONDITIONAL SEED SETS ARE MANDATORY.** *"The real seed set is J0711/J1713/J1909"* is a
**one-draw statement**; an estimator bootstrapping from that literal triple bootstraps from a set that
never recurs. Seed from **J1909-3744** and compute the seed set **per realisation**; **J0437 belongs in
any such set.**
[MEASURED: GEO §7.2 — and confirmed independently under real noise by FORGE §7.2]

**S3.3.5** **The selection function, and the trap in front of it.** Certification frequency correlates
with `1 − cos μ` to the nearest loud source **negatively and entirely through the fringe-breaking
evidence**: stratified within-pulsar `ρ(1−cos μ, dlnL) = −0.251` against `ρ(1−cos μ, K_counted) = +0.008`.
**The trials factor is BLIND to the loud sources** (at N_CW = 16, `dL` is set by whichever of the 16 has
the largest `f(1−cos μ)` — generically a faint one). The **marginal** correlation (ρ = −0.029) says
"geometry doesn't matter" and is **confounded by pulsar identity** — a trap. **Stratify.**
[MEASURED: GEO §4 — 4640 (pulsar, draw) cells, 180 certified; 15/18 pulsars have negative ρ]

**S3.3.6** `P(certify)` is **flat at 0.045–0.058** for `1 − cos μ ≲ 0.44` (μ ≲ 56°), then falls
monotonically to **0.009** in the top decile — a **5–6× suppression**. It does **not** rise as
`1 − cos μ → 0`: the near-alignment limit buys a wide fringe window but kills the pulsar-term response,
and the two cancel to a plateau. The ratio (mean `1−cos μ` when certifying / when not) is **< 1 for all
18 of 18** union pulsars, mean 0.50.
[MEASURED: GEO §4 deciles + §9.2]

**S3.3.7 The J1909 hemisphere anecdote — geometry, caught in the act.** J1909's non-certifying mean
`1 − cos μ = 1.0386` sits *just past* the hemisphere boundary (`1 − cos μ = 1 ⇔ cos μ < 0`): **the only
2 of 40 skies that break the array's most reliable pulsar are the two that put every loud source in the
opposite hemisphere from it.**
[MEASURED: GEO §9.2]

**S3.3.8** **Geometry, not weather, sets the yield.** Across 10 skies × 3 noise weathers the sky draw
dominates the variance (g03 gives 3–4 across all weathers; g05/g06/g08 give 0–1).
[MEASURED: FORGE §6 — consistent with GEO's selection function; the source of the standing ±0.2-class
sky-sampling error on every per-cell rate in IGNITE/IGNITE-2]

### S3.4 Pool vs selection — the two-stage frame

**S3.4.1** **Registration at truth is NECESSARY BUT NOT SUFFICIENT for certification.** *"The carrier set
is the pool of pulsars whose combs co-register, and certification then selects from that pool by prior
width and fringe-breaking margin."* Union-18 (GEO: certifies in ≥1 of 40 skies) and carrier-18 (P1:
registers at truth on the census sky) share **15**, Jaccard 0.714. **The equal size (18 = 18) is a
COINCIDENCE of two different quantities, not agreement.**
[MEASURED: GEO §9.1 — union-only: J0613, J1640, J2317 (the rarest certifiers); carrier-only: J1446−4701,
J1455−3330, J2145−0750 (register at truth, certify in 0 of 40)]

**S3.4.2** The carrier pool is **~18 loud-source-broken pulsars, not 3 and not 116**. At truth, true-reg
(`q > 0.5` & `k = 0`) = **18/116**, against an off-truth floor of 2–3; the census-3 register **3/3 at
truth, 0–1 off**.
[MEASURED: spec P1 result 4c]

---

## S4 — THE WALLS

### S4.1 The blind wall: three formalisms, one price

**S4.1.1 The needle EXISTS.** A unique joint-lnL registration needle sits **exactly at truth** at every
scale: coarse ±20° argmax at truth (next-best cell −102 nat); zoom ±2° argmax 0.00° from truth;
uniqueness patch ±0.05° an **interior 2-D maximum** (not a ridge or saddle). At exact truth with integers
fixed, `|grad| = 0` and the 24-loud-param Hessian is **strictly negative-definite**, eigenvalues
`[−5.96e11, −14.4]`. **The GPS-RTK / LAMBDA structure is real: the combs co-register ONLY at truth.**
[MEASURED: spec P1 + L2c — L2's "needle-absent" auto-label was an LBFGS artifact and is **WITHDRAWN**]

**S4.1.2 The needle is a CUSP, and its basin is curvature-microscopic.** Condition number **~4e10**;
registration curvature ~6e11 (the lever arm, S1.1.2); the quadratic bowl is valid only within
`δ ~ sqrt(2/6e11) ~ 1e-6` scaled. **Conditional re-solve pull-in radius < 1e-4 scaled** — from 1e-4, even
*with truth's integers*, damped Newton lands 0.046° out and 622 nat below truth.
[MEASURED: spec L2b/L2c, BANKED and float-independent]

**S4.1.3 The ladder does not span — this is the wall, and it is pure geometry.** Per-(pulsar, loud)
**sky** registration tolerance (the offset that shifts the pulsar-term phase by 1 rad): loosest
**1.85e-3**, median 3.81e-5, tightest 1.34e-6 scaled. **Pairs tolerating ≥ 0.05 (the F-stat float floor):
ZERO. Pairs tolerating ≥ 1e-2: ZERO.** The blind sky float lands at **0.05 (seed floor) to 0.21
(endpoint)** scaled. **WALL = 27× to 112× = 1.4–2.05 dex, with ZERO rungs in the gap.** The GPS-RTK
wide-lane cascade **cannot fix its first rung.**
[MEASURED: spec F2 ladder + Repair-3 float — all 348 pairs finite; trajectory-sensitive at ±0.2 dex]

**S4.1.4 Lever (i) confirmed: `tol ∝ 1/L`.** The three loosest rungs are the three **nearest** pulsars:
J0711−6830 (0.106 kpc, 1.85e-3), J1630+3734 (0.089 kpc, 1.39e-3), J0437−4715 (0.155 kpc, 1.02e-3), against
a median array pulsar at L = 1.38 kpc.
[MEASURED: spec F2 / DESIGN THEOREM lever (i)]

**S4.1.5 FORMALISM 1 — continuous optimisation: DEFEATED, and the mechanism was corrected.** The soft-EM
M-step drifts off truth and the census ceilings collapse to ~0 under all three weightings tried.
[MEASURED: spec P2 test-mstep]
[SUPERSEDED → the first mechanism given ("the soft objective's optimum is BROAD and DISPLACED; a
precision floor 100–500× coarser than the needle") was **WRONG**. The STEP-1a line scan shows HARD and
SOFT marginals **both peak sharply and exactly at truth** (argmax 0.0000 scaled; −8000 nat already at
±0.003). The real mechanism is a **CUSP ON A PLATEAU**: (a) the M-step holds `q` fixed through 40 Adam
steps and the fixed-`q` Q *is* displaced by the ~98 degenerate pulsars; (b) Adam falls off the cusp tip
onto the plateau. **Direction unchanged, mechanism fixed.**]

**S4.1.6 The plateau is a GAUGE CONSPIRACY, not absent physics.** In an all-snapped HARD surface every
one of the 116 combs is free to shift *together*, so for any displaced source there exists a globally
consistent wrong-fringe distance set that refits nearly as well. **Anchor the array and the conspiracy
breaks**: with other distances at truth, moving loud0's sky costs **−82 750 nat at 2°** (monotone to the
cusp), and `QTRUE` rises smoothly 4.9 (2°) → 9.6 (0.5°) → 20 (0.1°) → 24.45 (truth). **The information is
present but gauge-hidden in the simultaneous fringe freedom.**
[MEASURED: spec PROBE 1]

**S4.1.7 The loop closes at census-class anchoring.** HARD gap-drop (0.05° → 2°) grows monotonically with
anchors: k = 0 → 50, k = 1 → 50, **k = 3 → 132**, k = 6 → 172, k = 12 → 244, k = 24 → 324 nat; strict
monotonicity across the full gap closes at **k = 6**. And **`QTRUE` — the soft objective's driver — is
MONOTONE across the whole gap from k = 0** (6.25 at 2° → 24.45 at truth). The anchors that emerge *are*
the census set (idx 88/62/19 = J1909/J1713/J0711, q₀ = 1.00/0.99/0.95).
[MEASURED: spec PROBE 2 — the pre-registered FAVORABLE branch, which is why the soft path was then built
and tested]

**S4.1.8 FORMALISM 2 — discrete integer least-squares (LAMBDA/RTK): DEFEATED at BOTH gates.**
(1) **Integer neighbourhood**: the truth-blind float lands O(1) scaled off truth, so the bounded search
around the float's MAP integers (census −113/−143/−19; truth n = 0 for all 116) **never contains truth's
integers** — candidate sets (>1 % mass) span n ∈ [−182, 226] and do not include 0. (2) **Source pull-in**:
even GIVEN truth's integers, the basin is < 1e-4 scaled (S4.1.2) — 4+ orders inside any achievable float.
**GPS-RTK works because its code/pseudorange float lands within ONE carrier cycle; the PTA F-stat float
lands ~27–1300 pulsar-term cycles from the loosest baseline, and no fixing step spans that.**
[MEASURED: spec LAMBDA probe L0/L1/L2 + F2 EARNED VERDICT — NO-GO, re-derived independently of the
contaminated float]

**S4.1.9 The contamination asterisk, carried honestly.** The first NO-GO's float-dependent chain (L0
diverged float, L1 search-space spec, L2 blind search) carries a contamination asterisk: the M-step is
**numerically CHAOTIC** (multi-step FP-order amplification through a ~4e10-conditioned Adam — *not* the
NaN-poisoning first diagnosed, which was **tested and REFUTED**: the fix is a byte-identical no-op), and
the F-stat seeder was mis-seeding loud2 (fixed: sky-only NMS at 25°, 3/3 loud, rank-13/17.70° →
rank-11/11.88°). **UNAFFECTED and BANKED: L2b/L2c (fixed-integer pull-in) and F2 (pure geometry).** The
repaired float reaches the same wall.
[MEASURED: spec REPAIRS LOG 1–3 + CONTAMINATION ASTERISK]

**S4.1.10 FORMALISM 3 — the honest POSTERIOR: DEFEATED, and this is what makes the NO-GO
information-theoretic.** The fringe-marginalised posterior **does not concentrate at the needle — it
SMEARS across the plateau.** `ln Z_plateau − ln Z_needle = +14.19 nat` →
**f = Z_needle/(Z_needle + Z_plateau) = 6.9e-7.** The needle is **28 nat higher at its PEAK** but its
6-D sky volume (~e^−73) loses to the ±2° plateau by 14 nat. Census-3 `P(true)` at the posterior peak
(plateau MAP) = **0/0/0**.
[MEASURED: spec DELIVERABLE R — R1 `ln Z_needle = 405617.64` (quadrature) / 405617.84 (Laplace), agreeing
to 0.2 nat, 6/6 positive curvature, needle σ_sky 2.4–9.7e-6 scaled; R2 `ln Z_plateau = 405631.83 ± 0.053`
(2 seeds); all gates PASS including `exp(−m_p) == A2 P_true` to 3.7e-8]

**S4.1.11 The break-even, and the one thing that would buy it.** `f = 6.9e-7` ⇒ **break-even sky prior
θ* = 0.188° per loud source**: the needle wins only if the sky is pre-localised below ~0.19° — **32×
tighter than the F-stat blind floor (~6°)**. **An EM host supplies exactly this.** This is design-theorem
lever (iii), and it is why the programme pivoted to the targeted scenario.
[MEASURED: spec R3]

**S4.1.12 THE ONE PRICE.** All three formalisms fail on the **same** quantity — the lever arm
`2πL/dL ≈ 1.6e4` — expressed three ways: as a **cusp** too sharp for a continuous optimiser (S4.1.5), as
an **integer neighbourhood** too far for a bounded search (S4.1.8), and as a **volume ratio** too small
for the evidence (S4.1.10). *Certification == localisation == fixing the integers*, and the exchange rate
is the lever arm. **A LOCKED registration localises the source at mas class** (cusp width ~1e-5° ≈ tens of
mas) — which is precisely why the target needle is mas-narrow and the gap is the problem.
[MEASURED: spec P1 COROLLARY + TERMINAL VERDICT]

**S4.1.13 The wall expressed as a timescale (the pencil).** The Earth-term blind float reaches the F2
loosest rung (1.85e-3) only at **T ~ 11 kyr** (σ_sky ~ T^−1/2), or at fixed T at per-source
**SNR ~ 289** (`h ~ 10^−12.3`, ~27× the current loudest); the L2c pull-in at **T ~ 3.75 Myr**. **The wall
does not close on any observational horizon.**
[MEASURED: spec R5 / `PENCIL_t_crossing.md`, banked]

### S4.2 The targeted-circular wall: the chirp-mass wall

**S4.2.1 The premise that failed — `log10_mc` is a THIRD registration axis.** The B1 brief presumed one
remaining axis (frequency) once the sky is supplied. There are **two**. The model chirps the pulsar term
at the retarded time, so the pulsar-term phase depends on `log10_mc`; **F2's ladder used the chirp-free
phase, whose `∂/∂log10_mc` is IDENTICALLY ZERO — F2 was structurally blind to an mc lane and could not
have found one.**
[MEASURED: spec B1 STEP 1 — reproduction gate: the pilot's exact-phase autodiff reproduces F2's chirp-free
freq ladder to 3 digits (2.518e-2 / 1.035e-4 vs 2.52e-2 / 1.04e-4), so the disagreement is the chirp, not
the code]
[SUPERSEDED → **F2's SKY ladder and the Track B terminal verdict STAND** (both rest on the sky axis and on
L2b/L2c). **F2's FREQUENCY ladder is 0.72× optimistic** (median).]

**S4.2.2 The 8 source parameters split cleanly into three classes.** Per-axis **certification** tolerance
(E-step, one axis at a time, the delta at which census-3 `P_true` collapses): **sky 1e-5; log10_fgw 3e-5;
log10_mc 1e-3; extrinsics (cos_inc, log10_h, phase0, psi) NO collapse to 1e-2.** So: **4 harmless
extrinsics + 2 sky (supplied by the counterpart) + 2 registration axes (f, mc) a targeted pipeline must
still find.** This also *measures* the assumption R1/R2 asserted (extrinsic Laplace factors cancel in f).
[MEASURED: spec B1 pilot M2 — directly comparable to B0.2's ~1e-4, which perturbed all 8 at once and could
not attribute]

**S4.2.3 The Earth term cannot measure `fdot`.** Information gain over the prior: **162–389× for
log10_fgw** but **1.00–1.73× for log10_mc** (σ = 0.501/0.864/0.866 dex — *for two of the three loud
sources the posterior IS the prior*). `Δφ_E = π fdot T² ≈ 0.05 rad` for the loudest. Per-axis targeted
wall: freq **203×**, **mc 1730×**.
[MEASURED: spec B1 pilot M3 — HVP-assembled Earth-only Fisher, inverted WITH the generative priors, NEVER
pinv (pinv reports σ = 0 for an unconstrained axis)]
*(Independently reproduced by SIREN's N_seed = 0 gate to 1e-5–1e-3, and by SIREN §3.1: in **8 of 9** cells
the Earth-term posterior on `log10 Mc` **is** the prior, gain 1.000–1.051.)*

**S4.2.4 95–100 % of the pulsar-term phase wander is `mc`.** (J0711 95.0 %, J1713 100.0 %, J1909 99.9 %,
J0437 99.3 %.) Registration radius `R_a` (box shrink to reach 1 rad): max **2.14e-2** (J0711), median
2.04e-4. **IGNITION: 0 pulsars.** The loosest rung needs a **47×** box shrink; the binding census pulsar
(J1713) needs **3305×**. **Unlike sky, the obstruction is NOT rung spacing** — the ladder spans internally
(42 rungs, max gap 0.387 dex < F2's 0.7 criterion). **It is the mc BOX.**
[MEASURED: spec B1 STEP 1A]

**S4.2.5 The structure is BISTABILITY, and certification is SELF-REFERENTIAL.** Loop gain (conditional
(f, mc) Fisher at a fixed fringe subset — a masked likelihood, so **no optimiser and no L2c pull-in
problem**): **0.04** at the float, 0.12 (top-1), 0.70 (top-24), **1.34 (top-48)**, **4.50 (census-3
only)**. The registered state exists and strongly attracts (the P1 needle, seen from the (f, mc) side) —
**it is simply unreachable from the float.** ***~30 registered fringes are needed to measure the chirp
mass that lets any fringe be registered.***
[MEASURED: spec B1 STEP 1B — Asimov injected with the SAME mask ⇒ `H = −Fisher` exactly]

**S4.2.6 A pre-registered structural hypothesis was REFUTED and is recorded, not quietly dropped.**
The guess — *"loose rungs are loose BECAUSE they are blind to mc, so nothing can bootstrap"* — is **false**:
fixing J0711 alone cuts σ_mc **7.8×**; the census-3 alone cut it **1600×**. **Loose rungs carry mc
information through their AMPLITUDE/SNR, not only through `g_mc`.**
[MEASURED: spec B1 STEP 1B, KILLED HYPOTHESIS]

**S4.2.7 No conditioning tier ignites the cascade.** Tier A (sky only) best `R = 0.0214`; **Tier B (+ an
EM period, σ_P/P = 1e-3) `R = 0.0220` — an EM period buys NOTHING** (a 7× tighter f moves the best rung by
3 %, consistent with the 95–100 % mc share); Tier C (+ host `D_L`) **`R = 0.271`** — missing ignition by
**3.7×**. Tier C's mc box is set by the **ARRAY's own σ(log10_h)** via
`log10 h = (5/3)log10 Mc + (2/3)log10 f − log10 D_L`, i.e. **not by anything the counterpart supplies
beyond the redshift.**
[MEASURED: spec B1 STEP 1C — a pre-run prediction of "misses by 1.4×" was **WRONG**: it used the *median*
mc shrink; `R` is a quadrature dominated by the **worst** source (loud0, worst σ(log10_h)). Recorded.]

**S4.2.8 The targeted referendum: ~97 % of the posterior's evidence sits on the wrong-fringe plateau EVEN
WITH THE SKY EXACT.** Three tiers, frozen 4-seed protocol, R's count-once star-topology marginal verbatim:

| tier | conditioning | ln Z_needle | ln Z_box | d (nat) | **f** | ±2σ | gate |
|---|---|---|---|---|---|---|---|
| A | sky only | 405629.637 | 405632.017 | −2.380 | **0.0847** | ±0.131 | FAILED |
| B | + EM period | 405629.634 | 405632.619 | −2.986 | **0.0481** | ±0.0227 | FAILED |
| C | + host D_L | 405629.634 | 405632.734 | −3.101 | **0.0431** | ±0.0369 | FAILED |

**`ln Z_needle` is tier-independent to 0.003 nat** across three boxes and two independent bracket
algorithms. **The tier gradient is FLAT and mildly INVERTED** — f *falls* as conditioning tightens.
**Counterpart information does not help.**
[MEASURED: spec B1 STEP 2 / `b1_step2_table.npz`; all three FAIL the 0.95 bar by 4.4×–13× on the ±2σ band]
[DISPUTED: Tier C frozen at 4 seeds reads **f = 0.0323 ± 0.0134** (FAIL by 16.1×); the table's auto-ingest
used the completed 5-seed npz and reads **f = 0.0431 ± 0.0185**. **Both FAIL identically. The discrepancy
is recorded, not resolved by preference.**]

**S4.2.9 THE SATURATION MECHANISM — and it was caught only because a "logically foreclosed" arm was run
anyway.** Tier A's box is **e^10.1 ≈ 24 000×** larger in volume than Tier C's, yet `ln Z_box(A)` and
`ln Z_box(C)` agree to within the 0.86-nat sampler scatter (measured difference −1.02 ± 0.95 nat, i.e.
**consistent with EQUALITY**). `Z_box` is an integral, so `Z_box(A) ≥ Z_box(C)` is *required*; measuring
them equal means **the plateau's evidence has SATURATED — it is confined well inside Tier C's box.**
**Enlarging a prior box adds volume carrying negligible likelihood.**
[MEASURED: spec SATURATION NOTE + CONVENTION "logically-redundant measurements retain audit value" —
*saturation is invisible from one box BY CONSTRUCTION*]

**S4.2.10 The >20× bound (the corrected price), measured as a CURVE.** `Z_box` at
`λ_mc ∈ {1, 0.3, 0.12, 0.05}` (needle excision 0.0 % everywhere): `ln Z_box` = 405633.035 / 405631.754 /
405630.535 / 405628.910 → **f = 0.032 / 0.107 / 0.289 / 0.673.** Readings: **(a) the plateau's own
chirp-mass extent is ~0.02 dex** (a newly measured quantity — the width of the wrong-fringe plateau in Mc);
**(b) λ\* < 0.05 — a BOUND, not a value**; **(c) the deficit is > 20×**, i.e. **σ(log10_mc) must fall below
~0.003 dex.** *Even a 20×-shrunken mc box only reaches f = 0.673 — it still fails.*
[MEASURED: spec STEP 2C / `b1_breakeven_curve.npz`]
[SUPERSEDED → the deficit **8.29×** (and with it the two-term price "σ_h ×11.3 **AND** σ(log10 D_L) ≤ 1.0 %",
and WEAVE's `κ ≥ 8.3` threshold) was computed by extrapolating `Z_box ∝ λ_mc³` — **exactly across the range
over which `Z_box` is now shown to be insensitive.** It was **SUSPENDED** on discovering saturation and then
**corrected to >20×** by the response curve. **The suspension was vindicated: the real price is WORSE than
the proportionality implied.** *(A log-linear extrapolation of the last two points gives λ\* ~ 0.015,
deficit ~66×. **NOT reported as a result** — that is the very error the curve exists to correct.)*]

**S4.2.11 No combination of counterpart information delivers it.** σ(log10_mc) delivered by Tier C =
0.0364/0.0217/0.0205 dex. Setting **σ_h → 0** leaves a **FLOOR of 0.00301 dex**, set entirely by the
assumed 1 % host distance. At the corrected >20× deficit that floor is **at or above the requirement for
every loud source**: **strain alone cannot close it, and neither can any counterpart-supplied prior box.**
[MEASURED: spec "THE PRICE, REQUOTED OFF THE CURVE"]

**S4.2.12 THE HEADLINE MECHANISM, and it is the hinge of the forward story.** ***Conditioning the (f, mc)
PRIOR BOXES barely moves the evidence, because the plateau does not fill them. What moves the evidence is
LIKELIHOOD STRUCTURE — which the eccentric harmonic comb supplies and no prior box can.***
[MEASURED: spec B1 STEP 2 HEADLINE]

**S4.2.13 The last door, checked and closed.** A soft posterior-weighted cascade at `R ≈ 0.27` (spreading
over ~4 fringes rather than locking one) might still leak Mc. It does not: per-iteration σ_mc shrink
`[0.057, 0.306, 3.335, 1.428, 1.153]` — cumulative 3.77× but **NON-MONOTONE**, which is exactly why the
pre-registration demanded monotonicity — while the number of pulsars whose false-fringe mass **GREW**
climbs **54 → 70** (the soft analogue of the GPS wrong-fix). **Mechanism: σ_mc is ALREADY 1e-4–1e-3 dex at
iteration 0 — far narrower than the 0.0026–0.0044 needed. The local mode was never the problem.** Median
`q_max` is flat (0.067 → 0.070; ~16 effective fringes/pulsar, unchanged). ***The missing information is not
local curvature — it is WHICH of ~16 fringes, and that choice never sharpens. Local width ≠ global
concentration.***
[MEASURED: spec STEP 1D — pre-registered, one probe, one verdict: FAIL]

**S4.2.14 The needle is a thin SHELL, not a point.** The Hessian of `lnL_marg` at truth over (f, mc) has
**5/6 positive eigenvalues and one NEGATIVE** (−1.32e9), Richardson-stable (2.4e-2). It is **not a saddle**
(`lnL_marg` falls 159 nat at 1× base, 919 at the box edge along that eigenvector) but a **MICRO-DIP**:
`lnL_marg = lnL_ref + Σ_p m_p` with `m_p ≥ 0` growing as a pulsar de-registers, so a sub-fringe offset buys
more fringe entropy than it costs; the local max sits ~1e-5 scaled away, 0.12 nat higher. **CONSEQUENCE:
`Z_needle` is quadrature-only, NEVER Laplace, and negative eigenvalues must NEVER be clipped** — clipping
inflates `Z_needle`, i.e. **biases toward the answer that lets the campaign proceed.**
[MEASURED: spec "The needle is a thin SHELL"]

**S4.2.15 R is unaffected, and firmer.** R fixed `mc` at truth and argued its Laplace factor cancels; **for
mc that argument is WRONG** (needle σ ~1e-3 scaled vs plateau ~1.7). **The error is FAVORABLE**: including
mc as an active dim *shrinks* `Z_needle`, so **`f = 6.9e-7` STANDS and is FIRMER than reported.**
[MEASURED: spec B1 STEP 1 doc item (b) — a retroactive note that strengthens, and is recorded anyway]

**S4.2.16 PHYSICAL IDENTIFICATION — the frontier statement in its circular form.** ***The pulsar term is a
kyr-baseline TIMESTAMP of the source's phase. It cannot be read without the CLOCK RATE. The clock rate is
`fdot`, i.e. `Mc`. A 22-yr Earth term cannot measure `fdot`.*** The fringe index is set by the accumulated
chirp over `τ_p`, so the array must know `Mc` to place *any* fringe. **This is the E-track's
eccentric-harmonic mechanism in CIRCULAR form: the E-track and the targeted pipeline have MERGED — design
theorem lever (ii) is not an alternative to lever (iii), it is the ingredient lever (iii) is missing.**
[MEASURED: spec B1 STEP 1 doc item (c)]

### S4.3 The one-wall coincidence

**S4.3.1** Four independently derived thresholds land within a factor of ~1.5 of each other, **and nothing
in the repo explains why**:

| wall | quantity | value |
|---|---|---|
| blind (S4.1.3) | float ceiling / loosest sky rung | **27×** (1.4 dex; up to 112× at the chaotic endpoint) |
| targeted (S4.2.10) | required σ(log10_mc) improvement | **> 20×** |
| ATLAS relative bound | σ(log10_mc) improvement the comb must supply | **> 20×** |
| WEAVE self-clock | `κ ≥ 20` (`Δφ_E ≳ 1 rad` on the conditional Fisher) | **20** |

**STATUS: UNEXPLAINED.** These are *different objects* — a sky-localisation ratio, an evidence-derived
box-shrink bound, a σ-improvement factor, and a Fisher-ratio threshold. ATLAS explicitly warns that the
last three are **three distinct "20"s that earlier work conflated** and must be kept apart. Whether the
coincidence is (a) a genuine common origin in the lever arm `2πL/dL`, (b) an artifact of this array's
`L`-distribution and `T`, or (c) numerology, **is not settled by any banked measurement.**
[MEASURED: spec F2 / STEP 2C; ATLAS "THE QUALIFYING STATEMENT"]
[PENDING: QUILL — a first-principles treatment should either DERIVE the coincidence from the lever arm and
the array's `L`-distribution, or DEMOLISH it as three unrelated numbers that happen to be ~20. Until then,
**never quote the four as if they were one.**]

---

## S5 — THE CRITERION

### S5.1 The three layers (criterion-v2.1, CANONICAL)

**S5.1.1** The operative criterion:

```
DETECTION      dlnL_a > max( ln K_counted,a , floor(h, T, tol) )      [D1 family, D2 estimator]
CERTIFICATION  q_max,a > 0.9   (strict: > 0.99)   applied ONLY within detections
PURITY         NONE. Both candidate layers were tested and BOTH ARE REJECTED.
```

Layer 1 (`ln K`) is the **trials factor** — relative, per-pulsar. Layer 2 is the **absolute floor** —
the evidence a detection must carry regardless of how few fringes it had to beat. Layer 3 (`q_max`) is a
**quality bar on something already established as a detection**, not itself a detection statistic.
[MEASURED: spec CRITERION-v2.1 — gates `CW_transition/criterion_v2_gates.py` 8/8 PASS + census triple
bit-identical]

**S5.1.2 BINDING INVARIANT: certification is defined on `q_max`, NOT on `P_true`.** The cells where
`q_max > 0.9` but `P_true ≤ 0.9` **ARE the wrong-certifications** (2 in Arm A, 8 in Arm B). Scoring on
`P_true` would silently **define them out of existence**. Within detections, `q_max == P_true`.
[MEASURED: spec criterion-v1 — asserted in code (`trackB_criterion.py`)]

**S5.1.3 BINDING INVARIANT: GEO has no wrong-cert field BY CONSTRUCTION** (its criterion is defined on the
true fringe). **Not synthesised.**
[MEASURED: project_progress §10.0 footnote 1]

### S5.2 Why a two-layer gate is not enough — the trail

**S5.2.1** **The Bayesian bar (`P_true > 0.9`) is not a detection statistic at all.** It is a statement
about *where posterior mass sits among candidate fringes*, not about *whether there is anything there to
find*. **`nullN` — pure noise, NO CW in the data — certified 0.8 pulsars/realisation at the Bayesian bar**,
and 2 of them at P > 0.99.
[MEASURED: FORGE §3 / B1.2 — the finding that forced the criterion]

**S5.2.2** **The scrambled null FIRED, and the first hypothesis was REFUTED.** Bayesian certs: 2.2/real
(3 loud scrambled), **2.8/real (all 16 scrambled)**, 0.8/real (no CW at all). *Scrambling all 16 does not
reduce the count* — so the floor is **not** the faint sources staying coherent; it is **intrinsic to the
Bayesian criterion**: a prior-pinned floor (~0.8/real, pure noise) plus a **noise-lock excess** (→2.8/real)
when a wrong source model meets real CW data.
[MEASURED: FORGE §3]

**S5.2.3** **The two-layer gate (`dlnL > ln K` ⊕ `q > 0.9`) still lets the null fire — from small K.** A
tightly-EM-prior'd pulsar has so low a trials bar (**J0437: `ln K = 1.39`, the array minimum**) that a pure
noise fluctuation clears it. Two-layer null: **0.2–0.4/real**, and Arm-B's two-layer detections
(**0.13/real**) sit *at or below* it — **by count, the noisy pipeline does not detect above its own null.**
[MEASURED: FORGE §9.3/§9.4 — FORGE predicted the fix ("an absolute floor, ~8 nat") and correctly declined
to impose it, flagging it for the spec rather than silently changing the criterion mid-report]

**S5.2.4** **CONVENTION (adopted): confidence without a detection statistic is prior-pinning in disguise.**
Every confidence bar must sit downstream of a detection statistic **that can return zero**. A criterion
that cannot fire on a null is not a criterion. **Corollary: robustness to source error and vulnerability to
noise are THE SAME PROPERTY viewed from two sides.** Tiny `K` cuts both ways, always.
[MEASURED: spec CONVENTION — the J0437 double edge, measured in both directions: GEO detect-freq 0.65 at
zero noise → 0.10 under real noise; and the dominant residual-null false-alarm source]

### S5.3 Floors as Gumbel quantiles (D2)

**S5.3.1** **The floor is NOT a constant — it is loudness-relative.** `floor_fN ∝ h^1.66`
(per-(T,tier) fits span 1.5–1.7); `floor_fALL ∝ h^1.88` (span 1.7–2.0), measured in **every** baseline and
tier. **Mechanism, and it runs on data containing NO CW at all:** the E-step evaluates a model whose
pulsar-term amplitude ∝ h, and the per-fringe likelihood carries a **matched-filter cross term linear in the
MODEL amplitude** — so the null's fluctuations grow with the loudness of the **hypothesis**, not of the
data. With a *scrambled* source meeting loud real data the noise-lock grows ∝ h², which is why `fALL` scales
more steeply. Concretely the h-law lives in the **Gumbel scale**: `β` = 2.1–2.4 nat at h = −13.0, 4.2–7.0 at
−12.75, 13–24 at −12.5.
[MEASURED: IGNITE §4.1; spec D2.1]

**S5.3.2 CONSEQUENCE: the certified count is NON-MONOTONE in h.** T = 20 vlbi: **0.72/real at h = −13.25 →
0.38/real at h = −12.5** — a 10× louder source *lowers* the honest count. **Louder alone does not buy
certification; the bar rises almost as fast as the signal.**
[MEASURED: IGNITE §4.1]

**S5.3.3** **`9.01 nat` is RETIRED.** It was never a property of the pipeline; it was the **census-loudness
value of a function**.
[SUPERSEDED → criterion-v1's `DLNL_FLOOR = 9.01`, fitted as the smallest floor zeroing all 27 banked nulls,
with the binding cell a **`nullN` J1713+0747 fluctuation at `dlnL = 9.009` on data containing NO CW AT ALL**.
**Every criterion-v1 number stands exactly as measured** and is gated bit-identically (G1–G3); what changes
is the **scope**: they are census-loudness numbers and may not be extrapolated to other loudnesses,
baselines, or tolerances.]

**S5.3.4 The estimator had to change, and this is the part the data forced.** criterion-v1's floor is the
**maximum of N nulls**. The null offender statistic is itself a max over pulsars, so it lies in the **Gumbel
domain by construction** — and there **`sd(max_N) = 1.283·β`, INDEPENDENT of N**, while
`E[max_N] = μ + β·ln N` **creeps up without bound**. Measured: `sd(max_N)` = **8.91 / 8.68 / 8.79 / 8.74**
nat at N = 10/30/100/1000. **FLAT.** ***Banking more nulls does not stabilise a max-of-N floor — it inflates
it.*** Worse, max-of-N has **no fixed false-alarm rate**: it is implicitly the `1 − 1/(N+1)` quantile, so its
stringency was **an accident of how many nulls happened to be banked.**
[MEASURED: spec D2.2 / IGNITE G7 — IGNITE's own recommendation ("more nulls is the cheapest credibility
purchase") is **half right**, and this is the correction: **the estimator had to change first.**]

**S5.3.5 ADOPTED ESTIMATOR.** `floor(h, T, tol) = μ̂ + β̂·z(α)`, `z(α) = −ln(−ln(1−α))`, **α = 0.05**, via a
Gumbel (block-maximum) tail fit. Its scatter **does** shrink: `sd(floor̂) = 2.80·β/√N`. **Sizing (both rules
binding): N ≥ 100 nulls/cell** (scale-free: sd < 10 % of the floor at *any* loudness) and **N ≥ 150 at onset
cells** (absolute: sd < 1 nat). A 1-nat absolute target at the loud h = −12.5 cells would need N ≈ 2000–5000
and is **explicitly not adopted**. Cost: ~2 GPU-hours for 150 × 24 cells. ***The sizing was always cheap;
what was missing was the check that the estimator converges.***
[MEASURED: spec D2.2 — criterion-v1's max-of-27 was implicitly α ≈ 0.036, so α = 0.05 is an explicit,
N-independent version of a bar previously set by accident]

**S5.3.6 The tolerance hole is CLOSED, and it INVERTS.** The pure-noise floor is **flat-to-mildly-rising and
small** in registration tolerance: `fN` = **0.00 → 0.00 → 2.06 → 4.37 nat** at tol = 0/1/2/5. **It is the
TRUE-POSITIVE channel that dies of mis-registration** (true certs 0.14 → **0.00**/real by tol = 5), and **no
per-tol refit floor kills a single surviving true positive.** criterion-v1's *"the 9.01 floor kills the
correct `wrongpos` J0437 certification (dlnL = 4.41)"* pathology was an **artifact of applying a tol = 0
floor to a tol = 5 realisation** — calibrated at its own tolerance, the floor (4.37) sits *below* that
survivor.
[MEASURED: IGNITE §3 g3; spec D2.3 — the `fALL` spread {8.48, 14.03, 8.09, 4.37} is **sampling noise, not
tol dependence**: four independent 30-null redraws of one statistic]

**S5.3.7 Every margin in the programme is inside the calibration noise.** The floor's own scatter is **±5
nat per 30-null refit**. criterion-v1's **0.29-nat** margin and IGNITE's **0.01–2.0-nat** margins are
**annotated as WITHIN CALIBRATION NOISE and carry no evidential weight.** *A 0.29-nat gap against a ±5-nat
ruler is not a margin — it is a rounding error.*
[MEASURED: IGNITE §3 g2; spec D2.4]
[SUPERSEDED → *"the margin is thin (0.29 nat) and more nulls would tighten it"* — **the second half is
false** (S5.3.4). **No criterion-v1 verdict moves**: Arm B's largest `dlnL` is 8.0 against a floor whose
scatter band is 4–14 nat, so *"the honestly-gated noisy pipeline detects nothing at census loudness"* never
depended on the margin. **What dies is the margin's precision, not the result's sign.**]

### S5.4 The purity rejections, with anatomies

**S5.4.1 criterion-v1's purity property is DEAD above onset.** It claimed *"the gate does not merely thin
the count — it perfectly purifies what is left"* (0 wrong-certs, both arms, down from 2 and 8). **That was a
CENSUS-LOUDNESS ARTIFACT.** At the (−12.75, 30 yr, lit) onset cell, **23 of 50 realisations carry a wrong
certification** — the same noise-lock that raises the floor gives *wrong* fringes floor-beating gaps.
***Fringe correctness — the one discriminator that survived real noise at census loudness — degrades exactly
where the count turns on.*** **Any above-onset count without a purity number beside it is meaningless.**
[MEASURED: IGNITE §4.3 — the wrong certs concentrate in the wide-prior, data-driven pulsars (J1909, J1045,
J1603, J1713); **the anchor J0437 supplies 20/50 correct certs at the onset cell with 0 wrong**]

**S5.4.2 Above onset the workhorse FLIPS from J1909 to J0437.** Tiny `K` wins once `max(ln K, floor)` is
**floor-dominated** — the trials bar stops binding, and the anchor's genuine break rides a low bar.
[MEASURED: IGNITE §4.4]

**S5.4.3 D3 (the per-pulsar purity layer): PRE-REGISTERED, TESTED, REJECTED BY ITS OWN PRE-REGISTRATION.**
Scorecard: **(a) PASS** — 23/23 (100 %) of wrong certs killed, with the `R_det` control landing at *exactly*
its **pre-recorded 87 % ceiling** (which is what shows the co-registration idea is doing the work);
**(b) FAIL** — only **3 %** (lit, 2/77) and **67 %** (vlbi, 39/58) of TRUE certs survive, against a ≥90 %
bar; **(c)** 42/42 = **100 %** rejection of wrong-counterpart detections at the realisation level;
**(d) FAIL** — one scrambled-loop certification survives (J1909−3744, Δk = −4, `R_all` = 4.65 vs the 9.21
bar). **No threshold was tuned. No partial adoption was taken.**
[MEASURED: IGNITE2 §1; spec D3 pre-registration + criterion-v2.1 D3 verdict]

**S5.4.4 D3's ANATOMY: the REFERENCE is poisoned.** Above onset the array-wide fringe field is itself
poisoned — the leave-one-out reference `u_R` is dragged by confident wrong fringes **everywhere**, so a TRUE
cert (`u_a = 0`) **fails concordance with its own poisoned reference.** The veto measures *"this
realisation's fringe field is discordant"* (true of **every** realisation above onset) rather than *"this
pulsar disagrees with the others"*. **Structural, not a tuning artifact — and it fails hardest exactly where
it was needed most** (3 % survival at the wide-prior lit cell vs 67 % at the tight-prior vlbi cell).
[MEASURED: IGNITE2 §1.2]

**S5.4.5 D4 (the realisation-level salvage): DESIGNED, PRE-REGISTERED, TESTED, REJECTED.** D3 left exactly
one live lead — its (c) = 42/42. D4 promotes it to a gate: flag a realisation whose **detected set**
co-registers at a source *other* than the assumed counterpart, and veto every certification in it. The
statistic `S_det` was **chosen on the banked distributions before any condition was scored**, because it is
the only aggregate whose true-signal distribution **concentrates at the null value** (median **0.0** at both
onset cells). **VERDICT: (i) FAILS IN ALL EIGHT PRE-REGISTERED COMBINATIONS** (2 dk-conventions × 2
thresholds × 2 cells). Best catch at a ≤10 % false-flag rate: **90.3 %** against the **≥95 %** bar; the one
setting catching 97.5 % flags **44 %** of true-signal realisations against the ≤10 % bar.
[MEASURED: D4 §3; value gate: reproduces IGNITE-2's banked statistic for all **1089** candidates
(max |Δlog10 R| = 1.2e-10) — *"bank the statistic, not the verdict" paying for itself a second time:
D3's verdict could not have been re-examined, but D3's statistic could*]

**S5.4.6 D4's ANATOMY, and it is one statement: `S_det` is a `|Δk|` detector, and `|Δk|` is NOT the
difference between a right and a wrong counterpart.** **The misses:** every missed wrong-counterpart
realisation is a noise-lock **within ±1 fringe of truth** (median max|Δk| among detections = **1**, vs
**137** (lit) / **13** (vlbi) for the caught). **The limit case is decisive — one missed realisation has
Δk = 0: a wrong counterpart whose surviving detection sits on the TRUE fringe.** *The fringes co-register
because they are RIGHT; the SOURCE is wrong.* **The false flags:** at the lit onset cell **13 of 36 (36 %)**
detecting TRUE-signal realisations have an **impure detected set**, and the gate's false-flag rate there is
**36.1 %** — **these are the same number.** The gate faithfully measures the cell's own impurity and cannot
beat it. (VLBI cell: 12 % impure → **0 %** false flags.)
[MEASURED: D4 §4]

**S5.4.7 THE SCISSORS.** ***D3 failed because the REFERENCE was poisoned; D4 fails because the POPULATION IT
MUST PROTECT is itself poisoned. Same disease, one level up.*** Above onset a true-signal realisation and a
wrong-counterpart realisation **contain the same kind of object** — a confident noise-locked fringe — and a
geometry test cannot tell which of them is the counterpart. **Therefore: NO co-registration statistic can
close the D1 hole in general.**
[MEASURED: D4 §4 / project_progress §10.11]

**S5.4.8 The one genuine positive, reported in full.** **All three** scrambled-loop keepers ARE flagged by
the realisation-level form — **including the small-|Δk| J0437−4715 (Δk = −4) case that defeated the
per-pulsar statistic** (`R_all` = 4.65 → MISSED; `S_det` = 55.9 → FLAG); B1937+21 (Δk = +21) → 1728;
J0711−6830 (Δk = +231) → 3.2e5. ***The hole is closable on every instance this campaign holds — and no gate
that closes it survives condition (ii).*** That is the hole's status, exactly: **not** *"no statistic sees
these events"* **but** *"the statistic that sees them cannot distinguish them from the impurity the true
population already carries at the only cells where the count turns on."*
[MEASURED: D4 §5]

**S5.4.9 The purity number that travels beside every above-onset count, permanently: 14/50 realisations
carry a wrong certification at the lit onset cell under the fresh floors** (23/50 under the retired
max-of-10 floors).
[MEASURED: IGNITE2 §2; D4 §6.4]

### S5.5 The oracle-indexing caveat

**S5.5.1 CONVENTION (adopted): a statistic evaluated against truth is an ORACLE until its implementable form
is scored.** The fringe grid is indexed about the **EM-prior mean**, so D3's `dk = mapk − n_true_grid` is
referenced to the **TRUE** fringe — **which a real analysis does not know.** D4 therefore scored **both** the
ORACLE form and the IMPLEMENTABLE form (`dk = mapk`, prior-referenced, with the `(1−q_max)` factor
**dropped — forced by the change of reference, not tuned**).
[MEASURED: D4 §1.1; spec CONVENTION]

**S5.5.2 The implementable form is 2–4× WEAKER** (catch 25–52 % vs 43–97.5 %), **because `σ_EM/dL` is
O(150–800) fringes in the lit tier: the EM prior is wide enough to absorb almost any source displacement, so
a wrong counterpart does not look displaced relative to a prior that was never going to localise the fringe
anyway.**
[MEASURED: D4 §1.1]

**S5.5.3 THIS CAVEAT TRAVELS BACKWARD ONTO D3.** ***Every D3 number — (a) = 100 % and (c) = 42/42 included —
was computed in the ORACLE convention.*** **No co-registration number in this repo may be quoted as an
achievable power without its implementable-form value beside it.**
[MEASURED: D4 §1.1 warning block]

**S5.5.4 The constructive corollary: the gap CLOSES with `σ_d`.** D4-OBS is **1.6× stronger** in the VLBI
tier (51.6 %) than in the lit tier (32.9 %) — **the same lever RING identified and the same lever the onset
map rewards.**
[MEASURED: D4 §6.3 — see S7.2]

### S5.6 The null-family fork (D1)

**S5.6.1 DECISION (executed, with rationale, not reopened): the operative calibration family is `fN` —
counterpart-matched nulls** (pure noise / no CW in the data, recovery at the TRUE source position). **The
all-null family `fALL`** (adding wrong-counterpart scrambles) **is retained PERMANENTLY as the blind-robust
column and travels beside every `fN` number, in the docs and in the npz.**
[MEASURED: spec D1]

**S5.6.2 The rationale.** A **targeted** analysis faces exactly the counterpart-matched null: a real
counterpart exists **by construction**, and the false alarm it can actually suffer is **noise mimicking
fringe-breaking under the CORRECT model.** A sky-scramble null asks whether the pipeline can be fooled by a
source that is *not there* — a **blind-search** question the targeted analysis does not ask. Calibrating a
targeted criterion against a blind-search null **imports a bar the scenario never has to clear.**
[MEASURED: spec D1]

**S5.6.3 THE PRICE OF D1, RECORDED PERMANENTLY AND NEVER SUPPRESSED: under `fALL` there is NO ONSET ANYWHERE
in the modelled grid.** Best cell = **0.24 certifications/realisation, of which 0.22 correct** — against a
>1 bar, at **every one of the 24 cells**. The scrambled-source null's noise-lock grows ∝ h², so the
wrong-counterpart-robust floor rises **faster than the signal** and **closes the very window it was meant to
guard.** ***This is not a footnote. The targeted programme's onset exists BECAUSE OF the null family it is
calibrated against. Any onset number quoted outside this repo carries its `fALL` column or it is not
quoted.***
[MEASURED: IGNITE §4.2; spec D1]

**S5.6.4 What D1 gives up: the WRONG-COUNTERPART HOLE, and it is OPEN.** `fN` presumes the counterpart is
right, so it has **no defence against a wrong counterpart** — measured: IGNITE's scrambled-source loop
**detects in 2 of 5** realisations under the `fN` floor (`dlnL` up to ~15 nat > `fN` = 5.46); IGNITE-2's soft
loop, **6 of 10** scrambled realisations certify at some iteration, **3 keep one to the fixed point** (and
one of those clears even the `fALL` floor by 0.15 nat, inside that floor's own ±0.52 fit error — **even the
blind-robust column leaks this event**).
[MEASURED: IGNITE §5; IGNITE2 §3.1/§3.2]

**S5.6.5 The hole was to be bought back by geometry (D3/D4) rather than amplitude — and it cannot be.**
**STATUS: OPEN, PERMANENTLY STATED, and now known to be STRUCTURALLY UN-CLOSABLE BY CO-REGISTRATION** (S5.4.7).
[MEASURED: project_progress §10.6 item 5 / §10.12 item 3]
[PENDING: QUILL — the hole is currently *stated*, not *priced*. A first-principles criterion must either
supply a defence that is not a co-registration statistic, or convert the hole into an explicit
false-discovery-rate term that a paper can quote.]

### S5.7 What a first-principles treatment should derive

**S5.7.1** Everything in S5 is **empirical**: the floor is *fitted* per cell, its h-scaling is *measured*,
the purity layers were *guessed and tested*. **Nothing here is derived.** The criterion works, and nobody can
say from first principles why it has the shape it has.
[PENDING: QUILL — the first-principles criterion campaign]

**S5.7.2 What QUILL should deliver, in the order the programme needs it:**
1. **Derive `floor(h, T, tol)`** from the E-step's matched-filter cross term, rather than fitting `h^1.66`
   per cell. The mechanism is *stated* (S5.3.1) but never turned into a formula. *A derived floor would need
   no null banks at all.*
2. **Derive the `ln K` ⊕ floor combination** as a single properly-normalised statistic (currently `max(·,·)`
   of a relative bar and an absolute one — an engineering choice, not a likelihood-ratio test).
3. **Adjudicate the D1 fork on decision-theoretic grounds** rather than by scenario-assertion (S5.6.2),
   and **price the wrong-counterpart hole as an FDR term.**
4. **Explain or demolish the one-wall coincidence (S4.3.1).**
5. **State the criterion's operating characteristic** — a proper ROC over (h, T, tol, tier) — so that "onset"
   becomes a curve with an FPR attached rather than a threshold crossing.

---

## S6 — THE ONSET

### S6.1 The retraction

**S6.1.1 IGNITE measured an onset. IGNITE-2 RETRACTED it.** IGNITE's onset table (correct
certifications/realisation, `fN` calibration, **floors = max-of-10 nulls**):

| T | tier | over h = {−13.25, −13.0, −12.75, −12.5} | h\* (>1) |
|---|---|---|---|
| 15 | lit | 0.14, 0.06, 0.16, 0.48 | — |
| 15 | vlbi | 0.16, 0.22, 0.14, 0.52 | — |
| 20 | lit | 0.08, 0.72, 0.56, 0.54 | — |
| 20 | vlbi | 0.72, 0.72, 0.68, 0.38 | — |
| 30 | lit | 0.32, 0.96, **1.54**, 0.94 | **−12.75** |
| 30 | vlbi | **1.16**, 0.78, 0.98, 1.46 | **−13.25** |

[MEASURED: IGNITE §4.2 — 24 cells × 50 Arm-B signal realisations + 30 fresh nulls each]
[SUPERSEDED → **h\* DOES NOT SURVIVE THE D2 ESTIMATOR.**]

**S6.1.2 THE RETRACTION, with its numbers.** Fresh D2-sized floors (**N = 150 nullN, Gumbel MLE,
α = 0.05**) at the two pre-registered onset cells land **8 / 2 nat ABOVE** the max-of-10 floors IGNITE ran
under — **exactly as D2 predicts** (`E[max_10] = μ + β·ln 10` sits at the ~91st percentile with ±1.283β
scatter, while α = 0.05 is an explicit ~95th):

| cell | fresh `fN` floor (N = 150) | `fALL` (N = 270) | banked max-of-10 `fN` | corr certs/real: banked → **fresh** | wrong certs: banked → fresh |
|---|---|---|---|---|---|
| (−12.75, 30, lit) | **38.86 ± 1.47 nat** (μ 19.82, β 6.41) | 65.47 ± 2.26 | 30.89 | 1.54 → **0.92** | 23/50 → **14/50** |
| (−13.25, 30, vlbi) | **7.59 ± 0.48 nat** (μ 1.35, β 2.10) | 11.68 ± 0.52 | 5.46 | 1.16 → **0.54** | 5/50 → **3/50** |

***Re-cut under the honest floors, NEITHER pre-registered onset cell clears onset.*** **`h* = −12.75 /
−13.25` were partly artifacts of the retired floor estimator, and NO properly-calibrated onset exists
anywhere in the modelled box.**
[MEASURED: IGNITE2 §2 — the wrong-cert rate falls with the same stroke, **as it must: purity and count move
together**]
[SUPERSEDED → **the FIRST sentence stands; the SECOND is REFUTED.** SURFACE paid the D2 sizing at all
**108** cells of the extended box: IGNITE's two cells **remain below onset** (0.87 and 0.47 here, on
independent seeds and independent N = 100 null banks — reproducing 0.92/0.54 within the sky error), but
**59 other cells clear it.** *"No properly-calibrated onset exists anywhere in the modelled box"* was a
statement about **two** cells, generalised to twenty-four. **It is superseded, and the generalisation is
the part that was wrong.** → S6.3.2]

**S6.1.3 What SURVIVES the retraction, as measured:** the **relative** structure — onset is
**baseline-driven, not loudness-driven** (`T^{5/2}` fdot/coherence leverage beats the `h^{1.66}` floor race;
**louder alone does not**); **VLBI converts trials mass into detections** (ΣK 88 454 → 470); **no cell
reaches 3 correct certs/real anywhere**; and **under `fALL` the map never ignites at all.**
[MEASURED: IGNITE §4.2/§4.4; project_progress §10.8.2 retraction block]
[SUPERSEDED (in part) → **two of these four survive and two do not.** ✅ The baseline lever survives —
**but it has a CEILING inside the box** (S7.1.1). ✅ VLBI's trials conversion survives — **and it now has a
measured PRICE on the null side** (S7.2.8). ❌ *"No cell reaches 3 correct certs/real anywhere"* is
**FALSE**: SURFACE's 5+11 structure posts **9 cells with corr_lo > 3**, best **7.93/real**. ❌ *"Under
`fALL` the map never ignites at all"* is **FALSE**: **21 cells clear onset on `fALL`, best 2.57/real —
all 21 in the 5+11 structure.** **Both failures are the same failure: IGNITE never varied the population's
loudness STRUCTURE, and that is the axis that moves the count most.** → S6.3.2, S7.6.2]

**S6.1.4 CONVENTION NOW ENFORCED: every quoted onset carries its floor's N and its fit error, or it is not
quoted.** The other 22 cells rest on **10-null floors with ±2–18 nat fit errors** and **cannot support an
onset claim.**
[MEASURED: spec CONVENTION — adopted from IGNITE-2 §2]

**S6.1.5** Four cells still post > 1 under criterion-v2 (best: **1.32 at (−12.5, 30, vlbi)**) — but **every
one rests on an under-sized 10-null floor**, and the single datum we have on what proper sizing does (the lit
onset cell: 30.89 → 38.86 nat, 1.54 → 0.92 certs) says **those floors are biased low by about one of their
own sigmas.** *"Nothing measured contradicts the expectation that paying the sizing everywhere would close
the rest of the box."*
[MEASURED: IGNITE2 §4]
[SUPERSEDED → ***REFUTED. Paying the sizing everywhere OPENS the box.*** SURFACE paid it at all 108 cells
and **59 clear onset.** The expectation was not merely unmet — it pointed the wrong way.]
[~~DISPUTED (D-4)~~ → **CLOSED by SURFACE §7.** N = 100 floors resolve all four: **two RETRACT**
((−13.00, 30, lit) → **0.60** [0.47, 0.73]; (−12.75, 30, vlbi) → **0.73**), **one CONFIRMS**
((−12.50, 30, lit) → **1.13** [1.10, 1.23], floor 106.04 ± 4.62 — **the programme's first confirmed onset
cell**), **one stays MARGINAL** ((−12.50, 30, vlbi) → **1.17**, dies at floor + error). **All four
verdicts survive the criterion-v2.2 floor fix**, though the first cell's numbers moved (its Gumbel was
invalid: zero-fraction 0.27, floor restated 19.46 → 16.60 ± 1.60, count 0.37 → 0.60).
**And the MECHANISM was not what IGNITE-2's single datum suggested.** It is **not** "10-null floors are
biased low": at (−12.50, lit) the properly-sized floor came out **11 % LOWER** and the cell **SURVIVED**.
***A max-of-N floor is not biased in a fixed direction — it is an order statistic with ±1.283β of scatter
and no fixed false-alarm rate, so it lands wherever its ten draws put it.*** That is D2.2's argument, and
this is its first four-cell test.]

### S6.2 The two calibrated cells

**S6.2.1** **The programme currently holds exactly TWO properly-calibrated cells in the entire (h, T, tier,
tol) box:**

| cell | correct certs/real | floor (N, α, fit error) | wrong certs |
|---|---|---|---|
| (h = −12.75, T = 30 yr, lit) | **0.92** | 38.86 ± 1.47 nat (N = 150, α = 0.05) | 14/50 realisations |
| (h = −13.25, T = 30 yr, vlbi) | **0.54** | 7.59 ± 0.48 nat (N = 150, α = 0.05) | 3/50 realisations |

**Both are BELOW the >1 onset bar.**
[MEASURED: IGNITE2 §2, `ig2_floors.npz`]

**S6.2.2 A caveat that travels with the vlbi floor:** at that cell **45 % of `nullN` realisations have NO
offender** — a point mass at zero the Gumbel does not model. The fitted floor lands within 0.5 nat of the
empirical q95, so **the fit is serviceable, but the zero-fraction travels with the number.**
[MEASURED: IGNITE2 §2]
[SUPERSEDED → **this caveat was the edge of a DEFECT, and "serviceable" is not a property that
extrapolates.** ANCHOR §4 measured the zero-fraction across the box and found it reaches **57 %, 80 % and
93 %**. A Gumbel fitted to a point mass at zero is **dragged DOWN toward it**, understating the α = 0.05
bar by up to **2.8×** — **24σ and 12σ against its own quoted fit error.** *The fit error is not merely
wrong; it is confidently wrong: a Gumbel fitted to a 93 % point mass reports ±0.064 nat.* **And it errs in
the DANGEROUS direction** — detection is `dlnL > max(ln K, floor)`, so a floor that is too low is **too
PERMISSIVE.** **This cell's own floor is now RESTATED as the empirical q95 (7.06 ± 0.40, bootstrap), not
the Gumbel (7.59 ± 0.48) — but its count (0.54) and its verdict (BELOW onset) do not move.** → the
zero-fraction convention in the header; S6.5.]

**S6.2.3** Every per-cell rate carries a **±0.2-class sky-sampling error** (15 realisations/cell = 5 skies ×
3 weathers, and GEO says the sky draw dominates yield variance). The **dynamical** statements (S8) are
per-trajectory and do **not** share it.
[MEASURED: IGNITE2 §4.1]

### S6.3 The general surface

**S6.3.1** The onset surface is currently **two calibrated points and 22 uncalibrated ones**. The
`(h, T, tier, tol)` box has never been swept under a converging floor estimator.
[~~PENDING: SURFACE~~ → **DELIVERED.** 108 cells × (30 signal + 200 nulls) = **24 840 realisations**,
≈ 11 GPU-hours, every floor refit per cell at N = 100, α = 0.05. **Deliverables (i)–(v) supplied;
(vi) — the tolerance axis — was NOT run and remains open.** → S6.3.2]

**S6.3.2 THE ONSET EXISTS. `N_onset = 59` of 108 cells.** (MARGINAL 3, below 46 — a MARGINAL cell clears
the bar at its floor and dies at floor + its own error, and is not an onset.) The census-structure onset
sits at **h\* = −12.50** (**1.13** [1.10, 1.23] correct certs/real, floor **106.04 ± 4.62 nat**, N = 100,
α = 0.05, **zero-fraction 0.00**) — **one grid step LOUDER than IGNITE's retracted claim, not absent.**
***IGNITE's max-of-10 floors were drawn low enough to move the onset into the WRONG CELL — a subtler and
more instructive failure than "there is no onset".***
[MEASURED: SURFACE §3/§4, re-cut against the criterion-v2.2 floors; `recut_surface.npz`]

**S6.3.3 THE 59 IS A RE-CUT NUMBER, AND ITS STABILITY IS A COINCIDENCE.** The floor fix (S6.5) touched
**15 of the 108 cells**; **93 — including 57 of the 59 onsets — are bit-identical to the published
surface.** Of the touched cells, **two onsets DIED and two were BORN.** The pre-fix banked count was also
59. ***That the totals agree is a coincidence, not a confirmation. The number is stable; the map is not.***
[MEASURED: RECUT §1; gates A and B both 0.000e+00 — gate B proves the re-cut scorer IS the scorer that
produced the published surface]

**S6.3.4 `h*` IS NOT BOUNDED BELOW in 7 of the 18 frontier columns.** For those, SURFACE locates an onset
but **cannot locate its faint boundary** — the faint edge of the surface is **outside this grid**, and the
next grid is **fainter, not louder.** **Re-cut: still 7 of 18 — and NOT the same seven** (lost: lit 3+13
T = 40, lit 2+14 T = 50; gained: lit 3+13 T = 50, vlbi 2+14 T = 50; five unchanged). ***Quote the number;
never the same seven columns.*** *This is the number-not-names lesson (S3.3.1) arriving a second time, for
the same reason.*
[MEASURED: SURFACE §4 / RECUT §1.3]

**S6.3.5 A caveat on the frontier, stated because it is load-bearing for KINDLE.** Two of the faint-edge
onsets are **ISOLATED**: the cell certifies at h = −13.25 while **every louder h in the same column sits
below the bar.** This is **not** a re-cut artifact — the floor grows as ≈ `h^1.5–2`, so **a faint cell can
out-certify a loud one**, and the *published* surface already showed the pattern. **Whether an isolated
faint-edge onset is a frontier or a fluctuation is UNMEASURED.**
[MEASURED: RECUT §1.3 caveat] [PENDING: KINDLE item (ii)]

### S6.4 The real-data anchor

**S6.4.1** Every onset number in the programme is **simulated end to end**: simulated array, simulated noise
draws, injected sources, an extended baseline built by a **stated convention** (median-cadence TOAs, smooth-
basis timing-design extrapolation, zero-extended binary columns, span-scaled GP components) that is
**explicitly not a forecast of real future data** (no new backends, no DM events, no solar-wind realism).
**Nothing in the certification chain has ever touched a real TOA.**
[MEASURED: IGNITE §2 conventions + §7 caveats — the T-extension is gated as a strict no-op at dT = 0
(`lnL(truth|zero-noise) = 405413.512739`, bit-equal to the banked value)]

**S6.4.2** The one place a real number enters is the **prior**: `best_distances.txt` (frozen, canonical,
git-tracked), plus two hand-injected published composites (J0437 Reardon+2016, J1909 Reardon+2021).
[MEASURED: Anchor Census A0]

[~~PENDING: ANCHOR~~ → **RAN, AND ITS PREMISE WAS REFUTED AT TASK 1.** **The repo has no real residuals.**
The 116-pulsar feather set is a **MOCK** — telescope `AXIS`, a single 1440 MHz channel — and **its
`residuals` column IS the injected CW + CURN realisation `b20_cw_curn_r0`**, bit-identical (max|diff| = 0.0,
all 116). **The certification chain never reads it** (`data = inject_delay(θ_true) + NoiseDrawer.draw(seed)`),
so **no banked result depends on it** — ***but any future task that treats those residuals as "the data" is
measuring an injected CW.*** ANCHOR's brief was therefore amended in flight: **the anchor is REALISM, not
provenance** → S6.4.3. **The real-data question it was built to answer is NOT answered, and is re-tagged
[PENDING: REAL-ARRAY]** below.]
[MEASURED: ANCHOR §0, `anchor_data_forensics.npz`]

**S6.4.3 THE REALISM LADDER: CONSISTENT AT EVERY RUNG. THE ONSET SURFACE DOES NOT SHIFT.** Regenerate
noise-only nulls on the *same* array with noise the frozen analysis does **not** model — **R1** RN/GWB
spectral mis-specification, **R2** an unmodelled DM power-law GP, **R3** non-Gaussian tails — and the
realised noise rms rises **1.993 → 2.634 µs (×1.32)**. `NoiseDrawer.draw` is the **only** object replaced;
every other frozen object is untouched. ***The mismatch IS the experiment.*** Pooled over all 40
(rung × cell) comparisons: **18/40 cells "inflated" — a coin flip**; median Δq95 = **−0.18 nat**; **1/40**
significant under Mann-Whitney against **2 expected by chance.** **Verdict: CONSISTENT at every rung.**
[MEASURED: ANCHOR §3 — 7 200 realisations, T = 15 (the native, unextended array), N = 150 nullN per cell]

**S6.4.4 AND THE MECHANISM, WHICH IS A STRONGER STATEMENT THAN THE NULL RESULT.** Exactly 1 of 40 cells is
significant — and its **floor barely moves** (Δq95 = +0.36 ± 0.29) while **the rate at which the null
produces a candidate AT ALL more than doubles** (P(offender > 0): 0.200 → 0.453). ***Unmodelled noise makes
the null throw up MORE candidates, but not LOUDER ones.*** The floor is an **upper-tail** quantile, and the
upper tail is set by the **TEMPLATE, not the data**: the E-step's matched-filter cross term is linear in the
*model* amplitude, so the loudest fringe fluctuations scale with the hypothesis `h` — the `floor ∝ h^1.66`
mechanism (S5.3.1), which ANCHOR reproduces **from scratch at 1.67**, at a baseline where it was never
fitted.

> ***The criterion's floor is robust to noise mis-specification BECAUSE it is loudness-relative.*** The
> property S5.3.2 recorded as a **COST** — the bar rises almost as fast as the signal, so louder alone buys
> nothing — is **the same property that makes the floor immune to getting the noise model wrong.**
> **The tail is template-dominated; the body is noise-dominated.**

[MEASURED: ANCHOR §5]

**S6.4.5 THE OFFENDER ANATOMY IS SET BY THE TRIALS FACTOR, NOT THE NOISE BUDGET.** The same five pulsars
carry ~75 % of offenders in **every** rung, including the control; total offenders move **+6 %** from R0 to
R3 against a **+32 %** rise in noise rms. **The suspects did not show up:** J1824−2452A — the array's
reddest pulsar by far (feather χ²/N = 3243) — contributes 6 offenders in R0 and 13 in R3; B1937+21 never
appears. **Unmodelled DM and unmodelled red power do NOT route false alarms through the chromatically-dirty
pulsars.** What *does* move is **J0437−4715 (+31 under R3 — the largest single shift in the table)**, and
J0437 is **the array's smallest trials factor** (`ln K` = 1.39, the array minimum). ***The pulsar whose bar
is lowest is the one most exposed when you add noise the model does not know about.*** *That is the J0437
double edge (S5.2.4) — robustness to source error and vulnerability to noise are THE SAME PROPERTY —
showing up in a third independent place.*
[MEASURED: ANCHOR §6]

**S6.4.6 THE LADDER'S CEILING, stated because it BOUNDS S6.4.3.** **The array is single-frequency** — all
30 225 TOAs sit at 1440.0 MHz, so the chromatic factor `(f_ref/f_obs)²` is a **constant** (0.9450–0.9454)
and **chromaticity is UNIDENTIFIABLE BY CONSTRUCTION**: a DM GP enters the data as an extra *red* process,
and no analysis could separate it from an achromatic one. ***R2 tests UNMODELLED CHROMATIC-BAND POWER, NOT
CHROMATICITY.*** (ECORR is likewise degenerate with EQUAD here — one TOA per epoch.) **This is the ladder's
binding limitation, and only a real multi-frequency array lifts it.**
[MEASURED: ANCHOR §2.1]

[PENDING: **REAL-ARRAY** — port the criterion stack to a real PTA. ***Why it is a campaign, not a task:***
every prior in the criterion is keyed to the 116-pulsar mock — `best_distances.txt`, the per-pulsar
`ln K_counted`, ΣK, the lit/vlbi tiers, the geometry ensembles, the census, and the `NoiseDrawer`'s
hyperparameters. ***A real-array anchor is a RE-DERIVATION OF THE PRIOR STACK on a different array, not a
substitution of a residual vector.*** **On disk and verified (ANCHOR 2026-07-12):** NG 15 yr (66 psr,
615 294 TOAs, 16.03 yr, 1705 radio frequencies — loads in `discovery` unmodified); NG 20 yr (77 psr,
1 131 412 TOAs, 20.62 yr); MPTA DR3 (83 psr, par/tim only). **Deliverables:** (i) the prior stack
re-derived on the target array; (ii) the floor re-fit against that array's own noise; (iii) the programme's
first certification numbers that touch a real TOA; (iv) **it LIFTS ANCHOR's ceiling** — a multi-frequency
array makes chromaticity *identifiable*, so the DM channel can be tested as chromaticity rather than as
unmodelled red power (S6.4.6).]

### S6.5 The floor estimator has a validity domain — and outside it, it errs PERMISSIVE

**S6.5.1 THE DEFECT.** The offender statistic is **0.0 whenever a realisation has no cell passing layer 1 ⊕
layer 3.** At faint `h` that is *most* realisations — measured zero-fractions reach **57 %, 80 % and 93 %**.
**A Gumbel fitted to a point mass at zero is dragged DOWN toward it**, understating the α = 0.05 bar by up to
**2.8×** (0.845 fitted vs **2.395** empirical at h = −13.25 lit) — **24σ and 12σ against its own quoted fit
error.** ***And it errs in the DANGEROUS direction:*** detection is `dlnL > max(ln K, floor)`, so **a floor
that is too low is TOO PERMISSIVE — it lets pure-noise offenders through.** *Every faint-`h` cell in the
onset map was calibrated by this estimator.*
[MEASURED: ANCHOR §4 — found in the CONTROL arm of a campaign built to test something else]

**S6.5.2 THE CONVENTION (criterion-v2.2, binding).** ***The D2 Gumbel floor is valid ONLY where the nullN
zero-fraction is ≲ 20 %. Above that, quote the empirical (1−α) quantile with a BOOTSTRAP error, and bank the
zero-fraction beside it. The zero-fraction is a REQUIRED column, not a caveat.*** The onset **test** is
unchanged — `ONSET` iff the count at **floor + its own error** exceeds 1 — only the floor and the error
change.
[MEASURED: RECUT §0; adopted verbatim from ANCHOR §4]

**S6.5.3 WHAT IT COST, AND WHAT IT DID NOT.** **SURFACE: 15 of 108 cells touched, N_onset = 59 unchanged in
total and changed in membership** (S6.3.3). **CHORUS: ALL 26 of 26 cells fail the gate; 23 of 26 floors
rise; and the campaign's loudest headline dies** (S7.6.4). **ANCHOR: nothing retracted.** **IGNITE-2: both
calibrated cells keep their verdicts; the vlbi floor is restated** (S6.2.2). ***The fix bit exactly where the
provisional analysis predicted the exposure was — faint `h`, high zero-fraction — and nowhere else. 57 of the
59 onsets, the entire loud-cell result, both purity rejections, and the clock-sharing verdict never depended
on the estimator at all.***
[MEASURED: RECUT §6]

**S6.5.4 THE RULE THIS ESTABLISHES, and it is the durable part.** ***Anything of the form "cell A versus cell
B" — a difference, a ratio, an inversion, a peak location — must be RE-DERIVED from the corrected banks
before it is quoted. The counts moved, so the differences moved.*** **Two such claims were checked rather
than assumed after the fix, and BOTH moved**: SURFACE's `h*`-column membership (S6.3.4) and CHORUS's
trade-curve inversion (S7.6.5). *Neither was on anyone's list.*
[MEASURED: RECUT §6]

---

## S7 — THE LEVERS

*Ordered as IGNITE-2 measured them: T first, distances second, loudness third, eccentricity fourth,
geometry fifth — and the one premise nobody has tested, sixth.*

### S7.1 T — the strongest lever (`T^{5/2}`)

**S7.1.1** **`T^{5/2}` fdot/coherence leverage beats the `h^{1.66}` floor race. Louder alone does NOT.**
Onset appears at T = 30 yr and **nowhere at T = 15 or 20** anywhere in the box. This is the single clearest
lever statement in the programme.
[MEASURED: IGNITE §4.1/§4.2; project_progress §10.10(d)]
*(Note the honest asymmetry: T = 30 is where the retracted onset lived. Under fresh floors those cells read
0.92/0.54 — **the T-lever's DIRECTION survives the retraction; its sufficiency does not.**)*

**S7.1.1a THE STRONGEST LEVER HAS A CEILING, AND THE CEILING IS INSIDE THE BOX.** Over the 36
(h, tier, structure) columns of the extended grid: **T = 30 yr is optimal in 0 of 36.** T = 40 wins 19,
T = 50 wins 17 — **and the split is h-dependent, not random.** ***Loud cells PEAK at 40 yr and LOSE at 50***
(12 of 12 loud columns of the 3+13 and 5+11 structures), while faint cells keep gaining to 50 (6 of 6 loud
2+14 columns rise). **The mechanism is visible in the floors:** the counterpart-matched floor grows with
**data volume** as well as with loudness — the matched-filter cross term integrates more data — and between
40 and 50 yr **that growth OVERTAKES the `T^{5/2}` leverage** at loud `h`, while at faint `h` the floor is
still small enough that the leverage wins. **IGNITE's *"onset is baseline-driven; `T^{5/2}` beats the floor
race"* is TRUE UP TO A CEILING, and the ceiling is ~40 yr: past it, the floor race resumes and WINS.**
[MEASURED: SURFACE §9.2, re-cut — the 0/36, 12/12 and 6/6 counts are all re-derived on the corrected floors]
*(Convention caveat, load-bearing: the T = 50 cells extrapolate the timing model **35 yr** past the last real
TOA under a stated convention — median-cadence TOAs, smooth-basis design extrapolation, zero-extended binary
columns, span-scaled GP components. **The SIGN of the 40 → 50 fall is robust; the MAGNITUDE is a property of
the convention and must never be quoted as a forecast.**)*

**S7.1.2** The T-lever has an independent expression on the blind side: the float ceiling reaches the F2
loosest rung only at **T ~ 11 kyr** (S4.1.13). **The same lever, four orders of magnitude out of reach for
the blind problem, is merely "a longer PTA" for the targeted one** — which is the clearest single statement
of what the counterpart buys.
[MEASURED: spec R5 / PENCIL]

### S7.2 VLBI-tier distances — now TWO independent roles

**S7.2.1 ROLE 1 — VLBI buys DETECTIONS, by converting trials mass.** The tier reduces the union-18's trials
mass **ΣK = 88 454 → 470**, i.e. it converts trials into detections **at fixed `dlnL`**. *"VLBI's gain is
exactly where RING said it would be."*
[MEASURED: IGNITE §4.4 / gate IG4]

**S7.2.2 ROLE 2 — VLBI buys WRONG-COUNTERPART ROBUSTNESS, via the oracle/implementable gap.** D4-OBS (the
implementable form) is **1.6× stronger** in the VLBI tier (51.6 %) than in the lit tier (32.9 %), because the
gap closes with `σ_d`. ***Sub-3-pc distances are DOUBLY LOAD-BEARING: they buy detections AND they buy
wrong-counterpart robustness.*** (Role 2 is new with D4 and did not exist in the earlier design theorem.)
[MEASURED: D4 §6.3]

**S7.2.3 THE LADDER IS BINARY, NOT GRADED — and the Gaia tier is a NO-OP.** The coherence factor
`κ = exp(−½[2ω₀(1−cos μ)(kpc/c)σ_d]²)` demands **σ_d < 3.02 pc** for κ > ½ at `fgw = 1e-8`. Exactly **1 of 30**
ring pulsars qualifies (J0437, σ_d = 1 pc). **Gaia's factor-1.6 moves κ̄ from 0.0433 → 0.0550 and buys
1.0–1.8× in area while DEGRADING coverage.** ***A factor-1.6 distance improvement is worth nothing for CW
localisation at fgw = 1e-8.*** What matters is **crossing the coherence threshold — a VLBI / timing-parallax
regime, not a Gaia one.** (At `fgw = −9` the threshold relaxes 10× to σ_d < 30 pc, 5 of 30 qualify, and the
ladder becomes genuinely graded.)
[MEASURED: RING §0/§3.4/§5.4]

**S7.2.4 BAD DISTANCE PRIORS BIAS THE SKY — THEY DO NOT MERELY BROADEN IT.** At `fgw = −8` every non-exact
prior drives the sky MAP **3–6° off truth, INDEPENDENT of SNR**, while the 90 % area shrinks 4–17× per SNR
doubling. **Coverage therefore DEGRADES as the signal gets louder:** `inside90` = **0.90 → 0.50 → 0.00** at
SNR 5/10/20. ***Proven, not inferred***: the zero-noise control gives bias **2.73–5.28°** for every imperfect
tier and **exactly 0.0000°** for the exact tier, in all four configurations. **Mechanism, isolated: bias ∝
un-modelled pulsar-term power, and nothing else** (bias collapses to **0.033°** only when κ̄ reaches 0.290 —
a 54× reduction).
[MEASURED: RING §0/§3.5/§5.1 — 10 realisations/cell]

**S7.2.5 The worst regime is PARTIAL distance knowledge.** tier1 → tier2_k3 at SNR 10 shrinks area90 **6.5×**
while the zero-noise bias barely moves (2.73° → 3.07°) — the credible region contracts around an offset that
does not, so `inside90` falls **0.60 → 0.30**. **Tightening a prior that is still far from the coherence
threshold buys precision and pays for it in coverage.**
[MEASURED: RING §0]

**S7.2.6 RING's five stop-points travel with any future Q1 work** (cited, not restated):
**S-1** a coarse+zoom grid search **refines the wrong peak** (61 % of realisations; **100 % at SNR 20, in
every scenario**) — *a grid will not do*; **S-2** the tier ladder is binary at fgw = −8 and the Gaia tier is a
no-op; **S-3** the harness's timing-model prior is **internally inconsistent** (enterprise `1e-14·N_toa` vs
discovery `1e-14·N_toa/N_par`, factor ≈ 19) — harmless at fgw = −8 (`Var(z) = 1.17`, KS p = 0.73) but
**breaks fgw = −9** (`Var(z) = 6.42`, KS p = 0.0000): ***treat every `fgw ≲ −9` NOISY result from this harness
as uncalibrated***; **S-4** with a GWB, **even exact distances undercover** (`inside90` 0.40–0.50 vs nominal
0.90) — physics vs estimator error **not separated**, so **do not quote the scenario-C tier3 coverage
number**; **S-5** RING ran the feather σ_d, not the canonical `best_distances.txt` — **impact on conclusions:
none** (κ̄ 0.043 → 0.033, same coherent set), **impact on what it enables: large** — the canonical *means*
differ by **1.40σ on J0437** (a **0.55 rad** pulsar-term phase error), which is **a ready-to-run mis-centred-
prior arm** and the single highest-value follow-up RING points to.
[MEASURED: RING §7]

**S7.2.7** A structural limit of the RING harness, stated because it bounds every number above: **the prior
mean is ALWAYS exactly the true distance.** A "bad prior" there is **wide**, never **mis-centred**. Since the
width axis *alone* already produces a 3–6° SNR-independent bias, **a mis-centred prior can only be worse.**
[MEASURED: RING §2(a)/§2.1]

**S7.2.8 ROLE 3 — AND IT IS A COST. VLBI IS NOT FREE IN ONSET-MAP UNITS.** `floor(vlbi) − floor(lit)` at the
same `h`: **+2.9 ± 1.0 nat at h = −13.25, and nothing measurable anywhere else in the box** (the loud cells'
floor errors are 3–6 nat and cannot resolve a price of this size, so **none is quoted**). **It is a TIER
effect, not a realism effect** — the control arm already carries +2.59 ± 1.41, statistically
indistinguishable from the worst realistic rung's +2.76. **Realism does not change the price of VLBI.**
**The mechanism, and the sign is not a paradox:** VLBI shrinks `σ_d` → fewer fringes in the prior window →
smaller `K_counted` → **a LOWER trials bar** → **a pure-noise fluctuation clears layer 1 more easily** →
more and louder offenders → **a HIGHER absolute floor.** ***This is the J0437 double edge again, measured as
a TIER-level quantity: VLBI buys detections on the SIGNAL side (ΣK 88 454 → 470) and pays for them on the
NULL side.*** **This is the first number for what a VLBI campaign costs.**
[MEASURED: ANCHOR §7 — empirical q95, bootstrap errors, across all rungs]

### S7.3 Loudness — the SCOUT clock

**S7.3.1** **Loudness is a lottery, not a lever.** The floor rises almost as fast as the signal (S5.3.1), so
the certified count is **non-monotone in h** (S5.3.2).
[MEASURED: IGNITE §4.1]

**S7.3.2 SCOUT's census: NO named credible SMBHB candidate clears the operating point at its
literature-favoured mass.** OJ 287 (the best-studied binary, and the only secure one) is **~1.5 dex too
quiet** because its mass ratio is extreme (q ≈ 0.008): `log10 h = −15.86`. The three that come within 1 dex —
**3C 66B, PG 1302−102, PKS 2131−021** — do so **only at the loud edge of their mass range**, and that edge is
in each case **either already excluded by PTA non-detection** (3C 66B's Sudou mass; PKS 2131 at its
`Mc < 1.5e10` upper limit) **or the optimistic tail of a contested periodicity** (PG 1302).
[MEASURED: SCOUT §1/§2]

**S7.3.3 THE STRUCTURAL SCISSORS.** The `log10 h = −13.75` operating point sits **~0.35 dex ABOVE the NG15
sky-averaged CW upper limit (−14.1)**. **Any source genuinely at the operating point, in band, would already
be an individual-source DETECTION.** *"Certifiable and counterpart-identified" is asking for a source louder
than everything current PTAs have failed to detect.* **This is not a counterpart-availability problem; it is
a source-loudness problem.**
[MEASURED: SCOUT §2]

**S7.3.4 The population clock.** GWB-anchored population synthesis gives the mean number of individually
detectable SMBHBs as **N̄ ≲ 0.01 (current), ≲ 0.1 (any SMBH distribution consistent with the SGWB)**;
**O(0.1–1) in the SKA era.** **The joint (loud AND counterpart-identified) count is NOT suppressed by the
counterpart requirement** — EM astrometry beats the 22-arcsec tolerance by **≥10⁴×** for any VLBI/Gaia/optical
counterpart. ***The counterpart doesn't cost you; the loudness does.*** `N_joint ≈ 10⁻²…10⁻¹` (current).
[MEASURED: SCOUT §3/§4]

**S7.3.5 And IGNITE's onset is 0.5–1.0 dex LOUDER STILL than SCOUT's −13.75 operating point.**
[MEASURED: IGNITE §6 — see S11]

### S7.4 Eccentricity — the ATLAS corner and κ(e)

**S7.4.1 THE SELF-CLOCKING CORNER: eccentric (e ≳ 0.6), massive (Mc ≳ 10⁹ M⊙), HIGH IN THE BAND
(f_orb ≳ 10⁻⁸ Hz).** At `(f_orb = 10⁻⁸, Mc = 10⁹)` the comb self-clocks — σ(log10_mc) improves >20× — from
**e ≈ 0.58**; the threshold rises to **e ≈ 0.70–0.84** at lower Mc or `f_orb = 10⁻⁸·⁵`. ***Below
`f_orb = 10⁻⁸` the comb is buried in the red-noise/GWB band and NEVER self-clocks, at any e.*** **The first
source must live at the TOP of the band.**
[MEASURED: ATLAS M2 — 63 cells, e × Mc × f_orb; gates G0/G0-T1 (e = 0 stack == circular source to 7.14e-15)
and G1 (e → 0 reproduces SIREN's Earth-term σ to 1e-5–1e-3)]

**S7.4.2 THE QUALIFYING STATEMENT — the corner clears the RELATIVE bound, not the ABSOLUTE one.**
`κ_measured ≥ 20` by e ≈ 0.5–0.65, and marginal σ(log10_mc) reaches **0.008–0.02 dex** (a 40–115×
improvement) — which clears the **>20× relative** bound (S4.2.10) but **does NOT clear the ~0.003-dex
ABSOLUTE certification floor** (best in the valid tier: **0.0075 dex, 2.5× short**). **Two distinct criteria
that earlier work conflated. Keep them apart.**
[MEASURED: ATLAS M2 — and note this is a THIRD "20": WEAVE's `κ ≥ 20` self-clock threshold is *neither* of
the other two. See S4.3.1.]

**S7.4.3 κ IS FREQUENCY-DEPENDENT — which the white-noise `(n_eff/2)F(e)` scaling cannot capture.** At
`f_orb = 10⁻⁸` measured κ tracks the analytic (11.0 vs 5.6 at e = 0.5). At **10⁻⁸·⁵ it VASTLY EXCEEDS** it
(2216 vs 172 at e = 0.8) — the comb's higher harmonics reach further into the sensitive band. At **10⁻⁹ it is
BELOW analytic** at moderate e (0.92 vs 5.6) — the fundamental sits where red noise buries it, so spreading
power into the comb **reduces** the chirp Fisher until e climbs out (601 at e = 0.8).
[MEASURED: ATLAS M1]

**S7.4.4 THE THROTTLE (the honest surprise): eccentricity's value is CONDITIONING, not magnitude.** The
*conditional* chirp Fisher κ is enormous (~33 000× at e = 0.9) but the *marginal* σ(log10_mc) improves only
~40–115×: **the comb's chirp information is largely DEGENERATE with `e` and `f_orb`.** Only at high e does the
comb geometry (tooth spacing → f_orb, amplitude ratios → e) break the degeneracy.
[MEASURED: ATLAS M1]

**S7.4.5 The clock-cancellation ceiling BREAKS at high e — but ignition is marginal and shortest-lag-only.**
Broken honestly (Peters `e(t)` decay + 1PN periastron advance γ̇, RK4, autodiffed): at moderate e (0.3–0.65)
`R_rank3/R_scalar < 1` (the cancellation holds, as WEAVE predicted); at high e (0.8–0.9) it **BREAKS**, ratio
up to **41.6**. But ignition (`R_a ≥ 1`) is reached **only at `τ_a = 0.3 kyr`** (the nearest pulsars) — 3 cells
at f_orb = 10⁻⁹, 4 at 10⁻⁸·⁵, **max R = 3.53**; **at `τ_a ≥ 1 kyr`, NO ignition.** ***A refutation of the pure
cancellation at high e, plus a marginal, shortest-lag-only, CEILING ignition — not a clean null, and not a
robust cascade ignition.***
[MEASURED: ATLAS M3]

**S7.4.6 THE EOB-TIER VALIDITY LIMIT — the map cannot see the corner that matters.** The comb is a **toy tier**
(circular-kernel harmonic stack, fdot tie only). The `F(e)`-boosted `mc_n` makes the chirp term go negative —
the harmonic "coalesces" within the span — at the extreme e × Mc × f_orb corner; the clip is **binary** (a cell
is fully valid, or its whole comb coalesces → **TOY-TIER INVALID**, flagged, its κ **not** read as "not
self-clocking"). **5 of 63 cells flagged; 0 dropped.** ***Decisively: the cells that would clear the 0.003-dex
absolute floor (e ≳ 0.85) are EXACTLY the toy-invalid ones. The map's most important corner is the corner it
cannot see.***
[MEASURED: ATLAS "THE EOB-TIER VALIDITY LIMIT"]
[PENDING: EOB tier (arXiv:2511.19611) — **not a refinement; the load-bearing next step.** It must deliver
σ(log10_mc) at e ≳ 0.85 in the valid regime, and thereby settle whether the absolute floor is reachable at all.]

**S7.4.7 LEVER (ii) WAS REDEFINED MID-PROGRAMME, and the trail matters.** The design theorem's original
lever (ii) was **"wide lanes from eccentric harmonics"** — the *vernier* idea, that the orbital fundamental
supplies a registration lane `n_peak×` wider than the highest harmonic, populating the empty
`2e-3 → 0.05` band. **LANES tested it and it is a NO-GO.**
[MEASURED: LANES BRIDGE VERDICT — widest *usable* lane peaks at **8.6e-3 scaled at e = 0.9** (power-anchor),
still **5.8× short** of the 0.05 float ceiling; physically (residual-SNR anchor) the widest lane is
**1.85e-3–3.7e-3 — essentially the F2 rung itself**, 13–27× short. **No optimal-e window bridges.** The
negative is a **scissors**: (1) to sit the power-dominant harmonic at 27× needs e ≳ 0.85, but there the
fundamental falls below 5 % of peak and is unusable; (2) a PTA detects via *timing residuals*, so per-harmonic
`SNR² ∝ g(n,e)/n⁴` and **the SNR-dominant harmonic is ALWAYS n = 1–2, regardless of e** — the detection band
*is* the fundamental, so F2's rung already sits there and **there is no wider lane above it to climb to.**
***Eccentricity buys finer rungs (precision), not wider rungs (capture). It densifies the ladder Track B
already has; it does not extend it upward.***]
[MEASURED: spec B1 STEP 2 — lever (ii) was then **re-posed** as *likelihood structure* (κ), not lane width,
and in that form it is the ONLY thing that moves the targeted evidence (S4.2.12), which is what ATLAS then
measured. **Both statements are true of different objects. Do not quote LANES's NO-GO against ATLAS's corner,
or vice versa.**]

**S7.4.8** The eccentric-comb-as-distance-vernier angle is **NOVEL / UNCLAIMED** in the literature (Taylor+2016
explicitly excludes the pulsar term; the EOB paper maps the harmonic structure but states it cannot determine
the pulsar phase given poor distance knowledge; no one joins comb → distance). **Worth a footnote of priority
even though the toy says it does not bridge.**
[MEASURED: LANES "Literature"]

### S7.5 Geometry

**S7.5.1** **The sky draw dominates yield variance** (S3.3.8), the selection function is measured (S3.3.5–6),
and **the certified names are a measure-zero outcome** (S3.3.1). Geometry is a lever only in the sense that it
is a **variance** to be planned around — sky-conditional seed sets (S3.3.4) — not a knob anyone can turn.
[MEASURED: GEO; FORGE §6]

### S7.6 Mixed populations — THE ONE PREMISE NOBODY HAS TESTED

**S7.6.1** ***Every result in this repo is for a SINGLE-POPULATION source model. Nature supplies a MIXTURE.***
[MEASURED: project_progress §10.12 — flagged by IGNITE-2 as the queue head]
[~~PENDING: CHORUS~~ → **DELIVERED.** 26 mixture cells × 30 realisations + 4 000 nulls + 40 exact pairs +
30 soft loops. All three deliverables supplied. → S7.6.2–S7.6.6]

**S7.6.2 CERTIFICATION IS A PROPERTY OF THE POPULATION — measured, on TWO independent axes.**
- **LOUDNESS structure (SURFACE).** Promoting two of sixteen sources faint → loud (3+13 → 5+11) at fixed
  (h, T, tier) raises the count **up to 6.1×** (median **2.5×** across the 36 columns) — *super-linearly,
  and against a floor that itself rises 2–3× because the recovery model carries more amplitude.* Demoting
  one all but extinguishes certification. **Two loud sources is a dead population; five is a factory.**
  **The frontier moves by ≥ 0.75 dex.**
- **WAVEFORM structure (CHORUS).** Replacing ONE of the three loud members with an eccentric source at fixed
  census loudness moves the count **14.8× (lit) / 12.4× (vlbi)** at e = 0.7. ***Eccentric structure is the
  strongest single lever yet measured in the box.***

> ***Every single-source no-go in this repo is scoped to a population structure nature does not have to
> supply.*** *(Both levers are re-derived on the criterion-v2.2 floors, not inherited — S6.5.4. The two were
> measured on different floors, as they must be: each mix refits its own null.)*

[MEASURED: SURFACE §9.1 / CHORUS §1, re-cut]

**S7.6.3 BUT THE CLOCK IS NOT SHARED — and that was the campaign's reason to exist.** In **20 exact
seed-paired** realisations (e = 0 vs e = 0.7, sharing the noise seed AND the physical truth distances, each
banking the member0-inert-template rescore `dlnL_ct`), the certified count jumps in **every** pair (lit
0–3 → 1–15; vlbi 0–1 → 1–15). **But under the pre-registered attribution rule — ecc-attributed iff
`dlnL − dlnL_ct > 1` nat — ZERO of the ~120 lifted certifications across all 20 pairs are
circular-attributed.** Every pulsar the eccentric arm certifies is certified **through the eccentric
member's OWN comb template** — its harmonic pulsar terms — not through a lifted circular-member
registration. **Pre-registered verdict: MARGINAL (lit, 1/10 pairs above floor-refit noise) / ABSENT (vlbi,
0/10).** **And the joint-fit channel — where SIREN's lag-diversity mechanism (S9.2) would have to operate —
adds NOTHING:** all 20 signal soft-loop trajectories are **FLAT**; the eccentric-seeded loops hold their
large seed sets (18/6/2/2/6) exactly as the circular-seeded loops hold their sparse ones (0/0/0/0/1).
***The eccentric-seeded loop does not consolidate further than the circular one. Its entire advantage is in
its SEEDS.***

> ***Certification in a mixed population is a property of the SOURCE, not a shared array resource.***
> **The single-source no-gos are rescoped anyway — but through the comb's OWN pulsar terms, not through a
> shared clock.** *One moderately eccentric loud member certifies up to ~18 pulsars single-handedly where
> the entire circular population certifies none.*

[MEASURED: CHORUS §2 — **structural, and independent of the floor: the re-cut does not touch it**]

**S7.6.4 THE SWITCH-ON THRESHOLD IN e IS A MIXTURE PROPERTY — and CHORUS's published threshold is REFUTED.**
[SUPERSEDED → *"every eccentric mix clears the >1 bar; **a single e = 0.3 member suffices at either tier**"*
(counts 1.57 lit / 1.13 vlbi) — **FALSE.** Under the criterion-v2.2 floors the lit cell collapses to **0.70
— below the bar, REFUTED** — and the vlbi cell reads **1.03: MARGINAL, not confirmed** (it clears the bar at
the floor and **fails at floor + bootstrap error**, 0.60). The lit floor rose **7.39 → 11.30 nat: +53 %, a
6.2σ move against its own quoted fit error** — a Gumbel fitted to a **73 %-zero point mass.** ***This is the
single most expensive consequence of the floor defect in the whole repo.***]

> ### THE CORRECTED, BINDING EXTERNAL STATEMENT
> ### **With ONE eccentric member, the switch-on is at `e = 0.5`. With TWO OR MORE, it is at `e = 0.3` (CONFIRMED, both tiers).**
> ***The threshold is NOT a property of eccentricity alone — it depends on how many members carry it.***

| n_ecc | tier | e = 0.3 | e = 0.5 | e = 0.7 | switch-on |
|---|---|---|---|---|---|
| **1** | lit | 0.70 [0.43] **below** | 3.13 [2.70] CONFIRMED | 5.43 [4.90] CONFIRMED | **e = 0.5** |
| **1** | vlbi | 1.03 [0.60] MARGINAL | 2.27 [1.73] CONFIRMED | 5.77 [5.13] CONFIRMED | **e = 0.5** |
| 2 | lit / vlbi | 2.77 / 1.77 **CONFIRMED** | 4.90 / 3.97 | 5.47 / 4.10 | **e = 0.3** |
| 3 | lit / vlbi | 2.50 / 2.20 **CONFIRMED** | 5.83 / 4.50 | 4.07 / 5.07 | **e = 0.3** |

*(count at the adopted floor, [count at floor + bootstrap error]. CONFIRMED = clears the >1 bar at floor +
error; MARGINAL = clears it only at the floor.)*

**Note the direction this moves the programme:** the corrected single-member threshold (e ≳ 0.5) **pulls the
count criterion's requirement BACK TOWARD ATLAS's own self-clocking corner (e ≳ 0.58–0.6, S7.4.1)** — two
thresholds measured on **entirely different statistics**, which now nearly agree. *That is a harder
requirement than published CHORUS claimed, and a more coherent one.*
[MEASURED: CHORUS §1 / RECUT §2.1, `recut_chorus.npz`]

**S7.6.5 THE CAPACITY-VS-CLOCK TRADE CURVE — the deliverable exists; the CROSSING does not.**
[SUPERSEDED → *"the trade inverts between n_ecc = 2 and 3 at high e; the surface peaks at n_ecc = 2
(8.7/7.4); the capacity crossing sits at ~8–12 % band occupancy"* — **DEMOTED TO NOT-CLEAN.** Under the
corrected floors the inversion **flips status in 3 of the 8 (e, tier) combinations, and not in the same
direction**: published 5 of 8 inverted, re-cut **4 of 8 — but not the same four.** The surface no longer
peaks in the same place either: **vlbi peaks at n_ecc = 2, lit now peaks at n_ecc = 3.** **DO NOT QUOTE "the
trade inverts at n_ecc = 3" without re-deriving it from `recut_chorus.npz`. It is not refuted; it is no
longer a clean claim.** *This is a difference between two counts that BOTH moved — exactly the class S6.5.4
names.*]
**What SURVIVES is the mechanism, and it is floor-independent:** the binding cost at high occupancy is the
**`K_counted` trials term** (K_sum grows ~11× from m0 to m3e07) and the finer joint fringe grid — **not the
floor**, which *falls* at the n_ecc = 3 high-e cells while the count falls with it. ***The capacity ceiling
is real. WHERE it bites is now an open number.***
[MEASURED: CHORUS §3(3) / RECUT §2.3]

**S7.6.6 THE PRE-REGISTERED STOP FIRES — with the IGNITE-2 anatomy intact.** 2 of 10 scrambled-source
realisations certify at some iteration; **1 keeps** its certification to the fixed point (J1640+2224,
dlnL = 12.14 vs an 8.51 floor, qmax = 1.000, **Δk_oracle = −266** — a confident noise-lock under a scrambled
comb, **present at iteration 0 and never touched by the loop**). The other **SELF-CLEANS** — the M-step and
re-scored E-step DROP it, reproducing IGNITE-2's behaviour on the mixed problem. **No scrambled trajectory
grows; wrong-cert counts are flat in all 30/30 loop trajectories.** ***Every STOP event is inherited from the
criterion's seeds; none is generated by the loop.*** **The D1 hole travels unchanged, in its mixed-model
flavour.**
[MEASURED: CHORUS §2(c)]

### S7.7 The handoff that is not ours

**S7.7.1** Whether real SMBHB environmental evolution produces a **stochastic `df/f` above the P2-B
thresholds** (~1e-6…1e-5 over kyr lags, S2.2.6) is an **astrophysics** question (gas/stellar coupling) that
this project's methods **cannot answer**. **HANDOFF (Taylor/Farr), formulated and assigned outward, still
open.**
[MEASURED: project_progress PRONG-2 CLOSURE]

---

## S8 — THE LOOP

### S8.1 The hard-lock cascade

**S8.1.1 B1.3 (below onset): VERDICT (iii) DAMPED.** 12 Arm-B realisations with the source-fit channel wired
in (certified fringes fixed → re-fit (f, mc) on noisy data → shrink σ(log10_mc) → open the registration gate
→ re-certify). **0/12 realisations grew past their round-0 seed set**; median N_cert **1 → 1**; pooled
next-cycle gain **0.00**; within-loop σ(log10_mc) shrink **1.00× (saturates immediately)**. Two measured
mechanisms: (1) **local σ_mc is not the bottleneck** (the fit reaches 1.4e-4 dex, far below the 0.003-dex
bound, but that tight *local* width never becomes *global fringe concentration*); (2) **the seed set is
wrong-fringe-poisoned** — the census-sky realisation seeds {J1603, J1713, J1909}, **all three on the WRONG
fringe**. ***The loop cannot hot-wire its own clock.***
[MEASURED: FORGE §2]
[SUPERSEDED → **hard-lock-only. B1.3 was ALSO measured BELOW ONSET** — the loop was fed a seed set that,
under a gate the null cannot pass, **does not exist** (S3.2.4). A loop started from zero genuine detections
cannot compound, and that is arithmetic, not physics.]

**S8.1.2 IGNITE (above onset): a FOURTH MODE — a WRONG-CERTIFICATION CASCADE.** None of the three
pre-registered outcomes occurred:

| cell (h, T, tier; `fN` floor) | grew past seed | seeds → finals | wrong at fixed point | scrambled loop |
|---|---|---|---|---|
| −13.25, 30, vlbi; 5.46 | 2/10 | Σ15 → Σ28 | **26 of 28 wrong** | **2/5 DETECT (all wrong)** |
| −13.00, 30, vlbi; 15.55 | 1/10 | Σ5 → Σ6 | 5 of 6 wrong | — |
| −12.75, 30, lit; 30.89 | **6/10, incl. 2 runaways (3 → 116)** | Σ20 → **Σ359** | **356 of 359 wrong** | — |

**By raw count the loop "compounds" spectacularly (3 → 116 in three iterations); essentially every
certification it adds is on a FALSE fringe.** The genuine (on-true) certified count **never grows in any of
the 30 realisations**. ***The fit does not merely fail to compound — it DESTROYS the correct registration it
was fed*** (17 of 20 certs at the *quiet* fixed points are wrong post-fit, vs an 8–15 % wrong-rate in the raw
seeds). Once the source is mis-fit, the loop **is** a wrong source meeting loud real data — **the exact
configuration whose noise-lock sets the `fALL` floor** — and σ(log10_mc) collapses to ~1e-5 dex, opening the
registration gate for everything. ***Tight local width + wrong global registration = confident nonsense.***
Both pre-registered STOPs fired.
[MEASURED: IGNITE §5]
[SUPERSEDED → **hard-lock-only.**]

**S8.1.3 THE MECHANISM IS THE HARD LOCK, AND THE SPEC PREDICTED IT.** Both B1.3's and IGNITE's M-step **pins
each certified pulsar at its MAP fringe CENTRE** — up to half a fringe off the true within-fringe offset —
and fits the source against that delta. That mis-pin biases the (f, mc, extrinsics) gradient; one damped
Newton step moves the source materially (`src_mc_off` up to **1.6 dex in a single step**); the re-E-step at
the moved source re-registers the fringes; the poison compounds. **This is a GPS wrong-fix failure at loop
level. AND IT IS NOT THE M-STEP THE SPEC SPECIFIES:** spec §3 is explicit that the source update is **soft
and fringe-posterior-weighted** — `Q(θ) = Σ_p Σ_n q_p(n)·lnL(θ, L_p(n))` — i.e. it **marginalises** the fringe
uncertainty instead of committing to a fringe. ***The spec's own soft-fix discipline predicted exactly this
failure, and the implementation did not follow it.***
[MEASURED: spec §3 vs IGNITE §5 caveat; project_progress §10.8.3 block]

### S8.2 The soft loop — safe, consolidating, self-cleaning

**S8.2.1 IGNITE-2 ran the spec-§3 loop above onset. THE PREDICTION HELD: no cascade in 40/40 trajectories.**
(30 signal + 10 scrambled, `B1Marg`, **nothing ever hard-locked**, all 116 pulsars fringe-marginalised at
every step.) The wrong-cert count **NEVER GROWS** (against the hard lock's 3 → 116, 356/359 wrong); the
false-fringe mass `W` is **flat to ±1** (hard lock: W = 115–116 from the first iteration); `src_mc_off` stays
**< 1e-4 dex** at truth (hard lock: up to 1.6 dex in one step). **The soft M-step at the true source has ≈ zero
gradient and the line search REFUSES to walk off the peak** — 24 of 30 signal realisations terminate in 2
iterations with every proposed step rejected.
[MEASURED: IGNITE2 §3.1/§3.2]

**S8.2.2 σ(log10 mc) = 1e-5–1e-4 dex WITH `dk = 0`** — tight local width **AND** correct global registration,
**the exact inverse of the hard lock's "confident nonsense".**
[MEASURED: IGNITE2 §3.2]

**S8.2.3 Where the soft loop moves, it moves TOWARD the truth, both ways.**
- **Growth (2/30, both genuine, both `dk = 0`):** at the lit cell J1327−0755 crosses the 38.9-nat floor at
  iteration 2 of a smooth +27-nat `lnL_marg` climb and stays; at the vlbi cell J1713+0747 rises 10.7 → 12.9 nat
  over 9 accepted steps (+35 nat). ***The soft fit polishes the source against the fringe mixture and lifts
  genuine fringe-breaking evidence over the floor.*** **Real — but +1-class, in 2 of 30, with NO compounding.**
- **Self-cleaning (3 of 6 detecting scrambled realisations):** the M-step takes one large step (0.2–0.35
  scaled, +2000–3000 nat — *the scrambled source is nowhere near an optimum of loud real data*) and the
  re-scored E-step **DROPS the false certification** (2→0, 1→0, 1→0). **The hard lock did the opposite: its
  first step CREATED 116 detections.**
- **The one lost true cert is a MARGIN EVENT, not a wrong fix:** J1909 held at `dlnL = 8.97` against the
  7.59 floor — **a 1.4-nat margin**; two polish steps later the gap sits below the floor; **`dk = 0`
  throughout** — the pulsar never leaves the true fringe. *The D2 lesson (margins inside calibration scatter
  carry no weight) applied to a single object.*
[MEASURED: IGNITE2 §3.2]

**S8.2.4 FORMAL SUPERSESSION.** ***B1.3's DAMPED and IGNITE's CASCADE are HARD-LOCK-ONLY verdicts,
superseded-with-trail. The soft loop is spec §3's REFERENCE IMPLEMENTATION; the hard lock is retired.
Neither DAMPED nor CASCADE is a verdict on Hogg's iterated phase-up per se.***
[MEASURED: project_progress §10.8.3 / §10.10(c) / §10.12 item 1]

**S8.2.5 THE LOOP IS NOT THE PROBLEM. THE CRITERION IS.**
[MEASURED: IGNITE2 §3.3 — *"Hogg's iterated phase-up, implemented as the spec wrote it, is dynamically safe
and mildly consolidating — and it cannot rescue the criterion it is fed by."*]

### S8.3 The seeds problem

**S8.3.1 The pre-registered STOPs still fire — and EVERY failure is INHERITED FROM THE SEEDS; NONE is
generated by the loop.** At the lit cell, 4/15 signal realisations carry wrong certifications — J1545−4550
(Δk = −115, twice), J1525−5545 (−129), J1045−4509 (+133), J1125−6014 (−114) — **all at ITERATION 0**, i.e.
stage-1 impurity inherited as seeds, which the loop **neither amplifies nor removes** (their `dlnL` never
changes, because those realisations never accept a step). The 3 scrambled keepers are **small-|Δk|
noise-locks** (J0437 at Δk = −4, B1937+21 at +21, J0711 at +231) — **exactly D1's stated wrong-counterpart
hole, exactly the events the purity layer was meant to kill and, per S5.4, cannot.**
[MEASURED: IGNITE2 §3.2]

**S8.3.2 THE FRONTIER STATEMENT FOR THE LOOP: *the loop works given seeds. The modelled box supplies none.***
[MEASURED: project_progress §10.10(d)]

### S8.4 What is still untested about the loop

**S8.4.1** **The signal loops START AT THE TRUE SOURCE** (the conditional convention, as B1.3 and IGNITE's
hard loop did). *"Stable at truth where the hard lock cascaded from truth"* is the designed comparison —
**behaviour from a MIS-REGISTERED but unscrambled start is untested.**
[MEASURED: IGNITE2 §4.1 caveats]
[PENDING: ROBUST — the `tol > 0` trajectory through the **soft** loop. **Deliverables:** (i) does the soft
loop *recover* from a mis-registered start, *stay* mis-registered, or *walk further off*? (ii) the soft loop's
behaviour on the tolerance grid (tol = 0/1/2/5), scored under **per-tol refit floors** (S5.3.6); (iii) **D4
under a MIS-POSITIONED (rather than scrambled) counterpart** — the `tol > 0` axis of the discordance gate,
explicitly untested and open (D4 §7); (iv) whether the soft loop's *self-cleaning* property (S8.2.3) extends
to small counterpart offsets or is specific to gross scrambles.]

**S8.4.2** The M-step reuses the iteration-0 Hessian with a 3-candidate backtracking line search; steps are
accepted **only when `lnL_marg` + prior improves.** A "no-step" iteration is therefore a **measured property
of the objective**, not a solver failure — *but a richer optimiser could in principle move where this one
stayed* (it would have to fight an `lnL_marg` **decrease** to do so).
[MEASURED: IGNITE2 §4.1]

**S8.4.3** `σ(log10 mc)` as quoted by the loop is the **profile half-width of the M-step's own marginalised
objective**, **NOT a posterior credible width.** It must never be quoted as one.
[MEASURED: IGNITE2 §4.1]

**S8.4.4** Above-onset loop behaviour has been measured at **two properly-calibrated cells and nowhere else**
(S6.2.1).
[PENDING: SURFACE — the loop re-run at every cell that survives proper floor sizing. If the box closes
entirely (S6.3.1), the loop's above-onset behaviour becomes a statement about a regime nature does not supply,
and should be reported as such.]

---

## S9 — THE PAYOFF

### S9.1 The SIREN table — what certified pulsar terms are FOR

**S9.1.1 THE PAYOFF SENTENCE.** *"Conditional on `N_seed` certified pulsar terms, a single loud circular
SMBHB (per-source Earth-term SNR ≈ 33–54) is localised in luminosity distance to **σ(D_L)/D_L ≈ 6–12 %** for
`N_seed` = 3–5, against **332 %** from the Earth term alone, because the kyr-baseline pulsar terms measure the
chirp mass (σ(log10 Mc): **0.866 dex → 7e-4–0.03 dex**) and thereby break the chirp-mass/distance degeneracy.
**This is the same 10–30 % fractional-distance class that dark-siren H₀ programmes already treat as
cosmologically useful** — reached, in the nanohertz band, by three certified pulsar terms rather than by an
electromagnetic counterpart."*
[MEASURED: SIREN §4.2 — 189 exact Fishers; 9 cells (`log10_fgw` × `log10_mc`) × 7 seed sets × 3 arms]

**S9.1.2 The floor is the AMPLITUDE, not the chirp mass, and it is exact.** Once Mc is measured,
`σ(D_L)/D_L → ln10 · σ(log10 h)` — reproduced to **0.1–4 %** in every cell where Mc is measured. The crossover
is `σ(log10 Mc) < (3/5)·σ(log10 h) ≈ 0.02 dex`. ***σ(D_L)/D_L saturates; σ(log10 Mc) does not.*** **Beyond
`N_seed` ≈ 3–5 the siren improves only with SNR** — and `σ(log10 h)` is itself inflated **~3.5×** over
`1/(ln10·SNR)` by the **h–cos_inc–ψ degeneracy**. ***That degeneracy, not the chirp mass, is what a future
siren programme must attack.***
[MEASURED: SIREN §4.1/§8.4]

**S9.1.3 The sky collapses by 3–4 orders of magnitude** with seeds: ~1.6–4.3 deg² (Earth term) →
**~5e-4 deg² (≈ 1.7 arcmin², an ~80″ box)** with five seeds. **This IS the P1 needle seen from the amplitude
side**, and its scale agrees with B0.2's independently measured certification tolerance (~0.006° ≈ 22″).
`σ(log10 h)`, by contrast, is **flat in `N_seed`** — the amplitude never learns anything from the pulsar terms.
[MEASURED: SIREN §5]

### S9.2 Lag diversity, not count

**S9.2.1 "HOW MANY PULSARS?" IS THE WRONG QUESTION.** Adding a seed **at a lag you already have buys
nothing**: N1 → N2 adds J0437 (τ = 0.55 kyr) next to J1909 (0.69) — gain **1.02–1.07×, essentially nothing**.
Adding a seed at a **different** lag buys a lot, **longer OR shorter**: N3 → N5 adds the two *shortest* lags
(τ ≈ 0.22 kyr) and improves σ(log10 Mc) by **5.8×** while adding **+0.07 %** of mc lever. ***The gain is NOT
lever — it is breaking the Mc ↔ f_gw degeneracy.*** Because `∂Δφ_p/∂log10 f ∝ τ_p` and
`∂Δφ_p/∂log10 mc ∝ ḟ τ_p²`, their **ratio ∝ τ_p**: short-lag pulsars are near-pure **frequency** probes
(J1744: `g_f/g_mc` = 2400), and **pinning `f_gw` with them FREES the long lags to measure Mc.**
***Three well-chosen seeds beat five badly-chosen ones. Two seeds at the same lag are one seed.***
[MEASURED: SIREN §3.3/§6 — exact-phase autodiff of the phase `make_phase_connected_binary` actually builds]

**S9.2.2 Lag length dominates seed count.** **Three long-lag seeds beat FIVE frequency-ranked seeds in every
one of the nine cells** (by 1.005–1.50×), and beat the frequency-ranked *triple* by **5.8–7.9×**. Three
short-lag seeds are **3.3–4.0× WORSE** than the frequency-ranked triple. (Quadrature mc levers at
f = 10⁻⁸, Mc = 10⁹: shortest-3 = 1.891; freq-ranked-3 = 48.44; longest-3 = **2376** — a **49×** lever ratio
buying 7.7× in σ, sub-linear **because the long-lag set has a small lag SPREAD (0.128 dex): lever and
diversity are different resources and both are needed.**)
[MEASURED: SIREN §6]

### S9.3 The certification/siren tension — the campaign's structural finding

**S9.3.1 CERTIFIABILITY ∝ 1/τ; PAYOFF ∝ τ².** The registration tolerance is `tol ≈ 1/(2π f τ_p)`, so the
design theorem's *"wide lanes from nearby pulsars"* names the **NEAREST** pulsars — hence the **SHORTEST**
lags — while the chirp-mass lever goes as **τ²**. ***THE PULSARS THAT ARE EASIEST TO CERTIFY ARE THE WORST
CHIRP-MASS MEASURERS, AND THE RATIO GOES AS τ³.*** J0711 (τ = 0.220 kyr) has the array's loosest lane and an
mc lever of **0.408 rad/dex**; B1937+21 (τ = 7.768 kyr) has a lane **35× tighter** and a lever **4100×
larger**.
[MEASURED: SIREN §6.1]

**S9.3.2** ***The certification target list and the siren target list are DIFFERENT LISTS, and they
anti-correlate as τ³.*** The design-theorem target list **optimises the wrong objective** if the goal is the
standard siren. **The right objective is a LAG-DIVERSE set: short lags to certify and to pin `f_gw`, long lags
to carry Mc.**
[MEASURED: SIREN §6.1/§8.3 — new, measured, and it should change how a target sub-array is designed]

**S9.3.3 Cross-report coupling worth keeping.** RING says **only sub-3-pc (VLBI-class) σ_d matters**; SIREN's
headline arm B *assumes* **0.1 pc** seed distances — i.e. exactly that regime. The 0.1 pc premise is
**load-bearing for short-lag seeds** (arm C — one whole fringe of within-fringe ignorance — degrades N5 by
**5.0×**) but **irrelevant (0.2 %) for the long-lag triple.** *Mechanism: a short-lag pulsar's Mc contribution
rides on the linear phase `2πfτ_p`, which is exactly degenerate with `L_p`; a long-lag pulsar's rides on the
τ² chirp, which is not.*
[MEASURED: SIREN §7]

**S9.3.4 THE CRAMÉR-RAO CAVEAT — verbatim, and it governs every σ in S9.1–S9.3:**
> **These are Cramér–Rao bounds on zero-noise Asimov data with the fringe integers GIVEN.** A noisy
> realisation scatters; a full posterior with free integers can only be **wider**. Frozen GP hyperparameters
> (≤ 9 %, D5) and a single source (no confusion penalty) are the standing optimisms. The sky is free, which is
> the one place SIREN is *conservative* relative to the B1 targeted scenario.

**Every σ SIREN quotes is a LOWER BOUND.**
[MEASURED: SIREN §1/§8.5]

**S9.3.5 AND SIREN IS A GIVEN-SEEDS FORECAST.** Achievability is Track B's question and Track B's answer is
the information-theoretic NO-GO (`f = 6.9e-7`). **Under criterion-v1, Arm B certifies 0.000 seeds per
realisation — so SIREN's `N_seed` = 3–5 columns currently price a resource the noisy pipeline CANNOT DELIVER
AT ALL.** ***The payoff is real; the road to it is not through cold-start certification.***
[MEASURED: project_progress §10.3]

### S9.4 The eccentric self-siren — the door that needs no seeds

**S9.4.1 THE SOURCE IS ITS OWN SIREN.** An eccentric source at `(f_orb = 10⁻⁸, Mc ≳ 10⁹, e ≳ 0.58)` reaches
**σ(D_L)/D_L ≈ 12–14 %** — **the dark-siren-useful class — FROM ITS OWN EARTH TERM, WITH ZERO CERTIFIED PULSAR
TERMS.** SIREN reached the same class only with **3 certified seeds** (which the census recurs in **0 of 40**
skies, and which criterion-v1 says the noisy pipeline delivers **0.000** of). ***Eccentricity substitutes the
counterpart's own clock for the missing certified pulsar terms.***
[MEASURED: ATLAS M4 — 12.1 % at (10⁻⁸, 10⁹·⁵), 14.2 % at (10⁻⁸, 10⁹·⁰), against SIREN's circular N0 of
136 %/316 %]

**S9.4.2** It clears the **relative** bound only, so a residual **factor ~2–5 in σ_mc** remains for the EOB
tier to close (S7.4.6).
[MEASURED: ATLAS M4]

**S9.4.3** ***This is the expectation-value road, and it is the only one in the programme that does not depend
on a certification the population may never supply.***
[MEASURED: project_progress §10.8.4 / §10.12 item 4]

[DISPUTED: ATLAS's own consistency flag, **carried, not hidden**: the npz `e` column for the M4 rows holds
the **κ ≥ 20 min-e** (0.516/0.526/0.501) while the markdown labels it `e* (>20×)` (0.59/0.58/0.66); **σ_mc was
evaluated at the npz values.** *Do not silently reconcile.* A future session must decide which threshold the
M4 payoff numbers are quoted at, and re-derive if it is the other one.]

### S9.5 Imaging

**S9.5.1** The programme has calibrated the Boyle–Pen deconfusion bound (S2.1.3) — the capacity law of a PTA
*as an imaging telescope* — and has never used it as one. **The confusion-capacity law, the array-boost law
(`N_knee ≈ 0.40 · N_psr · T · Δf`, linear through 200 pulsars), and Hogg's spherical-harmonic framing
(angular resolution `min(λ/D, λ/[cT])` → sky structure to degree ~10) are three views of the same object, and
nobody has joined them.**
[MEASURED: ABACUS; P2-C; project_progress §3 (Hogg's remark)]

[PENDING: unnamed IMAGING campaign. **Deliverables:** (i) state the PTA's imaging capacity for the *real*
array (heterogeneous noise, real sky coverage) rather than the toy — the ABACUS mapping is exact at the
coefficient level but rests on the toy's homogeneous white noise and schematic amplitude normalisation;
(ii) connect the deconfusion bound to the **fringe-identification** layer, which Boyle–Pen never saw and which
S2.4.3 shows binds first for the real likelihood — *what does "characterise 2N/7 sources per bin" mean when
none of their pulsar distances is certified?*; (iii) decide whether the imaging frame is a **paper section** or
a **separate paper**; (iv) if the former, it belongs beside S2.1 as the outward-facing statement of the
confusion physics.]

---

## S10 — HONEST INFERENCE

### S10.1 What is already known about prong 3

**S10.1.1 The referendum settles what a sampler could and could not find.** The honestly fringe-marginalised
posterior **does not concentrate at truth** (`f = 6.9e-7`, S4.1.10). ***Prong 3 does not unpark on the data
alone*** — not because samplers are bad, but because **the marginal evidence itself prefers the wrong-fringe
re-registered plateau over truth.** The transient J1713 = 0.977 hint is real **per-pulsar** and **overwhelmed
in the aggregate evidence.** **Only more data, or external sky knowledge, helps.**
[MEASURED: spec R4 — this is the single most important input to any future sampler design]

**S10.1.2 The existing sampler's known failures, unchanged since the project's start.** The 3-phase annealed
snap-sampler locks distances reasonably but **struggles to find CW parameters**: **no global CW move, no
tempering swap** (single-chain annealing cannot recover a missed mode); the distance move can only migrate
±0.6 dL around a fresh prior draw (**~1/100s chance of hitting the right fringe**).
[MEASURED: project_progress §5 prong 3 + §7 issues 3–4]

**S10.1.3 The pipeline IS calibrated per-claim on signal-present data — and that was never the problem.** The
reliability curve tracks (claimed `q_max` 0.51 → 0.96 vs realised true-fringe fraction); **BH-FDR@0.05** gives
realised true fraction **1.000** (12 claims pass). ***The failure was never miscalibration — it was that a
calibrated confidence with no detection statistic under it still fires on pure noise.*** **Calibration and
detection are different questions.**
[MEASURED: FORGE §4 / B1.1]

**S10.1.4 A grid will not do.** RING S-1: a coarse+zoom sky grid **refines the wrong peak** in 61 % of
realisations and **100 % at SNR 20**, in every scenario, **because a louder signal narrows the fringes the
grid must resolve.** *"Any future Q1 sky search at `fgw ≈ −8` needs a sampler, or a two-stage
earth-then-pulsar-term search."*
[MEASURED: RING §3.4/§7 S-1]

**S10.1.5 And two coverage questions are already waiting for one.** RING S-4: with a GWB, **even exact
distances undercover** (`inside90` 0.40–0.50 vs nominal 0.90), and physics (a correlated background shifts
which fringe the MAP selects) has **not been separated** from estimator error (the "posterior" there is a
*profile* slice, not a marginal). **The scenario-C tier3 coverage number must not be quoted until a sampler
has been run.**
[MEASURED: RING §5.2/§7 S-4]

**S10.1.6** The soft-EM M-step is **numerically CHAOTIC** (cond ~4e10): *for any future use, replace Adam with
a Newton/trust-region step, or profile the extrinsics.* The soft loop (S8.2) already does the former; a
sampler must not inherit the former mistake.
[MEASURED: spec STATUS — TRACK B CLOSED]

**S10.1.7** The correct object for a sampler is **NOT** the conditional-on-truth ceiling. *"The right quantity
is the JOINT source+distance posterior (marginalising source-parameter uncertainty into the fringe posterior),
not the conditional-on-truth ceiling."*
[MEASURED: spec B0.2 IMPLICATION — the pivot that produced the whole joint-registration reformulation]

[PENDING: SAMPLER — the honest posterior machinery. **Deliverables:** (i) a sampler over the **joint**
(source × fringe-integer × within-fringe distance) posterior that does **not** condition on truth anywhere —
the star topology (spec §1) makes the fringe latents exactly factorisable given θ, which is the structure to
exploit; (ii) **SBC (simulation-based calibration)** on the fringe-integer posterior — the one calibration
test the programme has never run, and the only one that can validate `q_max` as a probability rather than as a
score; (iii) coverage of the **credible intervals**, not just the point verdicts, so that RING S-4 can be
closed; (iv) a global CW proposal + parallel tempering, which the existing sampler lacks; (v) a statement of
what a sampler *can* deliver given `f = 6.9e-7` — **the honest expectation is that it correctly reports a
diffuse fringe posterior, and a correct diffuse answer is a result, not a failure.**]

### S10.2 The realism ladder

**S10.2.1 What has been paid.** Real heterogeneous white noise + per-pulsar RN GP + HD-correlated GWB GP, all
**drawn** (not Asimov) in FORGE/IGNITE/IGNITE-2. **Noise is NOT a perturbation:** median per-pulsar CW rms
**458 ns** vs drawn noise **2005 ns** (white 1414, RN 805).
[MEASURED: project_progress B1 machinery gates, 8/8 PASS]

**S10.2.2 What has NOT been paid, in the order it matters:**
1. **Real data.** Nothing in the certification chain has touched a real TOA (S6.4.1). → [PENDING: ANCHOR]
2. **Mis-registered starts** through the soft loop (S8.4.1). → [PENDING: ROBUST]
3. **Mis-centred priors.** RING's harness can only make priors *wide*, never *mis-centred*; the canonical
   `best_distances.txt` means differ from the feather means by **1.40σ on J0437** (a **0.55 rad** pulsar-term
   phase error), and **the machinery supports this experiment with one dict swap** (S7.2.6, S-5). → [PENDING: ROBUST]
4. **Un-frozen GP hyperparameters.** Bounded at ≤ 9 % for the *amplitude* at ±0.5 dex (S2.4.4), but γ and the
   per-pulsar RN have never been marginalised. → [PENDING: ROBUST]
5. **The `fgw ≲ −9` regime**, where the RING harness is **uncalibrated** by a factor 6.42 in noise power
   (S7.2.6, S-3). → [PENDING: ROBUST]

[PENDING: ROBUST — the realism ladder. **Deliverables:** the five items above, each with a stated
*before → after* effect on the criterion-v2.1 count and purity, so that the programme can say which realism
gaps are ≤ 9 %-class (ignorable, like the GWB freeze) and which are verdict-moving (like the luminosity
function, and like Arm A → Arm B).]

### S10.3 End-to-end

**S10.3.1** No single artefact in this repo runs the whole chain. The campaigns are **fenced experiments on
banked machinery**, each answering one question with its own gates. **The chain — real data → noise model →
detection statistic → fringe posterior → certification → source fit → payoff — has never been executed as one
pipeline, on one dataset, with one set of conventions.**
[MEASURED: structural, from the file inventory: `stagec_*`, `trackB_*`, `hpc_harbor/{geo,ring,siren,atlas,
forge,ignite,ignite2}/` are seven separate drivers over shared primitives]

[PENDING: PIPELINE — the end-to-end campaign. **Deliverables:** (i) ONE driver that runs the full chain and
reproduces, as gates, the banked value anchors — **`lnL(truth | zero-noise) = 405413.512739`** (bit-equal
across H4c/R/GEO/SIREN/IGNITE), the census triple **(0.9534 / 0.9887 / 0.9984)**, and the criterion-v1 final
table (GEO 0.275 / Arm A 0.067 / Arm B 0.000 / nulls 0.000); (ii) one conventions file, so that the tier, the
null family, the floor estimator, the arm, and the tolerance are **stated once and inherited**, not re-decided
per campaign; (iii) a **provenance record per number** — which npz, which seed, which floor, which N — so that
a paper table is generated, never transcribed; (iv) the **cost** of the full chain, so that a reviewer's
"re-run it with X changed" is answerable in GPU-hours rather than in months.]

---

## S11 — THE FORECAST

### S11.1 The join: onset × ATLAS corner × SCOUT clock

**S11.1.1 THE FIRST SOURCE HAS FOUR REQUIREMENTS, AND THEY WERE FOUND ONE CAMPAIGN AT A TIME.** ATLAS: the
first source must be **eccentric (e ≳ 0.6), massive (Mc ≳ 10⁹), at the TOP of the band (f_orb ≳ 10⁻⁸)** to
self-clock. **IGNITE adds the fourth: NEAR.**
[MEASURED: IGNITE §6 — *the self-clocking corner and the certification-onset corner are THE SAME CORNER, plus
a distance cut*]

**S11.1.2 The join table** — maximum distance at which the first source can sit and still clear the (now
retracted) onset:

| f_orb (Hz) | log₁₀Mc | ATLAS min-e | `D_L ≤` at h\* = −12.75 (30 yr, lit) | `D_L ≤` at h\* = −13.25 (30 yr, VLBI) |
|---|---|---|---|---|
| 10⁻⁸ | 9.0 | 0.58 | **2.5 Mpc** | **7.8 Mpc** |
| 10⁻⁸ | 9.5 | 0.59 | **16.8 Mpc** | **53.1 Mpc** |
| 10⁻⁸ | 8.5 | 0.70 | 0.4 Mpc | 1.1 Mpc |
| 10⁻⁸·⁵ | 9.5 | 0.66 | 7.8 Mpc | 24.6 Mpc |
| 10⁻⁸·⁵ | 9.0 | 0.77 | 1.1 Mpc | 3.6 Mpc |
| 10⁻⁸·⁵ | 8.5 | 0.84 | 0.2 Mpc | 0.5 Mpc |

Only the `(Mc ≳ 10⁹·⁵, f_orb = 10⁻⁸)` corner clears onset **beyond the Virgo distance**, and only in the VLBI
tier at 30 yr; the reference `(10⁹, 10⁻⁸)` cell must sit **within ~8 Mpc** — an **M87-class host inside the
Virgo volume**.
[MEASURED: IGNITE §6 — strain convention verified against SIREN's banked D_L table]
[SUPERSEDED → **these D_L ceilings assume the onset that IGNITE-2 retracted (S6.1.2).** With no
properly-calibrated onset in the box, the requirement tightens from *"eccentric, massive, top-of-band, within
~8–53 Mpc"* to ***"louder or longer than anything in the modelled box."*** **The table is kept because its
STRUCTURE (the corner is the corner, plus a distance cut) survives; its NUMBERS are upper bounds on a bar that
moved up.**]
[SUPERSEDED AGAIN → **the onset EXISTS (S6.3.2), so the ceilings are recomputed — and they move in BOTH
directions at once.** This is the join as SURFACE re-measured it, against the census-structure onset
(h\* = −12.50) and against the 5+11 onset (h\* ≤ −13.25):]

| f_orb (Hz) | log₁₀Mc | ATLAS min-e | `D_L ≤` at **h\* = −12.50** (census 3+13, 30 yr, lit) | `D_L ≤` at **h\* = −13.25** (5+11, T ≥ 40) |
|---|---|---|---|---|
| 10⁻⁸ | 9.5 | 0.59 | **9.4 Mpc** | **53.1 Mpc** |
| 10⁻⁸ | 9.0 | 0.58 | 1.4 Mpc | 7.8 Mpc |
| 10⁻⁸ | 8.5 | 0.70 | 0.2 Mpc | 1.1 Mpc |
| 10⁻⁸·⁵ | 9.5 | 0.66 | 4.4 Mpc | 24.6 Mpc |
| 10⁻⁸·⁵ | 9.0 | 0.77 | 0.6 Mpc | 3.6 Mpc |
| 10⁻⁸·⁵ | 8.5 | 0.84 | 0.1 Mpc | 0.5 Mpc |

- **AGAINST the programme:** in the census structure the onset is **LOUDER** than IGNITE's retracted `h*`
  (−12.50 vs −12.75), so the reference `(10⁹·⁵, 10⁻⁸)` corner must sit inside **9.4 Mpc** — ***INSIDE the
  Virgo distance, not at 16.8 Mpc.*** SCOUT's population clock prices a source **1.25 dex quieter** than
  that. **On the census population, certification remains a variance play on the loud-nearby tail — and
  S7.1.1a caps the baseline lever at ~40 yr.**
- **FOR the programme, and this is the new thing:** ***the frontier is not a property of the SOURCE alone.***
  At 5+11 the onset reaches **h\* ≤ −13.25 (unbounded below by this grid)**, putting the same corner at
  **53 Mpc** — *outside Virgo, a **~180× larger volume**.* ***The single most valuable unknown in the
  forecast is therefore no longer "how loud is the loudest source" but "HOW MANY LOUD SOURCES ARE THERE" —
  and nobody has measured that number.***

[MEASURED: SURFACE §10 — strain convention verified against SIREN's banked DL table (Mc = 10⁹,
f_gw = 2×10⁻⁸ Hz at 7.77 Mpc → log₁₀h = −13.25). The h\* cells are both untouched by the floor fix
(zero-fraction 0.00), so these ceilings are re-cut-stable.]

**S11.1.3 THE POPULATION VERDICT.** Against SCOUT's clock (**N̄_detectable ≲ 0.01–0.1 current, O(0.1–1) SKA —
at −13.75-class loudness**), while IGNITE's onset is **0.5–1.0 dex LOUDER STILL** (and criterion-v2 moved that
bar **up**, not down):

> ***Nature does not supply a source above h\* at any epoch this campaign models
> (`N_joint(h > h*) ≪ N_joint(−13.75) ≲ 0.1`). The honest-certification programme is a VARIANCE PLAY on the
> loud-nearby tail, not an expectation-value plan.***

[MEASURED: IGNITE §6 / project_progress §10.8.4, strengthened by §10.10(d)]

**S11.1.4 THE EXPECTATION-VALUE ROAD.** ***It remains ATLAS M4's eccentric Earth-term standard siren, which
needs NO certified pulsar terms at all*** (S9.4).
[MEASURED: project_progress §10.12 item 4]

### S11.2 The frontier statement

**S11.2.1** ***The pulsar term is a kyr-baseline TIMESTAMP. It cannot be read without the CLOCK RATE. The
clock rate is `fdot`, i.e. `Mc`. A 22-yr Earth term cannot measure `fdot`.*** Everything else in this
programme is a consequence, a price, or an attempt to buy the clock somewhere else.
[MEASURED: spec B1 STEP 2 FRONTIER STATEMENT]

**S11.2.2** ***The loop works given seeds. The modelled box supplies none.*** (S8.3.2)

**S11.2.3** ***The certification bottleneck is the criterion's purity above onset, and it has NO geometric
fix.*** (S5.4.7)

### S11.3 The capabilities-to-science map

**S11.3.1** The programme's forward statement is **not** a prediction; it is a **map from capability to
science**. Each row is a lever whose price is measured and whose payoff is measured, and the honest posture is
that **no single row buys the certification programme by itself.**

| capability | what it buys | measured price / status | claim |
|---|---|---|---|
| **T → 40 yr** | the strongest lever — **and it CEILINGS here** | `T^{5/2}` beats the floor race up to ~40 yr; **T = 30 is optimal in 0 of 36 columns; loud cells LOSE at 50 (12/12)** | S7.1.1a, S6.3.2 |
| **T → 11 kyr** | the blind wall closes | **not an observational horizon** | S4.1.13 |
| **VLBI σ_d < 3 pc** | detections (ΣK 88 454 → 470) **AND** wrong-counterpart robustness (D4-OBS ×1.6) | a real VLBI campaign on the readable sub-array; **binary, not graded — Gaia is a NO-OP**; **and it COSTS +2.9 ± 1.0 nat of null floor at h = −13.25** | S7.2, S7.2.8 |
| **louder source (h)** | *almost nothing* — the floor rises as `h^1.72` | **the count is NON-MONOTONE in h**; and SCOUT says nature does not supply one | S5.3.2, S7.3 |
| **an EM counterpart (sky)** | removes the sky axis; **the break-even is 0.188°/source and a host supplies it** | free (**EM astrometry beats the tolerance by ≥10⁴×**) — **but it does NOT supply the clock** | S4.1.11, S4.2.7, S7.3.4 |
| **an EM period** | **NOTHING** (a 7× tighter `f` moves the best rung 3 %) | free, and worthless | S4.2.7 |
| **a host distance (1 %)** | σ(log10_mc) → 0.0205–0.0364 dex | still **>20× short**; the σ_h → 0 floor is 0.00301 dex | S4.2.11 |
| **eccentricity (e ≳ 0.6, Mc ≳ 10⁹, top-of-band)** | the **only** thing that moves the targeted evidence: **likelihood structure a prior box cannot supply**; and **σ(D_L)/D_L ≈ 12–14 % with ZERO seeds** | clears the *relative* >20× bound; **2.5× short of the absolute 0.003-dex floor**, and the corner that would clear it is **toy-invalid** | S4.2.12, S7.4, S9.4 |
| **the EOB tier** | the e ≳ 0.85 corner the toy cannot see | **still open** — *not a refinement, the load-bearing next step for the ABSOLUTE floor* | S7.4.6 |
| **ONE eccentric loud member (e ≳ 0.5)** | **the count switches ON at census loudness: 3.13 lit / 2.27 vlbi — the 3-pulsar siren threshold, exactly** | **MEASURED.** At e = 0.7 it reaches 5.4/5.8 — **14.8× / 12.4×** over circular. **But the clock is NOT shared: the gain is the comb's OWN pulsar terms** | S7.6.3, S7.6.4 |
| **TWO+ eccentric members (e ≳ 0.3)** | the same switch-on, at **one third the eccentricity** | **MEASURED, CONFIRMED both tiers** — *the threshold is a MIXTURE property, not a property of `e`* | S7.6.4 |
| **two more LOUD sources (3+13 → 5+11)** | **up to 6.1× the count (2.5× median), and the frontier moves ≥ 0.75 dex** — `h*` reaches the faint edge of the grid, putting the reference corner at **53 Mpc instead of 9.4** | **MEASURED.** Nobody has measured how many loud sources exist. ***This is now the single most valuable unknown in the forecast*** | S7.6.2, S11.1.2 |
| **3–5 certified pulsar terms** | **σ(D_L)/D_L = 6–12 %** — the dark-siren-useful class | the noisy pipeline currently delivers **0.000** of them; **the loop above onset has never been run on a low-impurity seed set** | S9.1, S9.3.5, PENDING: KINDLE |

**S11.3.2 The one-sentence forecast.**
[SUPERSEDED → *"…and the one premise that could overturn either statement — that certification is a property
of the POPULATION rather than of the source — has never been tested."* **It has now been tested, and the
answer is HALF YES AND HALF NO — and the half that is "no" is the half everyone expected to be "yes."**]

> ***Certification is a variance play on a loud, near, eccentric, top-of-band source that nature is not
> expected to supply; the eccentric Earth-term siren remains the expectation-value road and needs no
> certification at all; and the premise that could have overturned both — that certification is a property
> of the POPULATION rather than of the source — is TRUE, but not for the reason it was posed. The
> POPULATION's structure (how many loud members, how many eccentric ones) moves the count by 6× to 15× and
> moves the frontier by a factor of 180 in volume. The CLOCK, however, is NOT shared: an eccentric member
> certifies its own pulsar terms and drags the array's count up single-handedly. So the no-gos ARE rescoped
> — but by a population parameter nobody has measured, not by a mechanism anybody can exploit.***

[MEASURED: SURFACE §9.1/§10; CHORUS §1/§2 — both re-cut on the criterion-v2.2 floors]

---

## APPENDIX A — THE PENDING INVENTORY

**DELIVERED since criterion-v2.1** — struck here, kept for the trail:

| ~~tag~~ | delivered | where it now lives | what it did NOT supply |
|---|---|---|---|
| ~~**SURFACE**~~ | 108 cells, 24 840 realisations, ≈ 11 GPU-h | **S6.3.2–S6.3.5**, S7.1.1a, S7.6.2, S11.1.2 | **the `tol` axis was NOT swept** — it remains measured at one (h, T, tier) cell only → folded into **ROBUST** |
| ~~**CHORUS**~~ | 26 mixture cells, 4 000 nulls, 40 exact pairs, 30 loops | **S7.6.2–S7.6.6** | the trade-curve **crossing** (demoted, S7.6.5); mis-registered starts → **ROBUST** |
| ~~**ANCHOR**~~ | 7 200 realisations; realism ladder + the floor-estimator defect | **S6.4.3–S6.4.6**, **S6.5**, S7.2.8 | ***the REAL-DATA question — its premise was refuted at Task 1 (the repo has no real residuals). Re-tagged REAL-ARRAY.*** |

**OPEN:**

| tag | scope | first appears | what it must supply |
|---|---|---|---|
| **KINDLE** | the loop above onset, on the CORRECTED cells | S6.3.5, S9.3.5, S11.3.1 | **the gain contour on the corrected above-onset cells** — SURFACE reserved Pair B ((−12.75, 40, vlbi, 5+11) at 4.07/real and (−13.00, 40, vlbi, 5+11) at 3.57, **wrong-cert 0.07–0.13** — the first genuinely above-onset, low-impurity seed set the programme has ever had, and **both untouched by the floor fix**) and the soft loop was **never run on it**. **Two named questions must be answered BEFORE the contour is drawn, not after:** **(i) THE EXACTLY-1.000 MARGINALITY** — two cells ((−13.25, 40, vlbi, 3+13) and (−13.00, 40, vlbi, 3+13)) post a count of **precisely 1.000**: 30 correct certs in 30 realisations. **The count is quantised at 1/30, and the strict `> 1` bar lands exactly on a lattice point.** Both currently read *below onset* **on the strictness of an inequality.** ***Is that a measurement or a convention? It is a convention, it is undeclared, and it is doing work — a bar that can be decided by `>` versus `≥` is not yet a bar.*** **(ii) THE ISOLATED FAINT-EDGE ONSETS** (S6.3.5) — frontier, or fluctuation? |
| **REAL-ARRAY** | the real-data anchor *(re-tagged from ANCHOR)* | S6.4.2 | **a RE-DERIVATION OF THE PRIOR STACK on a real array, not a substitution of a residual vector.** On disk and verified: NG 15 yr (66 psr, 1705 frequencies), NG 20 yr (77 psr), MPTA DR3 (83 psr). Supplies: the prior stack re-derived; the floor re-fit against that array's own noise; **the programme's first certification numbers that touch a real TOA**; and **it lifts ANCHOR's single-frequency ceiling** (S6.4.6) |
| **QUILL** | first-principles criterion | S4.3.1, S5.6.5, S5.7 | derive `floor(h,T,tol)`; unify `ln K` ⊕ floor as one statistic; adjudicate D1 decision-theoretically; price the wrong-counterpart hole as an FDR term; explain or demolish the one-wall coincidence; state the criterion's ROC. **NEW: derive the floor's VALIDITY DOMAIN from first principles** — S6.5 bounds the D2 estimator empirically, at a zero-fraction of 20 %, and nothing explains why that is the number |
| **ROBUST** | the realism ladder | S8.4.1, S10.2 | soft loop from a mis-registered start; **the `tol` grid through the soft loop — SURFACE did not sweep it**; **D4 under a mis-positioned counterpart**; the mis-centred-prior arm (1.40σ on J0437); un-frozen GP γ + RN; the `fgw ≲ −9` recalibration |
| **SAMPLER** | honest posteriors + SBC | S1.3.3, S10.1 | joint (source × integer × distance) sampler with no conditioning on truth; **SBC on the fringe-integer posterior**; credible-interval coverage (closes RING S-4); global CW proposal + tempering |
| **PIPELINE** | end-to-end | S10.3 | one driver, one conventions file, provenance per number, cost of the full chain |
| **EOB tier** | the e ≳ 0.85 corner | S7.4.6, S11.3.1 | σ(log10_mc) where the toy comb coalesces — **the only place the absolute 0.003-dex floor could be cleared** |
| **IMAGING** *(unnamed)* | the PTA as an imaging telescope | S9.5 | the real-array capacity law; the deconfusion bound joined to the fringe-ID layer; paper-section vs separate-paper decision |
| *(ACCRE, CPU — a task, not a campaign)* | **the last piece of the floor fix** | S6.5.3 | **bank SURFACE's `fALL` offenders** and re-cut the 21-cell `fALL` ignition against the validity gate. The `fALL` zero-fractions were **never banked**, so the claim *"21 cells ignite under fALL"* **stands on the pre-fix estimator and cannot be checked from disk.** Proxy evidence says it is safe. **Proxy evidence is not a re-cut.** |
| *(handoff, not ours)* | environmental `df/f` | S2.2.6, S7.7 | whether real SMBHB environments produce stochastic `df/f` ≳ 1e-6–1e-5 over kyr lags — **Taylor/Farr** |

## APPENDIX B — THE DISPUTED INVENTORY

| # | dispute | where | how it gets resolved |
|---|---|---|---|
| D-1 | **ATLAS M4's `e` column**: npz holds the κ ≥ 20 min-e (0.516/0.526/0.501); the markdown labels it `e* (>20×)` (0.59/0.58/0.66); **σ_mc was evaluated at the npz values.** *Do not silently reconcile.* | S9.4 | decide which threshold the M4 payoff is quoted at; re-derive if it is the other one |
| D-2 | **Tier-C `f`**: frozen 4-seed protocol gives **0.0323 ± 0.0134**; the table's auto-ingest of the completed 5-seed npz gives **0.0431 ± 0.0185**. **Both FAIL identically.** | S4.2.8 | re-run the frozen protocol at a stated seed count and adopt one; the verdict does not depend on it |
| ~~D-3~~ | ~~**GEO's 0.275/draw** has never been re-cut under a properly-sized floor.~~ **CLOSED by SURFACE §8.** Applying a noisy-null floor to a noiseless statistic is a **CATEGORY ERROR**, not a mis-sizing. Zero-noise ceiling = **1.350 ± 0.82/draw** under the flat gate; **`0.275` RETIRED** (its sign was always safe; its value was never meaningful). And **GEO's count is IMPLEMENTABLE, not an oracle** (0 of 4640 cells have `dlnL < 0`, so `q_max ≡ P_true` at zero noise). | S3.2.2 | **RESOLVED 2026-07-12** |
| ~~D-4~~ | ~~**The four cells still posting > 1** are UNCALIBRATED.~~ **CLOSED by SURFACE §7.** N = 100 floors: **2 RETRACT, 1 CONFIRMS** ((−12.50, 30, lit) — the programme's first confirmed onset cell), **1 stays MARGINAL.** **All four verdicts survive the criterion-v2.2 floor fix.** And the mechanism is **not** "10-null floors are biased low" — at (−12.50, lit) the properly-sized floor came out **11 % LOWER** and the cell survived. *A max-of-N floor is an order statistic with no fixed false-alarm rate; it lands wherever its ten draws put it.* | S6.1.5 | **RESOLVED 2026-07-12/13** |
| D-5 | **The one-wall coincidence** (27× blind / >20× targeted / >20× ATLAS relative / κ ≥ 20 WEAVE): four *different objects* landing within ~1.5× of each other, with **no banked explanation**. ATLAS explicitly warns the last three are three distinct "20"s that earlier work conflated. | S4.3.1 | [PENDING: QUILL] item 4 |
| **D-6** | **The ONSET BAR IS AN UNDECLARED CONVENTION at the lattice point.** Two cells post a count of **exactly 1.000** (30 correct certs in 30 realisations). The count is **quantised at 1/30**; the test is **strict** (`> 1`); so both read *below onset* **on the strictness of an inequality that nobody ever declared.** `>` versus `≥` moves N_onset from 59 to 61. **Flagged, not resolved.** | S6.3.2, KINDLE(i) | [PENDING: KINDLE] item (i) — **declare the convention before the contour is drawn, not after** |
| **D-7** | **The `fALL` column was never re-cut** and its zero-fractions were never banked, so *"21 cells ignite under `fALL`"* (S6.1.3's refutation) **stands on the estimator that criterion-v2.2 bounds.** Proxy evidence (the matched `fN` zero-fractions of the 21 cells are 0.00–0.03) says it is **very probably safe.** ***Probably is not a re-cut.*** | S6.1.3, S6.5.3 | ACCRE CPU job — bank the `fALL` offenders; the banks already exist |
