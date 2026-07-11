# SIREN_payoff_chain.md — agent SIREN, ACCRE

Untracked working report. Cronus is canonical; this file never feeds back into
`project_progress.md`. No tracked file was edited, nothing was committed or pushed.

- Campaign: the **payoff Fisher** — Asimov, Fisher-level, no sampling.
- Grid: 1 loud circular source (`log10_h = -13.25`) × `log10_fgw ∈ {-8.5,-8.0,-7.7}`
  × `log10_mc ∈ {8.5,9.0,9.5}` × `N_seed ∈ {0,1,2,3,5}` + 2 lag-selected 3-sets × 3 arms
  = **9 cells × 7 seed sets × 3 arms = 189 exact Fishers.**
- Machinery: `trackB_b1_core` (`MaskedDelay`, `build_b1_amortised`), `trackB_estimator`
  (`phi_scale`/`phi_bounds`), `trackB_b1_pilot` (`dphi_exact`/`_grad_dphi`),
  `stagec_fisher`, `stagec_anchor_a2` — **imported and called verbatim.** Nothing
  reimplemented. HVP Hessian assembly is D1's production path.
- Jobs: gate/warm `12455052` (dgx03, 1h19m51s, exit 0, **ALL GATES PASS**);
  array `12456288_[0-2]` (3 tasks × 3 cells, all COMPLETED 0:0, 407–421 s).
- Artifacts: `SIREN_results/` — 9 lean cell npz + gate + geom + probe + summary + figure
  (**1.3 MB total**).

---

## 0. THE ANSWER IN ONE PARAGRAPH

**The payoff claim is correct, it is large, and it is not about how many pulsars you certify —
it is about the SPREAD of their lags.** With the Earth term alone the chirp mass is *not
measured at all*: in 8 of 9 cells the posterior on `log10 Mc` **is** the prior
(σ = 0.866 dex, information gain 1.00–1.05×). One certified pulsar term — J1909-3744, lag
τ_p = 0.69 kyr — takes σ(log10 Mc) to 0.33 dex at (f=10⁻⁸, Mc=10⁹) and to **0.013 dex** at
(f=10⁻⁷·⁷, Mc=10⁹·⁵); five take it to 0.025 and 7.3×10⁻⁴ dex. Propagated through the
amplitude relation, σ(D_L)/D_L falls from **332 %** (no seeds; the chirp-mass degeneracy,
unbroken) to **6–12 %** — squarely inside the 10–30 % class that dark-siren H₀ programmes
call useful. It then **stops falling**, because once σ(log10 Mc) < (3/5)·σ(log10 h) ≈ 0.02 dex
the distance error is entirely the *amplitude* error: the floor is exactly
`ln10 · σ(log10 h)`, reproduced to 0.1–4 % in every cell where Mc is measured. And the lag
structure dominates the seed count: **3 long-lag seeds (τ = 5.8–7.8 kyr) beat 5
frequency-ranked seeds** in every cell (by 1.0–1.5×) and beat the frequency-ranked *triple* by
**5.8–7.9×**, while the 3 shortest-lag seeds are **3.3–4× worse** than that same triple. The
mc lever arm goes as **τ²** and the competing frequency lever as **τ**, so short lags — the
nearest pulsars, which are exactly the design theorem's widest registration lanes — are the
*worst* chirp-mass measurers. **Certifiability ∝ 1/τ; payoff ∝ τ².** That tension is the
campaign's structural finding.

---

## 1. PROTOCOL, AND WHAT THESE NUMBERS ARE NOT

**GIVEN-SEEDS FORECAST.** Every number below is **conditional on N_seed certified pulsars.**
Whether that certification is achievable is a *separate, open question*, and the project's own
answer is currently negative: Track B's terminal verdict is that cold-start certification is an
**information-theoretic NO-GO** (f = Z_needle/(Z_needle+Z_plateau) = 6.9×10⁻⁷), and the census
ceilings are conditional-on-truth. Nothing in this report bears on achievability. This is a
*payoff* calculation: **if** you can certify N_seed pulsars, here is what you buy.

**Seed sets are SKY-CONDITIONAL.** Per `GEO_geometry_ensemble.md`, over 40 isotropic source-sky
draws J1909-3744 certifies in 95 % of skies, J0437-4715 in 80 %, J1713+0747 in 68 %,
J1744-1134 in 50 %, J0711-6830 in 40 %; the census's named triple recurs in **0 of 40**. I
therefore use GEO's **frequency-ranked** list as a nested sequence, never a fixed triple, and
report the lag-selected sets separately. A given sky realisation will hand you a different set;
read the N_seed columns as capability, not as a roster.

**Zero-noise (Asimov) data, so H = −Fisher exactly.** Data is injected at truth with the same
pulsar-term mask the model uses, so residual(truth) = 0 and the Hessian *is* minus the Fisher
(the D2 identity, and the STEP-1B loop-gain convention reused verbatim). These are therefore
**Cramér–Rao bounds**: the width of the central mode of the likelihood. A noisy realisation
scatters around them; a full posterior can only be *wider*.

**What "a seed" means, operationally.** `trackB_b1_core.MaskedDelay` puts a per-pulsar
pulsar-term switch in the likelihood as a runtime argument:
`delay_p = d_earth + m_p·(d_full − d_earth)`. A seed is `m_p = 1`. A non-seed is `m_p = 0` and
contributes only its Earth term — the Fisher-level statement of *"its fringe is unregistered,
so its pulsar-term phase carries no source information."* This is a **bracket, not an
identity**: a genuinely fringe-marginalised pulsar term retains a little fringe-averaged power.
The R POSTMORTEM (spec) shows that for **data-driven** pulsars — which is what the census-3 and
the GEO union are — the fringe information *evaporates* once joint registration fails
(P_true 0.12/0.014/0.019, exp H = 11–35 effective wrong fringes). So `m_p = 0` is the right
Fisher-level limit for an uncertified pulsar, and the **N_seed = 0 row IS the "marginalised"
arm** the brief asked for.

**Three arms, priced separately.** The brief's "certified = distance known to within-fringe,
~0.1 pc" is arm B, the headline.

| arm | seed pulsar terms | seed distance L_p | what it prices |
|---|---|---|---|
| **A** `exact-L` | ON | FIXED at truth | the ceiling |
| **B** `0.1 pc` | ON | free, N(L_true, 0.1 pc), Schur-marginalised | **the headline** |
| **C** `½-fringe` | ON | free, σ = dL_p/√12 (one whole fringe) | *sub-fringe* distance ignorance |
| — | OFF (`N_seed = 0`) | n/a | the marginalised / uncertified limit |

> **Arm C is NOT full fringe marginalisation.** It prices ignorance of where the pulsar sits
> *within* its fringe, with the integer still known. Full integer marginalisation is the
> N_seed=0 row (see above). Arm C is therefore an **upper bound** on what a partially-known
> distance supplies. I have labelled it `½-fringe`, not `fringe-marginalised`, for exactly
> this reason.

**The sky is FREE.** Unlike the B1 targeted scenario (sky supplied by an EM counterpart), every
Fisher here has `cos_gwtheta` and `gwphi` free with only the generative prior. SIREN's numbers
are therefore *conservative* relative to a counterpart-targeted analysis. All 8 source
parameters (sky ×2, `cos_inc`, `log10_mc`, `log10_fgw`, `log10_h`, `phase0`, `psi`) are free,
inverted **with the generative uniform priors and never with a pinv** — an axis the data does
not constrain has a zero Fisher eigenvalue, and pinv would silently report σ = 0 for it (the
exact `log10_mc` failure mode the STEP-1 pilot flagged).

**Inherited optimisms, unchanged.** GP hyperparameters (GWB, per-pulsar RN) frozen at truth
(D5 bounds the effect at ≤ 9 % on σ_L); linear timing model marginalised; single source, so no
source-confusion penalty.

**The fiducial source.** The census pop-draw's **loud0** (seed 3000): `cos_gwtheta = -0.61030`,
`gwphi = 4.23218`, `cos_inc = -0.60900`, `phase0 = 1.19323`, `psi = 1.50981`. Its drawn
parameters are **bit-identical** for `ncw=1` and `ncw=16` (gated, max|diff| = 0.000e+00), so
"the census's loud0 source, alone" is a well-defined object. Only `(log10_fgw, log10_mc)` are
swept and `log10_h` is pinned to the census loud value −13.25.

**Distances and lags are REAL.** `L_p = pe.pdist[0]` (the frozen EM prior means);
`τ_p = L_p (1 − cos μ_p)/c` under the loud0 sky.

---

## 2. GATES

| gate | result |
|---|---|
| **M3 reproduction** (ncw=16, mask=0, 24 loud params) vs the STEP-1 pilot | **PASS** (below) |
| Hessian symmetry residual, M3 | 4.83e-15 (pilot: 4.32e-15) |
| `ncw=1` vs `ncw=16` source-0 params | max|diff| **0.000e+00** |
| **FD curvature vs assembled Fisher** (`log10_mc`, `log10_h`, `log10_fgw`) | rel **3.7e-6 / 4.8e-6 / 8.8e-6** |
| arms A/B/C bit-identical at N_seed = 0 | **PASS** (0.8652942028931 ×3) |
| masked-off pulsar `∂lnL/∂L_p` (must be exactly 0) | **0.00e+00**, all cells |
| Hessian symmetry residual, campaign | max **1.01e-11**, and it occurs on the *least*-informative Hessian (f=−8.5, Mc=10⁸·⁵, Earth-only) |
| `(F+Π)` positive-definite; `‖A·A⁻¹ − I‖∞` | min eig 7.3e-6 > 0; residual ≤ 1.8e-11 |
| `lnL(truth | zero-noise)` | **405413.5127** (== H4c, == R, == GEO) |

### The N_seed = 0 gate, in full

The campaign's `N_seed=0` row must be the STEP-1 pilot's M3 Earth-term Fisher. Reproduced on
the pilot's own configuration (pop draw, N_CW=16, 3 loud sources, 24 free loud params,
inverted with the generative priors):

| quantity | SIREN loud0 / loud1 / loud2 | spec (`trackB_estimator_spec.md`) |
|---|---|---|
| σ(log10_fgw), scaled | 1.0076e-2 / 4.2012e-3 / 6.0919e-3 | 1.008e-2 / 4.201e-3 / 6.092e-3 |
| σ(log10_mc), scaled | 1.0020 / 1.7274 / 1.7310 | 1.002 / 1.727 / 1.731 |
| σ(log10_mc), dex | **0.5010 / 0.8637 / 0.8655** | **0.501 / 0.864 / 0.866** |
| σ(log10_h), dex | 0.0604 / 0.0359 / 0.0337 | 0.060 / 0.036 / 0.034 |
| prior σ, scaled | mc 1.7321, fgw 1.6358 | 1.732, 1.636 |

**Independent SNR cross-check** (not in the spec, added here). The conditional σ(log10 h) implies
per-loud-source Earth-term SNR = 1/(ln10·σ_cond) = **25.3 / 24.9 / 14.0**. Their quadrature is
38.2; adding the 13 faint sources (h = −14.25) in quadrature gives ≈ 39.2, against the directly
measured all-source Earth-term optimal SNR of **38.86**. And 25.3·√2 = 35.8 recovers P2-D's
banked full-model per-source SNR (~33.7 at h = −13.25) — the √2 being the pulsar term the
Earth-only fit does not see. The Fisher, the likelihood-based SNR, and the banked strain
reconciliation all agree.

---

## 3. THE HEADLINE TABLE — σ(log10 Mc) vs N_seed

Arm **B** (headline: seed distance known to 0.1 pc). Units: dex. Prior σ = 0.866 dex.
`N3S` = 3 shortest-lag, `N3L` = 3 longest-lag (union-18, under this sky).

| log10 f | log10 Mc | D_L (Mpc) | SNR | N=0 | N=1 | N=2 | N=3 | N=5 | N3S | N3L |
|---|---|---|---|---|---|---|---|---|---|---|
| −8.5 | 8.5 | 0.334 | 44.5 | 0.866 | 0.866 | 0.866 | 0.866 | 0.861 | 0.866 | 0.857 |
| −8.5 | 9.0 | 2.28 | 44.5 | 0.866 | 0.864 | 0.863 | 0.859 | 0.694 | 0.865 | 0.609 |
| −8.5 | 9.5 | 15.5 | 44.5 | 0.866 | 0.780 | 0.764 | 0.642 | 0.169 | 0.837 | 0.135 |
| −8.0 | 8.5 | 0.720 | 53.9 | 0.866 | 0.822 | 0.818 | 0.618 | 0.165 | 0.846 | 0.121 |
| −8.0 | 9.0 | 4.91 | 53.9 | 0.865 | 0.326 | 0.305 | 0.145 | **0.0251** | 0.480 | **0.0188** |
| −8.0 | 9.5 | 33.4 | 53.9 | 0.834 | 0.0633 | 0.0598 | 0.0216 | **3.87e-3** | 0.0803 | **3.19e-3** |
| −7.7 | 8.5 | 1.14 | 33.4 | 0.865 | 0.336 | 0.327 | 0.133 | 0.0252 | 0.488 | 0.0168 |
| −7.7 | 9.0 | 7.78 | 33.4 | 0.824 | 0.0611 | 0.0602 | 0.0198 | 3.86e-3 | 0.0780 | 2.82e-3 |
| −7.7 | 9.5 | 53.0 | 33.4 | **0.353** | 0.0131 | 0.0124 | 3.53e-3 | **7.28e-4** | 0.0129 | **6.12e-4** |

**Improvement ratio vs the N_seed = 0 row** (arm B):

| log10 f | log10 Mc | N=1 | N=2 | N=3 | N=5 | N3S | N3L |
|---|---|---|---|---|---|---|---|
| −8.5 | 8.5 | 1.00 | 1.00 | 1.00 | 1.01 | 1.00 | 1.01 |
| −8.5 | 9.0 | 1.00 | 1.00 | 1.01 | 1.25 | 1.00 | 1.42 |
| −8.5 | 9.5 | 1.11 | 1.13 | 1.35 | 5.13 | 1.04 | 6.42 |
| −8.0 | 8.5 | 1.05 | 1.06 | 1.40 | 5.26 | 1.02 | 7.15 |
| −8.0 | 9.0 | 2.66 | 2.84 | 5.98 | **34.5** | 1.80 | **46.1** |
| −8.0 | 9.5 | 13.2 | 13.9 | 38.7 | **215** | 10.4 | **262** |
| −7.7 | 8.5 | 2.57 | 2.65 | 6.51 | 34.3 | 1.77 | 51.5 |
| −7.7 | 9.0 | 13.5 | 13.7 | 41.7 | 214 | 10.6 | 292 |
| −7.7 | 9.5 | 27.0 | 28.6 | 99.9 | **485** | 27.4 | **577** |

### 3.1 The Earth term does not measure Mc at all

| log10 f | log10 Mc | Earth chirp phase π·ḟ·T² (rad) | σ(log10 Mc), N=0 | information gain over prior |
|---|---|---|---|---|
| −8.5 | 8.5 | 3.4e-5 | 0.8660 | 1.000 |
| −8.5 | 9.0 | 2.3e-4 | 0.8660 | 1.000 |
| −8.5 | 9.5 | 1.6e-3 | 0.8660 | 1.000 |
| −8.0 | 8.5 | 6.4e-4 | 0.8660 | 1.000 |
| −8.0 | 9.0 | 4.1e-3 | 0.8653 | 1.001 |
| −8.0 | 9.5 | 2.8e-2 | 0.8337 | 1.039 |
| −7.7 | 8.5 | 7.5e-3 | 0.8651 | 1.001 |
| −7.7 | 9.0 | 5.1e-2 | 0.8239 | 1.051 |
| −7.7 | 9.5 | **0.350** | **0.3530** | **2.453** |

In 8 of 9 cells the Earth-term posterior on `log10 Mc` **is** the prior. This is the STEP-1
pilot's `Δφ_E ≈ 0.05 rad` statement, now mapped: the Earth term begins to see the chirp only
when `π ḟ T²` approaches a radian, i.e. only for the most massive, highest-frequency source in
the grid.

> **On the field's "29 % within 1 dex" (arXiv:2502.16016).** I do *not* claim to reproduce that
> number, and our machinery cannot: the generative uniform prior on `log10 mc ∈ [7,10]` has
> σ = 0.866 dex, already inside 1 dex, so a "recovered within 1 dex" rate is a statement about
> the prior and the noise realisation, not about chirp information in the data. What our Fisher
> says is **sharper and strictly stronger**: for an Earth-term-era analysis of a source in
> most of this grid, the data contain *no chirp-mass information whatsoever* — the posterior is
> the prior to 0.1 %. That a 1-dex recovery then succeeds only ~29 % of the time is
> unsurprising; it is a prior-and-noise statistic. The 0.3-dex bar is the one that discriminates.

### 3.2 Where the bars fall (arm B)

First N_seed at which σ(log10 Mc) drops below each bar:

| log10 f | log10 Mc | < 1 dex | < 0.3 dex | < 0.1 dex | < 0.01 dex |
|---|---|---|---|---|---|
| −8.5 | 8.5 | 0 (prior) | never (N≤5) | never | never |
| −8.5 | 9.0 | 0 (prior) | never | never | never |
| −8.5 | 9.5 | 0 (prior) | 5 | never | never |
| −8.0 | 8.5 | 0 (prior) | 5 | never | never |
| −8.0 | 9.0 | 0 (prior) | 3 | 5 | never |
| −8.0 | 9.5 | 0 (prior) | 1 | 1 | 5 |
| −7.7 | 8.5 | 0 (prior) | 3 | 5 | never |
| −7.7 | 9.0 | 0 (prior) | 1 | 1 | 5 |
| −7.7 | 9.5 | 0 (prior) | 1 | 1 | 3 |

### 3.3 Saturation — it does not, and that is the wrong question

Marginal gain per added seed, σ ratio to the previous column (arm B):

| log10 f | log10 Mc | N0→N1 | N1→N2 | N2→N3 | N3→N5 |
|---|---|---|---|---|---|
| −8.0 | 9.0 | 2.66 | **1.07** | 2.11 | 5.77 |
| −8.0 | 9.5 | 13.2 | **1.06** | 2.77 | 5.57 |
| −7.7 | 9.0 | 13.5 | **1.02** | 3.05 | 5.12 |
| −7.7 | 9.5 | 27.0 | **1.06** | 3.50 | 4.86 |

σ(log10 Mc) has **not saturated by N_seed = 5** in any cell where the chirp is measurable — the
N3→N5 step still buys ~5×. But the gain per seed is wildly non-uniform, and the pattern names
the mechanism:

- **N1 → N2** adds J0437-4715 (τ = 0.55 kyr) to J1909-3744 (τ = 0.69 kyr) — a lag it already
  has. Gain: **1.02–1.07×. Essentially nothing.**
- **N2 → N3** adds J1713+0747 (τ = 1.22 kyr) — the longest so far. Gain 2.1–3.5×.
- **N3 → N5** adds J1744-1134 and J0711-6830 (τ ≈ 0.22 kyr) — the *shortest* so far. Gain
  4.9–5.8×.

Adding a seed at a lag you already have buys nothing. Adding a seed at a *different* lag buys a
lot — **whether longer or shorter.** The quantity that matters is **lag diversity.**

**The geometry proves it** (`siren_geom.npz`, exact-phase autodiff of the phase
`discovery.deterministic.make_phase_connected_binary` actually builds). The quadrature mc lever
`√Σ_p (∂Δφ_p/∂log10 mc)²` at (f = 10⁻⁸, Mc = 10⁹):

| set | τ range (kyr) | log10 τ spread (dex) | mc lever (rad/dex) |
|---|---|---|---|
| N3_freqrank | 0.55 – 1.22 | 0.346 | 48.444 |
| **N5_freqrank** | **0.22 – 1.22** | **0.743** | **48.480** |

The N3→N5 step adds **+0.07 %** of mc lever and **+0.40 dex** of lag spread — and improves
σ(log10 Mc) by **5.8×**. The gain is *not* lever. It is the breaking of the Mc ↔ f_gw
degeneracy: since `∂Δφ_p/∂log10 f ∝ τ_p` and `∂Δφ_p/∂log10 mc ∝ ḟ τ_p²`, their ratio is ∝ τ_p.
A short-lag pulsar is a nearly **pure frequency probe** (J1744: g_f/g_mc = 2400) while a
long-lag one is chirp-dominated (J1713: 120; B1937+21: 19). Pinning f_gw with the short lags
frees the long lags to measure Mc. **This is literally "Mc from the frequency ladder across
lags" — the payoff claim, measured.**

---

## 4. THE SIREN ROW — σ(D_L)/D_L

`log10 D_L = (5/3) log10 Mc + (2/3) log10 f_gw − log10 h + const`, which is *discovery's own*
amplitude relation (`dist = 2 mc^{5/3} w0^{2/3} / 10**log10_h`, `deterministic.py`, log10_h
branch). σ is propagated through the **full 8-parameter source covariance**, so the Mc–h
correlation is included, not assumed away. `σ(D_L)/D_L = ln10 · σ(log10 D_L)`. Arm B.

| log10 f | log10 Mc | D_L (Mpc) | N=0 | N=1 | N=2 | N=3 | N=5 | N3L |
|---|---|---|---|---|---|---|---|---|
| −8.5 | 8.5 | 0.334 | 3.32 | 3.32 | 3.32 | 3.32 | 3.30 | 3.29 |
| −8.5 | 9.0 | 2.28 | 3.32 | 3.32 | 3.32 | 3.30 | 2.66 | 2.35 |
| −8.5 | 9.5 | 15.5 | 3.32 | 3.00 | 2.93 | 2.47 | 0.648 | 0.538 |
| −8.0 | 8.5 | 0.720 | 3.32 | 3.16 | 3.14 | 2.37 | 0.635 | 0.468 |
| −8.0 | 9.0 | 4.91 | 3.32 | 1.26 | 1.18 | 0.559 | **0.118** | **0.0958** |
| −8.0 | 9.5 | 33.4 | 3.20 | 0.250 | 0.239 | **0.0988** | **0.0631** | **0.0657** |
| −7.7 | 8.5 | 1.14 | 3.32 | 1.30 | 1.27 | 0.523 | 0.141 | 0.124 |
| −7.7 | 9.0 | 7.78 | 3.16 | 0.262 | 0.260 | **0.124** | **0.106** | **0.104** |
| −7.7 | 9.5 | 53.0 | 1.36 | **0.114** | **0.113** | **0.102** | **0.0999** | **0.103** |

### 4.1 The floor is the amplitude, and it is exact

Once Mc is measured, `σ(D_L)/D_L → ln10 · σ(log10 h)`:

| log10 f | log10 Mc | σ(D_L)/D_L at N=5 | ln10 · σ(log10 h) | ratio |
|---|---|---|---|---|
| −8.0 | 9.5 | 0.0631 | 0.0606 | **1.041** |
| −7.7 | 9.0 | 0.1063 | 0.1045 | **1.017** |
| −7.7 | 9.5 | 0.0999 | 0.0998 | **1.001** |
| −8.0 | 9.0 | 0.1184 | 0.0655 | 1.81 (Mc not yet saturated) |
| −8.5 | 8.5 | 3.3042 | 0.0750 | 44.1 (Mc = prior) |

The crossover condition is `σ(log10 Mc) < (3/5)·σ(log10 h) ≈ 0.02 dex` — at which point the
`(5/3)σ_Mc` term drops below the `σ_h` term in quadrature. **σ(D_L)/D_L saturates, σ(log10 Mc)
does not.** So the operational saturation point for a standard siren is where Mc becomes
*sufficient*, not where it stops improving:

- (f = 10⁻⁷·⁷, Mc = 10⁹·⁵): **N_seed = 1 already saturates** (0.114 → 0.0999 from N1 to N5).
- (f = 10⁻⁷·⁷, Mc = 10⁹·⁰) and (f = 10⁻⁸, Mc = 10⁹·⁵): **N_seed = 3**.
- (f = 10⁻⁸, Mc = 10⁹·⁰): N_seed = 5, marginally.
- Low-f / low-Mc quadrant: **five seeds are not enough.**

**σ(log10 h) itself barely moves with N_seed** (0.0286 → 0.0263 at f = 10⁻⁸): the amplitude is
SNR-limited and is inflated ~3.5× over `1/(ln10·SNR)` by the h–cos_inc–ψ degeneracy. Improving
the siren beyond ~6 % therefore requires *more SNR*, not more seeds.

### 4.2 The sentence

> Conditional on N_seed certified pulsar terms, a single loud circular SMBHB
> (log10 h = −13.25, i.e. per-source Earth-term SNR ≈ 33–54 on this 116-pulsar array) is
> localised in luminosity distance to **σ(D_L)/D_L ≈ 6–12 %** for N_seed = 3–5, against
> **332 %** from the Earth term alone, because the kyr-baseline pulsar terms measure the chirp
> mass (σ(log10 Mc): 0.866 dex → 7×10⁻⁴–0.03 dex) and thereby break the chirp-mass/distance
> degeneracy. **This is the same 10–30 % fractional-distance class that dark-siren H₀
> programmes already treat as cosmologically useful** — reached, in the nanohertz band, by
> three certified pulsar terms rather than by an electromagnetic counterpart.

---

## 5. SKY AREA, AND σ(log10 h)

90 % credible sky area (deg²), arm B. Area = π·(−2 ln 0.1)·√det Cov[cos θ, φ]; the
(cos θ, φ) measure *is* solid angle.

| log10 f | log10 Mc | N=0 | N=1 | N=2 | N=3 | N=5 | σ(log10 h): N0 → N5 |
|---|---|---|---|---|---|---|---|
| −8.5 | 9.0 | 1.641 | 0.514 | 5.42e-2 | 1.91e-3 | 1.18e-3 | 0.0331 → 0.0325 |
| −8.0 | 9.0 | 1.633 | 0.440 | 4.54e-2 | 4.25e-3 | 5.18e-4 | 0.0286 → 0.0285 |
| −7.7 | 9.5 | 4.333 | 0.837 | 5.46e-2 | 2.59e-3 | 4.37e-4 | 0.0458 → 0.0433 |

The sky collapses by **3–4 orders of magnitude**, from ~1.6–4.3 deg² (Earth term) to
~5×10⁻⁴ deg² (≈ 1.7 arcmin², i.e. a ~80″ box) with five seeds. This *is* the P1 needle seen
from the amplitude side, and its scale agrees with B0.2's independently measured certification
tolerance (~0.006° ≈ 22″). σ(log10 h), by contrast, is flat in N_seed — the amplitude never
learns anything from the pulsar terms.

---

## 6. LAG LEVERAGE — the campaign's structural finding

Three seeds, chosen three ways (arm B, σ(log10 Mc), dex):

| log10 f | log10 Mc | N3 freq-ranked | N3 shortest-lag | N3 longest-lag | short/rank | **long/rank** |
|---|---|---|---|---|---|---|
| −8.5 | 9.5 | 0.642 | 0.837 | 0.135 | 0.77 | **4.76** |
| −8.0 | 8.5 | 0.618 | 0.846 | 0.121 | 0.73 | **5.10** |
| −8.0 | 9.0 | 0.145 | 0.480 | 0.0188 | 0.30 | **7.71** |
| −8.0 | 9.5 | 0.0216 | 0.0803 | 3.19e-3 | 0.27 | **6.77** |
| −7.7 | 8.5 | 0.133 | 0.488 | 0.0168 | 0.27 | **7.91** |
| −7.7 | 9.0 | 0.0198 | 0.0780 | 2.82e-3 | 0.25 | **6.99** |
| −7.7 | 9.5 | 3.53e-3 | 0.0129 | 6.12e-4 | 0.27 | **5.77** |

Seed sets under the loud0 sky:

- **freq-ranked-3** = J1909-3744 (0.692 kyr), J0437-4715 (0.549), J1713+0747 (1.218)
- **shortest-3** = J1545-4550 (0.098), J0711-6830 (0.220), J1744-1134 (0.225)
- **longest-3** = J0613-0200 (5.784), J2317+1439 (7.045), B1937+21 (7.768)

**Lag length dominates seed count.** Three long-lag seeds beat *five* frequency-ranked seeds in
**every one of the nine cells** (by 1.005–1.50×), and beat the frequency-ranked triple by
**5.8–7.9×** wherever the chirp is measurable. Three short-lag seeds are **3.3–4.0× worse** than
the frequency-ranked triple.

**Why**, measured (`siren_geom.npz`): the mc lever goes as **τ²** (Δφ_chirp = π ḟ τ²) and the
frequency lever as **τ** (Δφ_lin = 2π f τ). At (f = 10⁻⁸, Mc = 10⁹) the quadrature mc levers are

| set | lever (rad/dex) | ratio to freq-ranked-3 |
|---|---|---|
| N3 shortest-lag | 1.891 | **0.039** |
| N3 freq-ranked | 48.44 | 1 |
| N3 longest-lag | 2376 | **49.0** |

(σ improves sub-linearly in the lever — 49× lever buys 7.7× σ — because the long-lag set has a
*small* lag spread, 0.128 dex, so it separates Mc from f_gw poorly. Lever and diversity are
different resources and both are needed.)

### 6.1 The (T/τ)² tension, resolved: short lags LOSE

The brief anticipated that "nearby pulsars appear a **fourth** time if short lags win." **They
do not.** Short lags lose, decisively and in every cell. And this is not a null result — it is a
direct tension with the design theorem:

- **Certifiability.** The 1-radian registration tolerance is `tol ≈ 1/(2π f τ_p) ∝ 1/τ_p`. The
  design theorem's **lever (i)** ("wide lanes from nearby pulsars") names the three loosest F2
  rungs as J0711-6830 (0.106 kpc), J1630+3734 (0.089 kpc), J0437-4715 (0.155 kpc) — the
  **nearest** pulsars, hence the **shortest** lags.
- **Payoff.** The chirp-mass lever is `∂Δφ_p/∂log10 mc ∝ ḟ τ_p² ∝ τ_p²`.

**The pulsars that are easiest to certify are the worst chirp-mass measurers, and the ratio goes
as τ³.** J0711 (τ = 0.220 kyr) has the array's loosest registration lane and an mc lever of
0.408 rad/dex; B1937+21 (τ = 7.768 kyr) has a lane 35× tighter and a lever **4100× larger**.

Nearby pulsars therefore appear a fourth time — but on the **opposite side of the ledger**. The
design-theorem target list ("build the certification target list by crossing the readable-
pulsar-term sub-array WITH DISTANCE") optimises the wrong objective if the goal is the standard
siren. The right objective is a **lag-diverse** set: short lags to certify and to pin f_gw, long
lags to carry Mc. The frequency-ranked GEO set (0.22–1.22 kyr) already spans 0.74 dex of lag and
is a decent compromise; a set chosen for lag diversity would do substantially better.

---

## 7. HOW MUCH THE ANSWER LEANS ON THE CERTIFICATION PREMISE

σ(log10 Mc) at (f = 10⁻⁸, Mc = 10⁹), by arm:

| arm | N=0 | N=1 | N=2 | N=3 | N=5 | N3 longest-lag |
|---|---|---|---|---|---|---|
| A `exact-L` | 0.8653 | 0.3258 | 0.3049 | 0.1429 | 0.0160 | 0.01876 |
| **B `0.1 pc`** | 0.8653 | 0.3258 | 0.3051 | 0.1448 | **0.0251** | 0.01877 |
| C `½-fringe` | 0.8653 | 0.3259 | 0.3059 | 0.1533 | 0.1266 | 0.01880 |

Two readings, both worth keeping:

1. **Sub-fringe distance knowledge is load-bearing for the short-lag seeds and irrelevant for
   the long-lag ones.** Arm C degrades the N5 answer by **5.0×** (0.0251 → 0.1266) but the
   long-lag triple by **0.2 %**. Mechanism: a short-lag pulsar's Mc contribution rides on the
   *linear* phase 2πfτ_p, which is exactly degenerate with L_p (τ_p ∝ L_p); a long-lag pulsar's
   rides on the τ² chirp, which is not. **The 0.1 pc "within-fringe" premise is what makes the
   short-lag seeds useful at all.**
2. **The seeds' pulsar-term phases move by O(1) radian over ±1σ in Mc.** Conditional phase
   wander `max_p g_mc,p · σ(log10 Mc)` (rad), arm B:

   | log10 f | log10 Mc | N=3 | N=5 | N3 longest-lag |
   |---|---|---|---|---|
   | −8.0 | 9.0 | 6.57 | 1.14 | 31.4 |
   | −7.7 | 9.5 | 6.68 | 1.38 | 10.9 |

   These are **conditional** excursions (all other parameters held). Along the *profile*
   direction the Fisher's construction guarantees a 0.5-nat lnL drop at 1σ, i.e. a residual
   phase of ~1/SNR ≈ 0.02 rad — the compensating motion of (f_gw, sky, L_p) is real. But it is a
   **first-order** compensation, and for the long-lag set it must cancel tens of radians. Under
   the campaign's premise this is exactly right — with the fringe **integers fixed** the
   objective is a smooth quadratic in the source parameters (this is L2c's own finding, seen
   from the (Mc, f) side) — but it means the long-lag numbers lean hardest on certification. A
   full posterior with the integers free would show fringe-aliased secondary modes in (Mc, f);
   adding modes only widens the credible interval. **Every σ here is a lower bound.**

---

## 8. WHAT THIS DOES AND DOES NOT LICENSE

1. **The payoff chain is real and quantified.** Certified pulsar terms → Mc from the frequency
   ladder across kyr lags → D_L from the amplitude. Numbers: 0.866 dex → 7×10⁻⁴–0.03 dex on
   log10 Mc, 332 % → 6–12 % on D_L. **Conditional on N_seed.** Certification achievability is
   unchanged and remains a NO-GO from a cold start.
2. **"How many pulsars?" is the wrong question.** Lag *diversity* and lag *length* dominate
   count. Three well-chosen seeds beat five badly-chosen ones. Two seeds at the same lag are one
   seed.
3. **The certification target list and the siren target list are different lists**, and they
   anti-correlate as τ³. This is new, it is measured, and it should change how a target sub-array
   is designed (design theorem lever (i)).
4. **σ(D_L)/D_L is floored by the amplitude, not by Mc.** Beyond N_seed ≈ 3–5 the siren improves
   only with SNR. The floor is `ln10·σ(log10 h)`, reproduced to 0.1 %. σ(log10 h) is itself
   inflated ~3.5× over `1/(ln10·SNR)` by the h–cos_inc–ψ degeneracy — **that** degeneracy, not
   the chirp mass, is what a future siren programme must attack.
5. **These are Cramér–Rao bounds on zero-noise Asimov data with the fringe integers given.** A
   noisy realisation scatters; a full posterior with free integers can only be wider. Frozen GP
   hyperparameters (≤ 9 %, D5) and a single source (no confusion penalty) are the standing
   optimisms. The sky is free, which is the one place SIREN is *conservative* relative to the
   B1 targeted scenario.
6. **The grid is at fixed loudness, so the implied D_L is often unphysical.** At
   log10 h = −13.25 the cells span D_L = 0.33–53 Mpc; the four cells with D_L < 2 Mpc are
   astrophysically absurd for an SMBHB and are retained only because σ(D_L)/D_L depends on
   (f, Mc, h) and not on D_L separately. Read the grid as a map over *lever arm × SNR*, and
   read off the cell whose (f, Mc, SNR) you actually have. A source of the same Mc at 10× the
   distance has 10× lower SNR, and both σ(log10 Mc) and σ(log10 h) degrade ∝ 1/SNR once out of
   the prior.
7. **`log10_mc` does not enter the amplitude.** Verified in `discovery/deterministic.py`
   (`alpha = mc^{5/3}/(dist·ω^{1/3})` with `dist ∝ mc^{5/3}/h`, so mc cancels) and empirically:
   the optimal SNR is constant to 5 significant figures across the Mc column (44.46 / 44.46 /
   44.47). The chirp-mass information is *purely* an evolution/phase measurement — the spec's
   doc item (c), confirmed on the amplitude side.
8. **NOT measured here:** noise realisations (zero-noise throughout); source-sky redraws (one
   fiducial sky — GEO shows the seed set and the lags are sky-conditional, so the *lag spread*
   of a given N_seed set will vary from sky to sky and so will the payoff); eccentric harmonics
   (design-theorem lever (ii), which populates the same τ-lever channel); multi-source
   confusion; T-scaling (the Earth-term chirp phase ∝ T², so the N_seed=0 row moves with T while
   the seeded rows barely do — the *relative* payoff of seeds shrinks as T grows).

---

## 9. DIVERGENCES AND THINGS I STOPPED FOR

**Nothing unexpected occurred in the science.** Two bookkeeping items, reported and not patched:

1. **The "HELD list" does not exist in this checkout.** The brief instructs me to read
   "payoff-chain entries in the HELD list". `grep -rn "HELD"` over the repo returns exactly one
   hit — `trackB_estimator_spec.md:901`, *"Queue head = E-TRACK (spec in the HELD list)"* — and
   there is no such file in the working tree, in `git ls-files`, or in `git log --all`. Same
   shape as HARBOR's missing `PORTER_slurm_plan.md`. **If it exists on cronus, push it**; I
   proceeded from the payoff claim as stated in the brief plus the design theorem, R POSTMORTEM,
   and STEP-1 sections of the spec. I did not invent a payoff-chain entry.

2. **The rogue non-Slurm GPU process is still on dgx03, and it has moved cards.** Array task 2
   was allocated `GPU-ec71dc3d…` and found a foreign resident on it: pid 1393099,
   **51,286 MiB**, `env/deerdiff/bin/python`. HPC_SETUP §7.2 documents this pattern (a
   non-Slurm SSH session squatting cards Slurm still hands out) and §7.3 conclusion 4 establishes
   that contention changes **timing only, never values**. Task 2 ran 421 s vs 407/411 s for the
   clean tasks; its nine cells are used unreservedly. Note the squatted UUID differs from the
   three GEO recorded — the squatter restarts on different cards, as HPC_SETUP predicted. **This
   is still unreported to ACCRE admins.**

Also worth recording, though it is not a divergence: `Earth-only Asimov lnL(truth)` prints as
`405413.5127`, identical to the full-model value. That is not a coincidence and not a bug — with
residual(truth) = 0 the log-likelihood collapses to the (model-independent) log-determinant
terms. It is a free consistency check that the Earth-only injection is self-consistent.

---

## 10. EXECUTION RECORD

| item | value |
|---|---|
| lane | `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1` |
| gate/warm job | `12455052`, dgx03, **1h19m51s**, exit 0, ALL GATES PASS |
| array | `12456288_[0-2]`, 3 tasks × 3 cells, all COMPLETED 0:0, **407 / 411 / 421 s** |
| walltime requested | gate 3 h (2 cold compiles); array **30 min**/task (4× headroom) |
| **cold HVP compile, ncw=16 (24 cols)** | **2586 s** — paid once, in the gate job |
| **cold HVP compile, ncw=1 (26 cols, chunk 8)** | **1759 s** — paid once, in the gate job |
| warm per-process cost | build 162 s + graph materialisation **227 s** = 389 s fixed |
| science cost | **~0.4 s per Fisher**; 21 Fishers/task ⇒ ~15 s of the 411 s task |
| batching | 3 cells × 7 seed sets × 3 arms per task — fixed cost amortised 63× (HPC_SETUP §7.4) |
| idempotence / resume | each cell writes its own `siren_cell_fI_mcJ.npz` and is skipped if present |
| graph shapes | **one** `hvp` graph for the entire campaign: `sel` is fixed (8 source params + the 18 union distances); θ, data and mask are all runtime args |
| VRAM | comfortable; 26 columns in 4 chunks of 8, sized against the 30 GB squatted-card plan |
| disk | **1.3 MB total** (9 lean cell npz @ 20 kB + gate + geom + probe + summary + figure) |
| cache | 5286 → 5342 entries |

**Startup device log, per task (the convention).** Every job's first `[SIREN]` line is the GPU
UUID + foreign-resident query, after `sleep 12` to tolerate the previous tenant's CUDA-context
teardown ghost (HPC_SETUP §7.1). Three tasks walked three distinct A100 UUIDs; one carried the
rogue (§9.2).

### Files written (all untracked)

```
SIREN_payoff_chain.md                 this report
SIREN_results/siren_probe.npz         geometry + injection bookkeeping (CPU pre-flight)
SIREN_results/siren_gate_m3.npz       the M3 reproduction gate (H, Fs, sigmas, SNR)
SIREN_results/siren_geom.npz          exact-phase autodiff levers g_mc, g_f, dphi, tau, dL
SIREN_results/siren_cell_f{0..2}_mc{0..2}.npz
                                      per cell: 7 seed sets x 3 arms x
                                      {sig_log10_mc, sig_log10_fgw, sig_log10_h, sig_log10_DL,
                                       frac_DL, sky_deg2, cond, inv_resid, min_eig, sig_h_cond,
                                       sig_phys[8]} + snr, lnL, sym_resid, zero_off_max,
                                       seed names/lags, dL_union18
SIREN_results/siren_summary.npz       consolidated (3,3,7) arrays
SIREN_results/siren_payoff.png        3-panel: sigma(log10 Mc), sigma(D_L)/D_L, sky area vs N_seed
hpc_harbor/siren/siren_probe.py       CPU pre-flight (P1-P4)
hpc_harbor/siren/siren_payoff.py      the driver (gate | warm | run)
hpc_harbor/siren/siren_geom.py        CPU lever-arm diagnostic
hpc_harbor/siren/siren_consolidate.py tables + figure
hpc_harbor/siren/siren_gate.sbatch    warm-cache + value gate (3 h)
hpc_harbor/siren/siren_array.sbatch   3 x 30 min array
hpc_harbor/logs/siren_*.{out,err}     logs
hpc_harbor/logs/siren_consolidate.txt full tables
```

Nothing was committed. Nothing was pushed. No tracked file was edited.
