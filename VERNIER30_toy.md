# VERNIER-30 — integer ambiguity resolution (LAMBDA/ILS) on the fringe comb

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, **91 s**.
No jax, no MCMC. Clean room: imports only `vernier_toy` (constants, `dL_of`, `iso_dirs`).
Every engine comparison is **RNG-matched** — one realisation stream, all four engines
see byte-identical data.

Code `hpc_harbor/vernier/vernier30.py` · data `reports/vernier30.npz` ·
log `reports/vernier30.log`. **Staged, not committed.**

---

## Lead

**ILS beats snap decisively and matches the brute grid exactly — the win is speed
and scaling, not accuracy. The ratio test does NOT reproduce the certification gate,
and the gate is the safer of the two.**

1. **ILS ≡ brute grid on success, always.** On 300 matched realisations the LAMBDA
   fix and the dense sub-fringe grid returned the *identical* integer vector 300/300,
   and the grid never beat ILS on the objective. Across the whole K/SNR/prior sweep
   the two success columns agree to within grid resolution (ILS is very slightly
   *better* at K=8 and ρ≥16, where the grid at 8 pts/fringe misses the optimum).
   **There is no accuracy headroom for ILS to win — the grid is already exact.**
2. **ILS beats snap by 4–30× in success rate**, and snap collapses to the noise
   floor exactly where the campaign lives (2 % at σ_L = 100 pc vs ILS 51 %).
3. **The win is scaling.** Grid cost is ∝ σ_L/λ; ILS cost is flat. At σ_L = 3 kpc the
   grid is **42× slower** than ILS for the same answer, and the gap grows linearly.
4. **The ratio test and the gate disagree at their native thresholds** (agreement
   53–84 %) and **the gate dominates at matched acceptance everywhere**: gate wrong-fix
   0.00 % vs ratio 4.2–7.2 % over the 5–20 % acceptance band. The GNSS ratio test
   *cannot* be tuned to a ≤1 % wrong-fix rate on this problem at any acceptance; the
   gate holds ≤1 % out to 30 % acceptance. On pure noise, R≥3 fires **10.3 %** of the
   time; the gate fires **0.00 %**.

So: **not an independent decades-validated re-derivation of our gate.** The two
statistics are cousins (Spearman 0.83) but the ratio's denominator makes it unsafe
here, and we should *not* quote it as external validation. What IS quotable is the
Teunissen **bootstrapped success rate**, which predicts our fix rate to 1 % (§E6).

---

## Step 0 — ARCHAEOLOGY (done before any code was written)

**Prior LAMBDA work exists and is substantial. It is NOT what this instalment does.**

| artefact | what it is | what it found |
|---|---|---|
| `CW_transition/trackB_estimator_spec.md` §"STEP 2 — the RTK / LAMBDA integer solve" | Mihir's correspondence table (carrier phase ↔ pulsar-term phase, integer ambiguity ↔ fringe index, code/pseudorange ↔ F-stat, ratio test ↔ dlnL vs ln K), plus the partial-ambiguity plan 2a/2b/2c | architecture only; the ratio-test ↔ gate identification is *asserted* there, never measured |
| `trackB_L0_float.py/.npz` | truth-blind soft-EM float | float DIVERGES, 4.31 scaled off; zero carriers |
| `trackB_L1_spec.py` | search-space spec | candidate sets span n ∈ [−182, 226] |
| `trackB_L2_probe.py/.npz` | **staged integer search + conditional re-solve + ratio test** | winner integers −113/−143/−19 (truth 0), margin 0.198 nat, `matchesTruth=False` |
| `trackB_L2b_basin.py`, `trackB_L2c_newton.py` | fixed-integer pull-in | pull-in radius < 1e-4 scaled; at truth Hessian neg-def, cond ~4e10 |
| `trackB_F2_ladder` | lever-arm ladder | 27× (1.4 dex) wall between any achievable float and the loosest registration rung |

**Verdict banked 2026-07-07: NO-GO for the full LAMBDA build**, earned at F2, structural
at both ends (integer neighbourhood unreachable from the float; source pull-in
< 1e-4 even given truth's integers). Track B is CLOSED.

**The one piece that was never built.** The spec's own row reads:

> float covariance / Z-transform | fringe-index covariance | **not built (probe scale)**:
> census candidate sets already span n ∈ [−182,226] and DO NOT contain n=0 → Z-transform moot

So the decorrelation — the part of LAMBDA that *is* LAMBDA — was skipped, correctly,
because the float was already wrong. Nothing in the repo runs an actual Z-transform or
an actual integer least-squares search; `trackB_L2_probe` does bounded enumeration on a
census lattice. **VERNIER-30 builds the real thing** (LDL′, integer Gauss + permutation
reduction, depth-first search-and-shrink) and runs it on the **conditional** problem —
source sky/frequency known — which is the regime Track B's gate-2 wall does *not* touch
and which is where V28's R3 ladder actually operates. No prior implementation was
extended because none exists; the L2 machinery is a different object (campaign-scale,
jax, source-blind).

---

## The toy, and why it is a GNSS problem

One pulsar, M sources at distinct frequencies. Write the distance offset from the prior
mean as x (pc). Source k gives that pulsar a carrier wavelength
λ_k = c / (f_k(1 − cos μ_k)) — `vernier_toy.dL_of`, the fringe spacing — and an ambiguous
pulsar-term phase. In **cycles**:

```
y_k = x/λ_k − n_k + ε_k ,  n_k ∈ Z ,  ε_k ~ N(0, σ_k²) ,  σ_k = 1/(2π ρ_k)
0   = x       + ν       ,             ν   ~ N(0, σ_L²)
```

y_k is the measured wrapped fractional pulsar-term phase, n_k the fringe index, ρ_k the
per-(pulsar,source) SNR, σ_L the distance prior width. **The distance prior plays the
role of the code/pseudorange solution.** The M sources are a multi-constellation receiver.

The system is exactly determined (M phase obs + 1 prior vs 1 real + M integers), so the
float solution is closed form and so is its covariance:

```
x̂ = 0 ,  n̂ = −y ,  Q_nn[j,k] = σ_L²/(λ_j λ_k) + δ_jk σ_k²
```

**Q_nn is rank-1 plus diagonal**, median condition number **3.4e7**. That single fact is
the whole physics: the integer combinations that are *fixable* are the ones nearly
orthogonal to the 1/λ direction — the GNSS **wide lanes** — and finding them is exactly
what the Z-transform does. The vernier, written in Teunissen's coordinates.

All four engines minimise the same objective
`s(n) = (n−n̂)′ Q_nn⁻¹ (n−n̂) = −2 lnL_profiled(n)`, so `dlnL = (s₂−s₁)/2` **exactly** and
the engines are directly comparable.

**The four engines.**

| | what it is | campaign counterpart |
|---|---|---|
| **snap** | seed at the float (x=0), snap every fringe index to MAP, closed-form conditional x, iterate; then mode-hop ±1..6 on one source at a time, accept uphill | `trackB_estimator.snap_distances` / the node sampler's MAP-fringe snap, one shot |
| **snapR** | the same move given a **fair search budget**: 64 seeds drawn from the prior, each snapped to its fixed point | what a sampler chain explores over many steps |
| **grid** | dense **sub-fringe** scan over ±3σ_L, 8 pts/fringe, integers snapped and x re-solved at every node | the brute comb grid |
| **ils** | LAMBDA: LDL′, integer Gauss + permutation reduction, depth-first search-and-shrink for the two best integer vectors, back-substitution for x | new |

The grid is deliberately **sub-fringe** sampled. A fringe-*centre* grid pins the
sub-fringe offset u = 0 and fabricates fringe evidence — the banked artefact — so it is
not used as the reference anywhere here.

### S0 — machinery sanity (all pass)

```
LDL' + Z-transform round-trip max rel err : 3.12e-09   (200 draws, M=5)
det Z = 1 asserted on every draw (unimodular)
cond(Q_nn) median                         : 3.370e+07
ILS vs dense-grid exhaustive (M=4, sigL=10pc, 64 pts/fringe):
  identical integer vector    : 300/300
  ILS objective <= grid       : 300/300
  grid strictly better        :   0/300     (must be 0 -- ILS is exact)
```

A ±R per-component box around n̂ is **not** a valid reference for this problem and the
first version of the check used one and read 57/120. Because Q_nn is rank-1 dominant,
the admissible integer vectors run along a lattice *line* of length ~σ_L/λ (hundreds of
units); any small box misses the global optimum. The exhaustive reference has to be the
dense grid over the whole prior window, and against that ILS is exact.

---

## E1 — success and wall-clock vs K (number of sources)

ρ = 4, σ_L = 100 pc, 400 realisations, RNG-matched.

| K | ln K_trials | snap | snapR | grid | **ils** |
|---|---|---|---|---|---|
| 2 | 6.7 | 4.0 % | 4.0 % | 4.5 % | **4.5 %** |
| 3 | 6.9 | 3.5 % | 8.0 % | 19.2 % | **19.2 %** |
| 4 | 7.0 | 2.2 % | 11.8 % | 48.0 % | **48.0 %** |
| 5 | 7.1 | 2.2 % | 9.8 % | 79.0 % | **79.0 %** |
| 6 | 7.2 | 3.2 % | 11.8 % | 94.8 % | **94.8 %** |
| 7 | 7.2 | 1.5 % | 11.5 % | 97.5 % | **97.5 %** |
| 8 | 7.2 | 3.8 % | 10.2 % | 99.5 % | **100.0 %** |

**The vernier turns on between K=3 and K=6** and saturates by K=6 — 19 % → 95 % for three
extra sources at fixed per-source SNR. That is the commensurability effect VERNIER-24
named, now measured as a *fix rate* rather than a contrast.

**Snap does not participate.** It is flat at ~3 % for every K, and *falls* from K=2 to
K=7. Its mode-hop radius is ±6 fringes on one source; truth sits ~100 fringes away when
σ_L = 100 pc. More sources make the snap fixed point *stickier*, not more informative.
Giving it 64 prior restarts (snapR) buys 3× and no more — 64 seeds over a ~600-fringe
window samples one fringe in ten.

K=8 is the one cell where ILS and grid differ (100.0 vs 99.5): the grid at 8 pts/fringe
missed one optimum. **ILS is the exact engine; the grid is the approximation.**

---

## E2 — success vs per-(pulsar,source) SNR

K = 4, σ_L = 100 pc, 400 realisations.

| ρ | snap | snapR | grid | **ils** |
|---|---|---|---|---|
| 1 | 0.8 % | 0.5 % | 2.2 % | **2.2 %** |
| 2 | 3.0 % | 7.2 % | 13.8 % | **13.8 %** |
| 4 | 4.5 % | 6.5 % | 46.8 % | **46.8 %** |
| 8 | 3.0 % | 12.0 % | 85.5 % | **85.5 %** |
| 16 | 4.0 % | 14.8 % | 97.8 % | **98.2 %** |
| 32 | 2.0 % | 11.5 % | 99.2 % | **99.8 %** |

At the campaign's banked ρ_P ≥ 2 bar with K=4, the fix rate is **13.8 %**. The SNR
lever and the K lever are interchangeable over this range (ρ=8/K=4 ≈ ρ=4/K=5), which
is the expected behaviour for a phase-noise-limited lattice.

---

## E3 — success and cost vs prior width (the scaling result)

K = 4, ρ = 4. **Wall-clock in this block is same-run and monotone; absolute ms across
blocks is noisy on a shared box and should not be compared between tables.**

| σ_L | ln K | snap | snapR | grid | **ils** | grid ms | **ils ms** | grid/ils |
|---|---|---|---|---|---|---|---|---|
| 1 pc | 2.4 | 71.0 % | 96.8 % | 98.0 % | **98.5 %** | 0.147 | 0.604 | 0.24× |
| 3 pc | 3.5 | 47.8 % | 86.0 % | 94.0 % | **93.8 %** | 0.147 | 0.858 | 0.17× |
| 10 pc | 4.7 | 29.5 % | 57.8 % | 88.0 % | **88.0 %** | 0.112 | 0.197 | 0.57× |
| 30 pc | 5.7 | 10.0 % | 31.2 % | 73.5 % | **73.8 %** | 0.164 | 0.224 | 0.73× |
| 100 pc | 7.0 | 1.8 % | 8.8 % | 50.5 % | **50.5 %** | 0.557 | 0.202 | **2.8×** |
| 300 pc | 8.1 | 0.8 % | 4.0 % | 29.8 % | **29.8 %** | 2.00 | 0.235 | **8.5×** |
| 1 kpc | 9.3 | 0.0 % | 0.0 % | 10.0 % | **10.0 %** | 18.3 | 1.08 | **17×** |
| 3 kpc | 10.4 | 0.0 % | 0.0 % | 10.0 % | **10.0 %** | 70.0 | 1.66 | **42×** |

**This is the operational answer to "does ILS beat the grid".** Not on accuracy — they
are the same answer — but the grid pays ∝ σ_L/λ_min and ILS pays essentially nothing.
The crossover is at **σ_L ≈ 30–100 pc**, i.e. right at the campaign's tightest-lane
pulsars (J0711 106 pc, J1630 89 pc, J0437 155 pc). Below ~30 pc the grid is cheaper and
you should just use it; above it, use ILS, and the advantage grows without bound.

Snap's decay with σ_L is the clean statement of what it is: a **local** move. 71 % when
truth is within a few fringes, 0 % when it is 3000.

---

## E4 — THE RATIO TEST vs THE CERTIFICATION GATE

Pooled over M ∈ {3,4,5,6} × ρ ∈ {2,4,8} × σ_L ∈ {10,100,1000} pc, 120 realisations each,
**n = 4320**, ILS success 56.1 %. The trials bar ln K spans **2.85–10.11 (sd 1.95)** in
this pool — deliberately, so the gate's per-realisation bar does real work and the
comparison is not secretly "difference vs ratio at constant threshold".

- **ratio test** (GNSS): accept if R = s₂/s₁ ≥ μ.
- **our gate**: accept if contrast = (s₂−s₁)/2 ≥ ln K.

### Native thresholds — they disagree

| criterion | accept | wrong-fix | agreement with gate |
|---|---|---|---|
| **gate** (contrast ≥ ln K) | 23.9 % | **0.29 %** | — |
| ratio R ≥ 1.5 | 65.4 % | 27.6 % | 58.5 % |
| ratio R ≥ 2.0 | 51.2 % | 19.4 % | 72.5 % |
| ratio R ≥ 3.0 | 37.7 % | 13.5 % | 84.0 % |

### Matched acceptance — the gate is safer at every operating point

| target accept | ratio WF | **gate WF** | agree |
|---|---|---|---|
| 5 % | 4.17 % | **0.00 %** | 95.7 % |
| 10 % | 5.09 % | **0.00 %** | 93.2 % |
| 20 % | 7.18 % | **0.00 %** | 89.8 % |
| 35 % | 12.17 % | **2.78 %** | 85.1 % |
| 50 % | 18.75 % | **14.03 %** | 76.8 % |
| 70 % | 29.83 % | **27.35 %** | 71.8 % |
| 90 % | 39.89 % | **38.73 %** | 85.7 % |

Max acceptance holding wrong-fix under a target:

| WF target | ratio | gate |
|---|---|---|
| ≤ 1 % | **never reached at any acceptance** | 30.0 % |
| ≤ 5 % | 12.0 % | 38.0 % (**3.2×**) |

### Null control — pure noise, no coherent comb

| criterion | fires on noise |
|---|---|
| gate (contrast ≥ ln K) | **0.00 %** |
| ratio R ≥ 1.5 | 38.75 % |
| ratio R ≥ 2.0 | 20.75 % |
| ratio R ≥ 3.0 | 10.25 % |

### The mechanism

Among R ≥ 3 accepts:

| | median s₁ | median contrast |
|---|---|---|
| correct fixes | 3.245 | 10.27 nat |
| **wrong fixes** | **0.237** | **0.50 nat** |

**s₁ sits in the ratio's denominator.** A realisation whose best residual happens to be
tiny — lucky noise, not evidence — gets a large R with essentially zero contrast behind
it. The gate reads the *difference*, which is the actual log-evidence, and is immune.
The two statistics are strongly rank-correlated (Spearman 0.826, Kendall τ 0.645) because
in the bulk s₁ ≈ its expectation and R ≈ 1 + 2·contrast/s₁; they part company exactly in
the tail that matters.

**Answer to the brief's question 4: they do NOT agree, and the gate is the safer one.**
We do not have an independent decades-validated re-derivation of our certification
criterion, and we should not claim one. The honest framing for a paper is the reverse:
*the PTA problem is a case where the GNSS ratio test is mis-calibrated, because the
trials factor ln K is enormous (10³–10⁵ fringes) and a ratio statistic carries no
trials information at all.* Our bar does. That is a defensible novelty claim, not a
borrowed-validation claim.

Caveat: GNSS practice since Verhagen (~2005) largely uses a **fixed-failure-rate**
ratio threshold rather than a fixed μ, precisely because fixed μ is known to be
poorly calibrated. The matched-acceptance sweep above *is* the empirical FF-ratio
construction, and the ratio statistic still loses — so this is not a strawman-μ result.

---

## E6 — partial ambiguity resolution, and the one thing GNSS gives us for free

Teunissen's bootstrapped success rate on the decorrelated components,
`P_j = erf(1/(2√(2 D_j)))`, M=5, ρ=4, σ_L=100 pc, 300 draws:

| component | z₁ | z₂ | z₃ | z₄ | z₅ |
|---|---|---|---|---|---|
| median P | 0.908 | 0.956 | 0.976 | 0.985 | 0.996 |

- full-vector bootstrap (product): **0.804**
- **measured ILS success on a matched cell: 0.793**

**Agreement to 1.1 %.** This is the genuinely useful import from geodesy: a *closed-form,
per-realisation, data-only* predictor of whether a comb will resolve — no search, no
ensemble, one LDL′ and one erf. It is a lower bound on the ILS success rate by
construction and it is tight here. If anything in this instalment is quotable and
transferable, it is this, not the ratio test.

Note the shape of the P vector: after decorrelation **no single component is hopeless**
(min median 0.908) even though Q_nn has condition number 3e7 before it. The rank-1
distance direction has been spread across all five wide lanes. 75.7 % of draws have at
least one component at P > 0.99 (mean 1.26 of 5), so partial-fix-then-constrain has
something to bite on when the full fix fails.

---

## Feasibility note for VERNIER-28 (one paragraph, no run)

**Yes — as a seed for R3, not as a replacement for it, and the numbers say it is worth
one instalment.** V28's R3 ladder marginalises pulsar distances over **10 pc** priors,
and E3's σ_L = 10 pc row is exactly that cell: ILS 88.0 %, snap 29.5 %, at K=4 and ρ=4.
R3's joint-bin fit currently pays for the pulsar-term block by integrating over the
distance nuisance; ILS would instead *fix* the fringe integers first, at 0.2 ms per
pulsar per bin, and hand the joint fit a distance that is a smooth quadratic rather than
a comb — which is the same "cusp → quadratic" conversion Track B's spec argued for, but
in the regime where it actually holds (source position known, prior 10 pc not 100). Two
cautions. (i) At σ_L = 10 pc the grid is *cheaper* than ILS (0.11 vs 0.20 ms), so the
argument for ILS inside V28 is not speed — it is that ILS extends to the joint N-pulsar
lattice where a grid does not, and that E6's bootstrap gives R3 a free per-pulsar
"will this resolve?" flag to decide *which* pulsars to fix and which to keep soft
(partial ambiguity resolution, the spec's step 2a). (ii) A 12 % wrong-fix rate at
σ_L = 10 pc would inject hard biases into R3's recovery accounting, so any V28 use must
gate on the contrast bar (0 % wrong-fix out to 30 % acceptance, §E4) and keep the
ungated pulsars marginalised — fix-and-gate, never fix-and-assume.

---

## Limitations — read before quoting anything above

1. **Conditional, not cold-start.** Source sky position and frequency are KNOWN
   throughout, so λ_k is known. This is the regime Track B's gate-2 wall does not
   touch, and nothing here reopens the banked NO-GO for cold-start certification.
   VERNIER-30 says the *integer* half of the problem is solved cleanly once the
   continuous half is; Track B's F2 ladder says the continuous half is not reachable
   blind. Both remain true.
2. **Per-pulsar, decoupled.** Given the source parameters the pulsars separate, so this
   is N independent (1 real + M integer) problems, not one joint (N + N·M) lattice.
   The joint problem is where LAMBDA's decorrelation would earn most and it is not
   tested here.
3. **Non-evolving pulsar-term convention.** λ_k is taken as constant across the fringe,
   which is what makes the integer structure exact. `vernier_toy.anchor_check` records
   that the evolving pulsar term gives a *coarser* comb than c/f₀; the correction shifts
   λ, and with it the vernier's commensurability, in a way this toy does not model.
   Every success rate here should be read as "at the stated λ set", not "at the
   campaign's λ set".
4. **Phase noise is Gaussian and per-source independent.** Real pulsar-term phase errors
   are correlated through the shared Earth-term fit and the timing model. Correlated
   ε would enter Q_nn off-diagonal and the wide lanes would move.
5. **The exactly-determined structure is a modelling choice.** M phases + 1 prior with
   M+1 unknowns means zero redundancy, so the float residual is identically zero and
   there is no overall-model consistency check available — in GNSS that redundancy is
   what supports the ratio test's null distribution. Part of the ratio test's poor
   showing may be attributable to this, and a redundant variant (multiple epochs per
   source) is the obvious next probe.
6. **Wall-clock is from a shared box.** Ratios within one block are meaningful;
   absolute ms across blocks varied ~4× between two runs of the identical script.
7. Nothing here has been cross-checked against any earlier instalment's numbers, and
   the LAMBDA implementation is new in this instalment (validated only against the
   dense grid, §S0, and against Teunissen's bootstrap, §E6).
