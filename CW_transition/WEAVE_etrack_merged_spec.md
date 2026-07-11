# WEAVE — The E-track as the Eccentric Generalisation of the Chirp Cascade

**Agent:** WEAVE · **Date:** 2026-07-09 · Read-only + this file. CPU/numpy-scale checks only.
**Supersedes the scoping brief of** `LANES_eccentric_ladder.md` **as the E-track's objective**
(it does not overturn LANES' numbers; it re-uses them under a changed question).

**Sources:** `project_progress.md` §B1-STEP-1 (2026-07-09),
`CW_transition/trackB_estimator_spec.md` (design theorem; DELIVERABLE R; B1 STEP 1/1A/1B/1C),
`CW_transition/LANES_eccentric_ladder.md`. Anything not in those is tagged
**[DERIVED-WEAVE]** (my algebra + numpy, reproducible from this file) or **[UNSOURCED]**.

> **Doc flag.** The brief named `PROJECT_PROGRESS.md`. No such file exists; the progress doc is
> `project_progress.md` (lowercase, repo root, 135 KB). Read at lines 1780–1827 (B1 STEP 1 entry).

---

## 0. The merger in one paragraph

B1's `(f, mc)` pilot found that `log10_mc` is a **registration axis**, not a shape parameter:
`make_phase_connected_binary` builds the pulsar term at the retarded time with the full PN chirp
(`evolve_phase(tp)`, `discovery/deterministic.py:522`), so the pulsar-term phase depends on `Mc`.
Spec §STEP-1(c) states the physical identification verbatim:

> "The pulsar term is a kyr-baseline timestamp of the source's phase. What registers it is not the
> source's frequency but its frequency EVOLUTION: the fringe index is set by the accumulated chirp
> over `tau_p`, so the array must know `fdot` (hence `Mc`) to place any fringe. The 22-yr Earth term
> cannot measure `fdot` (`Delta_phi_E ~ 0.05 rad`). This is exactly the E-track's eccentric-harmonic
> mechanism in CIRCULAR form … **The E-track and the targeted pipeline have merged**: lever (ii) is
> not an alternative to lever (iii), it is the missing ingredient lever (iii) needs."

So the E-track's question is no longer "does eccentricity help?" It is: **the circular chirp is the
n = 2 special case of a harmonic clock; what does the full comb do to the cascade?** LANES asked
whether eccentric harmonics widen a *lane* (they do not — §2.4). The pilot shows the obstruction was
never lane width; it is the **mc box**. That is a different quantity, and it is what §2–§3 measure.

---

## 1. THE CIRCULAR BASELINE — the mc cascade as established

All numbers in this section are quoted verbatim from `trackB_estimator_spec.md` §"B1 STEP 1
(2026-07-09)". Zero-noise Asimov, truth-placed (LABELLED), pop config N_CW=16 / 116 psr,
`tau_p` median 3.52 kyr (min 0.019, max 86.0).

### 1.1 The registration ladders (pilot M1, 1-rad tolerance, SCALED)

| axis | loosest | median | tightest | # ≥ 1e-2 | # ≥ 1e-3 |
|---|---|---|---|---|---|
| sky | 1.852e-3 | 4.472e-5 | 3.847e-6 | 0 | 3 |
| log10_fgw (chirped) | 2.518e-2 | 1.435e-4 | 2.062e-5 | 2 | 27 |
| **log10_mc** | **6.047e+1** | **7.836e-4** | 1.345e-5 | 53 | 161 |
| freq, F2's chirp-free object | 2.520e-2 | 1.035e-4 | 5.578e-6 | — | — |

F2's ladder used the non-chirped phase `2 pi f L (1-cos mu)/c`, whose `d/dlog10_mc` is identically
zero: **F2 was structurally blind to an mc lane.** Its sky ladder is valid (reproduced to 3 digits);
its freq ladder is a median 0.72× optimistic. Per-axis certification tolerance (M2): sky 1e-5,
`log10_fgw` 3e-5, **`log10_mc` 1e-3**, extrinsics no collapse to 1e-2. The 8 source params split
cleanly into **4 harmless extrinsics + 2 sky (counterpart-supplied) + 2 registration axes (f, mc)**.

What the Earth term supplies (M3, HVP Hessian, inverted *with* the generative priors — never pinv):
`log10_fgw` info gain 162×/389×/269×; **`log10_mc` gain 1.73×/1.00×/1.00×** (σ = 0.501/0.864/0.866
dex — for loud1 and loud2 **the posterior IS the prior**). Targeted wall: freq 203×, **mc 1730×**.

### 1.2 Ignition and rung spacing (STEP 1A)

`R_a = 1 / sqrt( sum_k [(sig_f,k g_f,ak)^2 + (sig_mc,k g_mc,ak)^2] )`, register when `R_a >= 1`.

- **95–100 % of the pulsar-term phase wander is `log10_mc`** (J0711 95.0 %, J1713 100.0 %,
  J1909 99.9 %, J0437 99.3 %).
- max `R` = **2.143e-2** (J0711-6830, 0.106 kpc); median 2.036e-4; min 7.441e-6.
- **IGNITION: 0 pulsars at `R >= 1`.** Loosest rung needs a **47× box shrink**; the binding census
  pulsar J1713 needs **3305×**.
- The ladder **spans internally** — 42 rungs between float and census target, max `log10` gap
  0.387 dex < F2's 0.7 criterion. **Unlike sky, the failure is NOT rung spacing. It is the mc BOX.**

### 1.3 The loop-gain mechanism (STEP 1B)

Fixing a subset `S` of fringes and re-solving `(f, mc)` is a **masked likelihood** — `MaskedDelay`
makes the per-pulsar pulsar-term mask a runtime arg, so `F_S = -H[lnL(theta; data_S, mask_S)]` over
the 6 registration params is one call per subset. No optimiser ⇒ **no L2c pull-in problem.** Asimov
data injected with the same mask ⇒ residual(truth) = 0 ⇒ `H = -Fisher` exactly.

| subset fixed | sig_f | sig_mc | sig_f/tol | sig_mc/tol | next-rung R = GAIN |
|---|---|---|---|---|---|
| earth only (mask=0) | 1.592e-3 | 1.710e+0 | 53.1 | 1709.6 | **0.04** |
| top-1 (J0711) | 1.347e-3 | 2.204e-1 | 44.9 | 220.4 | 0.12 |
| top-3 | 1.013e-3 | 5.766e-2 | 33.8 | 57.7 | 0.15 |
| top-6 | 3.709e-4 | 7.720e-3 | 12.4 | 7.7 | 0.19 |
| top-12 | 2.720e-4 | 4.240e-3 | 9.1 | 4.2 | 0.36 |
| top-24 | 1.407e-4 | 1.196e-3 | 4.7 | 1.2 | 0.70 |
| **top-48** | 5.073e-5 | 2.663e-4 | 1.7 | 0.3 | **1.34** |
| top-116 | 1.363e-5 | 1.530e-5 | 0.5 | 0.0 | — |
| **census-3 only** | 1.403e-4 | 1.053e-3 | 4.7 | 1.1 | **4.50** |

**The loop.** A fixed loose rung is a **kyr-baseline chirp measurement** — its pulsar-term phase
`~ pi fdot tau_p^2` is a lever arm ~`(tau_p/T)^2 ~ 2.6e4` times the Earth term's. Fixing it collapses
`sig_mc`; a smaller mc box makes finer rungs reachable; those fix more chirp. **Loop gain is 0.04 at
the float, crosses 1 between 24 and 48 fixed fringes, and is 4.50 given the census-3.**

**KILLED HYPOTHESIS (recorded, not quietly dropped).** Pre-registered guess: since
`R_a ~ 1/(sig_mc g_mc,a)`, a pulsar registers exactly to the degree it is *blind* to mc, so no
fixing order bootstraps. **REFUTED**: fixing J0711 alone cuts `sig_mc` 7.8×; the census-3 alone cut
it 1600×. **Loose rungs carry mc information through their AMPLITUDE/SNR, not only through `g_mc`.**

**The real structure is BISTABILITY.** A self-consistent registered state exists and strongly
attracts (the P1 needle from the (f,mc) side); it is unreachable from the Earth-term float.
**Certification is SELF-REFERENTIAL: ~30 registered fringes are needed to measure the chirp mass
that lets any fringe be registered.**

### 1.4 Cascade from each conditioning tier (STEP 1C) — none ignites

| tier | conditioning | σ_f (scaled) | σ_mc (scaled) | best free R | ignites? |
|---|---|---|---|---|---|
| A | sky only | 0.0101/0.0042/0.0061 | 1.002/1.727/1.731 | 0.0214 | NO |
| B | + EM period, σ_P/P = 1e-3 | 0.00145 ×3 | unchanged | 0.0220 | NO |
| C | + host D_L (σ_mc 0.036/0.022/0.021 dex) | 0.00145 ×3 | 0.0727/0.0435/0.0409 | **0.271** | NO |

- **An EM period buys nothing** (0.0214 → 0.0220): a 7× tighter `f` moves the best rung 3 %.
  Consistent with the 95–100 % mc share.
- **Only the distance moves anything**, 12.7× against the 47× needed — **Tier C misses by 3.7×.**
  Tier C's mc box is set by `log10 h = (5/3) log10 Mc + (2/3) log10 f - log10 D_L + const`, i.e. by
  the **array's own** `sigma(log10_h)` (0.060/0.036/0.034 dex) plus `sigma(log10 D_L) = 0.005` — not
  by anything the counterpart supplies beyond the redshift.
- A pre-run prediction of "misses by 1.4×" was **WRONG**: it scaled `R` by the *median* mc shrink.
  `R` is a quadrature over the 3 loud sources and the loosest rung (J0711) is dominated by loud0,
  whose `σ_mc` shrinks only 13.8× because loud0 has the worst Earth-term `σ(log10_h)`.
  **The binding shrink is the WORST source's, not the median's.** This error class is the reason for
  the wrong-fix discipline below.

### 1.5 The wrong-fix discipline (MANDATORY for every E-track table)

Rationale, from the banked Track B failures:

- **Blind search fixes WRONG integers.** LAMBDA L2: winner integers −113/−143/−19 (truth 0),
  ratio-test margin 0.198 nat, `matchesTruth = False`. A confident-looking margin at a wrong fix.
- **`dlnL(n_true vs best-wrong)` = −13.6 / −10.9 / −8.4** at the wrong source: once the source is
  off, **the truth fringe sits BELOW the wrong fringe**. A cascade that fixes and never re-audits
  will report a healthy, self-consistent, wrong registered state.
- **Fixing is soft, and only at loose rungs.** Spec §2a: resolve integers only for concentrated
  carriers (`q_max > 0.9`, census-3 bootstrapped first); the ~98 degenerate pulsars stay
  marginalised — they carry Earth-term information, not fringe information. The validated M-step
  gate is `q_max > 0.5` to reassign a distance, else hold at the prior mean / last value.
- **The R POSTMORTEM separation.** *Prior-pinned* fringe information survives joint failure
  (J0437-4715: K=4, σ_EM = 0.25 pc, 88.7 % of plateau samples still MAP to the true fringe,
  P_true = 0.568, exp H = 3.41). *Data-driven* fringe information does not (census-3 J0711/J1713/
  J1909 certify at truth with P_true 0.953/0.989/0.998 but smear on the plateau to 0.117/0.014/0.019
  with 34.6/33.9/11.0 effective wrong solutions).

> **RULE.** Every cascade table in the E-track carries a **wrong-fix column**: at each rung, the
> fraction of fixed fringes whose integer ≠ `n_true`, and the `dlnL(n_true − n_fixed)` at the
> current source estimate. A rung that ignites with a non-zero wrong-fix count has not ignited.
> Report the wrong-fix column even when it is all zeros. **[Discipline mandated by the brief;
> the individual numbers above are sourced. The phrase "wrong-fix column" is the brief's, not the
> docs'; its docs-grounding is the B1 gate's "WRONG-certification count" and the L2 result.]**

---

## 2. THE ECCENTRIC EXTENSION — n clock hands on one timestamp

### 2.1 The setup

An eccentric binary emits a phase-locked Peters–Mathews comb at `f_n = n f_orb`, `n = 1, 2, 3, …`,
with power weights `g(n,e)` (Peters & Mathews 1963, Bessel form; LANES A3) satisfying
`sum_n g(n,e) = F(e)`, `F(e) = (1 + 73e²/24 + 37e⁴/96)/(1-e²)^{7/2}` — **numpy gate: exact to
4.7e-8 at e = 0.9, machine zero for e ≤ 0.7 [DERIVED-WEAVE].**

Each harmonic is *the same clock*, read by a different hand. All `n` share one `L_p`, one `Mc`, one
`e`, one `f_orb`. Over the pulsar-term lag `tau_p` the pulsar term timestamps **all of them at once**.

### 2.2 The phase derivative, analytically

At leading order in `fdot tau / f` (the pilot's own measured smallness: the chirp changes the exact
pterm-minus-earth GW phase by a median factor **0.955** vs the non-chirped `2 pi f tau`, i.e.
`fdot tau / f ≈ 0.09`):

```
Phi_n(t) - Phi_n(t-tau) = 2 pi n [ f_orb tau - (1/2) fdot_orb tau^2 + O(fddot tau^3) ]
```

Peters:  `fdot_orb = (96/5) (2 pi)^{8/3} (G Mc/c^3)^{5/3} f_orb^{11/3} F(e)`, so
`d fdot_orb / d log10 Mc = (5/3) ln10 · fdot_orb`. Hence

```
| d(Delta Phi_n) / d log10 Mc |  =  n · (5/3) ln10 · pi · fdot_orb · tau_p^2
tol_n(Mc)  =  1 rad / (that)  =  tol_1(Mc) / n ,      tol_1 ∝ 1 / (fdot_orb tau_p^2) ∝ 1/F(e)
```

**Answer to "how does harmonic n's mc-rung spacing scale?" — `tol_n = tol_1 / n`, and the whole
ladder shifts down (tighter) by `F(e)`.** Same `1/n` as LANES' frequency lane law, for the same
reason (lever arm ∝ emitting frequency); the new content is the `F(e)` factor and the fact that all
`n` are present simultaneously in one residual.

**numpy check at the pilot's own point [DERIVED-WEAVE].** From the pilot's `tau_p` median 3.52 kyr,
exact phase median 1.302e4 rad, chirp factor 0.955 ⇒ implied `f_gw = 1.953e-8 Hz`,
`fdot = 1.583e-20 Hz/s`. Circular (`n = 2` on `f_orb`, i.e. `Phi_gw`):

```
tol_mc = 1 / ( pi (5/3) ln10 · fdot · tau^2 ) = 4.247e-4 dex = 8.49e-4 scaled
pilot M1 median (measured, autodiff of the exact phase)  = 7.836e-4 scaled = 3.918e-4 dex
ratio 1.084   (8 % — leading-order expansion vs exact PN phase)
```

*(Scale factor 1 scaled unit = 0.5 dex for `log10_mc` [DERIVED-WEAVE]: `b1_referendum_tierA.npz`
`boxhalf` mc = 1.0 with `prov` = "generative prior U(8.5,9.5)" (half-width 0.5 dex); cross-checks
M3's `sigma_prior = 1.732` scaled = 0.866 dex.)*
Independent Peters cross-check at that `f_gw`: `log10_mc = 9.0` ⇒ `fdot = 3.13e-20`,
`tol = 2.15e-4 dex`; `log10_mc = 8.5` ⇒ `1.46e-3 dex`; `9.5` ⇒ `3.15e-5 dex` — the pilot's 3-loud
draw sits inside that spread. **The analytic law reproduces the measured ladder.**

### 2.3 Ladder densification — and the cancellation that governs it

Two consequences of `tol_n = tol_1/n`, pulling opposite ways.

**(a) The comb is a better chirp-mass meter.** Per-harmonic residual SNR obeys `SNR_n² ∝ g(n,e)/n⁴`
(LANES A4 — a PTA measures timing residuals `r ~ h/(2 pi f)`), while the `Mc`-derivative scales as
`n`. So per-harmonic **Mc-Fisher weight** `w_n ∝ SNR_n² · n² = g(n,e)/n²` — the `1/n⁴` residual tax
is bought back by one factor `n²`. Define the **eccentric clock boost** at fixed total SNR and fixed
`f_orb`:

```
kappa(e) = (n_eff/2) · F(e),      n_eff = sqrt( sum_n g/n^2  /  sum_n g/n^4 )
```

**[DERIVED-WEAVE, numpy]:**

| e | F(e) | n_peak(g) | n_eff | **κ(e)** | residual SNR(e)/SNR(0) at fixed (Mc, f_orb, D_L) |
|---|---|---|---|---|---|
| 0.0 | 1.000 | 2 | 2.000 | **1.00** | 1.000 |
| 0.3 | 1.776 | 3 | 2.095 | **1.86** | 0.985 |
| 0.5 | 4.884 | 4 | 2.294 | **5.60** | 0.957 |
| 0.6 | 10.228 | 6 | 2.462 | **12.59** | 0.938 |
| 0.7 | 27.266 | 10 | 2.708 | **36.92** | 0.915 |
| 0.75 | 51.145 | 13 | 2.878 | **73.60** | 0.901 |
| 0.8 | 110.902 | 18 | 3.099 | **171.8** | 0.887 |
| 0.9 | 1243.113 | 51 | 3.865 | **2402** | 0.854 |

Note the last column: at fixed `f_orb` the residual SNR is **flat to within 15 %** across all `e` —
the `F(e)` power boost and the `1/n⁴` residual tax cancel. Eccentricity is nearly free in SNR.
Taken alone, (a) looks decisive: the Earth-term chirp phase `Delta_phi_E ≈ pi fdot T² ≈ 0.05 rad`
(circular, M3) becomes `κ·0.05` — **≈ 1.8 rad at e = 0.7**. The Earth term starts measuring `fdot`.

**(b) The same `κ` multiplies the pulsar-term wander.** `g_mc,a ∝ n_eff F(e) fdot_orb tau_a²` — the
free rung's tolerance tightens by exactly the same `κ`. So with an mc box of width `σ_mc`, where
`σ_mc^{-2} = σ_prior^{-2} + κ² I_1` (`I_1` = circular Earth-term mc Fisher):

```
R_a(kappa) / R_a(1)  =  sqrt( sig_prior^-2 + kappa^2 I_1 ) / ( kappa · sqrt( sig_prior^-2 + I_1 ) )
                     =  sqrt( sig_prior^-2 / kappa^2 + I_1 ) / sqrt( sig_prior^-2 + I_1 )
```

which is **monotonically decreasing in `κ`, and never exceeds 1.** Numerically, using M3's measured
Earth-term mc info gains (loud0 1.73×, loud1/loud2 1.00×) **[DERIVED-WEAVE]:**

| e | κ | R(e)/R(circ), loud0 (gain 1.73) | R(e)/R(circ), loud1/2 (gain 1.00) |
|---|---|---|---|
| 0.0 | 1.00 | 1.000 | 1.000 |
| 0.3 | 1.86 | 0.873 | 0.537 |
| 0.5 | 5.60 | 0.823 | 0.179 |
| 0.7 | 36.9 | 0.816 | 0.027 |
| 0.9 | 2402 | 0.816 | 0.0004 |

Taking `κ → ∞` the ratio floors at `sqrt(I_1)/sqrt(σ_prior^{-2}+I_1)`, i.e. at `R_a` evaluated with
`σ_mc` supplied entirely by the Earth-term chirp. Substituting the scalings
(`I_1^{1/2} ∝ SNR·κ·fdot·T²`, `g_mc,a ∝ κ·fdot·tau_a²`), **`κ` and `fdot` both cancel**:

> ### **The clock-cancellation ceiling [DERIVED-WEAVE]**
> ```
> R_a^max  =  C · SNR · (T / tau_a)^2       (C = O(1) covariance penalty for f–fdot–phase0)
> ```
> **independent of `Mc`, `e`, and `f_orb`.** The clock-rate sensitivity that lets the Earth term
> measure `Mc` is the *same* sensitivity that makes the pulsar term wander. Whatever raises the
> meter raises the wander. The only surviving asymmetry is `T²` vs `tau_a²`.
>
> `R_a^max = 1` requires `tau_a <= T sqrt(SNR)` = **85 yr** (SNR 15, LANES A5), **156 yr** (SNR 50),
> **374 yr** (SNR 289 = the R5 pencil's per-source SNR for the sky wall).
> Array: `tau_p` median 3.52 kyr, max 86 kyr, **min 0.019 kyr = 19 yr**.

Consistency with the pilot: `R_max(tau = 3.52 kyr, SNR 15) = 5.9e-4` vs measured median
`R = 2.04e-4` (below the ceiling ✓); `R_max(tau ≈ 0.35 kyr)` = 5.9e-2 vs J0711's measured
`R = 2.14e-2` (below ✓). *(J0711's `tau ≈ 0.35 kyr` is my estimate from L = 0.106 kpc with
`(1-cos mu) ≈ 1` — **[UNSOURCED]**, the per-pulsar `tau_p` table is in `b1_pilot_m1.npz`, not read.)*
Tier C's `R = 0.271` legitimately exceeds `R_max`: its `σ_mc` comes from the **amplitude relation**
(`h`, `D_L`), an `fdot`-free constraint outside the cancellation. **That is precisely why STEP 1C
found "only the distance moves anything."**

**Pre-registered prediction, therefore, and the thing the map must falsify:** *at leading order in
`fdot tau/f`, and treating the comb as `Phi_n = n Phi_orb`, eccentricity CANNOT ignite the cascade —
it monotonically lowers `R`, worst for the sources whose Earth term is prior-limited (loud1/2, gain
1.00×, where `R ∝ 1/κ`).* Two further leading-order signs, both unfavourable:

- **Retarded-epoch asymmetry.** The pulsar term samples the source at `t - tau_p`, i.e. **earlier**,
  i.e. **more eccentric** (`e` decreases monotonically under Peters). So `κ_pterm > κ_earth` always:
  the wander is boosted by more than the meter. The residual asymmetry has a **definite unfavourable
  sign**. [DERIVED-WEAVE]
- **Vernier ceiling.** The joint fringe period of the comb is set by the **fundamental**,
  `dL_1 = c/(f_orb(1-cos mu))`; the circular source at the same *detection* frequency has
  `dL_2 = dL_1/2`. Because residual SNR pins the detection band at `n = 1–2` for all `e` (LANES A4;
  `n_eff ≤ 3.9` even at `e = 0.9`), the candidate-count reduction is **≤ 2×**, not `n_peak×`. Against
  the R-postmortem's prior-pinning threshold (`K ≲ 4`, J0437), 2× turns `K = 72/558/1264`
  (J1909/J0711/J1713) into `36/279/632`. **The comb does not prior-pin the census-3.** [DERIVED-WEAVE]

### 2.4 Where the cancellation breaks — the E-track's actual hypothesis

The cancellation is exact **only** if every harmonic's `Mc`-derivative is `n × (the fundamental's)`,
i.e. only if `Phi_n = n Phi_orb` and `Mc` enters solely through `fdot_orb`. It does not:

1. **`e(t)` evolves over `tau_p`.** `de/dt` is itself `∝ Mc^{5/3}`, so `d Phi_n/d log10 Mc` acquires a
   term `∂Phi_n/∂e · de/d log10 Mc` whose `n`-dependence runs through `∂ln g(n,e)/∂e` — **not** `∝ n`.
2. **Periastron advance.** `gamma_dot` (1PN) has a different `Mc`, `e`, `f_orb` dependence; the
   harmonic phases split as `n l + k gamma`, breaking the `∝ n` law at 1PN.
3. **Comb geometry is `fdot`-free.** Tooth *spacing* gives `f_orb` and amplitude *ratios* give `e`
   with no chirp involved. That is genuinely new information the scalar-`κ` argument cannot see.

So: the circular clock is **rank-1** in the evolution sector (only `fdot`, only `Mc`). The eccentric
clock is **rank-3** (`Mc`, `e`, `f_orb` via `fdot`, `edot`, `gammadot`), read by `n` hands with
linearly independent derivative rows. **The E-track's hypothesis is a RANK/CONDITIONING claim, not a
magnitude claim:** does the 3×3 evolution Fisher assembled over the pulsar-term lag beat the scalar
cancellation? §3 measures that. **[DERIVED-WEAVE framing; the ingredients (1)–(3) are standard
Peters/PN, the merger is the docs'.]**

### 2.5 Folding in LANES — what still stands, and what changes meaning

**The scissors stands, and it is not what the merged E-track relies on.** LANES:

1. *Power-anchor geometry.* Placing `n_peak` at F2's 1.85e-3 rung needs `e ≳ 0.85`; but then `n = 1`
   falls below 5 % of peak, so the usable band spans only ~5–14× in `n`, centred high. Widest usable
   lane tops out at **8.6e-3 scaled at e = 0.9** — **5.8× short** of the 0.05 float ceiling.
2. *Residual-SNR reality.* `SNR_n² ∝ g(n,e)/n⁴` ⇒ **the SNR-dominant harmonic is always `n = 1–2`,
   for every `e`.** The detection band *is* the fundamental. Physically-anchored widest lane
   **1.85e-3–3.7e-3** — i.e. essentially the F2 rung itself, 13–27× short.
3. `e ≈ 0.5` maximises usable lane *count* (5 rungs, min per-lane budget 2.9 nat at SNR_tot 15 —
   SNR is *not* the binding constraint), but every rung lands in the **already-populated fine band**
   (≤ 1.85e-3), not the empty coarse gap. **There is an SNR-optimal `e` (~0.5); there is no
   bridge-optimal `e`.**

**LANES verdict, restated: eccentricity buys FINER rungs (precision), not WIDER rungs (capture).**

**Why that is no longer fatal to the E-track — and what changed.** LANES was measuring against the
*sky* obstruction: a 27× gap between F2's loosest sky rung (1.85e-3) and the blind float (0.05),
where the failure mode **is** rung *reach*. The pilot's STEP 1A shows the `(f, mc)` obstruction is
categorically different: **the ladder spans internally** (42 rungs, max gap 0.387 dex < the 0.7
criterion) — *there is no gap to bridge*. The obstruction is the **mc box** (47× shrink needed).
So the E-track's target quantity changes from *lane width* to *box shrink*, i.e. from `tol` to
`σ_mc`. LANES' scissors constrains the former and is silent on the latter. §2.3 is the first
measurement of the latter — and it says the box shrink is cancelled by the wander boost.

**Stated explicitly, as the brief requires — the DEAD use vs the LIVE use of high harmonics:**

| | DEAD (LANES' bridging use) | LIVE (merged clock-hand use) |
|---|---|---|
| role of harmonic `n` | a **detection band** you migrate to, to get a wider lane at low `n` | a **clock hand** in the *same* coherent likelihood; never detected alone |
| requires | `SNR_n` individually detectable at high `n` | nothing — `n` enters through `w_n ∝ g/n²`, which decays as `1/n²`, not `SNR_n² ∝ g/n⁴` |
| killed by | scissors: `n_peak` needs `e ≳ 0.85`, but then `n=1` is < 5 % of peak; and residual SNR pins detection at `n = 1–2` regardless | **not** killed by the scissors — `n = 1–2` detection is *assumed*, not fought |
| killed by (merged) | — | the **clock cancellation** (§2.3b), unless the rank-3 structure (§2.4) breaks it |
| what it would buy | reach: a rung above 1.85e-3 | ignition: `σ_mc` shrink toward 47× |

**LANES Part 4 (evolution envelope) folds in as the certification ENDPOINT, and it is now binding.**
`K_eff ≈ pi/Delta_phi_cycle` (resolvable fundamental fringes before Earth↔pulsar de-coherence). Over
`tau_p = 3 kyr` the fundamental lane collapses (`K_eff → 1`) **only in the high-Mc / high-f_orb /
high-e corner** (`Mc ≳ 1e9 Msun`, `f_orb ~ 1e-8 Hz`, `e ≳ 0.7`); for the bulk (`Mc ≤ 1e8` or
`f_orb ≤ 3e-9`) `K_eff >> 1`. Two consequences under the merged reading:

- **Harmonic `n` de-coheres `n×` faster than the fundamental** (`Delta_phi_cycle → n Delta_phi_cycle`),
  so `K_eff,n ≈ K_eff,1/n`. The number of *usable clock hands* is bounded by
  `n_coh ≈ K_eff,1(Mc, e, f_orb)`. [DERIVED-WEAVE]
- Truncating the Fisher weight at `n ≤ n_coh` costs real information. **[DERIVED-WEAVE, numpy]**
  fraction of `W = sum_n g/n²` retained: at `e = 0.5`, `n_coh = 3` retains 54 %, `n_coh = 5` retains
  87 %; at `e = 0.7`, `n_coh = 5` retains 38 %, `n_coh = 10` retains 81 %; at `e = 0.9`, `n_coh = 10`
  retains only 13 %, `n_coh = 50` retains 78 %.

**The pincer.** Large `e` gives many hands (`n_peak = 51` at `e = 0.9`) but exactly the corner where
`K_eff → 1` kills them; small `e` keeps them coherent but there are `≈ 2` of them. The `(e, Mc, f_orb)`
map is therefore a **compact region problem**, not a monotone limit — which is why §3 must be a map
and not a scan.

---

## 3. THE MAP TO RUN — the revised E-track deliverable

**Object.** A map over `(e, Mc, f_orb)` of the *chirp cascade*, with the circular case as the
`e → 0` boundary. Three measured quantities per grid point, each with the §1.5 wrong-fix column.

### 3.1 The three quantities

**(M-i) CASCADE IGNITION.** Does the loosest mc rung stay reachable? Report `max_a R_a` and the
pulsar attaining it, at each conditioning tier A / B / C (§1.4), with `R_a` generalised to the comb
via the Fisher-weighted pulsar-term phase wander over `(Mc, e, f_orb)`:
`R_a = 1/sqrt( sum_{ij} sig_ij g_i,a g_j,a )`, `i,j ∈ {log10_mc, e, log10_f_orb}`.
Reference values to beat: circular tier A `R = 0.0214`, tier C `R = 0.271`; **ignition = `R ≥ 1`.**
*Pre-registered prediction (§2.3): `R` DECREASES monotonically with `e`, floor `0.816 × R_circ` for
loud0 and `∝ 1/κ` for loud1/2. A measured `R(e) > R(0)` anywhere falsifies the scalar-`κ` model and
localises the rank-3 break (§2.4) — that is the map's discovery channel.*

**(M-ii) CASCADE GAIN PER RUNG vs `e`.** Reproduce the STEP-1B loop-gain table (`MaskedDelay`,
`H = -Fisher` exactly at Asimov truth) as a function of `e`: gain at earth-only / top-1 / top-3 /
top-6 / top-12 / top-24 / top-48 / census-3. Report the **gain-crossing rung** (circular: between 24
and 48). *Prediction: `κ` cancels in the gain ratio too* (`gain = sqrt(sum_{p∈S} SNR_p² g_p²)/g_a`),
*so the crossing rung is `e`-invariant at leading order.* A measured shift is the same discovery
channel as (M-i), seen from the loop side. **This is the "do harmonics accelerate the Mc bootstrap?"
question, and it is answerable with no optimiser.**

**(M-iii) CERTIFICATION ENDPOINT.** `K_eff` at cascade convergence — LANES' Part-4 quantity, but
evaluated **at the registered state**, not at the float, and per harmonic (`K_eff,n ≈ K_eff,1/n`).
Overlay LANES' `K_eff = 1` and `K_eff = 3` contours (`lanes_kcontour.png`,
`Mc ∈ {1e8, 1e9, 5e9}`) with the `n_coh` truncation of `W(e)` from §2.5. Report `n_coh`, retained
`W` fraction, and the resulting `sigma_mc` at convergence. **The endpoint is where eccentricity is
expected to pay** (LANES: finer rungs, precision, the L2c pull-in end).

### 3.2 Grid

- `e ∈ {0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9}` — dense through 0.5–0.8 where `κ` crosses 10–100 and
  where LANES' lane count peaks (`e ≈ 0.5–0.65`) and collapses (`e ≳ 0.7`).
- `log10 Mc ∈ {8.0, 8.5, 9.0, 9.5}` — spans the generative prior `U(8.5, 9.5)` plus the LANES
  `K_eff` contour set.
- `f_orb ∈ {1e-9, 3e-9, 1e-8}` Hz — LANES' Part-4 axis; brackets the `K_eff → 1` corner and the
  pilot's implied `f_gw ≈ 1.95e-8 Hz` ⇒ `f_orb ≈ 1e-8 Hz` circular.
- Array, priors, seeds: the pop draw (N_CW = 16, 116 psr), `best_distances.txt` + LIT_INJECT,
  D4 seeds. Zero-noise Asimov, truth-placed, **LABELLED** (per the spec's language rule).

### 3.3 Machinery reuse — no re-derivation

**Toy tier — `MultiSourceDelay` harmonic stacks.** `lnL_GW_recovery_phase_connected.py:384`.
Its `__call__` sums `n_sources` copies of `make_phase_connected_binary` **sharing one
`{psr}_cw_p_dist`** (line 425) — i.e. one pulsar distance, `N` frequencies. **That is exactly a
harmonic stack.** Two reparametrisations are required and are the gate:

- **Amplitudes (Peters–Mathews):** `log10_h_n = log10_h + 0.5 log10 g(n,e)`.
- **Chirp tying (non-obvious — the stack's `log10_mc_n` must be *offset*, not shared).** The kernel
  is a *circular* chirper: at `f_n` it produces `fdot = (96/5) pi^{8/3} M_n^{5/3} f_n^{11/3}`. To make
  harmonic `n` chirp at the physical `n · fdot_orb` requires **[DERIVED-WEAVE]**
  ```
  log10_mc_n  =  log10_mc + (8/5) log10 2 + (3/5) log10 F(e) - (8/5) log10 n
  ```
  (check: `n = 2`, `e = 0` ⇒ `log10_mc_n = log10_mc` ✓ — the circular case is recovered exactly).
  **A naive shared `log10_mc` across the stack is WRONG and will silently mis-time the comb.**
- **Gate T1:** stack with `n ∈ {2}`, `e = 0` must reproduce `fisher_logL` on the circular pop draw to
  **0.00e+00** (the split's standard, spec §6.2 gate 0).
- **Gate T2:** the stack's `d Phi_n/d log10 Mc` (autodiff) must reproduce `n × ` the fundamental's to
  the leading-order accuracy of §2.2 (≈ 8 %) and reproduce the pilot's M1 mc ladder at `e = 0`
  to 3 digits (the pilot's own reproduction standard vs F2).
- **Limitation, stated:** the tying matches `fdot` only. `e(t)`, `gamma_dot`, and `fddot` are **not**
  reproduced — the toy tier is blind to precisely the rank-3 structure of §2.4, so **the toy tier can
  confirm the cancellation but cannot refute it.**

**Credible tier — EOB.** `arXiv:2511.19611` (Manzini & Babak, 2025-11): maps the rich harmonic
structure and **treats the pulsar term**, but "cannot determine the pulsar phase given poor knowledge
of pulsar distances" and does not link comb → distance (LANES literature review; the comb-as-vernier
angle is **novel / unclaimed**). Use its waveform for the `e(t)`, `gamma_dot`, `fddot` structure that
the toy tier omits. **The rank-3 measurement (M-i / M-ii's discovery channel) is only meaningful at
this tier.**

**Registration / Fisher machinery — all banked and gated, reuse verbatim:**
`MaskedDelay` (`trackB_b1_core.py:92`) for the loop-gain masked likelihood;
`trackB_split.Split` (47× E-step, 37 µs/candidate scorer, gates 0–4 all PASS);
`trackB_p1_map.batched_scan` as the **correctness oracle** (Amendment 6: no shared code);
`stagec_fisher.build_fisher_amortised`; `stagec_hessian` HVP; `stagec_anchor_a2.scan_all` /
`classify_full`. `B_CERT = 512`, fringe-centre evaluation, **never** the prominence peak finder.
XLA autotune OFF, jax persistent cache, x64, jug-gpu.
**Never `pinv`** for the Fisher inverse — it reports `σ = 0` for an unconstrained axis (M3's note);
invert **with** the generative priors.

### 3.4 Gates (pre-registered, in order; stop at the first failure)

| # | gate | pass condition | consequence of failure |
|---|---|---|---|
| **G0** | toy-stack fidelity | T1 = 0.00e+00; T2 ≤ 8 % / 3 digits | stack is mis-tied; fix before any map point |
| **G1** | circular boundary | the `e → 0` map column reproduces STEP 1A/1B: `max R = 2.143e-2`, gain 0.04 → 1.34 (top-48) → 4.50 (census-3) | the map is not measuring the pilot's object |
| **G2** | cancellation test (toy) | `R(e)/R(0)` matches the §2.3 prediction (0.816 floor for gain-1.73 sources; `1/κ` for gain-1.00) to 20 % | scalar-`κ` model wrong ⇒ **go to G3 immediately, this is the discovery** |
| **G3** | rank-3 break (EOB) | 3×3 evolution Fisher condition number and `max_a R_a` at tiers A/B/C | if `max R < 1` everywhere: E-track ignition **NO-GO**, report the map as a pricing instrument (§4) |
| **G4** | ignition | `max_a R_a ≥ 1` at some `(e, Mc, f_orb)` and some tier, **with wrong-fix count 0** | ignition with wrong fixes ≠ ignition |
| **G5** | endpoint | at any ignited point, `K_eff ≥ 3` for `n ≤ n_coh` and census P_true within 0.05 of the A2 ceilings (0.953/0.989/0.998) | ignition without certification |

**Cost discipline.** M-i and M-ii are **Fisher-only** — one masked-likelihood call per subset, no
optimiser, hence **no L2c pull-in problem** and no soft-EM (which is numerically CHAOTIC at
cond ~4e10 — spec's standing warning; if any optimiser is ever needed, replace Adam with a
Newton/trust-region step or profile the extrinsics). The whole toy map is a Fisher sweep.
**No silent caps:** if a grid point is dropped (compile OOM under GPU co-tenancy), log it.

### 3.5 What the STEP-2 tier referendum contributes — and its status

**[UNSOURCED-IN-SPEC / IN-FLIGHT — read from artefacts, not from a written-up section.]** The spec's
STEP-1 verdict says "STEP 2 (the (f,mc) referendum at tiers A/B/C) decides whether any tier is
legitimate conditioning for B1.0–1.4"; pre-registered threshold **`f ≥ 0.95`**. Artefacts present:

| tier | ln Z_needle (quad) | ln Z_box | **f** | verdict vs `f ≥ 0.95` | break-even |
|---|---|---|---|---|---|
| A (sky only) | 405631.862 | 405630.660 ± 1.684 | **0.769** | FAIL | shrink every active dim **1.222×**; mc box half → **0.611 dex** |
| C (+ host D_L) | 405630.368 | 405632.607 ± 0.504 | **0.096** | FAIL | `breakeven_lambda` 0.689 |

**Handle with care.** (1) These are **not** in `trackB_estimator_spec.md`; they are
`b1_referendum_tierA.npz` / `b1_referendum_tierC.npz` + `logs/b1_ref_tier{A,C}.log`, uncommitted.
(2) `logs/b1_ref_tierC.log` is **modified in the working tree and ends mid-quadrature** — tier C is
**PROVISIONAL**. (3) Tier A's `lnZn_lap` is **NaN** and `npos = 5` of 6 (one non-positive curvature
direction); the log records "a MICRO-DIP at truth inside a globally sharp peak … Z_needle is still a
well-defined local integral, but it is NOT a Laplace integral — quadrature only." (4) Tier C's `f` is
*lower* than tier A's despite a tighter box; that inversion is unexplained and is a reason to treat
both as provisional. **No tier B artefact found.**

**Why the E-track cares.** Tier A's break-even is an **mc-box** number: `1.222×` shrink, i.e. mc box
half-width `1.0 → 0.611` dex. §2.3 says an eccentric source at `e ≳ 0.6` supplies `κ ≥ 12.6` of
Earth-term chirp measurement — **an order of magnitude more shrink than tier A's break-even needs.**
So the merged prediction is a **split verdict**, and it is the map's headline:

> **Eccentricity is predicted to FLIP the tier-A (f, mc) referendum (a volume/inference statement) while
> LEAVING cascade ignition unchanged (a lever-arm/cancellation statement).** These are different
> questions about the same box, and the docs have kept them separate: `R` divides the box by the
> *pulsar-term* gradient; `f` compares the box against the *needle's* volume. `κ` cancels in the
> first and does not in the second.

Verifying that split is gate **G6** (add to §3.4): re-run the tier-A referendum machinery
(`trackB_b1_referendum.py`) with the mc box replaced by the eccentric Earth-term Fisher width
`σ_mc/κ(e)`; **pass = `f ≥ 0.95` at the smallest `e` on the grid that achieves it.** Report that `e`.

---

## 4. WHAT IT PRICES — the frontier framing

*Stated as functions of capability, never as verdicts. Nothing here is a claim about what the sky
contains; everything is a claim about what an array would be able to do with a given source.*

### 4.1 The self-clocking corner

A source is **self-clocking** when its own Earth term measures its chirp rate — when
`Delta_phi_E = pi fdot_gw T² ≳ 1 rad`. Circular, at the current loudest source: **0.05 rad** (M3), and
the Earth-term `log10_mc` information gain is **1.00–1.73×** — for two of the three loud sources the
posterior *is* the prior. The array is blind to `fdot`. Eccentricity multiplies that phase by `κ(e)`:

| `e` | `κ` | `Delta_phi_E` at loud0's `fdot` | Earth term measures `Mc`? | tier-A referendum (break-even 1.222×) |
|---|---|---|---|---|
| 0 (circular) | 1.0 | 0.05 rad | no (gain 1.00–1.73×) | `f = 0.769` FAIL |
| 0.5 | 5.6 | 0.28 rad | marginal | passes break-even (5.6× > 1.22×) — **predicted PASS** |
| 0.65 | ≈ 20 | ≈ 1 rad | **yes — threshold** | predicted PASS |
| 0.7 | 36.9 | 1.8 rad | yes, comfortably | predicted PASS |
| 0.8 | 172 | 8.6 rad | yes | predicted PASS; but `K_eff → 1` corner encroaches (LANES Part 4) |

**The corner: `e ≳ 0.65`, at the `(Mc, f_orb)` of the loud sources** (`Mc ~ 1e9 Msun`,
`f_orb ~ 1e-8 Hz`) — *and bounded above* by the LANES `K_eff → 1` collapse in the same
high-`Mc`/high-`f_orb`/`e ≳ 0.7` corner, which removes the clock hands it grants. **The map's payoff
region is a band, not a limit**, and its existence is the map's first deliverable. `Delta_phi_E ∝ T²`,
so the band's lower edge in `e` recedes as `T^{-2}·(SNR)`: a 2× longer baseline moves the
self-clocking threshold from `κ ≈ 20` to `κ ≈ 5`, i.e. from `e ≈ 0.65` to `e ≈ 0.5` (**[DERIVED-WEAVE]**,
same `Delta_phi_E` law).

### 4.2 What self-clocking buys, priced against the banked walls

| capability | circular (measured) | eccentric self-clocking (`e ≳ 0.65`, predicted) |
|---|---|---|
| `sigma(log10_mc)` from the source's own Earth term | 0.501 / 0.864 / 0.866 dex (≈ the prior) | shrink by `κ ≈ 20–37×` ⇒ ~0.02–0.04 dex |
| the mc box the targeted pipeline must search | 1.0 scaled (0.5 dex half-width, generative prior) | Fisher-width; **below tier A's 0.611 dex break-even** |
| tier-A `(f, mc)` referendum | `f = 0.769` (FAIL, threshold 0.95) [provisional] | predicted **PASS** — the targeted posterior concentrates |
| what conditioning tier B1 may legitimately assume | none yet established (STEP 2 all-FAIL so far) | **tier A (sky only)** — i.e. an EM counterpart's sky, nothing more |
| the external distance tier C must buy | `sigma(log10 D_L) = 0.005`, host galaxy | **not needed for the mc box** — the source supplies its own |
| cascade ignition (`max R ≥ 1`) | `R = 0.0214` (A) / `0.271` (C); needs 47× | **predicted unchanged** — clock cancellation (§2.3b) |
| what still gates certification | the mc box **and** the sky | **the sky alone** (lever iii) — and `tau_a` (lever i) |

**The exchange this prices.** Track B's terminal wall is a *sky* wall (F2: loosest sky rung
1.85e-3 scaled vs blind float 0.05–0.21; 1.4–2.05 dex, **zero rungs**), and deliverable R made it
information-theoretic (`f = 6.9e-7`, break-even sky prior **0.188 deg per loud source**, 32× tighter
than the F-stat blind floor ≈ 6 deg). Lever (iii) — an EM counterpart — supplies exactly that sky.
But the pilot then showed that removing the sky leaves **two** registration axes, and the mc box is
1730× too wide. **A circular first source therefore needs the counterpart's sky AND its host
distance, and even then misses ignition by 3.7×.** An eccentric first source with `e ≳ 0.65` needs
**only the sky**: it clocks itself. That is the price difference the map quotes — and, per §2.3, it
buys the *referendum*, not the *cascade*.

### 4.3 The three-lever ledger, re-drawn

The design theorem's levers, with the merger folded in:

- **(i) Wide lanes from nearby pulsars.** `tol ∝ 1/L_p`. J0711-6830 (0.106 kpc, 1.85e-3),
  J1630+3734 (0.089 kpc, 1.39e-3), J0437-4715 (0.155 kpc, 1.02e-3) vs median `L = 1.38` kpc.
  **Unchanged, and now the *only* ignition lever**: §2.3's ceiling `R_a^max = C·SNR·(T/tau_a)²` is a
  statement purely about `tau_a`, `T`, `SNR` — and `tau_a` is what a ~100 pc, well-timed target
  sub-array buys. `R_a^max = 1` needs `tau_a ≲ 85 yr` (SNR 15) to `374 yr` (SNR 289).
- **(ii) Eccentric harmonics.** **Re-scoped by this document.** *Not* a lane-widener (LANES, and the
  merged §2.3 confirms the sign). It is the **evolution/timestamp channel made measurable**: it
  supplies the `Mc` the pulsar term needs in order to be read, from the source's own Earth term.
  It prices out the mc axis; it does not move the lever arm.
- **(iii) Lowering the wall from above.** An EM counterpart removes the sky axis; `T` and `SNR`
  lower the float (R5 pencil: the Earth-term blind float reaches the F2 loosest rung only at
  `T ≈ 11 kyr`, or `SNR ≈ 289` at fixed `T`; the L2c pull-in at `T ≈ 3.75 Myr`).

> **Certification is feasible where lever (i)/(ii) raises a lane above the lever-(iii)-lowered float.**
> The merger's contribution to that sentence: **lever (ii) does not raise the lane — it removes the
> `Mc` that made the lane unreadable.** After (ii), lever (i)'s `tau_a` and lever (iii)'s sky are the
> whole remaining ledger, and the `R_a^max = C·SNR·(T/tau_a)²` ceiling says how far each must go.

### 4.4 Standing caveats

- The `κ` algebra is **leading-order in `fdot tau/f`** (measured 0.09 at the pilot's median) and
  assumes `Phi_n = n Phi_orb`. **§2.4 lists exactly the three effects that break it**; the map's
  discovery channel is their measurement, at the EOB tier.
- Every `R`, `κ`, and `Delta_phi_E` number here is **Fisher-scaling**, not likelihood-measured.
  The pilot's own `R` values were measured; mine are predictions **against** them.
- `SNR_tot = 15` is LANES' toy value (A5), not this array's measured per-source SNR.
- LANES A4 (`SNR² ∝ g/n⁴`, white noise) is **optimistic for high `n`**: real red noise / GWB rising
  to low `f` pushes the usable band lower, which **strengthens** the `n_eff ≈ 2–4` result used here.
- Zero-noise Asimov, truth-placed, throughout. Per the spec's language rule this is **LABELLED**
  diagnostic truth use; it is a *ceiling*, and noise can only make it worse.

---

## Sourcing ledger

**Verbatim from `trackB_estimator_spec.md`:** §1.1 M1/M2/M3 tables; §1.2 STEP-1A ignition, 47×,
3305×, 42 rungs / 0.387 dex; §1.3 STEP-1B loop-gain table + killed hypothesis + bistability;
§1.4 STEP-1C tiers; §1.5 L2 wrong-integer margins + R-postmortem table + §2a partial-fix rule;
§2.5 design theorem levers; §4 R referendum `f = 6.9e-7`, 0.188 deg, R5 pencil (11 kyr / SNR 289 /
3.75 Myr); F2 sky ladder; L2b/L2c pull-in and eigenvalues; §3.3 machinery and gates.
**Verbatim from `LANES_eccentric_ladder.md`:** scissors (both blades), Parts 2/3/4 tables, `K_eff`
corner, approximations A1–A6, EOB `arXiv:2511.19611` and the novelty finding.
**Verbatim from `project_progress.md`:1780–1827:** the STEP-1 summary and the merger sentence.

**[DERIVED-WEAVE]** (my algebra + `numpy`/`scipy`, all reproducible): the `d Phi_n/d log10 Mc` law
and its 1.084× agreement with the pilot's mc median; the 0.5-dex-per-scaled-unit conversion;
`sum_n g = F(e)` gate; `κ(e)`, `n_eff`, `W(e)` and the `n_coh` truncation table; the residual-SNR
flatness (0.85–1.00); the **clock-cancellation ceiling** `R_a^max = C·SNR·(T/tau_a)²` and the
`R(e)/R(0)` table; the retarded-epoch sign; the ≤2× vernier ceiling on `K`; the toy-stack
`log10_mc_n` tying formula; the `T^{-2}` recession of the self-clocking edge.

**[UNSOURCED]** and flagged in place: `PROJECT_PROGRESS.md` (file does not exist; used
`project_progress.md`); J0711's `tau_p ≈ 0.35 kyr` (estimated from `L`, not read from
`b1_pilot_m1.npz`); the phrase "wrong-fix column" (the brief's, not the docs'); the STEP-2 tier
referendum numbers (artefacts only, tier C mid-run, tier A Laplace NaN, no tier B, `f_C < f_A`
inversion unexplained); the constant `C` in `R_a^max` (an `O(1)` covariance penalty, not evaluated).

*(Agent WEAVE. Read-only; this file is the sole output. No computation beyond numpy-scale checks.)*
