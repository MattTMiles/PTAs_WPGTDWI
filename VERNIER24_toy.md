# VERNIER-24 — the clean-band wheel on the audited sky

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~5 min.
Clean room: imports only `vernier_toy` / `vernier6` / `vernier7` / `vernier16`.

Code `hpc_harbor/vernier/vernier24.py` · data `reports/vernier24.npz`.
Staged, not committed.

**Re-aimed mid-run** on the CENSUS_CLASH audit. The one-source-per-bin sky is a known
artefact and is **not reported as a bracket** — it is kept only as a control.

---

## Lead — realized separation collapses to the Earth-term value, and that closes VERNIER-23

**The counting bound is not reached, and the marginalised pulsar term buys nothing that
survives real competition.**

Realized frontier — largest M at which at least **half the bin** is recovered with errors
inside tolerance (σ_sky < half the typical separation **and** σ(A)/A < 0.5), distances
marginalised over their 10 pc priors:

| f | **half-bin M_max, marginalised** | **half-bin M_max, Earth-only** | top-source M_max | counting bound 3N/7 |
|---|---|---|---|---|
| 10 nHz | **5** | **5** | 17 | 50 |
| 15 nHz | **10** | **5** | 25 | 50 |
| 20 nHz | **10** | **10** | 25 | 50 |

*(N = 116 pulsars, ρ₁ fixed at 30 so loudness is not the limiter. Earth-only counting bound
2N/6 = 39; distances-known 4N/7 = 66.)*

**Realized marginalised M_max is 5–10 against a counting bound of 50 — a factor 5–10 short —
and it sits on top of the Earth-only realized value of 5–10.** The apparent 2× at 15 nHz is
one grid step on four realisations and is inside the scatter. **VERNIER-23's M_max = 17–20
was a counting bound that happened to coincide with its measured value at N = 40; at N = 116
against a real ladder the measured value does not follow the bound, and the marginalisation
gain the bound predicts is not realized.** That is the honest close on VERNIER-23.

**And at the census's own loudness it is worse: the half-bin frontier is 0 at every
frequency** — no bin in the band has half its sources recoverable, because ρ₁ = 4.4 / 2.4 /
2.1 at 10 / 15 / 20 nHz. Only the loudest source survives, and only up to 10 nHz.

| | 10 nHz | 15 nHz | 20 nHz |
|---|---|---|---|
| top-source M_max, census loudness | 17 | 10 | **0** |
| half-bin M_max, census loudness | **0** | **0** | **0** |

**Set against the audit's N_eff ≈ 8: the bins are marginally winnable on identifiability and
not winnable at all on sensitivity.** The binding constraint in-band is loudness, not
separability — the opposite of what VERNIER-23's Part 2 concluded on the artefact census, and
not the reason it gave.

The conditioning tells the same story without reference to any SNR: the smallest normalised
eigenvalue of the marginalised source block falls 7×10⁻⁴ → 3×10⁻⁵ → 3×10⁻⁷ → 2×10⁻⁸ →
3×10⁻¹¹ → 7×10⁻¹² across M = 2 → 50 at 10 nHz. **There is no knee; it degrades continuously,
which is why a rank-based M_max overstates what is recoverable.**

## The audited sky, rebuilt and checked

| f | N_bin (built) | N_bin (audit) | λ_fore (built) | λ_fore (audit) | N_eff | h₁ | implied log₁₀ Mc | ρ_array | ρ_P p90 |
|---|---|---|---|---|---|---|---|---|---|
| 10 nHz | 6.86×10⁴ | 6.9×10⁴ | **0.275** | 0.31 | 9.8 | 1.08e-15 | 9.36 | **4.43** | 0.86 |
| 15 nHz | 1.55×10⁴ | 1.6×10⁴ | **0.308** | 0.31 | 7.8 | 7.12e-16 | 9.18 | **2.43** | 0.44 |
| 20 nHz | 5.41×10³ | 5.4×10³ | **0.297** | 0.31 | 7.4 | 5.00e-16 | 9.03 | **2.09** | 0.32 |

Counts reproduce the audit exactly. **λ_fore comes out flat in N_bin across a 13× range in
occupancy — 0.275 at 6.9×10⁴ and 0.297 at 5.4×10³ — which is the structural
γ_eff = 1.56 < 2 behaviour the audit identified, and it is why no one-source bin exists at any
frequency.** The h(Mc, f, D) link gives a measured h²/bin slope of **−2.3333 against the
required −7/3**, so the census's γ_eff = 5.51 (h drawn independently of f) is not reproduced.

**One place the reconstruction only half-matches: N_eff.** The participation ratio
1/Σp² gives 7–10, the **lower** end of the audited 8–180. If the operative N_eff is nearer
180 the frontier comparison above gets worse, not better.

## Does the wheel turn on the corrected sky? It ignites, and it does not turn

| arm | turns | genuine | blend | ρ_array | ρ_P med | ρ_P p90 | **CERT** | true | scrambled null | depth | rank |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **B — physical sky** | 8.0 | **1.2** | 6.5 | 6.03 | 0.243 | 1.00 | **0** | 0 | 0 | +0.0017 | **5** |
| Am — zero-competitor control | 3.0 | 0.5 | 2.8 | 3.12 | 0.245 | 0.57 | **0** | 0 | 0 | +0.0006 | 4 |
| A — superseded artefact sky | 7.0 | 5.0 | 2.0 | 6.95 | 0.337 | 1.29 | **0** | 0 | 0 | +0.0144 | 5 |

**On the physical sky the wheel detects — but 6.5 of its 7.7 recoveries per run are blends,
and it certifies nothing.** Depth gained is +0.002 dex, rank 5.

**The one-source artefact sky is the only arm where the wheel looks like it works** (5 genuine
finds per run) — which is exactly why it must not be reported as a bracket.

### Why certification is zero, and it is not the ρ_P ≥ 2 bar

**41 % of Arm B's pulsar-term power sits in pulsars above VERNIER-1's ρ_P ≥ 2 bar, and zero
pulsars certify.** The bar was never the binding constraint here:

- **A single source cannot certify at any SNR.** Every fringe centre of one comb gives an
  identical template, so the fringe contrast is exactly zero up to the chirp envelope. The
  M = 1 control confirms it: detected 4/4, ρ_P max **5.46**, best fringe contrast **0.64 nat**,
  **0 certifications**.
- **What can certify is commensurability between several combs.** Accumulating sources across
  the band does build contrast — 0.00 → 0.01 → 0.17 → 0.17 → 0.18 nat over five turns in
  Arm A — **against a trials bar of ln K ≈ 4.8. It is 25× short after five sources.**

Better timing does not close it either:

| σ_toa | ρ_array | ρ_P p90 | best fringe contrast | bar (ln K) | CERT |
|---|---|---|---|---|---|
| 100 ns | 5.37 | 0.94 | 0.041 | 4.63 | 0 |
| 30 ns | 9.57 | 1.84 | 0.101 | 4.79 | 0 |
| **10 ns** | 10.71 | 1.81 | **1.015** | 4.98 | **0** |

### A methodological correction worth keeping

**The certification floor has to be the incoherent-comb null, not the no-comb null.** Built
the naive way — nulls with the coherent term removed — certification fires while the
source-level scrambled null fires at ~70 % of the same rate, because randomly offset combs
align by chance and the naive floor cannot see it. Rebuilt against nulls that carry the same
combs at random relative phase, the two separate cleanly (at ×100 loudness: 22 certifications,
**22 of them on the correct fringe**, against 3 in the scrambled null; the naive floor would
have passed 42). **The naive floor over-certifies by 1.9–3.7×, and every certification the
surplus buys is a chance alignment.**

## Loudness — what the bright end would have to weigh

| arm | gain | ρ_array | ρ_P med | ρ_P p90 | genuine | CERT | true | null | above ρ_P ≥ 2 | depth | rank |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **B** | ×1 | 8.80 | 0.42 | 1.49 | 1.0 | **0** | 0 | 0 | 0.41 | +0.000 | 4 |
| B | ×10 | 13.82 | 0.56 | 2.02 | 1.7 | 0 | 0 | 0 | 0.84 | +0.009 | 5 |
| B | ×100 | 26.99 | 1.17 | 4.60 | 1.0 | **0** | 0 | 0 | 0.97 | +0.017 | 5 |
| Am | ×30 | 24.00 | 0.87 | 4.01 | 7.0 | **1** | 1 | 0 | 0.92 | +0.115 | 7 |
| Am | ×100 | 44.84 | 1.51 | 6.39 | 7.0 | **3** | 3 | 0 | 0.98 | +0.380 | **18** |

**On the physical sky certification never fires, at any loudness up to ×100.** It fires only
in the zero-competitor control, and only from ×30 — where the bright end already holds 9× the
bin's measured GWB power and the run has stopped answering a physical question.

Loudness needed for ρ_P = 2 at the bin's own confusion, from the measured scaling:

| arm | median pulsar | p90 pulsar | best pulsar |
|---|---|---|---|
| **B** | ×1600 (+1.60 dex) | **×4.6 (+0.33 dex)** | already there |
| Am | ×154 (+1.09 dex) | ×10 (+0.50 dex) | ×2.3 (+0.18 dex) |

**The p90 pulsar needs +0.33 dex of bright-end amplitude to clear ρ_P = 2** — inside the
audit's own 0.25–0.67 dex uncertainty on the loudest in-bin strain. **So the ρ_P ≥ 2 bar is
genuinely marginal in-band. It is also not the constraint that matters**, per the section
above.

## VERNIER-23 Part 3, premise restored and run

Marginalised Fisher rank at the real array size (N = 116, σ_L = 10 pc, 14.75 nHz) is full to
M = 40–60 (rank/params 0.99–1.00), consistent with 3N/7 = 50 — **the counting bound is
confirmed at the real array size, and deliverable 1 shows it is not what is achieved.**

| | turns | runs finding anything | f_cert / f_marg | pulsar-term depth | total depth | rank |
|---|---|---|---|---|---|---|
| hard certification (v2.2) | 8.0 | 1/3 | 0.0000 | +0.0000 | +0.0059 | **5** |
| **marginalised pulsar terms** | 8.0 | 1/3 | **0.2053** | **+0.0559 dex** | **+0.0503** | **5** |

**Dropping the certification threshold entirely does exactly what VERNIER-23 predicted — the
marginalised pulsar term recovers f_eff = 0.21 where hard certification recovers 0.00, worth
+0.056 dex — and it does not move the rank.** Rank 5 against the target of 11, because the
chain never gets past its first genuine detection: only 1 run in 3 finds a genuine source at
all, and depth is dominated by what the search can subtract, not by what the pulsar terms add.

## The spin-off survives the re-cut

σ(log₁₀ Mc) from the pulsar-term line width, at the loudest in-bin source:

| f | line width × T | resolvable | **coherent competitors** | competitors as Gaussian | isolated (white only) | coh / iso |
|---|---|---|---|---|---|---|
| 10 nHz | 1.36 | yes | **0.0108** | 0.0210 | 0.0055 | 1.97 |
| 15 nHz | 6.41 | yes | **0.0054** | 0.0069 | 0.0034 | 1.60 |
| 20 nHz | 18.92 | yes | **0.0024** | 0.0033 | 0.0022 | 1.07 |

**The line is resolvable at all three frequencies (1.4–19 Fourier bins), and competition costs
only 1.1–2.0× — so 0.002–0.011 dex stands with N_eff ≈ 8 competitors present.** Against the
banked chirp-mass posterior at 0.50–0.87 dex that is still a 50–400× improvement, and it is
the most robust result in this instalment.

**Note the middle column: charging the competitors as an unresolved Gaussian sea is 1.3–1.9×
*worse* than modelling them.** That is the λ_eff suppression of VERNIER-21 doing real work —
a modelled interloper only injures the target through its template overlap.

---

## What this does not settle

1. **3–4 realisations per cell, and the ladder is heavy-tailed.** With γ_eff = 1.56 the drawn
   λ_fore ranges roughly 0.1–0.8 between realisations, so ρ₁ moves by ~2× draw to draw. The
   frontier numbers are also quantised to the M grid (2, 5, 10, 17, 25, 50). **Read 5 vs 10 as
   "one grid step", not as a measured difference.**
2. **Fisher errors are local.** σ_sky ≈ 10⁻⁴ rad is a within-fringe width, not a global
   uncertainty; the pulsar-term posterior is multi-modal in exactly the way the certification
   sections measure. The frontier is therefore an upper bound on realized separation.
3. **N_eff reproduces only the audit's lower end** (7–10 against 8–180). The whole frontier-vs-
   competition comparison is soft in that ratio.
4. **The search is Earth-term only.** A pulsar-term search is computationally unreachable
   (VERNIER-16: 5.3×10¹⁰ grid points), so the λ_eff suppression is measured on Earth-term
   templates, as in VERNIER-21/22.
5. **I regularised the 1/(1 + cos μ) antenna divergence at the series-wide 5×10⁻² clip.**
   VERNIER-16's `antenna_grid` clips at 10⁻⁶, which is harmless for an F-statistic but let one
   nearly antipodal pulsar carry a response of ~10⁵ and dominate every SNR and Fisher here.
   **Every `ρ_P max` column in this note is therefore ~an order of magnitude below the
   equivalent column in VERNIER-20/22 and is not directly comparable to them.**
6. **The certification statistic and its incoherent-comb floor are new in this instalment**
   and have not been cross-checked against any earlier result. The banked ρ_P ≥ 2 bar and this
   test disagree by construction (41 % of power above the bar, zero certifications); I believe
   the fringe-contrast test is the right one, but that is an argument, not a validation.
7. **Mc = 10⁹ M☉ fixed** for the chirp time throughout; the line-width numbers move with it.
8. **White noise plus a Gaussian sea, one polarisation model, no red noise, no spectral fit of
   the GWB, no eccentricity.**
9. **The loudness sweep leaves the measured GWB above ×1.1 (Arm B)** — those rows answer a
   counterfactual and are labelled with the bin fraction they imply.

## One-line verdict

**On the audited sky the deciding test comes back negative: realized half-bin separation with
the pulsar distances marginalised is M_max = 5–10 against a counting bound of 50, and it sits
on top of the Earth-only realized value of 5–10, so the marginalisation gain that VERNIER-23's
counting argument predicted — and that its N = 40 run appeared to confirm — is not realized
against a real ladder at N = 116; at the census's own loudness the half-bin frontier is 0 at
every frequency in the band, so the bins are marginally winnable on identifiability and not
winnable at all on sensitivity, with loudness the binding constraint rather than separability;
the wheel consequently ignites on the physical sky but does not turn (1.2 genuine recoveries
against 6.5 blends per run, +0.002 dex, rank 5, zero certifications), and certification fails
for a reason that is not the ρ_P ≥ 2 bar — 41 % of the pulsar-term power sits above that bar —
but because one source cannot break its own fringe degeneracy at any SNR and the
cross-frequency vernier reaches only 0.2 nat of fringe contrast against a 4.8 nat trials bar,
still 5× short at 10 ns timing; VERNIER-23's Part 3 runs cleanly now and behaves exactly as it
predicted, recovering f_eff = 0.21 and +0.056 dex where hard certification recovers nothing,
without moving the rank past 5 of a target 11; and the one thing that survives intact is the
spin-off — chirp mass to 0.002–0.011 dex from the line width with N_eff ≈ 8 competitors
present, degraded only 1.1–2.0× from the isolated case, against a banked posterior that is
prior-dominated at 0.50–0.87 dex.**
