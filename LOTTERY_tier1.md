# LOTTERY — tier 1: P(switch-on) across candidate eccentricity distributions

**Agent LOTTERY, ACCRE, 2026-07-17. CPU only (general lanes; no GPU touched).**
Staged, not committed. Banks/scripts listed in the add-list at the foot.

> **The question.** Nature supplies a *mixture* of binary eccentricities, not a single value.
> Given a candidate population distribution of eccentricity, what is the probability that a
> census-loudness mixed population **switches the certified count on**? And does the cheap
> analytic proxy that predicts it — the active-harmonic **channel budget** (S7.6.4a–c / MAGPIE
> J1) — agree with the externally-quoted **threshold rule** (RECUT §2.1: *any single member
> e ≥ 0.5, or any pair e ≥ 0.3*)? **Where they disagree is where tier 2 spends GPU.**

---

## 0. THE MODEL, AND WHY THE CHANNEL MAP IS EXACT (not interpolated)

The census is CHORUS's **3 + 13**: 3 loud members + 13 faint, all at census loudness
(h = −13.25, mc = 9, f_gw = −7.9, T = 30 yr). **Each loud member draws its eccentricity from
the candidate distribution; the 13 faint members are circular.** Two switch-on statistics are
scored on every draw:

**(1) CHANNEL BUDGET** — the *mechanism* isolated by the equal-κ contrast (S7.6.4a):
> `n_active = 13·(1) + Σ_{loud i} c(e_i)`,  **switch-on iff `n_active ≥ 30`.**

`c(e)` — the active-harmonic count of one member — is **not** a 5-point interpolation of the
README. It is computed from the **real ATLAS harmonic algebra**
(`atlas_harmonics.g_n`, the Peters–Mathews weights), active iff `g(n,e) ≥ 10⁻³·maxₙ g`, over the
production stack `n = 1..32`. This reproduces the `chorus_ch_README.md` census mapping and **every
banked `n_active` anchor exactly**:

| e | 0 | 0.15 | 0.20 | 0.25 | 0.28 | 0.30 | 0.35 | 0.45 | 0.50 | 0.55 | 0.70 | 0.80 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **c(e)** | **1** | 5 | 6 | 7 | 8 | **8** | 10 | 14 | **17** | 20 | **32** | **32** |

README wants 8 / 17 / 32 / 32 at e = 0.3 / 0.5 / 0.7 / 0.8 and 1 for circular — **all hit.**
`n_active`: m0 = 16, m1e03 = 23, m2e03 = 30, m3e03 = 37, m1e05 = 32, m1e07 = 47 — **all reproduced
to the integer.** The `N_HARM = 32` cap *is* the e = 0.8 magnitude truncation: `c(e)` saturates at
32 by e ≈ 0.65, so nothing above the ATLAS validity edge is counted.

**(2) THRESHOLD FORM** — the quoted external rule (RECUT §2.1):
> **switch-on iff** (any loud member `e_i ≥ 0.5`) **OR** (`≥ 2` loud members with `e_i ≥ 0.3`).

Both are per-draw booleans; `P(switch-on)` is the fraction over ≥ 2×10⁵ (named) / 5×10⁴ (mix)
draws, RNG seed 20260717.

---

## 1. THE HEADLINE STRUCTURAL RESULT — the threshold rule is a *conservative* proxy

> ### **On every distribution tested, threshold-ON ⟹ channel-ON. The disagreement is one-sided: the channel budget is uniformly the more permissive statistic.**
> Across the entire 21 × 17 mix grid, `P_chan − P_thr ∈ [0.000, 1.000]` — **never negative.**
> There is no cell where the quoted rule fires while the channel budget does not.

The mechanism is transparent: the quoted rule is a **hard cut on individual eccentricities**
(0.5 single, 0.3 pair), while the channel budget is **additive in `c(e)`**. Sub-threshold
eccentricities that individually trip neither cut still deposit channels — `c(0.25) = 7 > 1` — and
they accumulate. A population can therefore reach `n_active ≥ 30` with **no member at 0.5 and no
pair at 0.3.** The threshold rule is blind to exactly that regime; the channel budget is not.

**The blind spot is localised and sharp: `e_char ∈ [~0.17, 0.28)`** — below the 0.3 pair cut, so
the threshold rule reads **identically zero** there, while the channel budget rises with the
eccentric fraction (figure `LOTTERY_tier1_mix.png`, panel C). The maximal disagreement (ΔP = 1.00)
is a **whole population of e = 0.25 members**: `n_active = 13 + 3·7 = 34 ≥ 30` (channel: **ON,
certain**) while no member reaches 0.5 and no pair reaches 0.3 (threshold: **OFF, certain**).

---

## 2. P(switch-on) — the named distributions (N = 2×10⁵, ± binomial sd)

| distribution | family | `P_chan` | `P_thr` | \|disagree\| | ch-not-th | ⟨`n_active`⟩ |
|---|---|---|---|---|---|---|
| circ (δ at e = 0) | a | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 16.00 |
| lnN peak 0.1, w 0.1 | b | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 25.02 |
| lnN peak 0.1, w 0.2 | b | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 25.25 |
| **lnN peak 0.3, w 0.1** | b | **1.0000** | **0.5593** ± 0.0011 | **0.4407** | 0.4407 | 38.59 |
| **lnN peak 0.3, w 0.2** | b | **1.0000** | **0.6244** ± 0.0011 | **0.3756** | 0.3756 | 40.34 |
| lnN peak 0.5, w 0.1 | b | 1.0000 | 1.0000 | 0.0000 | 0.0000 | 65.96 |
| lnN peak 0.5, w 0.2 | b | 1.0000 | 1.0000 | 0.0000 | 0.0000 | 72.15 |
| **unif [0, 0.9]** | c | 0.9672 ± 0.0004 | 0.8889 ± 0.0007 | **0.0783** | 0.0783 | 63.98 |

*Log-normal parameterisation (Taylor/Farr plug-in): `e ~ LogNormal(μ, σ)` with `μ = ln(e_peak) + σ²`
so the **mode** is `e_peak`; `σ` is the log-space width `w`. `μ, σ` are the pluggable shape
parameters; draws clipped to [0, 0.9].*

**Reading it.**
- **The realistic astrophysical case is the disagreement case.** A population whose eccentricities
  *peak at 0.3* (exactly the low-e end nature is expected to populate) is the single largest
  disagreement among the named distributions: the channel budget says **certain switch-on**
  (three members near e = 0.3 → `n_active ≈ 38–40`), while the quoted rule gives only **0.56–0.62**
  because it needs *two* members to individually clear 0.3, and the log-normal scatter routinely
  drops one below. **The proxy and the rule diverge by ~0.4 exactly where the science lives.**
- Peaks at 0.1 are dead on both (⟨`n_active`⟩ ≈ 25 < 30); peaks at 0.5 are ON on both. The two
  statistics **only part company at low-to-moderate e** — precisely the interesting corner.
- Uniform [0, 0.9] switches on almost always on both (a fat tail past 0.5 dominates), with a small
  8 % residual disagreement from the sub-threshold draws.

---

## 3. The mix surface — the observer-facing panel (`LOTTERY_tier1_mix.png`)

Each loud member is eccentric (`e = e_char`) with probability `f_ecc`, else circular; the grid
sweeps `f_ecc ∈ [0,1]` × `e_char ∈ [0.1, 0.9]`, 5×10⁴ draws/cell.

- **Panel A** (channel budget): the `P = 0.5` contour **dips below e_char = 0.3** at high `f_ecc` —
  the channel budget fires on populations of purely sub-threshold members.
- **Panel B** (threshold rule): the contour **flattens at e_char ≈ 0.28** and the whole region
  below is identically zero — the rule cannot fire without an individual crossing.
- **Panel C** (disagreement `P_chan − P_thr`): a bright, sharply-bounded strip at
  `e_char ∈ [0.17, 0.28)`, brightening toward `f_ecc = 1`. **This strip is the tier-2 target.**

---

## 4. THE VALIDATION SET for tier 2 (selected; analytic P below, measured P to be filled)

Selected per the brief's priority order **(i) disagreement > (ii) P ≈ 0.5 > (iii) concordant**.
Because the mix is fixed-`e`, each point collapses to a small set of CHORUS-style cells
`(n_ecc, e, tier=lit, T=30)`; the decisive constituent cell is named so tier 2 refits **one** floor
per cell and scores ≥ 20 signal realisations against it. **`n_active ≥ 30` is the channel
prediction; the threshold-rule verdict and the two grades are what the GPU measures.**

| # | pri | distribution point | decisive cell | `n_active` | channel | threshold | analytic `P_chan` | analytic `P_thr` |
|---|---|---|---|---|---|---|---|---|
| **V1** | i | mix `f_ecc=1.00, e_char=0.25` | **m3 e=0.25** | **34** | **ON** | **OFF** | 1.000 | 0.000 |
| **V2** | i | mix `f_ecc=0.70, e_char=0.20` | **m3 e=0.20** | **31** | **ON** | **OFF** | 0.343 | 0.000 |
| **V3** | i | mix (near-knee) | **m2 e=0.28** | **30** | **ON** | **OFF** | — | 0.000 |
| **V4** | ii | knee, below | **m2 e=0.25** | **28** | OFF | OFF | — | 0.000 |
| **V5** | ii | knee, single-member | **m1 e=0.45** | **29** | OFF | OFF | — | 0.000 |
| **V6** | iii | concordant HIGH | **m3 e=0.70** | **109** | ON | ON | ~1.000 | ~1.000 |
| **V7** | iii | concordant LOW | **m1 e=0.25** | **22** | OFF | OFF | ~0.000 | ~0.000 |
| — | xcheck | banked concordance | m2 e=0.30 (CHORUS m2e03) | 30 | ON | ON | — | — |

**The decisive test is V1–V3:** cells the channel budget calls ON while the threshold rule calls
OFF. If their certified counts **clear the >1 bar at floor + bootstrap error**, the channel budget
is right and the quoted threshold rule is confirmed *too conservative* — and the whole
P(switch-on) surface is licensed at analytic cost. If they **fail**, the proxy is over-permissive
in the sub-threshold strip, and V4/V5 (n_active = 28, 29, just under the knee) bracket where it
turns over. **V6/V7 anchor the extremes.** (m2 e=0.30 = the banked CHORUS `m2e03`, CONFIRMED 2.77,
carried as a free concordance cross-check of the refit path.)

*Full mix candidate rankings and the P ≈ 0.5 concordant-steep cells (`f_ecc=0.50, e_char=0.35`;
`f_ecc=0.20, e_char=0.55`, both `P_chan ≈ P_thr ≈ 0.50`) are in the bank + selection log; the
steep-concordant pair is held as a spillover if the GPU budget allows, testing proxy calibration
where both statistics agree and are steep.*

### Tier-2 budget note (from CHORUS's measured per-realisation cost)
CHORUS stage-1 logs: `build_chorus_problem` ≈ 190–280 s (once per `n_ecc` shape, shared across all
cells at T = 30), first-compile per e-config ≈ 150–226 s (once per cell), **warm realisations
≈ 0.7–1.0 s.** A refit floor (≤ 100 nulls) + 20 signal reals per cell is therefore
**≈ 200 s compile + ≈ 120 s reals ≈ 0.1 GPU-hr/cell** once the three shape-builds (n_ecc 1/2/3,
≈ 0.2 GPU-hr) are paid. **7 cells ⇒ ≈ 1–3 GPU-hr — comfortably inside the ~60 GPU-hr cap**, no trim
needed. *(To be confirmed by the tier-2 smoke test before the wide submit; if warm reals prove
10× costlier the set still fits.)* Weekend general-queue H200/A100, fair-share, preemptible →
checkpoint per draw + resume drill precede the production array. **The reserved 4-card H200 share
and all SPARK-3 jobs are untouched.**

---

## 5. WHAT TIER 1 SETTLES

1. **`P(switch-on)` is computed for every candidate distribution** (named table §2 + the full mix
   surface §3), from a channel map that reproduces the census mapping and every banked `n_active`
   exactly — so the proxy is the campaign's own mechanism, not a new model.
2. **The two switch-on statistics disagree one-sidedly:** the quoted threshold rule is a
   *conservative* proxy; the channel budget fires on sub-threshold populations it misses.
3. **The disagreement is astrophysically located** — at the low-e peak nature is expected to
   populate (log-normal peak 0.3: ΔP ≈ 0.44) and in the mix strip `e_char ∈ [0.17, 0.28)`.
4. **The validation set is selected and budgeted** (§4). Tier 2 measures the certified-count grade
   at the disagreement cells; the analytic-vs-measured calibration curve is the deliverable that
   either licenses the whole surface at analytic cost or maps the proxy's failure direction.

> **TIER 2 EXECUTED (see `LOTTERY_tier2.md`).** The validation set was run on the general-queue
> H200 (0.74 GPU-hr). **Result: the channel-budget proxy does NOT track the certified-count grade.**
> On a common basis every tested mixture from `n_active = 22` to `34` clears the >1 bar, and the
> highest-channel cell (34) is only MARGINAL — because certification is gated by the **refit floor**
> (which rises with eccentric structure), a term the analytic proxy omits. The P(switch-on) surface
> above therefore stands as the *channel-budget statistic it is defined to be*, but is **not** a
> calibrated predictor of the GPU switch. Bonus finding: RECUT's "single e = 0.3 REFUTED" grade
> **flips to CONFIRMED on a second host** (near-bar grades are fragile at ~1 nat of floor error).

**STOP is at the foot of tier 2, not here** — tier 1 hands the validation set forward.
```
add-list (stage, do not commit):
  LOTTERY_tier1.md
  LOTTERY_tier1_mix.png
  reports/lottery_tier1.npz
  hpc_harbor/lottery/lottery_tier1.py
  hpc_harbor/lottery/lottery_tier1_figure.py
```
