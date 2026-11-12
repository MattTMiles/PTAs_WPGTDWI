# VERNIER-5 — the working point

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~2 min.
Clean room: imports only `vernier_toy` / `vernier2` / `vernier4` (my own anchored code).

Code `hpc_harbor/vernier/vernier5.py` · data `reports/vernier5.npz` ·
plot `reports/VERNIER5_workingpoint.png`. Staged, not committed.

---

## Lead

**I reproduce both papers' *arithmetic* exactly and neither paper's *claim* as a distance
measurement.** The within-fringe precision they report is cheap — it needs only
ρ_P ≈ 0.03–0.27 per pulsar. At that same SNR the fringe is not identified: P(unique
distance) ≤ 0.03. **Their sub-parsec numbers are within-fringe σ_L, which the programme's
own D3 established is necessary-but-not-sufficient.**

Per your instruction, that discrepancy is the deliverable and **I stopped before the
four-dial walk.**

**But the discrepancy cuts both ways, and the second half of this report is a correction to
VERNIER-4, not to the papers.** The literature threshold λ_fore ≥ 0.7 is not reachable by
any array — it is a property of the population's **bright-end counts slope γ**, and it needs
**γ ≤ 0.87** against Euclidean γ = 3. My census assumed γ = 3, which is exactly why VERNIER-4
found "never closes". **VERNIER-4's pessimistic column is an artefact of a missing mass
function**, and the knee you flagged at 25–40 nHz is the thing that would void it.

---

## 0. Convention check — no √(fT) error, and the numbers move *down*

VERNIER-4 built **both** sides of the ratio from per-bin residual power:
`Σ_k R_k² = 2 S_gwb(f)·Δf` and `R_1² = that/(3 N_c^(1/3))`. `h_c` never entered, so
`h_1/h_sea = R_1/√(ΣR_k²)` is like-with-like.

| f nHz | as built | via h_c | factor | √(fT) | λ_fore |
|---|---|---|---|---|---|
| 10 | 0.0902 | 0.0737 | **0.8165** | 2.18 | 0.0081 |
| 32 | 0.1837 | 0.1500 | **0.8165** | 3.89 | 0.0337 |
| 80 | 0.3215 | 0.2625 | **0.8165** | 6.15 | 0.1034 |

The conventions differ by a **constant 2/√6 = 0.8165**, not by √(fT). The √(fT) cancels
because source and background are expressed per the *same* Fourier bin. **Corrected
h_c-convention figures are 0.074 (10 nHz) and 0.150 (32 nHz)** — lower than reported, not
2.2–3.9× higher.

## 1. Targets A and B

Neither prompt specifies the source **strain**, **frequency**, or **distance-prior width**,
and the claimed precision depends on all three. So I inverted: what SNR does each claim
require? Using σ_L = dL/(2π ρ_P √N) — Fisher-exact for lnL = ρ²(cos(2π ΔL/dL) − 1), matching
the banked D3 scaling.

| target | f nHz | dL pc | ρ_array needed | **ρ_P per pulsar** | h₀ needed | **P(fringe ID)** |
|---|---|---|---|---|---|---|
| **A — Xiao** (85 psr, 20 yr, 50 ns, N=5) | 5 | 1.94 | 1.28 | 0.138 | 3.3e-17 | **0.00** |
| | 10 | 0.97 | 0.64 | 0.069 | 3.3e-17 | **0.00** |
| | 20 | 0.49 | 0.32 | 0.035 | 3.3e-17 | **0.01** |
| **B — Yu & Pan** (20 psr, 30 yr, 20 ns, N=8) | 5 | 1.94 | 1.22 | 0.273 | 2.1e-17 | **0.00** |
| | 10 | 0.97 | 0.61 | 0.137 | 2.1e-17 | **0.00** |
| | 20 | 0.49 | 0.31 | 0.068 | 2.1e-17 | **0.00** |

**The σ_L claim reproduces trivially** — sub-parsec within-fringe precision needs an array
SNR of order 1. That is the same result as the programme's D2 (116/116 pulsars below their
EM prior) and it is not in dispute.

**What it takes to actually identify the fringe** (1 pc prior, L* = 1 kpc):

| ρ_P per pulsar | 0.2 | 0.5 | 1.0 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|
| A (N=5) | 0.01 | 0.03 | 0.48 | **0.96** | 0.99 | 1.00 |
| B (N=8) | 0.00 | 0.02 | **0.94** | 1.00 | 1.00 | 1.00 |

- **A needs ρ_P ≳ 2**; its claim implies ~0.03–0.14 → **shortfall ~7× in SNR, ~44× in source power**
- **B needs ρ_P ≳ 1**; its claim implies ~0.07–0.27 → **shortfall ~3× in SNR, ~11× in power**

B does better than A because 8 sources beat 5 — consistent with VERNIER-1's K₉₀ = 3 at high
SNR degrading as SNR falls.

### Where I could be wrong about the papers

I want to be careful here, because three unspecified inputs each provide an escape:

1. **Prior width.** If their distance prior is *already* narrower than dL (~1 pc at 10 nHz),
   there is no fringe ambiguity to resolve and their claim stands as written. My test assumed
   a 1 pc prior, which is generous but still holds ~11 fringes at ±3σ (VERNIER-1, Toy 1a).
2. **Their sources may be louder than their claim implies.** My inversion pins ρ from the
   *reported* σ_L. If σ_L is limited by something other than the data, ρ is higher and the
   fringe may well be identified.
3. **Yu & Pan fix the chirp mass**, which VERNIER-2 showed is what breaks the Mc–L degeneracy
   in the pulsar-term frequency handle. I checked whether that rescues them: at their SNR,
   σ_L from the spectral line position is ~5 kpc. It does not.

**So the honest verdict is: the toy validates their precision arithmetic and flags an
unresolved question about fringe identification that their configurations, as specified to
me, do not appear to answer.** I have not read the papers; this is a comparison against the
numbers in your brief.

## 2. Field units — λ_fore

**Both targets set no GWB, so they run at λ_fore = 1.000 exactly** — the most favourable
value possible. That is the single largest difference between their configuration and a real
sky, and it is worth stating plainly before any comparison.

My realistic census (γ = 3):

| f nHz | 2 | 5 | 10 | 20 | 32 | 80 | 150 | 208 |
|---|---|---|---|---|---|---|---|---|
| N_c/bin | 2.5e7 | 8.7e5 | 6.9e4 | 5.4e3 | 965 | 33.5 | 3.34 | 1.01 |
| **λ_fore** | 0.0011 | 0.0035 | **0.0083** | 0.0199 | **0.0368** | 0.138 | 0.474 | 1.000 |
| vs 0.5 / 0.7 | fail | fail | fail | fail | fail | fail | fail | PASS |

(Exact Σ k^(−2/γ), not the 3N^(1/3) asymptote.) **λ_fore fails the 0.5 threshold at every
frequency in or near the band**, and only passes at 208 nHz — where N_c = 1.01, i.e. the
"population" is a single source and λ_fore = 1 trivially. That is not a knee; it is the
bin running out of sources. Between 32 and 150 nHz λ_fore climbs 0.037 → 0.474, but the
f^(−11/3) law is pure extrapolation there.

## 3. The knee — and the correction to VERNIER-4

**You were right that I answered the wrong question.** 208 nHz is the one-source-per-bin
frequency for the *full* population, and it is correct for that question — but the knee that
matters is where the *massive tail* truncates.

**The knee is a statement about the bright-end counts slope, not about N_c.** For
N(>h) ∝ h^−γ we get h_k² ∝ k^(−2/γ), so:

- **γ > 2**: Σ diverges, λ_fore → 0 as N grows. **No knee at any frequency.**
- **γ < 2**: Σ → ζ(2/γ), λ_fore → 1/ζ(2/γ), a *constant* independent of N_c.

| γ | 3.0 | 2.0 | 1.5 | 1.2 | 1.0 | **0.87** | 0.80 |
|---|---|---|---|---|---|---|---|
| λ_fore asymptote | 0 | 0 | 0.278 | 0.471 | 0.608 | **0.700** | 0.745 |

> **λ_fore ≥ 0.7 requires γ ≤ 0.87; λ_fore ≥ 0.5 requires γ ≤ 1.16. Euclidean is γ = 3.**

So the literature's 0.7 threshold **cannot be reached by a better array or a longer
baseline** — it is a property of the population's bright end and requires a hard
mass-function truncation. My census has none, at any frequency.

**What that does to VERNIER-4.** At λ_fore = 0.7 the loudest source is 0.84 of its bin's
amplitude — it *is* the bin. Detection then costs nothing extra, and VERNIER-4's "never
closes under correlated confusion" verdict is **void**, because there is no sea left to be
correlated with. VERNIER-4's "operates only above 13 nHz, closes only near 113 nHz" is an
artefact of γ = 3, not a result about the band.

| λ_fore | h₁/h_sea | implied γ |
|---|---|---|
| 0.034 (my census, 32 nHz) | 0.184 | 1.93 |
| 0.5 | 0.707 | 1.16 |
| 0.7 | 0.837 | 0.87 |

---

## What this settles, and what it does not

**Settles:** the convention question (no √(fT) error); that within-fringe σ_L is cheap and
reproduces both papers; that fringe identification at their implied SNR does not; and that
λ_fore ≥ 0.7 is a population property requiring γ ≤ 0.87.

**Does not settle:**

1. **The bright-end slope γ is now the single load-bearing number in this whole programme.**
   It decides λ_fore, which decides whether the wheel operates at all. I have not measured
   it — I assumed γ = 3 in VERNIER-4 and treated it as free here. **Measuring γ from
   GLACIER/CHORUS supersedes every other follow-up I have proposed.**
2. **I have not read either paper.** All three escape routes above (prior width, actual
   source strain, their inference method) could resolve the discrepancy in their favour.
   This should be checked against the actual texts before it is quoted.
3. **The ρ_P ≳ 1–2 fringe-ID threshold is my toy's**, from VERNIER-1, and carries all its
   caveats: Asimov, no GWB, no timing-model marginalisation, trials-factor threshold rather
   than a calibrated false-alarm rate.
4. **The four-dial walk was not run.** Stopped at the discrepancy as instructed.

## One-line verdict

The toy reproduces both published configurations' **precision** and neither one's **fringe
identification** at the SNR those precisions imply — but the more important finding is that
**λ_fore ≥ 0.7 demands a bright-end counts slope γ ≤ 0.87 against Euclidean 3**, which makes
γ, not the array, the quantity that decides whether the wheel can ever turn, and voids
VERNIER-4's pessimistic verdict.
