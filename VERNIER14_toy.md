# VERNIER-14 — reproduce them on their own terms

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~6 min.
Clean room: imports only `vernier_toy` / `vernier6`.

Code `hpc_harbor/vernier/vernier14.py`. Staged, not committed.

---

## Lead

**Yes, we reproduce them — comfortably — and their selection procedure holds its coverage.
Two of our three criticisms do not survive contact with their actual configuration.**

1. **VERNIER-5's failure to reproduce was our error, not theirs.** We inverted the SNR from
   their *claimed* precision. Their stated source gives ρ_P = 7–13 per pulsar, not the
   0.03–0.27 we inferred — a factor 40–100.
2. **The pairwise collapse does not survive their setup.** VERNIER-12's "pairs give zero
   contrast" is **configuration-specific** and must be restated that way.
3. **The min-over-84 selection does not overstate precision.** Coverage is **0.95 against a
   nominal 0.68** — it over-covers. **We drop that concern.**

---

## Test 1 — do we reproduce them? Yes, with room to spare

Computing the source from their stated parameters (Mc = 5×10⁹ M☉, d_L = 1 Gpc) rather than
inverting from their claim:

| f_GW | h₀ | residual | ρ_array (50 ns) | **ρ_P per pulsar** | σ_L (SNR-limited) |
|---|---|---|---|---|---|
| 3 nHz | 1.81e-15 | 95.9 ns | 116.6 | **12.65** | 0.018 pc |
| 10 nHz | 4.03e-15 | 64.2 ns | 78.0 | **8.47** | 0.008 pc |
| 15 nHz | 5.29e-15 | 56.1 ns | 68.2 | **7.40** | 0.006 pc |

**Sub-parsec is reached with ~100× margin.** Their claim of sub-parsec for 55–87% of pulsars
is *conservative* relative to what their configuration supports — the SNR-limited floor is
0.006–0.018 pc.

**Where VERNIER-5 went wrong, precisely.** It assumed σ_L was SNR-limited and inverted
ρ_P = dL/(2π σ_L √N) from the quoted 0.4–1 pc. That inversion is only valid if the quoted
precision *is* the SNR floor. It is not — theirs is limited by marginalisation over Mc, f and
d_L and by the pair procedure, which sit ~100× above the SNR floor. **The single assumption
that carried the whole VERNIER-5 discrepancy was treating their quoted number as
SNR-limited.** Their distance prior (10 pc, SKA-era timing parallax) is also *wider* than the
1 pc we generously assumed, so that was not the issue either.

Yu & Pan reproduce even more easily: ρ_P = 23–39 per pulsar, so their 0.4 pc at 1 kpc is
likewise well inside the configuration's reach.

## Test 2 — the pairwise collapse is configuration-specific

Their sky positions are fixed, so only frequency (and Mc) smear. With
σ_f = 1/(2πTρ_array), effective per-pulsar ρ after decoherence:

**Xiao (85 psr, 20 yr, 50 ns):**

| f | σ_f/f | per-pulsar | **PAIR (m=2)** | joint (m=85) |
|---|---|---|---|---|
| 3 nHz | 7.2e-4 | 4.75 | **7.75** | 12.50 |
| 10 nHz | 3.2e-4 | **0.95** | **2.84** | 8.25 |
| 15 nHz | 2.5e-4 | **0.42** | 1.77 | 7.15 |

**Yu & Pan (20 psr, 30 yr, 20 ns):** per-pulsar 12.7–31.8, pairs 17.0–35.1 — everything
works, the smear barely matters.

**Pairs clear the certification bar (~2) in Xiao's setup at 3 and 10 nHz, where per-pulsar
fails at 10 and 15 nHz.** So pairing is doing real work for them and it is sufficient.

**Why VERNIER-12 concluded otherwise:** it had ρ_P = 5 with σ_φ = 4.85 rad, giving pairs
5·exp(−3.43²/2) = 0.014. Xiao has ρ_P = 8.5 with σ_φ = 2.09 rad — **half the smear and 1.7×
the SNR** — which is the difference between exp(−5.9) and exp(−1.1).

> **Restating VERNIER-12 correctly: "pairs capture 0% of the gain" holds when σ_φ ≳ 4 rad. It
> does not hold when σ_φ ≲ 2 rad, which is the regime a loud source with fixed sky position
> produces. The m ≳ σ_φ² ≈ 24 threshold is right, but σ_φ is configuration-dependent and
> Xiao's is small enough that m = 2 suffices.** Their "full joint is computationally
> intractable" is therefore not a compromise that costs them the result.

## Test 3 — the selection effect: coverage holds, and we drop the concern

Simulating their procedure honestly — 84 pair posteriors per pulsar, adopt the narrowest:

| | quoted half-width | **actual error** | **coverage** (nominal 0.68) |
|---|---|---|---|
| **SELECTED (narrowest)** | 0.0189 pc | 0.0058 pc | **0.950** |
| typical (unselected) | 1.513 pc | 0.0075 pc | 0.983 |

**The selection does not produce undercoverage. It over-covers** — the quoted interval is
3.26× larger than the actual error.

**Why our criticism was wrong here.** Min-of-N produces selection bias when the *width
estimate* fluctuates independently of the *accuracy*. In this problem the posterior width is
dominated by whether the fringe degeneracy collapsed — a genuine, data-driven feature — and
when it collapses the estimate really is good. Selecting the narrowest therefore selects
genuinely-resolved cases, not lucky noise. **The banked convention "an order statistic is not
a threshold" is about detection thresholds; it does not automatically transfer to interval
coverage, and here it does not.**

**But the null shows a real limit:**

| scrambled null | quoted half-width | actual error | coverage |
|---|---|---|---|
| selected | **0.0190 pc** | **13.35 pc** | **0.000** |

On data with no coherent signal the procedure returns a 0.019 pc interval around an answer
that is 13 pc wrong. **The concern therefore shifts rather than vanishes: the precision is
not overstated when a source is genuinely present, but the procedure carries no internal
protection against being run on a non-detection.** It must be gated on a detection statistic
computed beforehand — which is a requirement on the pipeline, not a defect in the precision
they quote.

---

## What this does not settle

1. **I have still not read the papers** — this is their configuration as stated in your
   brief, and my reading of "adopt the smallest half-width" as min-over-partners.
2. **The pair model is a caricature**: partner quality enters as a single scalar drawn from
   U(0.3, 1), and the shared-parameter constraint scales as 1/√m. A real pair analysis
   marginalises the source parameters explicitly.
3. **Only frequency smear is modelled.** Chirp mass and luminosity distance are also
   marginalised in their analysis; d_L is degenerate with amplitude and should be benign, but
   Mc enters t_c and would add smear I have not included.
4. **Coverage was measured at one frequency (10 nHz) and one true distance (1 kpc)**, with
   120 realisations. The over-coverage is large enough (0.95 vs 0.68) that it is unlikely to
   be a sampling artefact, but it is not a stress test.
5. **Their 55–87% figure is a fraction across a real pulsar population** with heterogeneous
   distances and geometries; I checked the SNR floor at median geometry rather than
   reproducing the distribution.

## One-line verdict

**We reproduce both papers with ~100× margin once their actual source is used instead of an
SNR inverted from their claim — the VERNIER-5 discrepancy was entirely our assumption that
their quoted precision was SNR-limited; their pair-by-pair treatment is sufficient in their
configuration because a fixed sky and a loud source give σ_φ ≈ 2 rad rather than the ≈ 5 rad
that made VERNIER-12's pairs useless, so that conclusion is configuration-specific; and their
min-over-84 selection over-covers (0.95 against nominal 0.68) rather than overstating
precision, so we drop that criticism — with the single caveat that on a scrambled null the
same procedure returns a 0.019 pc interval around a 13 pc error, so it must be gated on a
prior detection statistic.**
