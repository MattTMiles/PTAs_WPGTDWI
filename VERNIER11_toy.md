# VERNIER-11 — population priors

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~4 min.
Clean room: imports only `vernier7` / `vernier8` for the measured census ladder.

Code `hpc_harbor/vernier/vernier11.py`. Staged, not committed.

---

## The deliverable: two numbers, and they say drop it

| | **yield gain** |
|---|---|
| **NUMBER 1 — signal** | **1.60×** |
| **NUMBER 2 — the NULL** | **5.98×** |

**The null gains 3.7× more than the signal.** This is not merely "comparable" — the
population prior is *anti-differential*: at a threshold held fixed, it raises the
false-alarm rate from 1.0% to 6.0% while raising the detection rate only from 0.348 to 0.557.

**And once the threshold is recalibrated on the null — which is mandatory — the gain is
0.88×, i.e. slightly negative.** Ranks recovered: 8 with a flat prior, 7 with the population
prior, at matched FPR = 0.01.

**Plainly: this is not real, and we should drop it.** It is not the cheapest sensitivity gain
in the programme; it is not a gain at all.

## Why — the mechanism is a pure Occam shift

For one source with everything known but amplitude, the matched-filter output is
x ~ N(ρ_true, 1) and the evidence is Z = ∫p(ρ)exp(ρx − ρ²/2)dρ. Everything about the prior
enters through p(ρ), so this is exact rather than a simulation artefact.

A narrower, better-placed prior raises lnB by removing the Occam penalty
ln(prior width / posterior width) ≈ ln(5 dex / 0.3 dex) ≈ 2.8 nats. **But that shift is
almost entirely data-independent — it applies whether or not a source is present.** The null's
lnB rises with it, the null distribution shifts right, and the threshold that preserves a
fixed false-alarm rate rises by the same amount. Median lnB on the null moves 0.856 → the
99th percentile moves 0.868 → 2.150. The gain is exactly cancelled.

The residual differential effect is small and, here, *negative*: the informed prior recovers
7 ranks vs the flat prior's 8.

## The naive LF prior is actively harmful

Using the census luminosity function literally as an amplitude prior — p(ln h) ∝ h^−(γ−1)
with γ = 4 — is **much worse than flat**, not better:

| prior | median lnB on signal | detection rate at fixed threshold |
|---|---|---|
| flat 5 dex | 7.07 | 0.900 |
| **naive LF, dN/dh ∝ h⁻⁴** | **0.001** | **0.187 (0.21×)** |

Because dN/dh ∝ h⁻⁴ is steeply falling, essentially all the prior mass sits at the faint end
— at strains far below detectability. The prior assigns ~zero weight to the amplitudes a
detectable source actually has, and the evidence integral collapses. **"Astrophysically
motivated" is not the same as "helpful for detection": the LF is a prior on a random member
of the population, not on the loudest one you are trying to detect.** The order-statistic
prior (used as "informed" above) is the correct object, and even it does not help.

## Misspecification — nothing to lose, so nothing is lost

| variant | threshold | FPR | ranks > 50% |
|---|---|---|---|
| shift −0.6 dex | 1.135 | 0.010 | 8 |
| shift −0.3 dex | 1.762 | 0.010 | 7 |
| **exact** | 2.150 | 0.010 | **8** |
| shift +0.3 dex | 2.013 | 0.010 | 8 |
| shift +0.6 dex | 1.247 | 0.010 | 7 |
| wrong slope −2 | 1.531 | 0.010 | 7 |
| wrong slope +2 | 1.793 | 0.010 | 7 |

**Robustness here is not reassuring — it is diagnostic.** The yield is 7–8 ranks whether the
prior is exactly right or 0.6 dex wrong with an inverted slope, because at fixed FPR the
prior is not doing any work in the first place. A prior that cannot be broken by
misspecification is a prior that was not contributing.

## Are the extra detections real?

At the *uncalibrated* threshold the informed prior finds sources the flat prior missed in
**20.9% of signal trials**, and a source genuinely is present in every one — so they are true
detections, not fabrications. **But they are bought at a 6× inflated false-alarm rate.** That
is the entire content of the "gain": moving the operating point down the ROC curve, which a
flat prior can do just as well by lowering its own threshold. Nothing about the population
prior improves the ROC curve itself.

## Against the 21× claim

I cannot reproduce anything like 21× as a *sensitivity* gain, and I do not think the two
numbers are measuring the same thing. A 21-of-114 vs 1-of-114 change in **catalogue flagging**
is a statement about how many marginal candidates sit near a significance cut — if the
candidate significances are piled up near threshold, a small shift moves many across it. That
is the source-density-near-threshold effect VERNIER-2 already identified, and it is a property
of the catalogue, not an improvement in the detector. **My result says the underlying
detection capability does not change; whether a published flagging count changes by 21× is a
separate question about that catalogue's significance distribution, which this toy cannot
address.**

---

## What this does not settle

1. **I tested an amplitude prior only.** Priors on chirp mass and frequency enter differently
   — they can genuinely reduce the search's trials factor (VERNIER-10 measured that channel:
   worth 0.079 dex). A population prior that shrinks the *parameter space searched* is not the
   same object as one that reweights amplitude, and I have not tested the former.
2. **Single-source, single-parameter model.** The real analysis is hierarchical across 114
   candidates with shared population hyperparameters; the common-mode cancellation could
   differ when the prior is itself being inferred from the data.
3. **The null here is pure Gaussian noise**, not the source-level scramble used elsewhere in
   the series — appropriate for this evidence calculation, but a different null.
4. **I have not read Goncharov et al.** This is a comparison against the numbers as stated in
   your brief.

## One-line verdict

**Drop it: the signal yield gain is 1.60× and the null's is 5.98×, so the population prior
raises false alarms nearly four times faster than detections, and at a properly recalibrated
false-alarm rate the gain is 0.88× — because the prior's benefit is a common-mode Occam
shift that the null absorbs entirely; the naive luminosity-function prior is worse still
(0.21× detection rate), since the LF describes a random population member rather than the
loudest source one is trying to detect.**
