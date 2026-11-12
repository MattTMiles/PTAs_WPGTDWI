# VERNIER-16 — does deconfusion feed the wheel?

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~4 min.
Clean room: imports only `vernier_toy`.

Code `hpc_harbor/vernier/vernier16.py` · data `reports/vernier16_part1.npz`.
Staged, not committed.

**Parts 2 and 3 were NOT run. Part 1 failed its gate.** Per the brief's instruction —
"if the spurious rate is high, Part 2 is meaningless and you should say so and stop."

---

## Lead — the search does not split sources, it blends them, and that is worse

| M injected (each at SNR 12) | recovered | genuine | **SPLIT** | **SPURIOUS** | split rate | **spur rate** | frac of truth found |
|---|---|---|---|---|---|---|---|
| **1** | 1.00 | 0.97 | 0.00 | 0.03 | 0.000 | **0.025** | **0.97** |
| 2 | 2.58 | 1.18 | 0.07 | 1.32 | 0.029 | **0.515** | 0.59 |
| 4 | 5.50 | 1.55 | 0.07 | 3.88 | 0.014 | **0.705** | 0.39 |
| 8 | 9.28 | 1.82 | 0.05 | 7.40 | 0.005 | **0.798** | 0.23 |
| 14 | 11.80 | 2.27 | 0.12 | 9.40 | 0.011 | **0.797** | 0.16 |
| 25 | 16.75 | 3.77 | 0.30 | 12.68 | 0.018 | **0.757** | 0.15 |

**Scrambled null (pure noise): 0.010 spurious per realisation.**

**The failure mode you were worried about does not happen. Splitting is ≤3% at every M.**
A different failure happens instead, and it is worse: **70–80% of recoveries sit at sky
positions that match no true source at all.** The greedy search fits one template to the
superposition of several sources and lands somewhere that belongs to none of them. A split at
least sits *near* a real source; a blend sits nowhere near one, so subtracting it injects a
fake signal at a wrong position — precisely the corruption the brief anticipated, arriving by
a different route.

**It is confusion, not noise.** The pure-noise null yields 0.010 spurious recoveries per
realisation; the M ≥ 4 cases yield 4–13. The spurious sources are made of real signal placed
wrongly.

**It sets in immediately.** At M = 2 — two well-detected sources — the spurious rate is
already 52%. There is no comfortable regime above one source.

**Therefore Part 2 and Part 3 were not run.** Measuring "residual confusion power after
subtracting the recovered sources" when 3 of every 4 recoveries are at fictitious positions
would produce a number with no meaning, and iterating that loop is, in the brief's own words,
a machine for manufacturing confident nonsense.

## Two design errors in my own experiment, both caught by controls

I report these because the first version of Part 1 gave a plausible-looking answer that was
an artefact:

1. **I injected sources below the detection threshold.** Drawing the two complex F-statistic
   coefficients from N(0,1) gave a median injected |signal| of 4.89 against a threshold
   corresponding to 3.39 — so a large fraction of "sources" were never detectable, and the
   search was scored on failing to find them. The **M = 1 control** exposed it: with a single
   source and therefore *zero possible confusion*, the spurious rate was 33%, which is
   impossible if the search works. Fixing the injection to a stated SNR gives M = 1 →
   genuine 0.97, spurious 0.025. **The search is sound; the experiment was not.**
2. **My matching tolerance was too loose.** At 0.35 rad a random sky position falls within
   tolerance of *some* true source 54% of the time at M = 25 — so "genuine" counts were
   inflated and "spurious" understated. Tightened to 0.15 rad, where the chance floor is
   4–13% and the measured spurious rate of 70–80% clears it by a wide margin.

Both corrections made the result *worse*, not better.

## And a feasibility wall found before any of this

Before building the search I checked whether the pulsar-term arm was even runnable. Using
VERNIER-10's decorrelation scale δθ = c/(2πf_P L_max):

| | grid points | compute for this sweep |
|---|---|---|
| Earth-only | 3.0e3 | 0.4 s |
| **pulsar-term** | **5.3e10** | **~0.2 years** |

**A brute-force pulsar-term sky search is roughly seven orders of magnitude beyond
feasible** — the same wall Xiao et al. describe as making the full joint analysis
"computationally intractable". So the VERNIER-15 result (8 → 20 separable sources) is an
*identifiability bound that no tractable search in this instalment can attain*, and the
comparison in the table above is Earth-term only.

---

## What this does and does not establish

**Establishes:** a greedy CLEAN-style search on Earth terms cannot deliver the deconfusion
chain. It recovers 15–23% of the true sources at M ≥ 6 and produces 3 spurious detections for
every genuine one. The `resolve → subtract → lower the floor → go deeper` chain breaks at the
first link.

**Does not establish that the chain is impossible.** This is one specific and rather weak
search:

1. **Greedy sequential subtraction is the worst reasonable algorithm for confusion.** Fitting
   sources *simultaneously* — or a trans-dimensional sampler that can add, remove and move
   components — should do substantially better, and VERNIER-15's rank result says the
   information is present. **The honest statement is "greedy CLEAN cannot do it", not
   "it cannot be done".**
2. **All sources were injected at equal SNR 12**, the hardest case for a greedy method
   (nothing dominates, so the first fit is always a blend). A realistic brightness ladder
   would let the loudest be peeled off cleanly first.
3. **Earth terms only**, for the feasibility reason above — so this says nothing about
   whether pulsar terms would improve *recovery*, only that a brute-force search for them is
   intractable.
4. **40 pulsars, one frequency bin, no red noise or GWB.**

## One-line verdict

**Part 1 fails its gate and Parts 2 and 3 were not run: the search does not split sources
(≤3% at every M) but it blends them — 70–80% of recoveries land at sky positions matching no
true source, rising from 52% at just two sources, against a pure-noise null of 0.01 per
realisation — so subtracting the recoveries would inject fake signal at fictitious positions
and any measured depth gain would be meaningless; two errors in my own setup (sub-threshold
injections, and a matching tolerance loose enough to call 54% of random positions genuine)
were caught by an M = 1 control and both made the answer worse; and the conclusion is
specifically that greedy CLEAN cannot deliver the deconfusion chain, not that a simultaneous
or trans-dimensional search could not.**
