# VERNIER-20 — run them together, honestly

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~10 min.
Clean room: imports only `vernier_toy` / `vernier16` / `vernier18` / `vernier19`.

Code `hpc_harbor/vernier/vernier20.py` · data `reports/vernier20_combined.npz`.
Staged, not committed.

---

## Lead — the mechanisms interfere. The 0.552 dex stack is refuted.

Simulated together in one self-consistent noise budget, 5.41–6.27 nHz (λ = 0.811, N = 5),
target = source 2:

| configuration | SNR | f_cert | **gain vs A** | rank |
|---|---|---|---|---|
| **A** baseline: uniform, Earth only, no sub | 2.794 | 0.000 | 0.0000 | 4 |
| **B** + heterogeneity | 2.104 | 0.000 | **−0.1231** | 4 |
| **C** + pulsar-term coherence | 2.122 | 0.000 | −0.1195 | 4 |
| **D** + one subtraction (all three) | 4.246 | 0.199 | **+0.1817** | **9** |
| (D but per-pulsar, not joint) | 3.849 | 0.000 | +0.1390 | 8 |
| **subtraction only, uniform, no pterm** | **4.858** | 0.000 | **+0.2402** | **11** |

> **VERNIER-19's stacked estimate: 0.552 dex → rank 28.**
> **Measured together: 0.1817 dex → rank 9.**
> **And the best configuration is subtraction ALONE on a uniform array: 0.240 dex → rank 11.**

**Adding all three mechanisms is worse than using one of them.** That was worth simulating.

### Why heterogeneity turns negative

**VERNIER-13 measured heterogeneity's benefit against WHITE noise. In a confusion-dominated
bin it evaporates and goes slightly negative** (−0.123 dex).

The reason is structural: confusion scales with the *same* per-pulsar weights as the signal.
A high-weight pulsar sees more signal *and* proportionally more confusion, so its
signal-to-confusion ratio is unchanged — it saturates at R_k²/conf, independent of weight.
Only the white-noise term fails to scale. So heterogeneity concentrates the array's
sensitivity into fewer effective measurements without improving any of them, and the effective
N falls. **VERNIER-13's f_cert = 0.80 and rank 11 stand for a white-noise-limited array; they
do not transfer to a confusion-limited bin.**

### One positive interaction, and it is the interesting one

**f_cert = 0.000 in configurations B and C, and 0.199 only after subtraction.** No pulsar
certifies while the confusion is present, because VERNIER-17's ceiling
ρ_P ≤ √(λ/(1−λ)) = 2.07 sits essentially on the bar of 2. **Removing the top source is what
lets certification happen at all.**

So the naive sum of parts (−0.123 + 0.004 + 0.240 = 0.121 dex) *understates* the measured
0.182, because subtraction enables coherence. The mechanisms are coupled in both directions —
one destructively, one constructively — and neither is captured by adding measured gains.

**Joint vs per-pulsar inference is worth +0.043 dex here** (0.182 vs 0.139), and the whole of
f_cert: per-pulsar treatment gives f_cert = 0.000 even after subtraction.

## Test 2 — joint position fitting does NOT sustain turn 2

| turn | n | greedy P(genuine) | **joint position P(genuine)** |
|---|---|---|---|
| 1 | 60 | 0.633 | **0.717** |
| 2 | 59 | 0.288 | **0.339** |

**Gate is 0.5. Turn 2 reaches 0.339 and fails.**

Your naive scaling predicted ~0.33 and that is what it does — the improvement is roughly
proportional (+13% at turn 1, +18% at turn 2) rather than transformative. **The loop still
manages one turn.** Joint position fitting is a real improvement and it is not enough to buy
a second turn.

## Test 3 — a bad subtraction costs almost nothing

| next-turn P(genuine) | |
|---|---|
| after a **clean** subtraction | 0.375 |
| after a **bad** subtraction (fake source injected at a wrong position) | 0.350 |
| **degradation** | **0.025** |

**This is the reassuring result of the instalment, and it settles the question the brief
posed.** A failed removal injects a fake negative source, but the next search is barely
affected — 2.5 percentage points. The reason is that a blended fit has a *small* amplitude
compared with the real sources (it is fitting a partial coincidence, not a source), so the
fake negative it leaves is weak relative to what remains.

> **So the loop STOPS; it does not actively damage the data.** That distinction matters: the
> mechanism is not ruled out, it simply does not compound. A pipeline that subtracts one source
> and stops is safe even when that subtraction is wrong 25–30% of the time.

---

## What this does not settle

1. **One bin, one target source (rank 2 of 5), one λ.** The interference result should be
   checked at other λ before it is generalised — though the mechanism (confusion scaling with
   the weights) is not λ-specific.
2. **The noise budget treats confusion as adding in quadrature per pulsar**, i.e. as
   incoherent between sources. Real interlopers are coherent signals with their own sky
   positions, which is what makes them subtractable in the first place — so this is the same
   Gaussian approximation VERNIER-17 flagged, and it is the weakest part of the model.
3. **Certification is a hard step at ρ_P ≥ 2.** f_cert jumping 0.000 → 0.199 across the
   subtraction is a step-function artefact of that bar; a graded criterion (VERNIER-6b's soft
   weights) would smooth it and probably raise the B/C configurations somewhat.
4. **Test 3 used a randomly-placed bad subtraction.** A real blend sits *between* true sources
   where the fit has more amplitude, so the true degradation is likely worse than 0.025 —
   this is a lower bound on the damage.
5. **40 pulsars, Earth terms only for the search**, white noise, no GWB.

## One-line verdict

**Simulated together the three mechanisms give 0.182 dex and rank 9, not the 0.552 dex and
rank 28 that adding them suggested — because heterogeneity, whose benefit VERNIER-13 measured
against white noise, goes NEGATIVE (−0.123 dex) in a confusion-dominated bin where the
confusion scales with the same weights as the signal, so that subtraction alone on a uniform
array (0.240 dex, rank 11) beats all three combined; the one constructive coupling is that
subtraction is what lets certification happen at all (f_cert 0.000 → 0.199); joint position
fitting lifts turn 2 from 0.288 to 0.339 and still fails the 0.5 gate, so the loop keeps its
single turn; and a bad subtraction costs only 0.025 in next-turn recovery, so the loop stops
rather than damaging the data — which leaves the mechanism unattractive but not disqualified.**
