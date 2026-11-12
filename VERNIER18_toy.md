# VERNIER-18 — deconfusion on a realistic ladder

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~8 min.
Clean room: imports only `vernier_toy` / `vernier16`.

Code `hpc_harbor/vernier/vernier18.py` · data `reports/vernier18.npz`,
`reports/vernier18_top1.npz`. Staged, not committed.

**Bookkeeping: VERNIER-13 DID run** — report `VERNIER13_toy.md`, staged since that session.
Its findings are appended at the end of this note since you did not receive them.

---

## Lead — the ladder helps, splitting never happens, and clean recovery reaches depth 1

**You were right that equal-SNR injection was the worst case.** Against VERNIER-16's 70–80%
spurious at equal strength, the measured ladders give **32–61%** — a real improvement. But the
usable result is not the spurious rate, it is the depth:

| bin | N | λ_fore | genuine | **split** | spurious | spur rate | **clean depth** |
|---|---|---|---|---|---|---|---|
| 2.00–2.85 | 156 | 0.263 | 5.16 | **0.04** | 2.48 | 0.323 | **0** |
| 2.85–3.71 | 46 | 0.579 | 1.32 | **0.00** | 2.00 | 0.602 | **0** (joint: 1) |
| 3.71–4.56 | 22 | 0.722 | 1.60 | **0.00** | 1.04 | 0.394 | **1** |
| 4.56–5.41 | 11 | 0.636 | 1.32 | **0.00** | 1.64 | 0.554 | **1** |
| 5.41–6.27 | 5 | 0.811 | 1.12 | **0.00** | 0.96 | 0.462 | **1** |

*(M = 1 control: genuine 1.00, spurious 0.00. Pure-noise null: 0.000 per realisation. Search
validated.)*

**Splitting is 0.00–0.04 across every bin — it does not happen.** The failure is blending,
exactly as VERNIER-16 found, and the realistic ladder reduces it without removing it.

**Clean recovery reaches depth 1**: the brightest source, and only it, in every bin with
λ ≥ 0.58. In the crowded 2.00–2.85 nHz bin (156 sources, λ = 0.263) even the top source is
not reliably recovered.

## And depth 1 is enough to matter

**Subtraction only needs the removed source to be real.** Testing the subtract-top-1 strategy
directly — is the single highest-F recovery trustworthy?

| bin | λ | N | **P(top-F is genuine)** | P(it is the brightest) | **depth gain if removed** |
|---|---|---|---|---|---|
| 2.00–2.85 | 0.263 | 156 | 0.683 | 0.133 | 0.066 dex |
| 2.85–3.71 | 0.579 | 46 | 0.717 | 0.533 | 0.188 dex |
| 3.71–4.56 | 0.722 | 22 | 0.767 | 0.717 | **0.278 dex** |
| 4.56–5.41 | 0.636 | 11 | 0.733 | 0.700 | 0.219 dex |
| **5.41–6.27** | **0.811** | **5** | **0.933** | **0.933** | **0.362 dex** |

Removing the top source removes a fraction λ of the bin's power, so everything else gains
0.5·log₁₀(1/(1−λ)).

> **In the 5.41–6.27 nHz bin the top recovery is genuine 93% of the time and removing it is
> worth 0.362 dex — larger than the entire pulsar-term coherence ceiling of 0.193–0.278 dex
> that has bounded every sensitivity route in this series.**

**That is the first mechanism found here that beats the ceiling, and it does so because it is
not a coherence argument at all** — no threshold or trials penalty applies, it is simply
removing power from the denominator.

**But the 7–23% error rate on that single subtraction is the whole risk.** Subtracting a
blended, fictitious source injects a fake negative signal. At λ = 0.811 the 6.7% failure rate
may be tolerable with verification; at λ = 0.263 the 32% failure rate is not, and the gain
there is only 0.066 dex anyway. **The usable window is exactly the high-λ bins — which is the
same window VERNIER-17 identified for certification (λ ≥ 0.70–0.80).** Both mechanisms live
or die in the same few bins.

## Simultaneous fitting: helps, does not rescue

Re-fitting all found amplitudes jointly each cycle (a CLEAN major cycle):

| bin | greedy spur rate | **joint spur rate** | greedy depth | joint depth |
|---|---|---|---|---|
| 2.00–2.85 | 0.323 | 0.395 | 0 | 0 |
| 2.85–3.71 | 0.602 | **0.400** | 0 | **1** |
| 3.71–4.56 | 0.394 | 0.379 | 1 | 1 |
| 4.56–5.41 | 0.554 | 0.606 | 1 | 1 |
| 5.41–6.27 | 0.462 | **0.383** | 1 | 1 |

**Joint amplitude refitting improves two bins and worsens two — it is not the fix.** It
recovers depth 1 in the 2.85–3.71 bin where greedy failed, and cuts the spurious rate there
from 0.60 to 0.40, but the pattern is not systematic at this sample size (25 realisations per
cell).

**This does not settle the algorithm question.** I refit *amplitudes* jointly at fixed
positions — the cheap linear part. A genuinely simultaneous fit would optimise the
*positions* jointly too, which is nonlinear and is where the blend bias actually lives. **So
"simultaneous fitting does not collapse the blend rate" is not established; what is
established is that the cheap half of it does not.**

---

## What this does not settle

1. **Top source at SNR 15** throughout. The reliability of the top-1 recovery scales with
   that; at lower SNR it degrades.
2. **The ladder is h_k² ∝ k^−α** fitted to reproduce λ_fore at the measured occupancy. That
   is one family of ladders consistent with the two numbers; the real census may not be a
   power law.
3. **40 pulsars, Earth terms only, one bin, white noise.** The pulsar-term arm remains
   intractable (VERNIER-16: 5.3e10 sky grid points).
4. **The iterated loop was not run.** Whether removing the top source lets the second be
   recovered cleanly — and so on — is VERNIER-16's Part 3 and needs the verification step
   built first, because a 7–23% subtraction error compounds.
5. **25 realisations per cell** for the greedy/joint comparison; differences below ~0.15 in
   rate are not resolved.

---

# Appendix — VERNIER-13, as reported at the time

**It ran.** Full report in `VERNIER13_toy.md`.

**Heterogeneity is what makes certification possible.** At matched total array sensitivity,
with a log-normal weight spread σ_log(w) = 0.40 dex (calibrated to give top-10 = 64% of
sensitivity, N_eff = 15 of 116 — consistent with "~10 pulsars dominate" and 23-of-68):

| array | ρ_array | **f_cert** | pulsars certifying | top-10 power share |
|---|---|---|---|---|
| **uniform** | 8 | **0.0000** | **0.0** | 8.6% |
| **heterogeneous** | 8 | **0.4001** | **3.0** | 64.2% |
| **uniform** | 20 | **0.0000** | 0.0 | 8.6% |
| **heterogeneous** | 20 | **0.7981** | 20.0 | 64.3% |

**Three pulsars certify and carry 40% of the array's pulsar-term power.** The scissors is an
artefact of assuming identical pulsars: uniform, every pulsar sits at 8/√116 = 0.743 against a
bar of 2 and *all* miss by the same amount; heterogeneous, the best few clear it comfortably.
**f_cert goes 0.000 → 0.400 with no change in detection capability at all.**

Fed back into the wheel: **rank 4 → 7 (ρ_array = 8), rank 4 → 8–9 (ρ_array = 20).** Revising
VERNIER-10's configuration D: f_cert = 1 assumed → rank 12; uniform (f_cert = 0) → **rank 6**;
heterogeneous (f_cert = 0.80) → **rank 11**. **VERNIER-10's rank 13 was reachable only on a
heterogeneous array; on a uniform one it is impossible, not merely optimistic.**

**Design question — a clean answer.** At fixed total observing effort split N ways,
ρ_a ∝ N^(−1/2) while ρ_array is **invariant**. Detection does not care how effort is split;
certification cares enormously:

| N_psr | 10 | **24** | 50 | 116 | 400 |
|---|---|---|---|---|---|
| f_cert (het) | 0.988 | 0.960 | 0.903 | 0.800 | 0.559 |
| rank | 9 | **9** | 9 | 8 | 7 |

**A small array of excellent pulsars beats a large array of mediocre ones at every N, with
the optimum at N ≈ 24** — the floor set by VERNIER-12's joint source-parameter fit
(m ≳ σ_φ² ≈ 24). Uniform arrays fail by a cliff (all certify below N = 100, none above);
heterogeneous ones degrade gracefully (f_cert = 0.56 still at N = 400). **This inverts the
usual PTA instinct: more pulsars is right for detection and the GWB, wrong for distance
certification at fixed effort.**

VERNIER-13's most leveraged caveat: f_cert is a step function of the ρ_P ≥ 2 bar, and
σ_log = 0.40 dex is calibrated rather than reconstructed from actual IPTA values (f_cert spans
0.12–0.55 across σ_log = 0.25–0.50).

## One-line verdict

**On the measured ladders the greedy search never splits a source (0.00–0.04) and blends far
less than at equal SNR (32–61% vs 70–80%), and clean recovery reaches depth 1 — the brightest
source only, in bins with λ ≥ 0.58 — which is enough, because removing just that source is
worth 0.5·log₁₀(1/(1−λ)) = 0.19–0.36 dex, exceeding the pulsar-term coherence ceiling that has
bounded every other route; the catch is that the single top recovery is genuine only 68–93% of
the time, so the usable window is the high-λ bins, which is the same window VERNIER-17 found
for certification — and joint amplitude refitting is not the fix, though joint *position*
fitting remains untested.**
