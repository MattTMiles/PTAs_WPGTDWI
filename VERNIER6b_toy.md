# VERNIER-6b — soft fringe weights, carried across the wheel

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~2 min.
Clean room: imports only `vernier_toy` / `vernier6` (my own anchored code).

Code `hpc_harbor/vernier/vernier6b.py` · data `reports/vernier6b.npz` ·
plot `reports/VERNIER6b_softweights.png`. Staged, not committed.

---

## Lead

**Soft weights do exactly what you predicted: every turn now gains something, and the loop
that VERNIER-3 reported as a flat zero is visibly moving. It is not enough to make the wheel
self-sustaining.** The wheel runs 4 turns instead of stalling at 1, K_eff falls by 1.6× while
the null stays flat, and the iteration converges *toward* the joint fit rather than away from
it — but the coherence it recovers is R ≈ 0.05, worth **0.013 dex per turn** against rungs
0.12 dex apart.

Settings — the **stalled** regime on purpose (ρ_P = 0.9 per pulsar per source, 30 pulsars,
30 yr, 20 ns, no GWB, e = 0.7 with 12 harmonics, 30 pc prior, 10 sources on a 0.12 dex
ladder).

| turn | n_det | **K_eff signal** | 16/84 | **K_eff null** | **w_true signal** | **w_true null** | R | gain dex |
|---|---|---|---|---|---|---|---|---|
| 1 | 1 | **196.0** | 159/240 | 224.5 | 0.0133 | 0.0073 | 0.001 | 0.0000 |
| 2 | 2 | 157.0 | 126/189 | 224.6 | 0.0179 | 0.0073 | 0.000 | 0.0007 |
| 3 | 3 | 139.8 | 112/152 | 222.5 | 0.0274 | 0.0073 | 0.018 | 0.0074 |
| 4 | 4 | **122.6** | 108/146 | 222.2 | **0.0256** | 0.0072 | 0.048 | **0.0131** |

Turns before stall: **median 4, range 3–4** — against VERNIER-3's **1 in 30/30 runs**.

---

## The metric works

**K_eff = exp(H(w)) is monotone, continuous, and non-zero exactly where the binary metric
read zero.** It falls 196 → 123 (**1.60×**) while the scrambled null sits flat at 222–225.
The distribution across pulsars is reported, not just the best: the 16/84 band tightens from
159/240 to 108/146, so the *whole array* is moving, not a lucky tail.

Nothing is ever collapsed. The pulsar's state is the full w(L); the gain comes from the
resultant length R = |∫ w e^{iφ_next} dL|, which is 0 for uniform weights and 1 for a delta.

## Anti-circularity — all three checks, and they pass

| check | result |
|---|---|
| **confidence** gain (K_eff ratio) | 1.60× |
| **truth** gain (w_true ratio) | **1.92×** |
| **null** confidence gain | **1.01×** (must be ~1) |

**Truth outpaces confidence (1.92× vs 1.60×), and the null does not concentrate at all.**
That is the signature we needed: the weights are moving toward the right fringe, not merely
becoming more certain. w_true ends at **0.0256 against a null of 0.0072 — 3.56× chance.**

Implementation of the three requirements:
1. **w is re-derived from the accumulated per-source likelihoods every turn.** No decision is
   propagated; only the *set* of detected sources carries forward. The per-(pulsar, source)
   noise is drawn once and cached, so re-deriving cannot re-roll the data in the wheel's
   favour.
2. **The scrambled null runs every turn** through the identical pipeline, scrambled at the
   *source* level (the bug VERNIER-6 found and fixed).
3. **w_true is reported separately from K_eff** in every row above.

## Is the iteration legitimate?

| | K_eff | w_true |
|---|---|---|
| joint fit, all 10 sources at once | **119.8** | 0.0265 |
| iterative wheel at stall (4 sources) | **122.6** | 0.0256 |
| null (joint) | 220.9 | — |

**The wheel reaches 95% of the joint fit's log-K_eff reduction. It converges toward the joint
solution and stops short of it; it does not drift away.** The shortfall is the sources it
never detected, not a bias in the ones it did.

**I should flag that the 95% flatters it.** On this 0.12 dex ladder, source 10 is 1.08 dex
below source 1, so its ρ² is ~0.007× — sources 5–10 carry almost no distance information.
The wheel gets 95% partly because there was little left to get. On a flatter ladder the gap
would be wider.

## The verdict on the currency change

**It converts a dead loop into a live one, and a live loop into a slow one.**

- VERNIER-3: certification gated at q > 0.9, nothing certified, f_cert = 0, gain = 0, stall at
  turn 1 in 30/30. **The loop could not start.**
- VERNIER-6b: every turn recovers a real coherence fraction, K_eff falls monotonically,
  4 turns. **The loop starts and inches.**
- But R climbs only to 0.048, so the gain is **0.0131 dex/turn** against ladder rungs of
  0.12 dex. It needs **~9× more gain** to buy one additional rung.

The reason is visible in the numbers: R ≈ 1/√K_eff to order of magnitude, and K_eff is still
123. **Coherence requires the distance pinned to better than the *next* source's fringe
spacing, and 123 effective fringes is nowhere near that.** Soft weights extract the
information that was being thrown away by the gate — which is real and was worth doing — but
the information was never large.

## What this does not settle

1. **The gain-per-turn is the whole question and it is small.** "A loop that inches can
   compound" is true in principle; here the inch is 0.013 dex and the step is 0.12 dex, so it
   compounds too slowly to reach the next rung before the sources run out.
2. **R = |∫ w e^{iφ}| assumes the pulsar term can be added coherently once the phase is
   known**, with no penalty for the amplitude and sky-position uncertainty that would come
   with it. That is optimistic.
3. **Still no GWB, sources pre-resolved, white noise only, sky positions known.** Same
   optimistic frame as VERNIER-6; realism was not added.
4. **The ladder is mine.** 0.12 dex spacing was chosen to put the wheel in the interesting
   regime. The conclusion "inches too small" is a statement about gain-vs-rung, and the rung
   is a choice — on a census with 0.02 dex rungs near rank 8 (VERNIER-3's measured ladder)
   0.013 dex/turn would be nearly enough, which is worth testing.
5. **K_eff is binned at a fixed reference fringe spacing** (median dL of source 0). It is
   comparable across turns by construction, but its absolute value depends on that choice.

## One-line verdict

**Soft weights are the right currency — they make every turn gain something, K_eff falls
1.6× while the null stays flat, truth outpaces confidence 1.92× vs 1.60×, and the iteration
converges to 95% of the joint fit — but the recovered coherence is R ≈ 0.05, worth 0.013 dex
per turn against 0.12 dex rungs, so the wheel now inches instead of stalling and still does
not sustain.**
