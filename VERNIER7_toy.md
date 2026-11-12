# VERNIER-7 — seed deep, on the measured ladder

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~4 min.
Clean room: imports only `vernier_toy` / `vernier6` / `vernier6b` (my own anchored code).

Code `hpc_harbor/vernier/vernier7.py` · data `reports/vernier7.npz` ·
plot `reports/VERNIER7_ignition.png`. Staged, not committed.

---

## Lead — ignition exists, and it is unreachable

**Your hypothesis is confirmed on its own terms and defeated by the other constraint.**

- **The sustain condition turns on at K_eff ≈ 100.** Gain/gap climbs through 1.0 there and
  reaches **2.1×** by K_eff = 60–100. Past that point the wheel *would* keep going unaided.
- **The wheel cannot get there.** Its own reach is **rank 5** — one rung past the unaided
  rank 4. Maximum gain measured anywhere is **0.0359 dex**; rank 6 needs **0.0527 dex** just
  to clear SNR 8.
- **And the decisive test: seeding to rank 16 puts K_eff at 96.9 — already below the
  ignition threshold — and the wheel still makes ZERO turns.** The sustain condition is met
  and nothing happens, because advancing one rung was never the binding constraint. Clearing
  SNR 8 is.

**So the wheel can neither start nor sustain, but for two different reasons, and the one
that stops it is not the one this session was designed to test.**

---

## 0. A ladder correction that matters

I built the first ladder as `cumsum(median gap)` and it was wrong. **Gaps are positively
skewed, so median-of-gaps does not sum to median-of-cumulative** — cumsum understates true
depth by **36% by rank 21**. Detectability is set by the cumulative, so the ladder must be
CUM_MEAS and the local gaps must be its differences.

| rank | 2 | 4 | 5 | 9 | 17 | 21 |
|---|---|---|---|---|---|---|
| cumsum(median gap) | 0.141 | 0.271 | 0.313 | 0.427 | 0.558 | 0.606 |
| **measured cumulative** | 0.141 | **0.335** | **0.399** | **0.565** | **0.758** | **0.825** |

Consistent local gaps (diff of cumulative) at ranks 1/4/8/16 are **0.141 / 0.064 / 0.033 /
0.021** — larger than the 0.142/0.043/0.023/0.013 in the brief, which are the median-gap
statistic. **The wheel therefore faces bigger rungs than the brief's figures imply**, and
this alone moved the answer: on the wrong (shallow) ladder the wheel reached rank 10; on the
right one it reaches rank 5.

## 1. Gain per turn vs K_eff — the acceleration is real

Detection threshold disabled so the descent continues past the stall (a diagnostic of the
gain law, not a claim the wheel gets there). 232 turn-samples over 8 realisations.

| K_eff bin | n | gain/turn | R | local gap | **gain/gap** | sustains? |
|---|---|---|---|---|---|---|
| 400–200 | 8 | 0.0046 | 0.003 | 0.0971 | 0.05 | no |
| 200–100 | 118 | 0.0172 | 0.091 | 0.0244 | 0.71 | no |
| **100–60** | 80 | **0.0301** | 0.109 | 0.0143 | **2.10** | **YES** |
| **60–40** | 26 | **0.0359** | 0.154 | 0.0179 | **2.01** | **YES** |

**Not a 4-point artefact.** Gain rises ~8× across the descent (0.0046 → 0.0359) as a smooth
function of K_eff, measured over the whole range rather than four consecutive turns. Both
curves do favour the wheel exactly as you predicted: **gain rises as distances sharpen while
gaps fall with depth, and they cross at K_eff ≈ 100.**

## 2. The two thresholds, and which one binds

| | condition | met? |
|---|---|---|
| **(a) SUSTAIN** | gain > local gap — advancing one rung | **met at K_eff < 100** |
| **(b) ABSOLUTE** | gain > dex needed to lift the next source over SNR 8 | **never met past rank 5** |

| rank | 4 | 5 | 6 | 8 | 16 | 20 |
|---|---|---|---|---|---|---|
| unaided SNR | 9.25 | 7.98 | 7.09 | 5.88 | 3.66 | 3.10 |
| **dex needed to reach 8** | 0 | 0.001 | **0.053** | 0.134 | 0.339 | 0.412 |
| max gain measured | — | — | 0.036 | 0.036 | 0.036 | 0.036 |

**The absolute requirement grows with depth (0.05 → 0.41 dex) while the gain saturates
around 0.036 dex.** The two curves diverge, so there is no depth at which the wheel escapes.

## 3. Seeding under soft weights — VERNIER-3's conclusion survives

This combination had never been run. It changes nothing.

| seed rank | final rank | turns | **own gain** | K_eff start → end | w_true end | null K_eff |
|---|---|---|---|---|---|---|
| 1 | 5 | 4 | **+4** | 202.2 → 140.1 | 0.0268 | 235.1 |
| 4 | 5 | 1 | **+1** | 140.3 → 140.1 | 0.0268 | 235.1 |
| 8 | 8 | 0 | **0** | 106.6 → 106.6 | 0.0286 | 234.8 |
| 16 | 16 | 0 | **0** | 96.9 → 96.9 | 0.0356 | 235.8 |

**Seeding buys nothing under soft weights either.** Seeds 1 and 4 converge on the same final
rank 5; seeds 8 and 16 add zero. VERNIER-3 measured this with the binary gate and concluded
total depth is invariant — **the conclusion is unchanged by the currency change**, which is
worth knowing because the gate was the obvious suspect.

Seed 16 is the sharpest cell in the whole session: **K_eff = 96.9, below the ignition
threshold, and still zero turns.**

## 4. Anti-circularity — passes, with an honest caveat

| seed | confidence × | truth × | null confidence × | verdict |
|---|---|---|---|---|
| 1 | 1.44 | **2.13** | 1.00 | truth outpaces confidence |
| 4 | 1.00 | 0.97 | 1.00 | (no movement) |
| 8 | 1.00 | 1.00 | 1.00 | (no movement) |
| 16 | 1.00 | 1.00 | 1.00 | (no movement) |

**Only seed 1 is a real test** — it is the only one that takes enough turns to move. There,
truth (2.13×) clearly outpaces confidence (1.44×) and the null does not concentrate
(1.00×). At seeds 4/8/16 the check is *trivially* satisfied because nothing happens; I am
not counting those as passes. w_true sits at 0.027–0.036 against a null of 0.0063 — **4–6×
chance at every seed**, so the weights that exist are on the right fringe.

---

## What this does not settle

1. **The gain law was measured with the threshold disabled.** The K_eff < 100 regime is
   reached only by forcing detections the wheel could not make. The gain/gap = 2.1 figure is
   real as a property of the likelihood, but it describes a state the wheel never occupies
   unaided.
2. **ρ_P = 0.9 per pulsar per source is a choice**, as is the 30-pulsar array. The absolute
   floor is set by the array, and VERNIER-4 already showed what array would move it (~10⁴
   pulsars at 10 ns under an optimistic confusion model).
3. **Ranks beyond 21 use an extrapolated gap** (0.30/k, fit to the measured tail). Nothing in
   the conclusion depends on it, since the wheel never gets past rank 5.
4. **Still no GWB, sources pre-resolved, white noise, known sky positions.** Same optimistic
   frame as VERNIER-6/6b; realism was not added.
5. **The code prints "measured reach rank 16"** — that is just seed 16's starting rank, not
   the wheel's own reach. The wheel's reach is **rank 5**. Corrected here; the npz carries
   both.

## One-line verdict

**Ignition is real and sits at K_eff ≈ 100 where gain/gap reaches 2.1 — but the wheel's own
reach is rank 5 against an absolute floor that demands 0.05–0.41 dex while gain saturates at
0.036, and seeding to rank 16 puts K_eff below the ignition threshold and still yields zero
turns, so the sustain condition is met and irrelevant: the binding constraint was never the
gap between rungs, it is whether the next source clears detection at all.**
