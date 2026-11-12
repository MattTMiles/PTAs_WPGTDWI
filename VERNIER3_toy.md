# VERNIER-3 — the seeded wheel

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~4 min compute.
Clean room preserved: imports only `vernier_toy` / `vernier2` (my own anchored code).

Code `hpc_harbor/vernier/vernier3.py` · data `reports/vernier3.npz` ·
plot `reports/VERNIER3_seeded.png`. Staged, not committed.

---

## Lead

**The engine is M2 — source subtraction — not M1 certification. At every seed N, and by
1.3–1.6×.** Running both is *superadditive*: 0.284 dex versus 0.090 (M1) + 0.137 (M2) = 0.227
summed, because subtraction lowers the noise that the pulsar terms are measured against, which
makes certification both easier and worth more.

| N seed | gain M1 | gain M2 | gain both | turns none / M1 / M2 / both | engine |
|---|---|---|---|---|---|
| 1 | 0.0904 | **0.1366** | 0.2838 | 4 / 6.5 / 7 / **15** | M2 |
| 2 | 0.1052 | **0.1329** | 0.2891 | 3.5 / 7.5 / 8.5 / **19** | M2 |
| 4 | 0.1114 | **0.1545** | 0.3077 | 0.5 / 4.5 / 7.5 / **15** | M2 |
| 8 | 0.1146 | **0.1715** | 0.3148 | 0 / 0 / 0 / **7** | M2 |
| 16 | 0.1451 | **0.1851** | 0.3437 | 0 / 0 / 0 / **4** | M2 |

**There is no critical N. The bulk-certification hypothesis is refuted.** f_cert does not step —
it ramps smoothly at ~0.2 per doubling of N across the whole range:

| N | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| f_cert | 0.209 | 0.409 | 0.626 | 0.806 | 0.942 |

Certification *is* a bulk event in the sense that one source certifies ~20 % of the array's
pulsar-term power at a stroke — but there is no threshold to cross. The largest single jump is
N = 2 → 4 (0.409 → 0.626), consistent with VERNIER-2's fractional K₉₀ = 2.32 at 1 pc, and it is
a factor 1.5, not 50. The expected "1/116 → 0.6" transition never happens because **N = 1
already delivers f_cert = 0.21, not 1/116** — a single source's comb is seen by the whole array
at once.

**The absolute floor binds, at every N, without exception.** This is the headline correction to
VERNIER-2.

---

## Gap versus gearing: the crossing, and why it does not matter

Measured gaps (20 000 realisations, 256 members, σ_logh = 0.6):

| rank k | gap(k) dex | cum. dex below rank 1 | SNR of rank k+1 | gap < M1 ceiling | gap < M1 realistic |
|---|---|---|---|---|---|
| 1 | **0.142** | 0.142 | 14.4 | YES | no |
| 2 | 0.074 | 0.255 | 11.1 | YES | YES |
| 4 | **0.043** | 0.397 | **8.02** | YES | YES |
| 8 | **0.023** | 0.562 | 5.48 | YES | YES |
| 16 | **0.013** | 0.758 | 3.49 | YES | YES |

Your rough expectation (0.20 / 0.06 / 0.03 / 0.02 at ranks 1 / 4 / 8 / 16) has the right shape
and is **~30 % high**; measured is 0.142 / 0.043 / 0.023 / 0.013.

- **Crossing: gap(k) falls below realistic M1 gearing (0.108 dex at f_cert = 0.45) at k = 2.**
  Below the M1 *ceiling* (0.193 dex) at k = 1.
- **Absolute floor: rank k+1 drops under SNR 8 at k = 4**, i.e. 0.399 dex below rank 1.

**So there is no window in which the relative gap is what binds.** The gap stops mattering at
rank 1–2 and the floor starts biting at rank 4. Past rank 2 the gaps are 0.01–0.05 dex — a
fifth of the gearing — and yet the wheel still stalls, because the sources are simply too faint.

At the stall, every seed classifies the same way:

| N seed | gap to next | gain delivered | SNR next | SNR no-gear | max conceivable reach | binds |
|---|---|---|---|---|---|---|
| 1 | 0.060 | 0.321 | 7.56 | 3.54 | 5.52 | **ABSOLUTE** |
| 2 | 0.094 | 0.328 | 7.69 | 3.70 | 5.76 | **ABSOLUTE** |
| 4 | 0.043 | 0.331 | 7.56 | 3.51 | 5.47 | **ABSOLUTE** |
| 8 | 0.041 | 0.318 | 7.38 | 3.40 | 5.30 | **ABSOLUTE** |
| 16 | 0.038 | 0.352 | 7.36 | 3.15 | 4.90 | **ABSOLUTE** |

The gain delivered (0.32–0.35 dex) exceeds the gap to the next source (0.04–0.09 dex) by 4–8×
in every case, and the wheel stalls anyway. **A small gap to an undetectable source is still a
stall, and that is the regime the census actually lives in.**

### Reconciling with VERNIER-2

VERNIER-2 concluded the brightness gap "binds first, and binds alone". That stands *for the
dial it walked* — an imposed uniform ladder of 0.15–1.2 dex. **But a real census does not have
uniform 0.5 dex rungs.** Natural gaps past rank 1 are 0.01–0.06 dex, an order of magnitude
inside the gearing. VERNIER-2's gap dial was exploring a regime that does not occur beyond
rank 1. **Census-realistic answer: the absolute detection floor binds, not the gap.** The
0.77–1.17 dex "conditioned rungs" figure is an artefact of conditioning on rank 0, not a
property of the luminosity function.

## How many turns each seed buys, and where it saturates

**Total sources known (seed + turns), median over 20 matched realisations:**

| arm | N=1 | N=2 | N=4 | N=8 | N=16 |
|---|---|---|---|---|---|
| none | 5 | 5.5 | 4.5 | 8 | 16 |
| M1 only | 7.5 | 9.5 | 8.5 | 8 | 16 |
| M2 only | 8 | 10.5 | 11.5 | 8 | 16 |
| **both** | **16** | **21** | **19** | **15** | **20** |

**Seeding does not extend reach.** With both mechanisms the wheel converges to a total of
**~15–21 sources (rank ≈ 20) regardless of whether it was handed 2 or 16.** Handing it 16
instead of 2 buys *nothing* — total goes 21 → 20, flat within scatter. The saturation depth is
set by the absolute floor: at rank 20 the no-gearing SNR is ~3.0, and 0.34 dex of gearing lifts
it to ~6.6, still under 8.

The turns-beyond-seed column falls with N purely because seeding consumes the bright end. The
invariant is the **total depth**, and it is fixed.

So the answer to "does the wheel self-sustain if handed N sources" is: **it runs a long way —
about 15–19 further detections with both mechanisms, far more than VERNIER-2's ~3 — but it
self-terminates at a fixed depth that no seed can move.**

---

## What this cannot settle — and the caveat that matters most

1. **M2's size is a direct function of two numbers I chose, and both are optimistic.**
   The analytic form is gain = ½ log₁₀[(1+C)/(1+C(1−f))], with C the confusion-to-white ratio
   and f the band-power fraction removed. My toy's measured 0.185 dex at N=16 matches this
   formula's 0.199 dex at C = 3, f = 0.83 — the mechanism is right, but the *magnitude* is set
   by C and f:

   | | f = 0.35 | f = 0.55 | f = 0.80 |
   |---|---|---|---|
   | C = 0.5 (white-dominated) | 0.027 | 0.044 | 0.067 |
   | C = 1 | 0.042 | 0.070 | 0.111 |
   | **C = 3 (used here)** | 0.066 | 0.115 | **0.199** |
   | C = 10 | 0.083 | 0.151 | 0.282 |

   **If the array is white-noise-limited rather than confusion-limited, M2 collapses to
   ~0.03–0.07 dex and M1 becomes the engine again.** The whole conclusion hinges on C.
2. **My population is very top-heavy.** At σ_logh = 0.6 the loudest single source holds **32 %**
   of band power and the top 16 hold **83 %**. A real GWB is built from many more faint
   sources, so f is likely smaller and M2 is likely overstated. **Measuring C and f from
   GLACIER/CHORUS is now the highest-value follow-up** — it decides which mechanism the
   programme should be building.
3. **Subtraction is assumed perfect.** Real CW subtraction leaves residuals from parameter
   error, and imperfect subtraction of a loud source can inject more power than it removes.
   There is no subtraction-error term here at all.
4. **Confusion is modelled as a flat, isotropic addition to the noise**, with no spectral shape
   and no cross-pulsar correlation — whereas real confusion from a few loud sources is neither.
5. **Everything from VERNIER-1/2 still applies**: Asimov (no noise realisations, so no
   false-alarm rate and no flicker), caricatured red noise, no timing-model marginalisation,
   detection by threshold rather than by search, and the marginalisation penalty
   ρ²_P,marg = max(0, ρ_P² − 2 ln K_a) is a trials-factor argument rather than a computed
   marginal likelihood.

## One-line verdict

**Subtraction, not certification, is the wheel's engine** (0.14–0.19 dex vs 0.09–0.15, and
superadditive together at ~0.31 dex); there is **no critical N** — f_cert ramps smoothly from
0.21 at N = 1 to 0.94 at N = 16; and the wheel **self-terminates at rank ≈ 20 regardless of
seed** because the **absolute detection floor**, not the brightness gap, is what binds — so the
next thing to measure is the confusion-to-white ratio C, which decides whether M2 is real.
