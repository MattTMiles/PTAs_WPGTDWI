# VERNIER-13 — heterogeneous arrays

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~5 min.
Clean room: imports only `vernier8` for the measured census ladder.

Staged, not committed.

---

## Lead — heterogeneity is what makes the wheel work at all

**At matched total array sensitivity, a uniform array certifies nothing and a realistic one
certifies most of its power.**

| array | ρ_array | **f_cert** | pulsars certifying | top-5 share | top-10 | top-20 |
|---|---|---|---|---|---|---|
| **uniform** | 8 | **0.0000** | **0.0** | 4.3% | 8.6% | 17.2% |
| heterogeneous σ=0.25 | 8 | 0.1168 | 1.3 | 27.3% | 40.4% | 57.3% |
| **heterogeneous σ=0.40** | 8 | **0.4001** | **3.0** | 49.3% | 64.2% | 78.9% |
| heterogeneous σ=0.50 | 8 | 0.5544 | 3.3 | 62.8% | 76.4% | 88.0% |
| **uniform** | 20 | **0.0000** | **0.0** | 4.3% | 8.6% | 17.2% |
| **heterogeneous σ=0.40** | 20 | **0.7981** | 20.0 | 49.3% | 64.3% | 79.0% |

**Three pulsars certify and they carry 40% of the array's pulsar-term power.** That is exactly
the asymmetry the brief predicted, and it is the number that matters: not the count, the
power share.

**The scissors is an artefact of assuming identical pulsars.** With a uniform array at
detection threshold, every pulsar sits at ρ_P = 8/√116 = 0.743 against a bar of 2 — *all of
them* miss, by the same amount, and f_cert is exactly 0. Give the array a realistic spread at
*identical total sensitivity* and the best few clear the bar comfortably while the mean still
misses. **f_cert goes 0.000 → 0.400 with no change in the array's detection capability
whatsoever.**

**Weight distribution used**: log-normal in log₁₀(w), w ∝ 1/σ_RMS, with σ_log = 0.40 dex as
fiducial. Calibrated to reproduce the known behaviour — it gives top-10 = 64% of sensitivity
and N_eff = 1/Σp² = 15 of 116, consistent with "~10 pulsars dominate" and with 23-of-68 style
findings. NANOGrav-15yr residual RMS spans ~5e-8 to ~1e-5 s (≈2.3 dex), so 0.40 dex is
conservative. σ = 0.25 and 0.50 are shown throughout as sensitivity.

## Fed back into the wheel

| array | f_cert | gain (dex) | **final rank** |
|---|---|---|---|
| uniform, ρ_array = 8 | 0.000 | 0.0000 | **4** |
| heterogeneous 0.40, ρ_array = 8 | 0.393 | 0.0969 | **7** |
| heterogeneous 0.40, ρ_array = 20 | 0.797 | 0.1652 | **8** |
| heterogeneous 0.50, ρ_array = 20 | 0.853 | 0.1732 | **9** |

**Revising VERNIER-10's configuration D with a realistic f_cert** (its 0.079 dex trials
benefit plus the signal term):

| assumption | total gain | rank |
|---|---|---|
| f_cert = 1 — what VERNIER-10 implicitly assumed | 0.2715 | **12** |
| uniform array, f_cert = 0 | 0.0787 | **6** |
| **heterogeneous, f_cert = 0.80** | **0.2443** | **11** |

**So VERNIER-10's rank 13 was only reachable because it assumed all pulsar-term power was
usable. On a uniform array that assumption is not merely optimistic, it is impossible
(f_cert = 0, rank 6). On a realistic heterogeneous array it is nearly right — rank 11.**
Heterogeneity is what rescues the earlier number, not what threatens it.

## The design question — and it is a clean answer

At **fixed total observing effort** split N ways, each pulsar gets T/N, so σ_a ∝ √(N/T),
w_a ∝ √(T/N), and:

> **ρ_a ∝ N^(−1/2) while ρ_array = √(Σρ_a²) is INVARIANT in N.**
> **Detection does not care how you split the effort. Certification cares enormously.**

| N_psr | ρ_a (uniform) | f_cert uniform | f_cert het | gain | rank |
|---|---|---|---|---|---|
| 10 | 6.33 | 1.00 | 0.988 | 0.1913 | 9 |
| **24** | 4.08 | 1.00 | 0.960 | **0.1876** | **9** |
| 50 | 2.83 | 1.00 | 0.903 | 0.1801 | 9 |
| 116 | 1.86 | **0.00** | 0.800 | 0.1656 | 8 |
| 400 | 1.00 | 0.00 | 0.559 | 0.1276 | 7 |

> **A small array of excellent pulsars beats a large array of mediocre ones — for
> certification — at every N, and the optimum is the smallest array you are otherwise allowed
> to have.**

The floor is **N ≈ 24**, set by VERNIER-12: below that the shared source-parameter smear is
not cancelled (the joint fit needs m ≳ σ_φ² ≈ 24 pulsars). So **the certification optimum is
N ≈ 24 excellent pulsars**, gain 0.188 dex, rank 9 — versus rank 8 at 116 and rank 7 at 400.

Two further observations worth having:

- **Uniform arrays fail by a cliff, heterogeneous ones degrade gracefully.** A uniform array
  certifies *everything* while ρ_array/√N ≥ 2 (i.e. N ≤ 100 at ρ_array = 20) and *nothing*
  above it. A heterogeneous array still has f_cert = 0.56 at N = 400, because its best pulsars
  stay above the bar long after the mean has fallen below.
- **This inverts the usual PTA instinct.** More pulsars is the right answer for detection and
  for the GWB, and it is the wrong answer for distance certification at fixed effort. The two
  goals genuinely conflict, which is the practical content of the scissors.

---

## What this does not settle

1. **σ_log = 0.40 dex is calibrated, not measured.** I did not reconstruct actual IPTA
   per-pulsar RMS values. f_cert at ρ_array = 8 spans 0.12 → 0.55 across σ_log = 0.25 → 0.50,
   so the headline is qualitatively robust and quantitatively soft.
2. **The ρ_P ≥ 2 bar is VERNIER-1's**, from a toy with its own caveats (Asimov, trials-factor
   threshold rather than a calibrated FPR). f_cert is a step function of that bar, so it is
   the single most leveraged number here.
3. **Fixed "observing effort" assumes σ ∝ t_obs^(−1/2)** — radiometer-limited. Real pulsar
   timing hits jitter and red-noise floors where extra integration buys nothing, which would
   *strengthen* the small-array conclusion (you cannot fix a bad pulsar with more time).
4. **The weight distribution is assumed independent of sky position.** In reality the best
   pulsars are not isotropically placed, which matters for both the antenna coupling and the
   sky trials factor.
5. **No per-turn wheel run here** — this note works through f_cert analytically and feeds it
   into the established gain law rather than re-running the soft-weight wheel, so there is no
   new w_true/null comparison in this instalment. The null discipline from VERNIER-8 onward is
   unchanged for the runs it applies to.

## One-line verdict

**Heterogeneity is not a complication, it is the thing that makes certification possible: at
identical total sensitivity a uniform 116-pulsar array certifies zero power while a realistic
one certifies 40–80%, because three pulsars carrying 40% of the weight clear a bar the mean
misses — which means VERNIER-10's rank 13 was reachable only on a heterogeneous array (rank
11) and impossible on a uniform one (rank 6); and since ρ_array is invariant under splitting
effort while ρ_P ∝ N^(−1/2), a small array of excellent pulsars beats a large mediocre one for
certification at every N, with the optimum at N ≈ 24 — the floor imposed by the joint
source-parameter fit.**
