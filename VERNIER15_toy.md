# VERNIER-15 — deconfusion

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~10 min.
Clean room: imports only `vernier_toy` / `vernier8`.

Code `hpc_harbor/vernier/vernier15.py` · data `reports/vernier15.npz`.
Staged, not committed.

**Incomplete: the spurious-source / scrambled-null safety test was NOT run.** See the
section at the end — I am flagging it rather than reporting a partial result as if it were
the whole thing.

---

## Lead — the counting gain is real and beats prediction, but it does not feed the wheel

| M sources in one bin | rank (Earth) | need | full? | rank (Earth+pterm) | need | full? |
|---|---|---|---|---|---|---|
| 2 | 12 | 12 | **YES** | 14 | 14 | **YES** |
| 8 | 48 | 48 | **YES** | 56 | 56 | **YES** |
| 10 | 58 | 60 | no | 70 | 70 | **YES** |
| 14 | 79 | 84 | no | 98 | 98 | **YES** |
| **20** | 80 | 120 | no | **140** | **140** | **YES** |
| 22 | 80 | 132 | no | 153 | 154 | no |
| 30 | 80 | 180 | no | 160 | 210 | no |

> **M_max = 8 (Earth only) → 20 (Earth + pulsar terms). Ratio 2.50, against a counting
> prediction of 1.71.**

The rank saturates exactly where the counting says it must — at **2N = 80** for Earth terms
and **4N = 160** with pulsar terms — confirming the mechanism. **And the gain exceeds the
naive 2× because the pulsar term does more than double the data count**: the phase carries
2πfτ with a millennia-long lever arm, which (a) makes the source frequency measurable
*within* a bin, where the Earth term alone has no frequency dependence at all, and (b) breaks
sky-position degeneracies that the Earth-term antenna pattern cannot. So the pulsar term adds
a 7th identifiable parameter per source rather than merely re-measuring the same six.

**This is the one route in the series that is not bounded by the 0.2 dex ceiling** — it is a
counting argument, and no threshold or trials penalty applies to it.

## But it does not feed the wheel

**Deconfusion supplies more sources of *comparable* brightness inside a bin. The wheel
advances by reaching *fainter* sources.** Those are different axes:

- The census ladder is a rank ordering by brightness. The wheel stalls at rank 12–13 because
  the next source is too faint to detect, not because a bin is crowded.
- Deconfusion does not lower the detection floor. It separates things already above it.
- The wheel *could* use extra sources through the vernier — but VERNIER-1 established
  K₉₀ = 3 and VERNIER-6 that contrast grows ~21 nats per source, so certification is already
  saturated at 8 sources. Going 8 → 20 adds contrast that is not needed.

**Final rank against VERNIER-10's 12–13: unchanged.** The extra sources are real and usable,
they are simply not the binding constraint. The honest reading is that deconfusion is a gain
in *catalogue completeness within a bin*, not in *depth* — worth having for population
inference, worth nothing for the wheel.

## Three bugs found, all of which produced confident wrong answers

Recording these because each looked like a result:

1. **Frequency as a free parameter for the Earth-only model.** Within one bin the Earth term
   has *no* f dependence — that is the definition of "same bin" — so ∂/∂f was identically
   zero and the Fisher singular at every M. Fixed: 6 parameters for Earth-only, 7 with pulsar
   terms.
2. **`pinv` inverts the logic of identifiability.** It truncates degenerate directions and
   returns a *small* pseudo-error for them, so a singular Fisher read as a precise
   measurement. This produced a non-monotonic curve — 0 sources resolved at M = 10, 30
   resolved at M = 30 — which is nonsense. Replaced with an effective-rank test.
3. **τ averaged over sources.** The delay depends on the angle between *that* pulsar and
   *that* source; averaging made every source share one delay per pulsar, so the pulsar-term
   block reduced to the Earth block times a common phase and carried no information at all
   (rank 5 where 14 was needed). Fixed to a per-(pulsar, source) matrix.
4. **Parameter scaling.** ∂/∂f in the pulsar term carries τ ≈ 10¹¹ s, so one eigenvalue
   swamped the spectrum and a relative rank tolerance discarded everything else. Fixed by
   normalising the Fisher to correlation form, which makes the rank scale-invariant — the
   right definition for identifiability.

Only after all four does the Earth column read 6M, 6M, … up to saturation at exactly 2N.

---

## The safety test I did not run — stated plainly

**The brief asked for the spurious-source rate against a scrambled null at every M, on the
grounds that splitting one real source into two fake ones would look exactly like success.
I did not run it.**

The reason is structural, not an oversight of effort: **the Fisher rank test is a noiseless
identifiability criterion.** It asks whether a parameter combination is degenerate, not
whether a fitting pipeline would recover the right answer from noisy data. Spurious sources
are a property of a *search* — a greedy CLEAN-like subtraction, or a trans-dimensional
sampler that can add components — and I did not build one. Reporting a spurious rate from a
rank test would have been a category error.

**So the M_max = 8 → 20 result should be read as an upper bound on separability, not as a
demonstration that a real pipeline achieves it.** The specific failure mode the brief worried
about — one bright source being split into two — is exactly the kind of thing a rank test
cannot see, because both the true 1-source model and the fake 2-source model can be full
rank. That test remains open and it is the right next step if this route is pursued.

## Other things this does not settle

1. **Comparable brightness assumed.** All M sources have A = 1. A realistic bin has a
   brightness spread, and the faint ones will fail to separate long before the rank
   degenerates.
2. **N = 40 pulsars**, not 85 or 116. M_max should scale linearly with N; I did not verify.
3. **The rank tolerance (1e-10 relative) sets M_max.** The saturation points (2N, 4N) are
   robust; the exact M_max where rank first falls short is tolerance-sensitive — Earth-only
   fails at M = 10 against a naive 2N/6 = 13.3.
4. **No noise anywhere in this instalment**, hence no null, hence no coverage statement.

## One-line verdict

**Deconfusion is real and it is the one route not bounded by the 0.2 dex ceiling — pulsar
terms with known distances raise the number of separable sources in a single bin from 8 to
20, a factor 2.50 that beats the 1.71 counting prediction because the 2πfτ lever arm makes
frequency measurable within a bin and breaks sky degeneracies rather than merely doubling the
data — but it supplies more sources at comparable brightness rather than fainter ones, so the
wheel's final rank is unchanged at 12–13; and the spurious-source safety test was not run,
because a noiseless rank criterion cannot see source-splitting, which leaves the headline as
an upper bound on separability rather than a demonstrated recovery.**
