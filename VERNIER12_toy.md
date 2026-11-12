# VERNIER-12 — pairwise, properly

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~12 min.
Clean room: imports only `vernier_toy` / `vernier6` / `vernier6b` / `vernier7` /
`vernier8` / `vernier9`.

Code `hpc_harbor/vernier/vernier12.py` · data `reports/vernier12_{mech,group,wheel}.npz`.
Staged, not committed.

---

## Bookkeeping first — and Q2 matters

**Q1. What was the "choice"?** It was **σ_f/f, the assumed fractional uncertainty in the
shared source frequency**, fixed at 1e-4 in VERNIER-6. That single number is what the 105
nats hinges on. Re-measured here across its range:

| σ_f/f | 0 | 1e-6 | 1e-5 | 3e-5 | **1e-4** | 3e-4 | 1e-3 |
|---|---|---|---|---|---|---|---|
| **joint − per-pulsar (nats)** | **0.00** | 0.02 | 1.82 | 15.5 | **100.5** | **112.9** | 10.5 |

(The fall at 1e-3 is both arms being destroyed, not the effect switching off.)

**Q2. Was pairwise active in 6b/7/8/9/10/11? No — it silently dropped out.**
`grep sigma_phi` returns **0 matches** in `vernier6b.py`, `vernier7.py`, `vernier9.py`. When
VERNIER-6b rewrote the likelihood for soft weights it never carried the smear term across,
and everything downstream inherited that.

**But the direction is the opposite of the worry.** Those runs were not "missing pairwise"
(pessimistic) — they were running with **no shared-parameter smear at all**, i.e. source
parameters known *exactly*, which is better than any joint treatment can achieve. **So every
wheel result from VERNIER-6b onward — including rank 13 — is an upper bound that silently
assumed perfect source parameters.** That needed saying before anything else, and it is now
quantified below.

## Mechanism test — passes exactly

The prediction was that the benefit must vanish when source parameters are known.

> **Gain at σ_f = 0 is 0.000 nats.** Exactly zero, to the printed precision.

It then rises monotonically with σ_f and reproduces VERNIER-6's figure (100.5 here vs 97.5
there). **The 105 nats is real and it is entirely the cancellation of correlated
shared-parameter smear** — nothing else is contributing.

**A realistic σ_f/f is right where the effect is largest.** For a CW of SNR ρ measured over T,
σ_f ≈ 1/(2πTρ):

| ρ | 8 | 20 | 50 |
|---|---|---|---|
| σ_f/f at 25 nHz | 8.4e-4 | **3.4e-4** | 1.3e-4 |

**At σ_f/f = 3e-4 the per-pulsar treatment gives contrast 0.00 — it does not work at all —
while the full joint gives 112.9 nats. Joint treatment is not an enhancement here, it is a
requirement.**

## "Pairwise" is a misnomer — pairs buy nothing

| group size m | 1 | **2** | 3 | 5 | 10 | 20 | 50 |
|---|---|---|---|---|---|---|---|
| σ_φ residual | 4.85 | 3.43 | 2.80 | 2.17 | 1.53 | 1.09 | 0.69 |
| contrast (nats) | 0.00 | **0.40** | 2.81 | 13.5 | 43.7 | 78.6 | **111.9** |
| **% of full-joint gain** | 0% | **0%** | 3% | 12% | 39% | 70% | 100% |

**Pairs capture 0% of the gain. Triples 3%. You need m ≈ 10 for a third of it and the full
array for all of it.**

The cutoff is analytic: the smear is σ_φ/√m, so the gain switches on only once
**m ≳ σ_φ² ≈ 24**. Below that the residual phase error still exceeds a radian and the comb is
still decohered. **There is a natural cutoff and it is the full array, not a pair.**

If Xiao et al.'s contribution is genuinely *pairwise* in the sense of two pulsars at a time,
this toy says that specific construction should not deliver the gain — the mechanism needs
the shared parameter constrained by many pulsars simultaneously. I have not read the paper
and cannot say whether "pairwise" there means the same thing.

## The wheel, with the smear finally present

At the VERNIER-9 working point (f_E = 20 nHz, ρ_P = 6, per-pulsar distances, 40 pulsars):

| treatment | **final rank** | w_true signal | w_true null |
|---|---|---|---|
| σ_f = 0 — *what 6b–10 assumed* | **10** | 1.0000 | 0.0000 |
| realistic 3e-4, **per-pulsar** | **5** | **0.0101** | 0.0061 |
| realistic 3e-4, **pairs (m=2)** | **6** | 0.0808 | 0.0058 |
| realistic 3e-4, **full joint** | **9** | 1.0000 | 0.0000 |

**Reading:**
- Adding realistic source-parameter uncertainty costs **5 ranks** (10 → 5) if handled
  per-pulsar. At that point w_true = 0.0101 against a null of 0.0061 — certification has
  essentially collapsed to chance.
- **Full joint treatment recovers 4 of those 5 ranks** (5 → 9) and restores w_true to 1.0000
  against a null of 0.0000.
- **Pairs recover one rank.**

So the net correction to the series is **one rank**: the wheel reaches rank 9 with realistic
source parameters properly handled, versus rank 10 with them assumed perfect. **The earlier
results were optimistic, but only slightly — provided the joint treatment is actually used.**
Without it they were optimistic by five ranks.

*(VERNIER-10's rank 13 is not directly comparable: that was configuration D with a threshold
lowered by known sky and frequency. This run uses the blind threshold of 8.)*

---

## What this does not settle

1. **σ_f ≈ 1/(2πTρ) is a Fisher estimate**, not a measured posterior width, and the whole
   effect scales with it. A factor 3 in σ_f moves the per-pulsar arm between "works" and
   "contrast zero".
2. **Only the frequency error is modelled.** Sky position and chirp mass are also shared and
   also smear correlated-ly; including them would make the per-pulsar arm worse still and the
   joint requirement stronger.
3. **The smear is modelled as a decoherence factor exp(−σ_φ²/2) on the envelope** — a
   Gaussian-marginalisation caricature, not an actual joint fit over a shared hyperparameter.
   A real joint fit would also pay a degrees-of-freedom cost this does not include.
4. **`group_m` is imposed, not achieved.** I assume fitting m pulsars jointly constrains the
   shared parameter m-fold. That is the right scaling for an unbiased joint estimator but
   assumes the group's pulsars actually carry independent information about f.
5. **Still no GWB, sources pre-resolved, white noise.** Same optimistic frame as VERNIER-6
   onward.

## One-line verdict

**The "choice" was σ_f/f = 1e-4; the smear term silently vanished from every wheel after
VERNIER-6, which made those runs optimistic rather than pessimistic because they implicitly
assumed exact source parameters; the mechanism test passes exactly (gain = 0.000 nats when
σ_f = 0), confirming the 105 nats is pure shared-parameter cancellation; but "pairwise" is
the wrong name — pairs capture 0% and triples 3%, with the gain switching on only at
m ≳ σ_φ² ≈ 24 pulsars fitted jointly — and with realistic σ_f/f = 3e-4 the per-pulsar wheel
collapses from rank 10 to rank 5 while the full-joint wheel holds rank 9, so joint treatment
is a requirement, not an enhancement, and the net correction to the series is one rank.**
