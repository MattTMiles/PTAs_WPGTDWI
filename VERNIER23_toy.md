# VERNIER-23 — drop the certification requirement

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~10 min.
Clean room: imports only `vernier_toy` / `vernier15` / `vernier4`.

Code `hpc_harbor/vernier/vernier23.py` · data `reports/vernier23_part1.npz`.
Staged, not committed.

---

## Lead — you were right, and it does not save the route

**Part 1: the separability gain survives marginalisation. Certification was never required.**

| M | Earth only | distances KNOWN | **marginalised (10 pc prior)** | marginalised (100 pc) |
|---|---|---|---|---|
| 8 | full | full | full | full |
| 12 | — | full | **full** | full |
| 17 | — | full | **full** | — |
| 20 | — | — | **full** | — |
| **M_max** | **8** | **17** | **17–20** | **16** |

*(counting prediction: Earth 2N/6 = 13, known 4N/7 = 23, marginalised 3N/7 = 17 — the
marginalised case is confirmed at 17–20.)*

**Marginalising the N_psr distances costs N_psr degrees of freedom out of 4N_psr data, leaving
3N_psr for sources.** The separable limit goes 8 → 17–20, essentially indistinguishable from
the known-distance case. **The per-pulsar ρ_P ≥ 2 bar never enters — it is an identifiability
question, not a per-pulsar detection question.** Your diagnosis was correct.

**Part 2 kills it anyway, for a different reason.**

The mechanism needs the pulsar-term line complex to be resolvably broad. The τ spread across a
0.3–5 kpc array is large (median 3302 yr, 5–95 % spread 5.6×10¹¹ s), but the *frequency* spread
it produces is not:

| f | Mc = 10^8.5 | 10^9.0 | **10^9.5** | resolvable (spread × T > 1)? |
|---|---|---|---|---|
| 3.0 nHz | 0.0 | 0.0 | 0.1 | no |
| **5.8 nHz** | 0.0 | 0.1 | **0.6** | **no** |
| 10 nHz | 0.1 | 0.7 | **2.4** | yes (heavy only) |
| 20 nHz | **1.4** | 4.6 | 7.2 | yes |

**At the operating band — 5.41–6.27 nHz, where every other mechanism works — the spread is
0.6 of a Fourier bin even for a 10^9.5 M☉ binary. The lever does not exist there.**

## The scissors, stated precisely

| f | PHOENIX (GLACIER) sources/bin | lever exists? |
|---|---|---|
| 2.4–5.8 nHz | **156, 46, 22, 11, 5** | **no** |
| 10–20 nHz | **1** | **yes** |
| above 32 | empty | yes |

**Where the bins are crowded the lever is absent; where the lever exists the bins already hold
one source. The mechanism and the problem it would solve do not occupy the same frequencies.**

**Part 3 was not run**, because its premise fails: there is nothing to separate at the
frequencies where separation would work. Running a subtraction chain on 1-source bins would
produce a number with no content.

## The one thing that could reopen it

**The two censuses disagree by 3–4 orders of magnitude exactly where it matters:**

| f | PHOENIX measured | VERNIER-4's f^(−11/3) model |
|---|---|---|
| 10 nHz | 1 | **6.9×10⁴** |
| 20 nHz | 1 | **5.4×10³** |
| 32 nHz | empty | 965 |

**If PHOENIX's measured GLACIER census is the operative one, this route is closed. If
VERNIER-4's theoretical population is closer to the truth, the lever and the confusion coexist
at 10–20 nHz and the route is open.** That disagreement — not the marginalisation question,
which is now settled positively — decides it, and this toy cannot adjudicate between a measured
catalogue and a population model.

## Spin-off worth having: chirp mass from the line width

The same calculation gives a **chirp-mass measurement with all distances marginalised**:

| f | Mc | line width | **σ(log₁₀ Mc)** |
|---|---|---|---|
| 5.8 nHz | 10^9.0 | 1.9e-10 Hz | **0.102 dex** |
| 5.8 nHz | 10^9.5 | 9.7e-10 Hz | **0.029 dex** |
| 20 nHz | 10^9.0 | 7.3e-9 Hz | **0.007 dex** |
| 20 nHz | 10^9.5 | 1.1e-8 Hz | 0.047 dex |

**Against the banked chirp-mass posterior, which is prior-dominated at σ = 0.50–0.87 dex, this
is a 5–100× improvement on a parameter the programme has never been able to measure.** It does
not need certified distances and it does not need the bins to be crowded — it needs one
detected source and the array's existing distance priors.

**Caveat: the 5.8 nHz rows use a width that Part 2 says is only 0.6 of a Fourier bin, i.e.
formally unresolved.** The Fisher error is still finite there because a sub-bin width is
measurable at high SNR, but it is an extrapolation past the resolution criterion and should be
checked before it is relied on. The 20 nHz rows sit in the resolvable regime.

---

## What this does not settle

1. **Rank is identifiability, not precision.** M_max = 17–20 says the marginalised Fisher is
   invertible, not that the parameters are *well* measured. VERNIER-16 showed a greedy search
   attains nothing like the rank limit.
2. **The 10 pc vs 100 pc prior changes M_max by 4** (20 → 16), so the result depends on the
   distance priors more than a pure counting argument suggests.
3. **f = 10 nHz was used for Part 1's Fisher**, which is in the regime where the lever exists.
   At 5.8 nHz the marginalised rank would degrade toward the Earth-only value — I did not
   measure that, and it is the number that would confirm Part 2's verdict from the Fisher side.
4. **40 pulsars, one bin, no noise realisation** — Part 1 is a noiseless rank calculation
   throughout.

## One-line verdict

**The certification requirement really was unnecessary — marginalising the pulsar distances
over their existing priors preserves the separability gain almost intact (8 → 17–20 against
17 with distances known), exactly as the counting argument 3N/7 predicts, so the per-pulsar
ρ_P ≥ 2 bar that killed the three previous routes never enters this one — but the route closes
anyway on frequency, because the pulsar-term line complex is only resolvably broad above
~10 nHz while PHOENIX's census puts every crowded bin below 6 nHz, so the lever and the
confusion never occupy the same band; the one thing that could reopen it is the 3–4 order
disagreement between PHOENIX's measured occupancy and VERNIER-4's f^(−11/3) model at 10–20 nHz,
which this toy cannot adjudicate; and the calculation yields a spin-off worth more than the
route itself — chirp mass to 0.03–0.10 dex from the line width alone, with all distances
marginalised, against a banked posterior that is prior-dominated at 0.50–0.87 dex.**
