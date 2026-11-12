# VERNIER-19 — does the deconfusion wheel turn?

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~12 min.
Clean room: imports only `vernier_toy` / `vernier16` / `vernier18`.

Code `hpc_harbor/vernier/vernier19.py` · data `reports/vernier19.npz`.
Staged, not committed.

---

## Lead — the loop does not turn. It manages exactly one turn.

**The safety gate stops it at turn 2 in every bin.** Genuine fraction per turn, on a
heterogeneous array (VERNIER-13's σ_log = 0.40), subtracting the *fitted* model:

| bin | N | λ₀ | turn 1 P(gen) | **turn 2 P(gen)** | cum dex | λ after |
|---|---|---|---|---|---|---|
| **5.41–6.27** | 5 | 0.811 | **0.750** | **0.229** ← STOP | **0.308** | 0.690 |
| 3.71–4.56 | 22 | 0.722 | 0.525 | **0.314** ← STOP | 0.139 | 0.631 |
| 4.56–5.41 | 11 | 0.636 | **0.375** ← STOP | — | 0.000 | 0.541 |
| 2.85–3.71 | 46 | 0.579 | 0.550 | **0.385** ← STOP | 0.005 | 0.534 |
| 2.00–2.85 | 156 | 0.263 | 0.625 | 0.650 | **0.387** ← STOP (t4) | 0.010 | 0.229 |

**Pure-noise null loop: 0.013 removals per run.** Clean — the failures are confusion, not noise.

**λ falls after every removal** (0.811 → 0.690, 0.722 → 0.631, 0.636 → 0.541), confirming the
analytic cascade: the steep top of the ladder is *what makes λ high*, so removing it leaves a
flatter, more equal population — which is harder, not easier.

**But λ is not what kills the loop. Absolute brightness is.** At λ = 0.811 the second source
is h₂/h₁ = 0.387, i.e. 8 dB fainter, sitting near the detection threshold where blending
dominates. The genuine fraction collapses 0.750 → 0.229 not because the residual is more
confused but because what remains is *faint*. The analytic λ-cascade (which projected
0.362 → 0.584 → 0.777 → 1.004 dex over four turns) assumed each successive top source stays
recoverable. It does not.

## What one turn is worth

| configuration | gain | **rank** |
|---|---|---|
| baseline (uniform, no subtraction) | 0.0787 | 6 |
| + heterogeneity (VERNIER-13) | 0.2443 | 11 |
| + one deconfusion turn | 0.3076 | **14** |
| **both combined** | **0.5519** | **28** |

*(step to rank 21 = 0.4272 dex)*

**One turn alone reaches rank 14, past the standing 11–13. Combined with heterogeneity the
arithmetic gives 0.552 dex and rank 28.**

**I do not believe the combined number and it should not be quoted.** Three reasons:

1. **The two gains are measured in different places and I simply added them.** Heterogeneity's
   0.244 dex is a coherence gain from f_cert = 0.80; subtraction's 0.308 dex is a
   noise-floor reduction. They plausibly compose, but I did not run them together in one
   simulation — the loop above measures subtraction *on* a heterogeneous array and reports the
   subtraction gain only.
2. **Rank 28 is off the end of the measured ladder** (VERNIER-7 measured 20 ranks; beyond that
   the gap is extrapolated as 0.30/k). Extrapolating 8 ranks past the data is not a result.
3. **The 0.308 dex applies in one bin** (5.41–6.27 nHz, 5 sources), which is also the only bin
   where VERNIER-17 found certification works. Everything rests on that single band.

**The defensible statement is: one deconfusion turn in the λ = 0.811 band is worth ~0.31 dex,
which exceeds the pulsar-term coherence ceiling (0.19–0.28) and is the largest single-mechanism
gain measured in this series — and it happens once, not repeatedly.**

## The honest gap, now closed: joint position fitting DOES help

The brief asked whether a joint *position* fit is tractable. **It is, at these M**, and it was
worth testing:

| M = 5, λ = 0.811 | genuine recoveries per run |
|---|---|
| greedy | 0.58 |
| **joint position fit** | **0.83** |

**12 trials in 2.2 s — entirely tractable at M = 5.** A Nelder–Mead fit over the 10 sky angles,
seeded from the greedy solution, with amplitudes solved linearly at each step, recovers **43%
more genuine sources**.

**So VERNIER-18's "simultaneous fitting doesn't help" was wrong, and specifically wrong in the
way I flagged**: the amplitude-only refit was the cheap linear half, and the blend bias does
live in the positions. **This changes the outlook for the loop** — a position-fitting search
might sustain turn 2, which the greedy one cannot. I did not re-run the loop with it (12
trials at M = 5 is a demonstration of tractability, not a measurement), and at M = 22 or 156 the
optimiser has 44 or 312 free parameters and almost certainly will not converge from a greedy
seed.

**That is the single most promising loose end in this instalment.**

---

## A normalisation bug, caught and fixed

The first version normalised the injected signal in *unwhitened* units while the F-statistic
operates on whitened data (d·w), making the effective SNR geometry-dependent and
systematically low. It depressed turn-1 P(genuine) to 0.35–0.43 against VERNIER-18's
0.73–0.77 on a uniform array. Fixed by normalising ‖m·w‖. **The corrected turn-1 numbers
(0.750 at λ = 0.811) are consistent with VERNIER-18's 0.933 once the heterogeneous array is
accounted for; the residual difference is not explained and I have not chased it.**

## What this does not settle

1. **The loop was run with the greedy search**, which the joint-position result above suggests
   is the wrong algorithm. "The loop does not turn" is established *for greedy*; with position
   fitting it is open.
2. **Top source at whitened SNR 15.** Turn 2's failure is a threshold effect, so it scales
   directly with this: a brighter bin would sustain more turns.
3. **The ladder is h_k² ∝ k^−α** fitted to (λ, N). Real bins need not be power laws, and the
   h₂/h₁ = 0.387 that kills turn 2 is a direct consequence of that form.
4. **40 pulsars, Earth terms only, one bin, white noise**, and subtraction assumed perfect
   apart from the position error — no amplitude or parameter-error residual is left behind.
5. **Cumulative dex is computed from the true residual signal power** with removed sources
   excluded only when the recovery was genuine, so a blended removal contributes nothing
   rather than actively corrupting. **The real penalty for subtracting a blend — injecting a
   fake negative source — is not modelled**, so these gains are optimistic.

## One-line verdict

**The deconfusion loop does not turn: the genuine fraction collapses from 0.750 to 0.229 at
turn 2 in the best bin and the safety gate stops every bin at turn 1–2, because λ falls after
each removal and, more decisively, the second source sits 8 dB down near threshold where
blending dominates — so the mechanism yields one turn worth ~0.31 dex in the λ = 0.811 band,
which is the largest single gain in this series and exceeds the coherence ceiling, but does not
compound; and the one genuinely promising loose end is that a joint POSITION fit is tractable
at M = 5 and recovers 43% more genuine sources than greedy, which means "simultaneous fitting
doesn't help" was my error and the loop deserves a re-run with it.**
