# VERNIER-21 — λ as a resolution parameter

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~10 min.
Clean room: imports only `vernier_toy` / `vernier16` / `vernier18` / `vernier20`.

Code `hpc_harbor/vernier/vernier21.py` · data `reports/vernier21_part1.npz`.
Staged, not committed.

---

## Lead — the hypothesis is right, and it still does not close the loop

**Part 1 confirms it, and by more than expected.** λ measured per *resolution element* rather
than per Fourier bin, from actual template overlaps:

| bin | N | λ₀ (raw) | **λ_eff Earth-only** | **λ_eff Earth+pterm** | array ceiling E | array ceiling E+P |
|---|---|---|---|---|---|---|
| 2.00–2.85 | 156 | 0.263 | **0.792** | **0.917** | 1.95 | 3.32 |
| 2.85–3.71 | 46 | 0.579 | 0.937 | 0.978 | 3.86 | 6.59 |
| 3.71–4.56 | 22 | 0.722 | 0.965 | **0.991** | 5.27 | 10.46 |
| 4.56–5.41 | 11 | 0.636 | 0.950 | 0.984 | 4.36 | 7.91 |
| 5.41–6.27 | 5 | 0.811 | 0.975 | **0.993** | 6.26 | 12.12 |

The pulsar term resolves frequency **1537× finer** than the Fourier bin (1/(2πτ) = 1.03e-12 Hz
vs 1/T = 1.58e-9 Hz at L = 1.5 kpc), and the measured λ_eff rises accordingly.

**But the larger surprise is the Earth-only column.** Even with no pulsar terms, λ_eff = 0.79–0.98
against raw λ of 0.26–0.81. **λ_fore as the field defines it — and as VERNIER-17 used it —
substantially overstates confusion**, because it counts all bin power as effective noise while
ignoring that interlopers at different sky positions barely project onto the target's template.
The pulsar term then roughly halves the residual leakage again.

**This is a material correction to VERNIER-17**, which concluded certification needs λ ≥ 0.80
and that most of the census fails. On λ_eff the whole census clears 0.79.

## Part 2 — it does not compound, and I nearly reported that it did

**I made an error and caught it.** My first Part 2 compared the ceiling √(λ_eff/(1−λ_eff))
against the certification bar of 2 and found every high-λ bin certifying at turn 1
(ceiling 6.03, 6.70). **That comparison is invalid.** λ_eff is computed from templates
*stacked over all pulsars*, so its ceiling is an **array-level** quantity. Certification
(VERNIER-1) is a **per-pulsar** test, and per-pulsar SNR ≈ array/√N_psr.

| bin | λ_eff (E+P) | array ceiling | **per-pulsar ceiling** | vs bar 2 |
|---|---|---|---|---|
| 2.00–2.85 | 0.917 | 3.32 | **0.524** | FAIL |
| 2.85–3.71 | 0.978 | 6.59 | 1.042 | FAIL |
| 3.71–4.56 | 0.991 | 10.44 | **1.650** | FAIL |
| 4.56–5.41 | 0.984 | 7.92 | 1.252 | FAIL |
| **5.41–6.27** | 0.993 | 12.09 | **1.911** | **FAIL** (just) |

**Every per-pulsar ceiling is below the bar.** The resolution gain is real and large at array
level; the √N_psr division puts certification out of reach. **The loop never starts, so there
is nothing to compound.**

This is the detection/certification scissors again, arriving through a third door. What array
would close it — needing array ceiling ≥ 2√N_psr:

| bin | array ceiling | certifies only if N_psr ≤ |
|---|---|---|
| 156 sources | 3.32 | **2.7** |
| 22 sources | 10.44 | **27.2** |
| 5 sources | 12.09 | **36.5** |

**A 22-source bin certifies only with ≤27 pulsars, a 5-source bin with ≤37.** That is the same
answer VERNIER-13 reached from the other side (optimum N ≈ 24), and it is consistent: small
arrays certify, large arrays detect.

**Caveat on Part 2 that matters:** my certification rule is binary (f_cert = 0 or 1), which
makes compounding impossible by construction. The 156-source bin sits at per-pulsar ceiling
0.524 and the 5-source bin at 1.911 — the latter is *just* below the bar, exactly where a
graded rule (VERNIER-6b's soft weights) might bootstrap. **"It does not compound" is
established for a hard bar and is untested for a soft one.**

## The untested corner — and it is the best configuration

| configuration | SNR | f_cert | gain | rank |
|---|---|---|---|---|
| A baseline: uniform, no sub, no pterm | 2.788 | 0.000 | 0.0000 | 4 |
| subtraction only, uniform | 4.826 | 0.000 | 0.2384 | 11 |
| **SUB + COHERENCE, UNIFORM** | **4.871** | 0.000 | **0.2424** | **11** |
| all three (heterogeneous) | 4.163 | 0.137 | 0.1742 | 9 |

**Subtraction + coherence on a uniform array is the best configuration measured: 0.2424 dex,
rank 11.** But coherence contributes only +0.004 dex over subtraction alone, because f_cert =
0.000 — no pulsar certifies, for the per-pulsar reason above. **The corner is worth testing and
it is essentially subtraction alone.** Heterogeneity remains the harmful ingredient (0.174 vs
0.242).

## Part 3 — not run

Whether a joint fit attains VERNIER-15's identifiability limit at M = 5/11/22 was **not
measured**. Part 2's failure removes the motivation — the resolution gain does not convert into
certification regardless of whether the limit is attainable — and I would rather say so than
report a number that changes nothing. VERNIER-20 established joint position fitting is
tractable at M = 5 (2.2 s) and gives +43% genuine recoveries; whether that reaches 20 separable
sources at M = 22 is open.

---

## What this does not settle

1. **λ_eff assumes the interloper leakage adds incoherently** (I sum h_i²·|overlap|²). Coherent
   addition would be worse; with many sources at random phases the incoherent sum is the right
   expectation, but it is an expectation.
2. **Overlaps are computed at fixed known source parameters.** A real search does not know the
   interlopers' sky positions, so it cannot exploit their orthogonality as cleanly as this
   calculation assumes. **λ_eff is therefore an upper bound on the achievable effective λ**, and
   the gap between it and raw λ_fore is where the real answer lies.
3. **Frequencies spread uniformly across one Fourier bin.** If sources cluster in frequency the
   resolution gain shrinks.
4. **The binary certification rule** (see Part 2 caveat) — the single most consequential
   modelling choice here.
5. **40 pulsars, one bin, white noise, no GWB.**

## One-line verdict

**λ is indeed a resolution parameter and not a sky parameter — measured per resolution element
it rises from 0.263–0.811 to 0.792–0.975 with Earth terms alone and 0.917–0.993 with certified
pulsar terms, which is a material correction to VERNIER-17's use of raw λ_fore — but the loop
still does not start, because λ_eff's ceiling is an ARRAY-level quantity while certification is
a PER-PULSAR test, and dividing by √N_psr leaves every bin below the bar (best case 1.911 at
5.41–6.27 nHz); the same arithmetic says a 22-source bin would certify with ≤27 pulsars and a
5-source bin with ≤37, agreeing with VERNIER-13's optimum of ≈24 from the opposite direction;
and the untested corner turns out to be the best configuration measured (subtraction +
coherence, uniform, 0.242 dex, rank 11) chiefly because it omits the harmful heterogeneity.**
