# Differential test pack — wheel demo ↔ JUG real-likelihood machine

Generated from the **corrected** machine (full-δμ distribution, F8 adopted).
Purpose: two independent implementations of the same physics, compared at matched
points. **Any disagreement outside tolerance is a finding, reported not
reconciled.**

## Comparison rule

Agreement within tolerance at **all five** configurations ⇒ the two machines
agree at matched points. This establishes *consistency*, not truth. Any single
disagreement outside tolerance ⇒ STOP-and-report, with the disagreeing quantity,
both values, and both implementations' definitions of it.

## Conventions (must match before any comparison)

| quantity | convention |
|---|---|
| κ (comb spacing) | κ = 2πf(1−cos μ)/c, units 1/pc |
| tooth spacing | λ = 2π/κ = c/[f(1−cos μ)] |
| σ_φ | radians, wrapped-normal |
| attenuation | R = ⟨cos δ⟩ = e^(−σ²/2) is **amplitude**; e^(−σ²) is its square (power) |
| evidence | **nats**, natural log throughout |
| K | K = 6σ_d/λ, consistent with a ±3σ_d prior window |
| c (pc/s) | 9.7159e−09 |

## C1 — single clean CW, white noise only

- one pulsar; one source at f = 20 nHz; (1−cos μ) = 0.5 ⇒ κ = 3.234 pc⁻¹, λ = 1.943 pc
- σ_φ = 0.05 rad (white noise only, no confusion)
- σ_d (prior) = 2.0 pc; L_true = +0.37 pc

| output | this machine | tolerance |
|---|---|---|
| σ_L | **0.00773 pc** | ±10 % |
| σ_L / (σ_φ/κ) | **1.0000** | ±0.10 |

## C2 — the fringe scan over one prior window

- as C1, but report the full likelihood over L ∈ [−3σ_d, +3σ_d]

| output | this machine | tolerance |
|---|---|---|
| number of teeth in the window (K) | **6.177** | ±0.05 |
| spacing between adjacent maxima | **1.943 pc** | ±1 % |
| ln L at truth − ln L at the adjacent tooth | **0** (degenerate at one source) | exact |

## C3 — two-source stack

- two sources sharing one pulsar, κ₁ = 3.234, κ₂ = 4.690 pc⁻¹ (incommensurate)
- σ_φ = 0.6 rad each; σ_d = 2.0 pc

| output | this machine | tolerance |
|---|---|---|
| margin, truth over best alias | **report from fig40_ladder.csv at s = 0.75: 7.89 nats** | ±1.5 nats |
| argmax on the true tooth | **yes** | must match |

## C4 — one confused pair sharing an element

- two sources in the same 1/T element; s = 0.5 ⇒ σ_φ = √((1−s)/s) = 1.0 rad
- otherwise as C3

| output | this machine | tolerance |
|---|---|---|
| fraction on the true tooth (24 realisations) | **0.79** | ±0.10 |
| median margin | **1.42 nats** | ±0.5 nats |

## C5 — the zero-signal null

- as C2 but with the source phases drawn uniform (no signal)

| output | this machine | tolerance |
|---|---|---|
| decode hit rate | **0.0427** | ±0.010 |
| predicted reference (prior mass in one tooth, 2·TOOTH_W/window) | **0.0417** | — |
| **NOTE** | the null reference is **prior mass, NOT 1/K** (=0.125). A machine reporting ~1/K is using a flat prior. | |

## Values CSV

`differential_test_pack.csv` carries the same numbers machine-readably.
