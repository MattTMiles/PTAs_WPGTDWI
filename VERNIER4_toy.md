# VERNIER-4 — realistic sky, generous array

**EXPLORATORY. NOT pre-registered.** Proof of principle, not a forecast. Toy model,
cronus, CPU/numpy, ~3 min. Clean room: imports only `vernier_toy` / `vernier2`.

Code `hpc_harbor/vernier/vernier4.py` · data `reports/vernier4.npz` ·
plot `reports/VERNIER4_realistic.png`. Staged, not committed.

---

## Lead

**The wheel does not survive a realistic population on any array we could build, and the
reason is not the one VERNIER-3 was worried about.**

- **Closure, optimistic confusion model: N_psr ≈ 10 000 at 10 ns** (in-band, SNR 9.8), or
  2 000 at 1 ns if you allow 113 nHz — above the band, where the population law is pure
  extrapolation.
- **Closure, pessimistic confusion model: never.** Max SNR of the loudest source is **0.32**
  at 50 000 pulsars and 1 ns timing. Not 0.32 short of threshold — 25× short.
- **M2 collapses on today's array, as suspected: 0.039 dex max** anywhere in the band at
  100 ns, versus the 0.199 dex VERNIER-3 assumed. **5× smaller. VERNIER-3's conclusion is
  inverted at 100 ns** — but it *recovers to 0.189 dex at 10 ns*, so M2 is array-dependent,
  not dead.
- **And a new blocker that neither previous session saw: M1 and detection pull in opposite
  directions in N_psr.** Certification is a *per-pulsar* statistic; detection is an *array*
  statistic. Building the array big enough to detect makes certification impossible.

**Which mechanism wins on a realistic sky: M2, but both are negligible.** At a closing array
(5 000 psr, 3 ns) the wheel manages 1–2 turns at 50–68 nHz with gains of 0.001–0.023 dex, and
**zero turns anywhere at or below 32 nHz**.

---

## Deliverable 1 — C(f) and f_rem(f), now derived rather than chosen

Census: `N_c(f) = 1000·(f/f_yr)^(−11/3)` per 1/T bin; strains from Euclidean volume
(`dN/dh ∝ h^−4` → `h_k² ∝ k^(−2/3)`); bright end normalised with **no free parameter** by
`Σ_k R_k² = 2 S_gwb(f)·Δf` with A_gwb = 2.4e-15. Only the top 300 per bin are handled
individually; the rest is a Gaussian sea. Cost is independent of N_c.

| f nHz | N_c / bin | C(100 ns) | C(30 ns) | C(10 ns) | f_rem(300) | R₁ ns |
|---|---|---|---|---|---|---|
| 2 | 2.5e7 | 1.0e4 | 1.1e5 | 1.0e6 | **0.023** | 34.1 |
| 5 | 8.7e5 | 189 | 2.1e3 | 1.9e4 | 0.070 | 8.20 |
| 10 | 6.9e4 | 9.36 | 104 | 936 | 0.164 | 2.79 |
| 15 | 1.6e4 | **1.61** | 17.9 | 161 | 0.268 | 1.49 |
| 20 | 5.4e3 | 0.464 | 5.16 | 46.4 | 0.381 | 0.949 |
| 32 | 965 | 0.061 | **0.673** | 6.06 | 0.678 | 0.457 |
| 80 | 33.5 | 0.001 | 0.013 | **0.114** | 1.000 | 0.110 |

**The scissors.** At fixed timing precision C ∝ f^(−13/3) falls while f_rem ∝ N_c^(−1/3) ∝
f^(+11/9) rises. Their product — which is what M2 gain depends on — goes as **f^(−28/9)**.
So M2 has exactly one operating window and it sits at the **GWB/white-noise crossover**:
below it the sea is enormous but 10⁵–10⁷ sources deep and therefore irreducible; above it a
few sources hold the power but there is no confusion left to remove.

**The crossover is not fixed** — better timing lowers S_white, raising C everywhere and
pushing the window up to where f_rem is large:

| σ_toa | crossover (C = 1) | peak M2 gain (subtract 300) | at f |
|---|---|---|---|
| 100 ns | 16.8 nHz | **0.0393 dex** | 15 nHz |
| 30 ns | 29.1 nHz | 0.0836 dex | 20 nHz |
| 10 ns | 48.4 nHz | **0.1891 dex** | 32 nHz |

So the honest statement is **not** "M2 is an artefact". It is: **M2 is worth 0.04 dex on a
100 ns array and 0.19 dex on a 10 ns array, and its operating frequency climbs from 15 to
32 nHz as timing improves.** VERNIER-3's assumed 0.199 dex was accidentally the *10 ns*
answer applied to a 100 ns array.

## Question 2 — where in the band

**Peak sensitivity sits at the crossover, and for a 100 ns array that is 13.5 nHz.** The
frequency of peak SNR is essentially independent of N_psr and set entirely by timing:

| σ_toa | 100 ns | 30 ns | 10 ns | 3 ns | 1 ns |
|---|---|---|---|---|---|
| peak-SNR frequency | **13.5 nHz** | 23.5 | 39.1 | 68.0 | 113.1 |

**This is a mechanism for ORACLE's "not one conditioned source lies above 10 nHz".**
Conditioning finds sources where the GWB is loudest — the bottom of the band — which is
precisely the regime where N_c is 10⁵–10⁷ per bin, the loudest source is 9 % of the sea
amplitude, and nothing can be resolved. **The conditioned sources sit exactly where the
wheel provably cannot operate.** The wheel's window is *above* 13 nHz; the conditioned
population is *below* 10 nHz. They do not overlap.

## Question 3 — what array closes the loop

Closure is defined as the loudest source in some bin reaching SNR 8 **unaided** — that is
the gate, because until one source is detected there is nothing to subtract and nothing to
certify. The answer brackets on one modelling choice:

| confusion model | in band (≤32 nHz) | anywhere (≤200 nHz) |
|---|---|---|
| **uncorrelated** (sea averages as 1/N_psr — optimistic) | **N = 10 000, σ = 10 ns** (SNR 9.8) | N = 2 000, σ = 1 ns (SNR 8.7 at 113 nHz) |
| **correlated** (sea is an N-independent floor — pessimistic) | **never** (max SNR 0.11) | **never** (max SNR 0.32 at N = 50 000, 1 ns) |

**The bracket is the result.** A real GWB is Hellings–Downs correlated between pulsars, so
it does *not* average away with more pulsars — the truth is much closer to the pessimistic
column than the optimistic one. The reason the two diverge so violently:

| f | N_c/bin | h₁ / h_sea |
|---|---|---|
| 10 nHz | 6.9e4 | 0.090 |
| 32 nHz | 965 | 0.184 |
| 80 nHz | 33.5 | 0.322 |
| 200 nHz | 1.16 | 0.563 |

**Until N_c ~ 1 the loudest source is below the sea it sits in**, so if the sea cannot be
averaged down it cannot be dug out at any N_psr. N_c = 1 only at **208 nHz**, far above the
band.

## Question 1 — which mechanism, and the new blocker

At a closing array (5 000 psr, 3 ns, 1 pc priors):

| f | none | M1 | M2 | both | engine |
|---|---|---|---|---|---|
| 10–32 nHz | 0 turns | 0 turns | 0 turns | 0 turns | — |
| 50 nHz | 1 turn | +0.0009 | +0.0059 | **+0.0144** | M2 |
| 68 nHz | 2 turns | +0.0000 | +0.0137 | **+0.0227** | M2 |

**M2 is the engine, confirming VERNIER-3's ordering — but both are 10× smaller than
VERNIER-3 measured, and nothing at all happens at or below 32 nHz.**

### The certification/detection scissors — the new blocker

Detection is an **array** statistic: ρ_array = √(Σ_a ρ_a²), growing as √N_psr. Certification
is a **per-pulsar** statistic: VERNIER-1 established a pulsar needs its *own* ρ_P ≳ 2 for its
comb to collapse. At fixed source strength ρ_a ≈ ρ_array/√N_psr, so **growing the array to
reach detection drives per-pulsar SNR down**:

| N_psr | per-psr ρ at detection (ρ_array = 8) | ρ_array needed to certify |
|---|---|---|
| 116 | 0.743 | 21.5 |
| 1 000 | 0.253 | 63.2 |
| 5 000 | **0.113** | **141** |
| 20 000 | 0.057 | 283 |

A 5 000-pulsar array detects at ρ = 8 but needs ρ ≈ 141 — a source **18× louder** — before
its distances certify. That is why f_cert collapses to **0.029** here versus 0.2–0.9 in
VERNIER-3's 116-pulsar toy. **M1 and detection pull in opposite directions in N_psr**, so
"build a bigger array" closes the detection gate and simultaneously slams the certification
gate. This was invisible in VERNIER-1–3 because all of them fixed N_psr ≈ 116.

---

## What this cannot settle

1. **The confusion-correlation bracket is the whole answer and I did not resolve it.**
   Whether the sea averages as 1/N_psr or not is the difference between "10 000 pulsars at
   10 ns" and "never". A proper treatment needs the Hellings–Downs covariance and a
   simultaneous CW + GWB fit, which this toy does not have. **Highest-value follow-up by a
   wide margin.**
2. **N_c(f_yr) = 1000 and A_gwb = 2.4e-15 are inputs, not derived.** N_c enters as N_c^(−1/3)
   in R₁ and as N_c^(−1/3) in f_rem, so a 10³ error in N_c is only a 10× error in each — the
   qualitative picture is robust, the numbers are not.
3. **`dN/dh ∝ h^−4` is the Euclidean/no-evolution limit.** Real populations have redshift
   evolution and a mass-function cutoff that steepen the bright end, which would *help*
   resolvability. I have not modelled it.
4. **Closure at 113 nHz is an extrapolation, not a prediction.** The f^(−11/3) law taken to
   113 nHz gives ~1 source per bin in the entire sky; such binaries are within ~10³ yr of
   merger. Treat the in-band row as the real answer.
5. **Confusion is modelled as an isotropic Gaussian addition to the noise** with no sky
   structure, so nothing here captures separating a loud source from the sea by sky position
   — which is the actual mechanism a real search would use, and the main reason the
   optimistic column is not absurd.
6. **Still Asimov, still 1 pc priors, still no timing-model marginalisation.** Everything
   from VERNIER-1–3's caveat lists carries over unchanged.

## One-line verdict

On a realistic sky the wheel needs **~10⁴ pulsars at 10 ns even under an optimistic confusion
model and never closes under a pessimistic one**; M2 beats M1 but both fall to ~0.02 dex; the
operating window is pinned to the GWB/white crossover at **13.5 nHz for a 100 ns array**,
which is *above* where conditioning puts its sources — and the new structural blocker is that
**certification needs per-pulsar SNR while detection needs array SNR, so the array size that
buys one destroys the other.**
