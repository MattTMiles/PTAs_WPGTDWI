# VERNIER — does the fringe degeneracy have a multiplicity lever?

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, ~6 min of compute.
Clean room: nothing from `glacier_loop`, `summit_*`, `stagec_*` or any campaign module is
imported. Geometry, waveform, match and fringe bookkeeping are implemented from the
statement of the problem so that agreement with banked numbers is evidence, not tautology.

Code `hpc_harbor/vernier/vernier_toy.py` · data `reports/vernier_toy.npz` ·
plots `reports/VERNIER_t1a_fringe_budget.png`, `..._t1c_vernier.png`,
`..._t1d_commensurate.png`, `..._t2_brightness_gap.png`. Staged, not committed.

---

## The answer, first

**Yes — source multiplicity breaks the fringe degeneracy, and it is a genuine mechanism
rather than a decoration. It takes 3 sources at a 1 pc prior and 5–6 at literature priors.
It then dies on the brightness gap, and the wheel cannot manufacture the gap it needs.**

| question | toy answer |
|---|---|
| how many sources for a unique distance? | **K = 3** (σ_EM = 1 pc), **K = 5–6** (σ_EM = 10–20 % of L), K = 1 (σ_EM = 0.1 pc — already sub-fringe) |
| how bright must source 2 be? | ρ₂ ≳ √(ln K_trials / 2) ≈ **1.1–2.0** in the *pulsar term* |
| as a census dex gap (ρ₁ = 10)? | **≤ 0.7–0.95 dex**, tier-dependent |
| GLACIER census gap rank0−rank1 | **0.77–1.17 dex** — straddles the bar at the top and fails below it |
| can the wheel bootstrap? | **No.** Self-gearing ceiling is **0.15 dex**; the gap to close is 0.8–1.2 dex |

The mechanism is *geometric*, and that is the encouraging part: it barely cares about SNR.
Per-source pulsar-term ρ = 2 still reaches P(unique) ≥ 0.9 by K = 5, and ρ = 20 only buys
K = 2. What the mechanism needs is **more distinct combs**, not louder ones — the comb
spacing dL = c/(f₀(1−cos μ)) depends on frequency *and* sky position, so K sources supply K
incommensurate rulers.

The failure is therefore not "multiplicity doesn't work". It is that **multiplicity has to be
delivered by the population and cannot be manufactured by certification.** Certifying every
pulsar in the array converts it from Earth-term-only to Earth-plus-pulsar-term coherent —
a factor √2, i.e. **0.15 dex** of array SNR, and that is the hard ceiling on one turn of the
wheel. A 0.8–1.2 dex gap to the next source is 5–8× more than one turn can pay. The wheel
has a mechanism and no bootstrap.

---

## Anchor — and one discrepancy worth more than the sweeps

Instructed to check the fringe spacing against the banked ~0.4 pc comb. With
`log10_fgw = −7.7` (the campaign's loud-source fiducial, `CW_transition/stagec_p2a.py:110`):

```
c = 9.715612e-09 pc/s   f0 = 1.9953e-08 Hz
dL = c/f0 = 0.4869 pc          banked (J0437, median geometry) = 0.489 pc      → 0.42%
```

**Anchor passes.** All five banked D3 spacings close on ordinary mid-sky geometry:

| pulsar | banked dL | implied 1−cos μ | implied μ |
|---|---|---|---|
| J1744−1134 | 0.361 | 1.349 | 110.4° |
| J1909−3744 | 0.424 | 1.148 | 98.5° |
| J1713+0747 | 0.448 | 1.087 | 95.0° |
| J0030+0451 | 0.485 | 1.004 | 90.2° |
| J0437−4715 | 0.489 | 0.996 | 89.8° |

**Second anchor, unplanned and stronger.** A single source certifies a distance in this toy
only when the prior is *already* sub-fringe — never by chirp-breaking:

```
K_trials in [1,2)  n=728   P(unique) = 1.000     <- prior already narrower than one fringe
K_trials in [2,4)  n=789   P(unique) = 0.001
K_trials >= 4      n=1483  P(unique) = 0.000
```

and against prior width, single source: σ = 1 pc → 0.230, 10 pc → 0.020, 50 pc → 0.007,
**235 pc (the array's median feather prior) → 0.000**. That independently reproduces D3's
`class (i) = 0/116` *and* D3's own reading that the rare "identifiable" pulsars were
"K_a-small edge cases … not physical chirp-breaking". Two implementations, one conclusion.

### The discrepancy: the banked comb assumes a pulsar term that has not evolved

`dL = c/f₀` pins the comb to the **Earth-term** frequency. That is only right if the pulsar
term is at f₀ too. At the campaign fiducial it is not:

```
L = 1 kpc, mu = 90 deg, log10 Mc = 9.5, f0 = 20 nHz
  tau     = 3262 yr          light-travel delay
  t_merge = 1015 yr          coalescence time
  tau/t_c = 3.19             the pulsar term is emitted BEFORE the source entered the band
  f_p/f0  = 0.584            ->  dL = 0.834 pc,  1.71x the banked value
```

Across the band, τ/t_c ≥ 1 for (f₀ ≥ 10 nHz, Mc ≥ 10⁹·⁵) and reaches 21.8 at
(20 nHz, 10¹⁰). Two consequences:

1. **The naive quadratic-phase model is out of validity there.** ΔΦ = 2πf₀τ − πḟτ² + 2πḟtτ
   needs |ḟτ| ≪ f₀; at the fiducial |ḟτ| *exceeds* f₀, the carrier frequency passes through
   zero and the model **fabricates fringe structure**. My first pass had exactly this and it
   produced spurious single-source certifications. The toy now uses the exact Peters
   solution, which is finite for all τ > 0. This is a cousin of the banked sub-fringe-offset
   artefact: a bookkeeping choice that invents fringe evidence.
2. **Combs are ~1.8× coarser than banked** for evolving pulsar terms (Toy 1a median ratio
   1.82×). That *helps* fringe identification, mildly — it does not change any verdict here,
   and every headline below was run in **both** conventions with the same K₉₀.

I could not determine from outside the clean room whether the campaign's pulsar term
evolves; D3's "chirp breaks the degeneracy only at the ~1e-3 lnL level" is consistent with a
**non-evolving** pulsar term. **Flagging, not resolving: worth one afternoon on the real
code, because if the campaign's pulsar term is static at Mc = 10⁹·⁵ then its comb spacing is
~1.7× too fine and every K_a in the D3 table is ~1.7× too large.**

Waveform validated against brute-force time averaging: static branch |Δm| ≤ 5e−13, exact
branch ≤ 3.7e−4 (the neglected in-window quadratic term — 0.04 nats at ρ = 10, versus
thresholds of 1.4–8 nats).

---

## Toy 1a — fringe budget

20 pulsars, isotropic sky, L ~ U(0.3, 5) kpc, one source at f₀ = 10^−7.7.
K₃ = 6σ/dL is the campaign's ±3σ convention; K₁ = 2σ/dL is the ±1σ count, reported because
"a couple of fringes at 1 pc" is a ±1σ intuition.

| tier | σ | median K₃ | median K₁ | K₃<1 | K₁<1 |
|---|---|---|---|---|---|
| lit10 | 10 % of L | 1521 | 507 | 0/20 | 0/20 |
| lit20 | 20 % of L | 3041 | 1014 | 0/20 | 0/20 |
| vlbi | 1 pc | **11.4** | 3.8 | 1/20 | 2/20 |
| subpc | 0.1 pc | **1.14** | 0.38 | 9/20 | 20/20 |

**The hypothesis is directionally right and ~1 dex optimistic at the key tier.** Literature
priors hold 10³–10⁴ fringes: hopeless, as expected. But 1 pc holds a **median of 11** fringes
at ±3σ (3.8 at ±1σ) — not "a couple". Only 0.1 pc is genuinely sub-fringe.

So **D2's saturation at 1 pc is not the prior running out of fringes** — at 1 pc there are
still ~10 to remove. On this toy the prior runs out near **0.1 pc**, which is below VLBI for
all but the nearest pulsars. That relocates the claim rather than refuting it: the mechanism
that must remove those last ~10 fringes is multiplicity, which is Toy 1c.

## Toy 1b — the clock handle H2 never reads the timestamp

In the exact chirp the pulsar-term carrier drifts with distance as
dfₚ/dτ = −(3/8) fₚ/(t_c + τ). Across **one comb period** (dτ = 1/fₚ) the carrier therefore
moves by Δf = (3/8)/(t_c + τ), independent of frequency. The clock separates *adjacent*
fringes when Δf·T ≥ 1, i.e.

> **T_obs ≥ (8/3)(t_c + τ)**

A closed-form NO-GO: τ is millennia, T_obs is decades. Met in **0 of 60** cells across
f₀ ∈ [3, 20] nHz × Mc ∈ [10⁸·⁵, 10¹⁰] × T ∈ {15, 30, 40} yr. Surviving fringes
n_surv = (8/3)(t_c+τ)/T runs from 6×10² to 1×10⁶; the **tightest clock-only window anywhere
in the grid is 237 pc**, still 237× looser than a VLBI prior.

**This corrects the brief's framing.** The stated condition |ḟτ| ≥ 1/T_obs is the *absolute*
one — "does the clock distinguish L from L = 0" — and it is met widely (24/60 cells in the
first pass). It does **not** identify a fringe. The condition that does is T ≥ (8/3)(t_c+τ),
and nothing in the band comes within two orders of magnitude. H2 truncates the comb to
~10²–10⁶ candidates. It never resolves one. **The clock is not a second handle for
identification; it is a weak envelope on the first.**

## Toy 1c — the vernier (the headline)

Noiseless Asimov dlnL = Σₖ ρₖ²(mₖ(L) − 1). A mode survives within ln(K_trials) of the peak
(the campaign's trials-factor criterion). 120 trials/cell. K₉₀ = smallest K with
P(unique) ≥ 0.9.

| tier | P(unique) vs K = 1…8 (evolving pterm, ρ = 10) | K₉₀ evolving | K₉₀ static |
|---|---|---|---|
| lit10 | 0.00 0.00 0.01 0.42 **0.93** 1.00 1.00 1.00 | 5 | 5 |
| lit20 | 0.00 0.00 0.00 0.35 0.89 **0.99** 1.00 1.00 | 6 | 6 |
| vlbi | 0.20 0.82 **0.97** 1.00 1.00 1.00 1.00 1.00 | **3** | 3 |
| subpc | 1.00 1.00 1.00 1.00 1.00 1.00 1.00 1.00 | 1 | 1 |

Median surviving modes at the vlbi tier collapse 4 → 1 → 1; at lit20, 1039 → 156 → 16 → 2 → 1.
**Both pulsar-term conventions give identical K₉₀**, so this result does not rest on the
discrepancy above.

**The mechanism is geometric, not photometric** — the most useful thing the toy says:

| per-source pulsar-term ρ | 1.0 | 2.0 | 3.0 | 5.0 | 10 | 20 |
|---|---|---|---|---|---|---|
| K₉₀ (vlbi tier) | 8 | **5** | 4 | 3 | 3 | **2** |

A factor of 20 in per-source SNR moves K₉₀ by a factor of 4. Sources are worth far more than
brightness. The K = 1 column (P ≈ 0.20–0.24) is entirely the already-sub-fringe tail — see
the second anchor — not chirp-breaking.

## Toy 1d — commensurate vs incommensurate: harmonics are a floor, not a lever

K independent circular sources (own f₀, own μ, own t_c) against n harmonics of **one**
eccentric binary (fₙ = n·f_orb, *same* μ, *same* t_c), total SNR matched so this is a
structure comparison and not a brightness comparison. Peters–Mathews g(n,e) sets the
harmonic weights.

Median surviving modes, lit20 tier:

| n or K | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | factor |
|---|---|---|---|---|---|---|---|---|---|
| K independent sources | 848 | 157 | 34 | 6 | 2 | 1 | 1 | 1 | **848×** |
| n harmonics, e = 0.1 | 544 | 942 | 452 | 474 | 496 | 550 | 552 | 470 | **1.2×** |
| n harmonics, e = 0.3 | 538 | 366 | 430 | 372 | 376 | 444 | 394 | 490 | **1.1×** |

**Prediction confirmed, and more sharply than stated.** Harmonics do not dilate the ambiguity
by ~n — they **floor** it. The reason is exact: harmonic n has dLₙ = dL₁/n, so every
harmonic's fringe set *contains* the fundamental's, and the intersection of all of them is
pinned at dL₁ forever. Adding harmonics sharpens each fringe but removes no alias. Eccentric
harmonics are **not** a bounded second-order help to identification; to first order they are
**no help at all**. Source multiplicity is the whole lever.

*Eccentricity stated properly*: modest e ∈ {0.1, 0.3}, harmonics propagated into the waveform
via g(n,e) with the correct one-binary chirp. **No result here depends on e being large** —
and the conclusion strengthens as e → 0 (at e = 0.1, 89 % of power is in n = 2 alone).

**Bug found and fixed en route.** My first pass gave harmonic n the chirp of an independent
binary at n·f_orb, i.e. ḟ ∝ f^{11/3}, inflating its clock handle by n^{8/3} — 250× at n = 8.
That fake clock drove the e = 0.3 harmonic arm all the way to 1 mode and would have produced
the opposite conclusion. One binary's harmonics share one ḟ_orb: ḟₙ = n·ḟ_orb.

## Toy 2 — the brightness gap, as a dex number

Source 1 at pulsar-term ρ₁ = 10; source 2 at ρ₁·10^−dex. Combs enter dlnL weighted by ρ²,
so a dex gap is a 10^(2·dex) *weight* gap.

**Analytic ceiling, independent of the toy.** Source 2 can demote an alias by at most
2ρ₂² nats, so it must clear the trials factor: **2ρ₂² > ln K_trials**.

| K_trials | ln K | required ρ₂ | max gap from ρ₁ = 10 |
|---|---|---|---|
| 4 | 1.39 | 0.83 | 1.08 dex |
| 12 (J0437, banked) | 2.48 | 1.11 | **0.95 dex** |
| 250 | 5.52 | 1.66 | 0.78 dex |
| 3000 | 8.01 | 2.00 | 0.70 dex |

Monte Carlo, dex at which P(source 2 cuts the ambiguity) falls to 0.5: **0.91 (vlbi)**,
**0.76 (lit20)** — consistent with the analytic bar.

**Census read-off across GLACIER's 0.77–1.17 dex:**

| tier | P(cuts) 0.77 → 1.17 dex | P(unique) 0.77 → 1.17 dex |
|---|---|---|
| vlbi | 0.59 → **0.03** | 0.47 → **0.25** |
| lit20 | 0.45 → 0.20 | 0.00 → 0.00 |

**The census gap straddles the bar and fails through it.** At the optimistic end (0.77 dex,
1 pc priors) rank 1 cuts the ambiguity about half the time; at the pessimistic end
(1.17 dex) it has essentially stopped contributing. With literature priors two sources never
suffice at any gap — consistent with Toy 1c needing K = 5–6 there.

**The dex gap is the wrong variable, and this is the actionable point.** What the physics
constrains is rank 1's **absolute pulsar-term SNR** (ρ₂ ≳ 1.1–2.0). A gap is admissible iff
ρ₁/10^dex clears that bar, so **a louder rank 0 buys a wider admissible gap**. Read rank 1's
absolute strain off the census, not the gap to rank 0.

## Toy 3 — the wheel stalls on turn 1, and it stalls on gearing

20 pulsars, all given 1 pc VLBI-tier priors (already the best real case), 6 sources on a
census-like 0.6–1.2 dex ladder, 30 realisations, detection at array SNR > 8. Full
per-iteration history kept; certification recomputed every turn, never collapsed to
"certified at any point".

| turn | med. detected | med. certified | med. array SNR | med. gain |
|---|---|---|---|---|
| 1 | 1 | 4 / 20 | 19.7 | 1.000 |
| 2 | 1 | 4 / 20 | **3.00** | 1.095 |
| 3–6 | — | — | stalled | — |

**Turns before stall: median 1, range 1–1, in 30/30 realisations. Cause: 30/30 "next source
below detection threshold."** Not a certification failure — a *detection* failure. The wheel
certifies what it can from source 1 and then cannot find source 2.

**The gearing arithmetic, which is the real result.** Certifying every pulsar in the array
converts it from Earth-term-only to Earth-plus-pulsar-term coherent: a factor √2 = **0.15 dex**
of array SNR, and that is the ceiling on one turn. The gap to the next source is 0.8–1.2 dex.
0.15 cannot close 0.8. **The wheel is under-geared by 5–8×, and no amount of certification
success changes that** — it is arithmetic about what a pulsar term is worth, not a fit.

The 4/20 certified at turn 1 are the already-sub-fringe geometries (second anchor), and they
are optimistic by construction: every pulsar here was handed a 1 pc prior, where the real
array's median feather prior is ~235 pc, at which single-source P(unique) = 0.000.

---

## What this toy cannot settle

1. **Whether the campaign's pulsar term evolves.** The τ/t_c ≥ 1 finding says the banked
   `dL = c/f₀` convention is only valid for a non-evolving pulsar term, and that at
   Mc = 10⁹·⁵ the comb may be ~1.7× coarser than banked. Resolving this needs the real
   waveform code, which the clean room excluded. **Highest-value follow-up.**
2. **Whether K = 3 survives real machinery.** No GWB, no timing-model marginalisation, no
   per-pulsar red noise, no source-parameter uncertainty, no sky/frequency covariance
   between sources. Every one of those *adds* confusion, so K₉₀ = 3 is a **lower bound** and
   the honest expectation is that the real number is larger.
3. **Whether the population supplies K = 3 at ρ₂ ≳ 2.** This is now the whole question, and
   it is a census question, not a fringe question. The toy says what is needed; only
   GLACIER/CHORUS can say whether the sky delivers it.
4. **Flicker.** Toy 3 recorded **zero** flicker events, but it is Asimov (noiseless) — and
   flicker in this programme is a noise-and-refit phenomenon. **Zero here is not evidence
   against flicker**; the toy simply cannot see it.
5. **Detection-threshold realism.** "Array SNR > 8" and the √(1 + N_cert/N_psr) coherent-gain
   model are caricatures. The 0.15 dex ceiling is robust (it is just √2), but the turn-1 stall
   depends on where the threshold sits relative to the ladder.
6. **Calibration.** Everything is Asimov mode-counting against a trials factor. There is no
   false-alarm rate here, so nothing in this report is a calibrated threshold and no count
   here is a certification in the programme's sense.

## One-line verdict

The wheel has a mechanism — **3 incommensurate combs, and SNR barely matters** — but no
bootstrap: certification pays 0.15 dex and the next source costs 0.8–1.2 dex. Stop asking
whether fringe identification is possible (it is, at K = 3) and start asking whether the sky
hands us three sources whose pulsar terms each reach ρ ≳ 2.
