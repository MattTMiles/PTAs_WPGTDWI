# VERNIER-25 — the wheel-friendly sky, specified; and the first-detection era

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, **~40 min** (not minutes —
the certification pass dominates and there are ~300 chains).
Clean room: imports only `vernier24` (and through it `vernier_toy` / `vernier6` / `vernier7` /
`vernier16`).

Code `hpc_harbor/vernier/vernier25.py` · data `reports/vernier25.npz`.
Staged, not committed.

Reuses VERNIER-24's corrected-sky construction, its realized M_max, and its incoherent-comb
certification floor.

---

## Lead — nothing ignites, on any sky and at any K, and the reason is a single number

**Ignition fraction is 0.00 in every cell of every arm: 0/20 at each of six slopes, 0/12 at
each of five sea suppressions, 0/10 at each of seven seed counts. γ\* does not exist inside
0.6–2.5 and K\* is not reached at K ≤ 16.**

**But the arms are not equally far away, and the thing they are short of is the same thing.**
Certification is gated by fringe contrast against a trials bar of ln K ≈ 4.6–5.3 nat, and
contrast is a function of the number of **accumulated, genuinely-recovered combs** — not of
any source's SNR, not of the slope, not of the sea:

| cell | combs accumulated | end-state contrast | bar | S/bar | best single chain |
|---|---|---|---|---|---|
| Arm G, γ = 1.2 | 2.5 | 0.031 | 4.58 | 0.007 | 0.15 |
| **Arm G, γ = 0.6** | 4.0 | 0.226 | 4.79 | 0.047 | 1.07 |
| Arm K, K = 5 | 5.0 | 0.278 | 4.75 | 0.059 | 1.60 |
| Arm K, K = 8 | 6.0 | 1.312 | 5.18 | 0.253 | 3.99 |
| Arm K, K = 12 | 6.5 | 1.201 | 5.34 | 0.225 | 2.48 |
| **Arm K, K = 16** | **8.0** | **2.724** | 5.32 | **0.512** | **5.93 — clears the bar** |

**Ordered by combs, the two arms fall on one curve.** Arm G's flattest sky (γ = 0.6) sits at
4 combs and Arm K's largest seeding at 8; the contrast each reaches is what that comb count
buys, and it does not matter which dial produced it. **At 8 combs the median chain is half way
to the bar and the best of six chains crosses it.**

**Extrapolating that curve, the bar is reached at ~10–11 accumulated combs.**

## The cross-arm answer: yes, one number — but it is not K, it is K_recovered

**γ\* and K\* are the same statement in different coordinates, and the coordinate both of them
are really tracking is the number of sources the chain GENUINELY RECOVERS at distinct
frequencies.** Neither raw availability nor raw seeding is that number:

| | sources available (K) | **genuinely recovered** | efficiency |
|---|---|---|---|
| Arm G, γ = 0.6 | 9.0 | **3.6** | 0.40 |
| Arm G, γ = 0.87 | 10.0 | **2.9** | 0.29 |
| Arm G, γ = 1.56 | 4.0 | **1.9** | 0.48 |
| Arm K, K = 8 | 8 seeded | **5.3** | 0.66 |
| Arm K, K = 16 | 16 seeded | **7.9** | 0.49 |

**The gap between the two columns is blends.** At γ = 0.87 ten sources clear the detection
floor and 9.2 of the chain's 12 turns land on blends rather than on any of them. **A flat
slope supplies detectable sources; it does not supply recovered ones**, because the F-statistic
peak in a 300-source bin sits between sources rather than on one.

**So the wheel's condition does collapse to one number — but seeding buys it about 1.6× more
efficiently than flattening the slope does, and no sky in the sweep supplies more than ~4.**

## Arm G — the gamma sweep

| γ | **ignite** | K | λ_fore | turns | genuine | blend | deepen | CERT | wrong | depth | rank | stall |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.60 | **0.00** | 9.0 | 0.874 | 11.0 | 3.6 | 7.5 | 0.00 | 0 | 0 | +0.121 | **7** | comb ambiguity |
| 0.87 | **0.00** | 10.0 | 0.637 | 12.0 | 2.9 | 9.2 | 0.00 | 0 | 0 | +0.073 | 6 | comb ambiguity |
| 1.20 | **0.00** | 8.5 | 0.541 | 12.5 | 2.6 | 9.8 | 0.05 | 0 | 0 | +0.056 | 5 | comb ambiguity |
| **1.56** | **0.00** | 4.0 | 0.321 | 10.0 | 1.9 | 8.6 | 0.00 | 0 | 0 | +0.005 | 5 | comb ambiguity |
| 2.00 | 0.00 | 0.5 | 0.125 | 7.0 | 0.6 | 6.2 | 0.00 | 0 | 0 | +0.000 | 4 | absolute floor |
| 2.50 | 0.00 | 0.0 | 0.044 | 4.0 | 0.1 | 3.5 | 0.00 | 0 | 0 | +0.000 | 4 | absolute floor |

**The detection channel has no positive feedback at all.** The `deepen` column is the fraction
of turns in which subtracting a source exposed something that had been *below the floor*, and
it is **0.00 everywhere except a single realisation at γ = 1.2**. On every sky tested, every
source the chain finds was already detectable before the chain started. **The wheel does not
lower its own floor.**

**Flattening the slope does buy depth, just not through the wheel:** γ = 0.6 reaches +0.121 dex
and rank 7 against +0.005 dex and rank 5 at γ = 1.56 — that is subtraction of a genuinely
brighter bright end, not a turning wheel.

### The boundary in four currencies

| γ | ignite | K | λ_fore | h₍₁₎ | A_eq/A₁₅ | N_eff | N_eff ≤ V24's M_max (5–10)? |
|---|---|---|---|---|---|---|---|
| 0.60 | 0.00 | 9.0 | 0.874 | 9.65e-16 | 1.66 | 1.3 | yes |
| 0.87 | 0.00 | 10.0 | 0.637 | 8.97e-16 | 1.54 | 1.8 | yes |
| 1.20 | 0.00 | 8.5 | 0.541 | 7.13e-16 | 1.23 | 3.3 | yes |
| **1.56** | 0.00 | 4.0 | 0.321 | 6.27e-16 | 1.08 | 5.7 | **yes** |
| 2.00 | 0.00 | 0.5 | 0.125 | 3.72e-16 | 0.64 | 36.5 | NO |
| 2.50 | 0.00 | 0.0 | 0.044 | 2.27e-16 | 0.39 | 183.9 | NO |

**Where γ = 1.56 sits — and which ingredient binds first.** It is short in exactly one:

- **N_eff = 5.7 is already INSIDE VERNIER-24's realized M_max of 5–10.** Separability is not
  the binding constraint at the true slope.
- **h₍₁₎ is 0.16 dex below the γ = 0.87 value** — a small deficit, and A_eq/A₁₅ = 1.08 says
  the true sky's loudest source is already at the effective loudness of an NG15-amplitude
  background.
- **K = 4 available, K_recovered = 1.9.** *That* is the deficit, and it is a factor of ~5
  against the ~10 combs the contrast curve needs.

**The binding ingredient is recovered multiplicity, and nothing else is close.**

## Arm G-FP — the final-parsec alternative, and why it has no lever

| sea kept | A_eq/A₁₅ | ignite | K | turns | genuine | CERT | depth | rank | stall |
|---|---|---|---|---|---|---|---|---|---|
| 1.00 | 1.00 | 0.00 | 4.5 | 11.5 | 1.8 | 0 | +0.010 | 5 | comb ambiguity |
| 0.50 | 0.98 | 0.00 | 4.5 | 11.5 | 1.7 | 0 | +0.006 | 5 | comb ambiguity |
| 0.20 | 0.97 | 0.00 | 6.0 | 12.0 | 2.2 | 0 | +0.023 | 5 | comb ambiguity |
| 0.05 | 0.96 | 0.00 | 8.0 | 14.0 | 2.2 | 0 | +0.016 | 5 | comb ambiguity |
| **0.00** | 0.96 | **0.00** | 6.0 | 14.0 | 1.9 | 0 | +0.006 | 5 | comb ambiguity |

**Deleting the faint sea entirely changes nothing, and the reason is structural: at γ = 1.56
the sea below the 300 explicitly-modelled sources holds only 2.3 % of the bin power.** A
Pareto with γ < 2 puts essentially everything in its bright end, so there is no faint sea left
to stall. (At γ = 2.0 the sea holds 41 % — but that sky has K = 0.5 and cannot detect anything
in the first place.)

**So the answer to "which mechanism does the wheel care about" is neither, and that is
informative rather than null:** the confusion that blocks the wheel at the true slope is not
an unresolved background at all. **It is the other resolvable sources in the same bin** — the
ranks 2–300 that a search turns into blends. A final-parsec cut that stalls small binaries
removes power the wheel was never blocked by. **The astrophysically quotable version is
therefore not "the wheel turns if the final parsec stalls small binaries" — it does not.**

## Arm K — the first-detection era

| K | A_eq/A₁₅ | ignite | turns | seeds hit | **beyond seed** | total known | CERT | wrong | null | depth | dex/turn | new det/turn | rank |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1.00 | 0.00 | 11.0 | 0.0 | 1.3 | 1.3 | 0 | 0 | 0 | +0.001 | +0.0001 | 0.12 | 5 |
| 2 | 1.65 | 0.00 | 12.0 | 0.7 | 1.5 | 2.2 | 0 | 0 | 0 | +0.015 | +0.0012 | 0.18 | 5 |
| 3 | 2.22 | 0.00 | 13.0 | 1.1 | 1.0 | 2.1 | 0 | 0 | 0 | +0.018 | +0.0014 | 0.16 | 5 |
| 5 | 1.99 | 0.00 | 12.5 | 2.6 | 0.9 | 3.5 | 0 | 0 | 0 | +0.034 | +0.0028 | 0.28 | 5 |
| 8 | 1.14 | 0.00 | 16.0 | 4.0 | 1.3 | 5.3 | 0 | 0 | 0 | +0.062 | +0.0038 | 0.33 | **6** |
| 12 | 1.60 | 0.00 | 16.0 | 5.7 | 1.4 | 7.1 | 0 | 0 | 0 | +0.053 | +0.0033 | 0.44 | 5 |
| 16 | 1.39 | 0.00 | 16.0 | 6.3 | 1.6 | 7.9 | 0.1 | **0.0** | 0 | +0.074 | +0.0046 | 0.49 | **6** |

**1. Does the vernier fire at K ≥ 3 against the real sea? No — but it is no longer flat.** Zero
certifications up to K = 12; **0.1 per run at K = 16, and every one of them on the correct
fringe (wrong = 0.0, scrambled null = 0.0)**. The untested middle between VERNIER-6's clean
corner and VERNIER-24's cold start turns out to be a slope, not a step: contrast climbs
0.00 → 0.04 → 0.18 → 0.56 → 0.71 (per-turn medians) and 0.28 → 1.31 → 2.72 (end-state) across
K = 5 → 16, against a bar that only creeps from 4.75 to 5.32.

**2. Seeding is NOT invariant on the honest sky — VERNIER-3's result does not survive the
census correction.** Total known = seed + turns:

| K | 0 | 2 | 3 | 5 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|---|
| **total known** | 1.3 | 2.2 | 2.1 | 3.5 | 5.3 | 7.1 | **7.9** |

**Invariance would be a flat row at 1.3; the measured row spans 1.3 → 7.9.** Seeding buys
almost exactly what it puts in and no more — *"beyond seed"* is pinned at **0.9–1.6 across the
whole range**, statistically flat. **So seeding buys total knowledge but buys no leverage: the
chain never converts a seed into an extra detection.** VERNIER-3's "seeding buys nothing" was
measured on the artefact census and is now wrong in its literal form and right in its spirit.

**3. Marginal value per turn rises, slowly, and never stalls in the interior.** +0.0001 dex and
0.12 new detections per turn at K = 0, rising monotonically to +0.0046 dex and 0.49 at K = 16.
**Nothing plateaus, so the chain is not hitting an interior limit** — it runs out of turns and
combs, not out of signal.

**4. Where it stalls, named: comb ambiguity, in every single Arm K cell.** Not the absolute
floor (which only binds at γ ≥ 2.0, where nothing is detectable), and not confusion (the
`conf_frac` diagnostic never selects it). **The chain detects; it cannot pin a fringe.**

**5. K\* is not reached at K ≤ 16, but it is now bracketed rather than open.** The best single
chain at K = 16 reaches 5.93 against a bar of 5.32 — **it crosses**. The median does not.
Reading the contrast-vs-combs curve, **the bar wants ~10–11 accumulated combs, which at the
measured recovery efficiency of ~0.5 means K\* ≈ 20–25 simultaneous detections.**

**Against the forecasts that motivated the arm — ~4 resolvable MeerKAT-class in 10 yr, 30–40
SKA in 30 yr — K\* ≈ 20–25 is out of reach this decade and inside the SKA number.**

## Controls

- **M = 1 control** (single bin, seeded, sea removed): detected 4/6, **0 certifications**,
  ρ_P max 3.00, fringe contrast **0.000 nat**. One comb cannot break its own degeneracy at any
  SNR, exactly as VERNIER-24 found.
- **Pure-noise null through the whole chain**: fired 1/6 against an expectation of 0.40 for
  10 bins at α = 0.05; **0 certifications, 0 ignitions**.
- **Per-turn source-level scrambled null**: `wrong` and `null` are 0.0 in every cell that
  certified at all, and a scrambled-null certification vetoes the deepening flag by
  construction, so **no ignition claim in this note could have been scrambled-dirty** — there
  are none to be dirty.
- **Chance-match floors**: 0.006 at M = 1, 0.055 at M = 10, 0.815 at M = 300. GENUINE is scored
  against the top-10 only, so the 0.055 floor is the operative one.

---

## What this does not settle

1. **The comb cap is mine and it binds at exactly the point that matters.** Certification
   carries at most **COMB_CAP = 8** combs (the cost is linear in the count). At K = 12 and
   K = 16 the chain has more recovered sources than that, so **the K = 12/16 contrasts are
   lower bounds, and the extrapolated K\* ≈ 20–25 is an upper bound on the requirement.** This
   is the single most leveraged assumption in the note and it should be relaxed before K\* is
   quoted anywhere.
2. **`max_turn = 16` truncates the largest cells.** At K = 16 the chain uses all 16 turns and
   still has seeds unrecovered (6.3 of 16 hit), so that row is turn-limited as well as
   comb-limited.
3. **λ_fore is a median of a very broad distribution.** At γ = 1.56 the IQR is 0.19–0.55 and a
   200-draw median has sd ≈ 0.02; the deliverable-0 value of 0.361 is a high draw against a
   long-run 0.317. Every γ-cell observable inherits that spread.
4. **20 realisations per γ, 12 per sea cell, 10 per K, 6 for the end-state supplement.** An
   ignition fraction of 0.00 out of 20 bounds the rate at ≲ 14 % (95 %), not at zero.
5. **The end-state contrast is max-over-turns of a per-pulsar max.** It is the most optimistic
   summary of the certification channel that is still honest, and it is the one the K\*
   extrapolation rests on; the per-turn medians in the Arm K table are 5–20× smaller.
6. **Detection is Earth-term only** (VERNIER-16/22's standing limitation), so the blend rate
   that eats the K → K_recovered efficiency is measured with the weaker search.
7. **A_equivalent is a loudness label, not a claim about a typical source** — it is defined
   from h₍₁₎ at fixed total power, and it is the only sense in which these cells are
   "louder" or "fainter" than NG15.
8. **The final-parsec arm suppresses only the sea below the top 300.** A stall cut that reached
   *into* the modelled bright end is a different experiment and is the one the null result
   above actually points at.
9. **Nothing here re-derives VERNIER-5's γ ≤ 0.87 threshold** — it is quoted as the target and
   the sweep straddles it; the wheel simply does not care where it falls.

## One-line verdict

**Nothing ignites — 0.00 ignition fraction in all eighteen cells across six bright-end slopes
from γ = 0.6 to 2.5, five degrees of final-parsec sea suppression, and seven seed counts from
K = 0 to 16 — so γ\* does not exist inside the swept range and K\* is not reached at K ≤ 16;
the detection channel has no positive feedback whatsoever (subtracting a source exposed
something previously below the floor in one realisation out of 240), so the whole question
lives in the certification channel, which is gated by fringe contrast against a ~5 nat trials
bar and is a function of one quantity only — the number of genuinely recovered combs at
distinct frequencies — on which Arm G and Arm K fall on a single curve, confirming that γ and
K are two clocks for one number; that number is not K but K_recovered, because a flat slope
supplies detectable sources and not recovered ones (10 available and 2.9 recovered at
γ = 0.87, nine turns in twelve landing on blends), and seeding supplies it ~1.6× more
efficiently while buying no leverage at all (total known rises 1.3 → 7.9 with the seed while
"beyond seed" stays pinned at 0.9–1.6, so VERNIER-3's invariance is broken in letter and
upheld in spirit); the final-parsec alternative has no lever because at γ < 2 the faint sea
holds 2.3 % of the bin and the confusion that actually blocks the wheel is the other
resolvable sources in the same bin, which is why the astrophysically quotable version is that
the wheel does NOT turn in a universe where the final parsec stalls small binaries; the
best-case cell (K = 16, eight accumulated combs) reaches an end-state contrast of 2.72 against
a bar of 5.32 with its best single chain crossing at 5.93 and every certification it makes
landing on the correct fringe, which brackets K\* at roughly 20-25 simultaneous detections —
inside the SKA 30-year forecast and out of reach of the MeerKAT-era one — subject to the
caveat that my own 8-comb cost cap binds at exactly the point where the mechanism starts to
work, so that bracket is the first thing to re-measure.**
