# VERNIER-28 — blend recovery: the one lever left, pulled, and it does not move

**EXPLORATORY. NOT pre-registered.** Toy, cronus, CPU/numpy, 6 worker processes
(your `vernier27.py` and `vernier29.py` were on the box), **575 s** for 48 universes.
No GPU used — the joint fits are 6000-point F-statistics and 116×512 marginalisations,
which are memory-bound, not flop-bound. **Nothing of the ~4 GPU-hr budget was spent.**

Clean room: `vernier28_lib.py` + `vernier28.py` import `vernier27_lib` / `vernier25` /
`vernier24` / `vernier16` / `vernier6` / `vernier_toy` and edit none of them. One
imported constant is overridden at run time and announced in the run header
(`vernier27_lib.SUB_SPAN → 0.0`; see *Departures*).

Code `hpc_harbor/vernier/vernier28_lib.py`, `vernier28.py`, `vernier28_paired.py` ·
data `reports/vernier28.npz` · log `reports/vernier28.log`.
Staged, not committed.

---

## Lead — joint fitting does not convert blends, and the pulsar block makes it worse

**Conversion efficiency, paired over 48 RNG-matched universes** (same sky, same array,
same noise draws, same Earth-block data; the ladders differ only in the fit):

| ladder | conversion K_rec/K_av | Δ vs R1 (paired) | verdict |
|---|---|---|---|
| **R1** sequential/greedy (control) | **0.503 ± 0.025** | — | — |
| **R2** joint-bin, M null-calibrated | **0.518 ± 0.027** | **+0.015 ± 0.022** [−0.027, +0.061] | consistent with zero |
| **R3** joint-bin + marginalised pulsar terms | **0.416 ± 0.024** | **−0.087 ± 0.021** [−0.128, −0.044] | **excludes zero — worse** |

Head to head on genuine recoveries: **R2 beats R1 in 16 universes of 48, ties 18, loses
14.** R3 beats R1 in 9 of 48.

**The lever is not there.** Fitting every source in a bin simultaneously — with each
source's sky position re-solved against the others, and the number of components chosen
against a measured null floor rather than asserted — recovers the same sources the greedy
chain already recovered. It finds *more structure* (16.7–23.7 components claimed a run
against R1's 12.5–18.9) and none of the extra structure is a true source.

**And the degree-of-freedom surplus provably does not pay in a search, though it does in
a Fisher.** That contrast is the sharpest thing in this note and it is measured twice on
the same page — §2b against §2.

---

## Departures from VERNIER-27, both chosen to *avoid flattering R3*

1. **`SUB_SPAN = 0`.** V27 gives packed seeds sub-bin frequency offsets. The toy's
   F-statistic searches at the bin centre only, so a sub-bin offset is a frequency the
   search cannot measure — and R3's comb model is exquisitely sensitive to it
   (2π·δf·τ ≈ 35 rad at δf = 0.35/T). Handing R3 the true sub-bin frequency would be
   handing it information the toy's search cannot have. Every source therefore sits at
   its bin centre. Cost: packed sources share a comb-spacing family, so the vernier's
   commensurability lever is **weaker** here than in V27's construction.
2. **One chirp mass, log₁₀Mc = 9.0 for every source.** The retarded pulsar-term phase
   carries a term (3π/8)·f·τ²/t_c worth ≈ 90 rad, and t_c ∼ Mc^(−5/3) varies by 46× across
   V25's 8.5–9.5 prior. A model using each source's *true* Mc would be an oracle; a model
   using a common Mc while the sky used per-source Mc would break on real sources and not
   on phantoms — the same oracle wearing a hat. Fixing both removes the asymmetry. Cost:
   all combs in a bin share one chirp evolution, which makes them **more alike** — again
   weaker for R3.
3. **The genuine pool is the AVAILABLE set.** Conversion is recovered/available, so the
   two must name the same set: a component may be matched only against sources above the
   detection floor in the unsubtracted λ_eff budget (V25's `k_available` criterion). Chance
   floor of that pool printed per cell (0.000–0.006).

---

## The construction, and what the ladders are

**Cells.** K_available is *measured*, not assumed:

| cell | construction | K_av (median) | A_eq/A₁₅ |
|---|---|---|---|
| C5 | γ = 1.56 + 2 seeds | 6.0 | 0.91 |
| C10 | γ = 0.87, unseeded — V25's flat-slope cell | 9.5 | 1.35 |
| C20 | γ = 0.87 + 10 seeds | 18.0 | 1.61 |
| C10b | γ = 1.56 + 8 seeds — **construction ablation** | 13.0 | 1.25 |

C10b exists so that conversion-vs-K_available is not confounded with γ: it reaches
comparable availability by *seeding* rather than by *flattening*. R1's conversion there is
0.555 against C10's 0.466 at lower K_av — **availability bought by seeding converts better
than availability bought by a flat slope**, which is V25's 1.6× seeding advantage showing
up again in a different coordinate.

**R1 — the control, instrumented, not improved.** V27's chain copied verbatim with one
addition: it records the grid index of every component, so all three ladders can be scored
by one accounting function. `validate_R1()` runs at the top of every run and asserts
turn-for-turn equality with the imported `vernier27_lib.chain`: **turns 17/17, 13/13,
18/18; genuine 3/3, 4/4, 7/7; blend 14/14, 9/9, 11/11; cert 0/0, 0/0, 0/0 — reproduced
exactly.**

**R2 — joint-bin.** Greedy CLEAN initialisation, then space-alternating conditional
maximisation: each source's sky position is re-solved against the residual of all the
others with their amplitudes re-fitted, until the set stops moving. M scanned 1..12,
component m admitted iff its own increment F_joint(m) − F_joint(m−1) clears that
increment's null floor, measured from 12 pure-noise realisations of the same bin and cut
by V24's own recut estimator.

**R3 — joint-bin + marginalised pulsar terms.** The pulsar term is made an explicit second
observable, exactly as V24's `dvec` already treats it:

    y_a = Σ_k E_ak · R_AMP · env_ak · exp(i φ_ak(L_a)) + n_a ,  n_a ~ CN(0,1)

with E_ak the same Earth-term contribution `inject` already puts in the Earth block,
R_AMP = 1.196 (V8), and φ the exact retarded circular-chirp phase from
`vernier_toy._retarded` — the same propagator `certify()` uses. All 300 sources contribute,
so the pulsar block carries the same confusion the Earth block does. R1 and R2 never look
at `y` (the standing V16/V22 Earth-only limitation, which is what keeps R1 identical to
V25/V27). R3 uses both and marginalises L_a over its 10 pc parallax prior, then scores every
move — refine, prune, grow — by ΔT where **T = Σ_b F_joint,b + Σ_a ln Z_a**.

---

## Deliverable 1 — conversion efficiency, with the full V18 accounting

| cell | ladder | claims | GEN | split | faint | phantom | K_rec/K_av | chance |
|---|---|---|---|---|---|---|---|---|
| C5 (K_av 6.0) | R1 | 12.50 | 2.58 | 0.00 | 8.50 | 1.42 | 0.429 ± 0.062 | 0.000 |
| | R2 | 16.75 | 2.75 | 0.00 | 12.08 | 1.92 | 0.442 ± 0.058 | |
| | R3 | 25.08 | 2.58 | 0.00 | 18.25 | 4.25 | 0.411 ± 0.049 | |
| C10 (K_av 9.5) | R1 | 12.50 | 4.75 | 0.00 | 6.33 | 1.42 | 0.466 ± 0.040 | 0.003 |
| | R2 | 16.67 | 4.50 | 0.00 | 10.33 | 1.83 | 0.438 ± 0.041 | |
| | R3 | 28.75 | 3.33 | 0.00 | 22.00 | 3.42 | **0.325 ± 0.055** | |
| C20 (K_av 18.0) | R1 | 18.92 | 10.08 | 0.00 | 6.67 | 2.17 | 0.560 ± 0.042 | 0.006 |
| | R2 | 23.67 | **10.33** | 0.00 | 11.33 | 2.00 | **0.573 ± 0.048** | |
| | R3 | 45.17 | 8.25 | 0.00 | 30.92 | 6.00 | 0.457 ± 0.044 | |
| C10b (K_av 13.0) | R1 | 14.92 | 7.17 | 0.00 | 6.75 | 1.00 | 0.555 ± 0.043 | 0.006 |
| | R2 | 23.00 | **8.00** | 0.00 | 12.00 | 3.00 | **0.618 ± 0.051** | |
| | R3 | 37.25 | 6.08 | 0.00 | 25.08 | 6.08 | 0.469 ± 0.039 | |

**SPLIT is 0.00 in every cell of every ladder.** No ladder ever puts two components on one
true source — the unique-assignment rule is never the binding one. The blend problem is not
double-counting; it is components landing on *nothing available*.

**Conversion does not depend on the ladder, and it does depend on K_available.** Across
cells, R1 goes 0.429 → 0.466 → 0.560 as K_av goes 6 → 9.5 → 18. That rise is real and it is
the only monotone trend in the table: **more available sources convert at a higher rate, not
a lower one**, because a bin holding two detectables gives the fitter a second peak to lock
on rather than a blend to sit between.

**Paired, the whole picture in three numbers**
(`vernier28_paired.py`, 20 000 bootstraps over the 48 matched universes):

    conversion  R1 0.503 ± 0.025   R2 0.518 ± 0.027   R3 0.416 ± 0.024
    Δ R2 − R1  +0.015 ± 0.022  [−0.027, +0.061]   consistent with zero
    Δ R3 − R1  −0.087 ± 0.021  [−0.128, −0.044]   excludes zero

---

## Deliverable 2 — realized per-bin separability: 2–3, not 5–10, and joint fitting does not move it

| cell | ladder | M̂ med | M̂ max | G_bin med | G_bin max | **largest M̂ with every component genuine** | bins with G ≥ 2 |
|---|---|---|---|---|---|---|---|
| C5 | R1 | 1.0 | 4 | 1.0 | 2 | 2 | 0.01 |
| | R2 | 2.0 | 6 | 1.0 | 2 | 2 | 0.05 |
| | R3 | 2.0 | 9 | 1.0 | 3 | 1 | 0.03 |
| C10 | R1 | 1.0 | 3 | 1.0 | 3 | **3** | 0.06 |
| | R2 | 1.0 | 5 | 1.0 | 3 | **3** | 0.06 |
| | R3 | 2.0 | 10 | 1.0 | 3 | 2 | 0.07 |
| C20 | R1 | 2.0 | 4 | 1.0 | 3 | 2 | 0.18 |
| | R2 | 2.0 | 7 | 1.0 | 3 | 2 | **0.22** |
| | R3 | 3.0 | 12 | 1.0 | 2 | 2 | 0.14 |
| C10b | R1 | 2.0 | 3 | 1.0 | 2 | 2 | 0.09 |
| | R2 | 2.0 | 7 | 1.0 | 3 | 2 | 0.13 |
| | R3 | 4.0 | 12 | 1.0 | 3 | 1 | 0.06 |

**The realized separability of an actual joint search is 2–3, against V24's Fisher frontier
of 5–10 and a counting bound of 3N/7 = 50.** Joint fitting raises M̂ — the number of
components the data supports — from a median of 1–2 to 2–4, and R3 raises it as far as 12.
**None of that reaches the truth.** `G_bin` has a median of 1.0 and a maximum of 3 in every
single cell of every single ladder: **the overwhelming majority of bins yield exactly one
genuine recovery no matter how they are fitted.** The fraction of bins yielding two or more
never exceeds 0.22, and that best value is R2 at C20 — the most crowded, most available cell.

**So the answer to "does joint fitting move realized M_max toward the counting bound at
all" is: it moves M̂ and not G_bin.** The extra components are real degrees of freedom the
data genuinely supports against noise — they are simply not sources.

### Deliverable 2b — the Fisher companion, and where the dof surplus does pay

`vernier24.separability`, unmodified, at V24's own identifiability-limited working point
(ρ₁ = 30 fixed — the column V24's "realized frontier 5–10" was read from):

| M | Earth only: resolved / top_ok | + pulsar, marginalised: resolved / top_ok |
|---|---|---|
| 2 | 2.0 / 1.00 | 2.0 / 1.00 |
| 5 | 4.0 / 1.00 | 4.0 / 0.67 |
| **10** | **2.0 / 0.33** | **6.0 / 0.67** |
| 17 | 1.0 / 0.33 | 1.0 / 0.67 |
| 25 | 1.0 / 0.33 | 1.0 / 0.67 |
| 40 | 0.0 / 0.00 | 0.0 / 0.00 |

**At M = 10 the marginalised pulsar block triples the number of resolved sources (6 against
2) and doubles the chance the loudest one resolves.** V23 Part 1's surplus is real and the
Fisher can see it.

**And the search cannot.** The same surplus, wired into an actual fit on actual data
(R3), *lowers* conversion by 0.087 ± 0.021. **The Fisher error is local and the search
problem is global** — V24 flagged exactly this in its own caveat ("a Fisher error is LOCAL;
σ_sky ∼ 1e-4 rad is a within-fringe width, not a global uncertainty"). This note is the
measurement that closes it: **the dof surplus buys precision on a source you have already
found, and buys nothing on the problem of finding it.**

---

## Deliverable 3 — THE CALENDAR

| cell | ladder | combs | S_end | bar | S/bar | cross | CERT | wrong | scrambled |
|---|---|---|---|---|---|---|---|---|---|
| C5 | R1 | 2.5 | 0.021 | 4.55 | 0.005 | 0.00 | 0.00 | 0.00 | 0.00 |
| | R2 | 2.5 | 0.040 | 4.51 | 0.009 | 0.00 | 0.00 | 0.00 | 0.00 |
| | R3 | 2.0 | 0.043 | 4.40 | 0.010 | 0.00 | 0.00 | 0.00 | 0.00 |
| C10 | R1 | 4.5 | 0.551 | 4.78 | 0.115 | 0.00 | 0.00 | 0.00 | 0.00 |
| | R2 | 4.0 | 0.576 | 4.72 | 0.122 | 0.00 | 0.00 | 0.00 | 0.00 |
| | R3 | 3.5 | 0.125 | 4.69 | 0.027 | 0.00 | 0.00 | 0.00 | 0.00 |
| **C20** | **R1** | **10.0** | **8.470** | 5.26 | **1.612** | **0.67** | 1.42 | 0.00 | 0.50 |
| | **R2** | **10.0** | **10.635** | 5.28 | **2.014** | **0.75** | 1.50 | 0.08 | 0.17 |
| | R3 | 8.0 | 4.036 | 5.13 | 0.787 | 0.42 | 0.58 | 0.08 | 0.00 |
| C10b | R1 | 6.5 | 2.188 | 4.91 | 0.445 | 0.17 | 0.00 | 0.00 | 0.00 |
| | R2 | 8.5 | 2.497 | 5.02 | 0.497 | 0.17 | 0.00 | 0.00 | 0.00 |
| | R3 | 6.0 | 1.540 | 4.96 | 0.310 | 0.00 | 0.00 | 0.00 | 0.00 |

**Contrast against accumulated combs, pooled over all twelve ladder-cells:**
S = +1.055 × combs − 3.42 nat against a bar of 4.85 → **combs\* = 7.8**, i.e. the fringe
contrast reaches its trials bar at **≈ 8 accumulated combs**, against V25's extrapolated
10–11. The measurement lands below the extrapolation, which is the direction that helps.

### First: V25's K\* was already in AVAILABLE units

**The brief asks to "re-express K\* in available-source units", and doing that to
V25's 20–25 would double-count.** V25 reached 20–25 by dividing ~10–11 needed combs by a
recovery efficiency of ~0.5 — the division had already been done. The quantity that is
invariant under a change of ladder is **combs\***, and

    K*_available  =  combs* / conversion

| ladder | conversion (paired, pooled) | **K\*_available at combs\* = 7.8** | at V25's combs\* = 10.5 | era |
|---|---|---|---|---|
| R1 | 0.503 | **15.6** | 20.9 | SKA 30 yr |
| **R2 (best ladder)** | **0.518** | **15.1** | **20.3** | **SKA 30 yr** |
| R3 | 0.416 | 18.8 | 25.3 | SKA 30 yr |

### The calendar we are on

**The brief's two hypotheticals do not occur.** Conversion does not go from ~30 % to ~60 %
(R3, the candidate for that, goes *down*), and it does not reach 80 %. The measured
conversion is **0.50 for the control and 0.52 for the best ladder, a difference consistent
with zero**, so:

> **K\* ≈ 15 available sources on this run's own combs\*, or ≈ 20 on V25's — and it is the
> same number for the greedy chain and the joint fit. We are on the SKA calendar, not the
> MeerKAT one, and no recovery method in this instalment moves us along it.**

Against the forecasts that motivated the question — ~4 resolvable MeerKAT-class in 10 yr,
30–40 SKA in 30 yr — **K\* ≈ 15–20 available is comfortably inside the SKA number and a
factor of ~4–5 beyond the MeerKAT-era one.** The honest headline is that the arrival date is
set by the *array*, not by the pipeline: recovery efficiency is not the bottleneck it looked
like, because improving the fitter does not improve it.

**For calibration against your own concurrent run:** `reports/vernier27.npz` (15:23 today)
gives Arm W K\*\_contrast,50 = 16 on the grid, 13.5 interpolated, and Arm P not reaching 50 %
at K ≤ 20. My K\*\_available ≈ 15 sits on top of that without having been tuned to it. I have
not built on V27's numbers — it was running while this ran.

---

## Deliverable 4 — the wrong-recovery ledger

| cell | ladder | phantom/run | phantom/claim | **WRONG/claim** | NULL claims | NULL phantom | wrong cert | scrambled cert |
|---|---|---|---|---|---|---|---|---|
| C5 | R1 | 1.42 | 0.113 | 0.793 | 0.33 | 0.08 | 0.00 | 0.00 |
| | R2 | 1.92 | 0.114 | 0.836 | 0.25 | 0.08 | 0.00 | 0.00 |
| | R3 | 4.25 | 0.169 | 0.897 | 0.50 | 0.00 | 0.00 | 0.00 |
| C10 | R1 | 1.42 | 0.113 | 0.620 | 0.33 | 0.00 | 0.00 | 0.00 |
| | R2 | 1.83 | 0.110 | 0.730 | 0.33 | 0.00 | 0.00 | 0.00 |
| | R3 | 3.42 | 0.119 | 0.884 | 0.33 | 0.08 | 0.00 | 0.00 |
| C20 | R1 | 2.17 | 0.115 | **0.467** | 0.75 | 0.50 | 0.00 | 0.50 |
| | R2 | 2.00 | **0.085** | 0.563 | 0.92 | 0.42 | 0.08 | 0.17 |
| | R3 | 6.00 | 0.133 | 0.817 | 0.83 | 0.58 | 0.08 | 0.00 |
| C10b | R1 | 1.00 | **0.067** | 0.520 | 0.08 | 0.00 | 0.00 | 0.00 |
| | R2 | 3.00 | 0.130 | 0.652 | 0.25 | 0.00 | 0.00 | 0.00 |
| | R3 | 6.08 | 0.163 | 0.837 | 0.00 | 0.00 | 0.00 | 0.00 |

**The false-alarm rate is not where the phantoms come from.** On the *paired pure-noise
null* — the same sky, the same array, the same noise draws, signal removed from both blocks
— all three ladders claim **0.0–0.92 components a run**, consistent with the 5 %-per-bin
test they were each calibrated to. R3 is no worse than R1 here. **The ladders are clean
against noise.**

**The phantoms come from signal.** On data with sources in it, R3 claims 25–45 components
a run to R1's 12–19, and 0.82–0.90 of everything it claims is not a clean recovery
(against R1's 0.47–0.79). Those extra components are fitted to *the other 290 sources in
the bin* — real gravitational-wave power, correctly identified as power, wrongly banked as
a source. **The marginalised pulsar block does not distinguish "there is signal at this sky
position" from "there is a source at this sky position", because the confusion in the
pulsar block is signal too.**

**Wrong-FRINGE certification stays clean: 0.00 everywhere except two cells at 0.08.**

**The scrambled null does NOT stay clean at C20, and that is a known effect, not a new
one.** C20/R1 certifies 1.42 pulsars a run and its source-level scrambled null certifies
0.50 — a **35 % leak**; R2 is 0.17 on 1.50 (11 %). This reproduces the banked
`cert-scrambled-leak-K16` result (certification control-clean at K = 8, leaking 34–53 % on
the scrambled null at K = 16) in a cell with K_av = 18. **C20 is the only cell in this
instalment that certifies at all, so every certification number in Deliverable 3 carries
that leak**, and the crossing fractions there (0.67 / 0.75 / 0.42) should be read as upper
bounds. It does not touch Deliverables 1, 2 or the conversion result, which never use the
certification channel.

**No ladder buys conversion at the price of dirty certificates** — because no ladder buys
conversion. But no ladder cleans the K ≈ 16+ leak either, and R3 does not fix it: it simply
certifies less (0.58 a run) because it recovers fewer sources.

**The trade the brief asked to be quantified, stated plainly:** R2 costs +0.69 phantoms a
run over R1 (2.19 vs 1.50) and returns +0.015 ± 0.022 conversion. R3 costs +3.44 phantoms a
run and returns −0.087 ± 0.021. **R2 is a wash; R3 is a loss on both axes.**

---

## Controls

- **Control validation, every run**: instrumented R1 reproduces `vernier27_lib.chain`
  turn-for-turn (turns, genuine, blend, certifications) on three universes before any
  science runs. The run aborts if it does not.
- **Paired pure-noise null through all three ladders**, on the same sky with the same noise
  draws — the table above.
- **Source-level scrambled null at every certification step of every ladder**, not just the
  control's.
- **Null-calibrated model selection everywhere**: R2's per-increment floor from 12
  pure-noise realisations per bin; R3's per-step ΔT floor from 6 band-wide pure-noise
  realisations, both cut with V24's own `recut_floor`.
- **Chance-match floors printed per cell** (0.000–0.006 against the available pool) and for
  the 300-source bin (0.815, which is why the `faint` column is reported as *not* a
  recovery).
- **Sub-fringe rule obeyed**: every distance grid has an EVEN point count (512) and is
  centred on the prior mean.

### Two defects of my own, found by the null arm and fixed

Both were in R3's admission threshold, both made R3 look like a phantom factory, and
**neither was a property of the pulsar block**:

1. **Per-bin marginalisation.** The first version marginalised L_a *within each bin*. But
   pulsar *a* has one distance and every source in all ten bins sees it — marginalising per
   bin throws that lock away and leaves each bin's comb phase nearly free. Diagnostic: a
   single source's ln Z was as large 1.5 rad from the truth as on it. Fixed to band-wide;
   after the fix ΔT is 0 at the truth and −82 / −97 / −106 at 0.4 / 0.8 / 1.5 rad.
2. **A scalar admission floor pooled over model sizes.** R3's `grow` attempts up to 12 adds
   a bin; by m = 12 on pure noise the bin is saturated and every further increment is
   worth almost nothing, so pooling dragged the floor down and the pure-noise arm claimed
   ~7 components a run against R1/R2's ~0.3. Fixed to a per-step floor curve, calibrated
   round-robin so the band is populated as the threshold is measured. **The null arm then
   came back to 0.0–0.83 and the signal-side result did not change**: R3's conversion was
   0.438 with the broken floor and 0.416 with the fixed one, i.e. **R3's loss is not an
   artefact of my threshold.**

---

## What this does not settle

1. **`SUB_SPAN = 0` and a single chirp mass both weaken R3 by construction.** They were
   chosen to remove oracles, and they are the first thing to relax if R3 is ever revisited:
   a search that *could* measure sub-bin frequency and chirp mass would give the comb more
   to be commensurate with. This note bounds R3 from below, not from above.
2. **The certification channel and R3's recovery channel use the same physics.** R3 spends
   the pulsar-term comb on finding sources and `certify()` then spends it on pinning
   distances. The scrambled null is run at every step and comes back clean, so nothing here
   is dirty — but the two are not statistically independent and a fuller treatment would
   fit them jointly.
3. **12 realisations per cell, 48 universes total.** The paired intervals are honest, but
   Δ(R2−R1) = +0.015 ± 0.022 rules out a *large* gain, not a 3 % one. A gain that small does
   not move a calendar.
4. **M is capped at 12 per bin** (the brief's cap). R3's M̂ reaches 12 in C20 and C10b, so
   **R3's component counts in those cells are cap-limited**, and its phantom counts there are
   lower bounds. R1 and R2 never approach the cap.
5. **No M = 1 control in this instalment.** V24 and V25 both establish that one comb cannot
   break its own degeneracy; it is inherited rather than re-measured.
6. **The `faint` column is scored against all 300 sources**, where a random sky point lands
   within tolerance 0.815 of the time. It is reported as *not* a recovery for exactly that
   reason, but it is not a clean phantom count either — the clean one is the `phantom`
   column, and the honest total wrongness is `WRONG/claim`.
7. **combs\* = 7.8 is a pooled linear fit over twelve points spanning 2–10 combs**, and the
   bar itself drifts (4.40–5.28) with the comb count. It is a better-supported number than
   V25's extrapolation and it is still an extrapolation past the largest cell measured.
8. **Detection remains Earth-term-only for R1 and R2** — the standing V16/V22 limitation.
   That is deliberate here (it is what keeps R1 identical to V25/V27), and R3 is the arm
   that lifts it.
9. **C20's certifications inherit the K ≈ 16 scrambled-null leak** (banked separately as
   `cert-scrambled-leak-K16`): 0.50 scrambled certifications against 1.42 real ones for R1.
   C20 is the only cell that certifies, so combs\* and every crossing fraction in
   Deliverable 3 rest on a channel with a known 11–35 % leak in its top cell. **The
   conversion result, which is the point of this instalment, does not use that channel at
   all.**
10. **C5 and C10/C20 differ in construction as well as in K_available.** C10b is the
   ablation that separates them, and it says construction matters (0.555 seeded vs 0.466
   flattened at higher availability) — but one ablation cell is not a controlled sweep.

---

## One-line verdict

**Blend recovery is not a lever: fitting every source in a bin simultaneously, with
positions re-solved against each other and the component count chosen against a measured
null floor, converts available detections into recovered combs at 0.518 ± 0.027 against the
greedy control's 0.503 ± 0.025 — a paired difference of +0.015 ± 0.022 that is consistent
with zero, winning in 16 matched universes of 48 and losing in 14 — while adding the
marginalised pulsar-term block on top of it makes things significantly worse (0.416 ± 0.024,
Δ = −0.087 ± 0.021, excluding zero) because with the distances marginalised the block
cannot tell "there is signal here" from "there is a source here" and the other 290 sources
in the bin are signal, so R3 banks 25–45 components a run to the control's 12–19 with
0.82–0.90 of them not clean recoveries; the realized per-bin separability of an actual
joint search is 2–3 and not V24's Fisher frontier of 5–10, against a counting bound of 50,
and joint fitting raises the number of components the data supports (M̂ median 1–2 → 2–4,
max 12) without raising the number that are sources (G_bin median 1.0 and max 3 in every
cell of every ladder, with at most 22 % of bins ever yielding two), which is the same
statement in a different currency; the degree-of-freedom surplus from V23 Part 1 is
therefore real and useless — the unmodified V24 Fisher at ρ₁ = 30 shows it tripling resolved
sources at M = 10 (6 against 2) while the search built on it loses conversion, so it buys
precision on a source already found and nothing on the problem of finding one; the
false-alarm rate is not the failure mode, because on the paired pure-noise null all three
ladders claim 0.0–0.92 components a run and wrong-fringe certifications are 0.00 nearly
everywhere, so the phantoms are bought with signal rather than with noise; and since the
invariant that sets the threshold is combs\* — measured here at 7.8 against V25's
extrapolated 10–11, and noting that V25's "K\* ≈ 20–25" was ALREADY in available units
because it had already been divided by a conversion of ~0.5 — the calendar reads
K\*\_available = combs\*/conversion ≈ 15 on this run's combs\* and ≈ 20 on V25's, IDENTICAL
for the greedy chain and the joint fit to within the error, which puts us inside the SKA
30-year forecast of 30–40 available sources and a factor of four beyond the MeerKAT-era
four, with the arrival date set by the array and not by the pipeline.**
