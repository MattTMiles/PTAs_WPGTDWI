# LOTTERY-2 — the REALISTIC LOUDNESS PRIOR, and what it does to `P_cond = 0`

**Agent LOTTERY-2, cronus, 2026-07-29. CPU only — no GPU touched, no ACCRE job submitted.**
Staged, not committed. Script `hpc_harbor/lottery/lottery2_realistic_pop.py`; bank
`reports/lottery2_realistic_pop.npz`. Add-list at the foot.

> **THE MISSION.** GLACIER Stage-2a banked **`P_cond(onset) = 0 EXACTLY`** — under the
> campaign's fixed-N power-conserving draw the brightest member of an NG15-normalised sky
> cannot exceed `log10 h = −13.867`, so the `−13.25` ladder rung is not merely unlikely, it
> is *outside the model class*. This report replaces that draw's hard ceiling with a
> **realistic loudness prior** and emits **P(sky above rung)** for the ladder's exact rungs
> {−13.87, −13.5, −13.25}, with the eccentricity joint carried at the Arm-A boundary.

---

## 0. WHAT THE CEILING ACTUALLY WAS — four suppressions, not three

`glacier_pop.draw_population` fixes `N = 256`, draws a bare lognormal strain *shape*
(`sigma_logh = 0.6`), throws the drawn masses and redshifts away, and then **renormalises**
so that `sum_i h_i^2 == A^2 (3/4)[x(flo) − x(fhi)]` with `log10 A = −14.6` **exactly**. The
ceiling `sqrt(target) = −13.8674` is what you get when one source carries 100 % of a
*constant* total. The brief named three suppressed tails; the arithmetic has four, and the
fourth is the one that does the work:

| | suppression | restored as |
|---|---|---|
| **T1** | `N` fixed at 256 | `N ~ Poisson(256)` |
| **T2** | loudness is a bare lognormal shape — masses/redshifts drawn then discarded, so there is no `N(>h) ~ h^-3` volume tail | a phenomenological SMBHB synthesis in (M, q, z) with the GW-driven residence time |
| **T3** | `log10 A` pinned at the point estimate | `log10 A ~ N(−14.6, sigma_A)` |
| **T-pin** | **the REALISED sky power is constrained to equal the ENSEMBLE mean** | the calibration moved to the ensemble: `lambda·<h^2> = P_target(A)` |

**T-pin is the ceiling.** T1 and T3 do not remove it — they only move it (T3) or leave it
(T1). §2's ablation measures exactly that.

### 0.1 The realistic prior, declared

Phenomenological merger-rate density (the standard form; **[VERIFY]** against
Middleton+2016 MNRAS 455 L72 / Chen, Sesana & Conselice 2019 MNRAS 488 401 / the NG15
astro-interpretation Agazie+2023 ApJL 952 L37):

```
d^2 n / (dlog10 M dz)  ~  (M/1e7)^(-alpha) exp(-M/Mstar) (1+z)^beta exp(-z/z0)
```

fiducial `alpha = 0`, `log10 Mstar = 8.0` (**TOTAL mass**), `beta = 2`, `z0 = 2`, all
**declared and swept, never fitted**. Total mass is pushed forward to chirp mass through a
declared `p(q) ~ q^0` on `[0.25, 1]` (`kappa(q) = q^(3/5)/(1+q)^(6/5)`). Band sources per
unit `ln f` carry the GW-driven residence time and the volume element:

```
dN/(dlog10 Mc dz dlnf) ~ Phi (dV_c/dz) (1+z)^-2 Mz^(-5/3) f^(-8/3)
h = (8/sqrt10) (G Mz)^(5/3) (pi f)^(2/3) / (c^4 d_L)
```

**[VERIFY]** the strain convention (Sesana, Vecchio & Colacino 2008 MNRAS 390 192): this is
the `h` whose quadrature sum over the band *is* `int h_c^2 dlnf` — i.e. the same `h` the
campaign pins, which is why the rungs mean the same thing on both sides of this report.
**[VERIFY]** `sigma_A`: NG15 at `gamma = 13/3` gives `A = 2.4 (+0.7/−0.6) e−15`
(Agazie+2023 ApJL 951 L8) → `log10 A = −14.62 +0.11/−0.12`. **The repo's declared
`SIGMA_NG15_DEX = 0.05` is ~2.2× narrower than NG15's own posterior.** Both are carried;
0.11 is the primary.

The rate normalisation `ndot0` is **not an input** — it is fixed by
`lambda·<h^2> = P_target(A)`, which is exactly how the field calibrates a synthesis to the
measured GWB. `ndot0` and `A` are degenerate; this removes the degeneracy the honest way.

### 0.2 Why the tail is exact rather than sampled

Band sources are a Poisson process, so the count above any `H` is Poisson with mean
`mu(H) = lambda·P(h >= H)` and

> **`P(brightest >= H) = 1 − exp(−mu(H))`** — exact at any depth, no MC noise.

`mu(H)` is quadrature on the `(log10 Mc, z)` grid with the `f`-integral in closed form, so
a `1e−8` tail is as accurate as a `1e−1` one. **The ≥ 1e5-draw MC (2×10⁵ draws) is run
beside it** as a cross-check and to carry the conditioning variants — and the two agree
wherever the MC can resolve anything at all (g2). *This is the honest division of labour:
no affordable MC resolves `3.6e−7`; the mission's "1e5+ draws" is met and is the control,
not the estimator.*

---

## 1. GATES

| gate | what it proves | result |
|---|---|---|
| **g0** | the incumbent draw is reproduced exactly | **PASS** — ceiling −13.867396 (banked −13.867396); p50/p90/p99.9 = −14.147/−13.978/−13.873 vs banked −14.147/−13.979/−13.874; `P(>= −13.9)` 0.0142 vs banked 0.0139 |
| **g1** | power identity `lambda·<h^2> == P_target` | **PASS** — relative residual **0.00e+00** |
| **g2** | the 2×10⁵-draw MC reproduces the analytic Poisson tail | **PASS** at every rung it can resolve (−14.98: 5.53e−2 vs 5.57e−2; −14.65: 5.65e−3 vs 5.37e−3; −14.41: 1.12e−3 vs 1.05e−3). The ladder rungs are reported **UNRESOLVABLE** (0.07–6.4 expected hits in 2×10⁵) — stated, never graded |
| **g3** | `mu(rung)` invariant to the integration box and grid | **PASS** on 5 of 6 graded rows (grid ×2, ×4: 3e−5; `z_lo` 1e−8: 2e−4; `z_max` 8: 8e−3; `mc_lo` 5.0: 1.6e−4). **`mc_lo 7.0` FAILS at 5.5e−2** — see below |
| **g4** | *scope disclosure, not pass/fail* | banked: **where** the sources above each rung live |

**g3's one failure is a real 5.5 % dependence, not a numerical one:** truncating the mass
function at `Mc = 10^7` removes 5.5 % of `mu` at the rungs. So ~5 % of the fiducial tail
comes from chirp masses below `10^7` — quoted, not hidden. Every other box edge is clean.

**g3 also caught a live bug and is the reason a default moved.** The first build used a
*linear* `z` grid and `z_lo = 1e−4`; g3 failed at 1.08 (i.e. `mu` moved by more than its own
value). The brightest-member tail is **near-field dominated**, so a linear grid starves it
of nodes and a `z` floor truncates it outright. The grid is now log-spaced with
`z_lo = 1e−7` and converged to < 0.1 %; the truncated box is retained as an **ungraded
regression row** (still fails at 0.93, as it must).

**g4 — the scope disclosure, and it governs how §3 may be worded.** Median `d_L` of the
sources that clear each rung, at the fiducial:

| rung | `d_L` 5/50/95 % [Mpc] | `log10 Mc` 5/50/95 % |
|---|---|---|
| −13.87 | 0.047 / **0.36** / 2.05 | 7.60 / 8.09 / 8.45 |
| −13.50 | 0.021 / **0.16** / 0.89 | 7.60 / 8.09 / 8.45 |
| −13.25 | 0.012 / **0.089** / 0.50 | 7.62 / 8.09 / 8.45 |
| −12.75 | 0.004 / **0.028** / 0.16 | 7.62 / 8.09 / 8.45 |

**At the fiducial mass function, a source at the ladder's rungs must sit inside the Local
Group.** That is not a defect of the method; it is what `Mstar = 10^8` *total* means. It is
also why the sweep (§6), not the fiducial, is the deliverable — at the sweep's loudest
corner the same sources sit at **3–120 Mpc**, which is an ordinary answer.

---

## 2. (a) P(brightest member >= rung) — THE DELIVERABLE

Band `(−8.7, −7.5)` (the banked scan's band), `sigma_A = 0.11`, 2×10⁵ draws.
Fixed-class ceiling = **−13.8674**.

| rung | σ_A above ceiling | **M0** incumbent | **M1** +Poisson N | **M2** +NG15 width | **M3** both | **M4 realistic** | odds |
|---|---|---|---|---|---|---|---|
| **−13.87** | −0.02 | 2.21e−4 | 1.75e−4 | 4.70e−2 | 4.70e−2 | **2.59e−5** | 1 in 38 600 |
| **−13.50** | +3.34 | **0** | **0** | 6.75e−6 | 6.76e−6 | **2.01e−6** | 1 in 497 000 |
| **−13.25** | +5.61 | **0** | **0** | 8.09e−11 | 8.13e−11 | **3.58e−7** | 1 in 2 800 000 |
| −13.90 *(rung 1 as realised)* | −0.30 | 1.44e−2 | 1.45e−2 | 6.90e−2 | 6.90e−2 | 3.19e−5 | 1 in 31 400 |
| −12.75 *(A-SKY robust line)* | +10.16 | **0** | **0** | 5.2e−27 | 5.2e−27 | **1.13e−8** | 1 in 8.9e7 |

**Four readings, in order of how much they change the campaign's language.**

1. **The ladder's top rung was never actually above the ceiling.** `−13.87` sits **0.0026
   dex BELOW** `−13.8674`, so it is reachable even in the incumbent class (M0: 2.2e−4).
   The GLACIER text realises rung 1 as a tail cut at `−13.9`, where M0 gives 1.4e−2 —
   consistent with the banked `P_cond = 1.4e−2`. **Only rungs 2 and 3 were ever
   `P_cond = 0`.**
2. **Poisson `N` does nothing.** M0 → M1 is flat at every rung and still exactly zero above
   the ceiling. At `lambda = 256` the fluctuation is ±6 % on a quantity the renormalisation
   immediately divides out. **T1 was not a load-bearing suppression.**
3. **The NG15 width is what *moves* the ceiling; the realistic loudness is what *removes*
   it.** M2 makes `−13.25` finite (8e−11) purely by letting `A` fluctuate up — a 5.6-σ
   excursion in `A` that drags the ceiling with it. The realistic prior gives **3.6e−7, a
   factor 4 400 larger**, because it does not need `A` to be wrong: it needs one nearby
   massive binary.
4. **The incumbent draw is not uniformly conservative.** Just below its ceiling it is
   *over-generous* — M2/M3 give 4.7e−2 at `−13.87` against the realistic prior's 2.6e−5, a
   factor 1 800 the other way, because pinning the total forces near-ceiling configurations
   to be common. **The fixed class is 3 orders too generous immediately below its ceiling
   and infinitely too strict immediately above it.** Neither error is small, and they point
   in opposite directions.

---

## 3. (b) THE `P_cond = 0` STATEMENT, RE-SCOPED — both results, side by side

> ### FIXED CLASS (unchanged, and still correct as stated)
> **`P_cond(−13.25) = 0 EXACTLY`.** Not small — outside the model class. A `−13.25` member
> alone carries ~16× the total NG15 band power; the hard ceiling is `−13.867`. This remains
> the right statement *about the power-conserving draw*, and GLACIER's framing sentence
> ("within the NG15-consistent class, no sky contains an onset-class source") is a true
> statement about that class.
>
> ### REALISTIC PRIOR (new)
> **`P(brightest >= −13.25) = 3.6e−7` — about 1 in 2.8 million skies.**
> **`P(brightest >= −13.5) = 2.0e−6` — 1 in 500 000.**
> **`P(brightest >= −13.87) = 2.6e−5` — 1 in 39 000.**
> Across the declared synthesis sweep (§6) the `−13.25` rung spans
> **2.2e−7 → 1.07e−4, i.e. 1 in 4.6 million to 1 in 9 400.**
>
> ### THE RE-SCOPED SENTENCE
> ***`−13.25` is RARE, not impossible. The zero was a property of the power-conserving
> draw — specifically of pinning the realised sky total to the ensemble mean — and not a
> property of nature. But "rare" here is 1 in 10⁴–10⁷, and it is bought by putting a
> massive binary unusually close: the odds and the required distance must always travel
> together.***

**The tension numbers on the ladder's figure axis need a correction that travels back to
GLACIER.** The banked axis reads **+7.8…8.0 σ at `−13.5`, +12.5 σ at `−13.25`** — those are
computed against the repo's declared `SIGMA_NG15_DEX = 0.05`. Against **NG15's own posterior
width (0.11 dex)** the same excursions are **+3.3 σ and +5.6 σ** (this report's
`sigma_above_ceiling` column reproduces them to 0.05 σ, an independent check of the banked
arithmetic). **The tension is real and large either way, but the printed number is a factor
2.2 sensitive to a width the repo declared rather than measured.** The ladder's conclusion
does not change; the axis label should.

### 3.1 The scope that must be quoted with the odds

The near-field domination (g4) is not a footnote. At the fiducial, a `−13.25` source sits at
a median `d_L` of **0.089 Mpc** — inside the Local Group, where a smooth comoving rate
density is an extrapolation and the true host count is a handful of galaxies, not a
continuum. Excluding the near field outright:

| rung | all d | d > 1 Mpc | d > 5 Mpc | d > 10 Mpc | d > 20 Mpc |
|---|---|---|---|---|---|
| −13.87 | 2.28e−5 | 3.96e−6 | 1.00e−7 | 6.03e−9 | 9.4e−11 |
| −13.50 | 1.77e−6 | 6.55e−8 | 1.88e−10 | 1.87e−12 | 1.9e−15 |
| −13.25 | 3.15e−7 | 2.57e−9 | 7.55e−13 | 1.11e−15 | 0 |

**At the fiducial mass function, excluding the Local Group costs the `−13.25` rung a factor
120, and excluding 5 Mpc costs it a factor 4×10⁵.** At the sweep's loudest corner
(`alpha = −0.5`, `log10 Mstar = 9.0` read as CHIRP) the same curve is nearly flat —
9.4e−5 → 9.3e−5 → 8.5e−5 at d > 0/1/5 Mpc — because there the sources are genuinely at
23 Mpc (median). ***So the honest form of the headline is conditional on the mass function:
with a `10^8`-total cutoff, `−13.25` is not "rare", it is "in the Local Group"; with a
`10^9`-chirp cutoff it is a plain 1-in-10⁴ sky at an ordinary distance.*** The rung
probability and the mass-function cutoff are the same statement.

### 3.2 Three ways to ask "is the sky NG15-consistent", and they differ

The MC carries the realised power, so the conditioning question can be posed three ways
(each at ±2 σ_A in power; C2/C3 keep 95.5 % of draws):

| rung | **C1** ensemble-calibrated | **C2** realised TOTAL pinned | **C3** realised BACKGROUND pinned |
|---|---|---|---|
| −13.87 | 2.50e−5 | 5.24e−6 | 2.62e−5 |
| −13.50 | 5.0e−6 | **0** | 5.24e−6 |
| −13.25 | 0 *(MC-unresolvable)* | **0** | 0 *(MC-unresolvable)* |

- **C2** (the realised total must itself sit inside the NG15 posterior) is the incumbent pin
  with a width — and it re-suppresses the tail by ~5× at `−13.87` and kills `−13.5`
  outright at MC resolution. **The incumbent draw is the `sigma_A → 0` limit of C2**, so C2
  is the fair like-for-like: it says the ceiling was *mostly* right if you insist the
  realised sky total equals the measurement.
- **C3** (the *unresolved* sum must sit inside the posterior) leaves the tail untouched.
  **This is the defensible reading: a single resolvable source is not HD-correlated, so it
  is not what NG15's cross-correlation amplitude measures.** A loud CW does not have to fit
  inside the GWB budget, because the GWB budget was measured with a statistic that does not
  see it.
- ***The gap between C2 and C3 is the whole re-scope.*** It is not about population
  synthesis at all; it is about which statistic the NG15 amplitude constrains.

---

## 4. (c) THE ECCENTRICITY JOINT AT THE ARM-A BOUNDARY

The Arm-A boundary `e*(h)` is **read from `reports/gen_armA_surface.npz::thr`**, never
retyped. At `(3+13, T = 30)` the measured box is `h ∈ {−13.0, −12.75, −12.5}` plus the
banked census column at `−13.25`:

| h | e* (3+13, T30) | provenance |
|---|---|---|
| −12.50 | 0.3 | measured |
| −12.75 | 0.3 | measured |
| −13.00 | 0.5 | measured |
| **−13.25** | **0.5** | banked CHORUS re-cut (census column) |
| **< −13.25** | **unmeasured** | **the box bottoms at −13.25** |

**Both fainter rungs sit below the measured box.** `e*` there is not known; monotonicity in
`h` says only `e* >= 0.5`, and it is entirely possible no `e` suffices. **Every joint at
`−13.87` and `−13.5` is therefore an UPPER BOUND**, flagged as such in the bank
(`joint_flag_key`).

Eccentricity prior: the LOTTERY tier-1 family `e ~ lnN(mode 0.3, width w)`, clipped [0, 0.9],
evaluated exactly (erfc) — at `w = 0.1` the `e >= 0.5` tail is 5.1 σ and a 2×10⁶-draw MC
returns a bare 0.0000, which would silently zero the joint.

| | `P(e >= 0.5)` | `P(e >= 0.3)` |
|---|---|---|
| w = 0.1 | 2.75e−7 | 0.540 |
| w = 0.2 | 9.28e−3 | 0.579 |
| w = 0.4 | 0.190 | 0.655 |

### 4.1 Single-member joints — SKY-MEDIAN-CLASS ONLY

`P(brightest >= rung AND e >= e*(h))`, `e ⊥ h` (see the caveat below):

| rung | e* used | w = 0.1 | w = 0.2 | w = 0.4 |
|---|---|---|---|---|
| −13.87 *(upper bound)* | 0.5 | 7.1e−12 | 2.4e−7 | 4.9e−6 |
| −13.50 *(upper bound)* | 0.5 | 5.5e−13 | 1.9e−8 | 3.8e−7 |
| **−13.25** *(banked census)* | **0.5** | **9.8e−14** | **3.3e−9** | **6.8e−8** |

> **THE A-SKY CAVEAT, CARRIED.** Every row above is a **single-member 3+13** claim, and
> GENERALISE §4b measured those to be a **sky lottery — 2–4 of 8 skies, the old `e = 0.5`
> edge at `h = −13.0` included (4/8)**. These joints are therefore quotable **only as
> sky-ensemble-median-class statements with the A-SKY band attached, never as per-sky
> probabilities.** A per-sky quote would multiply a ~10⁻⁷ population probability by a
> coin-flip nobody has marginalised, and the coin flip is the larger uncertainty.

### 4.2 The robust line — the one that survives sky marginalisation

A-SKY's sole surviving external sentence is **5+11 structure, `e >= 0.3`, `h >= −12.75`
(6/8 skies, median 1.77)**. Its `e*` is read from the bank as **0.3 (measured, `k = 5`)**.
Its loudness cost under the realistic prior:

- **necessary condition:** `P(brightest >= −12.75) = 1.13e−8` — **1 in 8.9×10⁷.**
- **with the 5+11 structure proxy** (`N >= 5` sources above `−12.75`, the Poisson count):
  `log10 P = −42.1`, and the joint with `P(e >= 0.3) ≈ 0.55` is `log10 P_joint ≈ −42.3`.
  At the sweep's loudest corner it is `log10 P = −29.7` — **still 30 orders down.**

> ### THE ROBUST LINE IS THE ROBUST *STATEMENT*, NOT A REACHABLE SKY.
> **The only eccentric switch-on claim that survives sky marginalisation requires a
> population the NG15-calibrated prior produces with probability ~10⁻⁴². The structure
> requirement, not the eccentricity and not the sky, is what is unreachable:** the
> necessary single-source condition alone costs 1 in 9×10⁷, and demanding *five* such
> sources costs another 34 orders. This is not a marginal call and no sweep corner rescues
> it. *The A-SKY survivor is the right thing to quote about the physics, and it describes a
> population that this prior says does not occur.*

**Independence caveat, stated because it bounds the joint.** `e` is drawn independently of
`(M, z, f)`, as in the incumbent draw. Physically the correlation runs the wrong way for
these numbers: an eccentric binary has a **shorter** residence time at fixed `f` (fewer in
band) and spreads its power across harmonics (**lower** `h` in the fundamental). **Both
effects suppress the eccentric-loud joint, so every number in §4 is an upper bound on that
count as well.**

---

## 5. WHAT THE `1e5+` DRAWS ACTUALLY BOUGHT

2×10⁵ draws, split-Poisson sampled (all sources above `h_split = −16.18`, `mu = 197`
explicit per draw; the faint remainder carried as its conditional-mean power with Poisson
fluctuation). The MC **confirms the analytic tail wherever it can resolve it** (g2, three
rungs at 1e−3…5e−2, agreement within 4 σ or 5 %) and **cannot resolve any ladder rung**
(0.07–6.4 expected hits). It is the control and the carrier of the C1/C2/C3 conditioning
masks; the Poisson-process identity is the estimator. *Reporting the MC's `0.000e+00` at
`−13.25` as the answer would have been the easy mistake here — it is 0 hits out of 2×10⁵ at
a true rate of 1 in 2.8 million.*

---

## 6. THE SENSITIVITY BAND — the deliverable is the band, not the fiducial

`P(brightest >= rung)` over `alpha ∈ {−0.5, 0, +0.5}` × `log10 Mstar ∈ {8.0, 8.5, 9.0}` ×
{Mstar read as TOTAL / as CHIRP}:

| rung | min | **fiducial** | max | spread |
|---|---|---|---|---|
| −13.87 | 1.57e−5 | **2.59e−5** | 7.80e−3 | **500×** |
| −13.50 | 1.22e−6 | **2.01e−6** | 6.02e−4 | **490×** |
| −13.25 | 2.17e−7 | **3.58e−7** | 1.07e−4 | **490×** |

The minimum is `(alpha = +0.5, Mstar = 10^8 total)`; the maximum is
`(alpha = −0.5, Mstar = 10^9 chirp)`. **The chirp-vs-total reading of `Mstar` alone is a
factor ~4.4 at every rung** — a convention, not a measurement, and the single largest
tractable uncertainty in this report after the mass cutoff itself. Two smaller axes:

- **strain-amplitude convention** (`8/sqrt10` vs `h0 = 2(GMz)^{5/3}…`): factor **0.79** on
  every rung. Small, and it moves everything together.
- **NG15 posterior width**: `sigma_A` = 0 / 0.05 / 0.11 / 0.15 gives `−13.25` =
  3.15e−7 / 3.23e−7 / 3.58e−7 / 3.99e−7. **A 27 % effect — the width matters for the
  *tension axis* (§3), not for the rung probability**, once the realistic loudness tail is
  present. That is the cleanest statement of what T3 was and was not doing.

---

## 7. WHAT THIS SETTLES, AND WHAT IT HANDS BACK

1. **`P_cond = 0` is re-scoped, not overturned.** It is exactly true of the power-conserving
   draw and false of nature. The replacement numbers are
   **2.6e−5 / 2.0e−6 / 3.6e−7** at `−13.87 / −13.5 / −13.25`, with a **~500× band** across
   the declared synthesis sweep. `−13.25` is **rare (1 in 10⁴–10⁷), not impossible.**
2. **The load-bearing suppression was the PIN, not `N` and not `A`.** Poisson `N` changes
   nothing; the `A` width moves the ceiling but still needs a 5.6-σ excursion; only dropping
   the realised-power constraint removes the ceiling. **And the pin errs in both directions**
   — 1 800× too generous just below the ceiling, infinitely too strict just above it.
3. **The re-scope's real content is a statistics question, not an astrophysics one** (§3.2):
   NG15's amplitude constrains the HD-correlated background (C3), not the realised total
   including a resolvable source (C2). The incumbent draw implements C2 at zero width.
4. **Every rung number must be quoted with its distance** (g4, §3.1). At the fiducial mass
   function the `−13.25` tail lives inside the Local Group and loses a factor 120 the moment
   the Local Group is excluded; at the loudest sweep corner it lives at 23 Mpc and loses
   nothing. **The rung probability and the mass-function cutoff are one statement.**
5. **The A-SKY robust line is unreachable under this prior by ~42 orders** (§4.2), driven by
   the 5+11 *structure* requirement, not by eccentricity or sky. Handed back as the sharpest
   consequence: the only sky-marginalisation-surviving switch-on claim describes a population
   an NG15-calibrated prior does not produce.
6. **A correction that travels to GLACIER's figure axis:** the banked tension numbers
   (+7.8…8.0 σ, +12.5 σ) are computed against a declared `sigma_A = 0.05`; against NG15's own
   `0.11` they are **+3.3 σ / +5.6 σ**. Independently reproduced here to 0.05 σ. The
   conclusion stands; the printed number should say which width it used.
7. **Two [VERIFY] flags that need a literature eye before any of this is quoted outside**:
   (i) the phenomenological `(alpha, Mstar, beta, z0)` and whether `Mstar` is total or chirp
   mass — a factor 4.4 per rung; (ii) the `8/sqrt10` strain convention — a factor 0.79.
   Both are declared in the module header and swept, not assumed away.
8. **Known limits, stated:** `e` is drawn independently of `(M, z, f)` (§4, upper bound);
   ~5 % of the fiducial tail comes from `Mc < 10^7` where the mass function is an
   extrapolation (g3's one failure); the eccentric residence-time and harmonic-power
   corrections are not modelled; nothing here is a GPU measurement — it is a prior, and the
   certification side of every joint is still the banked Arm-A/A-SKY grade.

```
add-list (stage, do not commit):
  LOTTERY2.md
  reports/lottery2_realistic_pop.npz
  hpc_harbor/lottery/lottery2_realistic_pop.py
```

**SIZE GLANCE.** `LOTTERY2.md` 24 KB · `reports/lottery2_realistic_pop.npz` **628 KB**
(32 keys, `savez_compressed`) · `hpc_harbor/lottery/lottery2_realistic_pop.py` 48 KB.
**Total staged ≈ 0.7 MB.** The bank was 4.6 MB on first write; the 2×10⁵-draw raw power
columns are banked as **histograms** (their only use, the C2/C3 masks, is already banked as
`cond_table`) and `hmax_mc` is float32 — no readback anyone needs was dropped. No
`*_results/` directory was created; there are no per-realisation banks, because there are no
realisations — this is CPU quadrature plus one MC, both re-derivable from the module in
minutes.

**STOP.** Mission complete: L1 population upgrade built and gated (g0–g4), L2 (a)/(b)/(c)
delivered with the sweep band, the A-SKY sky-lottery caveat carried on every single-member
joint and the structure-assisted survivor quoted as the robust line. No GPU used, no job
submitted, nothing committed. No further work without a fresh brief.
