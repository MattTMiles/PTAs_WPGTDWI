# VERNIER-29 — the VLBI targeting ranking (the Deller deliverable)

**EXPLORATORY. NOT pre-registered.** Toy GW sky, **real array**, cronus, CPU/numpy, **~5 min**.
Clean room: imports only `vernier24` / `vernier25` (and through them `vernier_toy` / `vernier6` /
`vernier16`), plus two real data products already tracked in this repo.

Code `hpc_harbor/vernier/vernier29.py` · data `reports/vernier29.npz`.
**Staged, not committed.**

Reuses VERNIER-24's corrected-sky construction and incoherent-comb certification floor,
VERNIER-25's Arm K first-detection seeding, VERNIER-13's heterogeneity treatment (here
replaced by the array's *measured* heterogeneity), and VERNIER-23 Part 3's soft channel.

---

## Lead — the ranking is a distance cut, not a timing cut, and it is short

A parallax campaign of precision σ_π buys a **distance** precision σ_d = σ_π d², because
π[arcsec] = 1/d[pc]. At the best precision anyone has actually demonstrated on an MSP
(**13 μas**, Ding et al. 2024 on J2222−0137, PX = 3.723 ± 0.013 mas — a number that is in
the catalogue file this run reads), sub-parsec is reachable only out to **277 pc**, and
sub-λ/2 only out to **158 pc**. Across the 116-pulsar array that is **6 pulsars** and
**3 pulsars** respectively. At a more typical 30 μas it is 3 and 1.

**So the VLBI team's own expectation of "2–3 more sub-parsec MSP distances next decade" is
not a funding constraint — it is reproduced here from the catalogue as a geometric one, and
no amount of telescope time moves it.** The ranking below is therefore, to first order,
a list sorted by distance, with timing weight breaking ties.

And the honest framing (D2/ORACLE, and the R referendum): **none of this buys ignition.**
VERNIER-24/25 put the ignition fraction at 0.00 on every sky at K ≤ 16 and the NO-GO is
information-theoretic. Everything below is priced in products.

---

## Deliverable 0 — the real array

116 pulsars with **real** sky unit vectors (tracked feather `pos_t`), **real** residual rms
(`reports/anchor_data_forensics.npz`), and **real** distance priors with the method that
produced each one (`CW_transition/anchor_a0_priors.npz` — `get_distance.py`'s canonical
output against the Cornell parallax catalogue). T = 20 yr, 1/T = 1.584 nHz, ten bins
10.00–24.26 nHz, λ_GW(15 nHz) = 0.648 pc.

Residual rms: median **756 ns**, range 196–40434 ns.

**VERNIER-13's assumed heterogeneity is confirmed by the real array**: the measured spread
is σ_log10(1/rms) = **0.370 dex** against the 0.400 dex the series has been assuming. That
is a free check on nine instalments' worth of arrays, and it passes.

Distance-prior category census — the catalogue values this run uses:

| category | N | median d (pc) | median σ_d (pc) | median σ_d/d | sub-pc today |
|---|---|---|---|---|---|
| VLBI | 19 | 1075.3 | 82.6 | 0.082 | 2 |
| timing parallax | 50 | 1125.3 | 219.2 | 0.208 | 0 |
| Ṗb Shklovskii | 5 | 930.0 | 197.6 | 0.200 | 0 |
| DM model (YMW16) | 12 | 2173.5 | 434.7 | 0.200 | 0 |
| DM default (no parallax-class measurement) | 30 | 1495.0 | 300.0 | 0.200 | 0 |

**Exactly two pulsars have σ_d < 1 pc today: J0437−4715 (0.11 pc) and J2222−0137 (0.97 pc).
Exactly one is inside λ/2 = 0.324 pc.**

J0437−4715 enters as an existing anchor at 0.11 pc, per the brief, and is excluded from the
candidate list. **Discrepancy, stated:** the canonical prior file in this repo still carries
Deller et al. (2008), PX = 6.396 ± 0.054 mas → d = 156.3 pc, σ_d = **1.32 pc**. The 0.11 pc
figure is the Reardon-era value. Deliverable 5 prices the anchor both ways.

---

## Deliverable 1 — what VLBI can actually reach

σ_d = σ_π d². The reachable set at 13 μas, sorted by reachable σ_d:

| # | pulsar | d (pc) | σ_d now (pc) | category now | σ_d reachable (pc) | gain |
|---|---|---|---|---|---|---|
| 1 | J1630+3734 | 89.0 | 24.00 | DM default | **0.103** | ×233 |
| 2 | J0711−6830 | 106.0 | 21.20 | DM model | **0.146** | ×145 |
| 3 | J0437−4715 | 156.3 | 0.11 | VLBI | 0.318 | ×0.3 *(anchor — already better)* |
| 4 | J1737−0811 | 206.0 | 41.20 | DM model | **0.552** | ×75 |
| 5 | J1017−7156 | 256.4 | 78.90 | timing PX | **0.855** | ×92 |
| 6 | J2222−0137 | 268.6 | 0.97 | VLBI | 0.938 | ×1.0 *(already there)* |
| 7 | J2124−3358 | 322.6 | 57.23 | timing PX | 1.353 | ×42 |
| 8 | J0030+0451 | 328.9 | 5.41 | VLBI | 1.407 | ×3.8 |
| 9 | J1744−1134 | 416.7 | 17.36 | timing PX | 2.257 | ×7.7 |
| 10 | J1045−4509 | 437.3 | 218.67 | timing PX | 2.487 | ×88 |
| 11 | J0613−0200 | 476.2 | 136.05 | timing PX | 2.948 | ×46 |

11 pulsars can be driven below 3 pc; **6 below 1 pc; 3 below λ/2**. Beyond ~480 pc nothing
in the array can be brought inside 3 pc at any astrometric precision on the table.

**Data caveat that bears directly on rows 1, 2 and 4.** J1630+3734's 89 pc, J0711−6830's
106 pc and J1737−0811's 206 pc are all *feather* `pdist` values with **no parallax-class
measurement behind them**. For J0711−6830 and J1737−0811 the repo's own `get_distance.py`
output disagrees with the feather by **8.0×** and **10.5×** (845.6 pc and 2161.6 pc from
DM_YMW16). J1630+3734 has no `best_distances.txt` entry at all. Ten of the array's
non-parallax pulsars disagree between the two files by more than 2×. **Three of the five
best targets are "nearby" only according to a distance nobody has measured — which is an
argument for observing them, but it means their rank is provisional until the parallax
exists.**

---

## Deliverables 2+3 — the ranked table

For each pulsar: force its prior to sub-λ, recompute all four products, and take the change.
Baseline and perturbed passes **share their noise draws**, so every Δ is a prior response
and not a fresh noise realisation. 24 sky realisations, 8 accumulated combs at array
SNR 12 on the corrected γ = 1.56 sky with its honest confusion sea.

Products: **(a)** certification — `pin` is the fraction of realisations in which the pulsar
goes from not-known to known; **(b)** soft-channel depth Δf_marg (VERNIER-23 Part 3, no
threshold at all); **(c)** localisation — fractional reduction in the 1-σ sky-ellipse area of
the loudest source; **(d)** fractional reduction in σ(log₁₀ M_c).

Baseline products at today's priors: f_marg = 0.199 (K=1) / 0.183 (K=8) / 0.233 (K=16);
f_cert = 0.021 at K=8; sky-ellipse area 5.76e−9 sr; σ(log₁₀ M_c) = 0.0169 dex.

### THE RANKING — achievable arm (σ_d → max(0.1 pc, 13 μas · d²))

| # | pulsar | cat now | d (pc) | σ_d now | σ_d new | w | **(b) Δf_marg K=8** | **(a) pin** | **(c) ΔA/A** | **(d) Δσ_Mc/σ** | S/N |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | J1216−6410 | timing PX | 939.7 | 185.02 | 11.48 | 2.18 | **+1.38e−2** | 0.00 | +1.2e−4 | +2.6e−2 | 4.5 |
| 2 | J0614−3329 | timing PX | 570.2 | 85.38 | 4.23 | 1.63 | +6.96e−3 | 0.04 | +3.0e−4 | +2.2e−3 | 4.9 |
| 3 | J1630+3734 | DM default | 89.0 | 24.00 | **0.103** | 1.75 | +6.52e−3 | **0.71** | **+9.7e−3** | +2.4e−4 | 4.7 |
| 4 | J1730−2304 | VLBI | 500.0 | 25.00 | 3.25 | 1.61 | +6.03e−3 | 0.29 | +4.3e−4 | +3.1e−3 | 3.5 |
| 5 | J1017−7156 | timing PX | 256.4 | 78.90 | **0.855** | 1.42 | +5.69e−3 | 0.33 | +1.5e−3 | +1.1e−3 | 5.6 |
| 6 | J2124−3358 | timing PX | 322.6 | 57.23 | 1.353 | 1.45 | +5.35e−3 | 0.17 | +8.8e−4 | +1.9e−3 | 4.8 |
| 7 | J0636−3044 | timing PX | 639.9 | 444.00 | 5.32 | 1.24 | +2.43e−3 | 0.00 | +2.3e−4 | +3.6e−3 | 2.3 |
| 8 | J1737−0811 | DM model | 206.0 | 41.20 | **0.552** | 0.64 | +1.66e−3 | 0.04 | +2.2e−3 | +6.4e−5 | 6.9 |
| 9 | J0125−2327 | timing PX | 1348.4 | 91.16 | 23.64 | 2.42 | +1.63e−3 | 0.00 | +4.9e−5 | +4.8e−2 | 1.5 |
| 10 | J2322−2650 | timing PX | 737.1 | 117.78 | 7.06 | 1.28 | +1.52e−3 | 0.00 | +1.8e−4 | +5.5e−3 | 1.8 |

**Is this a ranking, or is it noise?** Split-half over the 24 realisations:
Spearman(A, B) = **+0.400** across all 115 candidates, but top-10 overlap with the full
ranking is **9/10 in each half** and the top-10 median value/standard-error is **4.6**.
Only **11 of 115** candidates have S/N > 3. Read that as: **the top ~6–8 entries are real
and reproducible; ranks 9–10 are marginal (S/N 1.5–1.8); everything below the top ten is
noise** — because for most of the array the achievable improvement is nil (σ_π d² already
exceeds today's prior).

### One-line science case for the top five

1. **J1216−6410** — *the best-timed pulsar in the array with a genuinely bad prior.*
   w rank 3, σ_d 185 → 11.5 pc. Buys **soft-channel depth only**: the single largest
   Δf_marg at K = 8 (+1.4e−2, 27 % of what all ten anchors together deliver). It never
   becomes certifiable (pin 0.00) and it does nothing for localisation. Target it if you
   want the marginalised pulsar-term channel deepened, not if you want a pinned distance.
2. **J0614−3329** — *the same trade, one notch down.* σ_d 85 → 4.2 pc on a w-rank-13
   pulsar. Buys soft-channel depth (+7.0e−3). Cheapest large prior improvement in the array
   that is not distance-limited.
3. **J1630+3734** — *the only genuinely sub-parsec target in the array besides J0437.*
   89 pc, prior 24 pc → **0.10 pc**, a factor 233. It is the one pulsar that buys the
   **certification** product (pin **0.71**) and it is the **top localisation target by 4×**
   (ΔA/A = 9.7e−3). It buys essentially nothing for chirp mass (too nearby, no lever).
   **Caveat: its 89 pc is a feather `pdist` with no parallax-class measurement and no
   cross-check in the repo's own distance file.**
4. **J1730−2304** — *upgrade an existing VLBI parallax.* Ding et al. (2023),
   PX = 2.000 ± 0.100 mas → 500 ± 25 pc; reachable 3.25 pc. Buys soft depth (+6.0e−3) plus
   a real certification share (pin 0.29). The cleanest target in the list in the sense that
   its distance is already measured, so nothing about its rank is provisional.
5. **J1017−7156** — *sub-parsec, and the highest-confidence entry in the table.*
   256 pc, prior 78.9 → **0.855 pc**. Buys certification (pin 0.33), localisation
   (+1.5e−3, second best) and soft depth (+5.7e−3), at the best value/error of the top
   six (S/N 5.6).

### Driver decomposition — why the ranking looks like this

Rank correlation of marginal value against each candidate driver, over all 115 candidates:

| driver | achievable arm | oracle arm |
|---|---|---|
| 1 / distance | **+0.619** | −0.105 |
| current fringe count K_fr | **−0.608** | −0.003 |
| prior-width gain | +0.259 | −0.010 |
| pulsar-term power share | +0.214 | **+0.989** |
| timing weight w | +0.209 | **+0.990** |
| geometry ⟨dL_min⟩ | +0.038 | −0.095 |

**This is the whole result in one table.** In the *oracle* world — where any prior can be
driven to 0.1 pc for free — value is timing weight and nothing else (ρ = +0.990), and how
bad the prior was to begin with is irrelevant (ρ = −0.010). In the *achievable* world, the
astrometric d² wall dominates everything (ρ = +0.619 on 1/d), and timing weight falls to a
tie-breaker (ρ = +0.209). **Geometry does not enter at all (ρ = +0.04)** — and it should
not, because with an isotropic source sky every pulsar draws its μ from the same
distribution. Geometry only differentiates pulsars *once a specific CW is on the table*, at
which point this ranking should be recomputed against that source.

The oracle split-half is Spearman **+0.990** with 115/115 at S/N > 3, so the "value = timing
weight" statement is far better determined than the achievable ranking — which is another way
of saying the achievable ranking's noise comes from how few pulsars can actually be improved.

### The products do not agree with each other

Spearman between products, achievable arm, 115 candidates:

| | f_marg K=1 | f_marg K=8 | f_marg K=16 | bar K=8 | sky area | chirp M_c |
|---|---|---|---|---|---|---|
| f_marg K=1 | 1.000 | 0.260 | 0.356 | −0.017 | 0.347 | −0.362 |
| f_marg K=8 | 0.260 | 1.000 | 0.584 | 0.135 | 0.597 | −0.228 |
| f_marg K=16 | 0.356 | 0.584 | 1.000 | 0.251 | 0.682 | −0.284 |
| bar K=8 | −0.017 | 0.135 | 0.251 | 1.000 | 0.506 | −0.223 |
| sky area | 0.347 | 0.597 | 0.682 | 0.506 | 1.000 | **−0.588** |
| chirp M_c | −0.362 | −0.228 | −0.284 | −0.223 | −0.588 | 1.000 |

**Localisation and chirp mass are anti-correlated at ρ = −0.59.** Chirp mass is measured
from the pulsar-term line width, whose lever δ_a ∝ f τ_a / t_c grows with distance, so it
wants *far* pulsars whose fractional prior is bad. Localisation and certification want
*near* pulsars whose absolute prior can be driven inside a fringe. **There is no single
target list that serves both**, and a VLBI programme has to choose which product it is
buying. In the oracle arm the sign flips to +0.910, confirming the anti-correlation is
produced by the d² wall and not by the physics.

---

## Deliverable 4 — diminishing returns

**Structural warning first, because it changes what "saturation" means.** f_marg is *exactly*
separable — it is Σ_a R_a² pw_a / Σ pw with R_a depending on pulsar a's prior alone — so its
greedy curve *is* the sorted cumulative sum and **no anchor ever makes another anchor less
useful**. Any flattening there is a statement about the value distribution being
heavy-tailed. σ(log M_c) has a separable Fisher but σ = 1/√(ΣI_a) is not additive, so it gets
a true greedy (the naive cumulative sum overstates it ~2×). Only the localisation Fisher
genuinely couples pulsars, through the Schur complement.

| product | anchor #1 | total from 10 | 50 % at | 80 % at | 90 % at |
|---|---|---|---|---|---|
| f_marg (K=8), separable | +1.38e−2 | +5.16e−2 (**+28 %** of baseline 0.183) | #3 | #6 | #7 |
| f_marg (K=16), separable | +9.74e−3 | +5.87e−2 (+25 % of 0.233) | #4 | #8 | #9 |
| σ(log₁₀ M_c), true greedy | −8.0 % | **−19.8 %** (0.0169 → 0.0136 dex) | #2 | #5 | #7 |
| sky-ellipse area, **coupled** | −0.99 % | **−1.7 %** | #1 | **#3** | #5 |

Greedy order differs per product, as the anti-correlation demands. Chirp mass picks
J0406+3039, J0557+1551, J1902−5105, … (far, badly-primed, well-timed). Localisation picks
J1630+3734, J0711−6830, J1737−0811, J1017−7156, … (near, sub-parsec-reachable).

**Against the anchor literature's 3–6 saturation claim:** the *coupled* product reproduces it
exactly — localisation is 80 % done at anchor **3**, with marginal gains falling
#2/#1 = 0.25, #5/#1 = 0.056, #10/#1 = 0.024. **But the separable products do not saturate
there at all** — the soft channel is still returning 44 % of its remaining value between
anchors 4 and 10, and its curve flattens only because the *value distribution* is
heavy-tailed (one pulsar holds 27 % of the ten-anchor total), not because the anchors
interfere. **So "3–6 anchors is enough" is true for localisation and false for the soft
channel; the two claims have been conflated because nobody separated the coupled product
from the separable ones.** The correct statement for a programme proposal is: *for
localisation, three anchors capture the benefit; for pulsar-term depth and chirp mass, keep
going, and the reason to stop is the astrometric d² wall (only 11 pulsars are improvable at
all), not diminishing returns.*

The absolute sizes deserve equal emphasis. Ten anchors buy **+28 %** on soft-channel depth,
**−20 %** on σ(log M_c) — and **−1.7 %** on sky area. **The localisation case is real but
small.**

---

## Deliverable 5 — is J0437−4715's anchor status already banked?

| variant | σ_d (pc) | f_marg K=1 | f_marg K=8 | f_cert K=8 | area (sr) | σ(M_c) |
|---|---|---|---|---|---|---|
| anchor, brief (Reardon-era) | 0.110 | 0.19946 | 0.18334 | 0.02147 | 5.761e−9 | 0.0169 |
| anchor, repo catalogue (Deller 2008) | 1.320 | 0.19819 | 0.18305 | 0.01421 | 6.218e−9 | 0.0169 |
| **no anchor** (20 % DM prior) | 31.270 | 0.19814 | 0.17618 | 0.00150 | 6.234e−9 | 0.0169 |

Value already banked by the existing anchor: f_marg K=8 **+4.06 %**, sky area **−7.58 %**,
σ(M_c) −0.04 % (nil — J0437 is too nearby to carry a chirp lever).

**The headline is f_cert.** Today's certified pulsar-term power fraction at K = 8 is 0.0215.
Remove J0437's anchor and it falls to **0.0015** — **J0437 alone is ~93 % of every certified
distance the array has**, and the Deller-2008 value already delivers two thirds of that
(0.0142). Scored from its own no-anchor state against every other pulsar's achievable value,
J0437 would rank **#2 of 116**. **The community's one existing sub-parsec MSP parallax was
very close to the single best shot available, and its value is banked, not latent.**

Note also that a *further* J0437 campaign is worthless: 13 μas gives σ_d = 0.318 pc, which is
**worse** than the 0.11 pc it already has.

---

## Deliverable 6 — controls: is the sub-parsec certification real?

A 0.11 pc prior puts a pulsar in a window narrower than one fringe. In that regime both the
signal contrast and the incoherent-comb null contrast inflate, because "the best rival
fringe" is barely represented in the window — the sub-fringe artefact family this project
has already been bitten by. Audited three ways:

| variant | σ_d | K | f_cert | f_cert_true | **f_cert SCRAMBLED** | bar (J0437) | stat (J0437) |
|---|---|---|---|---|---|---|---|
| anchor 0.11 pc | 0.110 | 8 | 0.02101 | 0.02101 | **0.00000** | 35.30 | 37.22 |
| anchor 0.11 pc | 0.110 | 16 | 0.03574 | 0.03318 | **0.01574** | 51.39 | 105.10 |
| Deller 2008 | 1.320 | 8 | 0.01850 | 0.01850 | **0.00000** | 11.16 | 12.35 |
| Deller 2008 | 1.320 | 16 | 0.03205 | 0.02949 | **0.01536** | 23.75 | 38.05 |
| no anchor | 31.270 | 8 | 0.00046 | 0.00046 | **0.00000** | 7.14 | 1.56 |
| no anchor | 31.270 | 16 | 0.00683 | 0.00058 | **0.01020** | 16.12 | 7.41 |

**At K = 8 the channel is clean**: f_cert_true = f_cert exactly (every certified pulsar
picked the fringe containing its true distance) and the source-level scrambled null certifies
**nothing** in all three variants. The narrow-window inflation is real — J0437's bar rises
from 7.1 to 35.3 as its prior tightens — but the refit floor absorbs it exactly as the banked
lesson requires, and the statistic outruns it (37.2 > 35.3) rather than the floor failing.

**At K = 16 the channel is NOT clean.** The scrambled null certifies at 34–53 % of the signal
rate, f_cert_true < f_cert (some certified pulsars pick the wrong fringe), and in the
no-anchor row the scrambled column (0.0102) *exceeds* the genuine one (0.00058) by 18×. This
is **not** a null-count artefact: it persists unchanged at 32, 128 and 384 nulls. Sixteen
mutually-incoherent accumulated combs produce chance commensurability that the
incoherent-comb floor does not bound.

**Therefore: the (a) certification column is reportable at seed K = 8 and must not be used
at K = 16.** Products (b), (c) and (d) use no threshold at all and are unaffected.

---

## Two things this run fixes in the machinery it inherited

1. **VERNIER-24's `separability` builds per-pulsar timing weights and then never applies
   them.** The localisation Fisher here whitens the data vector by w before differentiating,
   so (c) responds to real timing quality.
2. **VERNIER-24's `chirp_precision` assumes τ_a is known**, i.e. that the distances are
   exact. The correct per-pulsar information after marginalising the prior is
   I_a = ((5/3)δ_a)² / (σ_f,a² + (δ_a σ_a/L_a)²), which reduces to V24's expression as
   σ_a → 0. With today's real priors that term costs a factor **3.5** in σ(log M_c) relative
   to the oracle — a pulsar with a 20 % DM distance contributes essentially nothing to chirp
   mass no matter how well it is timed. Baseline σ(log₁₀ M_c) is 0.0169 dex here, not the
   0.0023 dex the distance-free formula gives.

The K = 16 seeding also uses **16 distinct frequencies** (offset within-bin for the second
pass over the ten bins), so it does not repeat VERNIER-25's Arm K defect where K = 12 and
K = 16 collapsed onto the same ten-seed cell.

---

## What this does not settle

1. **The GW sky is a toy and the array is real.** Every Δ is a response of VERNIER-24's
   audited Pareto census (γ_eff = 1.56, N_bin from the audited counts, 300 explicit sources
   plus a Gaussian sea) seen by the true 116-pulsar geometry. The ranking inherits every
   caveat in VERNIER-24/25's sky and adds none of its own.
2. **Detections are seeded, not searched.** The K combs are placed at array SNR 12 at
   distinct frequencies, as VERNIER-25's Arm K does. The wheel's dominant failure mode —
   6.5 blends per 7.7 recoveries — is absent by construction, so every comb here is genuine.
3. **Three of the top eight targets rest on unmeasured distances** (§D1 caveat). If
   J1630+3734 is at 1 kpc rather than 89 pc it leaves the reachable set entirely and the
   localisation ranking loses its head.
4. **The absolute sky areas are not to be quoted.** 5.8e−9 sr corresponds to σ_sky ~ 4e−5
   rad, which is the known toy behaviour of a pulsar-term Fisher at 1–5 kpc with an
   OMC_MIN = 5e−2 antenna regularisation: fτ ≈ 1500 cycles makes the pulsar-term sky
   derivative enormous. Only the *fractional* change ΔA/A is meaningful, and even that is a
   Gaussian-approximation quantity.
5. **σ_π = 13 μas is a demonstrated best, not a survey-average.** It was achieved on a
   268 pc pulsar with a bright in-beam calibrator. Applying it uniformly across the array is
   optimistic; the 30 μas row is the realistic programme number and it halves the reachable
   set (6 → 3 sub-parsec).
6. **The achievable ranking's tail is noise** — split-half Spearman +0.40, only 11 of 115
   candidates at S/N > 3. Ranks 1–8 are reproducible; ranks 9–10 are marginal; ranks 11+
   should not be read as an ordering at all.
7. **Geometry is averaged over an isotropic source sky and therefore drops out.** For a
   *specific* detected CW the fringe spacing dL = c/(f(1−cos μ)) varies by 40× across the
   array and would re-order this list. This is a pre-detection targeting list.
8. **One fixed sky-realisation ensemble, no intrinsic pulsar red noise, no solar-system
   ephemeris error, no monopole/dipole**, and the pulsar term is placed in the same
   frequency bin as the Earth term at r_P = 1.196 (VERNIER-8), which is the series-wide
   simplification.
9. **No cost model.** Marginal value per pulsar is reported; hours per pulsar are not. Two
   targets with equal value and unequal calibrator availability rank equally here.
10. **The (a) column is audited only at K = 8** (§D6).

---

## One-line verdict

**The ranking exists, it is short, and it is a distance cut: because a parallax campaign
buys σ_d = σ_π d², the demonstrated 13 μas best precision reaches sub-parsec only inside
277 pc and sub-λ/2 only inside 158 pc, which is 6 and 3 of the 116 array pulsars
respectively — so the VLBI team's "2–3 more sub-parsec distances" is a geometric fact
reproduced here from the catalogue rather than a funding limit; value ranks on 1/d at
Spearman +0.62 with timing weight demoted to a tie-breaker at +0.21, and only in the
counterfactual oracle world where any prior is free does value become timing weight and
nothing else (+0.99); the top targets are J1216−6410 and J0614−3329 for soft-channel depth
alone (largest Δf_marg, never certifiable), J1630+3734 and J1017−7156 for certification and
localisation (the only genuinely sub-parsec-reachable pulsars besides the existing anchor),
and J1730−2304 as the one entry whose distance is already measured so its rank is not
provisional; localisation and chirp mass are anti-correlated at ρ = −0.59 because chirp mass
wants far pulsars with a long pulsar-term lever and localisation wants near ones that can be
driven inside a fringe, so no single list serves both and a programme must declare which
product it is buying; the diminishing-returns curve splits the literature's 3–6 saturation
claim in half — the coupled localisation Fisher reproduces it exactly (80 % of the benefit
by anchor 3, marginal #10/#1 = 0.024) while the exactly-separable soft channel does not
saturate at all and flattens only because one pulsar holds 27 % of the ten-anchor total, so
the real stopping rule is the astrometric d² wall and not interference between anchors; ten
anchors buy +28 % on soft-channel depth, −20 % on σ(log M_c) and only −1.7 % on sky area;
J0437−4715's existing anchor is banked, not latent — removing it drops the certified
pulsar-term power fraction from 0.0215 to 0.0015, so that one pulsar is ~93 % of every
certified distance the array has and would rank #2 of 116 if it were still a candidate,
while a further campaign on it would make its prior three times WORSE; and none of this
buys ignition, which VERNIER-24/25 and the R referendum have already settled at 0.00 on
every sky at K ≤ 16 — this is a products deliverable, and it is honest about the fact that
three of its top eight targets are "nearby" only according to a distance nobody has yet
measured.**
