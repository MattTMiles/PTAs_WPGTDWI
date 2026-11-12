# COMPASS — the isotropic control array (agent COMPASS, ACCRE)

**Session 2026-07-24 →. Venue: general H200 lane only (`batch_gpu`,
`taylor_group_acc`/`taylor_group_account_batch_gpu`, hgx03, 1-h backfill tasks under the
nodeupgrade drain); the reserved dgx share, the GLACIER pools and their spillover
untouched. Criterion-v2.2 floors throughout (Gumbel valid only at zero-fraction ≤ 20 %,
else empirical q95 ± bootstrap; zero-fraction a REQUIRED column; floors REFIT per unit —
per sky in arms I1/S1). Grades use the LOTTERY two-sided band rule. The scoring path is
GENERALISE's, reused verbatim (module hooks repointed), so gates g1/g3's licence carries.**

Scripts: `hpc_harbor/compass/compass.py` (+ `comp_gate.sbatch`, `comp_job.sbatch`).
Banks: `COMPASS_results/` (per-realisation lean npz with raw dlnL / lnK / qmax / mapk /
n_true_grid / ptrue columns + per-unit ledgers; `comp_isoarray_8100.npz`,
`comp_analysis.npz`). Real-side comparators: the banked A-SKY ledgers
(`GENERALISE_results/gen_ledger_AS_*.npz`) — **no real-array GPU respend**.

**MISSION.** Separate array-geometry effects from source-sky and physics effects, for
the first time: re-run the load-bearing boundary cells on a fake array that has the REAL
array's noise, distance and prior census but idealised (quasi-uniform) sky positions.
Whatever survives isotropisation is physics or source geometry; whatever changes was our
array's anisotropy all along.

---

## 0. THE BUILD — the fake array, and the audit that licenses it

116 pulsars at quasi-uniform positions: a Fibonacci-sphere (golden-angle) grid under a
seed-controlled random rotation, assignment by a RANK-PRESERVING map — the k-th best
real pulsar's noise bundle → the k-th position of a seed-controlled random enumeration
of the grid (`ISO_SEED_MAIN = 8100`). The **bundle travels intact**: TOA epochs/errors,
timing-model design matrix, the frozen per-pulsar red-noise assignment, the EM distance
(L0) and the lit/vlbi prior widths all stay attached to their pulsar; the sky position
(`pos/theta/phi`, on BOTH the enterprise and discovery objects) is the ONLY change,
applied by a runtime patch of `cw_helpers.load_pulsars` (the chorus-C7 pattern — no
tracked file touched) BEFORE any model build, so antenna patterns, pulsar-term
geometry, the HD-ORF GWB covariance and every fringe spacing follow consistently.

- **Quality rank** = descending TOA count (the harness flattens discovery-side TOA
  errors to 1 µs, so per-pulsar white-noise information ∝ n_toa; the real feather
  median TOA errors are banked beside it). Rank leader J1713+0747 (1 222 TOAs).
- **"Distances drawn from the real distance distribution" holds exactly** — the iso
  array's distance/prior census IS the real array's (bundles intact).
- **The audit** (`comp_isoarray_8100.npz`): nearest-neighbour separation
  min/median/max **16.50°/17.86°/18.25°** (real array: **0.95°/5.80°/27.97°**);
  position dipole |⟨n̂⟩| **0.0004** (real: **0.3956**). The real array's anisotropy —
  a 0.4-magnitude dipole and sub-degree pair clustering — is gone; the noise census
  (sorted n_toa / L0 / prior multisets) is bit-identical to the real array's.

**THE L_gwb VENUE-BANK HAZARD, CLOSED BEFORE IT FIRED.** The iso array's Phi_gwb has the
SAME shape (5336²) as the real array's, so GENERALISE's shape-keyed venue bank would
have been silently loaded with the WRONG (real-array) correlations. The compass loader
routes by iso state: iso OFF → the GENERALISE bank (`fp=f92c9e36b460d6f5`, required
bit-compatible by gate gA); iso ON → `comp_L_gwb_n5336_iso8100.npz`
(**fp=74ff91632ea8b407**, eigh at the pinned 8-thread convention on hgx03), keyed by the
assignment seed — arm S1's shuffled pairing changes the ORF and gets its own bank. The
run loop asserts the loaded provenance matches the unit's iso state before any
realisation is banked.

## 0.1 Gates (one small job, before any wide submit) — **ALL PASS** (job 12763272)

| gate | what it proves | result |
|---|---|---|
| gA (iso OFF) | the compass wrapper re-runs 5 nulls + 5 signal reals of the banked A-SKY survivor unit (AS_e03_h1275_k5_s4) with banked seeds and reproduces them | **PASS — strict bit-identity 10/10** (same host/threads/bank); certified+detected sets equal at the banked adopted floor; total certs 48 vs 48 |
| g0 (iso ON, the mission gate) | at zero signal the iso floor machinery reproduces Gumbel-validity and zero-fraction behaviour in class (one cell, N = 100 nulls: the survivor cell at GEO sky 4) | **PASS — zf 0.00, Gumbel valid, floor 145.48 ± 6.80 nat** — INSIDE the real cell's per-sky floor band (59.9–259.8), not even a level shift; + the uniformity audit above |

*(First gate submission (12763265) failed on a doubled path prefix in the real-ledger
template — fixed before any production realisation; its g0 nulls and the iso L_gwb bank
were banked and reused.)*

## 0.2 Budget header (before the wide submit)

Banked GENERALISE warm walls (0.51/0.66/0.83 s/real at n_ecc = 1/2/3, T = 30, this
venue): **I1 ≈ 1.9 + I2 ≈ 0.6 + S1 ≈ 0.5 ≈ 3.0 GPU-hr of compute** (A-SKY-class ×2,
inside the mission's 3–6 estimate) vs the **20 GPU-hr STOP bar → LAUNCH ALL ARMS, no
trims.** Worst-case walltime allocation 13 × 1-h backfill tasks = 13 GPU-hr, also under
the bar. Fan: I1 = 12763283 (8 tasks ×2-way packed), I2 = 12763284 (3), S1 = 12763285 (2).

SEEDS: I1 81/82.xM, I2 83/84.xM, S1 85/86.xM (+10M dist → 91–96.xM) — disjoint from
every banked range.

---

## 1. ARM 1 — THE ATTRIBUTION EXPERIMENT (the headline)

The four A-SKY boundary cells (incl. the survivor), 8 skies each (GEO draws 4–11,
**sky-paired with the real ensemble**, placements via the same `as_placement` map),
15 signal/sky, per-sky N = 100 floors. **Pre-registered reading, per cell:** the
sky-to-sky scatter ratio ρ = sd_sky(iso count)/sd_sky(real count). ρ < 0.5 →
**COLLAPSES** ("a uniform array de-lotteries the switch" — the lottery is array
anisotropy); ρ ≥ 0.75 → **PERSISTS** ("no array design escapes it" — intrinsic
source-geometry variance); between → MARGINAL (band). Primary = raw ratio (identical
protocol both sides, the 15-real sampling noise common-mode); a sampling-noise-
decomposed ratio is reported beside it, never substituted.

**COMPLETE** — 45/45 units (5 250 realisations; I1 32 units, I2 5, S1 8), all 13 fan
tasks COMPLETED in 7–11 min under backfill. The attribution table (count at adopted
floor, mean ± sky σ over the 8 paired skies; CONF = skies with count(floor+err) > 1;
real column = the banked A-SKY ledgers):

| cell | iso mean ± σ (median; CONF) | real mean ± σ (median; CONF) | **ρ = σ_iso/σ_real [16–84 %]** | floor swing iso / real | verdict |
|---|---|---|---|---|---|
| e0.3, h −12.75, 3+13 | 1.19 ± 0.85 (0.87; 2/8) | 1.33 ± 1.59 (0.90; 2/8) | **0.54 [0.28, 1.84]** | ×2.2 / ×2.6 | **MARGINAL (band)** |
| **e0.3, h −12.75, 5+11 (the survivor)** | **4.82 ± 2.69 (3.53; 8/8)** | 3.18 ± 3.02 (1.77; 6/8) | **0.89 [0.66, 1.18]** | ×2.2 / ×4.3 | **PERSISTS** |
| e0.3, h −12.50, 3+13 | 1.79 ± 1.53 (1.30; 5/8) | 1.43 ± 1.45 (0.90; 2/8) | **1.06 [0.40, 2.46]** | ×2.8 / ×3.4 | **PERSISTS** |
| e0.5, h −13.00, 3+13 (the old edge) | 1.61 ± 0.58 (1.73; 6/8) | 2.22 ± 2.51 (1.33; 4/8) | **0.23 [0.17, 0.40]** | ×1.6 / ×3.0 | **COLLAPSES** |

*(Sampling-noise-decomposed ratios 0.52 / 0.89 / 1.06 / 0.21 — indistinguishable from
the raw ones; the 15-real noise is common-mode as pre-registered. Wrong-cert rates
travel in the ledgers: iso 0.02–0.32/real vs real 0.09–0.23, same ordering — the
h = −12.5 cell is the wrong-cert-heavy one on both arrays. Zero-fractions 0.00 in 31/32
iso units; one old-edge sky at zf 0.21 adopted the empirical q95 — the v2.2 column
doing its work.)*

> ### VERDICT (pre-registered reading) — **THE LOTTERY IS NOT ONE THING.**
> **At the old single-member e = 0.5 edge, a uniform array DE-LOTTERIES the switch:
> ρ = 0.23 [0.17, 0.40], a 4× collapse of the sky-to-sky scatter — that cell's lottery
> was array anisotropy.** The paper sentence: *"a uniform array de-lotteries the
> e = 0.5 switch."*
> **At the e = 0.3 boundary cells the scatter is intrinsic: ρ = 0.89 [0.66, 1.18] at
> the survivor and 1.06 [0.40, 2.46] at h = −12.5 — no array design escapes it.** The
> paper sentence: *"the structure-assisted switch's sky band is source geometry, not
> array geometry."* (The faint 3+13 cell sits in the band, ρ = 0.54 [0.28, 1.84] —
> quotable only as MARGINAL at n = 8 skies.)
>
> **And the two HALVES of the criterion separate cleanly:** the FLOOR side of the
> lottery calms on the uniform array everywhere (per-cell floor sky-swing ×1.6–2.8 vs
> the real ×2.6–4.3 — the null tail's geometry sensitivity IS array anisotropy), while
> the COUNT side keeps its scatter wherever eccentric source geometry drives it.
>
> **The survivor gets STRONGER under isotropy: 8/8 skies CONFIRMED (was 6/8), mean
> 3.18 → 4.82, median 1.77 → 3.53.** The A-SKY external statement survives its
> geometry-control test with room to spare — its band is real, and it is the source's.

## 2. ARM 2 — THE CONSTANTS RE-MEASURED

(a) the channel-budget contrast (equal-κ, MAGPIE-J1 form: κ pinned at 2.65 while
channels run 23 → 30 → 37) at census loudness, iso vs the banked RECUT column
(0.70 / 2.77 / 2.50), plus m2/m3 at the survivor operating point; (b) the
union-ceiling analogue (same-protocol union certified set + per-pulsar frequencies,
iso vs real — **NOT** comparable to GEO's zero-noise union-18 or GALLERY's 30-readable
SNR_pterm census, which are different machineries); (c) one floor-law column (mean
per-sky floor vs h over the same e-mixture column as the real A-SKY readout — does
1.66 survive?).

### 2a. The channel-budget contrast — **the sharp ~30-channel flip is PARTLY OUR ARRAY**

Census loudness (h = −13.25, 3+13, lit, T = 30), κ pinned at 2.65 throughout,
`n_active` verified 23/30/37 from the banks:

| cell | channels | iso floor (v2.2) | iso count [pess] → grade | real re-cut |
|---|---|---|---|---|
| m1 e0.3 | 23 | 9.20 ± 0.89 (emp, zf 0.72) | 1.23 [0.97] → **MARGINAL** | 0.70 **below** |
| m2 e0.3 | 30 | 12.92 ± 0.64 (emp, zf 0.57) | 1.20 [1.00] → **MARGINAL** | 2.77 **CONFIRMED** |
| m3 e0.3 | 37 | 15.12 ± 1.36 (emp, zf 0.65) | 1.47 [1.23] → **CONFIRMED** | 2.50 **CONFIRMED** |

> **On the isotropic array the below→CONFIRMED FLIP between 23 and 30 channels
> disappears: the grades climb gradually (MARGINAL → MARGINAL → CONFIRMED), the
> single-member cell RISES (0.70 → 1.23) and the multi-member cells FALL (2.77 → 1.20,
> 2.50 → 1.47).** The channel budget remains monotonically helpful — the ordering
> survives — but the CRISPNESS of the ≥ 30 threshold (S7.6.4b) is partly a property of
> the real array's anisotropy, not of the physics alone. This lands exactly where Arm C
> already pointed ("necessary, not sufficient"); COMPASS adds: *the sharpness of the
> switch table's census-point flip does not survive geometry idealisation.*
> *(Caveat that travels: the real column is the banked dgx re-cut — cross-host class.
> g1 measured count-level host-robustness Δ = 0.000 at the neighbouring m1e05 cell, so
> the ±1.3–1.6-count moves here are far outside that systematic's reach.)*

At the survivor operating point (h = −12.75, 5+11): 23 ch (the I1 survivor) 4.82,
30 ch 4.07 [3.60], 37 ch 5.10 [4.50] — all CONFIRMED, Gumbel-valid floors ≈ 134–135
nat. **Above the structure-assisted onset the channel budget is saturated on the iso
array too; the 23-channel single member is already ON.**

### 2a.1 ADDENDUM (CPU harvest, 2026-07-25) — the carriers of the flip, decomposed

**The pre-posed question — is the real array's crisp 23 → 30 flip carried by a
concentrated readable set that iso disperses? — is answered NO, and INVERTED.** Per-
pulsar correct-certification counts over the 30 banked signal reals of each flip cell,
at each cell's adopted v2.2 floor (harvest reproduces every banked count/real exactly —
real 0.70/2.77/2.50, iso 1.23/1.20/1.47; readback-gated):

| cell (ch) | array | carriers | top-1 share | top-3 share | leaders (certs/30; mean exceedance over the bar, nat) |
|---|---|---|---|---|---|
| m1 (23) | real | 10 | 0.19 | **0.48** | J0711−6830 4 (2.0) · J1603−7202 3 (0.7) · J1713+0747 3 (2.4) |
| m1 (23) | iso | 13 | 0.32 | **0.57** | J1713+0747 12 (10.8) · J1909−3744 6 (4.7) · J0711−6830 3 (1.8) |
| m2 (30) | real | **26** | 0.13 | **0.33** | J1713+0747 11 (3.9) · J1545−4550 8 (3.0) · J1045−4509 8 (3.7) |
| m2 (30) | iso | **13** | 0.31 | **0.53** | J1713+0747 11 (18.0) · J1909−3744 4 (16.2) · J1603−7202 4 (8.8) |
| m3 (37) | real | 28 | 0.15 | 0.33 | J1713+0747 11 (5.7) · J0711−6830 8 (11.2) · J1909−3744 6 (16.8) |
| m3 (37) | iso | 13 | 0.23 | 0.48 | J1713+0747 10 (25.9) · J0711−6830 6 (9.1) · J1909−3744 5 (16.2) |

> **The real flip is a TAIL-RECRUITMENT event, and the tail is the anisotropy.**
> Crossing 23 → 30 channels on the real array the carrier pool BROADENS 10 → 26
> pulsars while the top-3 share FALLS 0.48 → 0.33: the flip is carried by a wide,
> shallow fringe of marginal carriers (per-cert margins ~3–4 nat) that the real
> array's clustered geometry places at favourable angles to the loud sources on lucky
> skies. **The uniform array never recruits that fringe: its carrier pool stays at 13
> in all three cells, MORE concentrated than the real array's (top-3 0.57/0.53/0.48),
> with the same few leaders carrying 3–5× deeper margins (J1713+0747 mean exceedance
> 10.8/18.0/25.9 nat vs 2.4/3.9/5.7 real).** Iso concentrates and deepens; real
> disperses and shallows — and the crisp census flip lives in the dispersal. This is
> the mechanism behind §2a's flattening, and the scope line the paper's channel-budget
> text needs: *the ≥ 30-channel switch's sharpness is a property of how an
> anisotropic array recruits marginal carriers, not of the channel arithmetic.*
> **A coherence worth one sentence:** on the uniform array the leaders ARE the
> best-timed pulsars (J1713/J1909/J0711 — ranks 0/1/3 of the build's n_toa metric),
> whereas the real m2 cell recruits J1545−4550, J1125−6014, … — S3.3.9's keeper
> ("certification is decoupled from timing precision") is an ANISOTROPIC-array
> statement; idealise the positions and certification re-couples to timing precision.
> *(Caveats: real side is the banked dgx re-cut (cross-host class, count-robust per
> g1); shares are over correct certs at the adopted floors; 30 reals/cell.)*

### 2b. The union-ceiling analogue — **the ceiling is a census property; the names are positional**

Same-protocol union over the 480 signal reals of each ensemble: **iso 75 vs real 70
pulsars** — isotropy buys ~7 % of union ceiling, nothing more. But the LEADERS reorder:
J1713+0747 tops both (115 vs 109 certs); **J1909−3744 (#3 on the real array, 72) falls
out of the iso top-8 entirely, and J0030+0451 (not in the real top-8) rises to #2
(81)**. Per-pulsar certification frequency is a POSITION property; the aggregate
ceiling is a NOISE-CENSUS property. (Scope: this union is criterion-v2.2 over these 4
cells × 8 skies — it must never be quoted against GEO's zero-noise union-18 or
GALLERY's 30-readable SNR_pterm census, which are different machineries on different
ensembles.)

### 2c. The floor-law column — **1.66 SURVIVES**

Mean per-sky floor over the same e-mixture column as the real A-SKY readout
(h = −13.0 (e05)/−12.75/−12.5 (e03), 3+13): iso 27.3 / 60.6 / 163.5 nat →
**floor ∝ h^1.55**; real 26.6 / 65.8 / 184.1 → **h^1.68**. Same absolute level, same
exponent class. **The floor-vs-loudness race is physics (the template-dominated null
tail), not array geometry — and it is one more axis on which the two arrays agree at
the ~10 % level while their sky-scatter behaviour differs by ×4.**

## 3. ARM 3 — REPORT-ONLY: the idealized-bound framing

**The isotropic onset at the tested cells is the geometry-best-case reference line for
any future-array sweep, and it says: uniformising TODAY'S array buys variance, not
onset.** With the real noise/distance/prior census held fixed and only the positions
idealised, the mean certified counts move 1.33 → 1.19, 3.18 → 4.82, 1.43 → 1.79,
2.22 → 1.61 — at most +52 % (the survivor), and DOWN at the old edge (the real array's
2.22 there was lottery wins, not capability) — while the mean floors track the real
ones to ~10 % (27.3/60.6/163.5 vs 26.6/65.8/184.1 nat) and the floor-law exponent is
unchanged. **No tested cell crosses a grade boundary in the mean; what a perfectly
uniform array buys at fixed census is the calm — floor sky-swings ×1.6–2.8 instead of
×2.6–4.3, and 8/8-sky confirmation where the real array manages 6/8.** A future-array
design that wants to MOVE the onset frontier must buy census (loudness reach, T,
distances, structure) — position optimisation alone converts a lottery into a
guarantee at the cells the census already reaches, and that is all it converts.

## 4. SCOPE LINES (travel with every number above)

- **The rank map is ONE choice of noise–geometry pairing — and it prices at ~0.1 sky-σ.**
  Arm S1 re-ran the survivor cell under a second, shuffled assignment (seed 8101, own
  L_gwb bank, fp banked): counts 5.19 ± 3.02 vs the main map's 4.82 ± 2.69
  (|Δmean| = 0.38 = **0.14 sky-σ**), 7/8 skies CONFIRMED (main: 8/8), floors 69–206 vs
  65–145 nat. Every Arm-1 conclusion is insensitive to the pairing choice at the
  resolution this campaign can see.
- The timing-model design matrix (incl. astrometry columns fit at the ORIGINAL
  position) travels with the noise bundle — the marginalisation structure is treated as
  a noise property, not a sky property. Second-order; stated, not hidden.
- Source skies are the banked GEO isotropic draws in BOTH ensembles — Arm 1 isolates
  the ARRAY's geometry at fixed source-sky protocol; it says nothing about anisotropic
  source populations.
- Iso-vs-real comparisons are same-venue, same-host, same-thread-convention
  (hgx03, 8 BLAS threads); the L_gwb provenance is banked per realisation.
- lit tier only; T = 30 only (the A-SKY pre-registration); fALL not measured
  (nullN-only floors, per the campaign spec).

## 5. GPU spend ledger + verdict block

**2.12 GPU-hr of walltime, sacct-verified, across all 15 COMPASS jobs of 2026-07-24**
(2 gate jobs + 13 fan tasks, all COMPLETED; the first gate's partial run was reused,
not wasted) vs the 20 GPU-hr STOP bar and the 3.0 GPU-hr header estimate — 2-way
packing and warm caches sold back the difference. All jobs on `batch_gpu`/hgx03 (H200)
via 1-h-class backfill under the nodeupgrade drain; the dgx share and every GLACIER
artefact untouched. On-disk: `COMPASS_results/` ≈ 5 250 realisation npz + 45 ledgers +
2 iso L_gwb banks (~230 MB each) — untracked, like every campaign results dir.

### VERDICT BLOCK (one line per arm)

- **BUILD/g0:** the isotropic control array exists and its floor machinery is in class
  (zf 0.00, Gumbel-valid, floor inside the real per-sky band); dipole 0.396 → 0.0004.
- **ARM 1:** the sky lottery separates — the e = 0.5 old edge's scatter was ARRAY
  anisotropy (ρ = 0.23, "a uniform array de-lotteries the switch"); the e = 0.3 cells'
  scatter is SOURCE geometry (ρ = 0.89/1.06, "no array design escapes it"); the floor
  half of the lottery is array, the count half at e0.3 is source; the survivor
  strengthens to 8/8 skies.
- **ARM 2:** the crisp ≥ 30-channel census flip is partly our array (grades flatten on
  iso, cells move both directions); the union ceiling is census (75 vs 70) while the
  leader names are positional (J1909 falls out, J0030 rises); floor ∝ h^1.55 vs 1.68 —
  the 1.66 law survives. **Addendum (§2a.1): the flip's mechanism is tail recruitment —
  the real array crosses 23 → 30 channels by broadening its carrier pool 10 → 26
  pulsars (top-3 share 0.48 → 0.33); the uniform array never recruits that fringe
  (13 carriers throughout, deeper margins) — and on it, certification re-couples to
  timing precision.**
- **ARM 3:** uniformising today's array buys calm, not onset — mean counts move ≤ 52 %
  (down at the old edge), floors track to ~10 %; position optimisation converts the
  lottery into a guarantee only where the census already reaches.
- **S1:** the rank-map pairing prices at 0.14 sky-σ — the build choice is not
  load-bearing.
