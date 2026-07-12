# IGNITE — the certification ONSET MAP + the above-onset Hogg loop

**Agent:** IGNITE · ACCRE · tag `criterion-v1` (`git describe --tags` → criterion-v1 @ db6ff05;
discovery `136b270f`, ee `f73b8e0`) · **Date:** 2026-07-12 · **PURE EXECUTION** (no commits,
no tracked-file edits).

**Scratch paths (host):** code `hpc_harbor/ignite/` (`ignite.py`, `ignite_consolidate.py`,
`ignite_analysis.py`, `ignite_figures.py`, 4 sbatch), results `IGNITE_results/` (2 110 lean
keyed npz + 3 figures ≈ 12 MB content; ~150 MB on-disk at panfs block granularity), logs
`hpc_harbor/logs/ignite_*`.

**Staged to `reports/` (sanctioned exception; Matt commits):** this report, the three
figures, `ignite_analysis.npz` (per-cell floors/onset/margins/tol tables),
`ignite_bank.npz` (the empirical basis: all 2 070 stage-0/1 realisations × 116 pulsars with
raw `dlnL_det/lnK/qmax/ptrue/on_true/mapk/n_true_grid` + cell metadata — refits every floor
in this report), the 35 `ig_loop*.npz` Stage-2 trajectories (raw finals + per-iteration
columns), and `ig_warm_T{20,30}.npz` (extension diagnostics). The full per-realisation tree
stays on scratch at the paths below; every staged number re-derives from `ignite_bank.npz`
alone. Original scratch paths: logs `hpc_harbor/logs/ignite_*`, container `/home/milesmt/soft/harbor/el9.sif`,
jax cache `/home/milesmt/.cache/jax_stagec_cache`. Lane
`-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1`.
Every summary npz carries the raw statistics — `dlnL_det`, `lnK`, `qmax`, `ptrue`, `on_true`,
`mapk`, `n_true_grid` — per the amended lean-npz convention. No verdict is banked without its
statistic.

---

## 0. THE ANSWER IN ONE PARAGRAPH

**The criterion-v1 floor is not a constant of the pipeline — it is loudness-relative.** The
per-cell refitted detection floor rises as `dlnL_floor ∝ h^{1.5–2.0}` (measured in every
baseline and tier, both null families), because the E-step's fringe statistic under a louder
*model* has proportionally larger null fluctuations — even on data containing **no CW at all**.
Making the source louder therefore raises the bar almost as fast as it raises the signal, and
the certified count is **non-monotone in h**. What does ignite the count is the **baseline**:
under the counterpart-matched (pure-noise) floor calibration, honest Arm-B certification first
exceeds 1 correct certification/realisation at **h\* = −12.75 for (T = 30 yr, literature
priors)** and **h\* = −13.25 for (T = 30 yr, VLBI tier)** — and **nowhere at T = 15 or 20 yr**
in the swept box (h up to −12.5). No cell reaches 3. Above census loudness the criterion-v1
purity property **collapses** (up to 23 wrong certifications per 50 realisations at the
onset cell): fringe correctness stops being automatic exactly where the count turns on. Under
the wrong-counterpart-robust (all-null) calibration the map never ignites anywhere — the
scrambled-source null's noise-lock grows ∝ h² and closes the window it was meant to guard.
The tolerance grid closes the K1 hole in the direction nobody expected: the pure-noise floor
is **small and flat-to-mildly-rising in registration tolerance** (0 → 4.4 nat at 5×10⁻⁴
scaled); it is the **true-positive channel that dies of mis-registration** (0.14 → 0.00
certs/real by tol = 5), not the null that inflates. And **Hogg's iterated phase-up does not
compound genuine certifications above onset — it compounds false ones**: at the three
above-onset cells the inherited source-fit channel amplifies seed impurity into a
wrong-certification cascade (up to 3 → 116 certified, 356/359 on false fringes at the
loudest cell), and the scrambled source run *through the loop* — the demonstration's own
null, never previously run — **detects (2/5)**. Both pre-registered STOPs fired; under the
wrong-counterpart-robust calibration the corollary is immediate (no seeds exist, loop INERT
by arithmetic). The campaign stopped there, per pre-registration.

---

## 1. STAGE 0 g1 — VALUE GATE (bit-identical, then five more gates)

`CW_transition/trackB_criterion.py` run against the banked `reports/` flats (code reused, not
reimplemented): binding invariants PASS; floor fit reproduces `floor_min = 9.0094 → 9.01` from
exactly the 9 banked offender cells (nullN J1713 9.009 binding); **THE FINAL TABLE reproduces
bit-identically: GEO 0.275 / Arm A 0.067 / Arm B 0.000 / all-nulls 0.000, wrong-certs 0.**

The IGNITE driver then passed its own value gates (job 12505333, ALL PASS, 464 s):

| gate | check | result |
|---|---|---|
| IG0 | criterion floor refit on banked flats | 9.0094 / 9 offenders / 27 nulls ✓ |
| IG1 | banked `flat_b10_g0_B_n900001` reproduced through the IGNITE path (h = −13.25, lit, T = 15) | `dlnL, lnK, qmax, ptrue, on_true, mapk, n_true` all **max\|diff\| = 0.0** |
| IG2 | banked `flat_wrongpos_g0_B_n920000` at tol = 0.5 (unit continuity: banked `tol_scale=5` ≡ 0.5×10⁻⁴-units) | **max\|diff\| = 0.0** |
| IG3 | loudness override live + monotone (same seeds, h = −12.5) | criterion-v1 detections 0 → **16** |
| IG4 | VLBI tier: union-18 `K_counted` | Σ 88 454 → **470**; runs finite |
| IG5 | T-extension path strict no-op at dT = 0 | `lnL(truth\|zero-noise)` = 405413.512739 **bit-equal**, = banked |

---

## 2. CONVENTIONS INTRODUCED (stated once, used everywhere)

- **Loudness axis.** The 3 loud sources' `log10_h` is overridden at runtime in θ, in **both**
  injection and recovery (the conditional pipeline knows the source to certification grade).
  Faint-13 stay at −14.25. dL/EV are strain-independent.
- **Prior tiers (binary, per RING).** `lit` = canonical literature priors. `vlbi` =
  σ_d → min(lit, 1 pc) on the **GEO union-18** (verified against `geo_summary.npz` at build
  time; improve-never-degrade, so J0437's 0.25 pc stands). The tier enters the FringeTables
  prior column, `K_counted`, **and the Arm-B truth draw**; the EV window stays the banked
  max-prior window in both tiers.
- **Baseline extension.** The fiducial array (true `getspan` = 22.15 yr; "the 15-yr array"
  throughout the campaign) **is** the T = 15 cell, unchanged. T = 20/30 append
  dT = T−15 yr of TOAs past each pulsar's own last TOA at its own median cadence (median
  frequency, modal backend flag, toaerr = the 1 µs load convention). Timing-design columns are
  extrapolated by per-column least squares on a smooth basis (polynomials ×
  annual/semiannual harmonics × 1/f² DM terms); columns fitting below R² = 0.99 — exactly the
  binary-orbital class (SINI/PB/A1/TASC/EPS/M2; day-period sinusoids far above the CW band;
  214 of 1 141 columns) — are zero-extended, i.e. their timing-fit leverage stays on the
  observed span. RN/GWB Fourier component counts scale with span (30→37→50, 14→17→23) so the
  noise-model f_max is span-invariant. `fdot` leverage ∼ T^{5/2} rides free. Injection
  parameters are drawn from the **fiducial** enterprise PTA (seed 3000 population), so the
  sources are identical across T by construction. T = 20: +14 275 TOAs (span 27.17 yr);
  T = 30: +42 835 (span 37.14 yr).
- **Tolerance units.** Registration offset = tol × 10⁻⁴ scaled (the spec-L286 certification
  tolerance), applied through FORGE's B1.4 `offset` mode; the banked `wrongpos` datum sits at
  tol = 0.5 in these units, reproduced bit-identically (IG2).
- **Per-cell nulls, two calibration families.** Every sweep cell banks ≥30 fresh nulls
  (10 nullN pure-noise/no-CW + 10 nullA all-16-scrambled + 10 nullL loud-scrambled, FORGE
  recipes) and refits its own floor — never inherited. The stage-1 decomposition (§4) forces
  the families apart: **fN** (nullN only — counterpart-matched: the conditional premise is a
  known source, the false alarm is noise mimicking fringe-breaking under the correct model)
  and **fALL** (all three — wrong-counterpart-robust, criterion-v1's original fit family).
- **Seeds** disjoint from every banked range (signal 1.0xx M, nulls 2.0xx M, tol 3.0xx M,
  loops 4.x M; dist_seed = noise_seed + 10⁷). Arm B throughout; Arm A retired. Truth
  (distances) redrawn per realisation.

---

## 3. STAGE 0 g2+g3 — NULL EXPANSION AND THE TOLERANCE GRID (the K1 hole)

**g2 — the floor's sampling scatter, measured.** 750 fresh nulls were banked (30 per sweep
cell + 90 on the tol grid). At the fiducial cell the fresh 30-null refit gives
**fALL = 8.48 nat** (banked 27-null value: 9.01) and **fN = 0.00** (no pure-noise cell even
passes layers 1+3 in 10 fresh draws — the banked 9.01 was set by a nullN J1713 fluctuation).
Because the scramble-family nulls are tolerance-independent, the four tol-grid refits are four
independent 30-null redraws of the same fALL statistic: **{8.48, 14.03, 8.09, 4.37} nat**.
The max-order-statistic the criterion floor rests on scatters by **±5 nat at fixed physics**
with 30-null samples. The criterion-v1 signal-side margin (0.29 nat) is **far inside this
sampling noise** — every onset-map margin below is quoted against its own cell's refit, and
they are all thin (0.01–2.0 nat).

**g3 — dlnL_floor(tol), the deliverable.** At the fiducial cell (h = −13.25, T = 15, lit):

| tol (×10⁻⁴ scaled) | fN (nullN refit) | fALL refit | true two-layer certs/real | max true dlnL | own-floor kills |
|---|---|---|---|---|---|
| 0 | 0.00 | 8.48 | 0.14 | 8.02 | 0/7 |
| 1 | 0.00 | 14.03 | 0.05 | 3.16 | 0/1 |
| 2 | 2.06 | 8.09 | 0.10 | 4.37 | 0/2 |
| 5 | 4.37 | 4.37 | 0.00 | — | 0/0 |

**Verdict: the counterpart-matched floor is FLAT-to-mildly-RISING in tol and small (≤4.4 nat);
the fALL spread is sampling noise, not tol dependence.** The K1 caveat inverts: what happens at
realistic registration error is not that the null inflates — it is that **the true-positive
channel dies of mis-registration** (true certs 0.14 → 0.00/real by tol = 5, exactly B1.4's
collapse), and no per-tol refit floor kills a single surviving true positive. The banked
9.01-kills-wrongpos-J0437 pathology was an artifact of applying a tol = 0 floor to a tol = 5
realisation; calibrated at its own tolerance the floor is smaller than every survivor.
(Fig. `IGNITE_results/ignite_onset_map.png`, panels e–f.)

**Resume discipline (the array's license).** Stage-0 task 12505409_0 was deliberately
`scancel`'d mid-run at 6 m 07 s with 29 npz banked; the resubmission logged
`already banked: 30; to run: 120` and completed — skip-on-exist resume proven on the
production SLURM path before any wide launch.

---

## 4. STAGE 1 — THE ONSET MAP

24 cells (h ∈ {−13.25, −13.0, −12.75, −12.5} × T ∈ {15, 20, 30} yr × {lit, vlbi}); per cell
50 Arm-B signal realisations (5 sky draws × 10 noise weathers, truth redrawn per realisation)
+ 30 fresh nulls; per-cell refit floors; criterion-v1 scoring with the cell's own floor.
1 920 realisations, 6 array tasks, 10–12 min each. **THE figure:
`IGNITE_results/ignite_onset_map.png`.**

### 4.1 The structural finding first — the floor scales with the model loudness

    dlnL_floor(fN)   ∝ h^{1.5–1.7}      (pure noise, NO CW in the data)
    dlnL_floor(fALL) ∝ h^{1.7–2.0}      (scrambled source meeting real CW data)

measured per (T, tier) across the full grid (fig. panel c). Mechanism: the E-step evaluates a
model whose pulsar-term amplitude ∝ h against the data; on pure noise the per-fringe
likelihood fluctuations grow with the model amplitude (matched-filter cross term), and with a
scrambled source meeting loud real data the noise-lock grows ∝ h². **The "absolute floor" of
criterion-v1 is loudness-relative**; 9.01 nat is its census-loudness value, nothing more.
Consequence: the certified count is **non-monotone in h** (e.g. T = 20 vlbi: 0.72/real at
−13.25 → 0.38/real at −12.5 — a 10× louder source *lowers* the honest count).

### 4.2 The onset table (correct certifications, counterpart-matched fN calibration)

| T | tier | corr certs/real over h = {−13.25, −13.0, −12.75, −12.5} | **h\* (>1)** | h\*\* (>3) |
|---|---|---|---|---|
| 15 | lit | 0.14, 0.06, 0.16, 0.48 | — | — |
| 15 | vlbi | 0.16, 0.22, 0.14, 0.52 | — | — |
| 20 | lit | 0.08, 0.72, 0.56, 0.54 | — | — |
| 20 | vlbi | 0.72, 0.72, 0.68, 0.38 | — | — |
| 30 | lit | 0.32, 0.96, **1.54**, 0.94 | **−12.75** | — |
| 30 | vlbi | **1.16**, 0.78, 0.98, 1.46 | **−13.25** | — |

**Onset is baseline-driven, not loudness-driven.** T^{5/2} fdot/coherence leverage beats the
h^{1.6} floor race; louder alone does not. T\* = 30 yr at both tiers; no cell reaches 3
correct certs/real anywhere in the box. Under the **fALL** calibration the map never ignites
(best cell 0.24/real): a criterion required to be robust to a wrong counterpart at these
loudnesses **closes the certification window entirely** — the choice of null family is now a
physics decision, not a technicality, and it gates every downstream claim.

**Null line and margins.** No null certifies under its own per-cell refit floor in any of the
24 cells (fit property, verified — the pre-registered STOP did not fire). The signal-side
margin (lowest surviving dlnL − floor) is **0.01–2.0 nat** across cells (loop cells: 0.03 nat
at (−13.25, 30, vlbi), 0.32 nat at (−12.75, 30, lit)) — the criterion-v1 thin-margin caveat
survives everywhere, against floors whose own sampling scatter is ±5 nat (§3).

### 4.3 Purity collapses above onset

Wrong certifications (criterion-v1 cert with `on_true = False`) per 50 realisations, fN
calibration: 0–2 at T = 15, up to 6 at T = 20, and **23/50 at the (−12.75, 30 yr, lit) onset
cell** (fig. panel d). criterion-v1's "the gate perfectly purifies what is left" was a
**census-loudness artifact**: above onset, wrong fringes acquire floor-beating gaps through
the same noise-lock that raises the floor. The wrong certs concentrate in the wide-prior,
data-driven pulsars (J1909, J1045, J1603, J1713 — never the anchor's clean channel: J0437
supplies 20/50 correct certs at the onset cell with 0 wrong). **Fringe correctness — the one
discriminator that survived real noise at census loudness — degrades exactly where the count
turns on.** Any above-onset certification claim needs a purity statement attached.

### 4.4 Per-pulsar: who arrives first, who second

At census loudness the first (and usually only) name is **J1909-3744** (wide prior, real
likelihood break). Above onset the workhorse flips to the **prior-pinned anchor J0437-4715**
(smallest K: 20/50 correct at both T = 30 onset cells), with J1909 second and the second-tier
names (J0711-6830, J1713+0747, J1045-4509, J1545-4550) arriving only at T = 30 or
h ≥ −12.75. It is **not** always J1909-first above onset — tiny K wins once the floor is the
binding constraint, because `max(lnK, floor)` is floor-dominated there and J0437's genuine
break rides a low trials bar. The VLBI tier's gain is concentrated exactly where RING said it
would be: it converts the union-18's trials mass (ΣK 88 454 → 470) into detections at fixed
dlnL — the whole (T = 30, vlbi, census-loudness) onset is this effect.

---

## 5. STAGE 2 — THE LOOP ABOVE ONSET. BOTH PRE-REGISTERED STOPS FIRED.

Job 12505592 (1 h 05 m; 18-param HVP cold compile 56 min, then 1.5–4 s/realisation).
Three above-onset cells × 10 Arm-B realisations + 5 scrambled-source realisations **through
the loop**, under the FULL criterion with each cell's own fN floor; source fit = sky fixed,
(f, mc, cos ι, h, φ₀, ψ) × 3 loud free (the inherited B1.3 machinery, extended per the brief);
soft-discipline W/W_grew columns every iteration. Trajectories:
`IGNITE_results/ignite_loop_trajectories.png`.

| cell (h, T, tier; fN floor) | grew past seed | seeds → finals | wrong at fixed point | scrambled loop |
|---|---|---|---|---|
| −13.25, 30, vlbi; 5.46 | 2/10 (3→7, 1→10) | Σ15 → Σ28 | **26 of 28 final certs wrong** | **2/5 DETECT (4 and 2 certs, all wrong)** |
| −13.00, 30, vlbi; 15.55 | 1/10 (1→2) | Σ5 → Σ6 | 5 of 6 wrong | — |
| −12.75, 30, lit; 30.89 | **6/10, incl. 2 full runaways (3→116, 3→116; also 2→79, 1→32, 5→8, 1→3)** | Σ20 → **Σ359** | **356 of 359 wrong** | — |

**VERDICT: none of the three pre-registered outcomes. The measured behaviour is a fourth
mode — a WRONG-CERTIFICATION CASCADE.** By raw count the loop "compounds" spectacularly at
the loudest cell (gain > 1 sustained ≥ 2 cycles in the runaway realisations, 3 → 116 in three
iterations); **essentially every certification it adds is on a false fringe** (356/359 at the
loudest cell, 26/28 even at the mission's expected (−13.25, 30 yr) cell). The genuine
(on-true) certified count never grows in any of the 30 realisations. The pre-registered STOP conditions
— *loop wrong-certs must stay 0* and *the scrambled loop must not manufacture detections* —
**both fired**, and per the mission the campaign stops here with the anatomy:

1. **The source-fit channel is the instability.** A certified pulsar is pinned at its MAP
   fringe *centre* (up to half a fringe off the true within-fringe offset); at onset loudness
   that mis-pin gives the (f, mc, extrinsics) gradient a large bias, and one damped Newton
   step moves the source materially — `src_mc_off` up to **1.6 dex** in a single step, 0.2 dex
   typically. The re-E-step at the moved source re-registers fringes: **17 of 20 final
   certifications at the quiet (no-growth) fixed points are wrong post-fit**, versus the
   ~8–15 % wrong-rate of the raw stage-1 seeds. The fit does not merely fail to compound —
   **it destroys the correct registration it was fed.**
2. **The runaway is the fALL mechanism, self-inflicted.** Once the source is mis-fit, the
   loop is a wrong source meeting loud real data — exactly the configuration whose noise-lock
   sets the fALL floor (~150 nat at this cell). After one poisoned step the E-step reports
   `n_detect = 116/116` above the 30.9-nat counterpart-matched floor, the fit's σ(log₁₀mc)
   collapses to ~10⁻⁵–0 dex (a huge Fisher at h = −12.75) so the registration gate R ≥ 1
   opens for everything, and the cascade certifies the whole array on false fringes within
   two more iterations. W = 115–116 from the first recorded iteration: the false-fringe mass
   saw it instantly. σ_mc "below the 0.003-dex self-clock bound" is reached trivially — and
   is meaningless: **tight local width + wrong global registration = confident nonsense**
   (B1.3's mechanism-1, weaponised above onset).
3. **The scrambled loop detects because the stage-1 family split says it must.** Its round-0
   detections (dlnL up to ~15 nat > fN = 5.46) are exactly the scramble-null events that
   separate fALL from fN. Under the counterpart-matched floor the loop has no defence against
   a wrong counterpart — that was the fALL floor's job.

**The pincer, closed.** Under the fALL (wrong-counterpart-robust) calibration, stage 1 finds
**no cell above onset anywhere in the box** — there is no seed set, and the loop is INERT by
arithmetic (no GPU time needed to know it). Under the fN (counterpart-matched) calibration
there IS an onset — and the loop above it is **unstable to its own fit channel and fails its
own null**. Either way, **Hogg's iterated phase-up does not compound genuine certifications
above onset with this machinery.** The fixed point never approaches the sky's GEO ceiling
(re-scored under criterion-v1: 0.275/draw at zero noise); it moves *away* from it, into
false-fringe territory. B1.3's DAMPED verdict, measured below onset, is superseded above
onset by something worse than damping.

**Caveat, stated honestly:** the instability is a property of the *inherited M-step*
(MAP-fringe-centre pinning + one damped Newton step + weak uniform prior), run exactly as B1.3
ran it (extended to extrinsics per the brief). A conservative M-step — within-fringe distance
refinement before the source step, trust-region on (f, mc), wrong-fix-aware soft weighting —
might tame the cascade; nothing here proves the *concept* of phase-up is dead, only that this
implementation of it, at these loudnesses, amplifies impurity instead of signal. That repair
is a design decision for cronus, not a fenced agent.

---

## 6. STAGE 3 — THE JOIN (onset × ATLAS self-clocking corner × SCOUT population clock)

Figure `IGNITE_results/ignite_join.png`; strain convention verified against SIREN's banked
DL table (Mc = 10⁹, f_gw = 2×10⁻⁸ Hz at 7.77 Mpc → log₁₀h = −13.25 exactly).

**The join table — the maximum distance at which the first source can sit and still be above
IGNITE's onset:**

| f_orb (Hz) | log₁₀Mc | ATLAS min-e (self-clock >20×) | D_L ≤ … at h\* = −12.75 (30 yr, lit) | D_L ≤ … at h\* = −13.25 (30 yr, VLBI) |
|---|---|---|---|---|
| 10⁻⁸ | 9.0 | 0.58 | **2.5 Mpc** | **7.8 Mpc** |
| 10⁻⁸ | 9.5 | 0.59 | **16.8 Mpc** | **53.1 Mpc** |
| 10⁻⁸ | 8.5 | 0.70 | 0.4 Mpc | 1.1 Mpc |
| 10⁻⁸·⁵ | 9.5 | 0.66 | 7.8 Mpc | 24.6 Mpc |
| 10⁻⁸·⁵ | 9.0 | 0.77 | 1.1 Mpc | 3.6 Mpc |
| 10⁻⁸·⁵ | 8.5 | 0.84 | 0.2 Mpc | 0.5 Mpc |

Read with ATLAS: the first source had to be *eccentric (e ≳ 0.6), massive (Mc ≳ 10⁹), at the
top of the band* to self-clock. IGNITE adds the fourth requirement: **near**. Only the
(Mc ≳ 10⁹·⁵, f_orb = 10⁻⁸) corner clears onset beyond the Virgo distance, and only in the
VLBI tier with a 30-yr baseline; the (10⁹, 10⁻⁸) reference cell must sit within ~8 Mpc even
then. The self-clocking corner and the certification-onset corner are the **same corner** —
loud, high-band, massive — plus a ≲50 Mpc distance cut.

**SCOUT's population clock.** SCOUT's census: no named credible SMBHB candidate clears even
h = −13.75 at its literature-favoured mass, the −13.75 operating point already sits above the
NG15 sky-averaged upper limit (−14.1), and the GWB-anchored expectation is
N̄_detectable ≲ 0.01–0.1 (current), O(0.1–1) (SKA era) — **at −13.75-class loudness**. IGNITE's
onset is 0.5–1.0 dex louder still (−13.25 VLBI / −12.75 lit, and only at 30 yr), i.e. an
Mc ≳ 10^9.5 binary at the top of the band within ~17–53 Mpc — an M87-class host inside the
Virgo volume. **When does nature supply a source above h\*? On expectation: not at any epoch
this campaign models — N_joint(h > h\*) ≪ N_joint(−13.75) ≲ 0.1.** The honest-certification
programme is therefore a *variance play on the loud-nearby tail* (SCOUT's jackpot caveat),
not an expectation-value plan; the expectation-value road remains the eccentric Earth-term
standard siren (ATLAS M4), which needs no certification at all.

---

## 7. WHAT IGNITE SETTLES, AND THE CAVEATS THAT TRAVEL

**Settled.**
1. **The criterion-v1 floor is loudness-relative** (∝ h^{1.5–2.0}, both null families, every
   baseline/tier). 9.01 nat is a census-loudness number. Any future loudness claim must refit
   per cell.
2. **The onset exists and is baseline-driven**: h\* = −12.75 (30 yr, lit), −13.25 (30 yr,
   VLBI); no onset at 15/20 yr in the box; no cell reaches 3 correct/real. Counterpart-matched
   calibration only — under wrong-counterpart-robust calibration there is **no onset
   anywhere**.
3. **Purity collapses above onset** (up to 23 wrong/50 real). Count and correctness must be
   quoted together from now on.
4. **The K1 tolerance hole is closed, inverted**: the pure-noise floor barely moves with
   registration tolerance (≤4.4 nat at 5×10⁻⁴); the true-positive channel is what dies.
5. **The floor's sampling scatter is ±5 nat per 30-null refit** — margins of 0.01–2 nat
   (including criterion-v1's 0.29) are inside the calibration noise. More nulls per cell is
   the single cheapest credibility purchase in the programme.
6. **Above onset the anchor leads**: J0437 (smallest K) overtakes J1909 once the floor
   dominates the trials bar.
7. **The above-onset loop question B1.3 left open is answered, and the answer is worse than
   DAMPED**: with the inherited machinery the loop is a wrong-certification cascade
   (§5) that fails its own scrambled null; genuine certifications never grow, and the fixed
   point moves *away* from the GEO ceiling. Both pre-registered STOPs fired. Any revival of
   the phase-up idea needs a redesigned M-step first (within-fringe refinement, trust
   region, wrong-fix-aware weighting) — flagged for cronus, not attempted here.

**Caveats.**
- Floors are max-statistics of 10 (fN) / 30 (fALL) fresh nulls per cell; individual floor
  values carry the ±5 nat scatter measured in §3. The *scaling law* and the *onset location*
  are robust against this (the h-exponent is fit across cells; onset cells exceed 1 by
  30–50 %); individual margins are not.
- The null-family split (fN vs fALL) is a **criterion-design decision this report surfaces
  but does not own**: fN presumes the counterpart is right (SCOUT/B1.4-style wrong-position
  events are then stage-2's scrambled-loop problem); fALL is criterion-v1's original family
  and closes the window. Cronus should adjudicate before any onset number is quoted outside.
- The T-extension is a stated convention (median-cadence TOAs, smooth-basis Mmat
  extrapolation, zero-extended binary columns, span-scaled GP components), gated as a strict
  no-op at dT = 0 — not a forecast of real future data (no new backends, no DM events, no
  solar-wind realism).
- Zero-noise GEO ceilings were not re-run at extended T; the onset map is entirely
  noisy-Arm-B.
- 5 sky draws × 10 weathers per cell; GEO showed the sky draw dominates yield variance, so
  the per-cell mean carries ~±0.2/real sky-sampling error (visible as the h-axis wiggle).

---

## 8. EXECUTION RECORD

| item | value |
|---|---|
| lane | `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1` |
| g1 (login, CPU, container) | criterion table bit-identical |
| gate | **12505333** ALL PASS, 464 s (IG0–IG5) |
| warm T=20 / T=30 | **12505408_20/_30**, 500/508 s, extension diagnostics banked (`ig_warm_T*.npz`) |
| stage-0 array + kill | **12505409_0** cancelled at 6 m 07 s (deliberate, 29 npz banked) |
| stage-0 resume | **12505435_0** COMPLETED 6 m 56 s — `already banked: 30; to run: 120` |
| stage-1 arrays | **12505485_0** (T15), **12505486_0-1** (T20), **12505487_0-1** (T30), all COMPLETED, 10–12 min each |
| stage-2 loop | **12505592** COMPLETED 1 h 05 m — 30 loop + 5 scrambled-loop realisations; 18-param HVP cold compile ≈56 min, then 1.5–4 s/realisation (first submit 12505591 died on an argparse leading-dash quirk; fixed with `--cells=`) |
| per-realisation wall (warm) | T=15: **0.5 s**; T=20/30 first-call 226/214 s (graph materialisation), warm ~1–2 s |
| disk | 2 110 lean npz + 3 figures in `IGNITE_results/` — ≈12 MB content (~150 MB on-disk, panfs block granularity); 4+ orders below the 3.1 TiB headroom |
| device log | every job logs GPU UUID + memory + foreign-process line as its first `[IGNITE]` line; squat lottery drew clean on every task |

**Nothing was committed. Nothing was pushed. No tracked file was edited.** `IGNITE_onset_map.md`,
`IGNITE_results/`, `hpc_harbor/ignite/` are untracked.
