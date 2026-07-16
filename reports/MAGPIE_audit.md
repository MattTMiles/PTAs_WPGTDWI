# MAGPIE — read-only reconnaissance audit

*Agent MAGPIE, ACCRE, 2026-07-15. Read-only: no new draws, no likelihood
evaluations, no tracked-file edits. All numbers below are cheap re-cuts /
re-reads of banks already on disk. Mission: find results we HAVE but never
LOOKED AT.*

`git pull --ff-only` → already up to date (HEAD d87db93). STORY.md read as the
claim skeleton. Scope: the campaign trees (GEO, RING, SIREN, ATLAS, FORGE,
IGNITE, IGNITE2, CHORUS, SURFACE, ANCHOR, RECUT, KINDLE, EMBER) + EXPLAINER.
The 41 357 npz on disk are dominated by the legacy `lnLs_*` per-cell scan trees
(pre-campaign); those are inventoried as buckets, not line-by-line.

---

## 1. FIGURE SWEEP

247 PNG on disk. Campaign + EXPLAINER figures (the ones tied to STORY claims):

| figure | shows | discussed in | flag |
|---|---|---|---|
| ATLAS_results/atlas_M2_contour_kappa.png | κ(e,mc) chirp-mass-wall contour | ATLAS_etrack_map.md | read |
| ATLAS_results/atlas_M3_ignition.png | ignition margins g_mc/g_e/g_forb | ATLAS_etrack_map.md | read |
| CHORUS_results/ch_pairs.png | mixed-pop pairwise count lifts | CHORUS_mixed_pop.md | read |
| CHORUS_results/ch_surface.png | count surface over (n_ecc,e,tier) | CHORUS_mixed_pop.md | read |
| CHORUS_results/ch_trade.png | the n_ecc trade / inversion | CHORUS_mixed_pop.md | read |
| EXPLAINER_results/F1..F14 (16 figs) | the narrative panels | reports/EXPLAINER.md | all read |
| FORGE_results/forge_figures.png | B1 conditional-noise loop, DAMPED | FORGE_b1_loop.md | read |
| GEO_results/geo_ensemble.png | certified-count distribution + selection fn | GEO_geometry_ensemble.md | read |
| IGNITE_results/ignite_onset_map.png | onset map over (h,T,tier) | IGNITE_onset_map.md | read |
| IGNITE_results/ignite_join.png, ignite_loop_trajectories.png | join panel + hard-loop traj | IGNITE_onset_map.md | read |
| IGNITE2_results/ig2_join_panel.png, ig2_soft_trajectories.png | soft-loop join + traj | IGNITE2_softloop.md | read |
| KINDLE_gain_contour.png | soft-loop gain contour (gain retired) | KINDLE_gain_contour.md | read |
| SIREN_results/siren_payoff.png | payoff chain, lag diversity | SIREN_payoff_chain.md | read |
| SURFACE_results/surface_onset.png, surface_join.png | onset surface (59 cells) + join | SURFACE_onset.md | read |
| **SURFACE_results/surface_floors.png** | **the per-cell null-floor figure** | **NOWHERE** | **UNREAD** |
| RING_results/ring_q1_fgw{-8,-9}_*.png (4) | area90 / inside90 / bias / detfrac | RING_q1_modernized.md (named, not embedded) | read (set) |

**UNREAD campaign figure: exactly one — `surface_floors.png`** (duplicated in
`reports/`). The SURFACE report cites the *table* `surface_floors_recut.npz` but
never the rendered floor figure; nobody has looked at the floor picture itself.
Note also RING references `fgw-9_{B,C}` which do **not exist on disk** (referenced-
but-missing, the inverse problem).

**Legacy / pre-campaign figure buckets (174 PNG, out of STORY scope):**
`lnLs_GWAmp*` (66), `lnLs_master_tests*` (49), `CW_transition` (17, prong-2 —
partly cited in project_progress.md/MANIFEST.md), `MAIN_PROJECT_QUESTIONS/Q1`
(16), `lnL_distance_scans` (24), `CW_node_sampling` (2). These are auto-generated
per-cell corner/diagnostic plots from the old scan pipeline; treat as pipeline
products, not un-mined results.

---

## 2. BANK SWEEP — banked columns no report section quotes

| campaign / bank | banked but un-quoted | note |
|---|---|---|
| **GEO** geo_dlnl_bank.npz | `feather_dlnL/K/P` **vs** `script_dlnL/K/P` (two independent pipelines) | report keeps both columns but never states the discrepancy → **38 cells with \|Δ\|>1 nat, max Δ=4.8 nat** |
| GEO geo_summary.npz | `omc_min_all`, `dL_pc` (per-psr), `freq_strict` (per-psr strict selection) | selection fn analysed on `omc_min_loud` only |
| **IGNITE** ignite_bank.npz | `weather` (2070 per-realisation noise levels), full `qmax`/`mapk`/`n_true_grid` cubes | weather never correlated with cert / wrong-cert |
| **IGNITE2** ig2_purity.npz | `sig_u_mag,sig_dk,sig_R_det/all,sig_nref_det` (809 wrong-cert rows) + `null_*` (280) | per-signal geometry aggregated only to a/b/c pass — the *anatomy* of wrong-certs is un-mined |
| **SIREN** siren_geom.npz | `lever_N{1,2,3,5}_freqrank`, `lever_N3_{shortlag,longlag}` | lag-lever grid; report quotes lag qualitatively, not the N-scaling levers |
| SIREN siren_probe.npz | `cosmu,tau_s,L0,sig_em` (116 per-psr) | per-pulsar siren geometry never joined to GEO selection |
| **ATLAS** atlas_kappa_forb*.npz | `kappa_analytic_trunc` vs `kappa_meas`, `sig_improve`, `cov_scaled` | measured-vs-analytic κ divergence not tabulated (see §4) |
| ATLAS atlas_ignition_forb*.npz | `tau_kyr` axis (0.3/1/3 kyr), `R_rank3` vs `R_scalar`, `ratio` | secondary τ axis swept, only τ=1 shown |
| **SURFACE** surface_analysis.npz | `corr99,wrong99,corr_fALL,wrong_fALL,corr_lo/hi,wrong_lo/hi` (99% arm + CIs), `laws` (18×5 power-law fits), `percert` (per-cell first-certifier identity) | report quotes point estimates; the **first-certifier identity is banked but never aggregated** (see §3) |
| **ANCHOR** anchor_data_forensics.npz | `chi2_over_N,resid_rms_s,span_yr,n_toa` (116 per-psr) | used once to declare "simulated"; per-psr distribution never joined to certification |
| ANCHOR anchor_ladder.npz | `mw_p` (Mann-Whitney), `d_gumbel`, `d_emp`, `psr_counts` (per-rung offender identity JSON) | offender pulsar identities banked, not surfaced |
| **KINDLE** kindle_g1_diagnosis.npz | `trans_keys/trans_counts` (12-way transition taxonomy) | gain retired; transition taxonomy never tabulated |
| **SAMPLER** ugrid_impact.npz | `CENTRE_ vs PROFILE_ vs MARGINAL_` qmax/mapk/gap arms | u-marginalisation impact banked, 3 arms, no report figure |
| **EMBER** (no summary bank) | **452 raw g-sweep partials, zero aggregation** | see §3 / §5 rank 2 |

---

## 3. CROSS-CAMPAIGN JOINS NEVER MADE

All keyed on the shared 116-pulsar spine or the shared e-axis — the banks are
already join-compatible; no join has been plotted.

**J1 — ATLAS κ(e) × CHORUS population lift** *(mission's own hypothesis)*.
Does single-source self-clocking κ **predict** the mixed-population count turn-on?
Teaser (banked columns, mc=9, forb=−9):

| e | ATLAS κ_meas | κ_analytic_trunc | CHORUS claim |
|---|---|---|---|
| 0.30 | **1.11** | 1.9 | "one e≥0.3 member turns the count on (10× lever)" |
| 0.50 | 0.92 | 5.6 | lever active |
| 0.65 | 3.80 | 20.7 | |
| 0.80 | 601 | 170 | |
| 0.90 | 16 857 | 1 750 | |

→ **κ is flat (~1) until e≈0.65, yet CHORUS's count lever is claimed to engage at
e=0.3.** On this teaser the single-source self-clocking κ does **NOT** explain the
population lift at low e — they are different mechanisms. This is the highest-value
join in the audit and is CPU-only to do properly.

**J2 — SURFACE first-certifier identity × ANCHOR offender identity.** SURFACE
`percert` aggregated over cells (never done before): J0437-4715 (690) ≫
J1713+0747 (491) > J0711-6830 (443) > J1909-3744 (419). ANCHOR `psr_counts`:
J0437 108–139/150 across every rung. **The same handful of pulsars carry both the
onset AND the null floor** — the onset may be a 3–4-pulsar phenomenon, not an array
one. Never stated.

**J3 — ANCHOR per-psr χ²/N × GEO certification frequency.** Spearman **ρ=+0.42
(p=3×10⁻⁶)** — pulsars that certify most often are the high-excess-power ones;
`resid_rms` is *uncorrelated* (ρ=0.036, p=0.7). Certification tracks excess power,
not timing precision. *(Caveat: spine is simulated CW+CURN, so high χ²/N partly =
more injected signal — but the RMS null result is still informative and unshown.)*

**J4 — GEO per-pulsar selection function × SIREN per-pulsar lag levers.** Both are
116-vectors on the same names (`geo.freq`/`omc_min` ↔ `siren_probe.cosmu/tau_s/L0`).
Does the geometry that selects a pulsar for certification also give it siren value?
Un-joined.

**J5 — IGNITE `weather` × wrong-cert rate.** Per-realisation noise level is banked
(2070 values) alongside every cert flag; the noise→wrong-cert link is never cut.

---

## 4. ANOMALY SCAN

- **ATLAS measured-vs-analytic κ divergence.** Non-monotonic and large: e=0.5
  κ_meas=0.92 while κ_analytic=5.6 (analytic *over*-predicts); e=0.9 κ_meas=16 857
  vs 1 750 (measured 10× analytic). The analytic truncation model breaks at both
  ends. Banked, never compared.
- **SURFACE null-floor blow-up.** `emp_q95` spans 5.7 → **8 740 nat** across cells;
  offender vectors reach 9 280 nat at high T/struct-k. Three orders of magnitude —
  the blind-robust floor is astronomically high in the high-T, high-k corner. Is
  this physical or an estimator artefact outside its validity domain (cf. S6.5)?
- **ANCHOR against-class pulsar.** J1824-2452A: χ²/N = **3 243** (magnetar-class
  excess) yet certification freq only 0.10 — huge power, rarely certifies. The one
  clear against-class member.
- **GEO pipeline disagreement.** feather vs script dlnL: 38/4640 cells disagree by
  >1 nat (max 4.8). Small mean (−0.003) so not a global bias, but the tail cells
  could flip a count. Never reconciled.
- **CHORUS fN/fALL inversion — NOT new.** The detect-floor-above-blind-floor
  inversion is systemic (19–22 of 26 cells) but is **already documented** in
  CHORUS_mixed_pop.md §1/§3 ("fALL inverts under the comb"). No *undocumented*
  second inversion case was found — SURFACE and IGNITE show no fN>fALL cells.
- **IGNITE2 b-arm fails** (`b_pass=False`, n_wrong_onset=23) and **KINDLE gain is
  degenerate** (NaN 90/140) — both already diagnosed in their reports.

---

## 5. RANKED TOP-10 (science-value per cost)

1. **ATLAS κ(e) × CHORUS lift — does self-clocking predict the population lever?**
   The teaser says no (κ flat to e≈0.65; CHORUS lever at e=0.3), which would mean
   the mixed-population payoff is a *distinct* mechanism from single-source κ — a
   structural claim about *why* eccentric members help. Both banks on disk, e-axis
   shared. **Cost: CPU-only** (one join script). *Top pick.*

2. **EMBER — aggregate the 452 orphaned partials.** An entire campaign (eccentric
   self-siren off-truth ladder, e=0.5, m=1e05, g∈{−1..3}, arms fix/map/truth) is
   fully banked with **no summary/analysis npz and no verdict** — only a
   pre-registration (EMBER_offtruth_ladder.md). The answer is on disk, unread.
   **Cost: CPU aggregation** (a surface.py-style consolidator). *Biggest "have-but-
   never-looked".*

3. **SURFACE first-certifier ranking (J2).** `percert` → the onset is carried by
   ~4 pulsars (J0437/J1713/J0711/J1909). If the onset is a handful-of-pulsars
   phenomenon it reframes the whole "array onset" claim. **Cost: CPU (essentially
   done above);** worth a figure + a paragraph.

4. **ANCHOR χ²/N × GEO cert-freq (J3), ρ=+0.42.** "Certifiers are high-excess-power,
   not low-RMS." Clean, significant, unshown. **Cost: CPU (done);** needs the
   simulated-spine caveat stated. Fold J4/J5 in as companion panels.

5. **GEO feather-vs-script reconciliation.** 38 cells disagree >1 nat. Could be a
   benign re-score artefact or a counting bug that moved the selection function.
   **Cost: CPU diff, at most one targeted re-score** to see which pipeline is right.

6. **IGNITE2 wrong-cert anatomy.** 809 signal + 280 null rows carry `u_mag`, `dk`,
   `nref`, `R` — the geometry of *why* a wrong cert fires, aggregated only to
   pass/fail. Mining it could yield the missing physical predictor of wrong-certs.
   **Cost: CPU.**

7. **SIREN lag levers × GEO selection (J4).** Do the pulsars geometry selects for
   certification also carry siren lag? Ties the census (S3) to the payoff (S9).
   **Cost: CPU.**

8. **SAMPLER u-marginalisation impact.** CENTRE/PROFILE/MARGINAL arms banked with
   qmax/gap; the effect of marginalising the fringe phase u on certification is a
   methods result never plotted. **Cost: CPU.**

9. **IGNITE weather × wrong-cert (J5).** Cheap correlation from banked columns; if
   noise weather predicts wrong-certs it's a free calibration knob. **Cost: CPU.**

10. **surface_floors.png + SURFACE 99%-arm/CI columns.** Look at the one UNREAD
    campaign figure and quote the banked `corr99/wrong99` + confidence intervals
    the report currently omits. **Cost: trivial** (open the png, quote 6 columns).

*(ATLAS κ_meas-vs-analytic divergence and the SURFACE floor blow-up are logged in
§4 as anomalies to explain, not separate work items — they fall out of items 1 and
10 respectively.)*

---

### One-line summary for Matt — top 3

1. **ATLAS κ(e) and CHORUS's population lever have never been plotted against each
   other, and the banked teaser says they disagree** — κ is flat until e≈0.65 but
   CHORUS's count turns on at e=0.3, so single-source self-clocking may not explain
   the mixed-population payoff (CPU-only to settle).
2. **EMBER is a whole campaign sitting unread** — 452 banked partials, no summary,
   no verdict, only its pre-registration; it just needs a consolidator run.
3. **The onset is carried by ~4 pulsars** (J0437/J1713/J0711/J1909 dominate SURFACE
   first-certs *and* ANCHOR floors), and **certification tracks excess power, not
   timing RMS** (ANCHOR χ²/N × GEO cert-freq ρ=+0.42, p=3e-6) — both joins are free
   and unmade.

---
---

# FOLLOW-UP — JOINS J1 & J3 EXECUTED (2026-07-15)

*Second brief: execute J1 and J3 only, CPU, read-only, one figure each. EMBER (#2)
left untouched (EMBER-2 is live). Both joins run on the CORRECTED banks. Figures:
`MAGPIE_J1_kappa_vs_channels.png`, `MAGPIE_J3_excess_power.png` (repo root,
untracked). The pre-RECUT teaser conclusions in §3/summary above are SUPERSEDED by
what follows — read these two pages, not the teaser.*

## J1 — κ vs the population lever → **VERDICT: MECHANISM SPLIT (refined)**

**The population count switch is set by the active-harmonic *channel budget*, not by
single-source κ.** Figure: `MAGPIE_J1_kappa_vs_channels.png`.

**Setup, corrected.** Counts + grades from `reports/ch_analysis_recut.npz` (the
post-RECUT bank; the e=0.3 single-member switch is refuted — current thresholds are
**e=0.5 single / e=0.3 two-member**, reproduced exactly here). κ from
`ATLAS_results/atlas_kappa_forb{0,1,2}.npz`. **Crucial correction to my own teaser:**
κ is *not* one curve — it is strongly f_orb-dependent, and CHORUS's members do **not**
sit at f_orb=1e−9. Per `hpc_harbor/chorus/ch_README.md` + the signal banks, CHORUS
census members are at **log f_gw=−7.9 ⇒ log f_orb≈−8.2** (bank `n_active_per_ecc`=17
for e=0.5 matches the README's harmonic count). So the relevant ATLAS rows are −8.5
and −8.0 (bracketing −8.2), where κ has **already** left ~1 by e=0.3 (κ≈2.7) and is
≈13 by e=0.5. My teaser compared to the wrong (−9) row.

**Why it is still not self-clocking — the decisive equal-κ contrast.** κ is a per-source
quantity; adding a *second* e=0.3 member does not change either source's κ. So m1e03,
m2e03, m3e03 all have **identical κ=2.65** yet:

| cell | e | κ@f_orb−8.2 | n_active (channels) | recut count | grade |
|---|---|---|---|---|---|
| m1e03 | 0.3 | **2.65** | 23 | 0.70 | below / MARGINAL |
| m2e03 | 0.3 | **2.65** | 30 | 2.77 | **CONFIRMED** |
| m3e03 | 0.3 | **2.65** | 37 | 2.50 | **CONFIRMED** |

κ held fixed, channel budget doubled, grade flips. κ **cannot** be the controlling
variable. And the channel threshold is clean: **every ON cell has n_active ≥ 30; every
OFF/marginal cell has n_active ≤ 23** — no κ analog (κ=2.65 appears in both ON and OFF
cells). Point-biserial grade-ON: r(n_active)=+0.53 vs r(κ)=+0.26. The per-source lift
is also **broadly distributed** — adding an eccentric member gives >0.5 nat to a median
of **47 pulsars**, with median max/sum≈0.15 (no single member/pulsar dominates) — the
registration-multiplication signature, not a concentrated Mc-tightening.

**What the banks support.** The count lift is **fringe-system multiplication** — the
comb's active harmonics acting as independent registration channels; the switch fires
at a fixed channel budget (n_active≈27), reached either by one high-e member (e=0.5→17
harmonics) or two low-e members (2×e=0.3→~16). Single-source Mc-information self-clocking
(κ) *is* real and *does* depart from 1 at CHORUS's f_orb by e=0.3, so a naïve
"count-turns-on-where-κ-departs" reading would wrongly call this CONSISTENT — but the
equal-κ contrast breaks the tie against it. **Caveat:** κ and n_active are collinear
along the e-axis (count ~ both: ρ=0.73 vs 0.81); they separate *only* on the
member-count axis at fixed e, which is exactly where this test lives.

## J3 — excess-power correlation, confound-checked → **VERDICT: ARTIFACT (signal-mediated)**

**The ρ=+0.42 is real but circular on the simulated spine; it does not transfer to real
data. The transferable finding is the *negative* one.** Figure: `MAGPIE_J3_excess_power.png`.

Base (reproduced): χ²/N × GEO cert-freq **ρ=+0.417 (p=3e−6)**; RMS × cert **+0.036
(p=0.7)**.

**(a) Noise-only control — passes.** In the IGNITE `no_cw` nulls (270 realisations, no
injected CW), per-pulsar **false**-cert frequency does **not** track spine χ²/N:
**ρ=−0.03 (p=0.8)**. So the effect is *not* a timing-structure artifact — high-χ²/N
pulsars have no intrinsic certification propensity absent signal.

**(b) Partial out precision & trials.** Controlling any one of σ_TOA / n_TOA / K leaves
it up (+0.36 / +0.25 / +0.36), but **jointly** controlling all three collapses it to
**+0.145 (p=0.12, n.s.)** at n=116.

**(c) Cross-checks.** Holds weaker on cert-freq_strict (+0.24, p=0.011) and on the
fALL/floor side — ANCHOR offender-count (+0.21, p=0.021).

**Name the path.** χ²/N = residual power ⁄ (N σ_TOA²). On the `b20_cw_curn_r0` spine the
residuals *contain the injected CW+CURN*, so χ²/N is elevated for exactly the pulsars
whose geometry/distance couples to the injected source — the same coupling that produces
the certifying dlnL. "Excess-power pulsars certify" ≈ "injected-signal-coupled pulsars
certify": expected, and it would **not** reproduce on a real array where χ²/N is not
dominated by one injected source. The noise-only pass (a) confirms it's *signal*-mediated,
not noise-structure — but signal here **is** the injection.

**The keeper (transferable, non-circular): certification is decoupled from timing
precision.** RMS×cert is null (+0.04) and joint precision/trials control removes
significance — so on a *real* array the loudest-coupled, not the best-timed, pulsars will
carry the first certifications. That is the part worth quoting; "excess power certifies"
should be retired as a target-selection claim on simulated data.

*(Correction logged: the §-summary line 3 phrasing "certification tracks excess power" is
demoted by J3 to signal-mediated/circular; the RMS-decoupling half stands.)*

**STOP.** Both joins complete; no further joins without a fresh brief.
