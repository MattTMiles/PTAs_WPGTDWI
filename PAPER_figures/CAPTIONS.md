# PAPER_figures — caption drafts, per-panel bank sources, intended sections

**Agent EASEL-5 · ACCRE · 2026-07-23.** Every number on every panel is read back off a banked
npz on this machine (CPU re-cuts only), except where a panel explicitly states otherwise on its
face; those exceptions are enumerated in `GAP_REPORT.md`. All counts follow criterion-v2.2
(Gumbel floor valid only at null zero-fraction ≤ 20%; empirical q95 + bootstrap above it;
certification on `q_max`, never `P_true`). Figures are drawn from the corrected (`_recut`)
banks; no superseded number is plotted. Style: `easel5_style.py` (Okabe–Ito palette, STIX
serif, PDF + 400-dpi PNG per figure).

---

## P1 — THE PROBLEM (`P1_the_problem`) · intended: Introduction / problem statement

**Caption draft.** The phase-connected CW likelihood is periodic in each pulsar's distance:
panel (a) shows a banked likelihood-vs-distance scan for J0437−4715, a comb of near-degenerate
fringes with the true distance (line) in no way preferred. (b) zooms to the electromagnetic
prior window (±1 pc, shaded): the measurement reduces to identifying *which* fringe the pulsar
occupies. (c) The census answer for the real array: the number of fringes inside each of the
116 pulsars' EM priors, median over 40 sky draws — the median pulsar has K ≈ 6,000 candidate
fringes, the best-anchored (J0437−4715) has 6, and the strict prior-alone anchor set (K ≤ 1) is
empty: **no free anchors exist**. The panel-(a/b) comb is a legacy single-source scan whose
teeth are coarser than the census's (the only banked lnL(D) curve in the repo); panel (c)
carries the quantitative claim.

**Sources.** (a,b) `lnLs_GWAmp_phase_connected/runD_3CW_test/J0437-4715_cw_p_dist/psrTerm/
runD_3CW_test_J0437-4715_cw_p_dist_psrTerm_0_0.npz`; prior `CW_transition/anchor_a0_priors.npz`.
(c) `reports/geo_dlnl_bank.npz` (`K_counted`, 40×116; median over draws). STORY S1.1, S3.1.

## P2 — THE WALLS (`P2_the_walls`) · intended: the blind and targeted no-gos

**Caption draft.** (a) The registration ladder: sky tolerance of all 348 pulsar–source pairs
(the offset that shifts a pulsar-term phase by 1 rad) against the measured blind-float band
(0.05–0.21 scaled). The gap is 27–112×, and **zero pairs sit inside it** — the wide-lane
cascade has no first rung. The fixed-integer pull-in basin (<10⁻⁴, dotted) is 4+ orders inside
any achievable float. (b) The same wall as a timescale, under both banked SNR scalings: the
blind float reaches the loosest rung only at T ≈ 11,000 yr (baseline scaling; 135 yr under the
optimistic T^{-3/2} bound), and the pull-in basin at 3.75 Myr (945 yr optimistic). (c) The
targeted referendum: even with the sky supplied exactly, posterior mass on the true fringe is
f = 0.085 (sky only) → 0.048 (+EM period) → 0.032 (+1% host distance; frozen 4-seed protocol,
adjudicated D-2; the 5-seed trail value 0.043 is the open marker) against the 0.95 bar — and
the inset shows the chirp-mass price: σ(log₁₀Mc) must improve by more than 20× (a bound, not a
value; even 20× reaches only f = 0.67). The scope line: panels (a,b) are pure geometry at zero
noise; panel (c) plots the canonical STORY values because their source npz exist only on cronus.

**Sources.** (a) `CW_transition/trackB_F2_ladder.npz` (`tol_sky`, 348); blind band + pull-in =
banked Track-B close-out numbers (spec F2/L2b/L2c). (b) `CW_transition/PENCIL_t_crossing.npz`
(`Tcross_loose_p05/p15`, `Tcross_l2c_p05/p15`). (c) STORY S4.2.8 (D-2 primary) + S4.2.10 —
`b1_step2_table.npz`/`b1_breakeven_curve.npz` **not in this checkout** (see GAP_REPORT §1).
STORY S4.1, S4.2.

## P3 — THE HONEST CRITERION (`P3_the_honest_criterion`) · intended: the criterion section

**Caption draft.** (a) Offender distributions at a representative cell (h = 10⁻¹³, T = 30 yr,
lit tier, census structure): 100 noise-only realisations (grey) against 30 signal realisations
(blue), re-scored by this figure's own gated scorer. The retired max-of-10 bar (13.5 nat)
yielded 1.03 certifications/realisation — a detection; the honest bar (empirical q95,
16.6 ± 1.6, N = 100) yields 0.60 and the detection is **retracted**. This cell's own Gumbel fit
(dotted, 19.5) is invalid: 27% of its nulls are silent. (b) Why the estimator has a validity
domain: the ratio of the Gumbel floor to the empirical q95 across all 108 SURFACE cells and the
48 ANCHOR-ladder cells, against the null zero-fraction. Below the 20% gate the ratio is ≈1;
above it the Gumbel is dragged toward the point mass at zero and errs PERMISSIVE, by up to 2.8×
at 93% zero-fraction — the criterion-v2.2 convention (empirical q95 above the gate) exists
because of this panel. (c) The realism ladder: per-cell floor shifts (Δq95 vs the R0 control)
under RN/GWB mis-specification (R1), an unmodelled DM GP (R2), and non-Gaussian tails (R3) are
consistent with zero (median −0.18 nat; 1 of 40 cells significant, ≈2 expected by chance) while
the realised noise rms rises ×1.32 — the floor's upper tail is template-dominated, so
mis-specifying the noise does not move it. (*The per-rung rms values are quoted from the
tracked ANCHOR report table; they are not npz-banked — GAP_REPORT §4.*)

**Sources.** (a) `SURFACE_results/sf_{nullN,sig}_h1300_T30_lit_k3_g*_n*.npz` (re-cut; scorer
gated to reproduce the banked count) + `reports/surface_analysis_recut.npz`
(`floor_adopted`, `fN`, `fN_zerofrac`, `corr`). (b) `reports/surface_analysis_recut.npz`
(`fN`, `fN_emp`, `fN_zerofrac`) + `reports/anchor_ladder.npz` (`gumbel`, `emp_q95`,
`zero_frac`). (c) `reports/anchor_ladder.npz` (`d_emp`, `mw_p`, rungs R1–R3); rms row from
`ANCHOR_realdata_null.md` §2 (prose). STORY S5.3, S6.5, S6.4.3.

## P4 — THE SWITCH (`P4_the_switch`) · intended: the eccentric switch-on

**Caption draft.** (a) Correct certifications per realisation at fixed census loudness as a
function of the eccentric members' eccentricity, for one, two, and three eccentric members in
both prior tiers, scored against per-mixture criterion-v2.2 floors; whiskers reach the count at
floor + bootstrap error (the CONFIRMED test). With one member the switch-on is at e = 0.5
(3.13 lit; e = 0.3 gives 0.70, below the bar); with two or more it is at e = 0.3 (2.77 lit /
1.77 vlbi, confirmed in both tiers). The all-circular population sits at 0.37. (b) The
mechanism is the total active-harmonic channel budget, not single-source self-clocking: at
fixed κ = 2.65 (the e = 0.3 members' clocking factor at their f_orb ≈ 10⁻⁸·²), growing the
channel budget 23 → 30 → 37 flips the grade — every ON cell has n_active ≥ 30 (line). (c)
LOTTERY's mixture forecast: P(switch-on) under the channel criterion over the (f_ecc, e_char)
population grid, with the conservative threshold rule (e ≥ 0.5 single / e ≥ 0.3 pair) as the
contour and the GPU-validated tier-2 cells overlaid: sub-rule cells DO certify, so the quoted
rule is a one-sided conservative subset of the channel criterion (and the channel proxy itself
is floor-gated, not a strict switch — LOTTERY tier 2).

**Sources.** (a) `reports/ch_analysis_recut.npz` (`surface_corr`, `surface_corr_lo`). (b)
`CHORUS_results/ch_analysis.npz` (`surface_nactive`; n_active convention S7.6.4c) ×
`ch_analysis_recut` counts; κ = 2.65 from `reports/MAGPIE_audit.md` J1 (ATLAS
`atlas_kappa_forb{1,2}.npz`). (c) `reports/lottery_tier1.npz` (`P_chan`, `P_thr`) +
`reports/lottery_tier2.npz` (`rows_json`, grades incl. the dgx-carried m1e03 'below' point —
the host-systematic band). STORY S7.6.4–S7.6.4c; LOTTERY tier1/2 + recut band.

## P5 — THE ONSET MAP (`P5_the_onset_map`) · intended: the onset surface

**Caption draft.** The complete corrected onset surface: correct certifications per realisation
over (h, T) for all six (structure × prior-tier) panels — 108 cells × 30 signal realisations,
each scored against its own criterion-v2.2 floor (N = 100 nulls, α = 0.05). The solid contour
is the onset boundary on the count at floor + error; the thin dashed contour is the bar at the
nominal floor, and the region between them is the floor-calibration error band. Hatched cells
are MARGINAL (clear the bar at the floor, die at floor + error; 3 cells — not onsets, D-6
strict-`>` convention). Dotted borders mark the 15 cells whose Gumbel failed the zero-fraction
gate and whose floors are the empirical q95 with bootstrap errors. N_onset = 59 of 108; the
census-structure onset sits at h* = −12.50 (T = 30, lit); T = 30 is optimal in 0 of 36 columns
(loud cells peak at 40 yr); promoting 3+13 → 5+11 moves the frontier to the faint edge of the
grid. Everything simulated end to end; the T = 50 row extrapolates the timing model under a
stated convention and is not a forecast.

**Sources.** `reports/surface_analysis_recut.npz` (`table` cols `corr`, `corr_lo`;
`verdicts`, `estimator`, `tiers`, `struct`; gates A/B = 0.000e+00). STORY S6.3.2–S6.3.5, S7.1.1a.

## P6 — THE TRANSITION (`P6_the_transition`) · intended: prong-2 anatomy

**Caption draft.** (a) Array-median marginal/conditional pulsar-distance information vs the
number of CW sources (band: 16–84%): 0.99 (N=1) → 0.87 (N≈100) → 0.53 (N≈400) → 0.08 (N=1000),
with the confusion knee at N ≈ 412 for the fiducial array. Inset: the knee scales as
0.40·N_psr^1.03 (linear through 200 pulsars), which calibrates Boyle & Pen's 2N/7 bound — the
2/5 ÷ 2/7 = 7/5 ratio is the price of the two source-sky parameters our metric supplies. (b)
The coherence law: distance information halves when the Earth–pulsar phase wander reaches
~1/SNR rad, in both banked functional forms; the two forms' 0.5-crossings differ by 0.30× and
the invariant is the location SNR²σ_φ² ≈ 1, beyond which σ_L → d_L σ_φ/2π, an SNR-independent
wander floor (absolute abscissa values are convention-dependent and never physical). (c) The
two transitions compound: the residual of R against the factorised product (commuting
reduction) is two-sided — a weak positive (de-confusion) lobe at N/N* < 1 and a dominant
negative lobe at N/N* ≳ 3, peak |residual| 0.202 — and the R = 0.5 contour bows toward the
origin, so a realistic background sits in the compounding region. The real likelihood never
reaches this knee at achievable N_CW; fringe identification binds first (S2.4.3) — this figure
bounds the regime, it does not describe the operating point.

**Sources.** (a) `CW_transition/prong2_results.npz` (`tot_n`, `tot_marg(+lo/hi)`, `tot_cond`,
`diag_T15_b3-20_knee`); inset `CW_transition/stagec_p2c_results.npz` (`A`, `b`, `n_psr`,
`knee_over_Nstar`). (b) `CW_transition/prong2_coherence.npz` (`x_grid`, `exp1_linear`,
`exp1_saturating`, `x_half_*`). (c) `CW_transition/prong2_Ntc_map.npz` (`resid_commuting`,
`R`, `x_grid`, `Y_grid`, `max_factor_dev_commuting`). STORY S2.1–S2.2.

## P7 — THE LOOP: SAFETY (`P7_the_loop_safety`) · intended: the loop section

**Caption draft.** The loop's final anatomy as the three-claim decomposition (S8.5.0; no
adjective crosses a campaign boundary). (a) Claim 1 — given certified seeds at/near truth the
soft (fringe-marginalised) loop HOLDS and self-cleans: all 70 banked IGNITE-2 + CHORUS loop
trajectories — 64 flat, 2 gain a genuine +1, 4 scrambled false alarms shed by the loop itself;
no trajectory cascades. Claim 2 (box) — from the honest Earth-term cold start (EMBER map arm,
N = 108 signal realisations) the loop is INERT: ΔN = 0.000 in all 9 cells, because the MAP
start lands ~a dex cold in log₁₀Mc; the EMBER truth arm shows the same machinery repairing
(ΔN up to +3) when already engaged. (b) Claim 3 — the danger mode: across the 112-realisation
scrambled census (every certification false by construction), 5 realisations manufacture false
certifications, and every one of them accepted M-step motion (5/5); motion under a wrong
counterpart, not start engagement, is the boundary. (c) The predictor comparison from the
banked contingency tables: motion sensitivity 1.00 (Fisher p = 0.0015) against engagement 0.60
(p = 0.33); the pre-registered joint (engaged AND mobile) boundary is REFUTED — it misses the 2
disengaged manufacturers. Scope: N = 112 with 5 events; the claim is the ordering, not a rate.

**Sources.** (a) `IGNITE2_results/ig_sloop*.npz` (40) + `CHORUS_results/ch_sloop*.npz` (30)
(`traj_n_cert_true`, `traj_wrong`); `EMBER_results/ember_analysis.npz` (`r_arm`, `r_dN`,
`r_scr`, `verdict`). (b,c) `EMBER_results/ember_predictors.npz` (`scr_*`, `cont_*`,
`fisher_p`, `logit_beta`, `motion_thr`). STORY S8.2, S8.5.0–S8.5.3.

## P8 — THE ARROWS (`P8_the_arrows`) · intended: the cascade arithmetic

**Caption draft.** (a) SPARK's certified-coherent detector on the 13-source faint reservoir:
per-state 2F values (points), signal medians (dotted), and null floors (bars; 1300 nulls each,
criterion-v2.2 gate passed). Certification is selective, not merely loud: from the decohered
state (floor 19.5, 4/13 clear) to the full-array ceiling (floor 7.5, 12/13 clear) the floor
falls 2.6× while the signal rises 2.96× — and the decisive state sC m1e05 recruits +2 with a
NEGATIVE signal gain (−0.20), purely on the floor drop. Honesty row: oracle-anchored with no
trials factor; the signal side is one noise draw; the plotted N_cert are the maximum over 30
banked realisations (means 5.5/3.2/2.8) — the robust readout is the relative one. (b) SPARK-3's
crossing ledger for the five crossing-eligible units under the optimistic Fisher bound and at
the true marginal width: 4/5 units survive at the marginal width, scrambled-clean
(EDGE-POSITIVE); the pessimistic (prior-width) bound gives zero crossings anywhere. (c) The
model-quality law from the same banks: the certification floor is 118 nat with the reservoir
unmodelled, 148 at truth, 276 under a conditional-width model, and 744 under a prior-width
model — a bad model does not merely fail to help, it inflates the null and suppresses
certification; this is the design law GLACIER carries.

**Sources.** (a) `SPARK_results/spark_s2c.npz` (`twoF_*`, `floor_*_floor`, `label_*`,
`corr_banked_*`). (b) `SPARK3_results/spark3_final_verdict.npz` (`unit`, `cross_opt`,
`cross_marg`, `survives`, `scrambled_clean`) + `spark3_ledger.npz` (`cross_pes`). (c)
`SPARK3_results/s3r_A_g3_r2_k8.npz` (`floor_un/tr/opt/pes`). STORY S8.6.3–S8.6.5; SPARK-3
§5.0/§5.3.

## P9 — THE PAYOFF (`P9_the_payoff`) · intended: the siren payoff

**Caption draft.** (a) Fractional luminosity-distance precision vs the number of certified
pulsar terms for the nine SIREN source cells (arm B, seed distances known to 0.1 pc): from 332%
with the Earth term alone to 6–12% at N_cert = 3–5 for the heavier pairs — the 10–30%
dark-siren-useful class (band) — because the kyr-baseline pulsar terms measure the chirp mass
(inset: σ(log₁₀Mc) falls from the prior 0.87 dex to 10⁻²–10⁻³ dex). The star is the eccentric
self-siren: 12% with ZERO certified pulsar terms (ATLAS M4, at the κ ≥ 20 self-clock threshold
e ≈ 0.50–0.53). (b) The same certified set localises the source: the 90% sky area collapses
3–4 orders of magnitude to below the one-galaxy line at 3–5 seeds. All σ are Cramér–Rao bounds
on zero-noise Asimov data with the fringe integers GIVEN — lower bounds, by the SIREN
convention. (c) The growth path on one axis, every rung a bank readback: 0.70 (a single e = 0.3
member — below the bar), 1.57 (the GEO strict-count mean over 40 skies), 7.93 (the best cell in
the corrected 108-cell box), 18 (pulsars certifying in ≥1 sky), 30 (the readable sub-array,
SNR_pterm > 1), 116 (the full array, never reached at any modelled loudness); the science needs
3 (siren-useful and host-unique) and saturates at 5.

**Sources.** (a,b) `SIREN_results/siren_summary.npz` (`B_0p1pc__frac_DL`, `__sig_log10_mc`,
`__sky_deg2`, `grid_mc`, `DL_Mpc`); star `reports/atlas_m2m4_summary.npz` (`m4`). (c)
`reports/ch_analysis_recut.npz` (m1e03_lit); `reports/geo_summary.npz` (`counts_strict`,
`freq`); `reports/surface_analysis_recut.npz` (max `corr`);
`CW_transition/stagec_p2b_results.npz` (`snr_pterm`). STORY S9.1–S9.4, S2.3.1, S3.2.1.

## P10 — GEOMETRY + WARNING (`P10_geometry_warning`) · intended: target selection / discussion

**Caption draft.** (a) The geometry lottery: the fraction of 40 isotropic source-sky draws in
which each pulsar certifies (census population, zero noise, honest gate). Only 18 of 116
pulsars ever certify; the other 98 never do, in any sky. The census names (red) are never
reproduced together (0 of 40 draws) — the count is plausibly robust, the names are not, and
sky-conditional seed sets are mandatory. (b,c) The RING warning: with today's distance priors
the sky MAP sits 4–6° from truth independent of SNR (b) while the claimed error box shrinks, so
the truth falls out of the 90% region — coverage collapses 90% → 50% → 0% as SNR rises 5 → 20
(c); with pc-class distances it stays honest at every loudness. Wrong distances move the
answer rather than blurring it: a precision distance campaign is bias prevention, not
enhancement.

**Sources.** (a) `reports/geo_summary.npz` (`freq`, `names`). (b,c)
`RING_results/npz/cell_fgw-8_A_snr{5,10,20}.npz` (`{feather,tier3}__map_offset_local_deg`,
`__inside90_local`; scenario A, 10 realisations/cell). STORY S3.3, S7.2.4–S7.2.5.

## P11 — ODDS + FORECAST (`P11_odds_forecast`) · intended: forecast section

**Caption draft.** (a–c) LOTTERY's observer panel: the probability that a random mixture of the
three loud members switches the certified count on, over the population's eccentric fraction
and characteristic eccentricity (50,000 mixture draws per cell), under the channel-budget
criterion (a) and the conservative quoted threshold rule (b), with their one-sided disagreement
(c) — the threshold rule never fires where the channel criterion would not. (d) Reserved
half-panel, PENDING: GLACIER's drain curve A_bg(iteration) — the fan is queued and no drain
bank exists yet; only the zero-noise Stage-0 precursor is banked. The layout is reserved so the
paper's forecast section composes without re-flowing when the fan lands.

**Sources.** (a–c) `reports/lottery_tier1.npz` (`P_chan`, `P_thr`, `P_dis`, `f_grid`,
`e_grid`, `conventions`). (d) status only: `GLACIER_capstone.md` (deliverable), queued jobs
observed PD; precursor `reports/glacier_g2_population.npz` (not plotted — its `drain_Abg` is a
zero placeholder). STORY S7.6, App-A GLACIER context; GAP_REPORT §6.

## P12 — THE CAVEATS TABLE (`P12_caveats_table`) · intended: discussion / scope-of-inference

**Caption draft.** The paper's claims against the domains in which each has been tested, built
from the tracked claim skeleton's scope lines: green = TESTED, naming the campaign that
measured it; amber = PENDING, naming the queued campaign; red = OPEN (the D1 wrong-counterpart
hole, for which both purity layers were pre-registered and rejected, and the real-data port —
no number in this paper has touched a real TOA). Grey dashes are not-applicable cells (e.g. the
zero-noise information-theoretic walls have no floor term to test). The sampler row-group
records that credible-interval coverage is TESTED at the banked N = 3 realisations (scenario C
vs control), SBC has never run, and the high-dimensional adaptation lesson (the chaotic Adam
M-step, retired for Newton steps) is noted.

**Sources.** STORY.md scope lines (S-references printed per row); this is the one figure whose
source is the tracked claim skeleton rather than an npz. Sampler N per
`SAMPLER_results/s4_scen{B,C}.npz` (bank supersedes the report prose; GAP_REPORT §5).

## P13 — THE HONEST POSTERIOR (`P13_honest_posterior_PENDING`) · intended: sampler section

**Caption draft (layout reserved).** The sampler figure the paper needs is two-regime by
design: (a) one above-onset realisation whose fringe-marginalised distance posterior
concentrates on the correct fringe (feasibility), (b) the same machinery at census loudness
reporting an honestly diffuse posterior (credibility — concentration tracks evidence, not
confidence), and (c) the per-parameter underestimate of linearised errors against sampled ones.
**None of the three is sourceable to a bank on this machine**: the SAMPLER posterior gates
(g1/g2/g3/SBC) did not run to completion on cronus (~7 GPU-hr/chain, host-RAM OOM) and no
chain or fringe-posterior npz exists in the repo; the σ_MLE-vs-σ_sampler comparison does not
exist in the SAMPLER record at all. The layout is reserved with each panel's missing artifact
stated; nothing is improvised. See GAP_REPORT §5.

**Sources.** `SAMPLER_dev.md` §3/§4 (gate status); absence verified by repo-wide search
(no `chain_*`/`sbc_*` npz).

## P14 — THE ARCHAEOLOGY (`P14_the_archaeology`) · intended: sampler close-out / why confidence ≠ evidence

**Caption draft.** Why previous samplers' confidence was not evidence. (a) Across the banked
loudness ladder (IGNITE signal bank, lit tier, T = 30), the median snap confidence q_max rises
0.16 → 0.94 while the purity of confident (q > 0.9) snaps falls 0.79 → 0.37: confidence rises
exactly where correctness collapses, because the same noise-lock that raises the floor gives
wrong fringes floor-beating gaps. (b) The joint quantity P(confident AND on the true fringe)
stays small (0.05–0.20) at every loudness. The MK9 corner — the regime in which delta-on-truth
sampling held (5 pulsars, one source at h = 10⁻¹², 1 μs white noise, credible-interval coverage
13/13) — sits ~56× louder than the realistic regime and is cited from its notebook-only
artifact (no npz bank; provenance printed on the panel), not plotted as data.

**Sources.** (a,b) `SAMPLER_results/archaeology_snap.npz` (`h`, `med_qmax`, `purity_q90`,
`frac_q90`, `p_confident_and_true`); raw re-derivation bank `reports/ignite_bank.npz`. MK9:
`CW_node_sampling/data_likelihood_sandbox_MK9_executed.ipynb` via `SAMPLER_dev.md` §8.2.
STORY S10; GAP_REPORT §5.

## P15 — CALIBRATION + THE SADDLE (`P15_calibration_saddle`) · intended: sampler section

**Caption draft.** (a) The S-4 coverage test, exactly as banked: the credible level at which
the truth sits under the 8-parameter loud-source sampler, for the control scenario (B,
white + RN: N = 2, both inside their 90% regions) and the full-GWB-marginalisation scenario
(C: N = 1, truth at the 94.3rd percentile, outside). The run was killed before the planned
ensemble; N is stated and nothing is extrapolated (the report prose quotes C values not present
in the npz — the bank is plotted, per the ARTIFACT READBACK convention). (b) The mechanism
panel: the banked Hessian condition numbers at truth (3.2×10⁸ in scenario C; up to 1.3×10⁹ in
B). The quoted eigenvalue spectrum — 3 of 8 directions negative, truth is a saddle of the
profiled objective, which is why mode-finding fails at truth even standing on it — was read off
a log absent from this checkout and is carried as provenance text, not drawn as data
(measured-vs-pending per the P12 convention).

**Sources.** `SAMPLER_results/s4_scenB.npz`, `s4_scenC.npz` (`c`, `inside90`, `cond`,
`seeds`). STORY S10.1.5 (RING S-4 context); GAP_REPORT §5.
