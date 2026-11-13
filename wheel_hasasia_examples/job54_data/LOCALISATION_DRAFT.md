# JOB 54 — CW localisation with a partly phased array (DRAFT)

## Five-line summary
1. Both anchors pass. All 116 pulsars phased, 10 nHz, SNR 10: 1.0e-5 deg^2 (Fisher) and 1.4e-5 deg^2 (sky likelihood); Taylor quotes ~1e-5. No pulsar phased: 330 deg^2 (Fisher) and 278 deg^2 (sky likelihood).
2. The gain from all-phased to none is 7.3–7.5 dex in area. A handful of phased pulsars buys almost none of it at SNR 10: with 3 the true 90% area shrinks by 1.1x, with 10 by 1.5x, with 30 by ~5–11x (interpolated). The collapse to the Fisher value happens between 40 and 80 phased pulsars.
3. The Fisher matrix says the opposite (3 phased: 460x smaller; 10: 10^4x). It is wrong there by factors 10^2–10^6: it measures one fringe peak, and the likelihood has hundreds to tens of thousands of equal peaks. Fisher and the sky likelihood agree (to within 30%) only once the truth peak wins: n_ph >= 80 at SNR 10, n_ph >= 20 at SNR 20.
4. The threshold is set by signal-to-noise per pulsar, not by a number of pulsars. At SNR 20 the collapse moves to n_ph ~ 10–20. In this equal-noise array each pulsar term has SNR ~ 0.65 at array SNR 10.
5. Distances must be good to a fraction of a wavelength. With every pulsar known to 0.3 pc the area is 3.5x the 0.01 pc value; at 1 pc the localisation is lost (64 deg^2); at 3 pc nothing is left (185 deg^2 against 278).

## Where the brief is contradicted by the code or the physics
- **Gaussian census priors in the Fisher matrix give a false n_ph = 0 area.** The brief asks for "every pulsar's distance with its Gaussian prior (add 1/sigma_L^2 on the diagonal)". Done literally, n_ph = 0 gives 0.42 deg^2 at SNR 10 and 0.33 deg^2 at SNR 20 (10 nHz; `job54_loc_fisher.py`, `__main__`, mode `census`). That is 800x below the Earth-term anchor and does not scale with SNR. Cause: census errors (1–9611 pc, median 232 pc) are many wavelengths (0.97 pc at 10 nHz). The pulsar-term phase wraps hundreds of times inside the prior. A linearised Gaussian cannot represent that; it treats a 20% distance prior as a constraint on (1 - cos mu). The correct statement is that an unphased pulsar-term phase is unknown. So in measure (a) unphased pulsars get no distance prior (sigma = infinity; `sig_vec(..., unphased="free")`). The literal version is kept in the CSV as measure `a_census`.
- **Measures (a) and (b) are not the same quantity.** (a) marginalises the six other source parameters. (b) holds them at truth. The like-for-like Fisher number is the conditional one (inverse of the 2x2 sky block of F), in the CSV as `a_cond`. It is 1.6–1.7x smaller than (a).
- **prong2's residual has no frequency evolution.** `/home/mattm/projects/HSYMT/CW_transition/prong2_transition.py:19-20`: "Pulsar-term amplitude == Earth-term amplitude (non-evolving), so distance enters ONLY / through the pulsar-term phase  -2π f (L/c)(1 - cos μ)." Lines 87–90: `tau_p = (L_kpc * KPC_OVER_C) * (1.0 - cos_mu)`, `Phi_e = 2.0 * jnp.pi * f * toas + phase0`, `Phi_p = Phi_e - 2.0 * jnp.pi * f * tau_p`. Line 98: `return zeta * (wave(Phi_e) - wave(Phi_p))`. It has antenna patterns (`antenna`, lines 60–74), inclination weights (lines 92–93) and amplitude zeta = h/(2 pi f) (line 86). It has no chirp mass. I re-implemented the residual with evolution (`job54_loc_core.py`, `residual` and `evol`). With the chirp mass set to 1e-6 Msun it equals `single_source_residual` to 4.0e-11 relative (50 random parameter sets; `job54_loc_test.py`, `test_residual_vs_prong2`). Tested.
- **Evolution at 1e9 Msun is not negligible in phase, and is in frequency only at 10 nHz.** For 1 kpc of retardation (`job54_loc_test.py`, `evolution_size`): 10 nHz: f_p/f = 0.9736, shift 0.26 nHz = 0.12 cycles over 15 yr (not resolved from the Earth term), extra phase -87 rad. 20 nHz: f_p/f = 0.866, shift 2.7 nHz = 1.3 cycles over the span (resolved), extra phase -945 rad. Drift inside the 15-yr span is ignored (phase error 0.002 rad at 10 nHz, 0.024 rad at 20 nHz).
- **The fringe-spacing estimate in the brief is right, and the grid cannot resolve the farthest pulsars over the whole n_ph = 0 region.** The drawn array contains a 49.9 kpc and two 15.7 kpc pulsars. Their fringes are 1.1e-3 deg and 3.5e-3 deg at 10 nHz. The n_ph = 0 99.9% region is ~1000 deg^2. See "Grid" for what was done and why the areas still hold.
- **Two processes not started by me were active in the same directory.** Another job (`job54_run.py`) shares `job54_data/` and `job54_logs/`; my files there are `fisher_*.np*`, `numbers.json`, `tables.md`, `tests.json`, `sky/`, and logs `fisher.log`, `test.log`, `post.log`, `sky_queue*.log`, `queue*.sh`, `queue.pid`. `git status` now shows my `job54_loc_*.py` and `talk_data/localisation.csv` as staged (`A`). I did not run `git add`; something else did. I left the index alone.

## Array used
- 116 pulsars at the sky positions of `prong2_transition.build_array(116, seed=777)` (isotropic). 15.0 yr span (4.7336e8 s), 7-day cadence (782 TOAs per pulsar, `build_toas`), white noise 100 ns per TOA, equal for all pulsars. No red noise, no timing-model fit.
- Distances: 116 (D, sigma) pairs drawn with replacement from `job45_t2_parallax_census.csv`, seed 54 (`job54_loc_core.build_all`). 64 distinct census rows. Distance min / median / max = 156 / 1346 / 49868 pc. Census error min / median / max = 1.0 / 232 / 9611 pc.
- Drawn pulsars with present error below 0.4 pc: 0. Below 1 pc: 0 (2 have exactly 1.0 pc). Below 3 pc: 7. (`job54_loc_post.tables`, key `census` in `numbers.json`.)
- Orderings: nearest-first = sort by distance. Random = `default_rng(5401).permutation(116)`. Both nested. The random order starts with far pulsars; the nearest-first set at n_ph = 10 reaches only 329 pc.
- Sources: 200, seed 5400 (`draw_sources`): isotropic sky, cos(inclination) uniform, polarisation and phase uniform. Chirp mass 1e9 Msun. Earth-term frequency 10 or 20 nHz. Amplitude scaled so the summed matched-filter SNR over all pulsars (Earth + pulsar terms) is 10 or 20 (`snr_unit_h`).

## Formulas
**Signal** (`job54_loc_core.residual`, `evol`). tau = L(1 - cos mu)/c. x = (256/5)(G Mc/c^3)^(5/3) (pi f)^(8/3) tau. f_p = f (1+x)^(-3/8). dPhi = 2 pi f tau [(1+x)^(5/8) - 1]/(5x/8). r(t) = zeta W(2 pi f t + phi0) - zeta (f/f_p)^(1/3) W(2 pi f_p t + phi0 - dPhi), with W(Phi) = F+ (1+cos^2 i) sin Phi - 2 Fx cos i cos Phi and zeta = h/(2 pi f).

**(a) Fisher** (`job54_loc_fisher.py`: `blocks`, `F8`, `area_from_F`). Parameters: theta, phi, ln f, ln h, phi0, inclination, psi, ln Mc, and one distance per pulsar (pc). Jacobians by `jax.jacfwd` of the residual in float64. Per pulsar a: A = Jg^T Jg / sigma^2, b = Jg^T jL / sigma^2, c = jL^T jL / sigma^2. Its distance is marginalised analytically: M_a = A - b b^T/(c + 1/sigma_La^2), computed as A_perp + b b^T p/(c(c+p)) with p = 1/sigma_La^2 and A_perp built from the Jacobian with the distance direction projected out. F = sum_a M_a (8x8). Phased: sigma_La = sigma_L. Unphased: sigma_La = infinity. Sigma = F^-1 after scaling F to unit diagonal.
A90 = 2 pi ln(10) sin(theta) sqrt(det Sigma[theta,phi]) steradian, times (180/pi)^2 for deg^2.
Stability at SNR 10, 200 sources x 2 frequencies, n_ph = 0, 4, all (`job54_loc_fisher.py __main__`, `job54_logs/fisher.log`): distances in kpc instead of pc changes the area by at most 8e-4; a naive (unprojected) Schur complement by at most 7e-3; a pseudo-inverse with threshold 1e-13 or 1e-10 by at most 2.5e-2 (median change 2e-15; the few-percent cases are near-face-on sources where inclination is nearly degenerate).

**(b) sky likelihood** (`job54_loc_core.parts_fn`, `job54_loc_sky.py`: `cases_main`, `run_pass`). Noiseless data: this is the expected log-likelihood, not a noise realisation. All source parameters except sky position are at truth. For pulsar a, data d = E0 + P0. The template at sky position n is E(n) + P(n), with P(n) = Re[c_P e^{i w_p t}]; c_P carries amplitude, antenna pattern and the phase phi0 - dPhi(n); w_p = 2 pi f_p(n).
- Phased pulsar (coherent, true distance): lnL_a = -|d - E(n) - P(n)|^2 / (2 sigma^2). Exact. Zero at the truth.
- Unphased pulsar: the pulsar-term phase is replaced by a free phase, uniform on [0, 2 pi). Define z = c_P sum_t (d - E(n))(t) e^{i w_p t}. Then lnL_a = -[ |d - E(n)|^2 + (1/2)|c_P|^2 N ] / (2 sigma^2) + ln I0(|z| / sigma^2). The Earth–pulsar cross term is inside z, because z uses d - E(n). Approximations: (i) the template norm drops the term (1/2) Re[c_P^2 sum_t e^{2 i w_p t}] (it is <= 3% of the norm at 10 nHz); (ii) the pulsar-term frequency and amplitude are evaluated at the catalogue (true) distance, only the phase is free. (ii) is weakest at 20 nHz for the worst-known pulsars (f_p shifts 0.27 nHz per 100 pc).
- Total: lnL(n) = sum_a inc_a(n) + sum_{a phased} [coh_a(n) - inc_a(n)]. Nested sets are cumulative sums, so one grid pass gives every n_ph and both orderings.
- Partial knowledge (`cases_partial`): the distance is marginalised over its Gaussian prior through the phase only. With cross term a cos(alpha) and sigma_phi = w_p sigma_L (1 - cos mu)/c: E[exp(a cos(alpha + dphi))] = sum_m I_m(a) exp(-m^2 sigma_phi^2 / 2) cos(m alpha), m = -24..24, I_m by Miller recurrence. It tends to the coherent form as sigma_L -> 0 and to I0 as sigma_L -> infinity. The 2-phi norm term is dropped here too.
- All inner products are closed-form sums of sinusoids on the regular TOA grid (Dirichlet kernels), not sums over TOAs.
Tests (`job54_loc_test.py`, `test_closed_form`, results in `job54_data/tests.json`): coherent lnL against a brute-force sum over TOAs at 41 sky points x 116 pulsars: max abs error 8e-11 (10 nHz SNR 10), 6e-10 (20 nHz SNR 20). Phase-marginalised form against exact numerical marginalisation over the phase: max abs error 0.003 (SNR 10, both frequencies), 0.048 (20 nHz, SNR 20). Partial form at 0.3 pc against numerical integration over distance: 0.021 / 0.008 / 0.069 in lnL.

**Grid** (`job54_loc_sky.run_pass`). Coordinates: longitude/latitude (u, v) in a frame rotated so the truth is at (0, 0); cell solid angle h^2 cos v. Steps: (1) whole sky at 0.5 deg with n_ph = 0; take the 99.9% set. (2) 512x512 over its bounding box, 99.9% set again, grown by one cell: this is the mask (median 990 deg^2 at 10 nHz SNR 10; the n_ph = 0 99% area is 584 deg^2). (3) Uniform fine grid over the mask with 3e7 cells (SNR 10) or 1e7 cells (SNR 20): h = 0.0037–0.0068 deg (10 nHz, SNR 10). (4) The 12x12 base cells around the truth are replaced by cells 64x finer (5.8e-5 to 1.1e-4 deg). (5) The 90% area is the area of the highest-density cells holding 90% of the probability, from area and mass histograms in lnL (bin 0.02). Islands: 8-connected components of that set on the base grid (on the zoom grid when the whole set is inside the zoom window).
Step against finest fringe (period = lambda_p / (L sin mu) at the truth; median over sources of period/h for the phased set):

| n_ph (nearest-first) | 1 | 2 | 3 | 4 | 6 | 10 | 20 | 40 | 80 | all |
|---|---|---|---|---|---|---|---|---|---|---|
| period / h, 10 nHz SNR 10 | 78 | 69 | 64 | 64 | 45 | 33 | 16 | 9.9 | 5.0 | 0.23 |
| period / h, 20 nHz SNR 10 | 49 | 41 | 40 | 40 | 27 | 21 | 10 | 6.1 | 3.1 | 0.14 |

Random ordering: period/h = 1.7 for n_ph = 1–3, 0.8 for 4–6, 0.23 from 10 on (10 nHz). So the base grid resolves the nearest-first sets up to n_ph = 80 (40 at 20 nHz) and does not resolve the random sets or the all-phased set away from the truth.
Rule used to accept an area (`job54_loc_post.load_b`): the base grid has >= 4 cells per finest fringe, OR >= 90% of the probability is inside the zoom window, OR the 90% set holds >= 1000 base cells. The third case is a broad set where an under-sampled fringe is sampled at effectively random phase; the area is then a sampling estimate and is right on average. The half-step test below supports this. Island counts are reported only when the base grid resolves the fringes or the set is entirely inside the zoom window; otherwise "n.r.".
Cases failing all three: 9 of 20 (10 nHz) and 10 of 20 (20 nHz) at SNR 20, random order, n_ph = 20. They are excluded from those medians. No other case failed.
The maximum of lnL over the grid minus its value at the truth was at most 2.0 (main runs) and 3.2 (partial runs); it is not exactly zero because the phase-marginalised likelihood need not peak at the truth. The truth was inside the 90% set in every case.

## Anchor checks
| | this work (a) Fisher, 200 sources | this work (b), 20 sources | reference |
|---|---|---|---|
| all phased, 0.01 pc, SNR 10, 10 nHz | 1.04e-5 [7.2e-6, 1.5e-5] deg^2 | 1.40e-5 [1.1e-5, 2.0e-5] deg^2 | Taylor ~1e-5 deg^2 |
| all phased, 20 nHz | 3.6e-6 [2.6e-6, 5.0e-6] | 3.6e-6 [2.9e-6, 5.6e-6] | |
| none phased, SNR 10, 10 nHz | 330 [147, 532] deg^2 | 278 [175, 377] deg^2 | tens to hundreds deg^2 |
| none phased, 20 nHz | 268 [139, 419] | 184 [149, 241] | |

Per source, (b)/(a), median [16th, 84th] (`job54_loc_post.tables`, keys `check_*`): all phased, 10 nHz SNR 10: 1.20 [1.06, 1.41]; 20 nHz: 1.03 [0.91, 1.12]. None phased: 0.77 [0.63, 1.06] and 0.59 [0.47, 0.92]. Both anchors pass, and (a) and (b) agree within a factor 2.
Against the like-for-like conditional Fisher: all phased (b)/(a_cond) = 1.93 [1.46, 3.22] at 10 nHz SNR 10, 1.30 at 20 nHz, 1.10 and 1.04 at SNR 20. With 80 nearest phased it is 1.07. The excess at "all" comes from the three pulsars beyond 15 kpc: Fisher credits their very fine fringes in full, but each has pulsar-term SNR below 1 and only modulates the likelihood. It goes away at SNR 20. This reading is an interpretation; it was not tested by removing those pulsars.
The n_ph = 0 area is at the top of the "tens to hundreds" range. Not checked against a specific published number. The 15-yr span, equal 100 ns noise and other parameters fixed (in (b)) all differ from any one paper.

## Main tables

**10 nHz, SNR 10.** 90% area in deg^2, median [16th, 84th percentile]. (a): 200 sources. (b): 20 sources; after the slash, the median number of islands (n.r. = not resolved).

| n_ph | (a) nearest | (a) random | (a_cond) nearest | (b) nearest / islands | (b) random / islands | (b) truth-island area, nearest (deg^2) | (b) cases unresolved near/rand |
|---|---|---|---|---|---|---|---|
| 0 | 330 [147, 532] | 330 [147, 532] | 204 [93.5, 376] | 278 [175, 377] / 1 | 278 [175, 377] / 1 | 278 | 0/0 |
| 1 | 40.8 [19.7, 105] | 86.1 [22, 255] | 7.05 [3.61, 24.7] | 278 [169, 366] / 6 | 271 [174, 375] / 3 | 277 | 0/0 |
| 2 | 4.89 [1.41, 18.8] | 14.2 [1.18, 47.3] | 0.678 [0.178, 3.74] | 260 [168, 358] / 20 | 255 [173, 356] / 86 | 240 | 0/0 |
| 3 | 0.714 [0.386, 3.9] | 0.22 [0.095, 1.25] | 0.179 [0.102, 0.891] | 253 [161, 358] / 168 | 238 [156, 325] / 6373 | 229 | 0/0 |
| 4 | 0.273 [0.135, 0.708] | 0.00681 [0.00246, 0.0455] | 0.112 [0.0768, 0.272] | 234 [149, 340] / 286 | 228 [153, 323] / n.r. | 210 | 0/0 |
| 6 | 0.129 [0.0705, 0.208] | 0.000498 [0.000234, 0.00165] | 0.0577 [0.0367, 0.0962] | 221 [144, 315] / 565 | 180 [152, 255] / n.r. | 188 | 0/0 |
| 10 | 0.0315 [0.0195, 0.05] | 0.00018 [6.25e-05, 0.000382] | 0.0215 [0.0151, 0.029] | 184 [109, 238] / 1341 | 157 [120, 221] / n.r. | 135 | 0/0 |
| 20 | 0.00381 [0.00287, 0.00595] | 6.92e-05 [3.31e-05, 0.00011] | 0.00276 [0.00214, 0.00443] | 115 [57.2, 172] / 7450 | 105 [75.2, 136] / n.r. | 39.9 | 0/0 |
| 40 | 0.000536 [0.000447, 0.00067] | 3.01e-05 [1.65e-05, 4.22e-05] | 0.000474 [0.000405, 0.000572] | 8.61 [0.0088, 59.5] / 9154 | 36.1 [19, 60.2] / n.r. | 0.00953 | 0/0 |
| 80 | 9.95e-05 [8.63e-05, 0.000116] | 1.24e-05 [8.07e-06, 1.85e-05] | 9.33e-05 [8.01e-05, 0.000108] | 0.000105 [8.83e-05, 0.000129] / 1 | 1.91e-05 [1.52e-05, 2.87e-05] / 1 | 0.000105 | 0/0 |
| all (116) | 1.04e-05 [7.15e-06, 1.49e-05] | 1.04e-05 [7.15e-06, 1.49e-05] | 5.97e-06 [3.78e-06, 1.04e-05] | 1.4e-05 [1.07e-05, 2.03e-05] / 1 | 1.4e-05 [1.07e-05, 2.03e-05] / 1 | 1.4e-05 | 0/0 |

**10 nHz, SNR 20.** 90% area in deg^2, median [16th, 84th percentile]. (a): 200 sources. (b): 20 sources; after the slash, the median number of islands (n.r. = not resolved).

| n_ph | (a) nearest | (a) random | (a_cond) nearest | (b) nearest / islands | (b) random / islands | (b) truth-island area, nearest (deg^2) | (b) cases unresolved near/rand |
|---|---|---|---|---|---|---|---|
| 0 | 82.4 [36.9, 133] | 82.4 [36.9, 133] | 51.1 [23.4, 94] | 79.1 [48.4, 111] / 1 | 79.1 [48.4, 111] / 1 | 79.1 | 0/0 |
| 1 | 10.2 [4.92, 26.3] | 21.5 [5.51, 63.8] | 1.76 [0.904, 6.18] | 53 [24.4, 97.5] / 10 | 50.6 [24.5, 78.3] / 1 | 33.9 | 0/0 |
| 2 | 1.22 [0.353, 4.71] | 3.55 [0.296, 11.8] | 0.17 [0.0447, 0.937] | 41.5 [18.5, 90.3] / 34 | 47.4 [16.7, 65.3] / 170 | 19 | 0/0 |
| 3 | 0.179 [0.0966, 0.976] | 0.0549 [0.0238, 0.313] | 0.0448 [0.0254, 0.223] | 37.4 [12.7, 86] / 182 | 21.7 [9.49, 55.5] / 7114 | 3.9 | 0/0 |
| 4 | 0.0684 [0.0339, 0.177] | 0.00171 [0.000616, 0.0114] | 0.028 [0.0193, 0.0682] | 19.9 [5.89, 43.8] / 300 | 17.3 [6.08, 25.8] / n.r. | 1.29 | 0/0 |
| 6 | 0.0323 [0.0176, 0.0521] | 0.000125 [5.86e-05, 0.000413] | 0.0145 [0.0092, 0.0241] | 10.5 [1.98, 31.3] / 331 | 5.93 [2.03, 19.6] / n.r. | 0.052 | 0/0 |
| 10 | 0.00791 [0.00489, 0.0125] | 4.5e-05 [1.57e-05, 9.56e-05] | 0.00538 [0.00379, 0.00727] | 0.981 [0.0712, 7.66] / 256 | 1.81 [0.559, 4.09] / n.r. | 0.0182 | 0/0 |
| 20 | 0.000954 [0.00072, 0.00149] | 1.73e-05 [8.3e-06, 2.76e-05] | 0.000691 [0.000537, 0.00111] | 0.000699 [0.000549, 0.0402] / 1 | 0.0597 [5.13e-05, 0.132] / 7 | 0.000699 | 0/9 |
| 40 | 0.000135 [0.000112, 0.000168] | 7.53e-06 [4.14e-06, 1.06e-05] | 0.000119 [0.000102, 0.000144] | 0.000128 [0.000105, 0.000154] / 1 | 5.32e-06 [3.14e-06, 1.02e-05] / 1 | 0.000128 | 0/0 |
| 80 | 2.49e-05 [2.16e-05, 2.9e-05] | 3.1e-06 [2.02e-06, 4.62e-06] | 2.34e-05 [2.01e-05, 2.71e-05] | 2.49e-05 [2.1e-05, 2.99e-05] / 1 | 2.09e-06 [1.27e-06, 3.67e-06] / 1 | 2.49e-05 | 0/0 |
| all (116) | 2.61e-06 [1.79e-06, 3.74e-06] | 2.61e-06 [1.79e-06, 3.74e-06] | 1.49e-06 [9.47e-07, 2.6e-06] | 1.8e-06 [1.14e-06, 3.11e-06] / 1 | 1.8e-06 [1.14e-06, 3.11e-06] / 1 | 1.8e-06 | 0/0 |

**20 nHz, SNR 10.** 90% area in deg^2, median [16th, 84th percentile]. (a): 200 sources. (b): 20 sources; after the slash, the median number of islands (n.r. = not resolved).

| n_ph | (a) nearest | (a) random | (a_cond) nearest | (b) nearest / islands | (b) random / islands | (b) truth-island area, nearest (deg^2) | (b) cases unresolved near/rand |
|---|---|---|---|---|---|---|---|
| 0 | 268 [139, 419] | 268 [139, 419] | 172 [88.6, 311] | 184 [149, 241] / 1 | 184 [149, 241] / 1 | 184 | 0/0 |
| 1 | 21.7 [10.3, 51.8] | 52.4 [11.9, 158] | 3.21 [1.68, 12.9] | 177 [143, 240] / 7 | 170 [136, 224] / 1 | 169 | 0/0 |
| 2 | 1.49 [0.405, 6.63] | 4.84 [0.585, 15.6] | 0.175 [0.0466, 1.14] | 169 [135, 228] / 64 | 168 [107, 208] / 378 | 167 | 0/0 |
| 3 | 0.206 [0.107, 1.11] | 0.063 [0.0313, 0.349] | 0.0454 [0.0274, 0.236] | 167 [128, 224] / 479 | 162 [103, 207] / 12586 | 160 | 0/0 |
| 4 | 0.0774 [0.0376, 0.183] | 0.00232 [0.000707, 0.0111] | 0.0298 [0.0197, 0.0721] | 153 [128, 213] / 646 | 149 [103, 193] / n.r. | 139 | 0/0 |
| 6 | 0.0378 [0.02, 0.0684] | 0.000155 [6.45e-05, 0.000425] | 0.0151 [0.00979, 0.0256] | 148 [111, 180] / 1504 | 130 [87.9, 186] / n.r. | 119 | 0/0 |
| 10 | 0.00905 [0.00518, 0.0142] | 5.63e-05 [1.5e-05, 0.000105] | 0.00571 [0.00394, 0.00774] | 122 [79.6, 166] / 3898 | 110 [71.8, 145] / n.r. | 76.5 | 0/0 |
| 20 | 0.00101 [0.000767, 0.00162] | 2.04e-05 [9.6e-06, 3.6e-05] | 0.000737 [0.000602, 0.0012] | 70.5 [38.9, 125] / 17417 | 64.6 [46.6, 98.4] / n.r. | 25.5 | 0/0 |
| 40 | 0.00015 [0.000126, 0.000183] | 9.62e-06 [5.74e-06, 1.44e-05] | 0.000133 [0.000117, 0.000156] | 13.4 [0.211, 40.9] / 35570 | 21.7 [12.6, 43] / n.r. | 0.00379 | 0/0 |
| 80 | 2.87e-05 [2.57e-05, 3.3e-05] | 4.29e-06 [2.99e-06, 6.33e-06] | 2.73e-05 [2.39e-05, 3.09e-05] | 3.09e-05 [2.64e-05, 3.4e-05] / 1 | 5.46e-06 [3.7e-06, 8.24e-06] / 1 | 3.09e-05 | 0/0 |
| all (116) | 3.59e-06 [2.64e-06, 4.97e-06] | 3.59e-06 [2.64e-06, 4.97e-06] | 2.41e-06 [1.69e-06, 3.88e-06] | 3.61e-06 [2.88e-06, 5.61e-06] / 1 | 3.61e-06 [2.88e-06, 5.61e-06] / 1 | 3.62e-06 | 0/0 |

**20 nHz, SNR 20.** 90% area in deg^2, median [16th, 84th percentile]. (a): 200 sources. (b): 20 sources; after the slash, the median number of islands (n.r. = not resolved).

| n_ph | (a) nearest | (a) random | (a_cond) nearest | (b) nearest / islands | (b) random / islands | (b) truth-island area, nearest (deg^2) | (b) cases unresolved near/rand |
|---|---|---|---|---|---|---|---|
| 0 | 67 [34.9, 105] | 67 [34.9, 105] | 42.9 [22.2, 77.7] | 34 [26.1, 43.4] / 1 | 34 [26.1, 43.4] / 1 | 34 | 0/0 |
| 1 | 5.42 [2.58, 13] | 13.1 [2.97, 39.4] | 0.806 [0.42, 3.23] | 28.2 [15.2, 41.3] / 10 | 23.5 [13.2, 36.1] / 1 | 27 | 0/0 |
| 2 | 0.375 [0.102, 1.66] | 1.21 [0.146, 3.9] | 0.0443 [0.0118, 0.286] | 22.6 [11.4, 35.5] / 80 | 20.6 [9.05, 27.6] / 360 | 13.6 | 0/0 |
| 3 | 0.0518 [0.0269, 0.279] | 0.0158 [0.00783, 0.0872] | 0.0115 [0.00688, 0.0595] | 15.1 [9.33, 33.9] / 388 | 11.7 [4.53, 22.2] / 11203 | 0.551 | 0/0 |
| 4 | 0.0194 [0.00949, 0.0458] | 0.000582 [0.000178, 0.00278] | 0.00751 [0.00498, 0.0183] | 10.8 [5.03, 17.9] / 431 | 6.2 [3.68, 9.88] / n.r. | 0.425 | 0/0 |
| 6 | 0.00954 [0.00505, 0.0173] | 3.89e-05 [1.62e-05, 0.000107] | 0.00382 [0.00248, 0.00649] | 3.62 [1.85, 15.2] / 570 | 2.27 [1.14, 6.52] / n.r. | 0.0151 | 0/0 |
| 10 | 0.00227 [0.0013, 0.0036] | 1.42e-05 [3.76e-06, 2.63e-05] | 0.00144 [0.000995, 0.00195] | 0.586 [0.167, 2.97] / 348 | 0.574 [0.241, 1.37] / n.r. | 0.00561 | 0/0 |
| 20 | 0.000254 [0.000193, 0.000409] | 5.12e-06 [2.41e-06, 9.09e-06] | 0.000186 [0.000152, 0.000302] | 0.000179 [0.000158, 0.0139] / 1 | 3.51e-05 [7.96e-06, 0.0452] / 4 | 0.000179 | 0/10 |
| 40 | 3.79e-05 [3.19e-05, 4.63e-05] | 2.42e-06 [1.44e-06, 3.61e-06] | 3.37e-05 [2.95e-05, 3.94e-05] | 3.69e-05 [3.05e-05, 4.05e-05] / 1 | 1.33e-06 [6.86e-07, 3.21e-06] / 1 | 3.69e-05 | 0/0 |
| 80 | 7.23e-06 [6.48e-06, 8.34e-06] | 1.08e-06 [7.54e-07, 1.59e-06] | 6.87e-06 [6.03e-06, 7.79e-06] | 7.33e-06 [6.32e-06, 8.07e-06] / 1 | 7.67e-07 [5.08e-07, 1.24e-06] / 1 | 7.34e-06 | 0/0 |
| all (116) | 9.02e-07 [6.64e-07, 1.25e-06] | 9.02e-07 [6.64e-07, 1.25e-06] | 6.06e-07 [4.25e-07, 9.77e-07] | 6.75e-07 [4.58e-07, 1.05e-06] / 1 | 6.75e-07 [4.58e-07, 1.05e-06] / 1 | 6.73e-07 | 0/0 |

**Partial knowledge, 10 nHz, SNR 10, nearest-first.** 90% area in deg^2, median [16th, 84th]. (a_lin) = linearised Fisher with a Gaussian distance prior sigma_L on the phased set (200 sources). (b) = distance marginalised over the prior (20 sources) / median islands.

| sigma_L (pc) | n_ph=4 (a_lin) | n_ph=4 (b) | n_ph=10 (a_lin) | n_ph=10 (b) | all (a_lin) | all (b) |
|---|---|---|---|---|---|---|
| 0.01 | 0.273 [0.135, 0.708] | 234 [149, 340] / 286 | 0.0315 [0.0195, 0.05] | 184 [109, 238] / 1341 | 1.04e-05 [7.15e-06, 1.49e-05] | 1.4e-05 [1.07e-05, 2.03e-05] / 1 |
| 0.1 | 0.288 [0.149, 0.731] | 235 [154, 343] / 210 | 0.0342 [0.021, 0.0532] | 189 [114, 248] / 1182 | 1.11e-05 [7.64e-06, 1.61e-05] | 1.65e-05 [1.32e-05, 2.38e-05] / 1 |
| 0.3 | 0.432 [0.225, 0.95] | 243 [171, 357] / 71 | 0.0511 [0.0314, 0.0779] | 213 [136, 278] / 619 | 1.56e-05 [1.11e-05, 2.39e-05] | 4.94e-05 [3.08e-05, 6.38e-05] / 1 |
| 1.0 | 1.74 [0.706, 4.23] | 261 [172, 366] / 4 | 0.182 [0.107, 0.329] | 254 [142, 337] / 102 | 5.83e-05 [3.65e-05, 9.27e-05] | 63.9 [30.7, 111] / n.r. |
| 3.0 | 8.98 [3.43, 21.9] | 275 [172, 366] / 1 | 0.956 [0.559, 1.93] | 263 [151, 363] / 7 | 0.000332 [0.000194, 0.000546] | 185 [81.3, 247] / n.r. |
| 10.0 | 49.4 [22.1, 117] | 278 [172, 366] / 1 | 7.43 [4.33, 17.2] | 271 [164, 366] / 1 | 0.00289 [0.00154, 0.00514] | 254 [131, 338] / n.r. |
In the partial table the 0.01 pc row and all (b) columns use the same 20 sources. For "all" at 1 pc and above, the base step (0.012–0.019 deg, 4e6 cells) does not resolve far-pulsar fringes; the area is accepted under the broad-set rule and is the least secure entry.

Source of every table number: (a) columns `job54_loc_fisher.py` (`blocks`, `F8`, `area_from_F`, `sig_vec`), file `job54_data/fisher_areas.npz`; (b) columns `job54_loc_sky.py` (`run_pass`, `cases_main`, `cases_partial`, `level`), files `job54_data/sky/*.npz`; medians and percentiles `job54_loc_post.py` (`tables`, `load_b`, `q`).

## Grid-resolution test
Source 2, 10 nHz, SNR 10, step halved (0.00641 -> 0.00320 deg; 3.09e7 -> 1.23e8 cells; `job54_loc_sky.run_pass` with `h_scale=0.5`, compared in `job54_loc_post.tables`, key `half_step`).
- 90% area: changes by at most 2.0e-4 (fraction) over all 22 cases (11 n_ph, both orderings). This includes the random-order cases whose fringes the base grid does not resolve. Requirement was < 10%.
- Islands, nearest-first: 1, 1, 2, 66, 309, 890 -> 891, 2160 -> 2162, 7966 -> 8022, 28608 -> 29184, 1, 1 for n_ph = 0…all. Stable to 2%.
- Islands, random order: change by factors 1.2–4. Not resolved, as flagged.
- Truth-island area at n_ph = 20: 70.6 -> 59.2 deg^2 (16% change). The set is near percolation there; treat truth-island areas at n_ph = 10–40 as rough.
Also: source 0 run at 2e6 and 3e7 cells (step 0.026 and 0.0068 deg) gave the same areas to 1% at every n_ph (logs only; first run not kept).
Only one source was run at half step. Not tested at 20 nHz or SNR 20.

## Answers
- **How much of the all-phased gain is reached.** Measure (b), SNR 10, 10 nHz (20 nHz in brackets), nearest-first: 3 phased give 1.10x (1.10x) = 0.6% (0.5%) of the 7.3 (7.7) dex; 10 give 1.5x (1.5x) = 2.5% (2.3%); 30 give ~11x (~7x) = 14% (11%), interpolated in log–log between 20 and 40, where the curve is steep, so uncertain by a factor of a few. Random order: 1.17x, 1.8x, ~5x. Measure (a) claims 460x (36%), 1.0e4x (54%) and 2.8e5x (73%) for the same three; those are the size of one peak, not the region. At SNR 20, 10 nHz, (b) nearest-first: 2.1x, 81x, ~3e5x (4%, 25%, 72%).
- **One region or many.** Many. At SNR 10 with 3 nearest pulsars phased the 90% set is 168 islands (10 nHz; 479 at 20 nHz), with 4 it is 286, with 10 it is 1341, with 20 it is 7450, and their total area is still 115–253 deg^2. All 20 sources had more than one island at n_ph = 3. A single small region appears only at n_ph = 80 (SNR 10) or n_ph = 20 (SNR 20).
- **What the expected supply buys.** Scaled to 116 pulsars the supply is 3.4–5.3 recovered and 0–0.5 certified distances (real priors), 5.7–9.7 recovered and 0.6–2.2 certified (Earth terms removed). At SNR 10 that is the n_ph = 1–10 range of measure (b): the area shrinks by 1.0–1.5x (278 -> 184–278 deg^2 at 10 nHz) and breaks into tens to ~1000 candidate patches; certified distances alone (<= 2) buy under 1.07x. At SNR 20 the same 3–10 distances buy 2–80x. Only the "every source known, Earth terms removed" case (39 certified, 58 recovered) reaches the edge of the collapse at SNR 10: n_ph = 40 gives 8.6 deg^2 median with a 16–84% range of 0.009–60 deg^2. And these distances must be good to ~0.1–0.3 pc.

## Pre-registered predictions
- **W1 both anchors pass: PASS.** All phased 1.04e-5 (a) and 1.40e-5 (b) deg^2 against ~1e-5; none phased 330 (a) and 278 (b) deg^2 against tens to hundreds. This holds only with unphased phases free in the Fisher matrix; with literal census priors the n_ph = 0 Fisher area is 0.42 deg^2 and the anchor fails by 800x.
- **W2 three random phased pulsars give a 5–100x smaller area and more than one island in most cases: FAIL on area, island part not decided for random order.** SNR 10: ratio of medians 1.17 (10 nHz) and 1.14 (20 nHz); per source 1.13 [1.03, 1.32]. Islands for the random set were resolved for only 1 of 20 sources (that one had many). For the nearest-first set, where the grid resolves the fringes, 20 of 20 sources have more than one island (median 168). At SNR 20 the ratio is 3.7 (10 nHz) and 2.9 (20 nHz): still below 5.
- **W3 (a) is more than 10x below (b) for n_ph <= 4 and they converge by n_ph ~ 20: first half mostly true, second half FAIL at SNR 10.** Median (b)/(a) per source, nearest-first, 10 nHz SNR 10: n_ph = 1: 4.4; 2: 33; 3: 183; 4: 789; 10: 6.0e3; 20: 2.2e4; 40: 1.8e4; 80: 1.00; all: 1.20. So at n_ph = 1 the gap is 4x, not 10x, and at 20 the gap is at its largest. Convergence is at n_ph = 80. At SNR 20 it is 0.80 at n_ph = 20: there the prediction holds.
- **W4 the all-phased area is not reached within 10x until n_ph >= 20: PASS, with margin.** Measure (b), SNR 10, nearest-first: first n_ph within 10x of all-phased is 80 (7.5x at 10 nHz; 8.6x at 20 nHz). Random order: also 80 (1.4x). Measure (a): 80 nearest-first, 20 random. At SNR 20, (b): not until all for nearest-first (14x at n_ph = 80 at 10 nHz, 11x at 20 nHz), 40 for random.
- **W5 with the phased set known to 3 pc the gain at n_ph = 10 is below 3x: PASS.** Measure (b), 20 sources, 10 nHz, SNR 10: gain 1.06 (ratio of medians), 1.02 [1.00, 1.15] per source. The linearised Fisher number is 345x and should not be used at 3 pc. At 0.01 pc the n_ph = 10 gain is itself only 1.5x, so this prediction passes for a weaker reason than intended.

## Limits
- One source, no background, white noise, equal noise in every pulsar. This is the gain available before confusion. Equal noise spreads the SNR evenly (each pulsar term has SNR ~ 0.65 at array SNR 10); a real array with a few dominant pulsars would move the collapse to fewer phased pulsars if those are the phased ones. Not tested.
- Measure (b) holds frequency, amplitude, phase, inclination, polarisation and chirp mass at truth, uses noiseless data, and counts islands of the expected likelihood. A real noise draw moves the peaks. Marginalising frequency and chirp mass would widen (b) further for small n_ph. Fisher (a) marginalises them and is 1.6–1.7x above the conditional value.
- Which pulsars are phased matters. In Fisher (a), 10 nHz SNR 10, over 30 random nested orderings the median area spans 0.053–6.6 deg^2 at n_ph = 3 (nearest-first 0.80), 2.1e-4–2.7e-3 deg^2 at n_ph = 10 (nearest-first 0.030) and 4.2e-5–2.0e-4 at n_ph = 30 (nearest-first 1.2e-3) (`job54_loc_post.tables`, key `spread_which_pulsars_fisher_f10_snr10`, every 4th source). Far pulsars have finer fringes, so Fisher prefers them. In measure (b) the ordering matters much less before the collapse (random/nearest area ratio 0.85–1.0 for n_ph <= 20 at SNR 10) and by 5x after it (n_ph = 80: 1.9e-5 against 1.05e-4 deg^2). Only one random ordering was run in (b).
- (b) used 20 of the 200 sources (indices 0–19), the same 20 for every frequency, SNR and n_ph. Partial knowledge was run at 10 nHz, SNR 10, nearest-first only.
- Island counts depend on the 90% threshold and on connectivity at the grid scale. They are an order-of-magnitude statement. Random-order island counts are not resolved.
- The n_ph = 30 values for (b) are interpolated, not run.
- The explanation of the collapse (the truth peak must outweigh ~10^5–10^7 fringe cells, which needs enough coherent SNR) is an interpretation consistent with the SNR 20 runs. It was not tested separately.

## Files
- Scripts: `job54_loc_core.py` (signal, array, closed-form likelihood), `job54_loc_test.py`, `job54_loc_fisher.py`, `job54_loc_sky.py`, `job54_loc_post.py`.
- `talk_data/localisation.csv`: measures `a` (Fisher, unphased phases free), `a_cond` (same, other parameters fixed), `a_census` (literal census Gaussian priors; do not use for small n_ph), `b`.
- `talk_data/skymaps.npz`: source index 4 (10 nHz, SNR 10). Rule: among the 20 sources, the one minimising the sum over the 11 n_ph of |log10 A90_b - log10 median A90_b| (nearest-first); score 0.77 dex summed, next best 1.05. Contents: expected log-likelihood minus its value at the truth; each 256x256 display cell is the MAXIMUM of the fine evaluation over that cell; float16; all 11 n_ph on one grid; axes are local lon/lat offsets (deg) from the truth, with RA/Dec arrays per cell; true RA/Dec; 90% areas, island counts and truth-island areas from the fine computation; zoomed maps (30x finer, block maximum of an 8x finer evaluation) for n_ph = 10, 20, 40, 80, all. 2.0 MB. See the `readme` string inside.
- Figures: `job54_fig1_area_vs_nph.png`, `job54_fig2_area_vs_sigmaL.png`.
