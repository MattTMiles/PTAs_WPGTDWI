# GALLERY — collaborator-update curation (everything since the last update)

*Read-only curation by agent GALLERY, 2026-07-08. Numbers verbatim from
`project_progress.md` / `trackB_estimator_spec.md` / the agent docs / cited npz.
Figure paths confirmed on disk; results with no figure flagged as such. The
previous update covered only the early prong-2 anatomy (the 2×3 confusion figure,
the A.2 coherence two-panel, the A.3 heatmap, the P2-A detect-vs-range plot).
Everything below is new since then.*

Currency note: "scaled" = source-parameter offset normalised by prior width. The
one number that ties the arc together is the **registration wall ≈ 27×**: the
loosest fringe-registration lane (1.85e-3 scaled sky) sits 27× below the best blind
source localisation (0.05 scaled). Certification becomes feasible when a lane is
raised above the float, or the float lowered onto a lane — the whole update is a map
of what moves those two numbers and when.

---

## 1. ANCHOR CENSUS (A0–A3) — do prior-alone "anchor" pulsars exist on the real array?

**Headline numbers.**
- **Strict anchor set (K ≤ 1, certifiable by the EM prior alone) = EMPTY** — zero
  pulsars under every prior, including hand-injected published composites. Only **one
  pulsar reaches K ≤ 3 at all**: J0437-4715, K = 3.07 at the Reardon+2016 composite
  (0.25 pc). K ≤ 10 list empty; smallest K from the canonical script = J0437 K_opt 11.88.
- **The 3 certifying pulsars** (honest Bayesian bar, realistic 3-loud+13-faint
  population, N_CW=16): **J0711-6830 P(true)=0.953, J1713+0747 P(true)=0.989,
  J1909-3744 P(true)=0.998** — all **DATA-DRIVEN** by the loud sources and
  **prior-independent** (they certify even under the feather prior).
- **Anchor-paper verdict:** the prior-alone-anchor assumption of arXiv:2603.28897
  **FAILS on the real array** — no pulsar's EM prior can phase-lock a fringe without
  help from the CW itself. Coverage: 69/116 parallax-class priors built (A0).

**Figure(s):** **NONE — flag.** All quantitative product is in npz only:
`anchor_a0_priors.npz` (A0 priors), `anchor_a1_Ktable.npz` (dual-prior K table),
`anchor_a2_results.npz` + `anchor_a2_diag.npz` (per-draw cause-tag audit). A
K-table / P(true)-onset bar figure would carry this block and does not exist.

**One sentence:** on today's array no pulsar's distance prior is tight enough to fix
a fringe alone; the only certifications come from loud sources, so the anchor list is
a *loud-source* list, not a *prior* list.

**Doc section:** `project_progress.md` §"Anchor Census — real EM priors vs the feather
priors (COMPLETE, 2026-07-03)".

---

## 1.5 P2-B / P2-C / P2-D — prong-2 close-out (post-previous-update)

**P2-B — coherence-axis grounding.**
- *Required-df/f threshold:* source-side stochastic df/f crosses the distance-info knee
  at **3.4e-5 (SNR=5), 8.4e-6 (SNR=20), 1.7e-6 (SNR=100)** over kyr pulsar-term lags;
  the real environmental magnitude is left explicitly OPEN → handoff (Taylor/Farr).
- *Readable sub-array (phase-floor rank):* of 116 pulsars, **30 can read their own
  pulsar term at realistic strain** (SNR_pterm > 1, i.e. σ_φ below the knee at 1 rad);
  array median SNR_pterm = 0.47, median σ_φ = 2.13 rad. Below-knee leaders:
  J0613-0200 (3.04), J0751+1807 (2.53), J1045-4509 (2.42), J1012+5307 (2.35),
  J0605+3757 / J0931-1902 / J1024-0719 (~2.1). Named chain: J1713+0747 σ_φ 1.26 /
  SNR_pterm 0.79; J0437-4715 0.74 / 1.35; J1909-3744 7.64 / 0.13.
- **Figure:** `CW_transition/p2b_coherence_grounding.png` (2 panels: SNR·σ_φ vs df/f
  knee-crossings; 116-pulsar σ_φ histogram). *One sentence:* red/DM noise is just the
  1/SNR measurement floor already inside the Fisher; only ~a quarter of the array can
  hear its own pulsar term, and the one genuinely new coherence knob is source-side df/f.

**P2-C — array-boost scaling law.**
- **knee/N\* = 0.40 · N_psr^1.03** (fit A = 0.3990, b = 1.0263); N_psr=116 → 52.75
  (reproduces the Stage-A ~52, gate PASS); **linear, NO saturation through 200 pulsars**.
  Forecasting form: N_knee ≈ 0.40 · N_psr · T · Δf.
- **Figure:** `CW_transition/p2c_array_boost.png`. *One sentence:* each doubling of the
  array pushes the confusion wall out proportionally — more pulsars buy resolvable-source
  count with no diminishing returns yet.

**P2-D — the two reconciliation numbers (verbatim).**
- *Item 1 (strain):* per-source optimal SNR (median/80 sources) — D4 equal-strain
  h=−13.75 = **10.7**; nb-03 loud h=−13 = 59.9; old prong-1 optimizer h=−12 = **599.5**.
  **Optimizer/D4 = 56.2× (= 10^1.75 exactly);** nb-03/D4 = 5.6×. The old "distances
  recoverable" successes ran ~56× louder than the realistic regime.
- *Item 2 (nb-01 verdict inversion):* single CW, h=−13, buggy peak-finder joint
  dlnL = **+296.09** ("truth preferred → separable") → corrected scanner **−108.35**
  (wrong fringe beats truth → NOT separable). **The bug flipped the sign of the verdict.**
- **Figure:** NONE — flag (P2-D is npz only: `stagec_p2d_item1.npz`,
  `stagec_p2d_item2.npz`). No figure expected, but noted so the +296→−108 flip and the
  56× are quoted as text, not illustrated.

**Doc section:** `project_progress.md` §"Prong-2 close-out (P2-A..D, 2026-07-02)".

---

## 2. TRACK A (F0–F2) — anisotropic-covariance distance channel

**Headline numbers.**
- Scaling exponent **α = 1.17 [1.09, 1.27]** (I_LL_aniso ∝ N^−α; shot-noise anisotropy,
  c_lm ∝ 1/√N). Slower than the resolved channel's N^−2 past the knee — a *formal*
  rate-crossover.
- **σ_L (within-fringe, best pulsar J0437-4715): 77.5 pc at N=300, 96.1 pc at N=1000**
  (iso baseline ~98.7 pc) — versus the resolved channel's ~2 mpc, **4–5 orders worse**.
- **Verdict: OUTCOME (ii) — channel real but formal; Farr's null effectively holds at
  realistic N.** By N≈1000–3000 the enhancement is 2–6% of the iso baseline and its
  16–84% band already includes zero; the channel inherits the same fringe-ID problem the
  Anchor Census showed fails.

**Figure:** `CW_transition/trackA_f2_scaling.png` (resolved channel vs covariance channel
on shared axes; the meaningful comparison is the decay slopes and physical σ_L, not the
absolute Fisher heights). *One sentence:* sky structure does carry distance information
with a clean, slowly-decaying exponent, but at realistic confusion N its σ_L is
astronomically worse than resolving sources, so anisotropy does not rescue the
confusion-regime route.

**Doc section:** `project_progress.md` §"Track A — anisotropic-confusion distance
information in the cross-pulsar covariance (F0–F2, 2026-07-03)".

---

## 3. TRACK B ARC (in sequence)

**B0.2 — the precision cliff.** Certification tolerance on the source params ≈ **1e-4
scaled (~0.006° sky)**, set by the pulsar-term lever arm 2πL/dL ≈ 1.6e4 rad; achievable
source precision (seed+EM) ≈ **0.5 scaled — 3–4 orders coarser**. F-stat seeding gate
PASS (3/3 loud recovered at 6.4/6.4/11.9°). The census ceilings are a conditional-on-truth
artifact, not an achievable measurement. *Section:* spec §"B0.2 RESULT".

**P1 — the registration needle.** A **unique joint-lnL registration needle sits exactly
at truth**; half-width ≈ **0.006° sky** (reproduces B0.2 as a needle); the census-3
register **3/3**; the needle is carried by **~18 loud-source-broken pulsars** — not just
3, not all 116. Cusp: −26 nat at one 0.0019° sky step, −126 nat at one 0.25-freq-bin step.
**Figure:** `CW_transition/trackB_p1_registration_needle.png` (6-panel registration map:
coarse / zoom1 / uniqueness-patch + the truth-placed cut_cos / cut_gwphi / freqzoom needle
cuts). *Section:* Anchor-Census STATUS UPDATE 2026-07-06 + spec §"P1 STATUS".

**Plateau / gauge-conspiracy probes.** PROBE 1: the flat plateau is a **gauge conspiracy**
(simultaneous fringe re-registration freedom), not absent physics — anchor the array and
moving loud0's sky costs **−82750 nat at 2°**, QTRUE rising smoothly 4.9→24.45. PROBE 2:
**the loop closes at census-class anchoring (k=6; strict monotone across the gap)** — HARD
gap-drop grows k=0→50, k=3→132, k=6→172, k=24→324 nat, and QTRUE is monotone 6.25→24.45
with **zero** anchors, so the *soft* objective is followable end-to-end.
**Figures:** `CW_transition/trackB_p2_anchorsweep.png` (anchor-fraction sweep, the loop
closing) and `CW_transition/trackB_p2_linescan.png`.

**P2 — continuous-methods floor.** The truth-blind soft-EM M-step converges only to
**~0.01–0.05 scaled — a precision floor 100–500× coarser** than the cusp/certification
tolerance; continuous ascent falls off the cusp onto the plateau. **Figure:**
`CW_transition/trackB_p2_linescan.png`. *Section:* running log 2026-07-07 "Track B P2
pipeline".

**LAMBDA probe + L2b/L2c basin (BANKED, float-independent).** The integer-least-squares
(GPS-RTK) cold start is **NO-GO**. Even with fringe integers fixed **at truth**, the
conditional source re-solve pull-in radius is **< 1e-4 scaled**; truth is a genuine sharp
neg-def max, Hessian eigenvalues **[−5.965e11, −14.43]**, condition number **~4e10**, cusp
curvature ~6e11. L2c damped-Newton from a 1e-4 start lands 0.046° / −622 nat (never locks).
*Section:* spec §"LAMBDA feasibility probe — RESULT" + §"F1/F2 AMENDMENT".

**F2 — lever-arm ladder.** Per-(pulsar,source) sky registration tolerance: **loosest
1.85e-3, median 3.81e-5, tightest 1.34e-6 scaled.** Pairs with tol ≥ 0.05: **0**;
≥ 1e-2: **0**; ≥ 1e-3: 3; ≥ 1e-4: 44. The whole ladder lives below any achievable blind
float. **Figure: MISSING — flag.** `trackB_F2_ladder.png` does **not** exist on disk; the
ladder spectrum lives in `trackB_F2_ladder.npz` only. This is the load-bearing "the ladder
does not span" result and has no figure. *Section:* spec §"F2 — the lever-arm ladder".

**The wall (final bracketed height).** WALL = repaired-float ceiling / loosest rung =
**0.05..0.21 / 1.85e-3 = 27×..112× = 1.4..2.05 dex, with ZERO rungs in the gap.** The
GPS-RTK wide-lane cascade cannot fix its first rung.

**Terminal verdict (quoted verbatim, incl. the PENDING-R line).**
> **NO-GO for COLD-START POINT ESTIMATION of CW distance certification on this array** —
> earned, bracketed FLOAT-INDEPENDENTLY … The census ceilings (J0711/J1713/J1909, P
> 0.953/0.989/0.998) remain REAL but CONDITIONAL — achievable only given source knowledge
> at the ~1e-3-scaled level THIS ARRAY CANNOT SELF-GENERATE.
>
> **TERMINAL CHARACTERIZATION PENDING deliverable R (posterior referendum):** whether the
> fringe-marginalised POSTERIOR concentrates at the needle (→ inference / prong-3 succeeds
> where optimisation fails; the NO-GO is only about POINT ESTIMATION) or smears across the
> plateau (→ the NO-GO deepens to INFORMATION-THEORETIC). Point estimation is settled
> NO-GO; the transient J1713=0.977 at a 0.05-scaled float is a tantalising hint that the
> POSTERIOR may carry more than the optimiser reaches — R decides it. R is the terminal
> deliverable.

**Doc section:** `CW_transition/trackB_estimator_spec.md` §§"TRACK B CLOSE-OUT" /
"TRACK B TERMINAL VERDICT (2026-07-07)".

---

## 4. FLEET (CPU-only scoping agents)

**PENCIL — when does the wall close from above?**
- Crossing times (baseline white-noise scaling, σ_sky ∝ T^−½): the float reaches the
  **loosest rung (1.85e-3) only at T ≈ 11 kyr (10,960 yr)**; the L2c pull-in (1e-4) at
  **T ≈ 3.75 Myr**. Buy it with strain instead: first rung needs **per-source SNR ≈ 289
  (h ≈ 10^−12.32)**, ~27× the current loudest source; L2c needs SNR ≈ 5,350. Optimistic
  (physically rejected) T³ bound only softens 11 kyr → 135 yr.
- **Figure:** `CW_transition/PENCIL_t_crossing.png` (σ_sky vs T, both rungs, the 15-yr
  anchor and float band, crossing markers). *One sentence:* the wall does close as data
  accumulate, but the first rung is ~11 kyr / SNR≈289 away — off any observational horizon
  for blind self-integration, though an EM host position is worth ~10⁴ yr of it instantly.

**LANES — do eccentric harmonics populate the gap?**
- **Scissors verdict: NO-GO** — does not span the 27× gap, no optimal-e window. Widest
  usable lane peaks at **~8.6e-3 scaled at e=0.9** (power-anchor), physically ~1.85e-3–3.7e-3
  (the F2 rung itself) — 5.8×–27× short of the 0.05 float ceiling. Two independent, each-fatal
  mechanisms: power-anchor geometry (fundamental drops below 5% of peak at high e) and
  residual-SNR reality (per-harmonic SNR² ∝ g(n,e)/n⁴ → SNR-dominant harmonic is always
  n=1–2, so the detection band *is* the fundamental). **Eccentricity buys precision, not
  capture.** Angle is novel/unclaimed in the literature.
- **Figures:** `CW_transition/lanes_ladder.png` (lane rungs vs e over the shaded gap) and
  `CW_transition/lanes_kcontour.png` (K_eff over (e, f_orb) for three Mc). *One sentence:*
  the same 1/n⁴ that could make the fundamental a loud low-registering lane also makes it
  the detection band, so there is no "detect high, register low" lever to climb.

**SCOUT — is there a loud *and* counterpart-identified source?**
- **No figure — flag** (literature census + closed-form arithmetic; the candidate table is
  the product, no plot). Joint counts: **N_joint ≈ 10^−2…10^−1 (current PTAs), O(0.1–1)
  in the SKA era**, set almost entirely by loudness (N̄_detectable ≲ 0.1), not the
  counterpart cut. No named credible SMBHB clears **h = 10^−13.75** at its real mass; 3C 66B
  / PG 1302 / PKS 2131 come within 1 dex only at an excluded or contested loud edge. The
  22-arcsec (0.006°) EM tolerance is met with **≥10⁴× margin** by any VLBI/Gaia counterpart —
  astrometry is never the binding constraint; loudness is. Operating point −13.75 sits ~0.35
  dex *above* the NG15 sky-averaged CW upper limit (−14.1). *One sentence:* the counterpart is
  free; the program lives or dies on the loud tail of the population, where the expectation is
  ≪ 1 source and the value is a variance play, not an expectation-value detection.

**ABACUS — reconciling the measured confusion law with Boyle & Pen (2012).**
- **RECONCILED.** Measured 2/5 (coefficient 0.40) vs Boyle–Pen **2/7** (0.2857); ratio
  **7/5 = 1.40, exact to <1%**, = the removal of Boyle–Pen's 2 source-sky parameters (7→5
  params/source, because the toy fixes the source sky and asks only whether distances
  survive). Residual ~14% in the raw count at 116 psr is soft-vs-hard threshold + finite-F,
  not a coefficient effect. **No figure — flag** (it calibrates the P2-C law;
  `CW_transition/p2c_array_boost.png` is the associated plot, ABACUS itself produced none).
  *One sentence:* our confusion-capacity measurement *calibrates and normalises* the
  Boyle–Pen deconfusion bound rather than contradicting it — same ∝ N_psr·F shape, coefficient
  fixed at 2/5 under the known-source-sky mapping.

**Doc sections:** `CW_transition/PENCIL_t_crossing.md`, `LANES_eccentric_ladder.md`,
`SCOUT_counterpart_census.md`, `ABACUS_boyle_pen.md`.

---

## 4b. SHOVEL — Q1 excavation verdict

**One line:** *Clean, reusable forward-Q1 harness (method healthy — a proper
`2·(logL_max−logL_null)` profile statistic, entirely free of the prominence-bug /
conditional-metric contamination); the specific Dec-2025 banked run's outputs are
SUPERSEDED (uncalibrated "SNR" ~10⁴× too loud, log10_fgw=−10.5 vs the design −8 which
removes the very phase-evolution mechanism, resolution-floored constant area90, n_real=1
scenario-A only) → **resume-and-modernize on the code, the banked figures do not
graduate**.* No figure. *Doc:* `MAIN_PROJECT_QUESTIONS/SHOVEL_q1_report.md`.

---

## 4c. THE DESIGN THEOREM (quoted as written)

> **Certification requires registration LANES below ~2e-3 scaled that a blind float can
> reach, OR a blind float above the loosest existing lane. tol ~ 1 / (2πf L_p (1−cos μ_ps)/c)
> per sky axis. Three constructive levers:**
> **(i) WIDE LANES FROM NEARBY PULSARS (measured).** tol ~ 1/L_p …
> **(ii) WIDE LANES FROM ECCENTRIC HARMONICS (E-track's sharpened job).** …
> **(iii) LOWERING THE WALL FROM ABOVE.** Better Earth-term localisation shrinks the float
> floor toward the lanes …
> Certification is feasible where lever (i)/(ii) raises a lane above the lever-(iii)-lowered
> float.

**Two tombstones (dead levers, with numbers).**
- **PENCIL (lever iii, blind integration):** first rung at **T ≈ 11 kyr / per-source
  SNR ≈ 289** — off any horizon.
- **LANES (lever ii, eccentric harmonics):** **scissors** — the 1/n⁴ that would make the
  fundamental a loud low-registering lane also makes it the detection band; widest usable
  lane ~8.6e-3, never reaches 0.05.

**Two survivors.**
- **Lever (i): the 1/L nearby-pulsar law.** The 3 loosest F2 rungs are the nearest pulsars —
  **J0711-6830 (0.106 kpc, tol 1.85e-3), J1630+3734 (0.089 kpc, 1.39e-3), J0437-4715
  (0.155 kpc, 1.02e-3)** — vs the median L=1.38 kpc; ~13× nearer → ~13× looser lane
  (J0437 holds the loosest rungs; already a census / anchor pulsar). A dedicated ~100 pc,
  well-timed sub-array raises the top lanes ~10× toward the float ceiling.
- **Lever (iii): the targeted (EM-counterpart) bypass.** A host position removes the sky axis
  entirely and clears the **22-arcsec (0.006°) tolerance with ≥10⁴× margin** (SCOUT); the
  gating clock is then SCOUT's **N_joint ≈ 10^−2…10^−1 now, O(0.1–1) SKA-era** — a loudness
  clock, not a localisation one.

**Doc section:** `CW_transition/trackB_estimator_spec.md` §"THE DESIGN THEOREM".

---

## 5. R — the posterior referendum (STILL RUNNING)

**Status: incomplete — R1 done, R2 in progress, R3/R4 not yet produced.** All foundation
gates A–E PASS (split == fisher_logL bit-exact; census exp(−m_p) reproduces the ceilings
0.953/0.989/0.998; additive/star-topology approximation error bounded ≤ 0.14 nat).

**Delivered so far (R1, banked in `trackB_R_znaddle.npz`):**
- **ln Z_needle = 405617.637 (quadrature) / 405617.844 (Laplace)**; needle binding σ (sharpest
  eig) = 2.886e-6 scaled; lnL_marg(truth) = 405686.34 = lnL_ref 405413.51 + Σ_p m_p 272.83.

**What it will still deliver (R2 → R3 → R4):**
- **R2** — ln Z_plateau by tempered SMC over the truth-centred ±2° box (needle excised),
  throughput ~570 ms/particle; live checkpoint `trackB_R_zplateau_s0_ckpt.npz`; will write
  `trackB_R_zplateau_summary.npz`.
- **R3** — the referendum number **f = Z_needle / (Z_needle + Z_plateau)**, its 1-σ band, the
  **break-even sky-prior radius θ\*** (deg per loud source at which Z_needle = Z_plateau) and
  the **break-even volume factor** V\*/V_2deg; writes `trackB_R_referendum_result.npz`.
- **R4** — pre-registered verdict branch: **f ≥ 0.95** → the data identify truth, Track B's
  failure is SEARCH-ONLY, guided/anchored inference (prong 3) unparks; **f ≤ 0.05** → NO-GO
  deepens to information-theoretic (only more data helps); **0.05 < f < 0.95** → intermediate,
  report f + the break-even sky radius, and B1 (noise) becomes decisive.

**No figure produced or expected.** Checkpoints/results all in npz under `CW_transition/`
(`trackB_R_znaddle.npz`, `trackB_R_znaddle_hess.npz`, `trackB_R_zplateau_s0_ckpt.npz`).
*Doc:* `trackB_R_referendum.py` header + running-log terminal-verdict PENDING-R line.

---

## 6. QUEUE AS PLANNED (each phrased as a forecast — what capability it prices/exploits)

- **B1 (targeted, noisy, scrambled null)** — prices *what noise costs the CONDITIONAL
  (source-localised) pipeline*, not cold-start certification: ≥20 realisations → certified
  set, P(true) distributions for the 3 ceiling pulsars, wrong-certification count, SBC
  reliability curve, BH-FDR@0.05. Exploits a supplied sky (counterpart / anchor regime); it
  is the price tag on certification *once the wall is bypassed*.
- **B1.5 forecast maps with prior tiers** — prices *how much better distance priors move the
  certification onset*: sweeps prior tiers (feather → EM composite → future tight) to map the
  certified-count surface, i.e. the payoff of each dex of prior improvement.
- **E-track surviving channel (timestamp / evolution-envelope)** — prices *eccentric
  harmonics for precision, not capture*: the surviving lever-(ii) job is to measure whether
  the timestamp/envelope (Peters-evolution) channel adds finer rungs at the L2c end, plus the
  two flagged escape hatches (an anomalously aligned 1−cos μ→0 pulsar with an unusually loose
  fundamental; a strain-limited non-residual regime).
- **Q1-modernized** — prices *the inverse computation done right*: real matched-filter SNR
  targeting + log10_fgw≈−8 + realistic TOA errors + off-grid localisation + N realisations +
  scenarios B/C, on the SHOVEL-certified harness. Buys dex on the float side (lever iii).
- **Inversion track (source localisation FROM certified distances)** — prices *the first loud
  detection turned into sky precision*: each certified distance lowers the blind-float floor
  toward the lanes, the direct lever-(iii) mechanism.
- **Payoff chain** — prices *the standard-siren endgame*: **localization → kyr-baseline chirp
  → chirp mass ℳ_c → luminosity distance D_L → standard sirens.** Named pain point to beat:
  the **29%-within-1-dex ℳ_c-recovery** rate (i.e. current chirp-mass recovery lands within a
  decade of truth only ~29% of the time) — the chain's weakest link, and the quantity a
  kyr-baseline localised chirp is meant to sharpen. *(This 29% figure is Matt's standing
  Mc-recovery number; it is not stated in the Track-B / prong-2 doc set curated here — flag
  its provenance for the letter.)*
- **Scintillation patch** — prices *a new tight-prior distance channel*: whether scintillation
  distances can add prior-alone anchors the parallax catalogue cannot (the Anchor Census found
  none), i.e. lever (i)/(iii) via a different distance estimator.
- **The A/B loop (one cycle):** **detection → host-ID** and **host → certification** are a
  single closed cycle — the first loud detection unlocks the counterpart catalog that supplies
  the sky that bypasses the wall that certifies the distance that (via the payoff chain)
  localises the next source. Every queue item above is one arc of this loop.

---

### Figure-existence ledger (confirmed on disk 2026-07-08)

| result | figure | on disk? |
|---|---|---|
| Anchor Census (A0–A3) | — | **NONE (flag)** |
| P2-B coherence grounding | `p2b_coherence_grounding.png` | yes |
| P2-C array boost | `p2c_array_boost.png` | yes |
| P2-D loose ends | — | **NONE (flag; npz only)** |
| Track A F2 scaling | `trackA_f2_scaling.png` | yes |
| Track B P1 needle | `trackB_p1_registration_needle.png` | yes |
| Track B P2 anchor sweep | `trackB_p2_anchorsweep.png` | yes |
| Track B P2 line scan | `trackB_p2_linescan.png` | yes |
| **Track B F2 ladder / the wall** | `trackB_F2_ladder.png` | **MISSING (flag; npz only)** |
| PENCIL crossing | `PENCIL_t_crossing.png` | yes |
| LANES ladder | `lanes_ladder.png` | yes |
| LANES K_eff contour | `lanes_kcontour.png` | yes |
| SCOUT census | — | **NONE (census table)** |
| ABACUS reconciliation | — | **NONE (calibrates p2c fig)** |
| R referendum | — | none yet (running; no figure expected) |

*All paths relative to `CW_transition/`.*
