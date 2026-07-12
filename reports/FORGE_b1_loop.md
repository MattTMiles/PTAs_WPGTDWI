# FORGE_b1_loop.md — agent FORGE, ACCRE (tag b1-v1)

*Untracked working report. PURE EXECUTION: no tracked file was edited, nothing committed or
pushed. Cronus is canonical; this file never feeds back into `project_progress.md`. All numbers
are from the A100 run recorded in §8.*

---

## 0. WHAT B1 IS (re-labelled per the Track B terminal verdict)

Blind AND targeted-circular **cold-start** certification are **CLOSED** — information-theoretic
NO-GO (`trackB_estimator_spec.md`: referendum f = 6.9e-7; the targeted (f, mc) wall; the
soft-cascade FAIL; the census "3 certified" pulsars are **CONDITIONAL** ceilings, achievable only
given source knowledge at the ~1e-3-scaled level this 15-yr array cannot self-generate). B1 does
**not** re-open them.

B1 measures the **CONDITIONAL pipeline under real noise** — conditional on the source knowledge a
self-clocking detection would supply (sky + f + mc known to certification grade) — and runs Hogg's
phase-up loop to its fixed point. **Every count below is a CONDITIONAL CEILING**, not an achievable
cold-start measurement.

The one genuinely open question (B1.3): does the loop's **source-fit channel** — certified pulsar
terms are kyr chirp measurements feeding Mc — **COMPOUND**? The static soft-cascade (STEP 1D)
FAILED, but it re-weighted a fixed fringe mixture and lacked this channel. B1.3 wires it: fix the
certified fringes → re-fit (f, mc) on the noisy data (smooth; fringes fixed; no L2c cusp) → shrink
σ(log10_mc) → open the registration gate R_a≥1 for the next rung → re-certify. **Gain > 1
sustained** = the loop hot-wires its own clock and the circular wall falls to bootstrapping.
**Damped** = E-track ignition is the only door, doubly confirmed.

Machinery (all banked + gated): `trackB_b1_core` (masked E-step, real 3-component noise draw,
arm-B truth draw, fringe posterior), `trackB_b1_ladder` (registration radius + conditional (f, mc)
Fisher), GEO sky-conditional seed sets. Driver `hpc_harbor/forge/forge_b1.py`; consolidation
`hpc_harbor/forge/forge_consolidate.py`; figures `FORGE_results/forge_figures.png`.

**Pre-registered value gate (FORGE gate, job 12491495).** The geometry-swap E-step reproduces the
banked census/GEO ceilings at zero noise to |Δ| ≤ 4e-4 (J0711 0.9534, J1713 0.9887, J1909 0.9984),
and geometry enters only through runtime θ (one build serves every sky). Scrambled-null and both
arms run finite. **ALL PASS.**

---

## 1. B1.0 HEADLINE — CONDITIONAL CEILINGS UNDER REAL NOISE (the ballgame)

30 realisations per arm = **10 sky draws (fiducial + 9 GEO redraws) × 3 noise weathers**; all three
noise components DRAWN (heterogeneous white + per-pulsar RN GP 30-comp + HD-correlated GWB GP
14-comp). Certification = q_max > 0.9 (and the strict > 0.99), B_CERT = 512. **CONDITIONAL CEILING:
the loud sources' sky + f + mc are held at certification grade; only the fringe integers are
inferred.**

| arm | certified count P>0.9 | P>0.99 | WRONG-cert (P>0.9) | WRONG-cert (P>0.99) |
|---|---|---|---|---|
| **A** — L_true = L0 (continuity vs census) | **2.87 ± 1.48** (med 3, range 1–7) | 0.87 ± 0.50 | 2 / 30 real (0.067/real) | 0 / 30 |
| **B** — n_true drawn (**HEADLINE**) | **1.43 ± 1.05** (med 1, range 0–4) | 0.37 ± 0.48 | **8 / 30 real (0.267/real)** | **3 / 30** |

Noise-weather statistic (lnL drop truth|clean → truth|noisy): median **15094.7**, range 14842–15133.

> **Read these counts WITH the §3 scrambled-null floor.** The Bayesian q_max>0.9 criterion carries a
> source-independent floor of ≈1 (prior-pinned) rising to ≈3 (noise-lock) under real noise. The
> absolute counts here are Bayesian ceilings *inclusive of that floor*; the A→B price and the
> per-pulsar prior-pinned-vs-data-driven split below are relative and unaffected by it.

**THE A→B PRICE — the headline.** Arm A under real noise reproduces the census count (**2.9 ≈ 3**);
noise alone barely dents the ceiling when truth sits at the prior mean. **Arm B — truth drawn off
the prior mean (the honest case) — HALVES the certified count (2.87 → 1.43) and QUADRUPLES the
wrong-certification rate (0.067 → 0.267 per realisation, and 3 confident P>0.99 wrong certs vs 0).**
Registration-from-the-prior-mean was worth a factor ~2 in yield and was suppressing essentially all
confident wrong certs. This is the first Track B number for truth *not* at the prior mean, under
real noise.

**PER-PULSAR, Arm B — PRIOR-PINNED survives, DATA-DRIVEN evaporates.** Median P(true fringe) and
certify-frequency across the 30 Arm-B realisations:

| pulsar | median P(true) | certify-freq (P>0.9) | on the TRUE fringe when certified | class |
|---|---|---|---|---|
| **J0437-4715** (anchor, K≤3) | **0.731** | 0.43 | **13/13 (100%)** | PRIOR-PINNED |
| J1909-3744 (census) | 0.156 | 0.47 | (supplies wrong certs) | data-driven |
| J1713+0747 (census) | 0.000 | 0.10 | (supplies wrong certs) | data-driven |
| J0711-6830 (census) | 0.003 | 0.07 | (supplies wrong certs) | data-driven |

The Arm-B wrong-certifications come **entirely from the data-driven census/loud-broken pulsars**
(J1909 ×3, J1713 ×3, J0711 ×1, J1603 ×1) — never from the anchor. **J0437-4715 certifies 13/30 and
is on the TRUE fringe every single time.** This is `trackB_estimator_spec.md`'s R-POSTMORTEM made
concrete under real noise: prior-pinned fringe information (tight composite EM prior → K≤4 candidate
fringes) is robust to truth-off-prior-mean and to noise; data-driven fringe information
(loud-source-broken) is not — it re-registers onto a wrong fringe once truth leaves L0. **The one
name that survives a realistic sky, noise, AND an off-prior-mean truth is J0437, the sole
Anchor-Census K≤3 pulsar — not any census-triple name.**

---

## 2. B1.3-ITERATED — THE HOGG PHASE-UP LOOP. VERDICT: (iii) DAMPED, doubly confirmed.

12 Arm-B realisations (fiducial + GEO redraws), the fixed-point loop: data-driven E-step seed →
conditional (f, mc) Fisher fit on the certified set (fringes fixed, noisy data) → registration gate
R_a≥1 → re-certify, ≤10 iterations. Primary axis = **N_cert growth** (does the certified set
compound?); σ(log10_mc) is the diagnostic.

| quantity | result |
|---|---|
| realisations that GREW past their round-0 seed set | **0 / 12** |
| median N_cert: seed → fixed point | **1 → 1** (every realisation converged in a single round) |
| pooled next-cycle gain > 1 | **0.00** (no second round of new certs anywhere) |
| σ(log10_mc): conditioning float → round-0 fit → within-loop shrink | ~0.025 dex → 0.039 dex (median) → **1.00× (saturates immediately)** |
| σ(log10_mc) below the 0.003-dex self-clock bound | 5 / 12 (locally, from the seed certs) |
| wrong-certification at the fixed point | median 1, **sum 12 / 12 realisations** |

**>>> B1.3 VERDICT: (iii) DAMPED. The source-fit channel does NOT compound.** N_cert never grows
past the data-driven seed set; the registration gate admits zero new pulsars in every one of the 12
skies. Two mechanisms, both measured:

1. **Local σ_mc is not the bottleneck** (STEP-1D, now confirmed *with* the fit channel the cascade
   lacked). The fit does tighten σ(log10_mc) — the census-sky seed reaches 1.4e-4 dex, well below
   the 0.003-dex bound — but that tight *local* width never becomes *global fringe concentration*.
   The uncertified pulsars' E-step q_max stays < 0.9 under real noise: the missing information is
   *which* of ~16 fringes, and the chirp-mass channel does not sharpen that choice. σ_mc saturates
   at the round-0 value (within-loop shrink 1.00×).
2. **Under real noise + Arm B the SEED SET is wrong-fringe-poisoned.** The census-sky realisation
   seeds {J1603, J1713, J1909} — **all three on the WRONG fringe** — so the source-fit channel
   bootstraps from FALSE fringes (biased Mc), and wrong-cert persists at every fixed point (12/12).
   Noise defeats the Asimov bistability (STEP-1B's gain 4.50 given census-3) precisely by
   contaminating the seed the loop would have compounded.

**The loop cannot hot-wire its own clock. E-track ignition (eccentric-harmonic likelihood structure,
κ ≥ the STEP-2C >20× bound) is the only door — doubly confirmed:** once statically (STEP-1D
soft-cascade) and now dynamically, under real noise, with the source-fit channel wired in.

Trajectory figure (N_cert and σ_mc vs iteration, all realisations overlaid):
`FORGE_results/forge_figures.png`.

---

## 3. B1.2 SCRAMBLED-NULL LINE — the null does NOT certify zero; anatomy resolves why

Pre-registered rule: *any confident cert in the null = STOP + anatomy.* The null **fired**, so here
is the anatomy (job 12491794, 5 realisations per variant, Arm-B distances):

| null variant | n_cert P>0.9 (per real) | total | P>0.99 total |
|---|---|---|---|
| **nullL** — 3 loud sources scrambled, 13 faint at truth | [1, 2, 2, 4, 2] → **2.2/real** | 11 | 4 |
| **nullA** — ALL 16 sources scrambled (honest Cornish-Sampson) | [2, 2, 2, 7, 1] → **2.8/real** | 14 | 4 |
| **nullN** — NO CW in the data, recovery at the TRUE source | [1, 2, 0, 0, 1] → **0.8/real** | 4 | **2** |

**My first hypothesis — that the loud-only scramble left the faint sources coherent — is REFUTED:
scrambling all 16 does not reduce the count (2.8 ≈ 2.2).** The floor is **intrinsic to the Bayesian
q_max > 0.9 criterion under real noise**, and it has two measured components:

1. **A prior-pinned floor (≈0.8/real, from `nullN`, pure noise + no CW).** Pulsars with a tight EM
   prior (K ≤ a few candidate fringes — J0437-class) concentrate q_max > 0.9 on the MAP fringe from
   the **prior tail alone**, independent of the source or even of any signal. **This is the SAME
   mechanism that makes J0437 the robust certifier in §1 and §5** — prior-pinning is a double-edged
   sword: robust to source errors, but it also fires in the null. Two of these are "confident"
   (P>0.99) from pure noise.
2. **A noise-lock excess (0.8 → 2.8/real) when a wrong source model meets real CW data.** Real-noise
   per-fringe likelihood peaks give a mis-placed source spurious fringe structure to lock onto,
   inflating the floor above the pure-prior level.

**SCRAMBLED-NULL LINE:** *the raw Bayesian certified count is NOT a clean detection statistic under
real noise — it carries a source-independent floor of ≈1 (prior-pinned) rising to ≈3 (scrambled
source + real CW). Separating source-driven from prior/noise-driven certifications requires the
DATA-DRIVEN flat criterion `dlnL_a > ln K_a` (the likelihood gap must beat the trials factor —
GEO's "flat count"), which this run's lean npz did not store.* This does **not** overturn the §1
A→B **price** or the §1/§5 **prior-pinned-vs-data-driven** result (both are relative / per-pulsar),
but the **absolute** B1.0 Bayesian counts (Arm A 2.87, Arm B 1.43) sit atop this floor and must be
quoted as Bayesian ceilings *inclusive of it*. Pre-registered response delivered: the STOP flag on
the raw Bayesian criterion is honoured with the mechanism; the disciplined follow-up is a re-cut on
the flat criterion (store `dlnL_a`, `ln K_a` per pulsar; one more array). **Flagged, not run.**

---

## 4. B1.1 CALIBRATION — the pipeline is CALIBRATED (not miscalibrated)

Reliability curve pooling all (pulsar, realisation) claims across the B1.0 ensemble — claimed q_max
band vs realized fraction landing on the TRUE fringe:

| claimed q_max band | realized fraction on true fringe | n |
|---|---|---|
| [0.50, 0.70) | 0.506 | 259 |
| [0.70, 0.80) | 0.675 | 83 |
| [0.80, 0.90) | 0.814 | 59 |
| [0.90, 0.95) | 0.884 | 43 |
| [0.95, 0.99) | 0.959 | 49 |
| [0.99, 1.00) | 0.919 | 37 |

The claimed q_max **tracks the realized true-fringe fraction across the whole range** (0.51 → 0.96)
**when a signal is present**, with only mild over-confidence in the topmost bin (0.919 realized at
claimed ≥0.99). **BH-FDR@0.05** over the pooled array-level claims: 12 claims pass, realized true
fraction among them = **1.000** (target ≥0.95). So the **per-claim posterior is calibrated on real
(signal-present) data** and the wrong-certs of §1 are its honest tail — *however*, §3 shows the same
q_max carries a source-independent false-alarm floor in the **null** (prior-pinning + noise-lock),
which the on-signal reliability curve cannot see. Pre-registered verdict axis: **(i) survives
calibrated per-claim on real data**, with Arm B **(ii) degrading honestly** vs Arm A; the null-floor
is a **STOP-flag on the absolute-count interpretation** (§3), resolved by the flat criterion, not a
per-claim miscalibration.

---

## 5. B1.4 WRONG-POSITION — PASS with mechanism

2 realisations, counterpart offset by **5× the certification tolerance**. Both certify **only
J0437-4715, and on the TRUE fringe.** Every loud-source-*dependent* (data-driven) certification
vanishes when the counterpart is wrong — the pipeline correctly fails to certify anything that
leans on the mis-placed source. The lone survivor is the **prior-pinned anchor**, which certifies
from its own tight EM prior *independent of the source position* — a correct certification, not a
false alarm. **B1.4 fails loud exactly where it should (data-driven certs → 0) and stays correct
exactly where it should (prior-pinned anchor → true fringe).** Same PRIOR-PINNED vs DATA-DRIVEN
split as §1.

---

## 6. B1.5-lite — CERTIFIED COUNT vs (sky draw × noise weather) [Arm B]

| sky | n001 | n002 | n003 |
|---|---|---|---|
| fiducial | 4 | 1 | 1 |
| g00 | 1 | 1 | 1 |
| g01 | 1 | 2 | 3 |
| g02 | 1 | 1 | 0 |
| g03 | 3 | 4 | 3 |
| g04 | 2 | 2 | 1 |
| g05 | 1 | 0 | 1 |
| g06 | 0 | 1 | 1 |
| g07 | 2 | 2 | 1 |
| g08 | 1 | 0 | 1 |

Empirical spread **1.43 ± 1.05, range 0–4**. The **sky draw dominates** the variance (g03 gives 3–4
across all weathers; g05/g06/g08 give 0–1) — geometry, not noise weather, sets the Arm-B yield,
consistent with GEO's geometry-driven selection function. The forecast maps will refine this
surface.

---

## 7. WHAT B1 SETTLES

1. **The conditional ceiling is real but modest and geometry/arm-dependent.** Under real noise with
   truth at the prior mean (Arm A) the census count reproduces (2.9 ≈ 3); with truth off the prior
   mean (Arm B, honest) it halves to 1.4 and grows a confident wrong-cert tail. Loudly: **CONDITIONAL
   CEILINGS**, and — per §3 — Bayesian ceilings that sit atop a ≈1–3 noise/prior false-alarm floor
   the flat criterion must subtract. Not achievable cold-start counts.
2. **Prior-pinned beats data-driven under every stress.** J0437-4715 (K≤3) is the single robust
   certifier across sky, noise, off-prior-mean truth, and a wrong counterpart; the census triple
   supplies the wrong-certs. Any seed set must be **sky-conditional and anchored on J0437 + J1909**,
   never the fixed census triple (0/40 GEO reproduction, and here the triple is the wrong-cert
   source).
3. **The phase-up loop does not self-clock (B1.3 DAMPED).** With the source-fit channel wired in and
   real noise, N_cert never compounds: local σ(log10_mc) tightens but never converts to global
   fringe concentration, and the seed set is wrong-fringe-poisoned under Arm-B noise. **E-track
   ignition is the only door, doubly confirmed.**
4. **The pipeline is calibrated per-claim on signal-present data** (B1.1 reliability + BH-FDR@0.05
   realized-true 1.000), so the wrong-certs are the honest tail of a calibrated posterior — **but**
   the scrambled null (§3) exposes a source-independent false-alarm floor (≈1–3) from prior-pinning +
   noise-lock that the raw Bayesian q_max>0.9 count cannot distinguish from a detection. The
   pre-registered STOP is honoured as an anatomy (§3): the DATA-DRIVEN flat criterion `dlnL > ln K`
   is required for an absolute count; the per-claim probabilities themselves are sound.

---

## 8. EXECUTION RECORD

| item | value |
|---|---|
| lane | `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1` |
| gate (v1, walltime kill) | 12491263 — FG1–3 PASS, FG4 hit the 60-min wall on a DOUBLE hvp cold-compile |
| gate (v2, fixed) | **12491495** — ALL PASS, 948 s; ceilings reproduce to 4e-4, single hvp+grad compile |
| production run | **12491565** — exit 0, 676 s; 74 cert + 12 iter realisations, one build amortised |
| null anatomy | **12491794** — nullL 2.2/real, nullA 2.8/real, nullN 0.8/real (§3); exit 0, 354 s |
| device | dgx03, A100-SXM4-80GB (drew clean / see logs); autotune off, x64, warm persistent cache |
| disk | 86 lean npz ≈ 0.4 MB total in `FORGE_results/` — 6+ orders below the 3.1 TiB group headroom |
| fix that mattered | the conditional (f,mc) HVP (`jit_hvp`, forward-over-reverse of the full logL) is a slow cold compile; it was being compiled TWICE (once wasted in a lazy grad rebuild) — collapsed to one, and the gate/run then fit walltime comfortably |

**Nothing was committed. Nothing was pushed. No tracked file was edited.** `FORGE_b1_loop.md` and
`FORGE_results/` are untracked; `hpc_harbor/forge/` is untracked.

---

## 9. TWO-LAYER RE-SCORE (RESCORE follow-up) — DETECTION `dlnL>ln K` ⊕ CERTIFICATION `P_true>0.9`

The §3 anatomy prescribed separating source-driven from prior/noise-driven certs with the flat
DETECTION gate. Done here. **No new realisations:** the banked (seed → data → certs) are
deterministic, so the B1.0 + null realisations were re-run only to *extract* the likelihood-only
fringe gap `dlnL` (best-minus-runner-up peak, prior-free) that the lean npz had not stored; the
Bayesian certs reproduce identically. GEO's zero-noise ceilings re-score purely from the banked
`geo_summary.npz` (`ident` = `dlnL>ln K`). Two layers: **DETECTION** = `dlnL_a > ln K_counted,a`;
**CERTIFICATION** = `q_max > 0.9` (and 0.99) *within* detections. Re-score job **12496241** (420 s);
tables from `flat_*.npz` (89) + banked GEO. Driver `forge_b1.py rescore`, consolidate
`forge_rescore.py`.

### 9.1 B1.0 under the two-layer gate — the flat gate KILLS every wrong-cert

| arm | Bayesian cert P>0.9 | DETECTIONS (`dlnL>lnK`) | **two-layer cert P>0.9** | **two-layer WRONG-cert** |
|---|---|---|---|---|
| A (L_true=L0) | 2.87 | 0.33 ± 0.54 | **0.33 ± 0.54** (P>0.99: 0.23) | **0 / 30** (was 2) |
| B (n_true drawn) | 1.43 | 0.17 ± 0.45 | **0.13 ± 0.43** (P>0.99: 0.07) | **0 / 30** (was 8) |

**The flat gate removes ~90 % of the Bayesian certs as prior/trials-driven and — decisively —
eliminates EVERY wrong-certification in both arms (2→0, 8→0).** What survives the two-layer gate on
real data is *always on the true fringe*. The honest genuine-detection ceiling under real noise is
tiny: **0.1–0.3 certified per realisation.**

### 9.2 GEO zero-noise ceilings, two-layer (banked, for continuity)

| | Bayesian P>0.9 | DETECTIONS (ident) | two-layer P>0.9 |
|---|---|---|---|
| GEO 40 draws | 4.50 (1–9) | 1.38 | **1.35 ± 0.82 (0–3)**; P>0.99 1.02 |

At zero noise the two-layer ceiling is **1.35/draw** (the Bayesian 4.5 was ~70 % prior/trials-driven).
Per-pulsar detect-freq: **J1909-3744 0.53, J0437-4715 0.65** survive as genuine likelihood detections;
**J0711-6830 0.05, J1713+0747 0.10** collapse (their Bayesian certs were prior/trials-driven). So the
real zero-noise detectors are J1909 (wide prior, large likelihood break) and J0437 (tiny K ⇒ low
trials bar AND a real break) — not the full census triple.

### 9.3 THE NULL LINE — reduced ~85 % but NOT zero → STOP + anatomy (as instructed)

| null bank | Bayesian cert (total) | DETECTIONS | **two-layer cert (total)** |
|---|---|---|---|
| null (loud, 910000, the original B1.2), 10 real | 22 (2.2/real) | 6 | **4** (0.4/real) |
| nullL (3 loud scrambled), 5 real | 11 | 3 | **2** |
| nullA (all 16 scrambled), 5 real | 14 | 2 | **2** |
| nullN (NO CW, true source), 5 real | 4 | 1 | **1** (P>0.99) |

**The two-layer gate cuts the null floor ~85 % (2.2–2.8 → 0.2–0.4 per realisation) but does NOT reach
zero.** Per the pre-registration this is a STOP; the anatomy:

- **The residual null certs are the SMALLEST-K pulsars.** Every residual null cert is J0437-4715
  (**ln K = 1.39**, the array minimum), J1909-3744 (ln K ≈ 4.2), or J1713 (a rare 1-in-580 noise
  fluctuation beating ln K = 6.75). For a pulsar with a tight EM prior the trials bar `ln K` is so
  low (J0437: 1.39 nat) that a **noise-induced likelihood fluctuation's `dlnL` exceeds it**. The very
  tiny-K that makes J0437 the robust "anchor" (§1, §5) is *also* what lets pure noise pass its flat
  gate. Prior-pinning's double edge, sharpened: it survives every gate, including on noise.
- **Consequence — the residual is an irreducible small-K statistical false-alarm, not a modelling
  bug.** `dlnL > ln K` is a per-pulsar likelihood-ratio test; with `ln K ≈ 1.4` its false-alarm
  probability under real noise is not negligible. The honest fix is an **absolute `dlnL` floor** in
  addition to the relative `ln K` (e.g. require `dlnL > max(ln K, ~8 nat)`), which the census/D3
  `ln K + 3` stricter variant already gestures at. Flagged for the estimator spec; **not imposed
  here** (it would be a new criterion, not a re-score).

### 9.4 The number that matters — real-data detections vs the null floor

Arm-B two-layer detections (**0.13/real**) sit **at or below the two-layer null false-alarm floor
(0.2–0.4/real)**, and both are the *same small-K pulsars* (J0437, J1909). **By count, the noisy
Arm-B conditional pipeline does not detect above its own null under a proper detection gate.** The
one thing that still separates signal from null: on **real data every two-layer cert is on the TRUE
fringe (0 wrong-cert)**, whereas null certs land on arbitrary fringes — so the discriminator that
survives real noise is *fringe correctness*, not *count excess*. Zero-noise GEO keeps 1.35 genuine
detections above a ~0 null; **real noise both adds a ~0.3 small-K false-alarm floor and smears the
genuine detections down into it** — noise, not the gate, is what erases the two-layer ceiling.

### 9.5 J0437 COLUMN — does the prior-pinned anchor survive the flat gate?

| | Bayesian cert-freq | DETECT-freq (`dlnL>lnK`) | two-layer cert-freq |
|---|---|---|---|
| GEO zero-noise | 0.80 | 0.65 | **0.62** (survives — genuine break, tiny K) |
| FORGE B1.0 Arm B, real noise | 0.43 | 0.10 | **0.10** (mostly prior-pinned; the 0.10 it does detect are on-true) |

**At zero noise J0437 SURVIVES the flat gate (0.62): its fringe is genuinely likelihood-broken, aided
by the array's smallest K.** Under real noise its genuine detection collapses to **0.10** — the other
0.33 of its Bayesian robustness was prior-pinning, real as robustness but NOT a data detection — and
the same tiny K makes it the dominant *residual-null* false-alarm source (§9.3). **Verdict: J0437's
robustness is genuine at zero noise and prior-pinned under real noise; tiny K cuts both ways.**

### 9.6 What the re-score changes and does not change

- **CHANGES the absolute counts and their meaning.** The clean genuine-detection ceiling is
  ~1.35/draw at zero noise and **~0.1–0.3/real under real noise** (vs the Bayesian 4.5 / 2.9 / 1.4);
  the Bayesian counts were ~70–90 % prior/trials-driven. **The flat gate ELIMINATES all
  wrong-certifications on real data (both arms → 0).**
- **DOES NOT change** the §2 B1.3 DAMPED verdict (that turns on N_cert growth, not the gate) or the
  §1 A→B *relative* price (Arm A two-layer 0.33 still > Arm B 0.13). It sharpens §1/§5: the
  prior-pinned anchor's robustness is genuine only at zero noise.
- **STOP flag honoured.** The null does not read exactly zero; the anatomy (small-K trials bar
  beaten by noise) is reported, with the pre-registered fix (absolute `dlnL` floor) flagged for the
  spec rather than silently applied.

Re-score execution: job **12496241** (dgx03, 420 s, deterministic; certs reproduce the §1 bank
exactly). 89 `flat_*.npz` ≈ 0.4 MB. Nothing committed; no tracked file edited.
