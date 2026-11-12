# JOB 49 — Physically correct data: rerun (Parts D, E, F)

**P4 holds; the paper's main claim survives.** In every arm, physical ones included, pulsars with K ≥ 60 are certified correctly 2–3% of the time (9–13% under the finest-tooth K convention). Pulsars with K ≥ 600 are certified correctly 0/79 (≤ 1.6% under the finest convention), and no K ≥ 600 pulsar passes the gate.

**Plain-language summary.**
1. The machine draws a Gaussian phase error instead of summing phasors (H1 holds). But its share s counts each other source once, while each source actually puts an Earth term **and** a pulsar term in the bin. Done as D0 specifies (both terms), the physical sum has the same small-angle phase variance the machine already uses.
2. **Physical arm per D0 (P2), matched filter:** 36.56 [34.29, 38.92] psr/sky @ 0.9083 [0.897, 0.918], against the booked 36.27 @ 0.9009. That is +0.8%, inside the interval. **The headline does not move.** The incumbent phase-only decoder does worse on the same physical data: 27.13 @ 0.886.
3. **Amplitude limit:** with the D0 physical data there is still no 0.90 crossing at A = 8.3 or 10.6 ×10⁻¹⁵ (max train purity 0.86 and 0.74). P3 fails.
4. If the data are physical *at the machine's own share* (P1, one phasor per source), the brief's expectation does appear: 39.69 @ 0.981, and a crossing at A = 8.3 (36.48 @ 0.909). That arm is internally inconsistent with the machine's Earth-term deprojection, which assumes two terms per source, so it is a bound, not a result.
5. **Heatmap (share as labelled):** physical data help where the source is weak. Along β = 0.2 the recovered fraction rises 2–3× at σ_d = 1–30 pc (1 pc: 0.24 → 0.44). Yet it stays ≤ 0.10 at σ_d ≥ 10 pc. The 90% contour moves right only at high β.
6. **Criterion refit:** the exact-coherence form (E-rice) describes the physical data with flat residuals against β. E-old leaves a ±0.08–0.10 trend. P5 holds.
7. **A defect in booked code, independent of jitter:** the booked gate sits on a smoothing-edge artifact (audit G1). With an edge-safe gate the headline's real 0.90 train crossing gives 39.69 @ 0.876 out-of-sample, below the contract.

---

## Arms and the switch

All arms run on **identical** skies, pulsars, priors, geometry, L_true and certification sets (common random numbers). The incumbent rng stream is consumed exactly as booked. Every physical datum comes from a separate seeded stream, `np.random.default_rng([49, seed0, sky, pulsar, rate·10⁶, mult, gauss])`. In incumbent mode the switch collector is **bit-identical** to the booked `collect36m` / `_fixed` at every configuration (13,071–14,402 records each; `identical: true` in every `job49_data/D2D3_*.json`).

| Arm | Datum per (pulsar, comb), units of the source amplitude | Decoder |
|---|---|---|
| **incumbent** | ε = κL + N(0, σ_inc²) (booked) | booked wrapped quadratic, weight 1/σ_inc² |
| **P1 D-inc / D-mf** | z = e^{iκL}e^{iη}[1 + Σ_j a_j e^{iθ_j} + CN(0, R) + CN(0, 2/ρ²)]. **One** phasor per other source, a_j² = D·c_j/Pt, so total confusion power = the incumbent's res·D/Pt | D-inc: booked decoder on arg(z). D-mf: Σ 2e^{−s_fz²/2}Re(z e^{−iκL})/P_n − L²/2σ_d², with P_n = res·D/Pt + 2/ρ² |
| **P2 D-inc / D-mf** (**the D0 arm**) | as P1 but **two** independent phasors per other source (Earth + pulsar term), and sea ×2 | same; P_n = 2·res·D/Pt + 2/ρ² |
| P2g (control) | P2 power, with the whole confusion as one complex Gaussian (no explicit phasors) | same |

Switch: `job49_common.collect36m_sw(..., jitter_model="incumbent")`. The default is incumbent; the alternatives are `"p1_inc"`, `"p1_mf"`, `"p2_inc"`, `"p2_mf"`, `"p2g_inc"`, `"p2g_mf"`. Every record carries every arm's outcome, and `set_arm(rec, arm)` selects one. The heatmap equivalent is `job49_D1.purity_cell_sw()`, which returns (incumbent, physical D-inc, physical D-mf).

### D0 — what the population does not supply, and how I supplied it

- **Per-pulsar phases** of every contributor: not in the population. I drew them U(0, 2π), independently per (pulsar, contributor, term).
- **The other sources' contributions:** `inv_rate`'s res is a sum over about 18,700 in-window sources (the ±100/T window covers the whole band). The **64 largest** contributions per comb are explicit phasors with their own sinc²-weighted amplitudes; they carry 30–73% of res in the sky-0 check. The remainder plus the sea is an exact-power complex Gaussian, which by the CLT is accurate for thousands of small terms. P2g replaces even the top 64 with a Gaussian, and the result barely changes (D2: 36.48 vs 36.56), so the explicit/Gaussian split is not driving anything.
- **The own Earth term** of the source: treated as removed, as the incumbent does (it never enters res).
- **Earth-term deprojection:** the incumbent's factor D = (1 − f_E E_K)(1 + K/(N−1)) is applied uniformly to the confusion power. That is the machine's own blended factor; a per-term treatment would remove only Earth phasors.
- **Pulsar-term frequencies:** placed at the source's Earth frequency f2, like the incumbent ledger. No chirp retardation, so the 15–40 nHz pile-up physics is not in either arm.
- **Measurement noise:** CN with total power 2/ρ², so its small-angle phase variance is exactly the incumbent's 1/ρ².
- **The fuzzy template term:** a template error, not a power share, so it is applied as a pure phase rotation η ~ N(0, δκσ_d), as in the incumbent.
- **Validation:** per-comb phasor sums reproduce `inv_rate`'s res to 4×10⁻¹⁶. The mean physical coherence matches the Rice × fuzzy prediction: P1 0.6149 vs 0.6163, P2 0.5059 vs 0.5059.

---

## D1 — heatmap

**Reproduction:** `job47c_grid.csv`: max |Δ| = **0**. `job47d_grid.csv`: max |Δ| = **0**, 0 n_real mismatches. This needed the rounded grid that 47d reads from the 47c CSV (audit G2).

**Physical arm on the heatmap:** the axis *defines* the share, so the datum is z = e^{iκL}(1 + C) with E|C|² = (1 − s)/s. C is each comb's **actual contributor mix** (top 64 phasors + Gaussian remainder) rescaled to that power. There is no measurement noise and no fuzzy term (47c has neither). Each arm gets the full 47d treatment: band cells re-measured at 1000 realisations, contour from the 3×3-mean field, cells shown raw.

Figure: `job49_D1_heatmaps.png`. Three maps on one colour scale with **β on the vertical axis** and the old σ_φ on the right axis, plus an overlay panel of the three 0.9 contours. Data: `job49_D1_grid.csv`, `job49_data/D1_rows.json`, `D1_grid.npz`.

**σ_d at which the 0.90 fraction is crossed** (dedicated rows, 40 log-spaced σ_d, 1000 realisations per cell):

| β | old σ_φ | incumbent | physical D-inc | physical D-mf | incumbent at σ_φ/√2 (P2 prediction) |
|---|---|---|---|---|---|
| 0.8 | 0.50 | 1.12 pc | 2.33 pc | **2.94 pc** | 3.11 pc |
| 0.5 | 1.00 | 0.226 | 0.223 | 0.343 | 0.519 |
| 0.3 | 1.53 | 0.144 | 0.158 | 0.176 | 0.178 |
| 0.2 | 2.00 | 0.141 | 0.147 | 0.157 | 0.143 |

At β ≤ 0.3 every arm's 0.9 crossing sits at σ_d ≈ 0.14–0.18 pc. That is the prior-only left edge: the 6σ_d window is narrower than the coarsest tooth. So these crossings measure the prior, not the data.

**Recovered fraction along β = 0.2 (σ_φ = 2.0)**, 4000 realisations per cell, ±1 binomial s.e.:

| σ_d [pc] | 0.1 | 0.3 | 1 | 3 | 10 | 30 | 100 |
|---|---|---|---|---|---|---|---|
| incumbent | 0.982±.002 | 0.595±.008 | 0.243±.007 | 0.096±.005 | 0.030±.003 | 0.010±.002 | 0.006±.001 |
| physical D-inc | 0.986±.002 | 0.707±.007 | 0.393±.008 | 0.206±.006 | 0.088±.004 | 0.035±.003 | 0.010±.002 |
| physical D-mf | 0.980±.002 | 0.749±.007 | 0.443±.008 | 0.231±.007 | 0.102±.005 | 0.035±.003 | 0.010±.002 |

**How to read D1 against the ledger point.** The heatmap's axis is a *labelled* share, and the physical arm takes that label as the true share. When machine combs are placed on the map through s = Pt/(Pt + res), their true share is s/(2 − s). Read that way, the incumbent map is already small-angle-correct, and the physical gain is only the Rice-shape gain at low β.

## D2 — headline (marginalised A, fuzzy on, full δμ, 100 train / 100 held-out skies)

Booked procedure: gate re-frozen on train skies by `book()`, applied blind to the held-out skies. Bootstrap per sky ID (per-ID counts), 400 resamples, 10/90% bands (the booked convention). Source: `job49_D2D3_table.csv` / `job49_data/D2D3_marg.json`.

| Arm | ungated purity | **yield @ purity (re-frozen gate)** | gate | yield @ purity at the **old threshold, unchanged** | diagnostic: edge-safe gate, 0.5–100% (not the booked procedure) |
|---|---|---|---|---|---|
| **incumbent** (reproduced) | 0.680 | **36.27 [34.02, 38.75] @ 0.9009 [0.890, 0.912]** | 58.6% | — | 39.69 [37.52, 41.94] @ 0.876, gate 66.1% |
| P1 D-inc | 0.750 | 38.72 [36.41, 41.17] @ 0.958 | 58.7%* | 38.60 @ 0.959 | 50.65 @ 0.883, 82.3% |
| P1 D-mf | 0.772 | 39.69 [37.35, 42.27] @ 0.981 | 58.8%* | 39.50 @ 0.981 | 52.81 @ 0.884, 85.7% |
| **P2 D-inc** | 0.648 | 27.13 [24.91, 29.71] @ 0.886 [0.874, 0.897] | 46.8% | 34.37 @ 0.854 | 27.13 @ 0.886, 46.8% |
| **P2 D-mf** | 0.691 | **36.56 [34.29, 38.92] @ 0.9083 [0.897, 0.918]** | 58.6%* | 36.57 @ 0.908 | 40.04 @ 0.892, 65.6% |
| P2g D-inc | 0.649 | 30.45 [28.13, 33.00] @ 0.888 | 51.0% | 34.81 @ 0.865 | 30.45 @ 0.888 |
| P2g D-mf | 0.689 | 36.48 [34.28, 38.91] @ 0.906 | 58.6%* | 36.46 @ 0.906 | 40.83 @ 0.882 |

\* A gate at 58.6–58.8% is the `book()` edge artifact (audit G1): train purity is still above 0.90 at the end of the 60% search range. For those arms the booked yield is **censored**, a lower bound on the yield at 0.90.

## D3 — amplitude limit (booked T3 harness, fixed count-scale)

All five booked rows reproduce exactly with untouched code (`_fixed` exec'd from job43c:790-829).

| A (10⁻¹⁵) | rate | incumbent (booked) | P1 D-inc | P1 D-mf | **P2 D-inc** | **P2 D-mf** |
|---|---|---|---|---|---|---|
| 3.7 | 0.33 | 42.35 @ 0.997* | 42.45 @ 1.000* | 42.44 @ 1.000* | 41.91 @ 0.987* | 42.24 @ 0.995* |
| 5.0 | 0.60 | 40.99 @ 0.971* | 41.93 @ 0.992* | 42.12 @ 0.997* | 38.97 @ 0.926* | 40.64 @ 0.963* |
| 6.4 | 1.00 | 31.57 @ 0.9035 | 39.15 @ 0.955* | 40.06 @ 0.977* | 2.84 @ 0.863 (gate 4.7%)† | 27.97 @ 0.893 |
| 8.3 | 1.70 | none (max 0.898) | 10.45 @ 0.913 (gate 17%) | **36.48 @ 0.909** | none (max 0.837) | none (max 0.863) |
| 10.6 | 2.70 | none (max 0.709) | none (0.840) | none (0.884) | none (0.656) | none (0.735) |

\* edge-censored (see D2). † `book()` takes the *first* train crossing, and P2 D-inc's train curve dips below 0.90 at a 4.7% gate before recovering. This is procedural fragility, not physics; the edge-safe diagnostic agrees (2.84). Diagnostic yields with the edge-safe gate are in `job49_data/analysis.json` → `diag_gate`. For example, at A = 3.7 every arm is 56–58 psr/sky @ 0.906–0.908.

**A 0.90 crossing exists:** incumbent up to A = 6.4. **P2 (D0) arms: up to 6.4 only; no new crossing.** P1 arms: up to 8.3. Nothing at 10.6.

## D4 — usable sources (200 headline skies)

From `job49_data/analysis.json` → D4. Inventory stage per sky, before the per-pulsar depro and fuzzy terms:

| per sky | mean | median |
|---|---|---|
| candidates (ρ ≥ 3, arr ≥ 3, 2–30 nHz, cap 1300) | 1160.7 | 1127 |
| pass the old cut σ < π (≡ β > 0.092 under σ² = (1−β)/β) | **16.6** | 14 |
| … after the booked top-12-by-nats truncation | 10.9 | 12 |
| pass the physical cut, one-term ledger | 1153.4 | 1126 |
| pass the physical cut, two-term ledger (D0) | 1102.7 | 1093 |
| … after top-12 truncation | 12.0 | 12 |

**The physical cut, and why.** The old cut σ < π admits a source whose *assumed* coherence e^{−σ²/2} is at least e^{−π²/2} = 0.0072 (v4 cell 15:10 calls it "positive information"). Under the Rice law the same coherence floor needs only γ ≥ 6.6×10⁻⁵. With physical data almost every candidate carries positive phase information. The old cut was a property of the Gaussian model, not of the data. Under physical data the binding limit on N_eff becomes the **top-12 truncation**, not readability.

Per certifiable pulsar, the decode uses 9.5 combs on average (median 12). Median incumbent σ over decode combs is 1.19 rad. Median comb share is s = 0.414 as the machine computes it; the true two-term share is 0.303. Job 40's "14.4 per sky, median 1.84 rad" is the inventory-stage count (here 16.6 / median 14), before depro.

Under P2 the incumbent and physical arms see the **same** usable set, because selection is a power-ledger quantity. A physical cut would only add sources. Re-running the inventory with nmax > 12 is outside this job's untouched-procedure scope.

## D5 — what decides the answer for the paper

Headline configuration, harness priors (synthetic px/pb; median σ_d = 1.09 pc; census widths are not wired into this harness). Fraction of certifiable pulsars whose argmax is the true tooth, with a one-sided 95% Clopper–Pearson upper limit and gated certifications. Source: `job49_data/D5_bothK.json`.

| subset | n | incumbent | P1 D-inc | P1 D-mf | **P2 D-inc** | **P2 D-mf** | gated (any arm) |
|---|---|---|---|---|---|---|---|
| K ≥ 60 (coarsest tooth, job43c:551 convention) | 1432 | 0.020 (UL 0.027) | 0.027 (0.035) | 0.027 (0.035) | 0.021 (0.028) | 0.021 (0.028) | 0 passed the gate |
| K ≥ 600 (coarsest) | 79 | 0/79 (UL 0.037) | 0/79 | 0/79 | 0/79 | 0/79 | 0 |
| K ≥ 60 (finest tooth, notebook lnKi convention) | 2869 | 0.091 (0.101) | 0.131 (0.142) | 0.133 (0.144) | 0.089 (0.098) | 0.096 (0.106) | 12/14, 14/14, 14/14, 3/6, 11/14 |
| K ≥ 600 (finest) | 620 | 0.011 (0.021) | 0.013 (0.023) | 0.016 (0.027) | 0.011 (0.021) | 0.016 (0.027) | 0 |

**Under every arm and both K conventions, K ≥ 60 stays below 20% and K ≥ 600 stays below 5%.** Physical data change these by about 4 points at most (finest K ≥ 60: 0.091 → 0.133). Realistic priors still fail.

---

## Part E — criterion refit

Data: per-pulsar argmax outcomes (headline configuration, K ≥ 2, n = 4532). The predictor is r = Σ_k f(β_k)/(2 ln K) per pulsar; the brief's N_eff·f(β) generalised to unequal shares.

- **E-old:** f = e^{−(1−β)/β}.
- **E-rice:** f = ⟨cos δφ⟩²(β), the exact Rice expression of §1a.
- **β_k by dataset:** for incumbent data, the generator's implied share 1/(1 + σ_inc²). For P1/P2 data, the datum's true share 1/(1 + m·res·D/Pt + 2/ρ²) with m = 1 or 2, and the fuzzy term as a coherence factor e^{−s_fz²} in both forms.
- **Fit:** logistic P(correct) = σ(a + b ln r) by maximum likelihood. c(p) is the r at which P = p. Collapse R² is computed over 12 quantile bins of r. Residuals are taken in 6 quantile bins of each pulsar's median comb β. Source: `job49_data/analysis.json` → E.

| data / decoder | form | c(0.5) | c(0.7) | c(0.8) | c(0.9) | c(0.95) | R² | lnL | residual vs β (low β → high β) |
|---|---|---|---|---|---|---|---|---|---|
| incumbent | **E-old** | 0.684 | 1.03 | 1.34 | 1.99 | 2.86 | 0.998 | **−1406.9** | −0.003 … −0.016 (flat) |
| incumbent | E-rice | 1.04 | 1.59 | 2.08 | 3.12 | 4.54 | 0.996 | −1408.8 | −0.018 … +0.008 (flat) |
| P1 D-inc | E-old | 0.506 | 0.734 | 0.931 | 1.33 | 1.85 | 0.998 | −1354.6 | +0.009 … −0.020 |
| P1 D-inc | **E-rice** | 0.691 | 1.04 | 1.35 | 2.00 | 2.87 | 0.998 | **−1347.7** | −0.023 … −0.001 |
| P1 D-mf | E-old | 0.463 | 0.646 | 0.799 | 1.10 | 1.47 | 1.000 | −1230.5 | +0.026 … −0.026 (trend) |
| P1 D-mf | **E-rice** | 0.627 | 0.904 | 1.14 | 1.62 | 2.23 | 0.999 | **−1209.3** | −0.007 … −0.005 (flat) |
| P2 D-inc | E-old | 0.416 | 0.675 | 0.918 | 1.46 | 2.23 | 0.997 | −1662.4 | **+0.065, +0.085 … −0.050, −0.077** |
| P2 D-inc | **E-rice** | 0.708 | 1.13 | 1.52 | 2.36 | 3.56 | 0.995 | **−1513.1** | −0.008 … −0.017 (flat) |
| P2 D-mf | E-old | 0.367 | 0.571 | 0.756 | 1.15 | 1.70 | 0.994 | −1593.5 | **+0.074, +0.103 … −0.055, −0.090** |
| P2 D-mf | **E-rice** | 0.626 | 0.959 | 1.26 | 1.89 | 2.76 | 0.993 | **−1428.1** | −0.003 … −0.024 (flat) |

**Which form describes which dataset.**
- **E-old** describes the incumbent data, as it must: they *are* wrapped-normal draws. ΔlnL is +1.9 in its favour.
- **E-rice** describes every physical dataset. ΔlnL is +6.9 (P1 D-inc), +21.2 (P1 MF), +149 (P2 D-inc) and +165 (P2 MF), and its residuals are flat across β.
- On the D0 data, E-old leaves a systematic trend: it under-predicts at low β by 7–10 points and over-predicts at high β by 5–9.
- Binned R² is high for all fits (≥ 0.993), so R² does not discriminate here. Use the likelihood and the β residuals.

c(0.9) here (1.1–2.4) is larger than the booked 1.01. That value came from the equal-share protocol; these are per-pulsar outcomes on the realistic harness, with unequal shares and the fuzzy term. The two are not the same quantity.

## Part F — what should not have changed

- **Pile-up counts.** They come from frequencies only: `f_src` (v4 cell 3) and `FP_ALL` (v4 cell 14:1586-1591, `FP_ALL = (f_src[None, :]**(-8/3) + (8/3)*_KCH_all[None, :]*_dt_all)**(-3/8)`). They are counted in cell 14 (0d) as `((f_src >= lo) & (f_src < hi)).sum()` and in `PAPER_figures/job46_partB_extract.py:35-37`. No phase is involved. **Recomputed once:** ≥ 50 nHz on the frozen sky (draw_population seed 3, pulsar 41): **Earth 232, pulsar 44**, unchanged (`job49_data/partF_counts.json`).
- **Share-versus-frequency figure** (`PAPER_figures/job46_partB_figure.py:11-14`): `P = h**2*f*T; … tot = np.bincount(b, weights=P, minlength=nb) + sea; return P/tot[b]`. It is a power ratio with no phases. It is unchanged by the jitter model. Two notes:
  - Its readability line `S_PI = 1/(1+np.pi**2)` (:10) *is* the σ² = (1−s)/s mapping. Under physical data that threshold has no special meaning (D4).
  - Its ledger is again one term per source, and only pulsar terms (sON).
- **Per-source shares** (`inv_rate`, job43c:498-505): pure power ledger. The incumbent and physical arms share the identical s and selection. The physical arm reads the shares; it never writes them.

---

## Pre-registered predictions

| | Prediction | Verdict |
|---|---|---|
| **P1** | Incumbent draws phase errors with variance (1−β)/β | **HOLDS.** job43c:749 and job47c:124. Part C: within 2 s.e. of WN((1−β)/β) at all β; in-situ KS p = 0.31. *Caveat:* β is the machine's one-term share, and with the two-term ledger that variance is the correct small-angle value. |
| **P2** | Heatmap 0.9 contour moves to larger σ_d at every β, by about √2 on the σ_φ axis | **PARTLY FAILS.** β = 0.8: yes, 1.12 → 2.94 pc (MF) against a √2 prediction of 3.11. β = 0.5: half the predicted move (0.226 → 0.343 vs 0.519), and D-inc does not move (0.223). β ≤ 0.3: the contour is pinned at the prior-only edge in all arms, so no test. Physical data do lift the β = 0.2 row 2–3× at 1–30 pc. |
| **P3** | Headline yield rises at 0.90 purity; a crossing appears at a new amplitude | **FAILS for the D0 arm (P2):** MF 36.56 vs 36.27, inside the band; no crossing at 8.3 or 10.6. *Holds only for P1* (one-term ledger): 39.69 @ 0.981 and a crossing at 8.3. P1 is inconsistent with the machine's own deprojection. |
| **P4** | K ≥ 60 < 20% and K ≥ 600 < 5% in the physical arm | **HOLDS** in every arm and both K conventions. K ≥ 60: 2.0–2.7% (coarsest), 8.9–13.3% (finest). K ≥ 600: 0/79 (UL 3.7%), 1.1–1.6% (finest). Nothing certified through the gate at K ≥ 600. |
| **P5** | E-rice fits the physical-arm data better than E-old | **HOLDS.** ΔlnL +7 to +165 in favour of E-rice; flat residuals vs β. E-old fits only the incumbent data. |

## Part G

**Not run.** Both reports are on disk, but the real-likelihood comparison (40 sky redraws × 116 pulsars through `CW_transition`) is a separate GPU campaign. Given your speed-first rule I did not start it unasked. Its first task would be to correct the √2 in `stagec_p2b.py:86` (audit A10) before comparing any σ_φ.

## Files

**Added** (all in `wheel_hasasia_examples/`, staged; nothing committed):
- Reports: `JOB49_STATUS.md`, `JOB49_AUDIT.md`, `JOB49_PHYSICAL.md`.
- `job49_common.py`: loader plus `collect36m_sw` (**switch `jitter_model`**, default `"incumbent"`), `book_sw`, `contributors`, `VARIANTS`/`ARMS`, `set_arm`.
- `job49_partC.py` → `job49_partC_table.csv`, `job49_partC_hist.png`, `job49_partC_verdict.json`.
- `job49_D2D3.py` (`marg` | `rate …`) → `job49_data/D2D3_*.json`, plus `rec_*.pkl` (120 MB, **not staged**; regenerate in about 4 min).
- `job49_D1.py` → `job49_D1_grid.csv`, `job49_D1_heatmaps.png`, `job49_data/D1_rows.json`, `job49_data/D1_grid.npz`.
- `job49_analysis.py` → `job49_D2D3_table.csv`, `job49_data/analysis.json`.
- `job49_data/D5_bothK.json`, `job49_data/partF_counts.json`.
- Logs: `job49_logs/*`, including `partB/` with the two reference scripts and their outputs.

**Changed:** none. No booked file was edited. The loader's only patch is in memory: it does not spawn the GPU gate worker, which affects speed only.

## Wall time

| Part | Time |
|---|---|
| A | ≈ 45 min reading |
| B | 32 s |
| C | 0.5 s |
| D2 + D3 | 4.2 min wall (3 processes in parallel; ≈ 65 s per configuration after a 58 s load) |
| D1 | 2.4 min (run twice: the first run exposed the 47d rounding) |
| D4/D5/E analysis | 1.0 min |
| F | 1 min |
| G | not run |
