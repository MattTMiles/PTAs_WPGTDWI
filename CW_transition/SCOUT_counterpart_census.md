# SCOUT — EM-Counterpart SMBHB Census vs the CW Certification Operating Point

*Agent SCOUT. CPU-only, literature census + closed-form arithmetic, no simulation.
Context read from `project_progress.md` (Stage C headline). No doc-write authority
beyond this file.*

Last updated: 2026-07-08

---

## 0. The question and the operating point

The project pivoted to **EM-targeted** CW sources: a confirmed counterpart hands us
the sky position, so the pulsar-term fringe / blind-localization wall (Stage C:
tens–hundreds deg² for a blind PTA search) is bypassed. What remains is a **joint**
population question:

> Does nature supply sources that are *simultaneously* (a) LOUD enough —
> per-source matched-filter SNR ~ 10 at 15 yr, i.e. **h ≈ 10⁻¹³·⁷⁵ = 1.78 × 10⁻¹⁴** —
> AND (b) counterpart-IDENTIFIED to the ~22 arcsec (≈0.006°) tolerance?

**Certification operating point:** `log₁₀ h = −13.75` (`h = 1.78 × 10⁻¹⁴`).

**Reference sensitivity anchor:** the NANOGrav 15-yr *individual-source* search sets a
sky-averaged 95% upper limit of **h₀ < 8 × 10⁻¹⁵ (log₁₀ h = −14.1)** at its most
sensitive frequency ~6 nHz [NG15 CW, Agazie+ 2023, ApJL 951 L50; arXiv:2306.16222].
**The operating point (−13.75) sits ~0.35 dex ABOVE this limit.** A source that loud,
in the sensitive band, would already be a *detection*, not a needle — this framing
tension recurs in the verdict below.

Strain arithmetic uses the standard circular-binary amplitude
`h₀ = 2 (G ℳ_c)^{5/3} (π f_gw,rest)^{2/3} / (c⁴ d_L)`, rest-frame chirp mass ℳ_c and
GW frequency `f_gw = 2 f_orb,rest`, flat ΛCDM (H₀=70, Ω_m=0.3) for d_L. Script:
`scratchpad/hcalc.py` (reproducible).

---

## 1. Candidate table

| Candidate | z | d_L (Mpc) | f_gw,rest (Hz) | ℳ_c range (M_⊙) | **log₁₀ h (implied)** | Astrometry class | Credibility status |
|-----------|---|-----------|----------------|-----------------|----------------------|------------------|--------------------|
| **OJ 287** | 0.306 | 1590 | 6.9e-9 (P_orb,rest≈9.2 yr) | ℳ_c≈1.0e9 (M₁=1.8e10,M₂=1.5e8, q≈0.008) | **−15.9** (lit. −15.8…−15.0) | radio **VLBI, mas** | Best-studied SMBHB; 12-yr outburst clock fits precessing-orbit model to days [Dey+2018; Valtonen]. Binary widely accepted but extreme q → GW-quiet. |
| **PKS 2131−021** | 1.283 | 8990 | 3.0e-8 (P_rest=2.082 yr) | 3e8 … 1.5e10 (UL) | **−17.1 … −14.2** | radio **VLBI, mas** | Strong coherent radio–optical sinusoid, ~2 yr [O'Neill+2022; Sheu+2024]. Compelling but Doppler-boost model favors *low* mass; NG15 targeted UL h₀<1e-14, ℳ_c<1.5e10. |
| **3C 66B** | 0.0213 | 93 | 6.2e-8 (P_obs=1.05 yr) | 6.9e8 (PPTA-DR3 UL) … 7e9 (Sudou'03) | **−14.3 (UL) … −12.6 (excluded)** | radio **VLBI, mas** | Nearest candidate. Original astrometric-wobble mass (Sudou 2003, total 5.4e10) **RULED OUT** by PTAs [NG11 Arzoumanian+2020; PPTA-DR3 2025 log₁₀h₀<−14.44, ℳ_c<6.9e8]. |
| **PG 1302−102** | 0.278 | 1430 | 1.6e-8 (P_rest=4.03 yr, P_obs=1884 d) | 3e8 … 5e9 | **−16.5 … −14.4** | optical/radio, **Gaia mas + VLBI** | Prototype CRTS periodic quasar [Graham+2015a]. Periodicity contested — extended baseline weakens it; red-noise-mimic risk [Liu+2018; Kovačević+2019]. |
| **Graham+ 2015 sample** | 0.2–1.5 | — | ~1e-8 (P_obs 1–5.9 yr) | ≲1e9 typical | **≲ −15** typ. | optical **~0.1–1″** (CRTS) | 111 candidates / 243k quasars. Population-level: most consistent with red noise; only a handful survive extended baselines [Vaughan+2016; Sesana+2018 PTA test]. |
| **Charisi+ 2016 sample** | — | — | ~3e-8 (P_obs few×100 d) | ≲1e9 | **≲ −15** typ. | optical **~0.1–1″** (PTF) | 33 candidates / 35k quasars, low q~0.01. Significance reassessed downward with null-signal templates [Robnik+2024, MNRAS 534 1609]. |
| *Post-2024:* SDSS J0252−0028 | 1.53 | — | ~2e-8 (P_rest≈1.7 yr) | ≲1e9 | ≲ −15 | optical mas (Gaia)+ZTF | ZTF/CRTS/Gaia ~4.6-cycle periodicity, 99.95% [Chen+2024]. Single new candidate; unconfirmed. |
| *Post-2024:* mpc-catalogue | var | — | var | var | mostly ≲ −15 | optical ~0.1–1″ | ~milli-pc-separation catalogue from long-baseline photometry [arXiv:2505.06656, 2025]; Gaia DR3 periodic-AGN search [arXiv:2505.16884, 2025]. Population candidates, individually low-credibility. |

**Arithmetic shown (representative, from `hcalc.py`):**
- OJ 287: ℳ_c=1.0e9, f=6.9e-9, d_L=1590 Mpc → log₁₀h = **−15.86** (matches lit. published
  1.6e-16…9.7e-16 = −15.8…−15.0). The extreme mass ratio (q≈0.008) suppresses ℳ_c despite
  the 1.8×10¹⁰ M_⊙ primary — it is GW-quiet.
- 3C 66B at the PTA-allowed mass (ℳ_c=6.9e8, f=6.2e-8, d_L=93 Mpc) → **−14.27**; at the
  original Sudou mass (ℳ_c≈7e9) → **−12.59** (but that mass is *excluded by PTA
  non-detection*).
- PG 1302 at optimistic ℳ_c=5e9 → **−14.42**; at fiducial ~2e9 → −15.08.
- PKS 2131 at the NG15 upper-limit mass ℳ_c=1.5e10 → **−14.23**; at Doppler-model-favored
  ~3e8 → −17.06.

---

## 2. THE VERDICT COLUMN — implied h vs the −13.75 operating point

| Candidate | best-case log₁₀ h | Δ from −13.75 | Clears −13.75? | Within 1 dex? | Catch |
|-----------|-------------------|---------------|----------------|---------------|-------|
| OJ 287 | −15.0 (lit. max) | −1.25 | **No** | **No** (~1.25 dex short) | GW-quiet: q≈0.008. Loud EM ≠ loud GW. |
| PKS 2131−021 | −14.2 (at ℳ_c=1.5e10 UL) | −0.45 | No | **Yes, at the EXCLUDED upper mass** | Fiducial (low-mass Doppler model) is −16 to −17: fails by 2–3 dex. |
| 3C 66B | −14.3 (PTA-allowed UL) | −0.55 | No | **Yes, at the UL** | Only clears (−12.6) at the mass PTAs already *ruled out*. |
| PG 1302−102 | −14.4 (ℳ_c=5e9) | −0.65 | No | **Yes, at optimistic mass** | Fiducial −15; periodicity itself contested. |
| Graham+/Charisi+ pop. | ≲ −15 typ. | ≲ −1.25 | No | No | Population dominated by red-noise mimics. |

**VERDICT SUMMARY.**
**No named credible SMBHB candidate clears h = 10⁻¹³·⁷⁵ at its literature-favored mass.**
The three that come within 1 dex — 3C 66B, PG 1302, PKS 2131 — do so *only* at the loud
edge of their mass range, and that loud edge is in each case either **already excluded by
PTA non-detection** (3C 66B Sudou mass; PKS 2131 at its ℳ_c<1.5e10 upper limit) or the
*optimistic* tail of a **contested** periodicity (PG 1302). OJ 287 — the single most
secure binary — is ~1.5 dex too quiet because its mass ratio is extreme.

**Structural reason (the scissors).** The operating point log₁₀h=−13.75 sits *above* the
NG15 sky-averaged CW upper limit (−14.1). Any source genuinely at the operating point,
in-band, would already be an individual-source *detection*. So "certifiable and
counterpart-identified" is asking for a source **louder than everything current PTAs have
failed to detect** — the census operating point is on the wrong side of the current
exclusion curve for essentially the whole sky. This is not a counterpart-availability
problem; it is a source-loudness problem.

---

## 3. Population arithmetic — joint (loud AND identified)

Order-of-magnitude, published numbers, assumptions stated.

**(A) How many resolvable single sources exist at all (current PTAs)?**
Population synthesis constrained to reproduce the measured GWB gives the mean number of
*individually detectable* SMBHBs as **N̄ ≲ 0.01 (current sensitivity), ≲ 0.1 (any SMBH
distribution consistent with the SGWB)** [Sato-Polito, Zaldarriaga & Kaiser 2024,
arXiv:2406.17010]. Non-detection of a CW is fully expected. Older, more optimistic
population work (Sesana+2009) suggested 5–10 sources above the background at
f_gw > 10⁻⁸ Hz in a 5-yr SKA set — but that predates the measured-amplitude
normalization and is now regarded as an upper edge.

**(B) SKA-era projection.** Detection probability of *any* single source rises to ~40%
(SKA1, 10 yr) / ~80% (SKA2) [Rosado, Sesana & Gair 2015, arXiv:1501.00127 / MNRAS]. So
O(1) resolved source is a next-decade (SKA) expectation, not a 15-yr one.

**(C) Joint fraction — loud AND counterpart-identified.**
Assumptions: (i) fraction of resolvable sources that are radio-loud/optically-periodic
enough to *have* an identified AGN counterpart ~ **0.1–1** (bright SMBHBs are AGN; the
counterpart is the easy part — see §4); (ii) counterpart astrometry always beats 22
arcsec (§4), so the localization-volume cut is ~1 (the counterpart *is* the position).
Then

> **N_joint ≈ N̄_detectable × f_counterpart ≈ (≲0.1) × (0.1–1) ≈ 10⁻² … 10⁻¹ (current);
> O(0.1–1) in the SKA era.**

The joint number is **not suppressed by the counterpart requirement** — it is set almost
entirely by N̄_detectable (loudness). Restated: the counterpart doesn't cost you; the
loudness does. **At current sensitivity the expected joint count is ≪ 1.** The pivot's
premise ("a counterpart bypasses the wall") is arithmetically correct *for the sky
position*, but does not manufacture a loud source where the population supplies none.

**Caveat.** These are GWB-anchored means; the variance is large and a single lucky nearby
massive binary (a "cosmic-variance jackpot") is not excluded — that is exactly what a
targeted EM list hedges against. The value of the EM-target program is therefore as a
**variance play on the loud tail**, not an expectation-value detection.

---

## 4. Astrometry requirements row (the EM column of the requirements table)

Tolerance: **22 arcsec = 0.006°.** Position-precision classes vs tolerance:

| Counterpart class | Typical σ_pos | Margin vs 22″ | Qualifies? | Examples |
|-------------------|---------------|---------------|------------|----------|
| Radio **VLBI** | ~0.1–1 **mas** | ~10⁴–10⁵× | **Trivially** | OJ 287, PKS 2131, 3C 66B |
| **Gaia** optical (QSO) | ~0.1–1 **mas** | ~10⁴–10⁵× | **Trivially** | most bright quasars |
| Optical survey centroid | ~0.1–1 **arcsec** | ~20–200× | **Trivially** | CRTS/PTF/ZTF periodic QSOs |
| Ground photometric (crowded/faint) | ~1–5 **arcsec** | ~4–20× | **Yes** | faint periodic-QSO tail |
| Large-PSF / unresolved | ~5–20 **arcsec** | ~1–4× | **Marginal** | rare, poorly-localized only |
| **Blind PTA CW** localization | ~tens–hundreds **deg²** | ~10⁻⁶× | **FAILS** | (this is the wall being bypassed) |

**Conclusion for the requirements table.** The 22-arcsec EM tolerance is met **with
enormous margin (≥10⁴×) by every credible AGN counterpart** — VLBI and Gaia are ~10⁴–10⁵×
tighter, and even a plain optical survey centroid is ~20–200× tighter. **EM astrometry is
never the binding constraint.** The binding constraint is entirely on the GW side:
source loudness (§2) and, once a position is supplied, the *within-fringe / fringe-ID*
distance physics of Stage C. The EM column should be recorded as **"trivially satisfied
(margin ≥10⁴×) for any VLBI/Gaia/optical counterpart; the requirement is that a
counterpart EXISTS and the source is loud, not that its position is precise."**

---

## 5. Cross-links & references

- Operating point / Stage C context: `project_progress.md` §"Stage C headline". The EM
  pivot is directly foreshadowed by arXiv:2512.10729 (*Two-Dimensional Pulsar Distance
  Inference from Nanohertz GW*, 2025) and the targeted-search machinery of
  Charisi/Taylor/Witt/Runnoe 2023 (arXiv:2304.03786, PRL 132 061401) — Earth-term-only
  approximation, ~100× faster targeted MM searches.
- VLBI-enabled CW localization in the SKA era: arXiv:2606.28721.
- NG15 individual sources: arXiv:2306.16222 (h₀<8e-15 @6nHz). NG15 targeted searches:
  arXiv:2508.16534. PPTA-DR3 CW: arXiv:2508.13944; 3C 66B multimessenger: arXiv:2508.20007.
- OJ 287: arXiv:1809.09473 (PTA detectability), Dey+2018. PKS 2131: O'Neill+2022
  (ApJL 926 L35), Sheu+2024 (arXiv:2407.09647), eccentricity arXiv:2606.05343.
- Populations: Graham+2015a (CRTS, 111); Charisi+2016 (PTF, 33); Robnik+2024
  (MNRAS 534 1609, PTF reassessment); Sesana+2018 (arXiv:1703.10611, PTA binary test);
  Chen+2024 (SDSS J0252, ZTF); arXiv:2505.06656, arXiv:2505.16884 (2025 catalogues).
- Population forecast: Sato-Polito+2024 (arXiv:2406.17010, N̄≲0.1); Rosado+2015
  (arXiv:1501.00127, SKA 40/80%); Sesana+2009.

**One-line takeaway for the requirements table:** EM astrometry is free (≥10⁴× margin);
the joint program lives or dies on the GW-loudness tail, where the population expectation
is N̄ ≲ 0.1 and no named credible candidate clears the −13.75 operating point at its real
mass.
