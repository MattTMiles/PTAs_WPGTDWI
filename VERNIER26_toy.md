# VERNIER-26 — does the loop sharpen the background measurement?

**EXPLORATORY. NOT pre-registered.** Toy model, cronus, CPU/numpy, **27 s**.
No MCMC: Fourier-domain quadratic estimators and a profiled power-law fit only.
**200 realisations** (the brief budgeted 20; at 0.13 s each I ran 10× and the
variance columns are the deliverable, so they wanted the samples). Every
comparison is **paired** — same sky, same noise, same array, different analysis.

Clean room: imports only `vernier24` (and through it `vernier_toy`).

Code `hpc_harbor/vernier/vernier26.py` · data `reports/vernier26.npz` ·
log `reports/vernier26.log`.
Staged, not committed.

Reuses VERNIER-24's Arm B corrected sky unchanged: Pareto γ_eff = 1.56 through the
h(Mc, f, D) link, N_bin from VERNIER-4's N_c ~ f^(-11/3), top 300 per bin modelled
individually. **One thing is new and it is what makes the question askable: the faint
sea is drawn HELLINGS-DOWNS CORRELATED across the array** (Cholesky of S_sea Γ_ab),
not as an uncorrelated per-pulsar PSD. Without that there is no HD estimator to sharpen.

---

## Lead — yes, and through M-A and M-C. M-B, the mechanism everyone would have guessed, is NEGATIVE.

**(iii) beats (ii) decisively at every R and by a widening margin, and the scrambled
control shows the gain is in the distances.** RMSE on the reconstructed background
amplitude, paired bootstrap over 200 realisations:

| R | (iii) vs (ii) | 95 % CI | p(gain ≤ 0) | **SCRAMBLED vs (ii)** |
|---|---|---|---|---|
| 4 | **+18.1 %** | [+13.9, +22.5] | < 1e-4 | −2.0 % |
| **8** | **+26.1 %** | [+21.1, +30.6] | < 1e-4 | −1.1 % |
| 16 | **+38.9 %** | [+34.8, +43.0] | < 1e-4 | +5.0 % |
| 25 | **+38.9 %** | [+34.5, +43.2] | < 1e-4 | +2.9 % |

**The hazard control passes: wrong distances buy −2 % to +5 % where right distances
buy +18 % to +39 %.** Modelling the pulsar term with a scrambled distance is worth
essentially nothing, which is the whole content of the CONDITIONED label.

**But the mechanism is not the one the brief's ordering implies.** Differencing the
three treatments:

| mechanism | where it is read | R = 4 | R = 8 | R = 16 | R = 25 |
|---|---|---|---|---|---|
| **M-A** bias, pulsar-term leak into CRN | AUTO estimator, (ii) → (iii) | **+0.118** | **+0.191** | **+0.283** | **+0.371** |
| **M-B** reassignment | RMSE of total, (i) → (ii) | **−4.9 %** | **−7.9 %** | **−31.8 %** | **−39.6 %** |
| **M-C** pulsar variance | HD variance floor, (ii) → (iii) | **+23.8 %** | **+27.1 %** | **+29.7 %** | **+31.1 %** |

**M-B is negative at every R.** Recovering a source and subtracting its Earth term
alone does not sharpen the background — it *degrades* it, by 5 % at R = 4 and by 40 %
at R = 25. **All of the gain comes from the pulsar terms**, split between a bias
channel (M-A) and a variance channel (M-C).

## M-C is real, and it is the unpublished one — flagged

**The Hellings-Downs variance floor falls 23.8 → 31.1 % as R goes 4 → 25, at ±1 %
on the mean, and the scrambled control moves it the WRONG way (−3.0 % to −5.8 %).**

| R | floor (i) | floor (ii) earth | **floor (iii) conditioned** | floor (iii-s) scrambled | **M-C gain** | HAZARD |
|---|---|---|---|---|---|---|
| 0 | 0.02424 | — | — | — | — | — |
| 4 | 0.02424 | 0.02284 | **0.01779** | 0.02354 | **+23.75 % ± 1.10** | −3.01 % ± 0.34 |
| 8 | 0.02424 | 0.02223 | **0.01653** | 0.02317 | **+27.12 % ± 1.01** | −4.23 % ± 0.39 |
| 16 | 0.02424 | 0.02164 | **0.01550** | 0.02284 | **+29.72 % ± 0.97** | −5.49 % ± 0.41 |
| 25 | 0.02424 | 0.02125 | **0.01491** | 0.02248 | **+31.12 % ± 0.94** | −5.75 % ± 0.37 |

*(floor = σ_A/A of the optimally weighted HD amplitude, computed from the per-pulsar
residual signal power.)*

**The logic is airtight by construction and that is the point.** Pulsar terms are
incoherent between pulsars, so they contribute **exactly nothing** to the mean of a
cross-correlation — no HD analysis can see them, and no amount of data removes them.
They sit in every pulsar's *auto* power, which is the variance normalisation of the
same estimator. So they are a pure variance floor, invisible to the statistic they
degrade. **Certified distances make them modellable, and the floor comes down by
a quarter to a third.**

**The size is set by how big they are, and they are big:** the median pulsar-term power
of the recovered sources is **97 % of the residual background at R = 4, rising to 169 %
at R = 25**. The loop's own catch is carrying a pulsar-term nuisance comparable to the
signal being measured.

Nobody has computed this. It is the one arrow in this instalment that is not a
re-derivation of something known.

## M-A is larger than M-C, and it is the known confound finally quantified

**Unmodelled sources inflate the common-red-noise amplitude by 49 % to 140 % in A²,
and Earth-term subtraction removes only two thirds of it.**

Bias in Â²/A²_true on the residual background, AUTO (CRN-like) estimator:

| R | (i) unmodelled | (ii) Earth only | **(iii) Earth + pulsar** | (iii-s) scrambled |
|---|---|---|---|---|
| 4 | **+0.492** | +0.103 | **−0.015** | +0.270 |
| 8 | **+0.758** | +0.168 | **−0.023** | +0.435 |
| 16 | **+1.114** | +0.251 | **−0.033** | +0.647 |
| 25 | **+1.400** | +0.327 | **−0.044** | +0.818 |

**Earth-only modelling leaves a +10 % to +33 % CRN bias that is entirely pulsar term,
and pulsar-term modelling drives it to −1.5 % to −4.4 %.** The scrambled control leaves
+27 % to +82 % — over half the bias survives a wrong distance, because a wrong-phase
template subtracts incoherently and puts back most of what it takes out.

**It shows in the slope too.** The CRN slope under Earth-only subtraction sits at
γ_auto = 4.16–4.20 against a true residual slope of 3.85 falling to 3.25 — **a tilt of
up to +0.95 in γ, purely from unmodelled pulsar terms.** Under (iii) the AUTO slope
tracks truth to 3.79/3.76, 3.62/3.58, 3.47/3.36, 3.35/3.17.

## Why M-B is negative — the finding I did not expect

Subtracting R sources in a bin projects out **2R complex directions** from a
116-pulsar vector, and those directions are Earth-term templates — precisely the
*coherent* part of the data. **The projection eats the background along with the
sources.** Two consequences, both in the tables:

1. **The cross-correlation detection weakens monotonically.** OS SNR runs
   12.78 → 7.70 → 5.66 → 3.78 → 2.91 across R = 0 → 25 under Earth-only subtraction.
   By R = 25 the residual background is barely detectable even though its amplitude is
   still well estimated.
2. **The reconstructed total gets worse, not better, under (ii).** RMSE 0.2020 at
   R = 0 → 0.2119 / 0.2180 / 0.2663 / 0.2822 at R = 4 / 8 / 16 / 25.

**A methodological point worth keeping.** Without the fitted-dof correction — referring
both estimators to the *projected* response Π Γ Πᴴ and Π diag(pw) Πᴴ rather than to Γ
and pw — the cross estimate of the residual background falls to **−1.00** (i.e. to zero)
by R = 25 and the loop appears to annihilate the background it was supposed to sharpen.
That is entirely an artefact of not accounting for the parameters spent. Every number
in this note is post-correction.

## The answer in the requested currency: % precision per recovered source

**Realised** RMSE on the reconstructed total amplitude, cross estimator, against the
fixed truth A = 2.4e-15:

| R | (i) | (ii) earth | **(iii) conditioned** | (iii-s) scrambled | **(iii) vs R = 0** | **per source** |
|---|---|---|---|---|---|---|
| 0 | 0.1010 | — | — | — | — | — |
| 4 | 0.1010 | 0.1060 | **0.0868** | 0.1076 | **+14.1 %** | **+3.5 %** |
| **8** | 0.1010 | 0.1090 | **0.0804** | 0.1101 | **+20.4 %** | **+2.6 %** |
| 16 | 0.1010 | 0.1332 | **0.0814** | 0.1265 | **+19.5 %** | **+1.2 %** |
| 25 | 0.1010 | 0.1411 | **0.0862** | 0.1366 | **+14.7 %** | **+0.6 %** |

On the **noise-limited floor** (what an optimally weighted analysis could reach):
+26.6 / +31.8 / +36.1 / +38.5 % at R = 4 / 8 / 16 / 25, i.e. **+6.7 / +4.0 / +2.3 /
+1.5 % per source**.

**Quote it as ~2.5 % per recovered source at R = 8, saturating.** The per-source number
falls by 6× from R = 4 to R = 25 because the loudest sources are recovered first and
the gain is proportional to the pulsar-term power removed. **Total realised gain peaks
near R = 8–16 at ~20 % and then turns over** — the projection cost of (ii) starts
outrunning the pulsar-term gain of (iii). The floor keeps improving because it does not
pay the projection cost.

## The slope

Truth γ = 13/3 = 4.333 for the total; the residual's true slope **genuinely flattens**
as the loop works, because the loudest sources sit in the lowest-frequency bins:
γ_true(residual) = 4.35 → 3.69 → 3.57 → 3.53 → 3.49 across R = 0 → 25.

**Treatment (iii) tracks the flattening to ±0.04** (3.80/3.76, 3.60/3.58, 3.40/3.36,
3.21/3.17 for estimate/truth). On the reconstructed **total**, where the answer should
stay at 4.333: (i) 4.37, (ii) 4.45–4.47, **(iii) 4.37–4.42**, scrambled 4.19–4.22.
**(iii) is the only treatment that neither steepens nor flattens the quoted total slope.**

## Controls

- **Pure-noise null** (no GWB, no sources), 200 realisations: Â²/A²_true = **−0.0021 ±
  0.0188** (auto) and **−0.0008 ± 0.0184** (cross); OS SNR **−0.010 ± 0.241**, s.e. of
  the mean 0.017. The mean is consistent with zero — the estimator does not manufacture
  a correlation. The spread is 0.24 rather than 1 because the weights assume a
  background the null does not contain, so the **SNR column is conservative by ~4× in
  the no-signal limit and is only ever compared across treatments.**
- **R = 0 rows** are identical across treatments by construction and are printed as such.
- **ORACLE hazard control** (S5.5.1): scrambled distances, everything else identical.
  Fails to reproduce M-C (−3.0 % to −5.8 % against +23.8 % to +31.1 %), fails to
  reproduce the RMSE gain (−2 % to +5 % against +18 % to +39 %), and leaves over half
  the M-A bias in place. **Treatment (iii) is CONDITIONED, not implementable, and the
  control is what licenses the claim.**
- **Measured HD-filter calibration** k = 1.0474 against the regularised antenna
  response, divided out. It is a constant, identical in every cell, so it cancels in
  every mechanism difference.

---

## What this does not settle

1. **CONDITIONED, not implementable.** (iii) uses oracle distances and oracle sky
   positions and frequencies for the recovered sources. VERNIER-24/25 say the loop
   certifies **zero** distances on the physical sky and needs ~10–11 accumulated combs
   (K\* ≈ 20–25 simultaneous detections) to certify any. **So this note prices a
   product the loop does not currently deliver.** It says what a certified distance is
   worth, not that one exists. The scrambled control bounds what a *wrong* certification
   would buy: nothing.
2. **The pulsar term is placed in the SAME frequency bin as the Earth term**, with
   amplitude ratio r_P = 1.196 (VERNIER-8) and a per-pulsar phase. In reality the pulsar
   term sits at an earlier, lower chirp frequency and can fall in a different bin. That
   simplification makes both the M-A leak and the M-C floor *easier* to model than they
   are, so both gains are optimistic.
3. **r_P = 1.196 breaks the textbook HD normalisation** (which assumes equal Earth and
   pulsar power). The Earth-coherent fraction is q = 1/(1 + r_P²) = 0.411, so the ORF
   used for both the sea and the filter is rescaled to 2q·HD off-diagonal. Using the
   textbook curve instead biases the cross estimator low by 18 %, which is a
   normalisation error and not a mechanism — but anyone comparing these numbers to a
   standard OS must apply the same rescaling.
4. **Source amplitudes are fitted; sky positions, frequencies and polarisation-basis
   geometry are not.** A real joint fit spends more parameters and pays a larger
   projection cost, so the M-B penalty reported here is a lower bound on the true one.
5. **Catalogue power is estimated as the coherent total ‖X_E â‖², not as a sum of
   per-source powers.** With several recovered sources in one bin at small angular
   separation the Earth-term columns go near-parallel and per-source powers diverge;
   the total is immune. The design matrices are SVD-truncated at rcond = 1e-6.
6. **The HD variance floor is a Gaussian-approximation quantity.** It comes out at
   σ_A/A = 0.024 while the *realised* RMSE is 0.101 — **4× larger**, because the sky is
   300 heavy-tailed point sources per bin and not a Gaussian field. The M-C column is
   therefore a statement about the noise-limited floor, and the realised RMSE column is
   the one to quote to an observer. They agree in direction and differ in size.
7. **The faint sea is only 8.2 % of the bin.** With γ_eff = 1.56 essentially all of the
   "background" is the unmodelled bright end, whose Γ-projection is an HD curve only in
   ensemble average. That is physically correct on VERNIER-24's audited sky and it is
   why the realised variance is 4× the Gaussian floor — but it means "background" here
   means "everything not yet resolved", which is not what a Gaussian-GWB analysis assumes.
8. **One fixed array across all realisations** (116 pulsars, 100 ns median, 0.40 dex
   heterogeneity, log-uniform 0.3–5 kpc). The ensemble is over universes, not arrays.
9. **Recovery is by array SNR rank, not by running the wheel.** R is a cut on truth,
   so blends, split detections and false recoveries — VERNIER-24/25's dominant failure
   mode, 6.5 blends per 7.7 recoveries — are absent by construction. **Every number here
   assumes the R recovered sources are genuine.**
10. **White noise plus the modelled sky. No intrinsic pulsar red noise, no solar-system
    ephemeris error, no monopole/dipole**, all of which are the real competitors to a CRN
    amplitude and would compete with the M-A bias channel specifically.
11. **Fixed-slope amplitude for the bias/variance columns.** The free two-parameter fit
    is reported for slope only: once the loudest sources are subtracted the residual is
    not a power law, and the A–γ degeneracy then moves A by tens of per cent for reasons
    unrelated to the mechanisms under test.

## One-line verdict

**The loop does sharpen the background, treatment (iii) beats treatment (ii) at every R
and by a widening margin — +18.1 % at R = 4, +26.1 % at R = 8, +38.9 % at R = 16 and 25
on the RMSE of the recovered amplitude, all at p < 1e-4 on a 200-realisation paired
bootstrap, against a scrambled-distance hazard control that buys −2 % to +5 % — but the
mechanism is not the obvious one: M-B, the reassignment of power from background to
foreground, is NEGATIVE at every R (−4.9 % to −39.6 %) because subtracting R sources
projects out 2R coherent directions and eats the background along with them, so the
entire gain lives in the pulsar terms and splits into a bias channel and a variance
channel; M-A is the larger of the two and is the known confound finally quantified —
unmodelled sources inflate the common-red-noise amplitude by +49 % to +140 % in A², an
Earth-only fit removes only two thirds of that and leaves a +10 % to +33 % residual
bias that is entirely pulsar term plus a slope tilt of up to +0.95 in γ, and modelling
the pulsar terms with certified distances drives it to −1.5 % to −4.4 % while a
scrambled distance leaves over half of it standing; M-C, the unpublished arrow, is real
and clean — pulsar terms contribute exactly nothing to the mean of a Hellings-Downs
cross-correlation and sit entirely in its variance normalisation, making them a floor
that is invisible to the statistic it degrades, and certified distances cut that floor
by +23.8 / +27.1 / +29.7 / +31.1 % at R = 4 / 8 / 16 / 25 with ±1 % errors while the
scrambled control moves it the wrong way by −3 % to −6 %, the size being set by a median
pulsar-term power that runs from 97 % to 169 % of the residual background it sits in;
quoted in the requested currency the gain is ~2.5 % of background-amplitude precision
per recovered source at R = 8 and saturates hard — +3.5 / +2.6 / +1.2 / +0.6 % per source
as R goes 4 → 25, with the realised total peaking near +20 % at R = 8-16 before the
projection cost turns it over, and the noise-limited floor continuing to +38.5 % because
it never pays that cost; and the whole result is CONDITIONED on oracle distances that
VERNIER-24/25 say the loop certifies for exactly zero pulsars today, so this prices a
product rather than delivering one.**
