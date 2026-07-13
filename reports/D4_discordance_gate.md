# D4 — the REALISATION-LEVEL DISCORDANCE GATE: designed, pre-registered, tested, **NOT ADOPTED**

**Agent:** cronus · **Date:** 2026-07-12 · tag `criterion-v2.1` · CPU-only, banked npz only, no GPU,
no new realisations. Code `CW_transition/criterion_v2_1_d4.py` (statistic + value gate + the
pre-registration text, machine-readable) and `CW_transition/run_d4_score.py` (scorer).
Empirical basis: `reports/ignite_bank.npz` (IGNITE Stage-1), the 540 `ig_fnull*` fresh nulls, the
40 `ig_sloop*`/`ig_sloopX*` soft-loop banks, `ig2_floors.npz`, `ig2_purity.npz`. Scores banked to
`reports/d4_score.npz`.

---

## 0. THE ANSWER IN ONE PARAGRAPH

**D4 fails its own pre-registration and is NOT adopted — and the anatomy is the sharpest statement
of the wrong-counterpart hole this campaign has produced.** IGNITE-2 left one live lead: the
co-registration statistic rejects wrong-counterpart *detections* perfectly at the realisation level
((c) = 42/42) even though the per-pulsar veto (D3) destroyed true certifications. D4 promotes that
lead to a gate: flag a realisation whose **detected set** co-registers at a source *other* than the
assumed counterpart, and veto every certification in it. The statistic is `S_det`, the
detected-set consensus displacement significance — chosen on the banked distributions before any
condition was scored, because it is the only aggregate whose true-signal distribution concentrates
at the null value (median **0.0** at both onset cells). **It fails condition (i) in all eight
pre-registered combinations** (2 dk-conventions × 2 thresholds × 2 cells; best catch at ≤10 %
false-flag = **90.3 %** against a 95 % bar), and the one setting that does catch 97.5 % flags
**44 %** of true-signal realisations against a 10 % bar. The two failures have **one** cause:
`S_det` is a `|Δk|` detector, and `|Δk|` is *not* the difference between a right and a wrong
counterpart. The misses are the wrong-counterpart realisations whose noise-lock lands within **±1
fringe** of truth (median max|Δk| = **1**; one at Δk = **0** exactly — the fringes co-register
because they *are* right, while the source is wrong); the false flags are the **36 %** of
true-signal realisations at the lit onset cell whose detected set already contains a wrong fringe.
**The gate is squeezed from both sides by the same quantity, and no threshold on it separates the
two populations.** What D4 *does* deliver is condition (iii): **all three** small-|Δk| scrambled-loop
keepers are flagged — including the J0437-4715 Δk = −4 case that the per-pulsar form missed
(`R_all` = 4.65 vs the 9.21 bar → `S_det` = 55.9). **The D1 hole is closable on every instance we
hold, and no gate that closes it survives condition (ii).** Separately and structurally: D4 is
scored in both a truth-referenced (**oracle**) and a prior-referenced (**implementable**) form, and
the implementable form is **2–4× weaker** — a caveat that applies retroactively to every D3 number.

---

## 1. THE STATISTIC, AND WHY THIS ONE (design justified on the banked distributions, before scoring)

D3 vetoed **per pulsar** with a leave-one-out chi-square `R_a` against a reference built from the
whole array, and failed on (b): above onset the reference `u_R` is itself poisoned by the array's
noise-locked MAP fringes, so a TRUE certification fails concordance *with its own poisoned
reference*. D4 asks the question one level up — **is this whole realisation's detected set
registered to a different source than the one we targeted?** — and vetoes the realisation, not the
pulsar.

    S_det = J^T I^-1 J   over the DETECTED set D
            I = sum_{b in D} g_b g_b^T / s2_b ,   J = sum_{b in D} g_b (2*pi*dk_b) / s2_b
          = chi2(u = 0) - min_u chi2(u)   on D

i.e. the generalised likelihood-ratio for *"the pulsars the data actually registered co-register at
a source **other** than the assumed counterpart."* It is also the **cheapest** aggregate available:
one 2×2 solve per realisation, with **no leave-one-out loop** — strictly cheaper than the D3 form it
replaces.

**The three candidate aggregates, scored on the banked distributions at the two onset cells, and the
choice made before any pre-registered condition was evaluated:**

| aggregate | true signal (must NOT flag) | wrong counterpart (must flag) | pure noise | verdict |
|---|---|---|---|---|
| `max/mean/min/frac` of per-candidate `R_all` | lit: median **min_R = 1.4e4** — even the *most concordant* detection of a true realisation is discordant against the 9.21 bar | 1.4e5–2.6e5 | 1.1e5 | **REJECTED** — inherits D3's poisoned reference in full |
| `S_ref` over **all 116** pulsars | lit median 2.2e4 | 4.0e5–5.0e5 | **4.4e5 — ABOVE the wrong-counterpart population** | **REJECTED** — measures "this fringe field is poisoned", which is true of *every* realisation above onset (IGNITE-2 §1.2) |
| **`S_det`** — consensus of the **detected** set | **median 0.0 at BOTH onset cells** | lit 3.3e4 · vlbi 1.8e2–2.1e2 | 4.1e3 / 6.0e1 | **SELECTED** — the only aggregate whose true-signal distribution concentrates at the null value |

**This inverts D3's variant ranking** (there, `R_all` carried the test and `R_det` was the control
with an 87 % ceiling), and IGNITE-2 §1.4 states exactly why: *"the `R_det` control kills 6/6 — but
for the degenerate reason that every detection is discordant with every other under a scrambled
source, which is exactly the realisation-level (not per-pulsar) signal."* **D4 is that sentence,
made a gate.** The detected set is the *clean* subset — the pulsars the data registered — and
building the reference from it is precisely what removes the poisoning that killed D3.

### 1.1 Two dk conventions — and the ORACLE / IMPLEMENTABLE gap (a new, structural caveat)

The fringe grid is indexed **−256 … 255 about the EM-prior mean** (measured on the bank:
per-pulsar `mean(n_true_grid) ≈ 0`, `std ≈ 150` fringes). So `k = 0` **is** the prior-mean fringe,
and two different statistics can be built from the same bank:

- **D4-ORACLE** — `dk = mapk − n_true_grid`, `s2 = (2π)²[1/6 + (1−q_max)(σ_EM/dL)²]`.
  **D3's exact convention**, reproduced bit-comparably (§2). It is referenced to the **TRUE** fringe,
  which a real analysis **does not know**. Its power is an **upper bound**.
- **D4-OBS** — `dk = mapk`, `s2 = (2π)²[1/6 + (σ_EM/dL)²]`.
  Referenced to the **EM-prior mean** fringe, which a real analysis **does** know. The prior's own
  distance error is then present in every `dk` regardless of how confident the fringe posterior is,
  so the `(1−q_max)` factor is **dropped — forced by the change of reference, not tuned.**

**Measured gap: the implementable form is 2–4× weaker** (catch 25–52 % vs the oracle's 43–97.5 %).
The mechanism is that `σ_EM/dL` is **O(150–800) fringes** in the lit tier: **the EM prior is wide
enough to absorb almost any source displacement, so a wrong counterpart does not look displaced
relative to a prior that was never going to localise the fringe anyway.** Consistent with RING
("only sub-3-pc, VLBI-class σ_d matters") — and indeed D4-OBS is 1.6× stronger at the VLBI cell
(51.6 %) than at the lit cell (32.9 %), the tight priors buying exactly the discrimination RING
predicted.

> ⚠️ **THIS CAVEAT TRAVELS BACKWARD ONTO D3.** Every D3 number — including (a) = 100 % and
> (c) = 42/42 — was computed in the **oracle** convention. The co-registration programme as a whole
> has been evaluated with a truth-referenced fringe origin, and its implementable form is materially
> weaker. **No co-registration number in this repo may be quoted as an achievable power without its
> implementable-form value beside it.**

---

## 2. THE VALUE GATE (passed before any scoring)

The D4 machinery is an independent cronus re-implementation of IGNITE-2's ACCRE-side purity scorer
(whose code is scratch, not in the repo). It is gated against the **banked statistic**, which the
lean-npz convention preserved:

    D4-GATE sig    n= 809   dk mismatch=0   max|dlog10 R_all|=2.6e-13   max|dlog10 R_det|=1.2e-10   PASS
    D4-GATE null   n= 280   dk mismatch=0   max|dlog10 R_all|=6.5e-14   max|dlog10 R_det|=1.1e-11   PASS

All 1 089 banked candidates reproduced. Two conventions were **recovered, not assumed**, and both
were forced by the gate: the phase variance is `(2π)²[1/6 + (1−q_max)(σ_EM/dL)²]` with **unscaled**
levers, and the rank handling projects onto the **reference set's** information subspace
(dof = rank of R — this is where the rank-1 `R_det` control lives). *This is the "bank the
statistic, not the verdict" convention paying for itself for the second time: D3's verdict could not
have been re-examined, but D3's statistic could.*

---

## 3. THE PRE-REGISTERED SCORECARD — **(i) FAILS IN ALL EIGHT COMBINATIONS**

Thresholds fixed before scoring, both reported: **T-CHI2** = the analytic χ²(2, p = 0.01) = 9.2103
bar inherited from D3 (no new dial); **T-CAL** = the (1−α) quantile of `S_det` at a **stated FPR
α = 0.10**, fitted on a **disjoint** true-signal calibration set (IGNITE Stage-1, 50 reals/cell) per
this project's own convention (*a calibrated threshold must state its false-alarm rate*), with
condition (ii) then scored **out of sample** on IGNITE-2's sloop realisations (disjoint seeds).
The χ²(2) null is *provably false* above onset — that is what killed D3 — so the empirical bar is
the honest form. Detection everywhere under the **fresh v2 Gumbel floors** (38.86 ± 1.47 /
7.59 ± 0.48 nat, α = 0.05, N = 150).

| cell | form | threshold | **(i) wrong-counterpart flagged** (bar ≥ 95 %) | **(ii) true-signal flagged, OOS** (bar ≤ 10 %) | adopt? |
|---|---|---|---|---|---|
| (−12.75, 30, lit) | ORACLE | T-CHI2 (9.21) | **77/79 = 97.5 %** ✅ | **4/9 = 44.4 %** ❌ (in-sample 36.1 %) | **NO** |
| (−12.75, 30, lit) | ORACLE | T-CAL (1.26e5) | 34/79 = 43.0 % ❌ | 2/9 = 22.2 % ❌ | **NO** |
| (−12.75, 30, lit) | OBS | T-CHI2 | 20/79 = 25.3 % ❌ | 0/9 = 0.0 % ✅ | **NO** |
| (−12.75, 30, lit) | OBS | T-CAL (4.70) | 26/79 = 32.9 % ❌ | 0/9 = 0.0 % ✅ | **NO** |
| (−13.25, 30, vlbi) | ORACLE | T-CHI2 | 52/62 = 83.9 % ❌ | 0/4 = 0.0 % ✅ | **NO** |
| (−13.25, 30, vlbi) | ORACLE | T-CAL (0.72) | **56/62 = 90.3 %** ❌ (best at ≤10 % FPR) | 0/4 = 0.0 % ✅ | **NO** |
| (−13.25, 30, vlbi) | OBS | T-CHI2 | 17/62 = 27.4 % ❌ | 0/4 = 0.0 % ✅ | **NO** |
| (−13.25, 30, vlbi) | OBS | T-CAL (4.43) | 32/62 = 51.6 % ❌ | 1/4 = 25.0 % ❌ | **NO** |

**ADOPTION: NO.** Adoption required **(i) AND (ii)**. No combination passes both, at either cell;
(ii) must hold at *both* onset cells and the lit cell never clears it in the oracle form. Per the
pre-registration: **no partial adoption, no post-hoc threshold tuning.** The wrong-counterpart hole
(D1) **remains OPEN and stated**, and the 14/50 (fresh-floor) wrong-certification rate continues to
travel beside every above-onset count.

**(iv) — the vacuity control.** 153/160 (lit) and 152/160 (vlbi) pure-noise realisations have **no
detection at all** → **no flag, by construction**; reported as such and *not as evidence*, exactly as
D3's (c) was. Of the 7 / 8 pure-noise realisations that **do** detect, the gate flags 14–100 %
depending on form and bar. That is *desirable* behaviour (it vetoes a pure-noise false certification)
and is recorded as anatomy, not as a condition.

---

## 4. THE ANATOMY — both failures are ONE statement

**`S_det` is a `|Δk|` detector. `|Δk|` is not the difference between a right and a wrong counterpart.**

**The misses (why (i) cannot reach 95 %).** Split the wrong-counterpart realisations by whether the
gate caught them (ORACLE / T-CHI2):

| cell | caught | missed |
|---|---|---|
| lit | 67 reals, median max\|Δk\| among detections = **137** | 2 reals, median max\|Δk\| = **1** |
| vlbi | 42 reals, median max\|Δk\| = **13** | 10 reals, median max\|Δk\| = **1** (one at **Δk = 0**) |

**Every miss is a noise-lock that landed within ±1 fringe of the true fringe.** A wrong fringe one
fringe from truth demands a source displacement of one fringe — which is *inside the co-registration
band by construction*. And the limit case is decisive: **at the VLBI cell one missed realisation has
max|Δk| = 0** — a wrong counterpart whose surviving detection sits on the **true** fringe. The
fringes co-register *because they are right*; the **source** is wrong. **A co-registration statistic
tests the fringes, not the counterpart, and when a wrong counterpart's noise-lock lands on the true
fringe there is nothing discordant to find — at any threshold, in any form.** This is not a tuning
failure. It is a proof that **no co-registration statistic can close the D1 hole in general.**

**The false flags (why (ii) cannot hold at the lit cell).** At the (−12.75, 30, lit) onset cell,
**13 of 36 (36 %)** detecting TRUE-signal realisations have an **impure detected set** — at least one
detection sitting on a wrong fringe (§10.8.3's noise-lock class). The ORACLE/T-CHI2 in-sample
false-flag rate is **36.1 %**. **These are the same number.** The gate is not misbehaving: it flags
realisations whose detected set is internally discordant, and at the lit cell **36 % of true
realisations genuinely have an internally discordant detected set.** At the VLBI cell only **3 of 24
(12 %)** are impure — and there the gate's false-flag rate is **0 %**. The false-flag rate is a
faithful *measurement of the cell's impurity*, and cannot be beaten by a better threshold.

> **THE SCISSORS.** To catch a wrong counterpart the bar must sit near |Δk| ≈ 1. To spare a true
> realisation the bar must sit above the displacement implied by the wrong fringes *already inside
> its own detected set* — which at the lit cell run to |Δk| = 25–395. **D3 failed because the
> reference was poisoned; D4 fails because the population it must protect is itself poisoned.
> Same disease, one level up.** The certification bottleneck is not the veto's design — it is that
> above onset, a true-signal realisation and a wrong-counterpart realisation *contain the same kind
> of object*, and a geometry test cannot tell which of them is the counterpart.

---

## 5. (iii) — THE D1 HOLE'S CLOSURE TEST: **the realisation-level form CATCHES what the per-pulsar form missed, 3/3**

The scrambled-loop **keepers** — the realisations that carry a false certification to the soft loop's
fixed point, the D1 hole in its most concrete form, and the exact events D3's scorecard (d) failed on:

| realisation | kept cert | Δk | D3 per-pulsar `R_all` (bar 9.21) | **D4 `S_det` ORACLE** | **D4 `S_det` OBS** |
|---|---|---|---|---|---|
| `sloopX` vlbi g0 | **J0437-4715** | **−4** | **4.65 → MISSED** (the (d) survivor) | **55.9 → FLAG** | **17.3 → FLAG** |
| `sloopX` vlbi g3 | B1937+21 | +21 | killed | **1 728 → FLAG** | **30.9 → FLAG** |
| `sloopX` lit g2 | J0711-6830 | +231 | killed | **3.2e5 → FLAG** | 5.9 → flagged only at the calibrated bar |

**All three keepers are flagged, in both forms, including the small-|Δk| J0437 case that defeated the
per-pulsar statistic.** The realisation-level aggregate sees what the leave-one-out form could not:
J0437's Δk = −4 is invisible *against a poisoned all-pulsar reference* (`R_all` = 4.65) but glaring
*against the other detected pulsars*, which under a scrambled source demand a wildly different source
(`S_det` = 55.9). **The D1 hole is closable, and the closure is demonstrated on every instance this
campaign holds.**

**And it changes nothing about adoption.** The gate that closes the hole flags 36–44 % of true
realisations at the lit onset cell. **Status of the D1 wrong-counterpart hole: OPEN — now with the
sharpest possible characterisation.** It is not that no statistic sees these events; it is that the
statistic which sees them cannot distinguish them from the impurity that the true-signal population
*already carries* at the only cells where the count turns on. The hole and the impurity are the same
physics (§10.8.3's noise-lock), and any gate keyed on that physics inherits both.

---

## 6. WHAT TRAVELS FORWARD

1. **criterion-v2.1 = criterion-v2 + the D4 rejection, recorded with anatomy.** The operative
   criterion is unchanged: DETECTION `dlnL > max(ln K, floor(h,T,tol))`, CERTIFICATION
   `q_max > 0.9`. **No purity layer, at either level.** D3 (per-pulsar) rejected by IGNITE-2; D4
   (realisation-level) rejected here. Both rejections are by their own pre-registration, with no
   threshold tuning at any point.
2. **The D1 wrong-counterpart hole is OPEN and is now known to be structurally un-closable by
   co-registration**: a noise-lock within ±1 fringe of truth (Δk = 0 in the limit) produces no
   discordance to detect, and a gate tuned to catch anything larger flags the true population's own
   impurity at the same rate.
3. **The ORACLE/IMPLEMENTABLE gap is a new, permanent caveat** on the whole co-registration
   programme, D3 included: every such number was computed against a truth-referenced fringe origin,
   and the implementable form is 2–4× weaker because the EM prior is 150–800 fringes wide.
   **Cross-link, and it is the constructive one:** the gap closes with `σ_d` — D4-OBS is 1.6×
   stronger in the VLBI tier — which is the *same* lever RING identified and the same lever the
   onset map rewards. **Sub-3-pc distances are now doubly load-bearing: they buy detections
   (IGNITE) and they buy wrong-counterpart robustness (D4-OBS).**
4. **The purity number travels beside every above-onset count**, permanently, as D3's rejection
   already required: 14/50 realisations carry a wrong certification at the lit onset cell under the
   fresh floors (23/50 under the retired max-of-10 floors).

## 7. CAVEATS THAT TRAVEL

- **The oracle/implementable caveat above is the governing one** and is stated in every table.
- Condition (ii)'s out-of-sample denominators are **small** (9 and 4 detecting sloop realisations);
  the in-sample calibration-set rates (36/24 realisations) are quoted beside them and agree in sign
  and magnitude. The verdict does not rest on the OOS numbers: (i) fails in all eight combinations,
  and the lit cell's (ii) failure is corroborated by the 36 % detected-set impurity, measured
  independently of the gate.
- The wrong-counterpart population is `nullA`/`nullL` (scrambled-source) + the scrambled soft loop.
  A *mis-positioned* (rather than scrambled) counterpart — the `tol > 0` axis — is **not** tested
  here; D4's behaviour under small counterpart offsets is **untested and open**.
- The 15-realisations-per-cell sky-sampling caveat (±0.2-class) travels from IGNITE-2 onto every
  sloop-derived number here.
- Fresh D2-sized floors exist at the **two onset cells only**; D4 was not scored anywhere else,
  and the other 22 cells' 10-null floors would not support it.
