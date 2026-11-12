# JOB 50 — What stops wide-prior pulsars from certifying (Parts B–E)

**Which conflicts were real, and what changes.**
- **Conflict 1 (retardation): real, and it matters.** The booked headline has no retardation. Adding it per pulsar cuts the median share of a decode pair from 0.27 to 0.07 and removes the 0.90 catalogue at the booked amplitude. **The booked 36.27 at 0.9009 does not survive.**
- **Conflict 2 (weak sources): the record is wrong under exact statistics.** With accurate templates, using every candidate instead of 12 raises the wide-prior (K ≥ 60) fraction from 0.25 to 0.54 without retardation, and from 0.035 to 0.24 with it.
- **Conflict 3 (template term): real, and it is what decides wide priors.** As coded (δκ·σ_d) it holds K ≥ 60 at 1–3% in every arm. With the true distance (δκ·L) nothing beyond the prior survives for anyone.
- **Code fact that contradicts the brief's framing:** the booked 601-point distance grid is itself too coarse for wide priors (step 1.5 pc at K ≥ 60 against a 0.25 pc window). Every wide-prior number is given on a fine grid as well.
- **The retardation size depends on a mass convention that two existing scripts disagree on** (see B2). Read the retardation-on numbers as "with Job 41's masses".

**Five-line summary for a reader with no context.**
1. A pulsar's distance is found by lining up the pulsar-term phases of several gravitational-wave sources. The wider the prior on the distance, the more wrong answers (K) must be rejected.
2. Wide-prior pulsars (K ≥ 60) are not stopped by confusion alone. With perfect knowledge of each source they would be decoded 24–54% of the time, against 0.2% by chance.
3. What stops them is **not knowing each source's sky position and frequency well enough**. That error is multiplied by the distance and scrambles the predicted phase.
4. Three further things each cost a factor: the pile-up of chirped pulsar terms, the machine's cut to 12 sources, and a search grid that is too coarse.
5. Ten times better timing helps a little (1% → 7%), and only through better templates.

Everything is paired: same 100 held-out skies (100–199), same 79 red-noise-clean pulsars per sky (7,900 pulsars), same random numbers across arms. Source: `job50_data/analysis.json`, `job50_logs/analysis.log`.

---

## Part B — the corrected arm `p2r_mf`

**Build** (`job50_common.collect50`). For each pulsar: every source's Earth term stays at f_E; every source's pulsar term, the target's included, moves to f_P = f_E(1 + τ/t_c)^(−3/8) with power scaled by (f_P/f_E)^(7/3). τ = L(1 − cos μ)/c uses this pulsar's distance and this source's angle. Sea, noise and decoder are as in `p2_mf`. The comb spacing uses f_P.

**B2 validation.**

| Check | Result |
|---|---|
| τ forced to 0 reproduces `p2_mf` | **Yes, exactly.** Flags equal on 15,800 of 15,800 pulsars; equal to Job 49's `p2_mf` on 13,923 of 13,923 booked records; booked gate gives **36.56 at 0.9083** |
| 232 Earth against 44 pulsar terms above 50 nHz | **Reproduced exactly** (pulsar 41, L = 2105 pc, cos μ ~ U(−0.9, 0.9)) *with the notebook's masses*. My function matches the notebook's `FP_ALL` to 1.4×10⁻¹⁵. |
| Same count with the masses `p2r` uses | **20**, not 44 (mean over 120 pulsars: 37 against 77) |
| Sources drifting more than one bin within the span | 1 of 132,514 decode pairs. Treated as monochromatic. |

**The mass discrepancy.** The notebook retards with source-frame chirp masses (`MCx.append(SRC_MC[ii])`, v4 cell 3). Job 41 retards with redshifted masses (`ts=_TSUN*mc*(1.0+z)/_MSUN`, job41:203), which is the observer-frame formula. `p2r` follows Job 41. So `p2r` retards harder than the code that produced 232/44. I did not run `p2r` with the notebook masses; that is the first sensitivity check to do.

**What retardation does to a decode pair** (20 held-out skies, 13,626 pairs; 10th / 50th / 90th percentile):

| | Off | On |
|---|---|---|
| true share of the bin | 0.094 / 0.269 / 0.600 | 0.011 / **0.071** / 0.274 |
| target's pulsar-term power relative to its Earth term | 1 | 0.10 / 0.39 / 0.88 |
| target's shift, in frequency bins | 0 | 0.6 / 3.6 / 10.5 |

87% of all sources move by more than 0.02 bin.

**Stated limits of the "every candidate" arms.** Each of the ~1,130 candidates gets its own noise draw, although they share about 45 frequency bins. Each counts all the others as confusion. This is an incoherent stack, not a joint fit; a joint fit could do better or worse. Only the 12 booked sources carry explicit phasor sums; the rest use an exact-power Gaussian (Job 49 showed this makes no difference).

---

## Part C — the 12 arms

Fraction decoded correctly, ungated, fine grid, with 5–95% bootstrap intervals over skies. K uses the coarsest tooth; the finest-tooth version is in the last column. "Info" is the median over pulsars of Σ(coherence × ρ)² over the sources used.

### Harness priors (K ≥ 60: 1,524 pulsars; K ≥ 600: 226; prior-only baseline 0.002 / 0.000 / 0.131)

| Retardation | Sources | Template | K ≥ 60 | K ≥ 600 | All | Sources used | Info | Finest K ≥ 60 / ≥ 600 |
|---|---|---|---|---|---|---|---|---|
| off | top 12 | oracle | 0.252 [0.212, 0.291] | 0.133 [0.089, 0.184] | 0.694 | 12 | 12.7 | 0.320 / 0.215 |
| off | top 12 | machine (= `p2_mf`) | 0.017 [0.011, 0.024] | 0.000 | 0.606 | 10 | 9.5 | 0.061 / 0.007 |
| off | top 12 | full | 0.001 | 0.000 | 0.142 | 12 | 0.0 | 0.005 / 0.001 |
| off | every | oracle | **0.543 [0.492, 0.584]** | **0.389 [0.327, 0.456]** | 0.835 | 1127 | 18.1 | 0.598 / 0.515 |
| off | every | machine | 0.028 [0.022, 0.036] | 0.000 | 0.688 | 1127 | 13.7 | 0.113 / 0.015 |
| off | every | full | 0.002 | 0.000 | 0.142 | 1127 | 0.0 | 0.004 / 0.001 |
| on | top 12 | oracle | 0.035 [0.026, 0.046] | 0.004 | 0.402 | 12 | 2.9 | 0.080 / 0.025 |
| on | top 12 | machine (= `p2r_mf`) | 0.005 [0.003, 0.008] | 0.000 | 0.377 | 10 | 2.2 | 0.030 / 0.003 |
| on | top 12 | full | 0.000 | 0.000 | 0.137 | 12 | 0.0 | 0.006 / 0.000 |
| on | every | oracle | **0.236 [0.203, 0.271]** | **0.133 [0.097, 0.172]** | 0.651 | 1127 | 7.1 | 0.329 / 0.215 |
| on | every | machine | 0.010 [0.006, 0.015] | 0.004 | 0.562 | 1127 | 5.4 | 0.065 / 0.009 |
| on | every | full | 0.001 | 0.000 | 0.135 | 1127 | 0.0 | 0.005 / 0.000 |

"All" here means all 79 red-noise-clean pulsars per sky, not only the ones the booked code certifies; that is why `p2_mf` reads 0.606 here and 0.691 in Job 49.

### Census priors (K ≥ 60: 6,578 of 7,900 pulsars; K ≥ 600: 2,475; prior-only 0.001 / 0.000 / 0.006)

Prior widths and distances are drawn as pairs from `job45_t2_parallax_census.csv`.

| Retardation | Sources | Template | K ≥ 60 | K ≥ 600 | All | Sources used | Info |
|---|---|---|---|---|---|---|---|
| on | top 12 | oracle | 0.039 [0.031, 0.047] | 0.013 | 0.080 | 12 | 3.2 |
| on | top 12 | machine | 0.003 [0.002, 0.005] | 0.001 | 0.030 | **2** | 0.2 |
| on | every | oracle | 0.224 [0.192, 0.252] | 0.124 [0.099, 0.149] | 0.277 | 1127 | 7.5 |
| on | every | machine | 0.008 [0.006, 0.010] | 0.001 | 0.050 | 1127 | 0.5 |
| off | every | oracle (for reference) | 0.480 [0.433, 0.524] | 0.387 | 0.523 | 1127 | 18.1 |

With real prior widths, 83% of pulsars have K ≥ 60, and the machine's template term leaves a median of 2 usable sources per pulsar.

### The booked grid, same arms

| Arm (harness) | K ≥ 60 booked grid → fine | K ≥ 600 booked → fine |
|---|---|---|
| off, top 12, oracle | 0.073 → 0.252 | 0.013 → 0.133 |
| off, every, oracle | 0.154 → 0.543 | 0.031 → 0.389 |
| on, every, oracle | 0.110 → 0.236 | 0.009 → 0.133 |
| any machine-template arm | unchanged within 0.01 | unchanged |

For 0.15% of pulsars (priors of several kpc) the fine grid hit its size cap and is coarser than 0.05 pc.

### The three answers

1. **Every candidate instead of 12, oracle templates:** yes, a large rise. K ≥ 60: 0.25 → 0.54 without retardation, 0.035 → 0.24 with it. **With the machine's template term:** barely. 0.017 → 0.028 without retardation, 0.005 → 0.010 with it.
2. **Retardation lowers it,** in every pairing: 0.54 → 0.24 (every, oracle), 0.25 → 0.035 (top 12, oracle), 0.028 → 0.010 (every, machine).
3. **Under the full template term nothing beyond the prior survives.** K ≥ 60: 0.000–0.002 against a prior-only 0.002. All pulsars: 0.135–0.142 against 0.131.

Figure: `job50_fraction_vs_K.png` (fraction correct against K for the 12 arms; left harness priors, right census priors).

### C2 — anchors

Number of the 79 pulsars per sky with K < 1 for a given source:

| Priors | Sources | Median per (sky, source) | Share of sources with ≥ 3 such pulsars |
|---|---|---|---|
| census | booked top 12 | 0 | 0.2% |
| census | every candidate | 0 | 1.0% |
| harness | booked top 12 | 1 | 36% |
| harness | every candidate | 14 | 98% |

Under the census priors almost no source has even three pulsars that could tie down its phase unambiguously. The harness priors are far more generous than the real ones.

---

## Part D — sensitivity

The census-law job scaled only the timing noise (`sig_ns = 100.0/boost`, job29:292), so that is what is scaled here; strains are untouched. Arm: `p2r_mf`, every candidate, harness priors, fine grid.

| n | Template | Candidates per pulsar | K ≥ 60 | K ≥ 600 | All | Info |
|---|---|---|---|---|---|---|
| 1 | oracle | 1161 | 0.236 [0.203, 0.271] | 0.133 | 0.651 | 7.1 |
| 3 | oracle | 1300 (cap) | 0.257 [0.222, 0.294] | 0.177 | 0.673 | 7.6 |
| 10 | oracle | 1300 (cap) | 0.257 [0.222, 0.294] | 0.177 | 0.673 | 7.6 |
| 1 | machine | 1161 | 0.010 [0.006, 0.015] | 0.004 | 0.562 | 5.4 |
| 3 | machine | 1300 (cap) | 0.028 [0.020, 0.036] | 0.000 | 0.599 | 5.9 |
| 10 | machine | 1300 (cap) | **0.072 [0.057, 0.088]** | 0.009 | 0.614 | 6.3 |

**Plainly:** certification of wide-prior pulsars does improve with sensitivity under exact statistics, but only through **better templates**. With oracle templates nothing moves beyond n = 3: the data are confusion-limited and the candidate list is capped at 1300. With the machine's template term K ≥ 60 rises sevenfold, from 1.0% to 7.2%, because the template error shrinks as 1/SNR. More candidates contribute almost nothing (1161 → 1300, and the cap binds).

Not run: the version that scales every strain. The census law did not use it.

---

## Part E — re-book

Gate, fixed in advance: 7-point mean without zero padding; scan 0.5–100%; last train crossing of 0.90; blind on held-out skies; no fallback. Booked certifiable set (13,923 records, 69.6 per sky). Intervals are 5–95% over skies.

| Arm | Yield (psr/sky) at held-out purity | Gate | Catalogue with K < 1 / K < 3 | Correct with K ≥ 10 |
|---|---|---|---|---|
| incumbent | 39.69 [36.95, 42.55] at 0.876 [0.859, 0.893] | 66.1% | 45% / 98% | 0.04 per sky |
| `p2_mf` (booked grid) | 40.04 [37.28, 42.89] at 0.892 [0.878, 0.906] | 65.6% | 45% / 98% | 0.03 per sky |
| `p2_mf` (fine grid) | 39.52 [36.75, 42.34] at 0.895 [0.881, 0.908] | 64.4% | 45% / 98% | 0.02 per sky |
| `p2r_mf`, top 12, machine template | **no 0.90 catalogue** (max train purity 0.861) | — | — | — |
| `p2r_mf`, every candidate, machine template | 17.46 [14.71, 20.57] at 0.881 [0.859, 0.901] | 34.8% | 45% / 99% | 0.00 per sky |

No arm reaches 0.90 on held-out skies at the marginalised amplitude; train crossings generalise to 0.88–0.90.

**Amplitude table, `p2r_mf`, every candidate, machine template (fine grid):**

| A (10⁻¹⁵) | Ungated purity | Result | Catalogue K < 3 |
|---|---|---|---|
| 3.7 | 0.800 | 57.15 [56.49, 57.79] at 0.904 [0.896, 0.913], gate 88% | 86% (2.2 correct per sky with K ≥ 10) |
| 5.0 | 0.728 | 49.24 [48.49, 49.99] at 0.898 [0.890, 0.906], gate 76% | 98% |
| 6.4 | 0.634 | **no 0.90 catalogue** (max train 0.858) | — |
| 8.3 | 0.529 | no 0.90 catalogue (0.760) | — |
| 10.6 | 0.443 | no 0.90 catalogue (0.593) | — |

With top 12 only, `p2r_mf` has no 0.90 catalogue at any amplitude (max train purity 0.804 at A = 3.7).

---

## Booked numbers: changed and confirmed

| Number | Old | New | Status |
|---|---|---|---|
| Headline yield at 0.90 | 36.27 at 0.9009 (booked gate) | **no 0.90 catalogue** with retardation and 12 sources; 17.5 at 0.881 with every candidate | **changes** |
| Amplitude limit | catalogue up to A = 6.4 | up to A = 5.0 (every candidate); none at any A with 12 sources | **changes** |
| `p2_mf` | 36.56 at 0.9083 | reproduced exactly; 39.5–40.0 at 0.89 under the fixed gate | confirmed (no retardation) |
| K ≥ 60 correct, machine as booked | 0.020–0.027 | 0.017 (`p2_mf`), 0.005 (`p2r_mf`), all 79 pulsars | confirmed in kind |
| K ≥ 600 | 0 of 79 | 0–1 of 226 with the machine template; **0.13–0.39 with oracle templates on a fine grid** | the *reason* changes |
| "Only ~4 sources usable; more never helps" | N_eff 4.2 | every candidate doubles to sevenfold the wide-prior fraction when templates are accurate | **does not hold under exact statistics** |
| "Does not improve with sensitivity" | 3.4 per sky at ×1, ×3, ×10 | K ≥ 60 rises 1.0% → 7.2% at ×10, through templates | **changes** |
| Pile-up 232 / 44 | — | reproduced exactly; becomes 232 / 20 with redshifted masses | confirmed, with a caveat |
| Heatmap (47c/47d) | — | untouched by this job; it has no retardation and no template term | confirmed as "no pile-up" |

## Pre-registered predictions

| | Prediction | Verdict |
|---|---|---|
| Q1 | τ = 0 reproduces `p2_mf` | **Holds**, exactly (13,923 / 13,923; 36.56 at 0.9083). |
| Q2 | Retardation lowers the median share and does not raise the wide-prior fraction | **Holds.** Share 0.269 → 0.071; K ≥ 60 falls in every pairing. |
| Q3 | Every candidate + oracle: K ≥ 60 above 0.5 | **Holds without retardation (0.543); fails with it (0.236).** It is not below 0.1 in either case, so the hypothesis is not refuted: weak sources do add up. |
| Q4 | Every candidate + machine template: K ≥ 60 below 0.1 | **Holds.** 0.028 (off), 0.010 (on). |
| Q5 | Full template term: indistinguishable from prior-only at K ≥ 60 | **Holds.** 0.000–0.002 against 0.002. |
| Q6 | ×10 sensitivity, machine template: K ≥ 60 rises but stays below 0.3 | **Holds.** 0.010 → 0.072. |
| Q7 | No arm reaches 0.90 held-out with more than a prior-dominated catalogue | **Holds.** At the marginalised amplitude no arm reaches 0.90 at all. At A = 3.7 the catalogue reaches 0.904 and is 86% K < 3. |

## Things in the prompt that the code shows to be wrong

- "100 held-out skies … K ≥ 600: 0 of 79" was on the booked certifiable set and the booked grid. On all red-noise-clean pulsars there are 226 such pulsars per 100 skies.
- The booked grid cannot decode wide priors (Part A addendum A6). The brief did not list it as a factor.
- The 232/44 count and Job 41's 18.7 come from two different mass conventions.
- "D12" is not in the repo.
- Under census priors the machine's own σ < π cut leaves a median of 2 sources per pulsar, not "about 4".

## What was not done

- `p2r` with the notebook's (source-frame) masses.
- Sensitivity with every strain scaled.
- A joint fit over candidates (the incoherent stack is a stand-in).
- Part F (next report).

## Files and wall time

Added, staged, nothing committed: `JOB50_A.md`, `JOB50_MAIN.md`, `JOB50_STATUS.md`, `job50_partA.py`, `job50_common.py`, `job50_run.py`, `job50_analysis.py`, `job50_fraction_vs_K.png`, `job50_data/*.json`, `job50_logs/*`. Run pickles (`job50_data/run_*.pkl`, `parts/`) are not staged; they regenerate with `job50_run.py`. No booked file was changed.

Wall time: Part A tests 1 min; production 37 min on 22 cores (a first attempt of about 40 min was lost because it saved nothing until the end); analysis 1 min.
