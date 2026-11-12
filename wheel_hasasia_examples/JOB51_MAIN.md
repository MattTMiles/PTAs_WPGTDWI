# JOB 51 — Is confusion the limit, or template accuracy?

## Five-line summary

1. **The oracle ceiling is real.** A proper joint fit over every candidate, with perfect templates, decodes 24% of K ≥ 60 pulsars at today's timing and 27% at 10× better timing (Job 50's stack: 24%). Prediction R1 (above 50%) is wrong.
2. **What sets the ceiling is the gravitational-wave power that is not in the candidate list**, not the candidates confusing each other and not timing noise. Candidates hold 44% of the band's power; of what is left, 93% is other sources and the sea, 7% is timing noise.
3. **Template accuracy has a hard floor.** Against the background, the array signal-to-noise ratio of the booked top 12 sources is 3.9 and does not move when timing noise goes to zero. Reaching a 1 rad template error for the median pair would need 6,900. No pair gets there.
4. **Most "candidates" are not detectable sources.** The candidate list is selected against timing noise only. Against the background, 7 of the 1,161 candidates per sky have an array ratio of 3 or more.
5. **Raw material for self-calibration is thin.** Per sky, 0.5 sources have at least 4 pulsars in which their pulsar term holds more than half its bin; 3.3 sources at a share above 0.2. With real (census) priors and K < 1, none.

## Where the prompt and the code disagree

- **"Candidates" are defined against timing noise only** (`job43c_rebooking.py:491-493`, quoted in Part 1). With the background counted they are mostly undetectable (Part 3). Every "every candidate" arm in Jobs 50 and 51 therefore assumes knowledge of ~1,150 sources per sky that the array could not find.
- **Job 50's stack number agreed with the joint fit partly by accident.** The Job 50 stack gives each candidate its own independent noise draw. The same stack run on bins that are built once and shared gives 2% at K ≥ 60, not 24% (Part 2). The joint fit recovers 24% by subtracting the candidates from each other.
- **The reference value 0.236 uses a marginalised amplitude; the joint fit cannot.** Job 50 scales "counts only" (other sources' power is multiplied by a rate, the target's is not: `job43c_rebooking.py:503`, `res=(...)*rate     # COUNTS ONLY`). A bin built once has no way to scale a source as confusion but not as signal. The joint arms use the fixed amplitude 6.4×10⁻¹⁵ (rate 1); the Job 50 stack rerun in the same setting gives 0.184, and that is the like-for-like comparison.
- **The prompt's 142 rad** (median δκ·L, Job 50 Part A) is for a different set of pairs. On the held-out skies it is 178 rad for booked-top-12 pairs and 87 rad for every-candidate pairs.

## Provenance of every number below

All code is uncommitted (see `JOB51_PROVENANCE.md`); files modified 2026-11-11.

| Numbers | Produced by | Raw output |
|---|---|---|
| Part 1 frames, `arr` | `job50_frames.py` | `job50_logs/frames.log`, `job50_data/frames.json` |
| Part 1 both mass conventions | `job50_common.collect50` (env `JOB50_MASS=z` or `src`) via `job50_run.py`; tables by `job50_analysis.py` | `job50_logs/analysis.log` (redshifted), `job50_logs/analysis_src.log` (source frame) |
| Part 2 | `job51_common.JointHook.psr` called from `collect50`; run by `job51_run.py` (`joint_*`) | `job51_data/run_joint_*.pkl`; `job51_logs/analysis.log` |
| Part 3, Part 4 | `job51_common.StatHook`; `job51_run.py` (`stat_*`, `dec_*`) | `job51_data/run_stat_*.pkl`, `run_dec_*.pkl`; `job51_logs/analysis.log` |
| All tables | `job51_analysis.py` | `job51_data/analysis.json` |

Intervals in square brackets are 95%, from resampling skies. "Per sky" pulsar counts refer to the 79 red-noise-clean pulsars of 120.

**Mass convention:** Parts 2, 3 and 4 use the redshifted chirp mass (the physically consistent one, Part 1). Only Part 1 was run under both.

**Rule 3 (reproduce first).** From Job 50's saved run: harness 0.236 / 0.133 / 0.651; census 0.224 / 0.124 / 0.277. Exact. In addition, `collect50` with the new hooks switched off returns records identical to Job 50's for sky 100 (79 pulsars, 6 keys compared; `job51_test.py`, `job51_logs/test.log`).

**Rule 6 (fine grid).** Step ≤ 0.05 pc everywhere except where the size cap (2²² points) binds: 5 of 7,900 pulsars with harness priors (step up to 41 pc), 0 of 1,580 with census priors.

---

## Part 1 — Loose ends from Job 50

**Frames.** Frequencies are observed-frame; the strain uses the redshifted chirp mass; the mass stored with each source is the source-frame one. Notebook cell 3, lines 9–12 and 100:

```
#   f_obs  : observed GW frequency -- the notebook's frequency axis.
#   strain : h = KH * Mc_z^(5/3) * f_obs^(2/3) / d_L, with REDSHIFTED chirp mass
SRC_AMP = KH*(SRC_MC*(1+SRC_Z))**(5/3)/SRC_DL                  # d_L + redshifted-mass form
```

The notebook's retardation (cell 14, lines 1589–1591) combines observed frequency and observer-time lag with the source-frame mass `mc_src`:

```
_KCH_all = (96/5)*np.pi**(8/3)*(G_SI*mc_src/C_SI**3)**(5/3)
FP_ALL = (f_src[None, :]**(-8/3) + (8/3)*_KCH_all[None, :]*_dt_all)**(-3/8)
```

Job 41 uses the redshifted mass (`job41_displaced_excess.py:206`, `ts=_TSUN*mc*(1.0+z)/_MSUN`). **The redshifted mass is the consistent one.** The notebook understates the chirp rate by (1+z)^(5/3): median 2.20 over its sky, 2.12 for the booked 12 sources per sky.

Tests: `SRC_AMP` equals the redshifted form (true) and not the source-frame form (false); `h_src` rebuilds from `mc_src·(1+z)` (true). Caveat: harness skies carry no mass column; (Mc, z) is recovered by amplitude lookup, exact for 79% of sources, an amplitude-neighbour's values for the rest.

**Retardation arms under both conventions** (harness priors, skies 100–199, 7,900 pulsars, fine grid, marginalised amplitude; fraction decoded correctly):

| Arm, retardation on | K ≥ 60 redshifted | K ≥ 60 source-frame | All, redshifted | All, source-frame |
|---|---|---|---|---|
| 12 sources, oracle | 0.035 | 0.073 | 0.402 | 0.499 |
| 12 sources, template as coded | 0.005 | 0.009 | 0.377 | 0.457 |
| every candidate, oracle | 0.236 | 0.314 | 0.651 | 0.720 |
| every candidate, template as coded | 0.010 | 0.015 | 0.562 | 0.620 |
| every candidate, full template | 0.001 | 0.001 | 0.135 | 0.137 |

Headline catalogue under the fixed (edge-safe) gate, pulsars per sky at purity:

| Arm | Redshifted (physical) | Source-frame (notebook) |
|---|---|---|
| 12 sources, template as coded | none reaches 0.90 (max 0.861) | 2.6 at 0.914 |
| every candidate, template as coded | 17.5 at 0.881 | 36.9 at 0.890 |

The notebook's convention nearly restores the booked count (36.3); the physical one halves it.

**How `arr` is computed: against timing noise only.** `job43c_rebooking.py:490-492`:

```
T=T_yr*YR_S; nt=T_yr*365.25/14.0
rp=(h2/(2*np.pi*f2))*np.sqrt(nt/2)/(100.0*1e-9)
arr=rp*np.sqrt(0.2*120)
```

100 ns white noise, 14-day cadence, no background, no red noise. Test: Part 3's `arr_noise`, built independently from bin powers, equals `arr` to 10⁻⁶ for every candidate (assertion in `StatHook.sky`).

---

## Part 2 — Is the oracle ceiling real?

**Arm `joint_oracle`** (retardation on, redshifted masses, fixed amplitude). For each pulsar the 47 Fourier bins k/T (T = 45 yr, up to 33 nHz) are built once: every source's Earth and pulsar term with amplitude leakage sinc((f − f_bin)T), the sea, one timing-noise draw. Candidate Earth terms are subtracted. The template at trial distance L is the sum of all candidates' pulsar-term phasors in each bin. Statistic as specified in the brief, with P_unmodelled = sea + non-candidate sources + timing noise.

Check that would fail if the fit were wrong: with the unmodelled part switched off, the joint fit decodes 79 of 79 pulsars on sky 100, including all 16 with K ≥ 60 (`job51_logs/test.log`).

Three statistics on the same skies, pulsars, distances and angles:

- **stack, own noise** — Job 50's method (each candidate has its own noise draw);
- **stack, shared bins** — same incoherent stack, but on the bins built once (candidate Earth terms left in);
- **joint fit** — on the same shared bins.

**Harness priors, skies 100–199, 7,900 pulsars (1,516 with K ≥ 60, 223 with K ≥ 600):**

| Timing | Statistic | K ≥ 60 | K ≥ 600 | All |
|---|---|---|---|---|
| ×1 | stack, own noise | 0.184 [0.165, 0.203] | 0.081 [0.045, 0.120] | 0.631 [0.618, 0.643] |
| ×1 | stack, shared bins | 0.020 [0.013, 0.028] | 0.009 [0.000, 0.023] | 0.547 [0.536, 0.558] |
| ×1 | **joint fit** | **0.238 [0.216, 0.260]** | 0.076 [0.037, 0.122] | **0.803 [0.794, 0.812]** |
| ×3 | stack, own noise | 0.214 [0.192, 0.235] | 0.126 [0.080, 0.174] | 0.655 [0.643, 0.666] |
| ×3 | stack, shared bins | 0.026 [0.019, 0.034] | 0.009 [0.000, 0.023] | 0.548 [0.537, 0.559] |
| ×3 | **joint fit** | **0.268 [0.244, 0.291]** | 0.108 [0.060, 0.152] | **0.813 [0.804, 0.822]** |
| ×10 | stack, own noise | 0.214 [0.192, 0.235] | 0.126 [0.080, 0.174] | 0.655 [0.643, 0.666] |
| ×10 | stack, shared bins | 0.026 [0.019, 0.034] | 0.009 [0.000, 0.023] | 0.548 [0.537, 0.559] |
| ×10 | **joint fit** | **0.270 [0.245, 0.294]** | 0.108 [0.060, 0.152] | **0.814 [0.804, 0.823]** |

Job 50 reference (stack, own noise, marginalised amplitude): 0.236 / 0.133 / 0.651.

**Census priors, skies 100–119 only (cost), 1,580 pulsars (1,318 with K ≥ 60, 452 with K ≥ 600):**

| Timing | Statistic | K ≥ 60 | K ≥ 600 | All |
|---|---|---|---|---|
| ×1 | stack, own noise | 0.174 [0.148, 0.200] | 0.064 [0.040, 0.086] | 0.225 [0.195, 0.253] |
| ×1 | stack, shared bins | 0.022 [0.014, 0.031] | 0.004 [0.000, 0.011] | 0.049 [0.039, 0.060] |
| ×1 | **joint fit** | **0.188 [0.162, 0.218]** | 0.080 [0.049, 0.110] | **0.271 [0.241, 0.308]** |
| ×3 | **joint fit** | **0.228 [0.197, 0.259]** | 0.095 [0.070, 0.122] | **0.300 [0.270, 0.332]** |
| ×10 | **joint fit** | **0.229 [0.198, 0.260]** | 0.095 [0.071, 0.121] | **0.301 [0.272, 0.332]** |

Job 50 reference (marginalised amplitude): 0.224 / 0.124 / 0.277 on 100 skies; 0.202 / 0.096 / 0.253 on these 20 skies. Census stacks at ×3 and ×10: own noise 0.178 / 0.077 / 0.233, shared bins 0.022–0.023 at K ≥ 60.

**Share of each bin's gravitational-wave power that the candidates account for** (harness, ×1; 1,127 candidates per pulsar, median):

| Measure | Value |
|---|---|
| median over bins (2–30 nHz) and pulsars | 0.046 |
| median over bins at or below 10 nHz | 0.457 |
| summed over the whole band | 0.444 |
| timing noise as a share of the unmodelled power, median bin | 0.071 (×1), 0.009 (×3), 0.001 (×10) |

Candidates sit almost entirely below about 10 nHz. Above about 12 nHz the bins hold essentially no candidate power.

**One sentence:** with perfect templates the limit is the unmodelled remainder — the joint fit is 12 times above the stack on the same data (0.238 against 0.020), so the candidates' mutual confusion is removable — but that remainder is other gravitational-wave sources and the sea, not timing noise, so the fraction barely rises with timing (0.238 → 0.270).

What else the table says:

- Across all pulsars the joint fit is clearly better than any stack: 0.80 against 0.63.
- The stack with its own noise draw per candidate overstates an incoherent stack by a factor of about nine at K ≥ 60. About 1,100 candidates share 47 bins, so their noise is not independent. Job 50's stack also multiplies confusion by the machine's cleaning factor (0.64–0.73 across the band), which the shared bins do not.

---

## Part 3 — Is there a floor on template accuracy?

Array signal-to-noise ratio of each candidate's Earth term, normalised as `arr` (√(0.2·120) × per-pulsar ratio, averaged over the 79 clean pulsars). Harness priors, skies 100–199, marginalised amplitude, redshifted masses, the ×1 candidate list (1,161 per sky; 10.8 booked sources per sky). Pairs are (pulsar, source): 855 per sky for the booked 12, 91,692 per sky for every candidate.

- `arr_noise`: as coded.
- `arr_bg`: noise per pulsar = timing + other sources' Earth terms + all pulsar terms in the bin + the sea.
- `arr_bg_opt`: noise = timing + pulsar terms in the bin + half the sea.

**Booked top 12:**

| Timing | Ratio | Median | 90th pct | δκ·L 10 / 50 / 90 % (rad) | Pairs < 1 rad | Pairs < π |
|---|---|---|---|---|---|---|
| ×1 | arr_noise | 35.1 | 428 | 6.7 / 178 / 926 | 0.009 | 0.052 |
| ×1 | arr_bg | 3.89 | 7.6 | 393 / 1,810 / 6,210 | 0.000 | 0.000 |
| ×1 | arr_bg_opt | 6.12 | 11.5 | 268 / 1,180 / 3,970 | 0.000 | 0.000 |
| ×10 | arr_noise | 351 | 4,280 | 0.67 / 17.8 / 92.6 | 0.126 | 0.224 |
| ×10 | arr_bg | 3.92 | 7.7 | 387 / 1,770 / 6,130 | 0.000 | 0.000 |
| ×10 | arr_bg_opt | 6.28 | 12.2 | 259 / 1,120 / 3,850 | 0.000 | 0.000 |
| ×100 | arr_noise | 3,510 | 42,800 | 0.07 / 1.78 / 9.3 | 0.383 | 0.642 |
| ×100 | arr_bg | 3.92 | 7.7 | 387 / 1,770 / 6,130 | 0.000 | 0.000 |
| ×100 | arr_bg_opt | 6.28 | 12.2 | 259 / 1,120 / 3,850 | 0.000 | 0.000 |
| zero | arr_noise | infinite | — | 0 | 1.000 | 1.000 |
| zero | arr_bg | 3.92 | 7.7 | 387 / 1,770 / 6,130 | 0.000 | 0.000 |
| zero | arr_bg_opt | 6.28 | 12.2 | 259 / 1,120 / 3,850 | 0.000 | 0.000 |

**Every candidate:**

| Timing | Ratio | Median | 90th pct | δκ·L 10 / 50 / 90 % (rad) | Pairs < 1 rad | Pairs < π |
|---|---|---|---|---|---|---|
| ×1 | arr_noise | 26.7 | 54.3 | 23.7 / 86.6 / 279 | 0.0001 | 0.002 |
| ×10 | arr_noise | 267 | 543 | 2.4 / 8.7 / 27.9 | 0.020 | 0.155 |
| ×100 | arr_noise | 2,670 | 5,430 | 0.24 / 0.87 / 2.8 | 0.560 | 0.923 |
| any | arr_bg | 0.127 | 0.417 | 4,560 / 17,000 / 49,800 | 0.000 | 0.000 |
| any | arr_bg_opt | 0.188 | 0.633 | 3,050 / 11,500 / 34,000 | 0.000 | 0.000 |

(`arr_bg` and `arr_bg_opt` for every candidate are the same to three figures at every timing level.)

**Array ratio needed for the median pair to reach 1 rad:** 6,931 for the booked top 12; 2,431 for every candidate.

**How many candidates are detectable at all** (array ratio ≥ 3, per sky): 1,161 of 1,161 with `arr_noise`; **7.4** with `arr_bg` (7.0 of the 10.8 booked); 12.8 with `arr_bg_opt` (9.8 of the booked).

**One sentence:** as timing improves without limit the template error stops at a floor set by the gravitational-wave power — median 1,770 rad for the booked top 12 (1,120 rad on the optimistic bound) — and the fraction of usable pairs at that floor is zero below π, for both source lists; the floor is already reached at today's 100 ns.

**Part D row rerun** (every candidate, full template term δκ·L_true, retardation on, fine grid, marginalised amplitude, skies 100–199; candidates per pulsar 1,127 at ×1, 1,300 at ×10 and zero; "zero" is timing noise ÷ 10⁶):

| δκ from | Timing | K ≥ 60 | K ≥ 600 | All |
|---|---|---|---|---|
| arr as coded (Job 50) | ×1 | 0.001 [0.000, 0.002] | 0.000 | 0.135 [0.127, 0.144] |
| arr_bg | ×1 | 0.002 [0.000, 0.004] | 0.000 | 0.131 [0.123, 0.138] |
| arr_bg_opt | ×1 | 0.002 [0.000, 0.004] | 0.000 | 0.131 [0.123, 0.138] |
| arr_bg | ×10 | 0.002 [0.000, 0.004] | 0.000 | 0.131 [0.123, 0.138] |
| arr_bg_opt | ×10 | 0.002 [0.000, 0.004] | 0.000 | 0.131 [0.123, 0.138] |
| arr_bg | zero | 0.002 [0.000, 0.004] | 0.000 | 0.131 [0.123, 0.138] |
| arr_bg_opt | zero | 0.002 [0.000, 0.004] | 0.000 | 0.131 [0.123, 0.138] |
| prior only (decode = prior centre) | — | 0.002 | 0.000 | 0.131 |

Every row equals the prior. With the template error tied to the background, no timing precision recovers anything.

---

## Part 4 — Raw material for self-calibration (counts only)

Share = pulsar-term power ÷ (pulsar-term power + all other gravitational-wave power in its bin), retardation on, no cleaning factor, redshifted masses, marginalised amplitude, 100 skies (100–199). Counts are pulsars out of 79 red-noise-clean per sky. Adding timing noise to the denominator changes every entry by at most 0.06.

**Number of pulsars per source in which its pulsar term exceeds the share threshold:**

| Priors | Source list | Share > 0.5: median / 90th / max | Share > 0.2 | Share > 0.1 |
|---|---|---|---|---|
| harness | booked top 12 (10.8 per sky) | 0 / 0 / 79 | 0 / 34 / 79 | 8 / 76 / 79 |
| harness | every candidate (1,161 per sky) | 0 / 0 / 79 | 0 / 0 / 79 | 0 / 0 / 79 |
| census | booked top 12 | 0 / 0 / 77 | 0 / 40 / 79 | 12 / 73 / 79 |
| census | every candidate | 0 / 0 / 77 | 0 / 0 / 79 | 0 / 0 / 79 |

**Number of sources per sky (mean over skies) with at least 4 / 10 / 30 such pulsars:**

| Priors | Pulsars counted | List | Share > 0.5 | Share > 0.2 | Share > 0.1 |
|---|---|---|---|---|---|
| harness | all | booked 12 | 0.48 / 0.25 / 0.14 | 3.27 / 2.29 / 1.18 | 6.48 / 5.00 / 2.97 |
| harness | all | every candidate | 0.48 / 0.25 / 0.14 | 3.37 / 2.35 / 1.20 | 7.48 / 5.81 / 3.35 |
| harness | K < 1 for that source | every candidate | 0.22 / 0.06 / 0.00 | 1.46 / 0.44 / 0.00 | 2.94 / 0.89 / 0.05 |
| harness | K < 10 | every candidate | 0.42 / 0.19 / 0.05 | 3.09 / 1.94 / 0.84 | 7.09 / 5.08 / 2.40 |
| census | all | booked 12 | 0.66 / 0.42 / 0.17 | 3.84 / 2.83 / 1.42 | 7.12 / 5.87 / 3.34 |
| census | all | every candidate | 0.66 / 0.42 / 0.17 | 3.94 / 2.90 / 1.44 | 8.26 / 6.72 / 3.81 |
| census | K < 1 for that source | every candidate | 0.00 / 0.00 / 0.00 | 0.00 / 0.00 / 0.00 | 0.00 / 0.00 / 0.00 |
| census | K < 10 | every candidate | 0.04 / 0.00 / 0.00 | 0.38 / 0.00 / 0.00 | 0.69 / 0.00 / 0.00 |
| census | K < 1 | booked 12 | 0.00 | 0.00 | 0.00 |
| census | K < 10 | booked 12 | 0.04 / 0.00 / 0.00 | 0.38 / 0.00 / 0.00 | 0.69 / 0.00 / 0.00 |

Reading: the sources with usable pulsar terms are almost all among the booked 12; going to every candidate adds about one source per sky at the 0.1 threshold and none at 0.5. A handful of sources per sky (about 3 at share > 0.2) dominate their bin in many pulsars. With census priors, almost none of those pulsars have a narrow enough prior to anchor the phase (largest count for any source: 3 pulsars at K < 1, 8 at K < 10).

---

## Verdicts on the predictions

| | Prediction | Result | Verdict |
|---|---|---|---|
| R1 | `joint_oracle` K ≥ 60 above 0.5 at ×1 | 0.238 [0.216, 0.260]; stack 0.184 same setting, 0.236 Job 50 | **Wrong.** Within 0.05 of the Job 50 stack. But the stated alternative reading is also not what the data show: mutual confusion is removed by the joint fit (shared-bin stack 0.020); the ceiling is the power outside the candidate list |
| R2 | zero timing: median `arr_bg` of booked 12 below 200; under 5% of pairs below 1 rad | 3.92; 0.000 | **Holds**, by a factor of 50 |
| R3 | full template, δκ from `arr_bg_opt`, zero timing: K ≥ 60 below 0.05 | 0.002 [0.000, 0.004] | **Holds** (equal to the prior) |
| R4 | none | counts in Part 4 | — |

## The three questions of §0

- **(a)** The oracle ceiling is real: a joint fit gives 24% (27% at ×10 timing), the same as the stack's figure, and what limits it is gravitational-wave power outside the candidate list, which no timing precision removes.
- **(b)** Yes, there is a floor: `arr` as coded counts timing noise only; with the background included the booked sources sit at an array ratio of 3.9 (6.3 on the optimistic bound) whatever the timing, against 6,900 needed, and no (pulsar, source) pair reaches a template error below π.
- **(c)** Barely: about 0.5 sources per sky dominate their bin (share > 0.5) in at least 4 pulsars, about 3 at a share above 0.2, and with real priors none of those pulsars has K < 1.

## Limits of these results

- **Joint template shape.** Only the phase varies with trial distance; each candidate's retarded frequency and amplitude are held at their true-distance values (as in Job 50). Each κ is moved to the nearest grid value (shift below half a grid step) so that data and template agree exactly.
- **Fixed amplitude** in Part 2 (see top). Parts 3 and 4 use Job 50's marginalised amplitude.
- **Census joint fit on 20 skies**, not 100.
- **`arr_bg` treats the background as independent noise in each pulsar** with the same √(0.2·120) combination as `arr`. The Earth-term background is correlated between pulsars, which is why `arr_bg_opt` is given as the optimistic bound. Neither is a full array analysis.
- **Part 3 statistics use the ×1 candidate list at every timing level**; the decode rows use each level's own list (as Job 50 Part D did).
- **Oracle depth.** The oracle knows the candidates only. A deeper list would raise the ceiling; 56% of the band's power is outside the present one.
- **Not run:** Parts 2–4 under the source-frame mass convention.
