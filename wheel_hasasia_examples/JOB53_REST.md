# JOB 53 — Parts 4 and 5: observing span, the six-worlds notebook

HEAD: `74cdf93` (branch `MM_playground`). All Job 52 and Job 53 files are uncommitted (staged), written 2026-11-11/12.

## Five-line summary

1. **Span is a clean argument of the collector, with one silent catch: a sky with fewer than two usable booked sources is dropped without notice.** At 15 years that removes 43 of 200 natural skies (the poorer ones); at 25 years 3; at 45 years none.
2. **Shorter spans are much worse.** N = 12 known, all pulsars: 0.022 at 15 yr, 0.032 at 25 yr, 0.046 at 45 yr (0.049 on the same 157 skies). Every candidate: 0.052, 0.097, 0.260.
3. **The reason is the share, not the number of bins.** The median known pulsar term holds 1.1% of what it competes with at 15 yr, 2.6% at 25 yr, 5.5% at 45 yr. Sources detectable against the background: 0.6, 1.9 and 4.7 per sky (2.2, 4.7, 10.2 with optimal combination).
4. **The six-worlds scoreboard quoted in the paper plan is not what the notebook contains.** The notebook stores, and reproduces exactly, 1.000 / 0.965 / 0.515 / 0.000 at K = 8.2. The quoted 1.000 / 0.985 / 0.540 / 0.000 at K = 10.1 is in `ACTIVE_WHEEL_REPORT.md` and describes an earlier version.
5. **The notebook is a self-contained toy** (no machine skies, phases drawn from a Gaussian, the √((1 − s)/s) relation). With the joint fit on complex data the same cases give 1.000 / 0.99 / 0.79 / 0.00; and the last zero is not a noise result: it is one source, a prior centred 1.42 pc from the truth, and nothing else.

## Where the prompt and the code disagree

- **"0.985 with red noise → 0.540 with a background … at K = 10.1"**: the notebook as it stands gives 0.965, 0.515 and K = 8.2 (summary point 4).
- **"Realistic loudness → 0.000"** is at K ≈ 55 for the one source that survives, not K = 10.1 or 8.2: the surviving source happens to have the finest tooth (0.22 pc).
- **Span**: the 45-year array is not only "generous" in bins; nothing can be said at 15 years for a fifth of the skies because the collector skips them.

---

## Part 4 — Observing span

### Is the span a clean argument?

Yes: `collect50(skies, T=45, …)` (`job50_common.py`, uncommitted). Everything that should follow the span does:

| What | Where | Follows T? |
|---|---|---|
| Frequency bins of the fit | `job53_common.py`: `nb = int(3.3e-8*Ts)+1`, `fb = np.arange(1, nb+1)/Ts` | yes: 16, 27, 47 bins |
| Number of observations (14-day cadence) and so every signal-to-noise against timing noise | `nt = T*365.25/14.0` | yes |
| Candidate list | same cut at each span; cap 1,300 | yes: 733, 999, 1,300 (cap binds only at 45 yr) |
| Booked list and the ranking by `arr_bg_opt` | rerun at each span (`stat15_nat`, `stat25_nat`) | yes: 3.3, 6.7, 11.6 booked sources per sky |
| Leakage between bins, the sea per bin | `sc = Ts/T_OBS`, `bsea = …*(T_OBS/Ts)` | yes |
| Census priors | fixed table | no (correct: they are today's errors) |
| Harness priors | `px ∝ sqrt(30/T)`, `pb = 3.0*(T/30.0)**(-2.5)` | yes, strongly; not used here |
| Population draw | drawn on 30-year bins, 0.2–100 nHz | no, and need not |
| Construction (`build_pop15`) | 30-year bins, share defined at 15 yr | not used on natural skies |

**Silent dependence.** `job50_common.py`:

```python
r = inv_rate(f2, h2, P0, sea, T, rate=rt)
if r is None: continue
```

`inv_rate` (`job43c_rebooking.py:506-507`: `good=sg<np.pi; cand,sg=cand[good],sg[good]` / `if len(cand)<2: return None`) returns nothing when fewer than two candidates have a phase scatter below π, and the collector then skips the whole sky. Test (`job53_diag_span.py` → `job53_logs/diag_span.log`): of the 200 natural skies, 43 are missing at 15 yr, 3 at 25 yr, 0 at 45 yr, and `inv_rate` returns `None` for exactly those 43 and 3. They are not empty skies (672 or more candidates pass the timing-noise cut on each). The joint fit does not need the booked list; the skip is inherited from the booked machine.

The skipped skies are the poorer ones: at 45 yr, where all are present, the 157 skies that survive at 15 yr give 0.049 for N = 12 and the other 43 give 0.036. So the 15-year figures below are biased upward.

Two things the machine ignores at every span, which matter more at short spans: the lowest one or two bins are in reality absorbed by the timing-model fit, and red-noise pulsars are simply excluded (79 of 120 kept) rather than modelled.

### Results (natural skies, census priors, reference; fraction right [95% interval])

| | T = 15 yr | T = 25 yr | T = 45 yr |
|---|---|---|---|
| Skies used (of 200) | 157 | 197 | 200 |
| Bins up to 33 nHz | 16 | 27 | 47 |
| Candidates per sky | 733 | 999 | 1,300 |
| Detectable against the background, `arr_bg` ≥ 3, per sky | 0.62 | 1.89 | 4.66 |
| The same, optimal combination, `arr_bg_opt` ≥ 3 | 2.15 | 4.69 | 10.15 |
| β_unmod of an N = 12 pair, 10th / **median** / 90th | 0.002 / **0.011** / 0.044 | 0.004 / **0.026** / 0.091 | 0.008 / **0.055** / 0.180 |
| Bins a known pulsar term moves, median | 0.57 | 1.22 | 2.73 |
| N = 7: all pulsars | 0.016 [0.014, 0.018] | 0.024 [0.022, 0.026] | 0.029 [0.026, 0.032] |
| N = 7: K ≥ 60 / K ≥ 600 | 0.004 / 0.001 | 0.006 / 0.002 | 0.007 / 0.001 |
| **N = 12: all pulsars** | **0.022 [0.020, 0.025]** | **0.032 [0.029, 0.034]** | **0.046 [0.043, 0.049]** |
| N = 12: K ≥ 60 / K ≥ 600 | 0.006 / 0.002 | 0.009 / 0.003 | 0.016 / 0.003 |
| N = 12: right per sky (of 79) | 1.8 | 2.5 | 3.6 |
| Every candidate: all pulsars | 0.052 [0.048, 0.056] | 0.097 [0.093, 0.102] | 0.260 [0.252, 0.269] |
| Every candidate: K ≥ 60 / K ≥ 600 | 0.020 / 0.007 | 0.051 / 0.020 | 0.194 / 0.080 |

Prior only: 0.006 / 0.001 / 0.000 at every span. K classes use the booked list at each span, so the class membership shifts (K ≥ 600: 6,322 of 12,403 pulsars at 15 yr; 4,638 of 15,800 at 45 yr).

**Same skies** (the 157 present at 15 yr), all pulsars:

| N known | 15 yr | 25 yr | 45 yr | 45 ÷ 15 |
|---|---|---|---|---|
| 7 | 0.016 | 0.025 | 0.031 | 1.9 |
| 12 | 0.022 | 0.033 | 0.049 | 2.2 |
| every candidate | 0.052 | 0.100 | 0.270 | 5.2 |

**N = 12 by the tooth count of the known set** (K from its coarsest tooth): K below 10: 0.29 (15 yr), 0.42 (25 yr), 0.49 (45 yr); K 10–100: 0.021, 0.037, 0.070; K above 100: 0.003, 0.005, 0.012. Span helps in every class.

Source: `job53_rest.py` → `job53_logs/rest.log`, `job53_data/rest.json`; runs `span15_nat`, `span25_nat`, `main_nat` in `job53_run.py`.

Reading: at 15 years, with the 12 best sources handed over perfectly, 2% of pulsars are right, against 0.6% from the prior alone, and that on the better four fifths of the skies. A 45-year array is where the 4.6% comes from.

---

## Part 5 — The six-worlds notebook

Source: `one_pulsar_six_worlds.ipynb` (file dated 2026-09-25), read cell by cell and re-executed unchanged in a scratch copy (16 s, no errors). Rebuild: `job53_sixworlds.py` → `job53_logs/sixworlds.log`.

### How the data are generated

Phases drawn from a real Gaussian; no complex sum. Cell 3, `surface`:

```python
eps=kk*Lt+(r.normal(0,ss) if jitter else 0.0)
resid=(kk[:,None]*GRID[None,:]-eps[:,None]+np.pi)%(2*np.pi)-np.pi
S=-(resid**2/(2*ss[:,None]**2)).sum(0)-GRID**2/(2*P["SIGMA_D"]**2)
```

Each source's datum is its true phase plus one Gaussian draw of width σ_φ. The decoder is the wrapped quadratic with the same σ_φ it was generated with, plus a Gaussian prior centred at zero. "Red noise" and "background" are not separate ingredients: each world is one value of σ_φ applied to all eight sources (0.30, 0.60, 1.00 rad).

### Which jitter relation

σ_φ = √((1 − s)/s), without the factor 2. Cell 16: `s_val=0.5; SIG_GWB=float(np.sqrt((1-s_val)/s_val))`. Cell 26 (realistic loudness): `return np.sqrt((Pw.sum()-Pw)/np.maximum(Pw,1e-30)+1.0/np.maximum(rho,1e-9)**2)`, then `keep=sig_struct<np.pi`.

The notebook's own phasor illustration (cell 14) contradicts the label it prints: recomputed with its seeds, the measured phase scatter is 0.245, 0.743, 1.142 rad at s = 0.9, 0.6, 0.3, against 0.333, 0.816, 1.528 from √((1 − s)/s) and 0.236, 0.577, 1.080 from √((1 − s)/(2s)). (That is Job 49's finding again.)

### Constructed or natural skies

Neither. The notebook imports only numpy and matplotlib and reads no file. The eight sources have frequencies drawn uniformly in 15–40 nHz and angles drawn until the coarsest tooth gives K between 8 and 12 (cell 1: `if 8.0 <= 6*P["SIGMA_D"]/lam.max() <= 12.0: break`). "Realistic loudness" is a hard-coded list, cell 26: `RHO=np.array([111.,6.3,2.1,1.5,1.0,1.0,1.0,1.0])`. The toy population in cell 32 feeds the later illustration cells only.

### Do the four numbers reproduce?

| World | Quoted in the paper plan | `ACTIVE_WHEEL_REPORT.md:6438-6442` | Notebook, stored | Notebook, re-executed |
|---|---|---|---|---|
| 1 clean | 1.000 | 1.000 | 1.000 | 1.000 |
| 2 + red noise | 0.985 | 0.985 | **0.965** | 0.965 |
| 3 + background | 0.540 | 0.540 | **0.515** | 0.515 |
| 4 + realistic loudness | 0.000 | 0.000 | 0.000 | 0.000 |
| K | 10.1 | 10.1 | **8.2** | 8.2 |

The notebook reproduces itself exactly (every purity, margin and bar; the regenerated `notebook_scoreboard.csv` is byte-identical). It does not reproduce the quoted numbers: the strings 0.985, 0.540 and 10.1 do not occur in it (test in `job53_sixworlds.py`). The report's table describes an earlier version of the notebook.

### The same four cases with the joint fit

Same eight sources (tooth spacings 0.22 to 1.47 pc; K = 8.2 from the coarsest, 55 from the finest), prior 2.0 pc centred at zero, 4,000 realisations. Data as a complex sum: source phasor plus complex Gaussian with power (1 − s)/s. Decoder: the joint fit on a fine grid, all eight sources kept with their weights.

| World | Notebook | Notebook's method rebuilt here | Joint fit, truth fixed at 1.42 pc | Joint fit, truth drawn from the prior | Notebook's method, truth drawn |
|---|---|---|---|---|---|
| 1 clean (σ_φ 0.30) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 2 + red noise (0.60) | 0.965 | 0.949 | **0.991** | 0.992 | 0.949 |
| 3 + background (1.00; s = 0.5) | 0.515 | 0.570 | **0.791** | 0.780 | 0.514 |
| 4 + realistic loudness | 0.000 | 0.000 | **0.001** | 0.133 | 0.123 |
| 5 real world (4 + red noise) | 0.000 | 0.000 | 0.001 | 0.130 | 0.116 |

Prior only: 0.000 with the truth at 1.42 pc; 0.098 with the truth drawn. (The rebuilt notebook column uses different random numbers from the notebook's 200 realisations; differences of 0.02 to 0.05 are sampling.)

Reading:

- **The middle of the scoreboard is pessimistic.** For a source holding half its bin, the notebook gives 0.52; the joint fit on complex data gives 0.79. Two causes: the Gaussian phase with √((1 − s)/s) overstates the scatter, and the wrapped quadratic is not the best decoder.
- **The final zero is real but is not a statistic.** In world 4 one source holds all the information (shares of the others are below 1%), so there is one comb with teeth every 0.22 pc and a prior centred 1.42 pc from the truth: the answer is the tooth nearest zero every time. With the truth drawn from the prior the same world gives 0.13, against 0.10 from the prior alone.
- **For the paper:** quote the notebook's own numbers (0.965, 0.515, K = 8.2) or the joint-fit column, and describe the last world as "one source, so the prior decides", not as purity collapsing to zero.

---

## Verdict

| | Prediction | Result | Verdict |
|---|---|---|---|
| V8 | At T = 15 yr the N = 12 all-pulsar fraction is lower than at 45 yr by at least a factor 1.5 | 0.022 against 0.046: factor 2.1 (2.2 on the same skies) | **Holds** |

## Limits

- 15-year results use 157 of 200 skies (the better ones).
- Part 5's joint-fit column is this job's rebuild of the notebook's cases, not a change to the notebook.
