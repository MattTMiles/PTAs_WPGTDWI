# JOB 53 — Part 0: what the skies are

**Every booked number in `ACTIVE_WHEEL_REPORT.md` since Job 16, the old booked headline, and every machine result of Jobs 49 to 52 used constructed skies. On natural skies the N = 7 to 12 result moves from 4.2–5.3% to 2.9–4.6% of pulsars right (census priors; 200 natural populations, `JOB53_MAIN.md`; on the first 20 natural skies, used in the tables below, 2.7–4.2%); the every-candidate result does not move (27.1% → 26.0%).**

HEAD: `74cdf93` (branch `MM_playground`). All Job 52 and Job 53 files are uncommitted (staged), written 2026-11-11/12.

## Five-line summary

1. **The machine's skies are the raw population draw with three sources per sky raised** so that each holds 52% of its 30-year bin; everything else in those three bins (about 4,100 sources and the bin's sea) is scaled down by about 0.75 in strain. Nothing else between the draw and the decoder changes a strain or a bin's power in the Job 52 configuration.
2. **There are 50 populations.** Sky j and sky j + 50 share the population and differ in everything else (pulsars, distances, priors, angles, noise). Job 52's census skies (100–119) are 20 distinct populations, so its census intervals stand; its harness intervals treated 100 skies as independent and are up to 16% too narrow.
3. **A natural mode now exists** (`sky="natural"` in `job50_common.collect50`): the raw draw, a new population for every sky. With the switch at "constructed" Job 52 reproduces flag for flag; at "natural" every strain and every bin power equals the raw draw.
4. **Natural skies have less power in the top sources** (top 3: 3.4% of the band against 9.0%; top 12: 9.3% against 19.1%) and **fewer sources detectable against the background** (4.5 against 6.5 per sky; 10.2 against 10.8 with optimal combination). The typical known pulsar term's share is unchanged (median 0.057 against 0.056); the upper tail shrinks (90th percentile 0.18 against 0.23).
5. **A second defect, found while testing: the machine looks up each source's chirp mass by its amplitude, and the lookup returns the wrong catalogue entry for 5% of sources and for all three raised sources.** On natural skies I use each source's own drawn mass; the effect on the decode fractions is 0.01 or less.

## Where the prompt and the code disagree

- **"Strain raised by up to ×4, median ×1.6" is right, but the larger change by count is downward**: on sky 100, 3 sources are raised and 4,141 are lowered (strain × 0.74 to 0.80), with the sea of those three bins.
- **The construction does not touch the band total or any other bin.** The "drain from the sea of other bins" branch never runs, because the brightness floor never binds (the report measured 0 of 150).
- **`j % 50` did not affect Job 52's main (census) result**: skies 100–119 map to populations 0–19, all distinct. It affects the harness runs (skies 100–199) and the old booked headline (200 skies).
- **The chirp-mass lookup** (summary point 5) is not in the prompt.

---

## 0.1 Construction

### The three functions, in full

`active_wheel_demo_v4.ipynb` cell 18, lines 11–15:

```python
def _pop18(j, s15=0.35):
    _sd, _f2, _h2 = SEED50[j][0], SEED50[j][1], SEED50[j][2]
    Pn, sea_n, _b, _o, _r, h_n = build_pop15(
        _f2, _h2, _sea50v[j], s_to_s30(s15), 3, seed=_sd)
    return _sd, _f2, h_n, Pn, sea_n
```

Cell 15, lines 12–16:

```python
def s_to_s30(s15, T_ref=15.0):
    """Dominance is defined at the T_ref element; convert to the T_OBS=30 yr
    element used by build_pop15 (ratio s/(1-s) scales linearly in T)."""
    r = s15/max(1 - s15, 1e-9)*(T_OBS/YR_S)/T_ref
    return r/(1 + r)
```

Cell 14, lines 4153–4199:

```python
_BAND15 = (2e-9, 1e-8)
def build_pop15(f2, h2, sea, s_tgt, n_bright, rho_tgt=4.0, seed=0):
    """Returns (f_l, P_l, sid, sea_new, bright_idx, ok) with the band total held."""
    rng = np.random.default_rng(950000 + seed)
    P = h2**2*f2*T_OBS
    sea = sea.copy()
    bmask = (bin_cent >= _BAND15[0]) & (bin_cent <= _BAND15[1])
    bidx = np.where(bmask)[0]
    binof = np.digitize(f2, bin_edges) - 1
    tot0 = P[(f2 >= _BAND15[0]) & (f2 <= _BAND15[1])].sum() + sea[bidx].sum()
    # rank candidate bins by existing top-source power; take n_bright of them
    cand = []
    for b in bidx:
        kk = np.where(binof == b)[0]
        if len(kk) == 0:
            continue
        cand.append((P[kk].max(), b, kk[np.argmax(P[kk])]))
    cand.sort(reverse=True)
    chosen = cand[:n_bright]
    bright = []
    dP = 0.0
    for _p, b, ktop in chosen:
        f0 = f2[ktop]
        Pneed = h_of_rhoP(rho_tgt, f0)**2*f0*T_OBS
        kk = np.where(binof == b)[0]
        Pel0 = P[kk].sum() + sea[b]
        Ptop = max(Pneed, s_tgt*Pel0)
        Pel = Ptop/s_tgt
        others = kk[kk != ktop]
        Po = P[others].sum() + sea[b]
        scale = ((1 - s_tgt)*Pel)/max(Po, 1e-300)
        P[others] *= scale
        sea[b] *= scale
        P[ktop] = Ptop
        dP += Pel - Pel0
        bright.append(ktop)
    # conserve the band total: drain dP from non-target bins' sea
    free = np.array([b for b in bidx if b not in [c[1] for c in chosen]],
                    dtype=int)
    pool = float(sea[free].sum()) if len(free) else 0.0
    ok = pool >= dP
    if pool > 0:
        sea[free] *= max(0.0, 1.0 - dP/pool)
    tot1 = P[(f2 >= _BAND15[0]) & (f2 <= _BAND15[1])].sum() + sea[bidx].sum()
    h_new = np.sqrt(np.maximum(P, 0)/(f2*T_OBS))
    # all arrays stay in the INPUT source order (bright indexes into it)
    return P, sea, np.array(bright, int), ok, tot1/max(tot0, 1e-300), h_new
```

The call that every job uses, `job50_common.py` (uncommitted; the line is unchanged since the committed version at `74cdf93`, where it is line 115):

```python
sd, f2, h2, P0, sea = _pop18(j % 50, 0.35)
```

### In plain words

- **The argument 0.35** is the share of its bin that the loudest source should hold *if the bin were 1/(15 yr) wide*. `s_to_s30` converts it to the 30-year bins the population is drawn on: the ratio share/(1 − share) doubles, giving **0.5185**.
- **The argument 3** is how many sources per sky are treated this way.
- **Which sources.** Take the 30-year bins between 2 and 10 nHz. For each, find its loudest source. Rank the bins by that source's power and take the top three. So the chosen sources are the three loudest sources in 2–10 nHz that sit in different bins.
- **What is done to them.** The source's power is set to 0.5185 of its bin's current total (discrete sources plus sea). Every other discrete source in that bin, and that bin's sea, is multiplied by one common factor so that together they hold the remaining 0.4815. The bin total is unchanged.
- **`Pneed`** is a brightness floor: the power at which the source's single-pulsar signal-to-noise against white timing noise (100 ns, 780 observations) equals `rho_tgt` = 4. If the floor were higher than 0.5185 of the bin, the source would be raised to the floor and the bin total would grow, with the excess taken from the sea of the other 2–10 nHz bins. The report measured that this never happens (`ACTIVE_WHEEL_REPORT.md:5776-5783`: "floor binds 0 / 150"), and the test below confirms that no other bin changes.
- **Net effect on one sky** (`job53_test.py`, test 3, sky 100 = population seed 3): 3 sources raised (strain × 1.47, 2.10, 1.96, at 3.0, 4.6, 9.1 nHz); 4,141 of 18,727 sources lowered (strain × 0.74 to 0.80); the sea changed in 3 of 94 bins. Across the 20 populations behind skies 100–119 the three boosts range from ×0.98 to ×3.98 in strain, median ×1.58 (`job52_logs/diag_pop.log`).

### Where it was introduced and why

- **Job 15** (`ACTIVE_WHEEL_REPORT.md:1308-1313`, "the s-ladder"): "Mechanism [DEMO-controlled, explicit]: bright-tail concentration at **fixed band hc²** — the element top is set to s·P_el, that element's other sources and its sea are rescaled to carry (1−s), and the brightness top-up is drained from non-target bins' sea." It was built as a dial, to measure how much dominance the scheme needs.
- **Job 16** (`ACTIVE_WHEEL_REPORT.md:1576-1579`, "the money run"): "s₁₅ = 0.35, T = 34 yr … Analytic T_ignite = 18.3(1−s)/s = 34.0 yr." The value 0.35 was chosen because it puts the predicted switch-on of the scheme exactly at the 34-year span of that run.
- The same section records that the natural value is lower (`ACTIVE_WHEEL_REPORT.md:1664`): "**Native s(T=15) = 0.137 (measured, not 0.25)**".
- `_pop18` (cell 18) froze the Job 16 setting as the default population, and Job 42A noted it (`ACTIVE_WHEEL_REPORT.md:5765`): "`_pop18` generates the population for **every booked number**."

So the constructed skies were a demonstration setting (dominance 0.35 where the model gives 0.137) that became the default.

### Every other step between the raw draw and what the decoder sees

| Step | Where | Changes a strain or a bin power? |
|---|---|---|
| Rejection of sources above the CW upper limit, inside the draw | cell 3:217-225 | Yes, in both sky types: a rejected source is redrawn. 1.5 rejections per sky (seeds 3–22), 1.85 per sky (natural skies 100–119) |
| Split between discrete sources and sea | cell 3:181-184, 202-205 | No power change. Sources above one tenth of the deepest detection bar are drawn one by one; the rest of each bin is a closed-form sum (the sea) |
| `build_pop15` | above | Yes: constructed skies only |
| "Counts only" rate | `job43c_rebooking.py:503`: `res=(c.sum()-P0[k]*(T/T_OBS)+sea[b]*(T_OBS/T))*rate     # COUNTS ONLY` | Scales all unmodelled power by `rate`. In Jobs 52 and 53 `rate=1.0` (`job52_run.py`, `BASE = dict(marg=False, rate=1.0, …)`), so it does nothing. The old booked headline drew it per sky (amplitude scatter 0.32 in the log, clipped to 0.2–3.5) |
| Candidate cut and cap | `job50_common.py`: `cand = np.where((arr_all >= 3.0) & (f2 >= 2e-9) & (f2 <= fmax_cand) & (rp_all >= 3.0))[0]`, then the loudest 1,300 | No power change; sets which sources can be "known". The cap binds on natural skies (1,300 candidates) and not on constructed ones (1,175 on average: lowering 4,000 sources pushes some below the cut) |
| Rescaling to the span | `sc = Ts/T_OBS`; `bsea = … *(T_OBS/Ts)` | Units only (30-year bins to 45-year bins) |
| Sea counted twice | `seab = 2.0*np.atleast_1d(bsea(fb))*rt` (`job52_common.py:73`) | Earth half plus pulsar half; same in both sky types |
| Pulsar-term amplitude | `pf_all = (fP_all/f2)**(7.0/3.0)` | Retardation; same in both sky types |
| Cleaning factor (`depro_sigma`) | `job50_common.py`, booked stream and Job 50 arms only | Not used by the joint fit of Jobs 51–53 |
| **Chirp-mass lookup** | `job50_common.py`, `_mcz_tc` | See below. Changes the retardation, not the power |

### The chirp-mass lookup (new finding)

`job50_common.py`, `_mcz_tc`:

```python
amp = h2/np.maximum(f2, 1e-30)**(2.0/3.0)
k = np.clip(np.searchsorted(_CACHE["AS"], amp), 0, len(_CACHE["AS"])-1); i = _CACHE["ordA"][k]
ts = TSUN*NS["SRC_MC"][i]*((1.0+NS["SRC_Z"][i]) if MASS == "z" else 1.0)/MSUN
```

The machine does not carry each source's mass. It recomputes the source's amplitude from its strain and frequency and looks for that amplitude in the catalogue of 2.8 million (mass, redshift, distance) samples. Two things go wrong:

- **Rounding.** The recomputed amplitude differs from the catalogue value in the last digit often enough that the search lands on the neighbouring entry, which is a different binary. Test (`job53_diag_mc.py` → `job53_logs/diag_mc.log`, 20 raw draws, 373,074 sources): wrong entry for **4.9% of sources and 3.3% of candidates**. For those, the chirp time used is wrong by a factor 0.1 to 3.9 (10th to 90th percentile; median 0.58).
- **The raised sources.** A source whose strain was raised by `build_pop15` is looked up at its new amplitude, so it is given the mass of a louder catalogue entry. Sky 100: redshifted chirp mass × 3.29, 1.59, 1.90 for the three raised sources; chirp time × 0.14, 0.46, 0.34. These three are always in the known set, and their pulsar terms were retarded two to seven times too fast in Jobs 50 to 52.

On natural skies I recover each source's own catalogue entry exactly (`job50_common.natural_pop`, by matching mass and distance among the neighbouring entries; every source matches) and use its mass. Effect on natural skies 100–119, reference, all / K ≥ 60 / K ≥ 600:

| N known | amplitude lookup (as Jobs 50–52) | drawn mass (used from here on) |
|---|---|---|
| 7 | 0.028 / 0.010 / 0.000 | 0.027 / 0.006 / 0.000 |
| 12 | 0.041 / 0.019 / 0.002 | 0.042 / 0.018 / 0.002 |
| every candidate | 0.268 / 0.196 / 0.095 | 0.258 / 0.188 / 0.077 |

The every-candidate figures move by about 0.01 (inside the interval); the small-N figures do not move.

## 0.2 Repeats

- **50 populations**: seeds 3 to 52 (`SWEEP_SEEDS = [3] + list(range(4, 16))`, cell 9, then `for _sd in range(16, 53)`, cell 14:1072).
- **Sky j and sky j + 50** use the same sources, strains and sea. They differ in the random stream that draws the pulsars (`rng = np.random.default_rng(seed0+sd+977*j)`): sky positions, distances, priors, the angle between each pulsar and each source, the true distance offsets, and every noise draw. The census pick also depends on j (`np.random.default_rng([50, 7, seed0, j])`).
- **Job 52's intervals** resampled skies as if independent (`job52_analysis.py`, bootstrap over sky number).
  - Census skies 100–119 are populations 0–19: no repeats, intervals correct.
  - Harness skies 100–199 are each population twice. Recomputed with populations resampled: N = 12, all pulsars [0.523, 0.547] against [0.524, 0.547]; every candidate, K ≥ 60 [0.212, 0.264] against [0.216, 0.260]. Widths change by −9% to +16%.
- From here on, intervals resample populations.

## 0.3 Natural mode

- **Switch:** `collect50(..., sky="constructed" | "natural", mass="exact" | "lookup")` in `job50_common.py`. Defaults leave Jobs 49–52 unchanged.
- **Natural sky j** is `draw_population(seed=5000 + j)`, nothing rescaled. Seeds 5000 + j were never used by the notebook (it uses 3 to 52), and every sky has its own. The draw's upper-limit rejection is kept and counted.
- **Pulsars and priors**: the census pick depends only on j, so natural sky j has the same pulsar distances and priors as constructed sky j. Sky angles and noise draws necessarily differ (the source lists differ).

**Tests** (`job53_test.py` → `job53_logs/test.log`):

| Test | Result |
|---|---|
| Constructed, census sky 100: Job 53 fit against Job 52's saved flags, 26 arms (reference at 12 values of N; Arm A at 4; reference at ×10 and zero noise; retardation off) | **79 of 79 pulsars equal in every arm**; largest difference in the fit margin 3 × 10⁻¹⁴ |
| Constructed, census skies 100–119, all / K ≥ 60 / K ≥ 600 | reference N = 12: 0.053 / 0.021 / 0.000; every candidate: 0.271 / 0.188 / 0.080; Arm A, N = 12: 0.098 / 0.054 / 0.011. **Exact** (`job53_skies.py` asserts it) |
| Natural sky 100: strain of each of 18,468 sources against the raw draw | largest difference **0** |
| Natural sky 100: power per bin (discrete plus sea) against the raw sum, 94 bins | largest difference **0** |

## 0.4 What changes

Census priors, skies 100–119, 1,580 pulsars each, 20 distinct populations each. Source: `job53_skies.py` → `job53_logs/skies.log`; fit `job53_common.Joint53`; runs `job53_run.py` (`main_con`, `main_nat`). All uncommitted, 2026-11-12.

**The skies**

| | Constructed | Natural |
|---|---|---|
| Sources per sky detectable against the background, `arr_bg` ≥ 3 | 6.5 (3 to 10) | **4.5** (2 to 7) |
| The same with optimal combination, `arr_bg_opt` ≥ 3 | 10.8 (7 to 17) | **10.2** (4 to 16) |
| Sources rejected by the upper limit, per sky | 1.5 | 1.85 (0 to 4) |
| Share of 2–30 nHz power held by the top 1 source | 0.019 | 0.012 |
| top 3 | 0.090 | **0.034** |
| top 7 | 0.172 | **0.068** |
| top 12 | 0.191 | **0.093** |
| top 30 | 0.238 | 0.155 |
| every candidate | 0.447 | 0.401 |
| Largest β_bin per pulsar, 10th / 50th / 90th percentile | 0.114 / 0.268 / 0.397 | 0.105 / 0.209 / 0.372 |
| β_unmod, N = 12 pairs, 10th / 50th / 90th percentile | 0.009 / 0.056 / 0.229 | 0.010 / **0.057** / **0.178** |
| β_unmod, 99th percentile and largest | 0.454, 0.624 | 0.364, 0.535 |
| Fraction of pairs with β_unmod above 0.3 | 0.062 | 0.030 |

(β_bin: the pulsar term's power over all gravitational-wave power in its bins. β_unmod: its power over itself plus everything unmodelled when 12 sources are known.)

**The decode** (fraction right [95% interval]; prior only: 0.007 / 0.003 / 0.000 constructed, 0.007 / 0.002 / 0.002 natural)

| Arm | N known | Sky | All pulsars | K ≥ 60 | K ≥ 600 | Right per sky (of 79) |
|---|---|---|---|---|---|---|
| reference | 7 | constructed | 0.042 [0.033, 0.051] | 0.014 [0.008, 0.022] | 0.000 | 3.3 |
| reference | 7 | **natural** | **0.027 [0.019, 0.035]** | 0.006 [0.002, 0.011] | 0.000 | 2.1 |
| reference | 12 | constructed | 0.053 [0.046, 0.060] | 0.021 [0.014, 0.030] | 0.000 | 4.2 |
| reference | 12 | **natural** | **0.042 [0.035, 0.051]** | 0.018 [0.014, 0.022] | 0.002 [0.000, 0.007] | 3.4 |
| reference | every candidate | constructed | 0.271 [0.242, 0.305] | 0.188 [0.163, 0.219] | 0.080 [0.050, 0.110] | 21.4 |
| reference | every candidate | **natural** | **0.258 [0.232, 0.283]** | 0.188 [0.165, 0.210] | 0.077 [0.049, 0.105] | 20.4 |
| Arm A | 7 | constructed | 0.066 [0.053, 0.080] | 0.031 [0.024, 0.038] | 0.007 [0.000, 0.014] | 5.2 |
| Arm A | 7 | **natural** | 0.042 [0.032, 0.053] | 0.020 [0.013, 0.028] | 0.005 [0.000, 0.012] | 3.4 |
| Arm A | 12 | constructed | 0.098 [0.084, 0.115] | 0.054 [0.041, 0.068] | 0.011 [0.002, 0.022] | 7.8 |
| Arm A | 12 | **natural** | 0.078 [0.061, 0.094] | 0.041 [0.030, 0.053] | 0.011 [0.004, 0.020] | 6.2 |
| Arm A | every candidate | constructed | 0.509 [0.472, 0.547] | 0.436 [0.395, 0.477] | 0.246 [0.194, 0.305] | 40.2 |
| Arm A | every candidate | **natural** | 0.509 [0.468, 0.548] | 0.450 [0.410, 0.488] | 0.274 [0.218, 0.333] | 40.2 |

Reading:

- **The construction helped the small-N result by a factor 1.3 to 1.6** (N = 7: 0.042 → 0.027; N = 12: 0.053 → 0.042). It did not create it.
- **The every-candidate ceiling is a property of the population, not of the construction** (0.271 → 0.258; Arm A 0.509 → 0.509).
- The conclusions of Job 52 stand on natural skies; the small-N numbers come down.

**Early verdicts (20 skies; final in `JOB53_MAIN.md`):** V1 holds so far (N = 12: 0.042, inside 0.015–0.045; `arr_bg_opt` ≥ 3: 10.2 per sky, **not** fewer than 8, so that half is wrong). V2 holds (0.258 against 0.271).

## What I will run next unless told otherwise

- **Natural skies: 200 distinct populations** (skies 100–299, seeds 5100–5299), each with its own drawn masses. One sky costs about 4 minutes of one core, so 200 fit well inside the two hours.
- Constructed skies stay at 100–119 for the comparison curve.
