# GEO_geometry_ensemble.md — agent GEO, ACCRE

Untracked working report. Cronus is canonical; this file never feeds back into
`project_progress.md`. No tracked file was edited, nothing was committed or pushed.

- Ensemble: **40 independent isotropic source-sky redraws** (seeds 700000–700039)
- Population: the census cell, **3 loud + 13 faint, N_CW=16**, canonical literature priors
  (`best_distances.txt` + `stagec_anchor_a2.py::LIT_INJECT`)
- Machinery: `stagec_anchor_a2.py` imported and called **verbatim** (`build_amo`, `prep_draw`,
  `make_vmap2`, `scan_all`, `classify_draw`, `classify_full`). Nothing reimplemented.
- Jobs: warm/gate `12452278` (dgx03, 10m07s, PASS); array `12452499_[0-7]` (8×5 draws,
  all COMPLETED 0:0, 5m46s–6m33s each)
- Artifacts: `GEO_results/geo_draw_{000..039}.npz`, `geo_summary.npz`, `geo_ensemble.png`
  (3.9 MB total)

---

## 0. THE ANSWER IN ONE PARAGRAPH

The standing caveat — *"the count is plausibly robust; the names are not"* — is **half right, and
the wrong half was the reassuring one.** The count is **not** 3: across 40 sky draws a realistic
population certifies **4.50 ± 1.48** pulsars (P>0.9), spanning **1–9**; the census's single draw
returned 3, below the ensemble mean. The names are **much** less robust than "not robust" conveys:
**the exact census triple {J0711, J1713, J1909} reproduces in 0 of 40 draws**, all three certify
together in only 6/40, and the union of everything that ever certifies is **18 pulsars**. One name
survives as near-universal (**J1909-3744, 38/40**); one name absent from the census's list
(**J0437-4715, 32/40**) outranks two of the three that are on it. Certification is genuinely
geometry-driven, and the lever is the **fringe-breaking evidence dlnL**, not the trials factor:
stratified within-pulsar, ρ(1−cos μ, dlnL) = **−0.251** while ρ(1−cos μ, K_counted) = **+0.008**.

---

## 1. PROTOCOL, AND WHAT THESE NUMBERS ARE NOT

**Zero-noise (Asimov) data; conditional-at-truth fringe scans.** Identical protocol to the banked
census. These are therefore **CONDITIONAL-ON-TRUTH CEILINGS**, computed at the exact truth source
parameters — *not* numbers achievable by a cold-start search (`project_progress.md` L1022-1034;
B0.2 puts the certification tolerance at ~0.006 deg sky, and P2 is testing whether the flat gap
between the F-stat basin and the needle is crossable at all). Everything below is a statement about
the ceiling's **geometry dependence**, and inherits the ceiling caveat unchanged.

**What was held fixed.** The population's loudness structure, exactly: the same 3 loud
(log10_h = −13.25) + 13 faint (−14.25) strains, the same 16 log10_fgw ∈ [−7.940, −7.554], same
log10_mc / cos_inc / phase0 / psi, the same frozen red noise and GWB, the same 116 pulsar-distance
truths. **Only the 16 source sky positions were redrawn**, isotropically
(cos_gwtheta ~ U(−1,1), gwphi ~ U(0,2π)), one draw per seed. Every draw's seed and positions are
saved in its npz (`seed`, `src_cos_gwtheta`, `src_gwphi`).

**Harness value gate (job 12452278).** Running the *unredrawn* census sky (seed 3000) through the
identical code path reproduces the banked cell exactly:

| | GEO | banked census |
|---|---|---|
| J0711-6830 P_true(lit) | **0.9534** | 0.9534 |
| J1713+0747 P_true(lit) | **0.9887** | 0.9887 |
| J1909-3744 P_true(lit) | **0.9984** | 0.9984 |
| population bayes count | **3** | 3 |
| population flat-ident count | **0** | 0 |
| lnL(truth, zero-noise) | **405413.512739** | 405413.512739 (H4c, spec L236/L637) |

`lnL(truth)` is bit-identical to the H4c value gate. The harness is the census's, not a lookalike.

**Cause-tag audit (required, and it fires zero times).** `classify_full` tags `geom_wide` when a
draw's dL exceeds 5× that draw's median dL — the "source nearly behind a pulsar, (1−cos μ)→0, dL
blows up, K<few by accident, pulsar-term SNR ≈ 0" artifact. Across all **4640 (pulsar, draw)
cells: 0 flagged, 0 `prior_cert` (K=0)**. Honest counts and raw counts are therefore *identical*
here, and nothing was excluded.

This is not luck, it is structural, and it is worth recording. The scanned fringe spacing is
`dL = min over ALL 16 sources` of c/[f(1−cos μ)]. For dL to blow up, **every one of the 16 sources**
must sit nearly behind the pulsar. At N_CW=16 that is effectively impossible. The geometry-wide
artifact is an N_CW ≲ 4 phenomenon — exactly where the census found it (the J0437 onset triple, its
dL = 12.3 pc antipodal draw at N_CW=1). **At the population cell the audit has nothing to remove.**
Reported separately as required, and the separate report is: none.

---

## 2. DELIVERABLE (a) — CERTIFIED-COUNT DISTRIBUTION

Over 40 draws, literature prior column, cause-tag audited (= unaudited, §1):

| criterion | mean ± sd | median | range | histogram (count : draws) |
|---|---|---|---|---|
| **Bayes P_true > 0.9** | **4.50 ± 1.48** | 5 | **1 – 9** | 1:1 2:1 3:10 4:5 5:15 6:6 7:1 9:1 |
| **strict P_true > 0.99** | **1.57 ± 0.98** | 2 | **0 – 4** | 0:5 1:15 2:13 3:6 4:1 |
| flat (dlnL > ln K_counted) | 1.38 ± — | 1 | 0 – 4 | — |

**The question as posed — "is it 3 ± 1 or 0–6?" — is answered: neither.** At P>0.9 it is
**4.5 ± 1.5 with a 1–9 span**; the distribution is broader than 3±1 and shifted above 3. The
census's single draw (3) sits at the **25th percentile** of the ensemble (12/40 draws certify ≤3).
Quoting "3" as the population's certified count understates the typical draw by ~1.5 pulsars and
conceals a factor-of-9 spread across the sky.

At the strict 0.99 bar the count is **1.6 ± 1.0**, and **5/40 draws certify nothing at all**. The
strict bar is where the "0–6" pessimism actually lives.

Note also that the **flat (ln K) bar is no longer identically zero**: it averages 1.4 pulsars per
draw (range 0–4), where the census's one draw gave 0. Even the conservative floor is a
geometry-dependent quantity, and the census drew a low value of it.

---

## 3. DELIVERABLE (b) — PER-PULSAR CERTIFICATION FREQUENCY

Union set (certifies in ≥1 draw at P>0.9): **18 pulsars**. Union at strict 0.99: 9 pulsars.
`med 1−cos μ` is to the *nearest loud source*.

| pulsar | n/40 | freq | strict freq | med P_true | med dL (pc) | med 1−cos μ (loud) | in census set? |
|---|---|---|---|---|---|---|---|
| **J1909-3744** | 38/40 | **0.95** | 0.68 | 0.996 | 0.228 | 0.464 | **YES** |
| **J0437-4715** | 32/40 | **0.80** | 0.35 | 0.973 | 0.230 | 0.398 | — |
| **J1713+0747** | 27/40 | **0.68** | 0.30 | 0.959 | 0.226 | 0.346 | **YES** |
| J1744-1134 | 20/40 | 0.50 | 0.07 | 0.902 | 0.235 | 0.460 | — |
| **J0711-6830** | 16/40 | **0.40** | 0.07 | 0.870 | 0.240 | 0.459 | **YES** |
| J1045-4509 | 12/40 | 0.30 | 0.03 | 0.721 | 0.236 | 0.327 | — |
| J0030+0451 | 8/40 | 0.20 | 0.03 | 0.647 | 0.219 | 0.355 | — |
| J1603-7202 | 4/40 | 0.10 | 0.03 | 0.646 | 0.233 | 0.566 | — |
| J1640+2224 | 4/40 | 0.10 | 0.00 | 0.449 | 0.220 | 0.429 | — |
| J1824-2452A | 4/40 | 0.10 | 0.00 | 0.572 | 0.225 | 0.540 | — |
| J2241-5236 | 3/40 | 0.07 | 0.00 | 0.503 | 0.226 | 0.461 | — |
| J1125-6014 | 3/40 | 0.07 | 0.00 | 0.617 | 0.241 | 0.408 | — |
| J0613-0200 | 2/40 | 0.05 | 0.00 | 0.367 | 0.231 | 0.397 | — |
| B1937+21 | 2/40 | 0.05 | 0.00 | 0.554 | 0.230 | 0.465 | — |
| J2317+1439 | 2/40 | 0.05 | 0.03 | 0.361 | 0.214 | 0.330 | — |
| J1012+5307 | 1/40 | 0.03 | 0.00 | 0.441 | 0.214 | 0.471 | — |
| J2222-0137 | 1/40 | 0.03 | 0.00 | 0.516 | 0.215 | 0.380 | — |
| J1545-4550 | 1/40 | 0.03 | 0.00 | 0.368 | 0.215 | 0.462 | — |

**Answer to the brief's framing: J1909 is an 18/20 pulsar (38/40 = 0.95). J0711 is a 6/20 pulsar
(16/40 = 0.40).** The two live at opposite ends of the same census sentence.

### Set stability against the census's one draw

| statistic | value |
|---|---|
| draws where the certified set **equals** the census set {J0711, J1713, J1909} | **0 / 40** |
| draws where all three census names certify (among others) | 6 / 40 (0.15) |
| draws where ≥1 census name **fails** to certify | **34 / 40** |
| Jaccard(draw set, census set) | 0.384 ± 0.132 |
| draws where J1909-3744 certifies | 38 / 40 |

**The census's named triple is not reproduced by a single one of 40 independent skies.** The
"names are not robust" caveat was correct in direction and badly understated in magnitude: it is
not that the list wobbles at the edges, it is that the list as published is a **measure-zero
outcome** of the sky draw. What survives redrawing is one name (J1909), a rough count, and a
selection function.

### The census draw was unlucky for J0437-4715

J0437 certifies in **32/40** redraws — more often than J1713 (27/40) and twice as often as J0711
(16/40) — yet it is absent from the census's population list, because the seed-3000 sky happened to
be one of its 8/40 failures. Independently corroborating: the census's own equal-strain N_CW=16
cell has J0437 certifying under the literature prior, and A1 gives it the array's smallest K
(K_opt = 11.88, K_lit = 3.07). **A conclusion that "the real seed set is J0711/J1713/J1909" would
have omitted the array's best-measured pulsar on the strength of one sky realisation.**

---

## 4. DELIVERABLE (c) — GEOMETRY CORRELATION / THE SELECTION FUNCTION

4640 (pulsar, draw) cells; 180 certified. `1−cos μ` is to the **nearest loud source**
(= |1 + Ω̂·p̂|, exactly the quantity in `compute_mode_spacing`; dL = c/[f(1−cos μ)]).

Certified cells sit at median 1−cos μ = **0.327**; uncertified at **0.444**.

### The marginal correlation is a trap

| Spearman ρ | value |
|---|---|
| ρ(1−cos μ_loud, certified) — **marginal** | **−0.029** |
| ρ(1−cos μ_loud, P_true) — marginal | −0.362 |
| ρ(dL_pc, certified) — marginal | −0.068 |

Read naively, the first line says geometry does not matter. It is **confounded by pulsar identity**:
whether a pulsar can certify at all is dominated by its σ_EM and timing precision, which are
uncorrelated with the sky draw, and that between-pulsar variance swamps the within-pulsar signal.
Stratifying — rank-normalising 1−cos μ *within each pulsar*, which removes every between-pulsar
difference — recovers it (18 pulsars with 1 ≤ n_cert < 40):

| stratified Spearman ρ (within-pulsar) | value | what it is |
|---|---|---|
| ρ(1−cos μ_loud, **certified**) | **−0.158** | the selection function |
| ρ(1−cos μ_loud, **P_true**) | **−0.323** | continuous version, same sign |
| ρ(1−cos μ_loud, **dlnL**) | **−0.251** | **fringe-BREAKING power** |
| ρ(1−cos μ_loud, **K_counted**) | **+0.008** | trials factor — **nil** |
| ρ(1−cos μ_loud, **dL_pc**) | −0.055 | fringe spacing — nil |

Per-pulsar: **15/18 have negative ρ**, median −0.083. The effect is consistent, not driven by
one or two pulsars.

### Mechanism, measured rather than asserted

At N_CW=16 the loud sources **do not set the fringe spacing**. dL is the minimum over all 16
sources, so it is fixed by whichever source has the largest f(1−cos μ) — generically a faint one.
Hence ρ(1−cos μ_loud, dL) = −0.055 and ρ(1−cos μ_loud, K_counted) = +0.008: the trials factor is
**blind** to where the loud sources are. Median dL is flat at 0.214–0.237 pc across every decile
of 1−cos μ.

What the loud-source geometry controls is **dlnL, the evidence separating the true fringe from the
best wrong one** (ρ = −0.251). Certification improves when a loud source lies *closer to the
pulsar's line of sight*. So the census's central claim — that the population's certifications are
**DATA-DRIVEN, carried by the loud sources** rather than by prior width or fringe-count luck — is
confirmed, and now has a measured geometric selection function attached to it.

### P(certify | 1−cos μ to nearest loud source), deciles

| decile | 1−cos μ range | n | n_cert | **P_cert** | med dL (pc) |
|---|---|---|---|---|---|
| 1 | [0.0001, 0.0735] | 464 | 24 | 0.052 | 0.237 |
| 2 | [0.0735, 0.1583] | 464 | 21 | 0.045 | 0.229 |
| 3 | [0.1583, 0.2381] | 464 | 21 | 0.045 | 0.229 |
| 4 | [0.2381, 0.3327] | 464 | 25 | 0.054 | 0.226 |
| 5 | [0.3327, 0.4377] | 464 | 27 | 0.058 | 0.224 |
| 6 | [0.4377, 0.5437] | 464 | 18 | 0.039 | 0.227 |
| 7 | [0.5437, 0.6627] | 464 | 11 | 0.024 | 0.235 |
| 8 | [0.6627, 0.8162] | 464 | 14 | 0.030 | 0.236 |
| 9 | [0.8162, 1.0563] | 464 | 15 | 0.032 | 0.226 |
| 10 | [1.0563, 1.9049] | 464 | **4** | **0.009** | 0.214 |

The selection function is **flat at P ≈ 0.045–0.058 for 1−cos μ ≲ 0.44** (μ ≲ 56°), then falls
monotonically, reaching **0.009 in the top decile** — a **5–6× suppression**. 1−cos μ > 1 means
cos μ < 0: **once the nearest loud source lies in the opposite sky hemisphere from the pulsar,
certification all but stops.** Note it does *not* rise as 1−cos μ → 0: the near-alignment limit
buys a wide fringe window but kills the pulsar-term response, and the two cancel to a plateau.

---

## 5. FIGURE

`GEO_results/geo_ensemble.png` — three panels:
(a) certified-count histogram, P>0.9 and strict P>0.99, with the census's single draw (3) marked;
(b) per-pulsar certification frequency, sorted, census names in red, strict overlay in black;
(c) the selection function P(certify) vs 1−cos μ to the nearest loud source (log x).

---

## 6. THE SENTENCE FOR THE PAPER

> Across 40 isotropic source-sky draws at fixed population loudness (3 loud + 13 faint, N_CW=16,
> literature priors, zero-noise conditional-at-truth ceilings), a realistic population certifies
> **4.5 ± 1.5 pulsars** (P_true > 0.9; range 1–9; **1.6 ± 1.0** at the strict 0.99 bar, range 0–4);
> the most frequent names are **J1909-3744 (95% of draws), J0437-4715 (80%), J1713+0747 (68%),
> J1744-1134 (50%) and J0711-6830 (40%)**, with a union of 18 pulsars certifying in at least one
> draw and the exact triple reported from a single draw reproduced in **none** of the 40;
> certification frequency correlates with (1−cos μ) to the nearest loud source **negatively and
> entirely through the fringe-breaking evidence** — stratified within-pulsar
> ρ(1−cos μ, dlnL) = −0.25 against ρ(1−cos μ, K_counted) = +0.01 — the selection function being
> flat at P ≈ 0.05 for μ ≲ 56° and suppressed 5–6× to P ≈ 0.01 once the nearest loud source lies in
> the opposite sky hemisphere.

---

## 7. WHAT THIS DOES AND DOES NOT LICENSE

1. **The census's count is a draw from a distribution, not a property of the array.** Any downstream
   text that says "the realistic population certifies 3 pulsars" should say **4.5 ± 1.5**, or
   explicitly say "3 in the fiducial sky draw".
2. **The census's names must not be used as a seed set.** `project_progress.md` L1013-1020 —
   *"the real seed set is ... the ~3 loud-source-broken pulsars that survive a realistic population
   (J0711/J1713/J1909)"* — is a **one-draw statement**. A sequential estimator bootstrapping from
   that literal triple is bootstrapping from a set that never recurs. What is defensible: seed from
   **J1909-3744** (0.95), and treat the seed set as sky-conditional, computed per realisation.
   J0437-4715 (0.80) belongs in any such set and is currently absent from it.
3. **The cause-tag audit is a no-op at N_CW=16 and must not be quietly dropped.** It fires zero
   times here for the structural reason in §1; it will fire again at low N_CW, where the census saw
   it. Keep it.
4. **Ceilings, not achievability.** Nothing here bears on P2 / reachability. If the flat gap proves
   uncrossable, this measures the geometry dependence of an unreachable ceiling.
5. **Not measured here:** noise realisations (zero-noise Asimov throughout), redraws of the loudness
   structure itself (strains and frequencies frozen by construction — this isolates geometry), and
   pulsar-distance-truth redraws. The 1–9 spread is the *geometry-only* spread; adding noise and
   loudness scatter can only widen it.
6. **Softmax-sampling caveat inherited unchanged.** Pulsars with wide windows (K_counted > B = 512)
   are under-sampled in the softmax denominator; per the census this makes P_true a safe
   over-estimate for the far fringes, not an inflated one. J1909 (σ_lit = 3 pc) remains fully
   sampled, which is reassuring for the one name that survives.

---

## 8. EXECUTION RECORD

| item | value |
|---|---|
| lane | `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1` |
| warm/gate job | `12452278`, dgx03, 10m07s, exit 0, gate PASS |
| array | `12452499_[0-7]`, 8 tasks × 5 draws, all COMPLETED 0:0, 332–379 s each |
| walltime requested | 30 min/task (6× headroom; `bf_resolution=1800` granularity) |
| batching | 5 draws per process — build 113 s + D=5 graph materialisation ~115 s (warm) amortised 5× |
| cold D=5 materialisation | 406 s, paid once in the warm job (cache 3583 → 4072 entries) |
| batched scan, D=5 | 23.7 s |
| disk | 3.9 MB total (40 lean npz + summary + png). Posteriors summarised, raw scans discarded. |
| VRAM | fit comfortably; D=5 chosen against the 30 GB squatted-card plan |

**Startup device log, per task (the convention).** 8 tasks walked 8 distinct A100 UUIDs. Three
landed on cards carrying `wut18`'s rogue non-Slurm PyTorch processes (~50 GB each) — precisely the
three UUIDs `1e10409a…`, `563c9580…`, `aa596b0f…` named in HPC_SETUP §7.2, still squatting. Per
§7.3 conclusion 4, contention changes **timing only, never values**; the contended tasks ran
374–379 s vs 332–350 s clean, and their npz are used unreservedly. **The rogue is still there and
still unreported to ACCRE admins.**

`sleep 12` before the first `nvidia-smi` tolerates the previous tenant's CUDA-context teardown
ghost (§7.1). All 8 device lines are the log's first `[GEO]` line, as required.

### Files written (all untracked)

```
GEO_geometry_ensemble.md              this report
GEO_results/geo_draw_{000..039}.npz   per draw: seed, 16 src positions, frozen h/f,
                                      omc_loud (116×3), omc_min_{loud,all}, dL_pc, marg_pc,
                                      P_true/dlnL/K_counted/ident/bayes/bayes_strict/
                                      prior_cert/geom_wide × {feather,script,lit}
GEO_results/geo_summary.npz           consolidated (40×116) arrays
GEO_results/geo_ensemble.png          the figure
hpc_harbor/geo/geo_ensemble.py        driver (imports stagec_anchor_a2 verbatim)
hpc_harbor/geo/geo_warm.sbatch        warm-cache + value-gate job
hpc_harbor/geo/geo_array.sbatch       8×5 array
hpc_harbor/geo/geo_consolidate.py     consolidation + figure
```

Nothing was committed. Nothing was pushed. No tracked file was edited.

---

## 9. FOLLOW-UP (no new draws; consolidation-only re-cut of the banked 40)

### 9.1 Union-18 vs the P1 needle carrier-18 — **distinct sets, equal size, 15/18 shared**

Carrier set reconstructed from `trackB_p1_{coarse,zoom1,patch}.npz` at the truth grid cell
(index 40, sky offset 0.00e+00 in all three): **true-reg = (qmax > 0.5) & (mapk == 0)**, which
returns exactly 18 and matches each file's stored `n_registered_true = 18` and
`registration_count = 18`. All three npz agree on the same 18 names. Pulsar indexing verified
against `census_idx` → {J0711-6830, J1713+0747, J1909-3744}, i.e. the P1 arrays share
`load_pulsars(116)` order with the GEO arrays.

**Answer: neither the same set nor a subset — distinct, with a large intersection.**

| | value |
|---|---|
| union-18 (GEO, certifies in ≥1 of 40 skies) | 18 |
| carrier-18 (P1, registers at truth) | 18 |
| **shared** | **15** |
| union-only | J0613-0200, J1640+2224, J2317+1439 |
| carrier-only | J1446-4701, J1455-3330, J2145-0750 |
| Jaccard | **0.714** |
| census-3 in both | yes |

The equal size (18 = 18) is a **coincidence of two different quantities** and should not be read
as agreement: carrier-18 is a *single-sky* count (the census's seed-3000 geometry, at truth),
while union-18 is the *union over 40 skies* — a strictly more permissive construction that still
lands on the same cardinality. The two disagree on 3 names each way.

The asymmetry is informative. All three **union-only** pulsars are the ensemble's rarest
certifiers (2–4 draws of 40, freq ≤ 0.10): they certify only on favourable geometry and do not
register on the census sky. All three **carrier-only** pulsars — J1446-4701, J1455-3330,
J2145-0750 — register at truth on the census sky yet certify in **0 of 40** draws. So
**registration at truth is necessary but not sufficient for certification**: the carrier set is
the pool of pulsars whose combs co-register, and certification then selects from that pool by
prior width and fringe-breaking margin. Consistent with the census's own framing (registration is
the geometry; certification adds the trials factor and the σ_EM window).

### 9.2 Per union-18 pulsar: mean (1−cos μ to nearest loud source), certifying vs not

Sorted by certification frequency. `carrier?` marks membership in the P1 carrier-18.

| pulsar | n/40 | mean 1−cos μ \| **certifying** | mean 1−cos μ \| **non-certifying** | ratio | carrier? |
|---|---|---|---|---|---|
| J1909-3744 | 38/40 | 0.5340 | 1.0386 | 0.51 | Y |
| J0437-4715 | 32/40 | 0.4439 | 0.5539 | 0.80 | Y |
| J1713+0747 | 27/40 | 0.3725 | 0.8181 | 0.46 | Y |
| J1744-1134 | 20/40 | 0.3238 | 0.7453 | 0.43 | Y |
| J0711-6830 | 16/40 | 0.4279 | 0.5266 | 0.81 | Y |
| J1045-4509 | 12/40 | 0.2112 | 0.4700 | 0.45 | Y |
| J0030+0451 | 8/40 | 0.1791 | 0.5457 | 0.33 | Y |
| J1603-7202 | 4/40 | 0.2493 | 0.5756 | 0.43 | Y |
| J1640+2224 | 4/40 | 0.2561 | 0.5544 | 0.46 | — |
| J1824-2452A | 4/40 | 0.3076 | 0.5826 | 0.53 | Y |
| J2241-5236 | 3/40 | 0.2310 | 0.5846 | 0.40 | Y |
| J1125-6014 | 3/40 | 0.1369 | 0.4709 | 0.29 | Y |
| J0613-0200 | 2/40 | 0.3448 | 0.4903 | 0.70 | — |
| B1937+21 | 2/40 | 0.2513 | 0.5142 | 0.49 | Y |
| J2317+1439 | 2/40 | 0.1017 | 0.4861 | 0.21 | — |
| J1012+5307 | 1/40 | 0.3813 | 0.4970 | 0.77 | Y |
| J2222-0137 | 1/40 | 0.2199 | 0.5089 | 0.43 | Y |
| J1545-4550 | 1/40 | 0.2382 | 0.5237 | 0.45 | Y |
| **pooled (180 / 540 cells)** | — | **0.3802** | **0.5413** | **0.70** | — |

**The ratio is < 1 for all 18 of 18 pulsars**, mean 0.50, with no exception — a stronger and
cleaner statement than the §4 median-based table (which gave 17/18) or the stratified
ρ = −0.158. Every union pulsar certifies on the draws where a loud source sits nearer its line
of sight. The effect is largest for the marginal certifiers (J2317+1439 at 0.21, J1125-6014 at
0.29) and mildest for the two most robust names (J0437 at 0.80, J0711 at 0.81), which is what a
saturating selection function predicts: a pulsar that certifies in 32/40 draws does so across most
of its geometry range, so its certifying and non-certifying draws differ less in 1−cos μ.

Note J1909-3744's non-certifying mean, 1.0386, sits **just past the 1−cos μ = 1 hemisphere
boundary** identified in §4: the only two skies (of 40) that break the array's most reliable
pulsar are the two that put every loud source in the opposite hemisphere from it.
