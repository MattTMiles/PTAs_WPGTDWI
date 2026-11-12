# VERIFY-1 — documentation verification sweep for the paper rebuild

**Read-only on the repo.** One compute item: re-reading banked `.npz` files (no
simulation). cronus, 2026-08-06. Two outputs staged, not committed:
`reports/VERIFY1_documentation_sweep.md` (this file) and
`reports/prong2_extracted_values.md`.

**Rule applied throughout:** verbatim quote or extracted number + `file:line`. Where the
answer is not on disk it is marked **NOT RESOLVABLE** and nothing is reconstructed from
memory.

---

# ACT 2

## A1. V23's data accounting — the "4 N_psr data" enumeration and the Part-1 construction

### (a) The 4 numbers per pulsar

The enumeration is stated in the source, not the report. `hpc_harbor/vernier/vernier15.py:10-13`:

> ```
> COUNTING PREDICTION
>   Per bin each pulsar supplies one complex Fourier coefficient = 2 real numbers.
>   Earth terms alone -> 2 N_psr measurements.
>   Because f_P != f_E, each pulsar's PULSAR term lands in its OWN separate bin, so
>   with known distances it is an extra independent complex datum -> 4 N_psr.
> ```

So the four reals per pulsar are:

| # | datum | where it comes from |
|---|---|---|
| 1 | Re(Earth-term coefficient) | `data_vector` first loop, `vernier15.py:77-78` |
| 2 | Im(Earth-term coefficient) | same |
| 3 | Re(pulsar-term coefficient) | `data_vector` second loop, `vernier15.py:80-82` |
| 4 | Im(pulsar-term coefficient) | same |

Confirmed structurally: `data_vector` returns `2·n_psr` complex entries when `taus is not
None`, and the Jacobian splits each into real+imag —
`J = np.zeros((2 * base.size, ncol))` with `np.concatenate([d.real, d.imag])`
(`vernier23.py:42`, `vernier23.py:50`). `2 × 2 n_psr = 4 n_psr` real rows.

V23's own restatement, `VERNIER23_toy.md:26-27` (verbatim):

> **Marginalising the N_psr distances costs N_psr degrees of freedom out of 4N_psr data, leaving
> 3N_psr for sources.**

and the three counting cases, `vernier23.py:13-17` (verbatim):

> ```
>   data = 4 N_psr (Earth block + pulsar-term block)
>   (a) Earth only      : M*6 <= 2 N        -> M_max = 2N/6
>   (b) distances known : M*7 <= 4 N        -> M_max = 4N/7
>   (c) marginalised    : M*7 + N <= 4 N    -> M_max = 3N/7
> ```

Per-source parameter counts: `NPAR_E = 6` (`th, ph, A, psi, cos iota, phase`) and
`NPAR_P = 7` (`... + f`), `vernier15.py:31-32`. The reason `f` is free only in the
pulsar-term case is stated at `vernier15.py:26-30`: within one bin the Earth term has no
frequency dependence, so including `f` "makes the Fisher singular by construction".

### (b) The exact Part-1 measurement construction

| aspect | answer | provenance |
|---|---|---|
| **how many sources** | M is *scanned*; the report displays M ∈ {8, 12, 17, 20} as table rows and returns M_max ∈ {8, 17, 20, 16}. The **actual M grid is NOT RESOLVABLE** — see (c). | `VERNIER23_toy.md:15-21`; `reports/vernier23_part1.npz` |
| **SNRs — equal or laddered** | **Equal, and in fact undefined.** Every source is drawn with amplitude literally `1.0`: `srcs.append([np.arccos(...), ..., 1.0, ...])`, `vernier15.py:108`. There is **no noise realisation at all** — the Fisher is `J.T @ J` with unit σ (`vernier15.py:99`, `vernier23.py:55`), so this is a pure rank/identifiability calculation, not an SNR-dependent one. `VERNIER23_toy.md:111-112` (verbatim): *"**40 pulsars, one bin, no noise realisation** — Part 1 is a noiseless rank calculation throughout."* |
| **single bin or several** | **Single bin for the Earth terms** — all M sources share one `f = 1e-8` Hz (10 nHz), the `f` argument threaded through `draw`/`fisher`/`data_vector` (`vernier15.py:102`, `vernier23.py:58`). The **pulsar terms each land in their own separate bin** by construction (`vernier15.py:12-13`). |
| **grid** | **There is no parameter grid.** Sources are random draws: sky isotropic (`np.arccos(rng.uniform(-1,1))`, `rng.uniform(0, 2π)`), `psi ~ U(0, π)`, `cos ι ~ U(−1, 1)`, `phase ~ U(0, 2π)` — `vernier15.py:106-109`. Pulsar distances `L ~ log-uniform(300, 5000) pc` (`vernier15.py:104`, redrawn at `vernier23.py:61`). Pulsar sky directions isotropic via `V.iso_dirs` (`vernier15.py:103`). |
| **n_psr** | **40** — but this appears **only in the report prose** (`VERNIER23_toy.md:111`), never in code or npz, because the driver is absent. |
| **prior widths** | σ_L = 10 pc and 100 pc, entering as `Fdd = F[ns:, ns:] + np.eye(nd) / sigma_L_pc ** 2` (`vernier23.py:64`). |
| **seeds / realisations per M** | **NOT RESOLVABLE** — no seed count is stored anywhere. |

### (c) The gap that limits everything in A1/A2

`hpc_harbor/vernier/vernier23.py` is **74 lines and contains no `__main__`, no driver, no
argv handling, and no npz write** (verified: `grep -n "__main__\|def part\|argv"` → 0
matches). It defines exactly two functions, `fisher_with_distances` and
`rank_marginalised`. A repo-wide search for any other file referencing `vernier23` or
`vernier23_part1` returns **nothing**.

`reports/vernier23_part1.npz` is **212 bytes** and holds exactly one array:

```
res  (4,) int64  = array([ 8, 17, 20, 16])
```

Nothing else. **The M grid, the seed count, the per-M "full / —" table body, and the
Earth-only and known-distance code paths are all off-disk.** The four scalars map
unambiguously onto the report's `M_max` row (Earth 8 / known 17 / marginalised-10pc 20 /
marginalised-100pc 16).

**Verdict: QUALIFY.** The counting argument and the 4-per-pulsar enumeration are fully
sourced and quotable. The Part-1 *construction* is quotable for everything above except
`n_psr = 40` (report prose only) and the M grid + seed count (unrecoverable). Any paper
sentence implying a scanned grid with stated resolution must be dropped or re-run.

---

## A2. The ordering anomaly — known 17 vs marginalised 17–20

### (a) First: "17–20" is not a measured interval

`reports/vernier23_part1.npz` stores **20**, a single integer, for the 10 pc marginalised
column. The report writes it as **"17–20"** (`VERNIER23_toy.md:21`, `:24`, `:118`). There
is **no 17 on disk for that column** — the lower end of the quoted range is the
*known-distance* column's value, re-used as a hedge. So one of the three candidate causes
you named — *range-vs-point reporting artefact* — **is present and confirmed**, but it is
a presentational overlay, not the anomaly itself. The anomaly on disk is the bare
**20 > 17**.

### (b) It is not a different estimator between columns

Both columns use the identical rank test: normalise the Fisher to correlation form by its
own diagonal, take `eigvalsh`, count eigenvalues above `tol * ev.max()` with
`tol = 1e-10`.

`vernier15.py:137-141` (known / Earth-only path, `rank_test`):
```python
d = np.sqrt(np.clip(np.diag(F), 1e-300, None))
F = F / np.outer(d, d)
ev = np.linalg.eigvalsh(F)
ev = ev[ev > 0]
rank = int(np.sum(ev > tol * ev.max())) if ev.size else 0
```

`vernier23.py:70-73` (marginalised path, `rank_marginalised`):
```python
d = np.sqrt(np.clip(np.diag(Fm), 1e-300, None))
Fm = Fm / np.outer(d, d)
ev = np.linalg.eigvalsh(Fm); ev = ev[ev > 0]
rank = int(np.sum(ev > tol * ev.max())) if ev.size else 0
```

Same estimator, same tolerance. **That candidate is ruled out.**

### (c) The named cause: the relative-tolerance denominator moves between columns

In exact arithmetic the ordering is **impossible**. The marginalised block is the Schur
complement `Fm = Fss − Fsd Fdd⁻¹ Fsdᵀ` (`vernier23.py:67`), with `Fdd` PSD, so
`Fm ⪯ Fss` in the Loewner order; hence `null(Fss) ⊆ null(Fm)` and
`rank(Fm) ≤ rank(Fss)`. Marginalisation cannot raise the rank. Full stop.

What breaks is that **neither column measures `rank(Fss)` or `rank(Fm)`** — each measures
a *thresholded* rank of a *differently rescaled* matrix:

- known column thresholds `D_ss^{−1/2} Fss D_ss^{−1/2}` at `1e-10 · λ_max(that matrix)`;
- marginalised column thresholds `D_m^{−1/2} Fm D_m^{−1/2}` at `1e-10 · λ_max(that matrix)`.

`D_ss ≠ D_m`, and the cutoff is *relative to each matrix's own largest eigenvalue*.
Subtracting the Schur term preferentially removes the best-determined (distance-coupled)
directions, which drives the normalised matrix toward the identity and **lowers
`λ_max`** — which **lowers the absolute cutoff** and lets more near-null eigenvalues be
counted as non-zero. The marginalised column is therefore graded against a *softer bar*
than the known column. This is the mechanism, and it is derivable from the code above
without re-running anything.

**Corroborating evidence, same file, no extra compute:** the two columns fail their own
counting predictions *in opposite directions*. With N = 40, `VERNIER23_toy.md:23-24`
predicts known `4N/7 = 23` and marginalised `3N/7 = 17`. Measured: known **17**
(undershoots by 6), marginalised **20** (overshoots by 3). A consistently-calibrated rank
threshold cannot undershoot one prediction and overshoot the other; a threshold whose
stringency floats with each matrix's own spectrum can, and does.

**Second contributing cause, also confirmed from code:** the columns are not run on the
same random realisation. `rank_marginalised` calls `D.draw(...)` — which internally draws
its own `L` at `vernier15.py:104` — and then **immediately discards it and redraws**
`L = 10 ** rng.uniform(np.log10(300.0), np.log10(5000.0), n_psr)` at `vernier23.py:61`.
So the geometry consumed from the shared `rng` stream differs between the two paths, and
with no stored seed count each `M_max` is, as far as disk can show, a **single
realisation**. Grid quantization cannot be assessed at all (the M grid is off-disk, §A1c)
— note only that `M_max = 16` for the 100 pc column is **not one of the table's displayed
rows** (8/12/17/20), so the scan was finer than the published table.

### (d) Verdict on the "known distances" column

**EXCISE the ordering; QUALIFY the column.**

The number 17 is not itself corrupt — it is produced by the same estimator as everything
else. What is not quotable is **the comparison**: "marginalised (17–20) ≥ known (17)" is a
thresholding artefact of a floating relative tolerance across two differently-scaled
matrices on two different draws, and it asserts something the data-processing inequality
forbids. The defensible reading of the four banked integers is the *weak* one, which is
also what V23's own verdict line actually needs:

> `VERNIER23_toy.md:116-118` — *"marginalising the pulsar distances over their existing priors
> **preserves the separability gain almost intact** (8 → 17–20 against 17 with distances known)"*

"Almost intact / indistinguishable" survives. "17–20" as an interval does not. **Quote it
as: Earth-only 8; with distances known and with distances marginalised, both ≈ 17–20 and
not separated by this estimator.** Do not print the two columns as distinct values, and do
not print marginalised above known. If the paper needs the columns separated, the run has
to be redone with a fixed absolute tolerance and a shared seed — which is a new
simulation, out of scope here.

---

## A3. Re-extracted [FIG-READ] values

**Delivered:** `reports/prong2_extracted_values.md` (staged).

**One scoping correction, stated up front.** No `[FIG-READ]` marker exists anywhere on
disk (repo-wide search across `*.md`, `*.py`, `*.txt`: zero hits), so the marker list
itself lives only in the external draft. I extracted every value the item names.

**Second scoping correction, material:**

1. **`CW_transition/prong2_results.npz` contains no 15-psr prototype.** All three
   prefixes carry `n_psr = 116` (`ps_n_psr`, `tot_n_psr`, `pop_n_psr`). The
   `diag_T10/T15/T20` labels are **Tspan in years**, not pulsar counts
   (`diag_T15_b3-20_T = 4.7336e8 s = 15 yr`). The prototype run was overwritten by the
   scaled production run. Its only surviving record is prose at
   `project_progress.md:141` (verbatim): *"`marg/cond` falls monotonically: ~0.95 (N=1) →
   ~0.53 (N=50) → ~0.18 (N=96)"*. **The 15-psr R(N) at sampled N is NOT RESOLVABLE.**
   The N_psr = 15 *knee* is resolvable: `stagec_p2c_results.npz` gives
   `n_psr = [15, 30, 60, 116, 200]`, `knee_over_Nstar = [6.47, 12.92, 26.73, 52.75, 91.45]`
   (matching `project_progress.md:1340`).
2. **A.2 and A.3 are not in `prong2_results.npz`** — they are in
   `prong2_coherence.npz` and `prong2_Ntc_map.npz` respectively. All requested values were
   extracted from those.

Headline extractions (full tables in the staged file):

| requested | extracted | key |
|---|---|---|
| R(N), 116-psr, both modes | 16-point table, `ps` ≡ `tot` to **2.804e-12** | `ps_*`, `tot_*`, `mode_indep_maxdev` |
| A.2 coherence crossing at R = 0.5 | **τ_p/t_c = 0.014417** (linear), **0.004389** (saturating), shift **0.304443×** | `x_half_linear`, `x_half_saturating`, `shape_shift` |
| A.3 raw max residual | **−0.330595**, at N/N* = **12.441**, Y = **0.2683** | `max_factor_dev`, recomputed location |
| A.3 commuting peak | **−0.201856**, at N/N* = **27.228**, Y = **0.1389** | `max_factor_dev_commuting` |
| A.3 commuting min (positive lobe) | **+0.195634**, at N/N* = 0.704, Y = 0.5179 | recomputed |

**Two caption-level corrections found while extracting** (both detailed in the staged file):

- *"weak positive lobe"* (`project_progress.md:304-305`) is true of the **raw** map
  (+0.012) but **not** of the commuting one (**+0.196**, i.e. 97 % of the negative peak's
  magnitude). The commuting reduction did not only shrink the negative lobe — it
  **manufactured a positive lobe 16× larger than the raw one**. **QUALIFY** — the COMPOUND
  verdict is untouched (the interior sign structure at N/N* ≳ 3 is unchanged), but the
  adjective must go.
- *"falls 0.99 → 0 across τ_p/t_c = 1e-3 → 1e3"* (`project_progress.md:252-253`) mixes the
  `t_c → ∞` gate value (**0.9931**, `ref_ratio`) with a grid endpoint. The grid's own first
  point is **0.926**. **QUALIFY.**

Also flagged: **three different seed counts** are in play and are cross-quoted in the docs
as if one run — Stage A production **30** (`ps_n_seeds`), Stage A.2 **24**
(`prong2_coherence.npz:n_seeds`), Stage A.3 **15** (`prong2_Ntc_map.npz:n_seeds`).

**Verdict: QUOTABLE** for the 116-psr production sweep, A.2 and A.3, from the staged table.
**NOT RESOLVABLE** for any 15-psr R(N).

---

## A4. The onset dip-and-recovery — diagnosis

**The described feature does not exist in any data on disk, and the doc it would come from
says the opposite.**

**(i) It is not in the 116-psr production curve.** `N* = 8.047`. Measured R at the
straddling grid points is **0.98907** (N = 6) and **0.98613** (N = 10) — nowhere near
0.74. R first reaches ≈ 0.74 at N ≈ 230, i.e. ~29 N*, and from there declines
monotonically. There is no 0.74 → 0.85 excursion anywhere in the 16-point series.

**(ii) The prototype's own record says monotonic.** `project_progress.md:141` (verbatim):

> `marg/cond` falls monotonically: ~0.95 (N=1) → ~0.53 (N=50) → ~0.18 (N=96).

A dip to 0.74 at N ≈ 8 followed by recovery to 0.85 is flatly inconsistent with that
sentence, and the prototype's underlying array is gone (§A3).

**(iii) What non-monotonicity actually exists, and what it is.** The production `ps` curve
has exactly two upward steps, and both are small and inside the seed band:

| step | ΔR | seed band at that N (`(cond_hi − cond_lo)/cond`) |
|---|---:|---:|
| N = 2 → 3 | **+0.0096** | 0.72 |
| N = 25 → 40 | **+0.0106** | 0.20 – 0.27 |

Full diff series: `[−0.0045, +0.0096, −0.0039, −0.0026, −0.0029, −0.0156, −0.0249,
+0.0106, −0.0390, −0.0435, −0.0679, −0.0984, −0.1766, −0.2377, −0.2088]`.

Both positive steps are **~1 % excursions sitting inside a 20–72 % realisation band at
30 seeds** — i.e. **seed scatter**, at roughly 1/20th of the band width, and an order of
magnitude smaller than the described 0.11 feature. They are not a bin-edge or window
artefact (the `rcond` sweep is a null: all three of 1e-8/1e-10/1e-12 give knee = **422.636**,
`rcond_knee_spread = 0.0`), and they are not real structure (they do not survive the
population arm — `pop` puts its upward steps at N = 16 and N = 251 instead, different N).

**Verdict: EXCISE.** The dip-and-recovery is not caption-safe. It cannot be sourced to any
banked array; the one document that describes the prototype states the curve is monotonic;
and the only non-monotonicity that does exist on disk is ~1 %, seed-scatter-sized, and
lands at different N in different power modes. If the feature is genuinely visible in the
prototype figure, the figure predates the surviving data and the claim would need the
prototype re-run to be quotable — a new simulation, out of scope.

---

# ACT 3

## B1. The criterion lineage, the v2.2 floor-validity wording, and the ln K combination rule

### (a) Lineage, one line each with provenance

All four are annotated git tags. Note that **v2.2's substance is not in its commit or tag
message** — both are one-liners; it lives in `project_progress.md §10.17`.

| step | tag → commit | date | what changed, one line |
|---|---|---|---|
| **v1** | `criterion-v1` → `db6ff05` | 2026-07-12 | Introduced the **absolute floor** as a third layer: `DETECTION dlnL > max(ln K_counted, 9.01 nat)`, the 9.01 **fitted to 27 banked nulls** as the smallest floor zeroing all null certifications — because FORGE §9's two-layer gate still let the null fire from small K (J0437, `ln K = 1.39`). |
| **v2** | `criterion-v2` → `fb277db` | 2026-07-12 | **Floor became a function**, `floor(h, T, tol)`, and the *estimator* changed: the constant 9.01 was a max-of-N order statistic with no fixed FPR and N-independent scatter (`sd(max_N) = 8.91/8.68/8.79/8.74` at N = 10/30/100/1000, **flat**), replaced by a **Gumbel (1−α) quantile at α = 0.05**. Plus D1 (counterpart-matched `fN` null family adopted) and D3 (purity layer **pre-registered, not adopted**). |
| **v2.1** | `criterion-v2.1` → `6bec3d6` | 2026-07-12 | **Onset RETRACTED**, the **soft loop adopted** as spec §3's reference implementation (superseding B1.3 DAMPED / IGNITE CASCADE as hard-lock-only verdicts), and **both purity layers rejected** — D3 by its own pre-registration, D4 (new) in all 8 pre-registered combinations. |
| **v2.2** | `criterion-v2.2` → `db2075a` | 2026-07-13 | Added the **zero-fraction floor-validity gate** bounding the D2 estimator's domain; **onset re-cut to N_onset = 59**; e = 0.3 single-member switch refuted; clock not shared. |

Tag subject lines, verbatim:

- v1: `criterion-v1: three-layer certification criterion, fitted floor, campaign consolidation`
- v2: `criterion-v2: null family, floor-as-function, purity layer pre-registered`
- v2.1: `criterion-v2.1: soft loop adopted, onset retracted, both purity layers rejected`
- v2.2 (annotated tag body, the **entire** body): `criterion-v2.2: floor-validity gate; onset re-cut to 59; e=0.3 switch refuted; clock not shared`
  — its **commit** message body is only `floor adopted, onset recut to 59`.

The scope statement, `project_progress.md:3146-3149` (verbatim):

> **criterion-v2.2 = criterion-v2.1 + the zero-fraction floor-validity gate + the corrected banks.**
> The three-layer criterion is **unchanged**; the *estimator* of layer 2 is now bounded, and its
> validity domain is stated. **No purity layer exists** (D3 and D4 both rejected, and the D1 hole is
> structurally un-closable by co-registration).

### (b) The v2.2 floor-validity wording, exact

`project_progress.md:3041-3046` (verbatim, including the emphasis as written):

> ### **CONVENTION ADOPTED (criterion-v2.2), BINDING:**
> ### **The D2 Gumbel floor is valid ONLY where the nullN zero-fraction is ≲ 20 %. Above that, quote the empirical (1−α) quantile with a BOOTSTRAP error, and bank the zero-fraction beside it. THE ZERO-FRACTION IS A REQUIRED COLUMN, NOT A CAVEAT.**
>
> The onset test itself is **unchanged** — `ONSET` iff the count at **floor + its own error** exceeds
> 1. Only the floor and the error change; above the gate the error is the bootstrap sd of the
> empirical quantile, so the onset test is made against a bootstrap error exactly as required.

Restated in the master at `STORY.md:1317-1318` and `STORY.md:26-28`; the failure direction
is at `project_progress.md:3036-3037` (verbatim): *"detection is `dlnL > max(ln K, floor)`,
so a floor that is too low is **too permissive** — it lets pure-noise offenders through."*

**Note for the paper:** the 20 % number is **empirical and unexplained**. `STORY.md:2761`
lists it as open work: *"**NEW: derive the floor's VALIDITY DOMAIN from first principles** —
S6.5 bounds the D2 estimator empirically, at a zero-fraction of 20 %, and nothing explains
why that is the number"*.

### (c) The combination rule for ln K_counted

**There is no combination rule — the trials term is applied per pulsar, element-wise, and
is never summed into the detection bar.**

The bar, `CW_transition/trackB_estimator_spec.md:1218-1223` (verbatim):

> For pulsar `a`, with `dlnL_a` = the likelihood-only fringe gap (best minus runner-up peak,
> prior-free), `K_counted,a` = the counted candidate fringes, `q_max,a` = the E-step fringe
> posterior's modal mass:
>
>     DETECTION      dlnL_a > max( ln K_counted,a , DLNL_FLOOR )      DLNL_FLOOR = 9.01 nat
>     CERTIFICATION  q_max,a > 0.9   (strict: > 0.99)   applied ONLY within detections

and `trackB_estimator_spec.md:1225-1226` (verbatim):

> - **Layer 1, `ln K` — the trials factor.** The fringe-breaking evidence must beat the number of
>   fringes it chose among. Relative, per-pulsar. This is GEO's `ident` / "flat count".

The implementation, `CW_transition/trackB_criterion.py:72-74`:

```python
def detect(dlnL, lnK, floor):
    return dlnL > np.maximum(lnK, floor)
```

`np.maximum` is element-wise over the pulsar axis. **Per-pulsar counts never combine**;
the array-level count is the number of pulsars independently clearing their own bar.

`K_sum` exists but is a **diagnostic aggregate only** — it never enters the bar. What it
tracks is how the *per-pulsar* `K_counted` grows with the source list, via a shared
minimum: `reports/SPARK_launch_criterion.md:131-133` (verbatim):

> **lnK growth per added source** (exact, CPU — `dL` is the MIN over the source list of the mode
> spacing (`forge_b1.apply_geometry:91-94`), so a longer list can only shrink `dL` → more fringes
> in the prior window → `K_counted` grows):

| n_src | K_sum | lnK_med | dL_med |
|---|---|---|---|
| 1 | 348 330 | 7.164 | 7.795e-04 |
| 3 (the loud) | 581 900 | 8.183 | 3.814e-04 |
| 4 (+1 recruited) | 622 702 | 8.291 | 3.642e-04 |
| 16 (the full census) | 993 132 | 8.761 | 2.167e-04 |

> **Growing the list 3 → 16 costs +0.578 nat** on the median pulsar's trials bar (K_sum ×1.707;
> the first recruit alone costs +0.108 nat).

**Verdict: QUOTABLE.** Lineage, floor-validity wording and the combination rule are all
fully sourced. Two caveats travel with them: v2.2's commit/tag messages carry no substance
(cite `project_progress.md §10.17` instead), and the 20 % gate is empirical with no
derivation.

---

## B2. EMBER / basin engagement

**The item's premise is inverted by the source.** The EMBER table does not find that the
predictor is basin engagement — it finds that **engagement is not significant and MOTION
is the boundary**. Quoting `reports/EMBER_offtruth_ladder.md:354-359` in full, the
contingency table on all 112 scrambled realisations (target = manufacturing, ΔN > 0):

| predictor | sensitivity | specificity | Fisher p (1-sided) | verdict |
|---|---|---|---|---|
| ENGAGEMENT C0 ≥ 1 | 0.60 (3/5) | 0.60 | **0.333 (n.s.)** | misses the 2 disengaged manufacturers |
| **MOTION** (moved, step > 0) | **1.00 (5/5)** | 0.72 | **0.002** | **in-sample NECESSARY** (0/5 manufacturers had zero motion) |
| JOINT engaged **AND** mobile | 0.60 | 0.87 | 0.024 | **NOT the boundary** — sensitivity collapses to the engaged subset |
| JOINT engaged **OR** mobile | 1.00 | 0.48 | 0.044 | catches all, but no better than motion alone and far less specific |

`reports/EMBER_offtruth_ladder.md:361-366` (verbatim):

> Joint logistic (standardised, L2 = 1 for the quasi-separation): **β_motion = +0.357 > β_engagement
> = +0.308** — **motion dominates.** Necessary/sufficient: engagement is **neither** necessary (2/5
> disengaged) **nor** sufficient (43/46 engaged starts were safe); motion is **necessary in-sample**
> (all 5 moved) but **not sufficient** (27/32 mobile starts were safe). **The boundary is MOTION —
> the loop manufacturing a false cert requires it to take a real accepted M-step; the hypothesised
> "engaged AND mobile" joint is refuted (sensitivity 0.60).**

and the corrected statement, `:373-374` (verbatim):

> **The corrected two-condition boundary is therefore
> (wrong counterpart) × (motion) — not (wrong counterpart) × (engagement).**

### The quotable rates

**Denominators, exact** (`:325-327`): **112 scrambled realisations**; *"No signal sits at
the assumed counterpart, so **every certification here is false**: `C0` is seed-inherited,
**ΔN is loop-made**. **5 of 112 manufactured** (grew the false-cert count, ΔN > 0)."*

| class | count | rate | provenance |
|---|---|---|---|
| manufactured (ΔN > 0), all starts | **5 / 112** | 4.5 % | `:327` |
| engaged starts that stayed clean | **43 / 46** | 93.5 % | `:363` |
| engaged starts that manufactured | **3 / 46** | 6.5 % | derived from `:356`, `:363` |
| disengaged starts that manufactured | **2 / 66** | 3.0 % | derived (5 − 3 of 112 − 46) |
| mobile starts that stayed clean | **27 / 32** | 84.4 % | `:363` |
| mobile starts that manufactured | **5 / 32** | 15.6 % | `:363` |
| motionless starts that manufactured | **0** | 0 % | `:357` |
| self-cleaning (ΔN < 0) | seed 9510200, **ΔN = −4** | — | `:376-378` |

**By start class (arm), signal side** — `:308-312`, verbatim table:

| arm | what it is | ΔN across the 9 cells | loop-grown wrong certs | accepted motion |
|---|---|---|---|---|
| **map** (b) | honest Earth-term MAP (~0.74–1.58 dex cold in mc) | **0.000 in all 9** | **0** | ~0 (inert) |
| **truth** | on-truth, already engaged | **+0.000 … +0.417** (repair in 4 cells) | 0 (and **cleaned 1**) | 0.8–2.8 steps |
| **fix** (c) | truth + 1e-3 dex on mc | **0.000 in all 9** | 0 | ~0 (inert) |

with the mechanism at `:304-306` (verbatim): *"the honest start is **disengaged and
motionless** — in-basin fraction 0.00–0.42, median accepted M-step ≈ 0 (`steps` column
0.00 in 7 of 9 map cells) — so there is nothing for the additive gain to add."*

The D1 framing, `:437-441` (verbatim): *"Across the 324 signal realisations the loop
carries **90 inherited wrong certs** ... arm totals inherited 8 / 28 / 107 for truth / map
/ fix. It **grows 0** of them and **cleans 1**"*.

**Scope that must travel with any quote** (`:380-384`, verbatim):

> **Scope (travels with §4).** N = 112, **5 events** — small-N, one-sided Fisher, logistic under
> L2 shrinkage. Motion is a *process* quantity: its "predictive" edge partly reflects that it sits
> closer to the manufacturing mechanism than a t = 0 flag can. The claim is not "motion causes
> manufacturing at rate X"; it is the **ordering** — motion separates the 5 events cleanly where
> engagement does not, and every manufacturer moved.

**Verdict: QUALIFY.** The table is quotable exactly as printed, but **the headline must be
rewritten**: the predictor is **motion, not basin engagement** — engagement is n.s. at
p = 0.333. And it rests on **5 events**; the claim is an ordering, not a rate. Any paper
sentence of the form "the predictor is basin engagement" is contradicted by its own source.

---

## B3. The e = 0.3 switch refutation

### (a) The v2.2-era wording

The canonical statement, `reports/RECUT_floors.md:127-129` (verbatim):

> ### **THE HEADLINE: "a single e = 0.3 member switches the count on" is REFUTED.**
> **lit collapses from 1.57 to 0.70 — below the bar. vlbi survives at 1.03, but only marginally:
> it does not clear its own floor's bootstrap error (0.60).**

The evidence table, `reports/RECUT_floors.md:131-136`:

| the switch-on cells | zero-f | floor (Gumbel → **adopted**) | ratio | count (old → **re-cut**) | @floor+bs | margin vs >1 | verdict |
|---|---|---|---|---|---|---|---|
| **m1 e=0.3 lit** | 0.73 | 7.39 ± 0.63 → **11.30 ± 1.02** | 1.53× (6.2σ) | 1.57 → **0.70** | 0.43 | **−0.30** | **below — REFUTED** |
| **m1 e=0.3 vlbi** | 0.48 | 10.58 ± 0.83 → **10.78 ± 1.54** | 1.02× (0.2σ) | 1.13 → **1.03** | 0.60 | **+0.03** | **MARGINAL — not confirmed** |
| m1 e=0.5 lit | 0.81 | 4.61 ± 0.41 → **8.58 ± 0.87** | 1.86× | 3.60 → **3.13** | 2.70 | +2.13 | **CONFIRMED** |
| m1 e=0.5 vlbi | 0.51 | 9.55 ± 0.76 → **9.87 ± 0.91** | 1.03× | 2.43 → **2.27** | 1.73 | +1.27 | **CONFIRMED** |

The corrected, externally-quotable replacement, `reports/RECUT_floors.md:145-150` (verbatim):

> ### **With ONE eccentric member: the switch-on is at e = 0.5.**
> ### **With TWO OR MORE: it is at e = 0.3 (confirmed, both tiers).**
> **The docs' interim binding to e = 0.5 STANDS — and it was the right call.**
>
> The threshold is **not a property of eccentricity alone**; it depends on how many eccentric
> members carry it. That is the corrected, quotable statement.

Restated in the masters at `project_progress.md:3162-3163` and `STORY.md:1610`; and in
`reports/EMBER_offtruth_ladder.md:433-435` (verbatim): *"the chorus onset is **e = 0.5
single-member (C1)** or **e = 0.3 two-member (C2, `m2e03`)** — **never** e = 0.3
single-member (that reading is refuted)."*

### (b) The refutation is itself contested — and this is post-v2.2

`LOTTERY_tier2.md:92-97` (verbatim):

> ### **THE CONSEQUENTIAL FINDING (for the CHORUS / RECUT owners).**
> **RECUT's "single e = 0.3, lit, REFUTED" — called "the single most expensive consequence of the
> floor defect in the whole repo" — sits on the floor's own error knife-edge. It is a below/CONFIRMED
> coin-flip at the ~1-nat level, NOT a settled refutation.** m2 e=0.30 (well above the bar) is solid;
> the e = 0.3 *single-member* result is not. This does not touch CHORUS's structural results (clock
> not shared, the 10× lever, the two-regime channel statement) — only the one near-bar grade.

The measurement behind it, `LOTTERY_tier2.md:80-90`:

| cell | floor (banked → refit) | count (banked → refit) | grade |
|---|---|---|---|
| **m2 e=0.30** | 10.05 → **9.93** (Δ −0.12) | 2.77 → **2.80** (Δ +0.03) | CONFIRMED → CONFIRMED ✓ |
| **m1 e=0.30** | 11.30 → **9.91** (Δ **−1.40**) | 0.70 → **1.20** (Δ +0.50) | **below → CONFIRMED** ✗ |

> **m1 e=0.30 ... FLIPS.** Its refit floor drops 1.40 nat across hosts and the count rises through the bar:
> **below → CONFIRMED.** The 1.40-nat move is ≈ 1 σ of the floor's own combined bootstrap error
> (0.86 ⊕ 1.02 ≈ 1.33) — i.e. **re-estimating this floor from a fresh 100-null draw, on a different
> host, moves it by about its own sampling error, and that is enough to reverse the verdict.**

and the disposition, `LOTTERY_tier2.md:148-150` (verbatim): *"**The near-bar certified-count
grade is fragile at ~1 nat** — floor sampling/host error flips below↔CONFIRMED. **RECUT's
single-e=0.3-REFUTED headline flips on a second host** and should be re-examined against
floor error (m2 e=0.30 and all structural CHORUS results are unaffected)."*

**Verdict: QUALIFY (split).**

- **QUOTABLE:** the *positive* statement — **e = 0.5 for one eccentric member, e = 0.3 for
  two or more; the threshold is a mixture property, not a property of eccentricity alone.**
  The two-member cell (`m2 e=0.30`) is well above the bar and reproduces across hosts to
  0.12 nat / 0.03 count.
- **DO NOT QUOTE as settled:** *"a single e = 0.3 member is REFUTED."* It is a ~1-nat
  coin-flip that reverses on a second host. If the paper needs it, it must be stated as
  "not confirmed" with the host-flip and the ±1 nat floor error attached — never as
  "refuted" unqualified.

---

## B4. BASELINE_field

**NOT RESOLVABLE.**

Searched exhaustively, all negative:

- No file named `BASELINE_field*` anywhere in the repo (`find -iname "*baseline*"` → no
  match outside `.git`).
- No case-insensitive occurrence of the string `baseline_field` in any tracked or
  untracked file (`grep -ril "baseline_field" .` → no match).
- No such file in git history: `git log --all --diff-filter=A -- '*BASELINE*'` → empty, so
  it was never added and later deleted.
- The token `BASELINE` occurs in exactly four places, none of them a field-search arm:
  `hpc_harbor/vernier/vernier28.py:23` and `vernier28_lib.py:11` (V27/V25's *sequential
  greedy per-source recovery* chain, labelled arm `R1 BASELINE`),
  `hpc_harbor/vernier/vernier29.py:695` (a print of "BASELINE PRODUCTS" under today's
  priors), and `CW_transition/WEAVE_etrack_merged_spec.md:38` ("THE CIRCULAR BASELINE — the
  mc cascade as established").
- None of the target figures is on disk: no `9.3`/`9.5` confident-IDs-per-realisation
  number, no `92–93 %` wrong rate, and no purity pair `0.78 / 0.96` appears in any `.md`
  in the repo in that context. (The 0.78 hits that exist are unrelated — RECUT/SURFACE
  zero-fractions and per-source counts.)

**Consequence:** the sentence *"9.3–9.5 confident IDs/realisation at 92–93 % wrong vs
criterion ~1 at purity 0.78/0.96"* has **no on-disk backing** — not the construction
(whether the field arm is an F_e-statistic and which implementation), not the two cells'
configs, and not the purity/count figures. **Do not carry it into the paper without the
run being re-staged.** Per the instruction, nothing is reconstructed from memory here.

---

## B5. SIREN c3 — the q ≥ 0.99 gate

**NOT RESOLVABLE.**

- **There is no LEDGER master on disk.** No file matches `*ledger*` (case-insensitive,
  excluding `.git`). The string `LEDGER` appears in only three documents, none of which is
  a claim-numbered master: `reports/SPARK_launch_criterion.md:113` ("THE LEDGER" — a
  five-row `s0`/`sC_*`/`sMAX` state table for SPARK's cascade arithmetic, unrelated to
  SIREN), `GLACIER_capstone.md:222` and `:562` (GLACIER bookkeeping), and
  `VERNIER2_toy.md:94`, which references the object by name but does not contain it:
  *"This is a mechanism for LEDGER-C3's banked "reach does not ..."*. That single dangling
  reference confirms a claim-numbered LEDGER **existed**; it is not in the repo.
- **The specific numbers are absent.** `3.94` occurs once in the whole repo and is
  unrelated (`lnL_distance_scans/04_gradient_distance_recovery.md:2987`,
  `3.947881e-09`). No `q ≥ 0.99` / `q = 0.99` gate table exists. No "more seeds = worse"
  statement exists.
- **What IS on disk for SIREN** (`reports/SIREN_payoff_chain.md`) is adjacent but not the
  same claim. §7 "HOW MUCH THE ANSWER LEANS ON THE CERTIFICATION PREMISE" (`:413-421`)
  prices the *distance*-knowledge premise, not a `q` gate — three arms A `exact-L` / B
  `0.1 pc` / C `½-fringe`, with arm C degrading the N = 5 answer by **5.0×**
  (0.0251 → 0.1266 dex) and the long-lag triple by **0.2 %** (`:426-428`). The nearest
  thing to "more seeds = worse" is `:309` (verbatim): *"(f = 10⁻⁷·⁷, Mc = 10⁹·⁵):
  **N_seed = 1 already saturates** (0.114 → 0.0999 from N1 to N5)"* and `:316`: *"the siren
  beyond ~6 % therefore requires *more SNR*, not more seeds"* — which is **saturation, not
  degradation**, and carries no `q` dependence at all.

**Consequence:** the inflation figures (≤ 1.04× at q = 0.99, 3.94× at q = 0.90) and the
"more seeds = worse" phrasing cannot be confirmed, and the ×5.0 / "more SNR not more
seeds" numbers that *are* banked are **different claims about a different premise** — they
must not be substituted for the missing ones. Per the instruction, nothing is reconstructed
from memory.

---

# SUMMARY OF VERDICTS

| item | verdict | one-line reason |
|---|---|---|
| **A1** data accounting | **QUALIFY** | Enumeration + construction fully sourced; `n_psr = 40` is prose-only and the M grid / seed count are off-disk (V23 has no driver; npz is 4 integers). |
| **A2** ordering anomaly | **EXCISE** the ordering, **QUALIFY** the column | "17–20" is an author hedge over a stored `20`; `20 > 17` violates `rank(Fm) ≤ rank(Fss)` and is caused by a *relative* rank tolerance rescaled per matrix, on two different draws. Quote as "both ≈ 17–20, not separated". |
| **A3** [FIG-READ] re-extraction | **QUOTABLE** (116-psr, A.2, A.3); **NOT RESOLVABLE** (15-psr) | Staged table replaces every eyeballed number; the prototype run was overwritten and A.2/A.3 live in different npz. |
| **A4** dip-and-recovery | **EXCISE** | Not in any banked array; the prototype's own record says "falls monotonically"; the only real non-monotonicity is ~1 %, seed-scatter-sized, at different N per mode. |
| **B1** criterion lineage | **QUOTABLE** | Four tags, wording and per-pulsar `np.maximum` rule all verbatim-sourced; cite §10.17 for v2.2, not its commit message. |
| **B2** EMBER predictor | **QUALIFY** | Table is exact, but the predictor is **motion** (p = 0.002), **not** engagement (p = 0.333, n.s.); 5 events total. |
| **B3** e = 0.3 refutation | **QUALIFY (split)** | The mixture statement (0.5 single / 0.3 two-member) is solid; the *single-member REFUTED* headline flips on a second host at ~1 nat. |
| **B4** BASELINE_field | **NOT RESOLVABLE** | No such file, ever (not in git history); none of the quoted figures exists on disk. |
| **B5** SIREN c3 q ≥ 0.99 | **NOT RESOLVABLE** | No LEDGER master on disk (one dangling `LEDGER-C3` reference only); `3.94` and the q-gate table absent; SIREN §7 prices a different premise. |
