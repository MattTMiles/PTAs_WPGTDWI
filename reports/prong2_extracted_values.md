# prong2 — extracted values (VERIFY-1, ACT-2 item A3)

**Read-only re-extraction.** Every number below is read directly from the banked `.npz`
files; nothing is re-simulated and nothing is eyeballed off a figure. Generated
2026-08-06 on cronus.

Sources:

| file | what it holds |
|---|---|
| `CW_transition/prong2_results.npz` | Stage A production sweep (R(N), knee, rcond, σ_L) |
| `CW_transition/prong2_coherence.npz` | Stage A.2 coherence transition (Exp 1, Exp 2) |
| `CW_transition/prong2_Ntc_map.npz` | Stage A.3 N × t_c map (R, residuals, commuting recut) |

Throughout, **R ≡ marginal / conditional distance information** (`*_marg / *_cond`).

---

## 1. Run configuration (all three prefixes in `prong2_results.npz`)

Read from `ps_*`, `tot_*`, `pop_*` scalars — **identical across all three modes**:

| quantity | value | key |
|---|---|---|
| N_psr | **116** | `ps_n_psr` |
| seeds | **30** | `ps_n_seeds` |
| T | 473 364 000 s = **15.00 yr** | `ps_T` |
| cadence | **7.0 d** (weekly) | `ps_cadence_days` |
| band | **3.0 – 20.0 nHz** | `ps_f_band` |
| σ (white noise) | **1e-7 s** = 100 ns | `ps_sigma` |
| N* = T·Δf | **8.047188** | `ps_N_star` / `ps_Tspan_bins` |
| N_max | **1000** | `ps_n_max` |
| population mode | k_loud = **3**, h_loud/h_faint = **10.0** | `pop_k_loud`, `pop_h_ratio` |

> **There is no 15-pulsar run in this file.** All three prefixes are `n_psr = 116`. The
> `diag_T10/T15/T20` labels are **Tspan in years**, not pulsar counts
> (`diag_T15_b3-20_T = 4.7336e8 s`). See §5.

---

## 2. R(N) — production 116-psr sweep, both power modes

`ps` = fixed per-source strain · `tot` = fixed total power · `pop` = 3-loud population.
Sampled N are the 16 log-grid points in `ps_n`.

| N | R (`ps`, fixed per-source) | R (`tot`, fixed total) | R (`pop`, 3 loud) |
|---:|---:|---:|---:|
| 1 | 0.99054 | 0.99054 | 0.99054 |
| 2 | 0.98603 | 0.98603 | 0.98603 |
| 3 | 0.99559 | 0.99559 | 0.99559 |
| 4 | 0.99171 | 0.99171 | 0.99488 |
| 6 | 0.98907 | 0.98907 | 0.98566 |
| 10 | 0.98613 | 0.98613 | 0.97258 |
| 16 | 0.97057 | 0.97057 | 0.98472 |
| 25 | 0.94570 | 0.94570 | 0.96003 |
| 40 | 0.95633 | 0.95633 | 0.94564 |
| 63 | 0.91733 | 0.91733 | 0.92782 |
| 100 | 0.87388 | 0.87388 | 0.86012 |
| 158 | 0.80597 | 0.80597 | 0.78090 |
| 251 | 0.70755 | 0.70755 | 0.72867 |
| 398 | 0.53097 | 0.53097 | 0.53031 |
| 631 | 0.29330 | 0.29330 | 0.30349 |
| 1000 | **0.08449** | **0.08449** | 0.07458 |

**Mode-independence** (`ps` vs `tot`): `max |ΔR| = 2.804e-12` (`mode_indep_maxdev`),
i.e. identical to machine precision at every sampled N — **not** merely "0.000".
`mode_indep_break_N = 1000` (no break found in the grid).

**Population departure**: first |ΔR| > 0.02 vs the uniform curve at `pop_depart_N = 158`.

### Derived scalars (production)

| quantity | value | key |
|---|---|---|
| 0.5-knee in N | **422.636** | `knee_N` |
| knee / N* | **52.52** | derived |
| knee, rcond ∈ {1e-10, 1e-8, 1e-12} | **422.636 / 422.636 / 422.636** | `rcond_knees` |
| rcond knee spread | **0.0** (exactly zero, all three identical) | `rcond_knee_spread` |
| production rcond | **1e-10** | `prod_rcond` |
| SNR at the knee (median, N=423) | **195.076** | `snr_knee_median`, `snr_knee_N` |
| J0437-like L | **0.15679 kpc** = 156.79 pc | `j_L_kpc` |
| EM prior line | **0.25 pc** | `em_sigma_pc` |
| N where σ_L crosses sub-pc | **NaN** (never crossed in N ≤ 1000) | `sigma_subpc_N` |
| N where σ_L crosses the EM prior | **NaN** (never crossed) | `sigma_em_N` |

---

## 3. A.2 — the coherence crossing (`prong2_coherence.npz`)

N = 1 (zero confusion), 116 psr, **24 seeds** (`n_seeds`, *not* 30 — this differs from the
Stage-A production run).

| quantity | exact banked value | key |
|---|---|---|
| **0.5-crossing, linear σ_φ² form** | **τ_p/t_c = 0.014417** | `x_half_linear` |
| **0.5-crossing, saturating form** | **τ_p/t_c = 0.004389** | `x_half_saturating` |
| shape shift (sat / lin) | **0.304443×** | `shape_shift` |
| coherent reference (t_c → ∞) | **0.993121** | `ref_ratio` |
| gate (i): t_c → ∞ reproduces it to | \|0.993121 − 0.993118\| = **3.0e-6** | `ref_ratio` vs `coh_inf` |
| reference τ_p | **1.048531e11 s** | `ref_tau` |
| mid-transition t_c used for Exp 2 | **7.272778e12 s** | `t_c_mid` |

**Exp 1 curve endpoints (grid, not extrapolated):** at τ_p/t_c = 1e-3 the linear form gives
**0.92619**; at τ_p/t_c = 1e3 it gives **1.7127e-5**.

> **Caption correction.** The prose "falls 0.99 → 0 across τ_p/t_c = 1e-3 → 1e3" conflates
> two different objects: **0.9931** is the `t_c → ∞` gate value (`ref_ratio`), *not* a grid
> point. The first actual grid point (τ_p/t_c = 1e-3) is **0.926**. Quote either
> "0.93 → 1.7e-5 across the grid" or "0.9931 in the coherent limit", not the mixture.

**Exp 2 — T-scaling (falsified pre-registered guess):**

| T (yr) | coherence R | confusion knee |
|---:|---:|---:|
| 10 | **0.615444** | 287.755 |
| 15 | **0.516349** | 412.179 |
| 20 | **0.447762** | 555.111 |
| 30 | **0.349726** | 755.983 |

Fitted form `1/(1 + T/T0)` with **T0 = 16.074 yr** (`coh_T0_yr`) — the docs' "1/(1 + T/16 yr)".

---

## 4. A.3 — factorisation residuals (`prong2_Ntc_map.npz`)

Grid: 14 × 15 (N × Y). N* = **4.260276**, band **3–12 nHz**, 116 psr, **15 seeds**
(`n_seeds` — again different from both runs above). SNR² = 85736.16, τ_med = 1.017768e11 s,
t_c_50 = 1.327637e13 s.

Residual ≡ `R − R(N,0)·R(1,Y)`.

| reduction | extremum | value | location |
|---|---|---:|---|
| **raw** (`resid`, median-of-ratios) | **peak (min)** | **−0.330595** | N/N* = **12.441** (N = 53), Y = **0.2683** |
| raw | max (positive lobe) | **+0.012031** | N/N* = 0.704 (N = 3), Y = 26.827 |
| **commuting-reduced** (`resid_commuting`) | **peak (min)** | **−0.201856** | N/N* = **27.228** (N = 116), Y = **0.1389** |
| commuting | **max (positive lobe)** | **+0.195634** | N/N* = **0.704** (N = 3), Y = **0.5179** |

Stored scalars: `max_factor_dev = −0.3305952`, `max_factor_dev_commuting = −0.2018560`.
Both match the recomputation exactly.

**Fraction of grid cells with residual < −0.1:** raw **0.371** → commuting **0.148**.

> **Two corrections that matter for the caption.**
> 1. The **raw** peak sits at N/N* = 12.4, Y = 0.27; the **commuting** peak sits at
>    N/N* = 27.2, Y = 0.139. The docs quote the commuting location (correct) but the two
>    are not the same cell — do not attach "27, 0.14" to the 0.331 number.
> 2. The commuting positive lobe is **+0.196**, i.e. **97 % of the negative peak's
>    magnitude**. Calling it a "weak positive lobe" is only true of the *raw* map (+0.012).
>    The commuting reduction did not merely shrink the negative lobe — it **manufactured a
>    positive lobe 16× larger than the raw one**. The COMPOUND verdict rests on the sign
>    structure of the interior (N/N* ≳ 3), which is unchanged; the adjective is not safe.

**1-D edge slices** (banked, used as the gates that anchor the 2-D map):

- `R_conf` (Y → 0): 0.98543 → 0.37855 across N = 1 → 256.
- `R_coh` (N = 1): 0.98543 → 0.010162 across Y = 0.01 → 100; crosses 0.5 between
  Y = 0.5179 (**0.66226**) and Y = 1.0 (**0.49574**).

---

## 5. What this file does NOT contain — and where the missing numbers went

1. **No 15-psr prototype curve.** `prong2_results.npz` was overwritten by the scaled
   116-psr production run. The only surviving record of the prototype is prose:
   `project_progress.md:141` — *"`marg/cond` falls monotonically: ~0.95 (N=1) → ~0.53
   (N=50) → ~0.18 (N=96)"* — which is itself approximate ("~") and carries no sampled-N
   grid. **Any 15-psr R(N) figure cannot be re-derived from disk.**
2. **The N_psr array-size sweep** (`stagec_p2c_results.npz`) holds only knee/N* per array
   size — `n_psr = [15, 30, 60, 116, 200]`, `knee_over_Nstar = [6.47, 12.92, 26.73,
   52.75, 91.45]` — **not** the underlying R(N) curves. So the N_psr = 15 knee is
   quotable; the N_psr = 15 *curve shape* is not.
3. **Three different seed counts** are in play and must not be cross-quoted as one run:
   Stage A production **30**, Stage A.2 coherence **24**, Stage A.3 map **15**.

---

## 6. Replacement table for eyeballed figure values

| eyeballed claim | replace with | source |
|---|---|---|
| "R ≈ 0.99 at N = 1" | **0.99054** | `ps_marg/ps_cond[0]` |
| "R ≈ 0.87 at N ≈ 100" | **0.87388** (at N = 100 exactly) | `ps_*[10]` |
| "R ≈ 0.53 at N ≈ 400" | **0.53097** (at N = 398) | `ps_*[13]` |
| "R ≈ 0.08 at N = 1000" | **0.08449** | `ps_*[15]` |
| "knee at N ≈ 410" | **422.636** | `knee_N` |
| "knee/N* ≈ 52" | **52.52** | derived from `knee_N`, `ps_N_star` |
| "rcond spread 0 %" | **exactly 0.0** — all three rcond give 422.636 | `rcond_knees` |
| "modes agree to 0.000" | **2.804e-12** | `mode_indep_maxdev` |
| "0.5-crossing 0.0144 / 0.0044, 0.30× shift" | **0.014417 / 0.004389 / 0.304443** | `prong2_coherence.npz` |
| "peak residual 0.331 → 0.202" | **−0.330595 → −0.201856** (signed; different cells) | `prong2_Ntc_map.npz` |
| "weak positive lobe" | **+0.196 commuting** (vs +0.012 raw) — not weak | `resid_commuting` |
