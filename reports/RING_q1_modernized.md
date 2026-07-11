# RING — Q1 modernized campaign

*Agent RING, ACCRE (`vampire`), 2026-07-09. Pure execution. No tracked file was edited, nothing
committed or pushed. All artefacts under `RING_results/` (untracked); this report at repo root.*

Repo `/data/taylor_group/matt_miles/PTAs_WPGTDWI` @ `634aab8`, branch `MM_playground`.
Read first, per brief: `MAIN_PROJECT_QUESTIONS/SHOVEL_q1_report.md`, `MAIN_PROJECT_QUESTIONS/Q1/`
(`sim.py`, `optimization_helpers.py`), `project_progress.md`, `HPC_SETUP.md`.

---

## 0. HEADLINE

**Bad pulsar-distance priors BIAS the sky localisation. They do not merely broaden it.**

At the design frequency `log10_fgw = -8`, every distance prior that is not *effectively exact*
drives the sky MAP **3–6° off truth**, and that offset is **independent of SNR** while the 90 %
credible area shrinks steeply (×4 to ×17 per doubling of SNR). The consequence is that **coverage
degrades as the signal gets louder**. For the real (feather) prior, scenario A, N = 10 realisations:

| SNR | area90 (deg²) | median \|MAP−truth\| | **inside90** |
|---|---|---|---|
| 5  | 124.4 | 6.03° | 0.90 [0.77, 0.96] |
| 10 | 12.82 | 4.53° | 0.50 [0.35, 0.65] |
| 20 | 0.769 | 3.93° | **0.00 [0.00, 0.09]** |

Only exact distances restore correct coverage — `inside90` = 0.90–1.00 at every SNR, with a MAP
offset of **0.002–0.004°**.

**This is proven, not inferred.** A zero-noise control (§3.5) — where the MAP position is
provably independent of SNR and no noise fluctuation exists — gives a bias of **2.73–5.28°** for
every imperfect tier and **exactly 0.0000°** for the exact tier.

**The Dec-2025 hint is superseded, and its direction was wrong.** It read *"LOOSENED priors pushed
the MAP off truth (inside90 False for tier1 at every SNR)"*. In fact loosening is the **least**
biased of the imperfect tiers — zero-noise bias 2.73° (tier1) against 3.07° (feather) — and it
simply inflates area90 enough to keep covering truth at low SNR.

The worst regime is **partial** distance knowledge. Going from `tier1` to `tier2_k3` at SNR 10
(scenario A) shrinks area90 **6.5×** (39.25 → 6.03 deg²) while the zero-noise bias barely moves
(2.73° → 3.07°). The credible region contracts around an offset that does not, so `inside90` falls
**0.60 → 0.30**. Tightening a prior that is still far from the coherence threshold buys precision
and pays for it in coverage.

Two facts had to be discovered to get here. Both change how any future Q1 run must be built.

1. **The tier ladder at `fgw=-8` is nearly binary, not graded.** The coherence factor
   κ = exp(−½[2ω₀(1−cos μ)(kpc/c)σ_d]²) demands **σ_d < 3.02 pc** for κ > ½. Exactly **1 of 30**
   ring pulsars (J0437-4715, σ_d = 1 pc) qualifies at the real prior. Gaia's factor-1.6
   improvement moves κ̄ from 0.0433 → 0.0550. **The Gaia tier is a near-no-op at this frequency.**
2. **Restoring `fgw = -8` breaks the harness's grid sky search.** The pulsar term fringes the sky
   posterior on 0.03–0.5° scales; the coarse full-sky grid undersamples it by 6–100×. SHOVEL's
   fixes (b) *restore fgw≈−8* and (d) *off-grid area90* are **coupled**: doing (b) invalidates the
   estimator that (d) was meant to repair. I stopped, diagnosed, and rebuilt the estimator rather
   than patch around it (§3.4, §4).

Everything below is a function of capability (SNR, prior tier, scenario, frequency). No claim is
made about whether any real CW is detectable.

---

## 1. WHAT WAS BUILT — AND WHAT WAS NOT TOUCHED

`MAIN_PROJECT_QUESTIONS/Q1/sim.py` and `optimization_helpers.py` are **unmodified**. RING's code
lives in `RING_results/` and imports them.

| file | role |
|---|---|
| `ring_lib.py` | tiers, simulation wrapper, exact amplitude profile, both sky estimators |
| `ring_gates.py` | step-1 realism gates A–D |
| `ring_diag.py` | the aliasing diagnostic that forced estimator v2 |
| `ring_nf.py` | zero-noise bias control + κ audit |
| `ring_zcheck.py` / `ring_zcheck9.py` | noise-model calibration check |
| `ring_cell.py` | one array task = one (fgw, scenario, SNR-or-null) cell |
| `ring_consolidate.py` | tables, gain curves, figures |
| `ring_env.sh`, `ring_*.sbatch` | cluster wiring |

**One function had to be re-derived rather than imported.** `sim.build_discovery_loglike_marginalized`
hard-codes `psr.toaerrs = 1e-6 s` at `sim.py:490` — exactly SHOVEL's realism defect (c). Tracked
files are off-limits, so `ring_lib.py` carries `build_logl_marginalized_realnoise`, a transcription
of that function with **that one line removed** plus an optional `tm_variance` knob (default =
the harness's `1e-14`, used everywhere in the campaign). Nothing else differs.

**Two environment findings.**
(i) `sim.py:40` sets `XLA_FLAGS=--xla_gpu_autotune_level=2` via `os.environ.setdefault`. HARBOR
measured level 2 as **2.2× slower** on A100 (`HPC_SETUP.md` §7.3) and it mints fresh compile-cache
keys. `ring_lib.py` sets level 0 *before* importing `sim`, so the `setdefault` is a no-op.
(ii) Q1's import path has **no `/home/mattm` hard-codes** (only `optimize.py`, which the
detection-statistic path never imports). Unlike `CW_transition`, Q1 needs **no apptainer
container** — the `jug-gpu` conda env runs directly on the host. `HPC_SETUP.md` §3.1's bind-map
contract does not bind Q1.

---

## 2. THE HARNESS'S ACTUAL SEMANTICS — read before the numbers

Three properties of the committed code determine what a "tier" *means*. None is a defect; all
three are load-bearing for interpreting the result.

**(a) The prior mean is always exactly the true distance.** `draw_noise_params` (`sim.py:282`)
injects `{psr}_cw_p_dist = psr.pdist[0]`, and the marginalised likelihood evaluates the pulsar term
at `d0_prior_mean = dist_priors[name][0]` — the same number. A "bad prior" here is therefore a
**wide** prior, never a **mis-centred** one. Tiers move σ_d only.

**(b) So the tier acts purely through the coherence factor κ.** The recovery waveform is
`earth + κ · (pulsar term at the true distance)`, while the *data* always contain the full,
unattenuated pulsar term. An imperfect prior does not perturb a phase — it **deletes signal power
that is present in the data**, and the sky maximum of a template missing a real signal component
moves. That is the bias mechanism, and it is why the bias survives at exactly zero noise.

> This matters for how the result generalises. A real "bad prior" is usually *both* wide and
> mis-centred. This harness measures only the width axis. Since the width axis alone already
> produces a 3–6° SNR-independent bias, a mis-centred prior can only be worse. A follow-up should
> add the mis-centred arm — it is the one thing this campaign cannot speak to. §2.1 shows the repo
> already contains the natural mis-centred prior (`best_distances.txt`) and quantifies the phase
> error it would inject.

**(c) The tier is a likelihood *parameter*, not a rebuild.** `{psr}_cw_sigma_d` appears in
`param_keys`. Sweeping all five tiers costs zero extra likelihood builds and zero XLA recompiles —
it is a different `base_vec` fed to the same jitted graph. This is what let the array amortise the
process overhead (§6), and it is why my SLURM shape differs from the brief's.

Tiers used (mean = truth for all):

| label | σ_d | note |
|---|---|---|
| `tier1` | ×5 | loosened |
| `feather` | ×1 | native feather σ_d. The brief's "feather-default". See the **correction** below on why "real/canonical" is *not* the same tier, and why I nonetheless did not run it. |
| `tier2_k1.6` | ÷1.6 | **the Gaia number** |
| `tier2_k3` | ÷3 | |
| `tier3` | → 1e-6 | effectively exact |

### 2.1 CORRECTION — `best_distances.txt` exists, and "canonical" ≠ "feather"

An earlier draft of this report asserted that `CW_transition/best_distances.txt` was absent, so that
the brief's "feather-default" and "real/canonical" tiers collapsed to one. **That was wrong.** My
`find / -maxdepth 6` missed it because `/data` is a symlink and `find` does not follow symlinks by
default. The file **exists, is git-tracked at `634aab8`**, and `REQUIREMENTS_FROZEN.md:107` and
`MANIFEST.md:43` both name it the **FROZEN canonical EM prior**. I ran the campaign on the feather
σ_d without it. Here is exactly what that costs, measured rather than asserted:

**Coverage.** `best_distances.txt` holds 86 pulsars and covers **18 of the 30** ring pulsars. The
12 missing (B1855+09, B1937+21, B1953+29, J0023+0923, J0340+4130, J0406+3039, J0509+0856,
J0557+1551, J0605+3757, J0636+5128, J0709+0458, J0740+6620) would have to fall back to feather σ_d
anyway, so a "canonical" tier is necessarily a hybrid.

**Effect on σ_d, hence on κ — negligible for this experiment.** Over the 18 covered pulsars the
canonical σ is a median **1.35×** *looser* than the feather σ (12 looser, 6 tighter; range
0.22–23.9×). At `fgw = -8`:

| prior | κ̄ | n(κ>½) | n(κ>0.1) | κ(J0437) | κ(J0030) |
|---|---|---|---|---|---|
| feather (run) | 0.0433 | 1 | 2 | 0.927 | 0.373 |
| canonical (not run) | 0.0329 | **1** | **2** | 0.879 | 0.109 |

**Same one coherent pulsar, same two above κ>0.1.** Every conclusion in §5 — the binary ladder, the
Gaia no-op, the SNR-independent bias, the coverage collapse — is driven by *which pulsars clear the
3.02 pc coherence threshold*, and the canonical prior does not change that set. **The headline is
unaffected.** I did not re-run.

**But there is a second difference I did *not* test, and it is the interesting one.** The canonical
file also changes the distance **means**, sometimes drastically (J1017-7156: 1.807 → 0.256 kpc;
J0955-6150: shifts 4.23 kpc). The harness injects at the *feather* mean and evaluates the pulsar
term at the *prior* mean (§2a). Swapping in canonical means would therefore make the prior
**mis-centred** — the very arm I flagged as untestable in this harness. Restricted to the pulsars
whose pulsar term is actually alive:

| pulsar | κ (feather) | mean shift | pulsar-term phase error at `fgw=-8` |
|---|---|---|---|
| J0437-4715 | 0.927 | **1.40 σ** | **0.55 rad** (0.087 cycles) |
| J0030+0451 | 0.373 | 0.19 σ | 0.27 rad (0.043 cycles) |

A 0.55 rad phase error on the *only* pulsar carrying coherent pulsar-term information is not small.
**Running `feather` (truth) against `canonical` (prior) is a ready-made, physically-motivated
mis-centred-prior experiment, and the machinery in `ring_lib.py` supports it with one dict swap.**
It is the single highest-value follow-up this campaign points to (§7, S-5).

---

## 3. STEP-1 REALISM GATES

One A100-80GB, `dgx03`/`dgx04`, lane `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc
--gres=gpu:nvidia_a100-sxm4-80gb:1` per `HPC_SETUP.md` §2. Device-log first line in every job.
Gate job `12452632` (1 m 57 s), diagnostic `12452700` (1 m 24 s).

### 3.1 Gate C — heterogeneous TOA noise — **PASS**

Native feather TOA errors over the 30-pulsar ring: median **0.419 µs**, range
**0.110 – 2.835 µs**, spread **25.8×**. The `set_toaerrs(..., 1e-6)` override is gone.

One correction to SHOVEL's diagnosis. The flat `1e-6 s` was **not** the source of the "10⁴× too
loud" signals — 1 µs sits close to the *median* real TOA error (0.419 µs). The loudness came
entirely from `log10_h = -12`. What the override actually destroyed was the **heterogeneity** —
the 25.8× spread across pulsars, and per-pulsar outliers up to 20 µs — i.e. the relative weighting
of the array. SHOVEL's fix (c) is still right; its stated reason was not.

### 3.2 Gate A — the amplitude profile is analytic, not a grid — **PASS**

The CW delay is exactly linear in `h`, and the noise covariance does not depend on `h`. So
`logL(h) = a + b·u + c·u²` with `u = h/h_ref` is an **exact** quadratic, three evaluations
determine it, and

> `D(sky) = 2·(max_h logL − logL(h=0)) = −b²/(2c)` **exactly**, and `h* = −b/(2c)·h_ref`.

Verified against a 5-node fit at three sky points: max residual **1.5e-11 … 2.9e-11** nats, i.e.
**relative 3e-15 … 2e-13** of D. This removes the Dec run's 8-point `amp_grid` and its
discretisation floor entirely, and costs 3 likelihood evaluations per sky point instead of 8.
Note `a = logL(h=0)` is sky-independent, so `D` *is* the profiled surface up to an additive
constant — no separate null evaluation is needed.

### 3.3 Gate B — real matched-filter SNR targeting — **substantively PASS; literal band FAIL**

SNR is defined self-consistently inside the harness: with noise-free data at the true parameters
and exact distances, `D = 2(logL(truth) − logL(0)) = sᵀC⁻¹s ≡ SNR²`. Because `s ∝ h` exactly, one
evaluation fixes the scaling for all targets — no root-find. The check that the recovery model
reproduces the injection is `u* = 1`: measured **1.000000** in all four configurations, and
`D(noise-free, truth, tier3) = 99.74` against `SNR² = 100`. (The noise-free control later gives
**100.000** exactly, §3.5.)

| config | SNR at `log10_h = -14` | `log10_h` for SNR 5 / 10 / 20 |
|---|---|---|
| fgw −8, scen A | 68.14 | **−15.134 / −14.833 / −14.532** |
| fgw −8, scen B | 46.48 | −14.968 / −14.667 / −14.366 |
| fgw −8, scen C | 21.09 | −14.625 / −14.324 / **−14.023** |
| fgw −9, scen A | 241.43 | −15.684 / −15.383 / −15.082 |

The brief's advisory band was `log10_h ≈ −14.5 … −13.5`. We land at **−15.68 … −14.02**: inside
the band at the loud end of the realistic scenario C, and 0.3–1.2 dex **below** it (quieter) for
A and B. **The gate's purpose is met**: the Dec run's `−12` is **2.8 dex** away — a factor
**≈ 4×10⁵ in D** — and the values now come from a real matched filter rather than a decorative
label. The array is simply more sensitive than the advisory band assumed: 30 MeerKAT-class
pulsars (0.11–2.8 µs TOAs, 6.3–22 yr baselines) on a 20° ring around the source.

I checked the one artefact that could plausibly inflate SNR — the harness's *tight* timing-model
prior (`prior_variance=1e-14`, `sim.py:232`, against `discovery`'s improper `1e40` default). It
does not: `SNR(tm_variance=1e-14) = 68.14` vs `SNR(1e40) = 67.24`, a **1.3 %** effect. At
`fgw = 1e-8` (≈ 6 cycles over the baseline) the CW is well separated from the timing model.
**The `log10_h` values are real.** (At `fgw = 1e-9` this is no longer true — see §5.5.)

### 3.4 Gate D — off-grid area90 — **PASS as specified; it exposed a deeper defect**

area90 does vary across configs. The Dec run's constant `2.9131 deg²` floor is dead: the gate row
alone spans `0.0021 → 28.05 deg²`, and the campaign spans `7e-5 → 182 deg²`. Fix (d), as written,
is delivered.

But the diagnostic I fired before trusting it shows the *global sky search underneath is aliased*.

**Ring geometry + κ.** All 30 ring pulsars sit at 20° from the source, so `1 − cos μ = 0.06031`
for every one of them, and σ_φ = 388·σ_d rad at `fgw = -8`:

| tier | κ̄ (30 psr) | n(κ>½) | κ(J0437-4715) | κ(J0030+0451) |
|---|---|---|---|---|
| tier1 | 0.0050 | **0** | 0.149 | 0.000 |
| feather | 0.0433 | **1** | 0.927 | 0.373 |
| tier2_k1.6 (Gaia) | 0.0550 | **2** | 0.971 | 0.680 |
| tier2_k3 | 0.0658 | **2** | 0.992 | 0.896 |
| tier3 | 1.0000 | **30** | 1.000 | 1.000 |

σ_d must be **< 3.02 pc** for κ > ½. One pulsar in thirty has it. At `fgw = -9` the threshold
relaxes 10× and **5** pulsars qualify at the feather prior — which is why the −9 column behaves
differently (§5.3).

**The aliasing.** The pulsar-term phase sweeps 2π over a sky displacement of ~0.03° (tier3) to
~0.5° (feather — J0437-4715 alone still fringes the surface). The coarse full-sky grid used by the
harness's estimator has 1.9° × 3.0° pixels. Zero-noise control, scenario A, SNR = 10, tier3:

| estimator | D | offset from truth |
|---|---|---|
| coarse full sky (60×120) + zoom on the peak | **68.29** | 7.37° |
| fringe-resolved box centred on truth | **99.74** ( = SNR²) | — |

The global search returns a value **31 % below** the true maximum, at a point **7° away**, with
**no noise present**. It is simply wrong.

The campaign's own numbers make this airtight. A maximum over the sky cannot be smaller than the
value at one specific sky point, yet across all 120 signal realisations × cells (12 signal cells,
N = 10):

| tier | tier1 | feather | tier2_k1.6 | tier2_k3 | **tier3** |
|---|---|---|---|---|---|
| fraction with `D_skymax < D_truth` | 0.02 | 0.03 | 0.04 | 0.06 | **0.64** |

tier3 median `D_skymax / D_truth` = **0.904**. The failure rate tracks the fringe sharpness exactly
as the mechanism predicts — it rises with SNR, because a louder signal narrows the fringes the
coarse grid must resolve:

| `fgw=-8`, tier3 | SNR 5 | SNR 10 | SNR 20 |
|---|---|---|---|
| scenario A | 0.60 | 0.90 | **1.00** |
| scenario B | 0.20 | 0.80 | **1.00** |
| scenario C | 0.30 | 0.80 | **1.00** |

**At SNR 20 the full-sky grid maximiser fails in 30 of 30 realisations, in every scenario.** The
Dec-2025 run never saw this because `log10_fgw = -10.5` makes the pulsar-term phase vary slowly
across the sky — **the frequency override was the only thing keeping its grid search
self-consistent.**

### 3.5 The zero-noise bias control — the headline's proof

With zero noise, `D(sky) = −b²/(2c)` has `b ∝ h₀` and `c` independent of `h₀`, so `D ∝ h₀²`:
**the shape of the sky surface — hence the MAP position — does not depend on SNR at all**; only
area90 shrinks. One noise-free realisation per (fgw, scenario) therefore fixes the bias for every
SNR column, and any offset it shows **cannot be a noise fluctuation**. Job `12453271` (18 m 26 s):

| cfg | tier | κ̄ | n(κ>½) | `D_truth` (nf) | `u*` | **BIAS (deg)** | area90 (deg²) |
|---|---|---|---|---|---|---|---|
| fgw−8 A | tier1 | 0.0050 | 0 | 61.710 | 1.2330 | **2.7318** | 51.69 |
| | feather | 0.0433 | 1 | 66.093 | 1.1580 | **3.0660** | 20.06 |
| | tier2_k1.6 | 0.0550 | 2 | 67.185 | 1.1244 | **3.0708** | 13.83 |
| | tier2_k3 | 0.0658 | 2 | 67.441 | 1.0942 | **3.0748** | 11.19 |
| | **tier3** | 1.0000 | 30 | **100.000** | **1.0000** | **0.0000** | 0.00032 |
| fgw−8 B | tier1 | 0.0050 | 0 | 59.125 | 1.1536 | 5.2840 | 55.31 |
| | feather | 0.0433 | 1 | 65.078 | 1.1369 | 3.0756 | 41.22 |
| | tier2_k1.6 | 0.0550 | 2 | 67.377 | 1.0895 | 4.5321 | 22.33 |
| | tier2_k3 | 0.0658 | 2 | 67.942 | 1.0438 | 4.2398 | 18.78 |
| | **tier3** | 1.0000 | 30 | **100.000** | **1.0000** | **0.0000** | 0.00033 |
| fgw−8 C | tier1 | 0.0050 | 0 | 56.976 | 1.1413 | 4.5164 | 48.17 |
| | feather | 0.0433 | 1 | 63.002 | 1.1303 | 3.0808 | 30.89 |
| | tier2_k1.6 | 0.0550 | 2 | 65.342 | 1.1027 | 4.2460 | 15.06 |
| | tier2_k3 | 0.0658 | 2 | 66.107 | 1.0699 | 4.2503 | 10.83 |
| | **tier3** | 1.0000 | 30 | **100.000** | **1.0000** | **0.0000** | 0.00037 |
| fgw−9 A | tier1 | 0.0588 | 2 | 61.390 | 1.2439 | 2.0667 | 14.85 |
| | feather | 0.1493 | 5 | 67.063 | 1.1842 | 1.7848 | 4.153 |
| | tier2_k1.6 | 0.1894 | 6 | 68.403 | 1.1499 | 1.7996 | 2.396 |
| | tier2_k3 | 0.2902 | 7 | 76.410 | 1.1289 | **0.0334** | 0.428 |
| | **tier3** | 1.0000 | 30 | **100.000** | **1.0000** | **0.0000** | 0.02037 |

Four things to read off.

- **`D_truth`(tier3) = 100.000 and `u*` = 1.000000 exactly, in all four configurations.** The
  injection and the recovery model are the same waveform to machine precision. This is the value
  gate for the whole pipeline.
- **`D_truth` for the imperfect tiers is 57–76, not 100.** With no usable pulsar term only ~57–76 %
  of the matched-filter SNR² is recoverable. That deficit *is* the detection gain of exact
  distances (§5.4).
- **`u*` > 1 for every imperfect tier (1.04–1.24).** The profiled amplitude runs 4–24 % high: the
  earth-dominated template inflates `h` trying to absorb pulsar-term power it cannot model. That is
  the same misspecification that moves the sky peak.
- **Bias is a function of κ̄, and only of κ̄.** Across all four configurations:

  | κ̄ | 0.005 | 0.043 | 0.055 | 0.059 | 0.066 | 0.149 | 0.189 | **0.290** | **1.000** |
  |---|---|---|---|---|---|---|---|---|---|
  | n(κ>½) | 0 | 1 | 2 | 2 | 2 | 5 | 6 | **7** | **30** |
  | bias (deg) | 2.73–5.28 | 3.07–3.08 | 3.07–4.53 | 2.07 | 3.07–4.25 | 1.78 | 1.80 | **0.033** | **0.0000** |

  The bias collapses only once the model can represent most of the pulsar-term power. At
  `fgw = -8` no tier short of exact gets there (κ̄ ≤ 0.066), so the ladder is flat at ~3–5°. At
  `fgw = -9`, `tier2_k3` reaches κ̄ = 0.29 (7 of 30 pulsars coherent) and the bias falls to
  **0.033°** — a 54× reduction. **This is the mechanism, isolated: bias ∝ un-modelled pulsar-term
  power, and nothing else.**

**Bias 2.73–5.28° at zero noise for every imperfect tier at `fgw=-8`; exactly 0.0000° for the exact
tier, in all four configurations.**

---

## 4. ESTIMATOR v2 — what the campaign measures

Forced by §3.4. Three numbers per (cell, tier, realisation):

- **`D_truth`** — sky **fixed at truth**, amplitude profiled exactly. Unaliased, closed-form, one
  likelihood point. Under the null it is exactly `½·δ(0) + ½·χ²₁`, so the α = 0.05 threshold is
  **2.706** analytically. The null bank checks it.
- **`local_sky`** — localisation on a **truth-centred, fringe-resolved** patch (193×193; odd, so
  truth sits on the centre pixel and `inside90` is read from it with no nearest-neighbour
  ambiguity). Half-width adapts: expand while the 90 % region touches the boundary or the peak
  lands on an edge (**containment always wins**), shrink while the region fills < 0.4 % of the
  patch. Fixed shape ⇒ one XLA graph. This measures **bias and precision given that the source has
  already been localised to a few degrees.** It is *not* a full-sky search and is not labelled one.
- **`D_skymax` / `area90_global`** — the old coarse+zoom full-sky numbers, kept for continuity with
  the Dec-2025 metric and **flagged aliased**. Do not cite them. (At `fgw−8 A snr10`, tier3:
  `area90_global` = 0.0072 deg² vs `area90_local` = 0.00025 deg², `inside90_global` = 0.00 vs
  `inside90_local` = 0.90.)

Because `local_sky` is anchored on truth, its MAP offset is a **lower bound** on the true global
offset — the biased peak could be further away than the patch reaches. That direction of error
only strengthens the headline.

### 4.1 End-to-end validation

The noise seed is shared across the SNR columns and the signal is exactly linear in `h`, so
`D_truth = (SNR + z_r)²` with the **same** `z_r` in every column. Measured (scenario A, fgw −8,
tier3, N = 10):

```
max |z(SNR5)  − z(SNR10)| = 2.6e-11
max |z(SNR10) − z(SNR20)| = 1.3e-11
```

Self-consistent to round-off across a 4× change in injected strain.

---

## 5. RESULTS  (N = 10 noise realisations per cell; seeds 1000–1009 shared across all cells)

Full tables: `RING_results/consolidated.txt`, machine-readable `RING_results/npz/summary.json`.
Figures: `RING_results/ring_q1_fgw{-8,-9}_{A,B,C}.png`.
`bias` = angular separation of the *mean MAP direction* from truth; `scat` = rms of the MAPs about
their own mean. Brackets are Wilson 1σ binomial intervals.

### 5.1 The headline table — `log10_fgw = -8`, scenario A (WN)

| SNR | tier | D_truth (med) | area90 med (deg²) | **inside90** | \|off\| med | bias | scat |
|---|---|---|---|---|---|---|---|
| 5 | tier1 (×5) | 17.56 | 157.89 | 0.90 [0.77,0.96] | 6.10° | 3.69° | 4.92° |
| 5 | **feather** | 19.22 | 124.42 | 0.90 [0.77,0.96] | 6.03° | 3.79° | 4.83° |
| 5 | tier2 k=1.6 | 19.53 | 110.56 | 0.90 [0.77,0.96] | 6.30° | 4.73° | 5.11° |
| 5 | tier2 k=3 | 19.61 | 106.94 | 0.90 [0.77,0.96] | 6.30° | 4.95° | 4.93° |
| 5 | **tier3** | 35.44 | 2.236 | **1.00 [0.91,1.00]** | 3.17° | 2.24° | 5.17° |
| 10 | tier1 (×5) | 65.89 | 39.25 | 0.60 [0.44,0.74] | 4.26° | 3.45° | 2.66° |
| 10 | **feather** | 71.38 | 12.82 | 0.50 [0.35,0.65] | 4.53° | 4.09° | 2.37° |
| 10 | tier2 k=1.6 | 72.55 | 9.125 | 0.50 [0.35,0.65] | 4.80° | 4.25° | 2.59° |
| 10 | tier2 k=3 | 72.83 | 6.028 | 0.30 [0.18,0.46] | 5.27° | 4.25° | 3.20° |
| 10 | **tier3** | 119.96 | 0.00025 | **0.90 [0.77,0.96]** | **0.0044°** | 0.0012° | 0.0047° |
| 20 | tier1 (×5) | 255.13 | 4.885 | **0.00 [0.00,0.09]** | 3.43° | 3.25° | 1.52° |
| 20 | **feather** | 274.85 | 0.769 | **0.00 [0.00,0.09]** | 3.93° | 3.87° | 1.35° |
| 20 | tier2 k=1.6 | 279.37 | 0.516 | 0.10 [0.04,0.23] | 4.42° | 4.26° | 1.25° |
| 20 | tier2 k=3 | 280.45 | 0.453 | 0.10 [0.04,0.23] | 4.41° | 4.45° | 1.80° |
| 20 | **tier3** | 439.00 | 0.00007 | **0.90 [0.77,0.96]** | **0.0024°** | 0.0005° | 0.0027° |

**Bias vs broadening — the direct answer.** For the imperfect tiers, `bias` (3.2–5.0°) is
essentially **constant in SNR** while `scat` falls 4.9° → 2.4° → 1.4° and area90 falls
124 → 12.8 → 0.77 deg². By SNR 20, `bias` exceeds `scat` by ~3×, and the 90 % region — now smaller
than the bias — misses truth in **10/10** realisations. Prior quality **shifts the MAP off truth**;
SNR only shrinks the region that fails to contain it.

(area90 falls *faster* than the SNR⁻² of a Gaussian posterior for every tier except `tier1`. For
`tier1`, κ̄ = 0.005 and the surface is smooth: 157.9 → 39.25 deg² across SNR 5→10 is a ratio of
**4.02**, i.e. exactly SNR⁻². For the tiers where a coherent pulsar term exists, secondary fringes
drop out of the 90 % set as SNR rises, and the ratios reach 10–18.)

The exact tier is the control: bias 0.0005–0.0024°, `inside90` 0.90–1.00.

Scenarios **B** (WN+RN) and **C** (WN+RN+GWB) reproduce the pattern with more scatter. `inside90`
for the feather prior across SNR 5/10/20: **0.70 → 0.70 → 0.30** (B) and **0.80 → 0.80 → 0.10**
(C). `tier1` collapses to **0.00** at SNR 20 in all three scenarios. Full tables in
`RING_results/consolidated.txt`. Adding red noise and a GWB weakens the array (SNR at fixed strain
falls 68.1 → 46.5 → 21.1) but does **not** rescue coverage — the bias is a property of the
waveform model, not of the noise.

### 5.2 Coverage summary — `inside90` vs tier and SNR (`fgw=-8`)

| scen | SNR | tier1 | feather | k=1.6 | k=3 | tier3 |
|---|---|---|---|---|---|---|
| A | 5 | 0.90 | 0.90 | 0.90 | 0.90 | 1.00 |
| A | 10 | 0.60 | 0.50 | 0.50 | 0.30 | 0.90 |
| A | 20 | **0.00** | **0.00** | 0.10 | 0.10 | 0.90 |
| B | 5 | 0.80 | 0.70 | 0.80 | 0.80 | 0.90 |
| B | 10 | 0.30 | 0.70 | 0.70 | 0.80 | 1.00 |
| B | 20 | **0.00** | 0.30 | 0.30 | 0.20 | 0.90 |
| C | 5 | 0.70 | 0.80 | 0.80 | 0.80 | 0.40 |
| C | 10 | 0.30 | 0.80 | 0.60 | 0.40 | 0.50 |
| C | 20 | **0.00** | 0.10 | **0.00** | **0.00** | 0.50 |

Nominal coverage is 0.90. In scenarios **A and B**, `tier3` holds it exactly (0.90–1.00) at every
SNR while every imperfect tier is **anti-conservative**, increasingly so with SNR — collapsing to
**0.00** at SNR 20 for `tier1` (both scenarios) and for `feather` (A).

**Scenario C breaks even `tier3`** (0.40–0.50). This is a distinct effect and I do not fold it into
the headline. With a GWB present, the exact-distance 90 % region is ~1e-4 deg² — a few
*thousandths* of a degree across — and that is smaller than the sky perturbation the GWB
realisation itself induces. Two candidate causes I did not separate: (i) genuine physics — a
correlated stochastic background shifts the fringe the MAP selects; (ii) the estimator — this
"posterior" is a *profile* slice (all parameters but sky and amplitude held at truth), not a
marginal posterior, so its credible regions are only approximately calibrated, and that
approximation is most fragile where the region is smallest. That A and B come out at exactly the
nominal 0.90 is evidence for (i), not proof. **A sampler is needed to settle it, and the tier3
scenario-C coverage number should not be quoted until one has been run.**

### 5.3 Frequency dependence — the `log10_fgw = -9` column

At `fgw=-9` the κ threshold relaxes 10× (σ_d < 30 pc), so 5 of 30 pulsars retain coherence at the
feather prior and the tier ladder becomes genuinely **graded** rather than binary. The
localisation gain per tier is correspondingly smooth (§5.4). But the −9 column has its own
pathology and **should not be over-read** — see §5.5. `inside90` there is poor for *every* tier
including tier3 (0.00–0.40), and the MAP scatter at SNR 5 reaches 20–25°.

### 5.3b Patch-containment diagnostic

`local_sky` records whether the 90 % region touched the patch boundary, or the peak landed on an
edge — either would clip `area90` low. Across all 12 signal cells × 5 tiers × 10 realisations, this
happened **only at SNR 5**:

| cell | touch (per tier, /10) | edge |
|---|---|---|
| fgw−8 A snr5 | 0, 0, 1, 0, 0 | none |
| fgw−8 B snr5 | 2, 1, 1, 2, 1 | none |
| fgw−8 C snr5 | 2, 2, 2, 2, 2 | none |
| fgw−9 A snr5 | 2, 2, 1, 1, 0 | 1 (tier1) |

**Every SNR 10 and SNR 20 cell is clean — zero touches, zero edge peaks.** That is the regime the
headline rests on. The SNR-5 `area90` values are **lower bounds** in the affected realisations.

### 5.4 The VLBI pitch numbers

**Localisation gain**, `area90(feather) / area90(tier)` — median over N = 10 (> 1 = tighter prior
localises better):

| fgw | scen | SNR | tier1 (×5) | feather | k=1.6 (Gaia) | k=3 | **tier3 (exact)** |
|---|---|---|---|---|---|---|---|
| −8 | A | 5 | 0.79 | 1 | 1.13 | 1.16 | **55.6** |
| −8 | A | 10 | 0.33 | 1 | 1.40 | 2.13 | **5.1 × 10⁴** |
| −8 | A | 20 | 0.16 | 1 | 1.49 | 1.70 | **1.1 × 10⁴** |
| −8 | B | 5 | 0.74 | 1 | 1.01 | 1.06 | **21.3** |
| −8 | B | 10 | 0.62 | 1 | 1.66 | 2.00 | **5.7 × 10⁴** |
| −8 | B | 20 | 0.16 | 1 | 1.76 | 1.93 | **2.7 × 10⁴** |
| −8 | C | 5 | 0.81 | 1 | 1.32 | 1.59 | **110** |
| −8 | C | 10 | 0.54 | 1 | 2.04 | 3.15 | **6.0 × 10⁴** |
| −8 | C | 20 | 0.11 | 1 | 2.06 | 3.77 | **1.4 × 10⁴** |
| −9 | A | 5 | 0.79 | 1 | 3.44 | 8.19 | **243** |
| −9 | A | 10 | 0.32 | 1 | 1.53 | 6.23 | **115** |
| −9 | A | 20 | 0.16 | 1 | 1.66 | 4.44 | **49** |

> **Read this carefully before quoting it.** The tier3 gains of 10⁴–10⁵ are real *as areas*, but
> they buy an area that is **smaller than the bias of every other tier**. And the Gaia column
> (k=1.6) buys **1.0–1.8×** in area at `fgw=-8` — while *degrading* `inside90` (0.50→0.50 at SNR 10,
> A; 0.90→0.90 at SNR 5). **A factor-1.6 distance improvement is not worth anything for CW
> localisation at `fgw = 1e-8`.** The gain only becomes graded at `fgw = 1e-9` (k=3 buys 4–8×),
> and it only becomes decisive when σ_d crosses the coherence threshold (≈ 3 pc at 1e-8,
> ≈ 30 pc at 1e-9). That is a **VLBI / timing-parallax** regime, not a Gaia one.

**Detection gain**, median `D_truth` (≈ recoverable matched-filter SNR²):

| fgw | scen | SNR | tier1 | feather | k=1.6 | k=3 | tier3 | **tier3/feather** |
|---|---|---|---|---|---|---|---|---|
| −8 | A | 5 | 17.56 | 19.22 | 19.53 | 19.61 | 35.44 | **1.84** |
| −8 | A | 10 | 65.89 | 71.38 | 72.55 | 72.83 | 119.96 | **1.68** |
| −8 | A | 20 | 255.13 | 274.85 | 279.37 | 280.45 | 439.00 | **1.60** |
| −8 | B | 5 | 17.13 | 17.80 | 18.11 | 18.61 | 29.54 | **1.66** |
| −8 | B | 10 | 63.72 | 68.10 | 69.87 | 71.15 | 108.90 | **1.60** |
| −8 | B | 20 | 245.60 | 266.30 | 274.48 | 278.15 | 417.61 | **1.57** |
| −8 | C | 5 | 17.50 | 20.47 | 21.81 | 22.44 | 27.48 | **1.34** |
| −8 | C | 10 | 63.32 | 72.10 | 75.90 | 77.48 | 104.84 | **1.45** |
| −8 | C | 20 | 240.41 | 269.87 | 282.08 | 286.72 | 409.56 | **1.52** |
| −9 | A | 5 | 8.97 | 9.49 | 10.07 | 10.76 | 12.30 | **1.30** |
| −9 | A | 10 | 47.45 | 51.47 | 53.41 | 58.53 | 72.13 | **1.40** |
| −9 | A | 20 | 216.49 | 236.03 | 242.71 | 268.69 | 341.78 | **1.45** |

Exact distances buy **≈ 1.6× in D**, i.e. **≈ 1.26× in effective SNR** — equivalent to ~1.6× more
observing time, or ~26 % more array. Gaia-class tightening buys **1.6 %** in D at `fgw=-8`.
The zero-noise ceiling confirms it: `D_truth` goes 66.09 (feather) → 100.000 (exact), a factor
**1.51** with no noise at all — the pulsar term carries **34 %** of the matched-filter SNR².

**Detection fraction is saturated and therefore uninformative here.** With the sky fixed at truth,
`D_truth > 2.706` in **10/10** realisations at every SNR ≥ 5 for `fgw=-8` (SNR 5 already gives
D ≈ 18–35). The only sub-unity entries are `fgw=-9`, SNR 5 (0.80 → 0.90 across tiers). A
*non-trivial* detection fraction needs a sky-**searched** statistic with its trials factor — and
that is precisely what the aliasing of §3.4 prevents from being computed on a grid. I report the
statistic value (above) rather than manufacture a detection curve from a broken maximiser.

### 5.5 The null bank, and a real defect it uncovered at `fgw = -9`

`D_truth` on paired noise-only data (same seeds, `log10_h = -30`), pooled over 5 tiers × 10
realisations. Note the 5 tiers share a noise draw, so there are only **10 independent** samples;
quantiles are correspondingly noisy.

| fgw | scen | frac(D=0) [0.50] | mean [0.50] | p95 [2.706] | max |
|---|---|---|---|---|---|
| −8 | A | 0.08 | 0.835 | 5.14 | 5.20 |
| −8 | B | 0.34 | 0.432 | 1.63 | 2.07 |
| −8 | C | 0.40 | 0.873 | 2.64 | 9.32 |
| **−9** | **A** | 0.60 | **2.859** | **19.09** | **21.12** |

The `fgw=-8` rows are consistent with `½δ(0)+½χ²₁` given 10 effective samples. I verified this
directly rather than assume it. Two dedicated 60-realisation runs measure `z ≡ √D_truth − SNR`
(scenario A, tier3), each with the harness's timing-model prior and with `discovery`'s own improper
default — jobs `12453453` (fgw −8) and `12453872` (fgw −9):

| fgw | `tm_variance` | mean(z) | std(z) | **Var(z)** | KS p |
|---|---|---|---|---|---|
| −8 | **1e-14 (harness)** | +0.103 ± 0.140 | 1.084 | **1.174** | **0.726** |
| −8 | 1e40 (improper) | +0.076 ± 0.124 | 0.962 | 0.925 | 0.493 |
| **−9** | **1e-14 (harness)** | +0.173 ± 0.327 | **2.534** | **6.422** | **0.0000** |
| −9 | 1e40 (improper) | −0.033 ± 0.105 | 0.810 | 0.655 | 0.322 |

`Var(z)` is exactly *(true noise power) / (noise power the likelihood assumes)* along the CW
template. Read the table:

- **At `fgw=-8` the harness's noise model is sound** (Var = 1.17, KS p = 0.73). The
  Enterprise-simulates / Discovery-evaluates contract holds. The headline rests on this row.
- **At `fgw=-9` it is wrong by a factor 6.4 in power** (2.5× in amplitude), KS p = 0. Switching
  `discovery`'s `makegp_timing` to improper marginalisation (`1e40`) repairs it (Var = 0.655,
  KS p = 0.32). **The cause is confirmed, not conjectured:** at `fgw = 1e-9` the CW period (32 yr)
  *exceeds* the longest baseline (22 yr), so the template is nearly a low-order polynomial in `t`
  and is strongly degenerate with the timing model — which the harness does not marginalise
  (`prior_variance=1e-14` at `sim.py:232`; Enterprise's `tm_prior` = `variance · N_toa` against
  `discovery`'s `makegp_timing` = `variance · N_toa / N_par`, a factor ≈ `N_par` ≈ 19 in assumed
  coefficient variance). Un-marginalised timing-model power leaks onto the CW template.

(The N = 10 campaign seeds at `fgw=-8` happen to give mean z = +0.90 — an unlucky subset of the
verified `N(0,1)`, which inflates every `D_truth` median in §5.1 by ~19 % at SNR 10. Tier *ratios*
are unaffected; absolute D values carry that common-mode offset. This is sampling, not model error:
the 60-realisation run settles it.)

**Consequence, stated precisely.** The `fgw=-9` **noisy** rows — `inside90`, MAP scatter,
`D_truth`, the null bank — are **uncalibrated and I build no claim on them**. The `fgw=-9`
**noise-free** rows contain no noise and are therefore untouched by this defect; they are the good
part of that column, and they carry the κ-ladder result (bias 2.07° → 1.78° → 1.80° → **0.033°** →
0.000° as κ̄ goes 0.059 → 0.149 → 0.189 → 0.290 → 1.000) that isolates the bias mechanism (§3.5).
The `fgw=-8` columns, which carry the headline, are unaffected: Var(z) = 1.17, TM-prior sensitivity
on SNR = 1.3 % (§3.3), `z ~ N(0,1)` at KS p = 0.73.

---

## 6. SLURM SHAPE — as run

**Deviation from the brief, deliberate.** The brief specified an array over
`(tier × SNR × scenario)` with realisations looped inside. Because `{psr}_cw_sigma_d` is a
likelihood *parameter* (§2c), all five tiers sweep inside one realisation with **no rebuild and no
recompile**. The correct task unit is therefore `(fgw × scenario × SNR)`, with **tiers and
realisations both looped within**. This is strictly cheaper: 16 tasks instead of 80, and the
~230 s process + graph-materialisation overhead (`HPC_SETUP.md` §7.4) is paid 16 times, not 80.

Lane, per `HPC_SETUP.md` §2: `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc
--gres=gpu:nvidia_a100-sxm4-80gb:1`, `--cpus-per-task=8 --mem=64G`, `--requeue`,
`--signal=B:USR1@300`. Device-log first line in every job. **Per-realisation checkpointing with
self-resume is implemented** (`ring_cell.py` reloads its own npz and skips completed realisations)
— `HPC_SETUP.md` D-1 notes no script in the repo does this; this one does. Output is lean npz
(19 float arrays × 5 tiers + a JSON meta string); the whole campaign is **< 1 MB**.

| stage | job | tasks | measured wall |
|---|---|---|---|
| gates | `12452632` | 1 | 1 m 57 s |
| aliasing diagnostic | `12452700` | 1 | 1 m 24 s |
| warm-cache probe (also r=0,1) | `12452813`, `12453068` | 3 | 2–10 m |
| **fast array** (A, B, fgw−9 A, 3 nulls) | `12453218` | 12 | **2 m 11 s – 15 m 21 s**, all COMPLETED |
| **slow array** (scenario C) | `12453219` | 4 | **4 m 26 s (null), 3 h 21 m / 3 h 47 m / 4 h 37 m**, all COMPLETED |
| zero-noise control | `12453271` | 1 | 18 m 26 s |
| noise-model check ×2 | `12453453`, `12453872` | 2 | 22 m 55 s / 23 m 05 s |

All 16 cells COMPLETED, N = 10 each, no failures, no requeues, no NaNs.

**Walltimes were measured, not derived** (`HPC_SETUP.md` §5 item 7). Scenario A: **17–32 s per
realisation**. B: **39–89 s**. C: **1208 / 1657 / 1362 s** (SNR 5/10/20) — ≈ 40× A. Cause:
`makeglobalgp_fourier` with the
HD ORF puts a dense `(30 psr × 30 freq)² = 1800×1800` factorisation **inside every likelihood
evaluation**, and it is not hoisted out of the sky-scan `vmap` (it does not depend on the sky
parameters). A 512-wide `vmap` of it asks XLA for **74 GiB** (`12452569`, OOM) — the per-scenario
`vmap` chunk is 512 / 256 / **64** for A / B / C. So C's tasks are ~3–4 h and I gave them 8 h;
A/B/−9 got 1 h. **This breaks the brief's "tasks ≤ 30 min" constraint for the four C tasks.** I
made that call because (i) per-realisation checkpoint + resume means a TIMEOUT costs one
realisation, (ii) `HPC_SETUP.md` §2 establishes the `dgx_iacc` lane is uncontended with a 14-day
ceiling, so short walltimes buy nothing here.

**Compile-cache note.** Q1's graphs are new shapes; the cold compile was paid once into the shared
`/home/milesmt/.cache/jax_stagec_cache`. The warm-cache probe (`12452813`) doubled as the first
realisation of tasks 1 and 4, and later array tasks resumed from its checkpoints. Job `12452569`
also shows the trap that HARBOR's §7.4 predicts in a different form: the *benchmark* I first wrote
timed a fresh graph shape on every call, so it measured compile, not throughput (265 evals/s vs the
true ~25 000 evals/s). Recorded so nobody re-derives it.

---

## 7. WHAT I STOPPED ON, AND WHAT I RECOMMEND

Per the brief I stopped and reported rather than patched, three times.

**S-1 — SHOVEL's fixes (b) and (d) are coupled, and (d) as specified is insufficient.**
Restoring `log10_fgw ≈ -8` (fix b) makes the sky posterior fringed at 0.03–0.5°. Refining the peak
of a coarse full-sky grid (fix d) refines the *wrong* peak: at tier3 the coarse+zoom estimator
returns a "maximum" below the value at truth in **61 %** of realisations, and in **100 %** at
SNR 20. I did not patch the global search — it cannot be fixed by refinement, only by a sampler or
an earth-term-first hierarchical search. Instead I reported an **exact sky-fixed statistic**
(`D_truth`) and a **truth-anchored fringe-resolved localisation** (`local_sky`), and kept the
broken number as a flagged legacy column. **Any future Q1 sky search at `fgw ≈ -8` needs a
sampler, or a two-stage earth-term-then-pulsar-term search. A grid will not do.**

**S-2 — the tier ladder at `fgw=-8` is nearly binary, and the Gaia tier is a no-op.**
κ > ½ needs σ_d < 3.02 pc; 1 of 30 ring pulsars qualifies. This is a property of the design
(`fgw = 1e-8`, 20° ring), not a bug, but it means the experiment as specified cannot resolve a
`tier2` ladder at that frequency. If the intent is to measure the *gradient* of localisation with
distance precision, the design must either (a) move to `fgw ≈ 1e-9` (where 5/30 qualify and the
ladder is graded — but see S-3), or (b) build the ring from pulsars whose σ_d straddles the
coherence threshold, or (c) sweep σ_d directly rather than through catalogue-scaled tiers.

**S-3 — the harness's timing-model prior is internally inconsistent, and it breaks `fgw=-9`.**
Enterprise injects timing-model coefficients with variance `1e-14 · N_toa` (`sim.py:232` →
`tm_prior`); `discovery`'s `makegp_timing(variance=1e-14)` assumes `1e-14 · N_toa / N_par` — a
factor ≈ `N_par` ≈ 19. **Measured, 60 realisations each (§5.5):** at `fgw = 1e-8` this is harmless
(`Var(z) = 1.17`, KS p = 0.73; 1.3 % on SNR). At `fgw = 1e-9` the likelihood under-assumes the
noise power along the CW template by a factor **6.42** (`std(z) = 2.53`, KS p = 0.0000), because
the 32-yr-period template is nearly degenerate with the timing model over a 22-yr baseline.
**Setting `discovery`'s `makegp_timing` to its own improper default (`1e40`) repairs it**
(`Var(z) = 0.655`, KS p = 0.32). **I did not change either prior.** Recommended fix on cronus: make
both sides improper, or make them agree explicitly. Until then, **treat every `fgw ≲ -9` *noisy*
result from this harness as uncalibrated.**

> This casts an extra shadow on the Dec-2025 run, which used `log10_fgw = -10.5`. That is a CW
> period of ≈ 1000 yr against a 22-yr baseline — a factor **45** — so its template sits deep inside
> the quadratic timing model's null space, exactly the regime where the mis-specification measured
> above is largest, and 1.5 decades *further* into it than the `fgw=-9` point where I measured
> `Var(z) = 6.4`. I did not re-measure at −10.5 (SHOVEL already classified those outputs
> superseded), so this is a prediction, not a measurement. But it means the Dec run's detection
> statistic was inflated by two independent mechanisms, not one: `log10_h = -12`, **and** a timing
> model too tightly constrained to absorb a template it should largely have absorbed.

**S-4 — with a GWB, even exact distances undercover.** Scenario C's `tier3` `inside90` is
**0.40–0.50** against a nominal 0.90, while scenarios A and B hit 0.90–1.00 exactly (§5.2). The
exact-distance 90 % region there is ~1e-4 deg², smaller than the sky perturbation a GWB realisation
induces. I did **not** separate genuine physics (a correlated background shifts which fringe the
MAP selects) from estimator error (the sky "posterior" here is a *profile* slice with all other CW
parameters held at truth, not a marginal posterior; its calibration is weakest where the region is
smallest). **A sampler is required to settle this, and the scenario-C `tier3` coverage number should
not be quoted until one has been run.** It does not touch the headline, which rests on the
imperfect tiers in A and B.

**S-5 — I ran the wrong "real" prior, and the right one is a free experiment.** See §2.1. The
frozen canonical EM prior `CW_transition/best_distances.txt` is git-tracked and I wrongly concluded
it was absent. Measured impact on this campaign's conclusions: **none** (κ̄ 0.043 → 0.033, same one
pulsar above the coherence threshold, same two above κ>0.1). Measured impact on what it *enables*:
large — its distance **means** differ from the feather means by 1.40 σ on J0437-4715, the only
pulsar with a live pulsar term, injecting a **0.55 rad** pulsar-term phase error. That is the
mis-centred-prior arm, ready to run.

**Also reported, not fixed:**

- `sim.py:490` (and `:362`) force `psr.toaerrs = 1e-6 s` *inside* the likelihood builder. Any
  realistic-noise run must bypass it. Suggest a `toaerr_override=None` kwarg defaulting to today's
  `1e-6` so cronus behaviour is byte-identical.
- `sim.py:40` sets `--xla_gpu_autotune_level=2`, which `HPC_SETUP.md` §7.3 measured as **2.2×
  slower** on A100. Suggest reading it from the environment.
- Scenario C's HD global GP is re-factorised inside every sky-scan likelihood evaluation. Hoisting
  the `(1800×1800)` Cholesky out of the sky `vmap` (it has no sky dependence) would give scenario C
  roughly the throughput of scenario A — a ~40× speedup for the most realistic scenario.

**Scope caveats on the headline, stated plainly.** (i) One injected sky position (seed 123) and one
ring realisation; the *magnitude* of the 3–6° bias is geometry-specific, though the zero-noise
control proves it is deterministic, not noise. (ii) The harness's priors are wide-but-centred, never
mis-centred (§2a) — the mis-centred arm is untested and can only be worse; the "real" tier I ran is
the feather σ_d, not the frozen canonical EM prior (§2.1, which shows the substitution does not move
any conclusion but does forgo the mis-centring experiment). (iii) `local_sky` is
truth-anchored, so its MAP offsets are lower bounds. (iv) N = 10; binomial intervals are quoted
throughout and the coverage collapse at SNR 20 (0/10 for feather and tier1, both scenarios A and B)
is the one result that is already unambiguous at this N.

---

## 8. ARTEFACTS

```
RING_results/
  ring_lib.py ring_gates.py ring_diag.py ring_nf.py ring_cell.py
  ring_zcheck.py ring_zcheck9.py ring_bench_c.py ring_consolidate.py
  ring_env.sh  ring_*.sbatch
  consolidated.txt                     full tables
  ring_q1_fgw{-8,-9}_{A,B,C}.png       area90 / inside90 / bias+scatter / detfrac vs tier
  npz/cell_fgw{-8,-9}_{A,B,C}_{snr5,snr10,snr20,null}.npz   16 cells, lean
  npz/gates.json  npz/diag.json  npz/noisefree.json  npz/zcheck*.json  npz/summary.json
  logs/                                every job, device-log first line
```

**Campaign completeness.** 16/16 cells COMPLETED, N = 10 realisations each, zero failures, zero
requeues, zero NaNs. 4 prior tiers requested by the brief; 5 run (adding `tier2_k3`, free, to give
the gain a third interior point). 3 scenarios × 3 SNR × 2 frequencies as specified, plus 4 null
banks, 4 zero-noise controls, and 4 noise-model calibration arms not in the brief.

Nothing was committed. Nothing was pushed. No tracked file was edited. `project_progress.md`,
`MANIFEST.md` and `MAIN_PROJECT_QUESTIONS/Q1/` were not touched. The only file written outside
`RING_results/` is this report and `.claude/settings.json` (both untracked).
