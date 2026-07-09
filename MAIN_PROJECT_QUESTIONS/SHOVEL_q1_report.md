# SHOVEL — Q1 excavation + b20 provenance

*Read-only audit by agent SHOVEL. No reruns, no deletions, no edits to any file
except this one. Context read: `project_progress.md`, `CW_transition/trackB_estimator_spec.md`,
`MANIFEST.md`, plus Q1 tree + result pickles + b20 pickles.*

Repo: `/home/mattm/projects/HSYMT`. Subject: `MAIN_PROJECT_QUESTIONS/Q1/`.

---

## TASK 1 — Q1 excavation

### (a) Inventory

Single banked run: **`results/run_20251218_103717/`** — timestamp = **2025-12-18 10:37:17**
(the `run_YYYYmmdd_HHMMSS` dir name is `datetime.now()` at save time; all file mtimes are
2026-05-04, a later checkout/touch, so the dir name is the true compute date).

**Scaffolding (2026-05-04 mtimes; code):**
- `sim.py` (25 KB) — load feathers, build ring geometry, Enterprise PTA sim → Discovery
  likelihood; `build_discovery_loglike_marginalized()` (analytic distance marginalisation).
- `optimize.py` (6 KB) — MAP ascent / SVI / LBFGS helpers (NOT used by the banked run).
- `optimization_helpers.py` (8.6 KB) — `compute_detection_statistic_optimized()`,
  `LikelihoodCache`, cached pulsar loading. **This is what produced the numbers.**
- `run_q1_map.py` (CLI, one scenario/tier) + `test_optimization.py`.
- Notebooks: `run_q1_det_stat.ipynb`, `run_q1_det_stat_fast.ipynb` (**the run driver**),
  `run_q1_det_stat_fast_scenB.ipynb`, `run_q1_map.ipynb`, `run_q1_sample.ipynb`.
- Docs: `CLAUDE.md`, `README_Q1.md`, `OPTIMIZATION_{GUIDE,SUMMARY}.md`.

**Outputs in the run dir:**
- Data: `metadata.pkl`, `results.pkl` (21 rows), `ensemble_D.pkl`, `injected_positions.pkl`,
  `posterior_cache.pkl` (5.2 MB).
- Figures: `dstat_{combined,per_snr,trends}{,_no_tier1}.{png,pdf}`,
  `sky_area_vs_tier.{png,pdf}`, `sky_maps/skymap_{real,tier1,tier2_k_*,tier3}.{png,pdf}`.
- Root of Q1/: `D_vs_k_snr{5,10,20}.png`.

### (b) What was computed

**Question (as coded):** Q1 forward direction — *how do CW detection statistic & sky
localisation change with the pulsar-distance prior width?* Swept distance-prior tiers × nominal
SNR, single realisation.

**Config (`metadata.pkl`):** `scenario='A'` (WN only), `n_pulsars=30`, `seed=123`,
`snrs=[5,10,20]`, `tiers=[real, tier1, tier2_k=1.6, tier2_k=2, tier2_k=3, tier2_k=5, tier3]`,
sky grid `n_cos=n_phi=120`, `n_amp=8`, `amp_span=1.5`, **`log10_fgw_override=-10.5`**.

**Tiers** (`apply_distance_prior_tier`): real = native feather σ_d; tier1 = σ×5 (loosen);
tier2_k = σ/k (tighten, k∈{1.6,2,3,5}); tier3 = σ→1e-6 (≈exact).

**Metric machinery:** likelihood = `build_discovery_loglike_marginalized` (distances integrated
out analytically via coherence factor κ = exp(−½[2ω₀(1−cosμ)(kpc/c)σ_d]²)). Detection statistic
= brute-force profile over a (cos_gwtheta, phi, log10_h) grid, max over amplitude:
`D = 2·(logL_max − logL_null)` (`optimization_helpers.py:156`). Sky area = 90%-mass credible
region on the profiled logL grid; `inside90` = is the true pixel inside it.

### (c) ERA-DATE vs the correction timeline — **CLEAN of the bug era, but realism-broken**

- **`find_best_wrong_mode_in_prior` / prominence peak-finder?** **NO.** Grepped all Q1 `.py`
  and the driver notebook — zero references. Q1 lives in a separate code path
  (`sim.py` + `optimization_helpers.py` + `discovery`), never imports `CW_lnL_check/cw_helpers.py`.
  It is **immune to the prominence-floor-0.5 bug** (`project_progress.md` D3/D5, P2-D).
- **Conditional joint-ΔlnL metric (`compute_joint_best_wrong_in_prior`)?** **NO.** Q1 uses a
  proper profile-likelihood `2·(logL_max − logL_null)` with an explicit no-signal null —
  the *marginal, prominence-free* style the tracker recommends over the deprecated conditional
  metric (`project_progress.md:340-354`). **Metric methodology is HEALTHY.**
- **Injection strain / the 56× realism issue?** **This is where Q1 fails — and worse than 56×.**
  1. **No SNR calibration.** "SNR" is a decorative label. The only place it enters:
     `cw_draws[0]["log10_h"] += np.log10(max(snr,1)/10.0)` on top of a hardcoded
     `log10_h=-12.0` anchored at "snr=10" (`sim.py:131`, notebook `simulate_data`). No
     matched-filter SNR is ever computed (README/CLAUDE both list SNR targeting as an
     unimplemented "next step"). Result: labels 5/10/20 just scale strain ×0.5/×1/×2.
  2. **Signals are ~10⁴× too loud.** With TOA errors force-set to 1e-6 s (`set_toaerrs`,
     `sim.py:172,362,490`) the banked `D_stat` ranges **1.7e10 → 2.8e11** (`results.pkl`),
     i.e. effective optimal SNR ~10⁴–10⁵, not 5–20. The other prong's documented "56×
     too loud" (`project_progress.md:795,1320`) is a mild version of the same disease here.
  3. **Wrong frequency kills the mechanism.** Design calls for `log10_fgw ≈ -8`
     ("evolution happens but isn't insane", `main_project_questions_readme.md:25`, `CLAUDE.md:27`).
     The run overrode to **`log10_fgw=-10.5`** → ω₀ tiny → pulsar-term phase barely evolves →
     κ≈1 regardless of σ_d → **the distance-degeneracy-breaking effect that Q1 exists to
     measure is absent.** Consequence visible in the numbers: D changes only ~3% across
     real→tier3 at fixed SNR (1.73e10 → 1.78e10 at SNR5).

### (d) Headline numbers (with context)

From `results.pkl` (single realisation, seed 123, scenario A, log10_fgw=-10.5):

| SNR | tier | D_stat | inside90 | area90 (deg²) |
|----|------|--------|----------|---------------|
| 5  | real | 1.735e10 | ✔ | 2.9131 |
| 5  | tier1 (σ×5) | 1.567e10 | **✘** | 2.9131 |
| 5  | tier3 (exact) | 1.781e10 | ✔ | 2.9131 |
| 10 | real | 6.940e10 | ✔ | 2.9131 |
| 10 | tier1 | 6.269e10 | **✘** | 2.9131 |
| 10 | tier3 | 7.123e10 | ✔ | 2.9131 |
| 20 | real | 2.776e11 | ✔ | 2.9131 |
| 20 | tier1 | 2.508e11 | **✘** | 2.9131 |
| 20 | tier3 | 2.849e11 | ✔ | 2.9131 |

Patterns:
- **`area90_deg2` is IDENTICAL (2.913138990849465) in all 21 rows.** The 1e-6 s TOA forcing
  makes the profiled posterior a near-delta spike → the 90% region floors at grid resolution
  → the localisation metric is **resolution-limited and uninformative** (it cannot move with
  SNR or tier). `sky_area_vs_tier.png` is a set of flat lines.
- **D grows ~×16 from SNR5→SNR20** — consistent with D∝strain², i.e. purely the label's
  amplitude scaling, not a physical detection curve.
- **D barely moves across tiers (~3%)** and **non-monotonically**: loosening (tier1) lowers D,
  tightening/exact raises it slightly — directionally sensible but tiny, because the low
  frequency suppresses the distance term.
- **`inside90` is False only for tier1, at every SNR.** The one qualitatively interesting
  signal: loosening the distance prior pushes the MAP >1 pixel off truth. But on a
  resolution-floored area with N=1, it is a coarse yes/no, not a calibrated statement.
- **`ensemble_D.pkl`**: D-vs-k sweep (k = σ_d scale, 0.05→3.0, 25 pts) is **n_real=1** despite
  the "ensemble" name; the ±σ bands in `D_vs_k_snr*.png` are kernel-smoothing artifacts, not
  realisation scatter. Curve shows a shallow (~7%) peak near k≈0.66 then decline.
- **Scenario coverage:** only **A (WN)** was run. B (WN+RN) / C (WN+RN+GWB) never executed
  (a `_scenB` notebook exists but produced no banked run). Design wanted A→B→C.

### (e) VERDICT — **MIXED**

- **HEALTHY / resume-and-modernize (the code + method):** the scaffolding
  (`sim.py`, `optimize.py`, `optimization_helpers.py`), the analytic distance-marginalised
  likelihood, the tier-application machinery, and the **`2·(logL_max−logL_null)` profile-grid
  detection statistic** are methodologically sound and **entirely free of the prominence-bug /
  conditional-metric contamination**. This is a reusable forward-Q1 harness. → The queued
  inverse-Q1 backlog item can **resume on this code**, not restart.
- **CONTAMINATED / classify-superseded (the banked numbers + figures in
  `run_20251218_103717/`):** not citable, because of four compounding realism defects —
  (1) uncalibrated "SNR" (labels ≠ physical SNR; signals ~10⁴× too loud),
  (2) `log10_fgw=-10.5` vs design −8, which removes the phase-evolution mechanism the study
  measures (tier sensitivity washed to ~3%),
  (3) resolution-floored constant `area90` (localisation metric dead),
  (4) `n_real=1` and scenario A only.

**One-line:** *Clean, reusable harness; the specific Dec-2025 run's outputs are superseded.*
The MANIFEST's "register Q1 as a new §; reconcile with inverse-Q1 backlog" stands — with the
note that only the **code** graduates; the **banked figures do not**.

**Modernize checklist before any citable rerun:** implement real matched-filter SNR targeting
(root-find log10_h → target SNR); restore `log10_fgw≈-8`; use realistic per-pulsar TOA errors
instead of 1e-6 s; refine sky localisation off-grid (Hessian/zoom) so area90 is not
pixel-floored; run N realisations; execute scenarios B and C.

---

## TASK 2 — b20 provenance

Files (both git-tracked, `data_products/`):
- `b20_cw_curn_r0.pkl` — 20.0 MB, mtime **2025-10-04**, committed `e8ca312` 2025-10-23.
- `b20_cw_curn_r0_w_flags.pkl` — 21.0 MB, mtime **2025-10-14** (matches the frozen feather-set
  date in MANIFEST), newer numpy (`numpy._core`), `_flags` populated.

**Structure (opcode walk, no full unpickle — safe):** top level is a **Python `list`**
(`EMPTY_LIST` first op) of **116 objects** (`NEWOBJ ×116`), each an
**`enterprise.pulsar_edited.Tempo2Pulsar`** (custom fork class). 2090 `REDUCE` / 1974 `BUILD`
(numpy arrays + object state). Per-object attributes seen: `name`, `_toas`, `_stoas`,
`_residuals`, `_toaerrs`, `_flags`, `_designmatrix`, `planetssb`, `_pos`, `pdist`. 113 distinct
`Jhhmm±ddmm` pulsar names embedded (regex undercount of the 116 instances by a few edge cases).

**Reading of the name "b20":** **a build/batch tag, NOT "20 pulsars."** The list holds the
**full 116-pulsar array** — the same array as the 116 `data_products/*.feather` files.
`cw` = injected continuous-wave source; `curn` = common uncorrelated red noise; `r0` =
realisation 0.

**Best provenance guess (HIGH confidence):** a **bundled Enterprise-simulated dataset of the
full 116-pulsar array with an injected CW + CURN, realisation 0**, generated Oct 2025 during the
`hessian_check_enterprise_edited` lnL investigations (the commit that introduced it,
`45882ad`, is that workstream). `_w_flags` is a re-dump one build later (2025-10-14) that
carries backend/system `_flags` and was serialised with a newer numpy. Consistent with
MANIFEST §E ("20-pulsar CW+CURN bundled sim") **except the count: it is 116, not 20** — the
"b20" tag misleads.

**No doc ties any banked result to these pickles** (confirmed via MANIFEST §B/§E). Provenance is
inferred from structure + git, not from a doc reference. **Nothing deleted or modified.**

**VERDICT (b20): BENIGN DATA ARTIFACT, provenance now identified** — keep. Suggest labelling it
in the tracker as "full-116-array CW+CURN Enterprise sim, r0, Oct-2025" and flagging that the
`b20` prefix does not mean 20 pulsars.
