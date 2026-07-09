# Speedup brief — Track B P1/P2 compute hot-path (for a helper agent)

**Goal:** make the Track B fringe/registration experiments ~100–1000× faster. Everything
below runs in the `jug-gpu` conda env (`/home/mattm/miniforge3/envs/jug-gpu/bin/python`,
jax 0.10.1) on a SHARED RTX 4090 (keep `XLA_FLAGS=--xla_gpu_autotune_level=0`; a co-tenant
`kyleg` spikes GPU memory and OTHER-user autotune profiling OOMs the compile). `discovery`
is editable at `/home/mattm/soft/discovery` (jax port; `src/discovery/matrix.py`).

## 1. The one hot primitive (this is the whole cost)

`fisher_logL(theta_arr, data_tuple)` in `CW_transition/stagec_fisher.py` (~line 177),
returned by `build_fisher_amortised(...)`. It is the FULL 116-pulsar array marginalised
CW log-likelihood (white + frozen GWB GP + frozen per-pulsar red-noise GP marginalised,
Woodbury). `theta` = [116 distances] ++ [8·N_CW source params] (n_theta=244 for N_CW=16).
Measured: **0.544 ms per eval** warm; ~465 s one-time compile (cached at
`~/.cache/jax_stagec_cache`). `fisher_logL_batched = jit(vmap(fisher_logL,(0,None)))`.

## 2. Why every experiment is slow

The E-step (`batched_scan` in `CW_transition/trackB_p1_map.py`) evaluates, for ONE source
position, each pulsar's distance swept over B=512 fringe centres, others fixed:
**116 × 512 = 59392 full-array evals ≈ 32 s per source position.** Every test is
N_source_positions × 32 s:
- P1 map: 169 trials/stage × ~44 min/stage.
- probes / line scan / gate: dozens–hundreds of source positions × 32 s each.
- The soft M-step (`make_soft_mstep` in `trackB_p2_pipeline.py`) is a scan over ~npsr·KTOP
  mixture configs, one full-array BACKWARD (jax.grad) per config, ×40 Adam steps ×30 EM iters.

All of this is bound by the per-fringe full-array eval. Only ONE pulsar's distance changes
between adjacent fringe evals, yet the whole array likelihood is recomputed.

## 3. The fix — low-rank QuickCW split (algebra ALREADY VERIFIED vs the code)

See `CW_transition/trackB_estimator_spec.md` §6 / §6.1 for the full derivation (traced through
`discovery/matrix.py:make_kernelterms_vary`, line ~2631). Summary:

`fisher_logL` is a TWO-LEVEL Woodbury. Distances enter at exactly two PER-PULSAR-LOCAL places:
```
fisher_logL(theta) = Σ_p a_p(L_p)  +  ½[ b(theta)ᵀ M⁻¹ b(theta) − ldP − logdet_M ]
  a_p(L_p) = −½( yᵀN⁻¹y_p − FtNmy_rn,p·cf_p⁻¹FtNmy_rn,p ) [+ frozen consts]   (per-pulsar scalar)
  b_p(L_p) = TtNmy_gwb,p − TtNmF_p·(cf_p⁻¹ FtNmy_rn,p)                        (n_gwb vector)
  M ≡ cf_cached = matrix_factor(Φ_gwb⁻¹ + blockdiag(c))    -- FROZEN (GWB HD-ORF coupling)
```
`c` (terms[2]) is distance-INDEPENDENT → `M` built once. Sweeping pulsar p's fringe changes
ONLY `a_p` (scalar) and `b_p` (n_gwb vector); every other pulsar's block and `M` are fixed.
So with `u ≡ M⁻¹ b_base` (one solve per source position) and the per-pulsar blocks `(M⁻¹)_pp`
(**precompute ONCE per experiment — M is frozen**), each fringe eval is a rank-block update:
```
Q(L_p)  = Q_base + 2·Δb_p·u_p + Δb_pᵀ (M⁻¹)_pp Δb_p ,   Δb_p = b_p(L_p) − b_p(base)
lnL(L_p) = lnL_base + [a_p(L_p) − a_p(base)] + ½[Q(L_p) − Q_base]
```
Per-fringe cost collapses from a full (npsr·n_gwb) triangular solve + 116-pulsar delays to
ONE single-pulsar kernelterm (`yprod_p`, matrix.py ~2646) + two n_gwb dot products — ~1000×.
**The chirp is retained** (b_p re-evaluates pulsar p's residual incl. frequency evolution over
L/c — the ~1e-3-nat fringe-breaking signal), only the array coupling is amortised.

### What to build
A function `estep_split(P, source_params) -> LNL[pulsar, fringe]` that reproduces
`batched_scan(P,[theta])[0]` **bit-exact to ~1e-8** but via the rank-block update:
1. once per experiment: extract/pre-solve `(M⁻¹)_pp` for all p (or only the carrier subset).
2. per source position: one array kernelterm for the base `a_base,b_base`; `u=M⁻¹ b_base`.
3. per (p, fringe): recompute `a_p(L_p), b_p(L_p)` (single-pulsar `yprod_p`), apply the update.

Needs reaching into `discovery`: `fml.vsm` exposes `Ns` (per-pulsar noise solves), `Fs`
(rednoise basis), the gwb basis `Ts`, and `P_var` (common RN); `matrix_factor/solve` are in
`discovery.matrix`. The per-pulsar `yprod` closure (matrix.py:2646) already returns
`(ytNmy_p, FtNmy_rn_p, TtNmy_gwb_p)` — reuse it rather than re-deriving.

### Gate (mandatory before trusting it)
`estep_split` vs `batched_scan` on ~8 random source positions, max|Δ| < 1e-8, on the pop
config (N_CW=16, 116 psr). Reuse `trackB_p1_map.batched_scan` as the reference.

## 4. Where it plugs in (all callers of the E-step)
`trackB_p1_map.batched_scan` (P1 map), `trackB_p2_probe.py` / `trackB_p2_probe2.py`
(probes), `trackB_p2_linescan.py` (line scan), `trackB_p2_pipeline.py` (soft-EM gate +
M-step). Swapping `batched_scan` for `estep_split` accelerates every one. The soft M-step's
grad also benefits: with a_p,b_p closed-form, the source-gradient of each mixture term is a
small analytic expression instead of a full-array reverse-mode pass.

## 5. Lower-value but easy wins
- `build_problem` rebuilds `amo_full` AND `amo_earth` (~2 compiles); P1/probes only need
  `amo_full`. Skipping `amo_earth` saves one compile per script.
- Batched source-position axis: many tests evaluate one source position at a time
  (`batched_scan([theta])`); batching M positions into one vmapped call improves occupancy.
- KTOP / #fringes B: B=512 per pulsar is the science default; the SEARCH phases can use
  fewer fringe candidates (the window is dominated by a handful of near-fringes) — but keep
  B=512 for any CERTIFICATION eval.

## 6. Constraints / gotchas
- x64 required (`jax_enable_x64`). Zero-noise Asimov: residual(truth)=0 ⇒ H=−Fisher exactly.
- The M-step grad OOMs (32 GiB) if you vmap-all mixture configs; use a scan/accumulate (grad
  of a sum = sum of grads) — already done in `trackB_p2_pipeline.make_soft_mstep`.
- Persistent compile cache is on; first compile of any new shape is 100s–500s (one-time).
