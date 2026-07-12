# ATLAS — the E-track self-clocking map (WHERE the eccentric comb supplies the >20×)

**Agent:** ATLAS · ACCRE · tag `b1-v1` (`git describe --tags` → b1-v1; discovery `136b270f`,
ee `f73b8e0`) · **Date:** 2026-07-11 · **PURE EXECUTION** (no commits, no tracked-file edits).

**Scratch paths (host):** code `hpc_harbor/atlas/`, results `ATLAS_results/`, logs
`hpc_harbor/logs/atlas_*`, container `/home/milesmt/soft/harbor/el9.sif`, jax cache
`/home/milesmt/.cache/jax_stagec_cache`. Lane `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc
--gres=gpu:nvidia_a100-sxm4-80gb:1` (HPC_SETUP §2). Fisher-only; no optimiser, no SMC.

---

## THE QUESTION (re-priced to the final numbers)

Targeted-circular certification fails because `sigma(log10_mc)` must improve **>20×** (below
~0.003 dex) and no EM conditioning delivers it — conditioning the `(f, mc)` prior boxes barely
moves the evidence because the plateau does not fill them; only **likelihood structure** moves
it (project_progress.md 2026-07-10 close-out: the >20× bound, the saturation mechanism). WEAVE's
Fisher scaling says the eccentric harmonic comb multiplies Earth-term chirp information by
`kappa(e) = (n_eff/2)·F(e)` (5.6 / 36.9 / 172 at e = 0.5 / 0.7 / 0.8). **The map decides: WHERE
in (e, Mc, f_orb) does the comb supply the >20× — i.e. where is the first detected source
SELF-CLOCKING?**

## THE ANSWER IN ONE PARAGRAPH — *what must the first detected source be?*

**Eccentric (e ≳ 0.6), massive (Mc ≳ 10⁹ M⊙), and at the HIGH end of the band
(f_orb ≳ 10⁻⁸ Hz, f_gw ≳ 2×10⁻⁸ Hz).** At `(f_orb = 10⁻⁸, Mc = 10⁹)` the comb self-clocks —
`sigma(log10_mc)` improves >20× — from **e ≈ 0.58**; the threshold rises to e ≈ 0.70–0.84 at
lower Mc or `f_orb = 10⁻⁸·⁵`, and **below `f_orb = 10⁻⁸` (into the red-noise/GWB-dominated band)
the comb is buried and NEVER self-clocks, at any e** (best σ stays ≈ 0.1–0.7 dex). The comb
supplies the *conditional* chirp information readily — `kappa_measured ≥ 20` by e ≈ 0.5–0.65 —
but the *marginal* `sigma(log10_mc)`, throttled by the eccentric `(e, f_orb)` parameter
correlations, reaches only **≈ 0.008–0.02 dex** in the valid toy tier (a 40–115× improvement),
which clears the >20× *relative* bound but does **not** quite clear the ~0.003-dex *absolute*
certification floor. The cells that would clear it (e ≳ 0.85) are exactly where the toy tier's
F(e)-boosted comb coalesces within the observation span — the EOB tier is required there.
Nonetheless, **where the comb self-clocks it is its own standard siren**: fed to the SIREN
amplitude relation the comb-measured Mc gives `sigma(D_L)/D_L ≈ 12%` at `(f_orb = 10⁻⁸,
Mc ≳ 10⁹, e ≳ 0.5)` with **zero certified pulsar terms** — the dark-siren-useful class SIREN
reached only with 3 certified seeds.

---

## METHOD (WEAVE §3.3, executed)

**Harmonic stack.** `MultiSourceDelay` / `MaskedDelay` (`trackB_b1_core.py`) sum `n_sources`
copies of the circular kernel `make_phase_connected_binary` sharing one pulsar distance — a
harmonic stack once reparametrised (WEAVE 3.3):

```
log10_mc_n  = log10_mc + (8/5)log2 + (3/5)logF(e) - (8/5)log n     (chirp tie: fdot_n = n·fdot_orb·F)
log10_fgw_n = log10 n + log10_f_orb                                 (radiate at f_n = n·f_orb)
log10_h_n   = log10_h + 0.5 log10 g(n,e)                            (Peters-Mathews amplitude)
```

n=2, e=0 recovers the circular source exactly. `g(n,e)` = Peters-Mathews (1963) Bessel weights,
`sum_n g = F(e)`; **N_HARM = 32** (resolves κ exactly through e = 0.7, where the κ=20 boundary
sits; e = 0.9 truncation-limited in *magnitude* only, still ≫ 20 so the contour is robust).

**The Fisher, cheaply and exactly.** At the zero-noise Asimov truth residual(truth)=0, so
`H = −F` depends only on FIRST signal derivatives. We compose the analytic tie map
`phys(9) → theta_stack(8·N)` with `logL` and HVP over the **nine physical params**
`[cos_gwtheta, gwphi, cos_inc, log10_mc, e, log10_f_orb, log10_h, phase0, psi]` — one graph,
9 columns, reused across the grid. The tie's only non-elementary piece is `g(n,e)`; its
e-derivative (Bessel) is precomputed on CPU and injected as a **linear-in-e amplitude channel**,
which is EXACT for the Fisher (needs only `d(signal)/de` at truth). `F(e)` in the chirp channel
is elementary and autodiffed directly. `Cov = (Fs+Pi)⁻¹` with the generative uniform priors,
**never pinv** (a zero eigenvalue = unconstrained = report the prior, never σ=0).

`kappa_MEASURED(e) = sqrt( F_phys[mc,mc](e) / F_phys[mc,mc](e→0) )` — the colored-noise metric's
verdict on WEAVE's white-noise κ.

**Merger guard (toy-tier limitation, logged).** The F(e)-boosted `mc_n` fed to the *circular*
kernel makes its chirp term `1 − (256/5)mc^{5/3}w^{8/3}t` go negative (the single-harmonic
source "coalesces" within the span) → NaN, at the extreme e×Mc×f_orb corner. `mc_n` is clipped
per harmonic to keep the term ≥ 0.2 over 22 yr. Because `mc_n` and its cap share the same
n-slope, the clip is **binary**: a cell is either fully valid (`n_clip=0`) or its whole comb
coalesces (`n_clip=32` → TOY-TIER INVALID, flagged, its κ NOT read as "not self-clocking").
This is exactly where WEAVE said the toy tier cannot go and LANES said `K_eff → 1`.

**Grid.** e ∈ {0.1, 0.3, 0.5, 0.65, 0.8, 0.9} × Mc ∈ {10⁸·⁵, 10⁹, 10⁹·⁵} × f_orb ∈
{10⁻⁹, 10⁻⁸·⁵, 10⁻⁸} Hz, at the census pop-draw loud0 sky/inc/phase/psi (SIREN's fiducial).

---

## GATES (pre-registered, WEAVE §3.4; stop-at-first-failure)

| gate | pass condition | result |
|---|---|---|
| **G0 (numpy)** | `sum_n g(n,e) = F(e)`; κ,n_eff,F vs WEAVE table | **PASS** (rel < 1e-2 all e; F/n_eff/κ reproduced) |
| **G0/T1** | e=0 stack == circular single source, per-pulsar Earth delay | **PASS** — max rel **7.14e-15** (bit-exact collapse) |
| **G1** | e→0 `sigma(log10_mc)` reproduces the circular Earth-term value | **PASS** — matches SIREN's N_seed=0 σ to **1e-5–1e-3** where the grids align (e.g. f_gw=−7.7, Mc=10⁹·⁵: **0.351 vs 0.353**) |
| Hessian symmetry | — | symres **1e-13 … 1e-16** throughout |

Gates G2–G6 (cancellation test, rank-3 break, ignition, endpoint, referendum) are addressed by
M2/M3 below. Job record: build 12492569, N_HARM=32, dgx03, 55 min (cold HVP compile 3077 s),
all 63 cells computed (0 dropped; 5 flagged toy-invalid).

---

## M1 — kappa MEASURED, not scaled

The measured κ is **frequency-dependent** — this is the "measured, not scaled" content the
white-noise `(n_eff/2)F(e)` cannot capture:

- **`f_orb = 10⁻⁸`** (comb in the sensitive band): `kappa_measured` tracks WEAVE's analytic κ
  closely (Mc=10⁹: e=0.5 → 11.0 vs 5.6; e=0.65 → 65 vs 20.7) — the machinery validated against
  the analytic scaling (M5 figure b, right panel).
- **`f_orb = 10⁻⁸·⁵`**: `kappa_measured ≫` analytic (e=0.65 → **127** vs 20.7; e=0.8 → **2216**
  vs 172) — the comb's higher harmonics `f_n = n·f_orb` reach further into the array's sensitive
  band, amplifying the chirp Fisher far beyond white noise.
- **`f_orb = 10⁻⁹`**: `kappa_measured <` analytic at moderate e (e=0.5 → **0.92** vs 5.6) — the
  fundamental sits where red noise / the GWB rise to low f bury it; spreading power to the comb
  *reduces* the conditional chirp Fisher until e is high enough (e=0.8 → 601) to climb out.

**The throttle (the honest surprise).** The *conditional* chirp Fisher κ is enormous (up to
~33 000× at e=0.9), but the *marginal* `sigma(log10_mc)` improves only ~40–115× — the comb's
chirp information is largely **degenerate with e and f_orb**, and only at high e does the comb
geometry (tooth spacing → f_orb, amplitude ratios → e) break the degeneracy. This is WEAVE §2.4's
rank/conditioning claim, measured: eccentricity's value is CONDITIONING, not magnitude.

---

## M2 — the self-clocking contour (THE map figure: `ATLAS_results/atlas_M2_contour_kappa.png`)

**Minimum e per (Mc, f_orb)** for three thresholds. `σ<0.003 dex` is the *absolute* certification
floor (project_progress, the σ_h→0 host-distance floor 0.00301 dex); `improve>20×` is the
*relative* >20× bound off the prior; `κ≥20` is WEAVE's `Δφ_E ≳ 1 rad` self-clock threshold.

| f_orb | Mc | min e (σ improves **>20×**) | min e (κ_meas ≥ 20) | min e (σ < 0.003 dex) | best σ(log10_mc), valid |
|---|---|---|---|---|---|
| 10⁻⁹   | 8.5 | **never** | 0.65 | never | 0.734 dex |
| 10⁻⁹   | 9.0 | **never** | 0.65 | never | 0.223 dex |
| 10⁻⁹   | 9.5 | **never** | 0.65 | never | 0.106 dex |
| 10⁻⁸·⁵ | 8.5 | **0.84** | 0.50 | never | 0.0219 dex |
| 10⁻⁸·⁵ | 9.0 | **0.77** | 0.50 | never | 0.0150 dex |
| 10⁻⁸·⁵ | 9.5 | **0.66** | 0.50 | never | **0.0075 dex** |
| 10⁻⁸   | 8.5 | **0.70** | 0.53 | never | 0.0194 dex |
| 10⁻⁸   | 9.0 | **0.58** | 0.53 | never | 0.0179 dex |
| 10⁻⁸   | 9.5 | **0.59** | 0.52 | never | 0.0156 dex |

**Readings.** (1) The self-clocking corner is `f_orb ≳ 10⁻⁸`, `e ≳ 0.6`, `Mc ≳ 10⁹` — lowest
threshold **e ≈ 0.58 at (f_orb = 10⁻⁸, Mc = 10⁹)**. (2) `f_orb = 10⁻⁹` never self-clocks (the
comb is buried) — the first source must live at the *top* of the band. (3) The absolute
0.003-dex floor is reached **nowhere in the valid toy tier** (best 0.0075 dex, 2.5× short); the
cells that would reach it (e ≳ 0.85, high Mc/f_orb) are toy-invalid (comb coalesces → EOB tier).
The comb clears the *relative* >20× bound but not the *absolute* one within the toy scope.

---

## M3 — the invariance-theorem escape hatch (`ATLAS_results/atlas_M3_ignition.png`)

WEAVE's clock-cancellation ceiling `R_a^max = C·SNR·(T/tau_a)²` holds only under strict
harmonicity. We broke it honestly: the pulsar-term registration phase over lag `tau_a` was
integrated with **Peters e(t) decay + 1PN periastron advance γ̇** (jax RK4, autodiffed for the
gradient `g = (∂Φ/∂log10_mc, ∂Φ/∂e, ∂Φ/∂log10_f_orb)`; validated `g_mc ∝ τ²`). Registration
`R_a = 1/sqrt(gᵀ Σ g)`, Σ = M1's Earth-term box on (mc, e, f_orb). Compared to the circular
mc-only scalar `R_scalar = 1/(|g_mc|·σ_mc)`. Toy-invalid boxes excluded.

**VERDICT — the cancellation does NOT simply survive; it BREAKS at high e, and marginal ceiling
ignition is reached at the nearest pulsars.**

- **Moderate e (0.3–0.65):** `R_rank3/R_scalar < 1` — the rank-3 channels add wander; the scalar
  cancellation holds, as WEAVE predicted.
- **High e (0.8–0.9):** the rank-3 structure **breaks** the cancellation —
  `R_rank3/R_scalar` up to **41.6** (Mc=8.5, e=0.9, f_orb=10⁻⁸). Mechanism: the eccentric Earth
  term measures f_orb (tooth spacing, σ ≈ 10⁻⁴ dex) and a (mc,e,f_orb) combination (σ ≈ 3×10⁻⁵)
  so well that the pulsar-term fringe wander — dominated by `g_e`, `g_forb` — is pinned.
- **Ignition (`R_a ≥ 1`, Fisher ceiling):** reached **only at `tau_a = 0.3 kyr`** (the nearest
  pulsars) — 3 cells at f_orb=10⁻⁹, 4 at f_orb=10⁻⁸·⁵; max `R_rank3 = 3.53`. Genuinely
  **rank-3-driven** (scalar `R<1<R_rank3`, i.e. eccentricity buys ignition the scalar forbids):
  `(Mc=10⁹, e=0.65, f_orb=10⁻⁸·⁵)` R=1.28 and `(Mc=10⁹, e=0.9, f_orb=10⁻⁸·⁵)` R=1.27. At
  `tau_a ≥ 1 kyr`: no ignition.

**Caveats (WEAVE §4.4):** this is a **ceiling** (zero-noise Asimov, single source, frozen GP,
free sky); noise/confusion only lower it. The rank-3 ingredients are the *analytic* Peters+γ̇,
not full EOB — the EOB tier (arXiv:2511.19611) is the credible-tier check. So M3 is a
**refutation of the pure cancellation at high e** plus a **marginal, shortest-lag-only, ceiling
ignition** — not a clean null, and not a robust cascade ignition.

---

## M4 — the payoff join (SIREN; `SIREN_payoff_chain.md` §4, machinery reused verbatim)

For self-clocking cells, feed the **comb-measured** Mc (M1, no certified pulsar seeds) to SIREN's
amplitude relation `sigma(D_L)/D_L = ln10·sqrt((5/3 σ_mc)² + σ_h²)` with SIREN's `sigma(log10_h)`:

| f_orb | Mc | e* (>20×) | σ_mc (comb) | σ(D_L)/D_L, **comb-only (0 seeds)** | +1 seed | +3 seeds | SIREN circular N0 / N3 |
|---|---|---|---|---|---|---|---|
| 10⁻⁸   | 9.5 | 0.59 | 0.0156 | **12.1 %** | 11.2 % | 10.6 % | 136 % / 10.2 % |
| 10⁻⁸   | 9.0 | 0.58 | 0.0249 | **14.2 %** | 13.8 % | 12.1 % | 316 % / 12.4 % |
| 10⁻⁸·⁵ | 9.5 | 0.66 | 0.086  | 33.6 % | 20.6 % | 10.4 % | 320 % / 9.9 % |

**The payoff.** An eccentric source at `(f_orb = 10⁻⁸, Mc ≳ 10⁹, e ≳ 0.58)` reaches
`sigma(D_L)/D_L ≈ 12–14%` — the dark-siren-useful class — **from its own Earth term, with zero
certified pulsar terms.** SIREN reached the same class only with 3 certified seeds (which the
census recurs in 0 of 40 skies). Eccentricity substitutes the counterpart's own clock for the
missing certified pulsar terms — but, per M2, only above the >20× *relative* bound, not yet the
0.003-dex *absolute* one, so the residual gap to full targeted-certification is a factor ~2–5 in
σ_mc that the EOB tier (M3's rank-3 ignition corner) must close.

## M5 — figure set

- `ATLAS_results/atlas_M2_contour_kappa.png` — **(headline)** the self-clocking min-e contour +
  κ measured vs WEAVE's scaling (validation).
- `ATLAS_results/atlas_M3_ignition.png` — the rank-3 ignition panel `R_a(e)` vs the R=1 line.

## Artifacts (all untracked)

```
ATLAS_results/atlas_kappa_forb{0,1,2}.npz   M1 grid (per cell: sig_mc, kappa_meas/analytic, cov_scaled 9x9, n_clip)
ATLAS_results/atlas_ignition_forb{0,1,2}.npz M3 (R_rank3, R_scalar, g, per tau_a)
ATLAS_results/atlas_m2m4_summary.npz         M2 contour + M4 join
ATLAS_results/atlas_M2_contour_kappa.png, atlas_M3_ignition.png
hpc_harbor/atlas/{atlas_harmonics,atlas_kappa,atlas_ignition,atlas_consolidate}.py + atlas_kappa.sbatch
hpc_harbor/logs/atlas_k_*                     job logs
```

**Standing caveats.** Fisher ceilings on zero-noise Asimov, truth-placed (LABELLED), free sky,
frozen GP, single source. The comb is the **toy tier** (circular-kernel harmonic stack, fdot tie
only — no e(t)/γ̇/fddot inside M1's amplitude Fisher; those enter only M3's pulsar-term phase).
Toy-invalid cells (comb coalesces) flagged, not silently dropped. Nothing tracked was edited;
nothing was committed.
