"""PRONG2-N: the N_psr dependence of the confusion-knee offset.

Stage A measured the confusion knee at ~52*N* (N* = T*Delta_f) for ONE array size,
N_psr = 116.  P2-C (stagec_p2c.py) fitted a power law over N_psr = 15/30/60/116/200 and
got knee/N* = 0.399 * N_psr^1.026 -- but it never reached an array size at which the two
candidate COUNTING LAWS separate, and it reported no genuine seed scatter (its
`bootstrap_knee_scatter` re-ran the *same* seed batch every time, so it always returned
~0).  This module supplies both.

The counting laws
-----------------
Per 1/T frequency bin an array of N_psr pulsars can resolve

    M_max = (d - m) * N_psr / p

sources, where p is the number of per-source nuisance parameters the fit must absorb.
Multiplying by the number of bins N* = T*Delta_f gives the predicted knee, so the
"offset"  knee/N*  is predicted to be

    (a) N_p / 2      = 0.5000 * N_p   -- p=2, d=2, m=1: THIS runner's parameterisation
                                         (phase0 + log10h are the only per-source
                                         nuisances marginalised in fisher_blocks)
    (b) 3 * N_p / 7  = 0.4286 * N_p   -- the full 7-parameter CW source
    (c) 2 * N_p / 5  = 0.4000 * N_p   -- the coefficient ABACUS_boyle_pen.md banks and
                                         STORY S2.1.3 asserts: Boyle & Pen's 2/7 with
                                         their 2 source-sky parameters removed (7 -> 5)
    (d) 2 * N_p / 7  = 0.2857 * N_p   -- raw Boyle & Pen (quoted in PAPER_figures/p06)

All four are LINEAR in N_p, so the log-log slope is 1 under every one of them; the slope
tests the counting law's FORM, the coefficient tests WHICH law.  At N_p = 116 the (a)/(b)
predictions are 58.0 and 49.7 -- straddling the measured 52.75, hence unresolved.  At
N_p = 232 they are 116.0 and 99.4, a 16.7% split, and (c) predicts 92.8.

Note on (c): ABACUS identifies "our coefficient 0.40 = 2/5" with P2-C's power-law
AMPLITUDE A = 0.3990.  A is only the coefficient if the exponent is exactly 1; P2-C fitted
b = 1.026, and the LOCAL coefficient at N_psr=116 is 52.75/116 = 0.4547, not 0.40.  So the
232 point tests the 2/5 identification directly: under (c) the offset must be 92.8, under
"local coefficient stays ~0.455" it must be ~106.  Whether b is really > 1 is the same
question as whether the counting law's linear FORM holds.

What this runs
--------------
Production Stage-A config, unchanged: T = 15 yr, band 3-20 nHz, 100 ns TOAs, 7 d cadence,
fixed_persource, isotropic sky, rcond = 1e-10, 30 source seeds, at

    N_psr in {30, 58, 116, 232}      (58 = 116/2, 232 = 116*2)

For each N_psr we build PER-SEED marg/cond curves (not just the seed-median curve), so the
knee carries a real seed scatter and the fitted slope carries a real error bar.

Implementation notes
--------------------
* CHUNKED Fisher assembly.  prong2_transition.fisher_blocks vmaps over the whole array at
  once, materialising H = (N_psr, n_toa, 2*n_max); at N_psr=232, n_max=1486 that is 4.4 GB
  before intermediates, which is what OOM'd P2-C's 200-pulsar run under a co-tenant.  We
  accumulate F_ss = sum_a H_a^T H_a over pulsar chunks instead: identical arithmetic
  (gated below to 1e-14), peak 3.1 GiB.
* PER-CONFIG n_max.  Padding makes the result exactly n_max-independent (padded sources are
  masked to zero and their null directions are projected out by the eigh pseudo-inverse),
  so we size n_max to ~1.75x the expected knee per N_psr rather than paying N_psr=232's
  cost at N_psr=30.  Gated directly (G2) rather than assumed.
* TWO-STAGE N grid.  A coarse log grid locates the knee, then a fine log grid (~7% spacing)
  around it resolves the crossing.  The 16.7% (a)/(b) split needs the interpolation
  systematic well under one grid step, which P2-C's single 34%-spaced grid did not give.

Run:  ~/miniforge3/envs/jug-gpu/bin/python stagec_p2n_arraylaw.py
"""

import os

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys
import time

import numpy as np

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import prong2_transition as p2

import jax
import jax.numpy as jnp

YR = p2.YR
REPO = "/home/mattm/projects/HSYMT/CW_transition"
OUT = f"{REPO}/stagec_p2n_results.npz"
FIG = f"{REPO}/p2n_array_law.png"

# ---- production Stage-A config (identical to P2-C / the T15_b3-20 diagnostic row) ----
T = 15 * YR
F_BAND = (3e-9, 2.0e-8)
SIGMA = 1e-7                 # 100 ns
CADENCE_D = 7.0
MODE = "fixed_persource"
LOG10_H = -14.3
RCOND = 1e-10
N_SEEDS = 30
SEED0 = 10                   # same RNG convention as run_sweep
SKY_SEED = 777               # build_array default; arrays are nested (prefix) in N_psr

N_PSR_LIST = [30, 58, 116, 232]

# P2-C anchor used only to SIZE the grids (never to fit); 52.75/116 = 0.4547
COEF_ANCHOR = 52.75 / 116.0
P2C_K116 = 52.75

# candidate counting-law coefficients for knee/N* = c * N_psr
LAWS = {"N_p/2   (p=2, d-m=1)": 0.5,
        "3N_p/7  (7-param)": 3.0 / 7.0,
        "2N_p/5  (ABACUS/BP-sky-fixed)": 2.0 / 5.0,
        "2N_p/7  (Boyle&Pen raw)": 2.0 / 7.0}


# ----------------------------------------------------------------------
# Chunked Fisher assembly (arithmetically identical to p2.fisher_blocks)
# ----------------------------------------------------------------------
@jax.jit
def _chunk_blocks(toas, p_hats_c, L_c, srcs, mask, sigma):
    g, H = jax.vmap(lambda ph, L: p2.analytic_psr_grads(toas, ph, L, srcs, mask))(
        p_hats_c, L_c)
    inv = 1.0 / sigma ** 2
    return (jnp.einsum("pt,pt->p", g, g) * inv,
            jnp.einsum("pt,pta->pa", g, H) * inv,
            jnp.einsum("pta,ptb->ab", H, H) * inv)


def fisher_blocks_chunked(toas, p_hats, L_arr, srcs, mask, sigma, chunk):
    """cond, F_Ls, F_ss with F_ss accumulated over pulsar chunks."""
    cs, fls, F_ss = [], [], None
    for i in range(0, p_hats.shape[0], chunk):
        c, fl, fs = _chunk_blocks(toas, p_hats[i:i + chunk], L_arr[i:i + chunk],
                                  srcs, mask, sigma)
        cs.append(c)
        fls.append(fl)
        F_ss = fs if F_ss is None else F_ss + fs
    return jnp.concatenate(cs), jnp.concatenate(fls), F_ss


def _chunk_for(n_psr, target=32):
    """Largest divisor of n_psr that is <= target (keeps a single compiled chunk shape)."""
    for c in range(min(target, n_psr), 0, -1):
        if n_psr % c == 0:
            return c
    return 1


# ----------------------------------------------------------------------
# Per-seed sweep
# ----------------------------------------------------------------------
def per_seed_ratio(n_values, n_psr, n_max, p_hats, toas, chunk, n_seeds=N_SEEDS,
                   seed0=SEED0, sky_tag=""):
    """marg/cond for every (N, seed).  Returns (len(n_values), n_seeds).

    Per seed the ratio is median_over_pulsars(marg) / median_over_pulsars(cond) --
    the same ratio-of-medians run_sweep forms, just kept un-collapsed over seeds.
    """
    L_arr = jnp.ones(n_psr)
    R = np.zeros((len(n_values), n_seeds))
    for i, n_src in enumerate(n_values):
        t0 = time.time()
        for s in range(n_seeds):
            rng = np.random.default_rng(seed0 + 1000 * s + n_src)
            srcs = p2.draw_sources(n_src, rng, F_BAND, LOG10_H, mode=MODE)
            pad, mask = p2.pad_sources(srcs, n_max)
            blocks = fisher_blocks_chunked(toas, p_hats, L_arr, jnp.array(pad),
                                           jnp.array(mask), SIGMA, chunk)
            marg = p2.marg_from_blocks(*blocks, RCOND)
            R[i, s] = float(jnp.median(marg)) / float(jnp.median(blocks[0]))
        print(f"    [{n_psr}{sky_tag}] N={n_src:5d}  R(med over seeds)={np.median(R[i]):.3f} "
              f"[{np.percentile(R[i],16):.3f},{np.percentile(R[i],84):.3f}]  "
              f"({time.time()-t0:.1f}s)", flush=True)
    return R


def knee_localfit(n, R, level=0.5, lo=0.25, hi=0.75):
    """Knee from a local linear fit of R vs log N over the transition band.

    p2.crossing_N takes the FIRST grid point that drops below `level`, so on a fine grid a
    single noisy point can move it by a whole step.  Regressing R on log N across all
    points with lo < R < hi and solving for R = level averages that noise down, which is
    what the 16.7% coefficient split at N_psr=232 needs.  Falls back to crossing_N if the
    transition band is too thinly sampled.
    """
    n = np.asarray(n, float)
    R = np.asarray(R, float)
    m = (R > lo) & (R < hi) & np.isfinite(R)
    if m.sum() < 3:
        return p2.crossing_N(n, R, level)
    slope, icpt = np.polyfit(np.log(n[m]), R[m], 1)
    if slope >= 0:
        return np.nan
    return float(np.exp((level - icpt) / slope))


def knees_from(n_values, R):
    """Seed-median-curve knee (P2-C convention) + one knee per seed."""
    med = np.median(R, axis=1)
    knee_med = p2.crossing_N(n_values, med, 0.5)
    per_seed = np.array([p2.crossing_N(n_values, R[:, s], 0.5) for s in range(R.shape[1])])
    # non-monotone diagnostic: does any seed curve climb back above 0.5 after crossing?
    recross = 0
    for s in range(R.shape[1]):
        below = np.where(R[:, s] < 0.5)[0]
        if len(below) and np.any(R[below[0]:, s] > 0.5):
            recross += 1
    per_seed_lf = np.array([knee_localfit(n_values, R[:, s]) for s in range(R.shape[1])])
    return knee_med, per_seed, per_seed_lf, recross


def run_one(n_psr, p_hats=None, sky_tag="", n_seeds=N_SEEDS, n_fine=18, knee_hint=None):
    """Two-stage (coarse locate -> fine refine) knee measurement at one array size.

    knee_hint skips the coarse stage (used by the sky-realisation check, where the main
    run has already located the knee for this N_psr).
    """
    N_star = p2.confusion_onset(T, F_BAND)
    n_hi = int(np.ceil(1.75 * COEF_ANCHOR * N_star * n_psr))
    chunk = _chunk_for(n_psr)
    toas = p2.build_toas(T, CADENCE_D)
    if p_hats is None:
        p_hats = p2.build_array(n_psr, SKY_SEED)
    print(f"\n=== N_psr={n_psr}{sky_tag}  n_max={n_hi}  chunk={chunk}  "
          f"n_seeds={n_seeds} ===", flush=True)

    if knee_hint is None:
        coarse = np.unique(np.round(np.logspace(0, np.log10(n_hi), 14)).astype(int))
        R_c = per_seed_ratio(coarse, n_psr, n_hi, p_hats, toas, chunk, n_seeds,
                             sky_tag=sky_tag)
        knee_c = p2.crossing_N(coarse, np.median(R_c, axis=1), 0.5)
        print(f"    -> coarse knee = {knee_c:.1f}", flush=True)
    else:
        coarse = np.array([], dtype=int)
        R_c = np.zeros((0, n_seeds))
        knee_c = float(knee_hint)

    fine = np.unique(np.round(np.logspace(np.log10(max(2.0, 0.60 * knee_c)),
                                          np.log10(min(float(n_hi), 1.70 * knee_c)),
                                          n_fine)).astype(int))
    fine = np.array([n for n in fine if n not in set(coarse.tolist())])
    R_f = per_seed_ratio(fine, n_psr, n_hi, p_hats, toas, chunk, n_seeds, sky_tag=sky_tag)

    n_all = np.concatenate([coarse, fine])
    order = np.argsort(n_all)
    n_all = n_all[order]
    R_all = np.vstack([R_c, R_f])[order]

    knee_med, per_seed, per_seed_lf, recross = knees_from(n_all, R_all)
    knee_lf = knee_localfit(n_all, np.median(R_all, axis=1))
    good = np.isfinite(per_seed)
    good_lf = np.isfinite(per_seed_lf)
    step = float(np.median(np.diff(np.log(fine)))) if len(fine) > 1 else np.nan
    print(f"    -> knee(median curve, first-crossing) = {knee_med:7.1f}   "
          f"knee/N* = {knee_med/N_star:7.2f}\n"
          f"    -> knee(median curve, local fit)      = {knee_lf:7.1f}   "
          f"knee/N* = {knee_lf/N_star:7.2f}\n"
          f"    -> per-seed local fit = {np.median(per_seed_lf[good_lf]):7.1f} "
          f"[16,84] = [{np.percentile(per_seed_lf[good_lf],16):.1f},"
          f"{np.percentile(per_seed_lf[good_lf],84):.1f}]  "
          f"({good_lf.sum()}/{n_seeds} ok; first-crossing {good.sum()}/{n_seeds}, "
          f"{recross} re-crossed)\n"
          f"    -> fine-grid log step = {step:.4f} ({100*(np.exp(step)-1):.1f}%)",
          flush=True)
    return dict(n_psr=n_psr, n_values=n_all, R=R_all, n_max=n_hi, N_star=N_star,
                knee_med=knee_med, knee_lf=knee_lf, knee_seeds=per_seed,
                knee_seeds_lf=per_seed_lf, recross=recross, fine_log_step=step)


# ----------------------------------------------------------------------
# Gates
# ----------------------------------------------------------------------
def gate_chunking(toas):
    """G1: chunked assembly == monolithic p2.fisher_blocks."""
    n_psr, n_max, n_src = 58, 200, 120
    p_hats = p2.build_array(n_psr, SKY_SEED)
    L_arr = jnp.ones(n_psr)
    srcs = p2.draw_sources(n_src, np.random.default_rng(3), F_BAND, LOG10_H, mode=MODE)
    pad, mask = p2.pad_sources(srcs, n_max)
    pad, mask = jnp.array(pad), jnp.array(mask)
    a = fisher_blocks_chunked(toas, p_hats, L_arr, pad, mask, SIGMA, _chunk_for(n_psr))
    b = p2.fisher_blocks(toas, p_hats, L_arr, pad, mask, n_psr, SIGMA)
    rels = [float(jnp.max(jnp.abs(x - y)) / jnp.max(jnp.abs(y))) for x, y in zip(a, b)]
    ma = p2.marg_from_blocks(*a, RCOND)
    mb = p2.marg_from_blocks(*b, RCOND)
    r_marg = float(jnp.max(jnp.abs(ma - mb)) / jnp.max(jnp.abs(mb)))
    print(f"[G1 chunk-vs-monolithic] cond/F_Ls/F_ss rel = "
          f"{rels[0]:.1e}/{rels[1]:.1e}/{rels[2]:.1e}, marg rel = {r_marg:.1e}", flush=True)
    assert max(rels + [r_marg]) < 1e-10, "chunked assembly != monolithic"
    return max(rels + [r_marg])


def gate_nmax(toas):
    """G2: marg/cond is independent of the padding size n_max (per-config n_max is safe)."""
    n_psr = 116
    p_hats = p2.build_array(n_psr, SKY_SEED)
    L_arr = jnp.ones(n_psr)
    worst = 0.0
    for n_src in (100, 400, 700):
        rng_state = np.random.default_rng(SEED0 + n_src)
        srcs = p2.draw_sources(n_src, rng_state, F_BAND, LOG10_H, mode=MODE)
        vals = []
        for n_max in (743, 1486):
            pad, mask = p2.pad_sources(srcs, n_max)
            bl = fisher_blocks_chunked(toas, p_hats, L_arr, jnp.array(pad),
                                       jnp.array(mask), SIGMA, _chunk_for(n_psr))
            marg = p2.marg_from_blocks(*bl, RCOND)
            vals.append(float(jnp.median(marg)) / float(jnp.median(bl[0])))
        rel = abs(vals[1] - vals[0]) / vals[0]
        worst = max(worst, rel)
        print(f"[G2 n_max independence] N={n_src:4d}: R(743)={vals[0]:.6f} "
              f"R(1486)={vals[1]:.6f}  rel={rel:.2e}", flush=True)
    assert worst < 1e-6, f"marg/cond depends on n_max ({worst:.2e})"
    return worst


# ----------------------------------------------------------------------
def main():
    t_start = time.time()
    toas = p2.build_toas(T, CADENCE_D)
    N_star = p2.confusion_onset(T, F_BAND)
    print("=== PRONG2-N: knee offset vs N_psr ===", flush=True)
    print(f"config: T=15 yr, band {F_BAND[0]*1e9:.0f}-{F_BAND[1]*1e9:.0f} nHz, "
          f"sigma={SIGMA*1e9:.0f} ns, cadence {CADENCE_D:.0f} d, n_toa={len(toas)}, "
          f"mode={MODE}, rcond={RCOND:.0e}, n_seeds={N_SEEDS}", flush=True)
    print(f"N* = T*Delta_f = {N_star:.4f}\n", flush=True)

    g1 = gate_chunking(toas)
    g2 = gate_nmax(toas)
    print("", flush=True)

    results = [run_one(n_psr) for n_psr in N_PSR_LIST]

    npsr = np.array([r["n_psr"] for r in results], float)
    knee_xing = np.array([r["knee_med"] for r in results])   # P2-C convention
    knee_med = np.array([r["knee_lf"] for r in results])     # local fit (primary)
    off_xing = knee_xing / N_star
    off_med = knee_med / N_star                              # knee/N* (primary)
    off_seed = np.array([r["knee_seeds_lf"] for r in results], float) / N_star

    off_seed_lo = np.nanpercentile(off_seed, 16, axis=1)
    off_seed_hi = np.nanpercentile(off_seed, 84, axis=1)
    off_seed_med = np.nanmedian(off_seed, axis=1)
    n_ok = np.isfinite(off_seed).sum(axis=1)
    # SE of the offset: sigma_log / sqrt(n)   (log space; the fit is log-log)
    sig_log = np.nanstd(np.log(off_seed), axis=1)
    se_log = sig_log / np.sqrt(n_ok)

    # ---- GATE G3: N_psr=116 reproduces P2-C / Stage A ----
    i116 = N_PSR_LIST.index(116)
    k116 = off_xing[i116]
    dev116 = abs(k116 - P2C_K116) / P2C_K116
    tol = max(3 * se_log[i116], 0.05)
    print(f"\n[G3 116-reproduction] knee/N* = {k116:.2f} (first-crossing, P2-C convention) "
          f"vs P2-C {P2C_K116:.2f} (Stage A ~52): rel dev {100*dev116:.1f}%, "
          f"tol {100*tol:.1f}% -> {'PASS' if dev116 < tol else 'CHECK'}", flush=True)

    # ---- power-law fit knee/N* = A * N_psr^b, weighted by the per-N_psr log SE ----
    x = np.log(npsr)
    y = np.log(off_med)
    w = 1.0 / se_log ** 2
    W = np.diag(w)
    X = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ (X.T @ W @ y)
    logA, b = beta
    A = float(np.exp(logA))
    sb = float(np.sqrt(cov[1, 1]))
    chi2 = float((y - X @ beta) @ W @ (y - X @ beta))
    dof = len(x) - 2

    print(f"\n================ PRONG2-N SUMMARY ================", flush=True)
    print(f"{'N_psr':>6s} {'knee':>8s} {'knee/N*':>9s} {'[16,84] seeds':>20s} "
          f"{'c=off/N_psr':>12s} {'off(1st-x)':>11s} {'ok':>7s}", flush=True)
    for i, r in enumerate(results):
        print(f"{int(npsr[i]):>6d} {knee_med[i]:>8.1f} {off_med[i]:>9.2f} "
              f"{'[%.2f, %.2f]' % (off_seed_lo[i], off_seed_hi[i]):>20s} "
              f"{off_med[i]/npsr[i]:>12.4f} {off_xing[i]:>11.2f} "
              f"{n_ok[i]:>3d}/{N_SEEDS}", flush=True)

    print(f"\npower law: knee/N* = {A:.4f} * N_psr^{b:.4f}  (slope err {sb:.4f}; "
          f"chi2/dof = {chi2:.2f}/{dof})", flush=True)
    z_slope = (b - 1.0) / sb
    print(f"  counting-law prediction is slope = 1 (all three laws are linear): "
          f"measured {b:.3f} +- {sb:.3f}  -> {z_slope:+.1f} sigma "
          f"({'CONSISTENT' if abs(z_slope) < 3 else 'INCONSISTENT'} with linear)",
          flush=True)

    # ---- coefficient discrimination at N_psr = 232 ----
    i232 = N_PSR_LIST.index(232)
    off232 = off_med[i232]
    se232 = off232 * se_log[i232]                    # (1) source-seed SE
    est232 = abs(off_med[i232] - off_xing[i232])     # (2) knee-estimator systematic
    print(f"\n---- coefficient test at N_psr = 232 ----", flush=True)
    print(f"measured knee/N* = {off232:.2f}   error budget: seed SE {se232:.2f}, "
          f"estimator (local fit vs first crossing) {est232:.2f}; "
          f"per-seed 16/84 = [{off_seed_lo[i232]:.1f}, {off_seed_hi[i232]:.1f}]", flush=True)
    print("  (the per-seed spread is small because all seeds share ONE pulsar array; the "
          "sky-realisation term below is the dominant statistical error)", flush=True)

    # ---- sky-realisation systematic at N_psr = 232 (single sky per size otherwise) ----
    print(f"\n---- sky-realisation check at N_psr = 232 (2 extra sky seeds, 10 source "
          f"seeds) ----", flush=True)
    sky_offsets = [off232]
    for sky in (778, 779):
        ph = p2.build_array(232, sky)
        rr = run_one(232, p_hats=ph, sky_tag=f" sky{sky}", n_seeds=10, n_fine=12,
                     knee_hint=knee_med[i232])
        sky_offsets.append(rr["knee_lf"] / N_star)
    sky_offsets = np.array(sky_offsets)
    sky_spread = float(np.std(sky_offsets, ddof=1))
    print(f"[sky] knee/N* over sky seeds {SKY_SEED},778,779 = "
          f"{np.round(sky_offsets, 2).tolist()}  -> spread (sd) {sky_spread:.2f} "
          f"({100*sky_spread/off232:.1f}%)", flush=True)

    tot232 = float(np.sqrt(se232 ** 2 + est232 ** 2 + sky_spread ** 2))
    print(f"\n---- N_psr = 232 VERDICT ----", flush=True)
    print(f"knee/N* = {off232:.2f} +- {tot232:.2f}  (seed {se232:.2f} + estimator "
          f"{est232:.2f} + sky {sky_spread:.2f}, in quadrature)", flush=True)
    zs = {}
    for name, c in LAWS.items():
        pred = c * 232.0
        z = (off232 - pred) / tot232
        zs[name] = z
        print(f"  {name:>22s}: predicts {pred:7.2f}  -> measured is "
              f"{100*(off232-pred)/pred:+6.1f}% ({z:+.1f} sigma)", flush=True)
    best = min(zs, key=lambda k: abs(zs[k]))
    two = ("N_p/2   (p=2, d-m=1)", "3N_p/7  (7-param)")
    if abs(zs[best]) < 3:
        others = [n for n in LAWS if n != best and abs(zs[n]) >= 3]
        verdict = (f"SUPPORTS {best} (|z|={abs(zs[best]):.1f})"
                   + (f"; rejects {', '.join(others)}" if others else
                      "; but does NOT reject the alternatives"))
        if all(abs(zs[n]) < 3 for n in two):
            verdict = (f"the two coefficients still do NOT separate "
                       f"(|z| = {abs(zs[two[0]]):.1f} and {abs(zs[two[1]]):.1f}); "
                       f"measured c = {off232/232:.4f}")
    else:
        verdict = (f"REJECTS ALL THREE at >3 sigma (closest {best}, "
                   f"|z|={abs(zs[best]):.1f}); measured c = {off232/232:.4f} "
                   f"lies between 3/7={3/7:.4f} and 1/2=0.5000")
    print(f"  VERDICT: {verdict}", flush=True)

    np.savez(OUT,
             n_psr=npsr, knee=knee_med, knee_over_Nstar=off_med, N_star=N_star,
             knee_firstcrossing=knee_xing, knee_over_Nstar_firstcrossing=off_xing,
             knee_seeds=np.array([r["knee_seeds"] for r in results]),
             knee_seeds_localfit=np.array([r["knee_seeds_lf"] for r in results]),
             offset_seed_median=off_seed_med, offset_seed_lo=off_seed_lo,
             offset_seed_hi=off_seed_hi, offset_se_log=se_log, n_seeds_ok=n_ok,
             A=A, b=b, b_err=sb, chi2=chi2, dof=dof,
             coef=off_med / npsr, coef_232=off232 / 232.0,
             off232=off232, se232=se232, sky_offsets=sky_offsets,
             sky_spread=sky_spread, tot_err_232=tot232, est_err_232=est232,
             law_names=np.array(list(LAWS.keys())),
             law_coefs=np.array(list(LAWS.values())),
             law_z=np.array([zs[k] for k in LAWS]),
             verdict=np.array(verdict),
             gate_chunk=g1, gate_nmax=g2, gate_116=k116, p2c_k116=P2C_K116,
             n_max_used=np.array([r["n_max"] for r in results]),
             fine_log_step=np.array([r["fine_log_step"] for r in results]),
             recross=np.array([r["recross"] for r in results]),
             config=np.array(f"T=15yr band={F_BAND} sigma={SIGMA} cadence={CADENCE_D} "
                             f"mode={MODE} rcond={RCOND} n_seeds={N_SEEDS}"),
             **{f"R_{int(npsr[i])}": results[i]["R"] for i in range(len(results))},
             **{f"n_{int(npsr[i])}": results[i]["n_values"] for i in range(len(results))})
    print(f"\nsaved {OUT}", flush=True)

    _figure(npsr, off_med, off_seed_lo, off_seed_hi, A, b, sb, results, N_star, off232,
            tot232)
    print(f"[timing] TOTAL {time.time()-t_start:.1f}s", flush=True)
    return 0


def _figure(npsr, off, lo, hi, A, b, sb, results, N_star, off232, tot232):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.0))

    a = ax[0]
    for r in results:
        med = np.median(r["R"], axis=1)
        a.plot(r["n_values"], med, "-o", ms=3, lw=1.3, label=f"$N_p$={r['n_psr']}")
        a.axvline(r["knee_lf"], ls=":", lw=1, alpha=.5,
                  color=a.lines[-1].get_color())
    a.axhline(0.5, color="k", lw=.7, ls="--", alpha=.6)
    a.set_xscale("log")
    a.set_xlabel("number of CW sources $N$")
    a.set_ylabel("marginal / conditional information")
    a.set_ylim(0, 1.02)
    a.set_title("Confusion knee moves right with array size", fontsize=10)
    a.grid(alpha=.25, which="both")
    a.legend(fontsize=8)

    a = ax[1]
    err = np.vstack([off - lo, hi - off])
    a.errorbar(npsr, off, yerr=err, fmt="o", ms=8, color="#2c7fb8", capsize=3,
               label="measured knee/$N^*$ (16/84 over seeds)")
    xs = np.logspace(np.log10(25), np.log10(280), 100)
    a.plot(xs, A * xs ** b, "-", color="#d95f0e", lw=1.6,
           label=f"fit ${A:.3f}\\,N_p^{{{b:.3f}\\pm{sb:.3f}}}$")
    for (nm, c), col, ls in zip(LAWS.items(),
                                ("#333333", "#7b3294", "#c51b7d", "#4d9221"),
                                ("--", "-.", (0, (3, 1, 1, 1)), ":")):
        a.plot(xs, c * xs, color=col, ls=ls, lw=1.2, alpha=.85, label=nm.strip())
    a.errorbar([232], [off232], yerr=[[tot232], [tot232]],
               fmt="none", ecolor="k", capsize=5, lw=1.8, zorder=5,
               label="232: full error budget")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("$N_{\\rm psr}$")
    a.set_ylabel("knee / $N^*$   (array-resolution offset)")
    a.set_xticks(N_PSR_LIST)
    a.set_xticklabels(N_PSR_LIST)
    a.set_title("Counting-law test: slope and coefficient", fontsize=10)
    a.grid(alpha=.25, which="both")
    a.legend(fontsize=8, loc="upper left")

    fig.suptitle(f"PRONG2-N: confusion-knee offset vs array size "
                 f"(T=15 yr, 3-20 nHz, $N^*$={N_star:.2f}, {N_SEEDS} seeds)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130, bbox_inches="tight")
    print(f"saved {FIG}", flush=True)


def replot():
    """Rebuild the figure (and the chi2-inflated slope error) from the saved npz."""
    d = np.load(OUT, allow_pickle=True)
    npsr = d["n_psr"]
    N_star = float(d["N_star"])
    off = d["knee_over_Nstar"]
    results = [dict(n_psr=int(npsr[i]), R=d[f"R_{int(npsr[i])}"],
                    n_values=d[f"n_{int(npsr[i])}"], knee_lf=float(d["knee"][i]))
               for i in range(len(npsr))]
    A, b, sb = float(d["A"]), float(d["b"]), float(d["b_err"])
    chi2, dof = float(d["chi2"]), int(d["dof"])
    # the 4 points scatter by more than the source-seed SE alone (chi2/dof > 1), so the
    # naive slope error is optimistic; inflate it by sqrt(chi2/dof) as usual.
    infl = max(1.0, np.sqrt(chi2 / dof))
    sb_inf = sb * infl
    print(f"power law: knee/N* = {A:.4f} * N_psr^{b:.4f}", flush=True)
    print(f"  slope err {sb:.4f} (seed-SE only) -> {sb_inf:.4f} after chi2/dof="
          f"{chi2/dof:.2f} inflation; slope-1 = {(b-1)/sb_inf:+.1f} sigma", flush=True)
    # the local coefficient c(N_p) = A*N_p^(b-1) is NOT constant if b != 1
    for c_name, c_val in LAWS.items():
        if b != 1.0:
            n_cross = float(np.exp(np.log(c_val / A) / (b - 1.0)))
            print(f"  c(N_p) passes {c_name.strip()} (c={c_val:.4f}) at N_psr = "
                  f"{n_cross:.3g}", flush=True)
    out = {k: d[k] for k in d.files}
    out["b_err_inflated"] = sb_inf
    out["chi2_over_dof"] = chi2 / dof
    np.savez(OUT, **out)
    _figure(npsr, off, d["offset_seed_lo"], d["offset_seed_hi"], A, b, sb_inf, results,
            N_star, float(d["off232"]), float(d["tot_err_232"]))
    print(f"re-saved {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "plot":
        sys.exit(replot())
    sys.exit(main())
