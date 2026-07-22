"""FORGE-G2 gates: SMASK as a JITTED RUNTIME ARGUMENT (TURBINE's proposal).

WHAT CHANGED (trackB_b1_core / trackB_b1). The FORGE-G feed weight `smask` used to be
closed over inside the jitted evaluators via `_params` at trace time -- every pattern
change invalidated and recompiled all 117 evaluators (116 per-pulsar + 1 batched base;
TURBINE measured 83-800 s per frontier update, 17-32 GPU-hr over the GLACIER campaign).
It is now a runtime argument of every jitted evaluator: a pattern change is a plain
argument change (no cache clear, no recompile); only switching None <-> array retraces
once per direction. The baked-in-constant path is retired. MaskedDelay additionally
gained the TURBINE grad guard (DOUBLE-WHERE): carried-source waveform INPUTS are
sanitised to SMASK_SAFE before evaluation, so autodiff never touches a carried source's
actual (possibly non-evaluable) params.

GATES (this script; Bg5a/b/c re-run separately via trackB_b1_softsource.py gate):
  SG1  pattern bit-exactness vs the BAKED incumbent: every (pattern, theta) value equals
       the pre-edit capture (b1_smask_baked_ref.npz) at max|diff| = 0.0. This is
       TURBINE's bit-exactness check AND certifies the double-where primal is inert.
  SG2  mechanism identity: a freshly-jitted closure that CAPTURES the pattern as a
       constant (the incumbent bake mechanism) == the runtime-arg call, 0.0.
  SG3  warm flip: repeated pattern flips on compiled evaluators, median < 10 ms,
       jit cache sizes constant (no recompile).
  SG4  default-inert: smask None -> bit-identical to the incumbent B1Marg path
       (fast-full default AND masked path), 0.0 vs the pre-edit capture.
  SG5  gradient at a carried->fed boundary: with a carried source at NON-EVALUABLE
       params, grad(theta) and grad(smask) are finite, the carried block's grad is
       exactly 0, and the fed-block grads are invariant to the carried junk; the
       single-where NEGATIVE CONTROL reproduces the NaN the guard exists to stop.
  SG6  full-stack B1Marg: lnmarg/posterior at loud-only == pre-edit capture (0.0),
       junk-reservoir invariance (Bg5c) preserved, pattern flips leave the batched
       evaluator caches untouched.

Env: jug-gpu, XLA autotune off, x64. Matt commits; never commit from here.
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time
import numpy as np

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import trackB_b1_core as C
import trackB_b1 as B1
import trackB_estimator as TE
import jax
import jax.numpy as jnp

CWT = "/home/mattm/projects/HSYMT/CW_transition"
REF = f"{CWT}/b1_smask_baked_ref.npz"
NCW = 16
N_LOUD = C.N_LOUD
FLIP_MS = 10.0                        # gate 3 budget per warm pattern flip


def cache_size(fn):
    try:
        return int(fn._cache_size())
    except Exception:
        return -1                     # API absent -> timing-only evidence


def patterns():
    r01 = np.random.default_rng(42).integers(0, 2, NCW).astype(float)
    r01[0] = 1.0; r01[NCW - 1] = 0.0
    return dict(
        none=None,
        ones=np.ones(NCW),
        loud=np.array([1.0] * N_LOUD + [0.0] * (NCW - N_LOUD)),
        zeros=np.zeros(NCW),
        rand01=r01,
        frac=np.array([1.0, 0.6, 0.3, 0.5] + [0.0] * (NCW - 4)),
    )


def main():
    print("=== FORGE-G2 RUNTIME-SMASK GATES ===", flush=True)
    ok = True
    have_ref = os.path.exists(REF)
    ref = dict(np.load(REF)) if have_ref else None
    if not have_ref:
        print(f"  !! no baked reference at {REF} -- SG1/SG4/SG6 vs-incumbent checks SKIPPED "
              f"(mechanism-identity SG2 still certifies the runtime path)", flush=True)

    P = C.build_b1_problem()
    sp = P["sp"]; amo = P["amo"]; tt = P["theta_truth"]; nd = P["n_dist"]
    npsr = P["npsr"]
    ones_mask = jnp.ones(npsr)
    data0 = amo["inject_delay"](jnp.asarray(tt), ones_mask)

    # the capture's exact deterministic theta set
    rng = np.random.default_rng(0)
    scale = TE.phi_scale(P); plo, phi_hi = TE.phi_bounds(P)
    tt_pert = tt.copy()
    tt_pert[nd:] = np.clip(tt_pert[nd:] + rng.normal(size=len(tt) - nd) * scale * 0.05,
                           plo, phi_hi)
    tt_junk = tt.copy()
    src_j = tt_junk[nd:].reshape(NCW, 8)
    src_j[N_LOUD:] = np.random.default_rng(31337).normal(size=(NCW - N_LOUD, 8))
    src_j[N_LOUD, 0] = 1.5
    tt_junk[nd:] = src_j.reshape(-1)
    THETAS = dict(tt=tt, pert=tt_pert, junk=tt_junk)
    PATS = patterns()

    # ---- SG1: pattern bit-exactness vs the baked incumbent capture ----
    print("\n-- SG1: runtime-smask == BAKED incumbent, every pattern, max|diff| = 0.0 --", flush=True)
    w1 = 0.0
    if have_ref:
        for pk, w in PATS.items():
            sp.set_smask(w)
            a, b = sp._ppab(jnp.asarray(tt), data0, ones_mask, sp.smask)
            da = float(np.max(np.abs(np.asarray(a) - ref[f"{pk}_a"])))
            db = float(np.max(np.abs(np.asarray(b) - ref[f"{pk}_b"])))
            w1 = max(w1, da, db)
            for thk, th in THETAS.items():
                if f"{pk}_lnl_{thk}" not in ref:
                    continue
                d = abs(sp.lnL(th, data0, ones_mask) - float(ref[f"{pk}_lnl_{thk}"]))
                w1 = max(w1, d)
            print(f"   [{pk:6s}] max|d a|={da:.3e}  max|d b|={db:.3e}", flush=True)
        sp.set_smask(None)
        g1 = w1 == 0.0; ok &= g1
        print(f"   [SG1] max|diff| over patterns x thetas = {w1:.3e} -> "
              f"{'PASS' if g1 else 'FAIL'} (== 0.0)", flush=True)
    else:
        print("   [SG1] SKIPPED (no reference)", flush=True)

    # ---- SG2: mechanism identity -- constant-captured (baked) jit == runtime arg ----
    print("\n-- SG2: freshly-BAKED closure (constant capture) == runtime-arg call, 0.0 --", flush=True)
    w2 = 0.0
    for pk in ("loud", "frac", "rand01"):
        wj = jnp.asarray(PATS[pk])
        baked = jax.jit(lambda t, d, m, _w=wj: sp._per_pulsar_ab_impl(t, d, m, _w))
        a0, b0 = baked(jnp.asarray(tt), data0, ones_mask)
        sp.set_smask(PATS[pk])
        a1, b1 = sp._ppab(jnp.asarray(tt), data0, ones_mask, sp.smask)
        d = max(float(np.max(np.abs(np.asarray(a0) - np.asarray(a1)))),
                float(np.max(np.abs(np.asarray(b0) - np.asarray(b1)))))
        w2 = max(w2, d)
        print(f"   [{pk:6s}] max|baked - runtime| = {d:.3e}", flush=True)
    g2 = w2 == 0.0; ok &= g2
    print(f"   [SG2] -> {'PASS' if g2 else 'FAIL'} (== 0.0)", flush=True)

    # ---- SG3: warm pattern flip < 10 ms, no recompile ----
    print("\n-- SG3: warm pattern flip (compiled evaluators, alternating patterns) --", flush=True)
    p0 = int(P["census_idx"][0])
    sp.AI = np.asarray(P["AI"])
    fn = sp._pulsar_ab_fn(p0)
    tb = jnp.asarray(tt); EVp = jnp.asarray(P["EV_truth"][p0])
    wA, wB = PATS["loud"], PATS["rand01"]
    for w in (wA, wB):                             # warm both patterns once
        sp.set_smask(w)
        jax.block_until_ready(fn(tb, EVp, data0, ones_mask, sp.smask))
        jax.block_until_ready(sp._ppab(tb, data0, ones_mask, sp.smask))
    cs_fn0, cs_pp0 = cache_size(fn), cache_size(sp._ppab)
    # fixed-pattern baseline (no flip)
    base = []
    for _ in range(20):
        t0 = time.perf_counter()
        jax.block_until_ready(fn(tb, EVp, data0, ones_mask, sp.smask))
        base.append(time.perf_counter() - t0)
    # alternating flips: set_smask + evaluator call
    flips = []
    for k in range(20):
        w = wA if k % 2 == 0 else wB
        t0 = time.perf_counter()
        sp.set_smask(w)
        jax.block_until_ready(fn(tb, EVp, data0, ones_mask, sp.smask))
        flips.append(time.perf_counter() - t0)
    # full 116-pulsar base evaluator under alternating patterns (reported, not gated on 10 ms)
    ppab_flip = []
    for k in range(6):
        sp.set_smask(wA if k % 2 == 0 else wB)
        t0 = time.perf_counter()
        jax.block_until_ready(sp._ppab(tb, data0, ones_mask, sp.smask))
        ppab_flip.append(time.perf_counter() - t0)
    cs_fn1, cs_pp1 = cache_size(fn), cache_size(sp._ppab)
    sp.set_smask(None)
    med_b = float(np.median(base) * 1e3); med_f = float(np.median(flips) * 1e3)
    no_recompile = (cs_fn0 == cs_fn1) and (cs_pp0 == cs_pp1)
    g3 = (med_f < FLIP_MS) and (cs_fn1 == -1 or no_recompile)
    ok &= g3
    print(f"   per-pulsar evaluator: warm fixed {med_b:.2f} ms, warm FLIP {med_f:.2f} ms "
          f"(median of 20; gate < {FLIP_MS:.0f} ms)", flush=True)
    print(f"   full-array _ppab under alternating patterns: median {np.median(ppab_flip)*1e3:.0f} ms "
          f"(reference; a baked flip re-COMPILED here)", flush=True)
    print(f"   jit cache sizes across flips: per-pulsar {cs_fn0}->{cs_fn1}, _ppab {cs_pp0}->{cs_pp1} "
          f"({'unchanged -> NO RECOMPILE' if no_recompile else 'CHANGED' if cs_fn1 != -1 else 'API n/a'})",
          flush=True)
    if have_ref:
        inc = ", ".join(f"{pk} {float(ref[f'{pk}_t_flip']):.0f}s" for pk in ("ones", "loud", "frac"))
        print(f"   incumbent BAKED flip cost on this box (pre-edit capture): {inc}", flush=True)
    print(f"   [SG3] -> {'PASS' if g3 else 'FAIL'}", flush=True)

    # ---- SG4: default-inert -- smask None bit-identical to the incumbent ----
    print("\n-- SG4: default-inert (smask None == incumbent, bit-identical) --", flush=True)
    if have_ref:
        w4 = 0.0
        a, b = sp._ppab(jnp.asarray(tt), data0, ones_mask, None)
        w4 = max(w4, float(np.max(np.abs(np.asarray(a) - ref["none_a"]))),
                 float(np.max(np.abs(np.asarray(b) - ref["none_b"]))))
        for thk in ("tt", "pert"):
            w4 = max(w4, abs(sp.lnL(THETAS[thk], data0, ones_mask) - float(ref[f"none_lnl_{thk}"])))
        g4 = w4 == 0.0; ok &= g4
        print(f"   split level: max|diff| = {w4:.3e} -> {'PASS' if g4 else 'FAIL'} (== 0.0)", flush=True)
    else:
        g4 = True
        print("   split level: SKIPPED (no reference)", flush=True)

    # ---- SG5: gradient at a carried->fed boundary (the TURBINE guard) ----
    print("\n-- SG5: grad through where() with a carried NON-EVALUABLE source --", flush=True)
    msd = amo["internals"]["msds"][p0]
    theta_g = tt.copy()
    gsrc = theta_g[nd:].reshape(NCW, 8)
    gsrc[N_LOUD] = [1.5, 0.0, 0.0, 30.0, -2.0, -13.7, 0.0, 0.0]   # arccos domain + absurd mc/fgw
    theta_g[nd:] = gsrc.reshape(-1)
    theta_g2 = theta_g.copy()                                     # different junk, same carried slot
    g2src = theta_g2[nd:].reshape(NCW, 8)
    g2src[N_LOUD] = [2.5, 1.0, 3.0, 45.0, -1.0, -13.7, 1.0, 1.0]
    theta_g2[nd:] = g2src.reshape(-1)
    w_loud = jnp.asarray(PATS["loud"])

    def f_new(theta, w):
        params = sp._params(theta, data0, ones_mask, w)
        return jnp.sum(msd(params) ** 2)

    def f_single_where(theta, w):
        """The PRE-GUARD masked delay (single where, inputs NOT sanitised) -- negative control."""
        params = sp._params(theta, data0, ones_mask, None)
        cw = {k: jnp.stack([params[f"cw_{k}{s}"] for s in msd._suf]) for k in msd._G}
        L = params[msd._pdist_key]; m = params[C.PMASK][msd.ipsr]
        args = (msd.psr.toas, msd.psr.pos, cw["cos_gwtheta"], cw["gwphi"], cw["cos_inc"],
                cw["log10_mc"], cw["log10_fgw"], cw["log10_h"], cw["phase0"], cw["psi"])
        full = msd._full(*args, L); earth = msd._earth(*args)
        gate = (w > 0.0)[:, jnp.newaxis]; wcol = w[:, jnp.newaxis]
        d_full = jnp.sum(jnp.where(gate, wcol * full, 0.0), axis=0)
        d_earth = jnp.sum(jnp.where(gate, wcol * earth, 0.0), axis=0)
        return jnp.sum((d_earth + m * (d_full - d_earth)) ** 2)

    tj = jnp.asarray(theta_g)
    v_new = float(f_new(tj, w_loud)); v_ctl = float(f_single_where(tj, w_loud))
    gt = np.asarray(jax.grad(f_new, 0)(tj, w_loud))
    gw = np.asarray(jax.grad(f_new, 1)(tj, w_loud))
    gt2 = np.asarray(jax.grad(f_new, 0)(jnp.asarray(theta_g2), w_loud))
    gt_ctl = np.asarray(jax.grad(f_single_where, 0)(tj, w_loud))
    carried = np.zeros(len(tt), bool)
    for s in range(NCW):
        if PATS["loud"][s] == 0.0:
            carried[nd + 8 * s: nd + 8 * (s + 1)] = True
    b_fin = bool(np.all(np.isfinite(gt)) and np.all(np.isfinite(gw)))
    b_zero = float(np.max(np.abs(gt[carried]))) == 0.0
    b_inv = float(np.max(np.abs(gt[~carried] - gt2[~carried]))) == 0.0
    b_wcarr = float(np.max(np.abs(gw[np.asarray(PATS['loud']) == 0.0]))) == 0.0
    b_ctl = not bool(np.all(np.isfinite(gt_ctl)))
    b_prim = np.isfinite(v_new) and np.isfinite(v_ctl) and v_new == v_ctl
    print(f"   primal: double-where {v_new:.6e} == single-where {v_ctl:.6e}: {b_prim}", flush=True)
    print(f"   grad(theta) finite: {b_fin}; carried-block grad == 0: {b_zero}; "
          f"fed-block grad invariant to carried junk: {b_inv}", flush=True)
    print(f"   grad(smask) finite + zero on carried (w_s == 0 boundary): "
          f"{bool(np.all(np.isfinite(gw)))} / {b_wcarr}", flush=True)
    print(f"   NEGATIVE CONTROL single-where grad has NaN: {b_ctl} "
          f"(n nonfinite = {int(np.sum(~np.isfinite(gt_ctl)))})", flush=True)
    g5 = b_fin and b_zero and b_inv and b_wcarr and b_ctl and b_prim
    ok &= g5
    print(f"   [SG5] -> {'PASS' if g5 else 'FAIL'}", flush=True)

    # ---- SG6: full-stack B1Marg vs incumbent capture + flip stability ----
    print("\n-- SG6: B1Marg full stack (default-inert, loud-only == capture, flips free) --", flush=True)
    rng2 = np.random.default_rng(7)
    src0 = tt[nd:].copy()
    srcs = np.vstack([src0] + [np.clip(src0 + rng2.normal(size=src0.size) * scale * 0.02,
                                       plo, phi_hi) for _ in range(3)])
    E = B1.B1Marg(P, data0, P["mask_one"])
    w6 = 0.0; g6 = True
    lm_ff = E.lnmarg(srcs)
    if have_ref:
        w6 = max(w6, float(np.max(np.abs(lm_ff - ref["marg_none_lm"]))))
        print(f"   default (fast-full, no smask): max|lnmarg - incumbent| = {w6:.3e}", flush=True)
    E.sp.enable_fast_full(False); E.all_on = False; E._abb = {}
    E.set_smask(None)
    lm_mn = E.lnmarg(srcs)
    if have_ref:
        d = float(np.max(np.abs(lm_mn - ref["marg_none_masked_lm"])))
        w6 = max(w6, d)
        print(f"   masked path, smask None:       max|lnmarg - incumbent| = {d:.3e}", flush=True)
    t0 = time.time()
    E.set_smask(PATS["loud"])
    lm_loud = E.lnmarg(srcs)
    t_first = time.time() - t0
    post_l, LNL_l, _ = E.posterior(src0)
    lm_junk = E.lnmarg(src_j.reshape(-1)[None])
    if have_ref:
        d1 = float(np.max(np.abs(lm_loud - ref["marg_loud_lm"])))
        d2 = float(np.max(np.abs(post_l["q_max"] - ref["marg_loud_qmax"])))
        d3 = float(np.max(np.abs(LNL_l - ref["marg_loud_LNL"])))
        d4 = float(np.max(np.abs(lm_junk - ref["marg_loud_lm_junk"])))
        w6 = max(w6, d1, d2, d3, d4)
        print(f"   loud-only: |d lnmarg|={d1:.3e} |d q_max|={d2:.3e} |d LNL|={d3:.3e} "
              f"|d junk-res|={d4:.3e}  (None->array structure retrace paid once: {t_first:.0f}s)",
              flush=True)
    cs0 = [cache_size(E._abb[p]) for p in sorted(E._abb)][:8]
    t0 = time.time(); E.set_smask(PATS["frac"]); lm_frac = E.lnmarg(srcs); t_frac = time.time() - t0
    t0 = time.time(); E.set_smask(PATS["rand01"]); _ = E.lnmarg(srcs); t_r01 = time.time() - t0
    t0 = time.time(); E.set_smask(PATS["loud"]); lm_loud2 = E.lnmarg(srcs); t_re = time.time() - t0
    cs1 = [cache_size(E._abb[p]) for p in sorted(E._abb)][:8]
    d_re = float(np.max(np.abs(lm_loud2 - lm_loud)))
    b_cache = cs0 == cs1
    b_finite = bool(np.all(np.isfinite(lm_frac)))
    g6 = (w6 == 0.0) and (d_re == 0.0) and b_cache and b_finite
    ok &= g6
    print(f"   pattern flips on WARM evaluators: frac {t_frac:.1f}s, rand01 {t_r01:.1f}s, "
          f"back-to-loud {t_re:.1f}s (vs first-structure {t_first:.0f}s); "
          f"loud repeat |diff| = {d_re:.3e}", flush=True)
    print(f"   batched-evaluator cache sizes unchanged across flips: {b_cache} ({cs0[:4]}...)", flush=True)
    print(f"   [SG6] -> {'PASS' if g6 else 'FAIL'}", flush=True)
    E.set_smask(None)

    print(f"\n=== FORGE-G2 RUNTIME-SMASK GATES: {'ALL PASS' if ok else 'FAIL'} ===", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    sys.exit(main())
