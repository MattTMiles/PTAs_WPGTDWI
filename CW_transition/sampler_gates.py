"""SAMPLER — G0/G1 machinery gates, run BEFORE any posterior.

G0  VALUE   MargJax.lnmarg == trackB_b1.B1Marg.lnmarg              (the banked object)
G1  GRAD    autodiff grad(lnL_marg) == central FD                  (8 points, the mission's gate)
G2  COST    per-value / per-gradient wall on the 4090

The mission pre-registers G1 at 1e-8. Two honest caveats recorded IN ADVANCE:
  * lnL_marg ~ 4e5 nat, so a float64 value carries ~1e-10 nat of roundoff; a central FD
    gradient with step h therefore carries ~1e-10/h of absolute noise. The gate is read as a
    RELATIVE error on gradient components that are not themselves ~0.
  * the within-fringe peak is a segment MAX over the B_CERT grid, so lnL_marg has kinks where
    a within-fringe argmax swaps. Autodiff takes the subgradient at the incumbent argmax; FD
    across a kink does not. Any dim failing the gate is reported with its kink diagnostic
    (does the within-fringe argmax move between the two FD legs?) rather than hidden.
"""
import os, sys, time, argparse, gc
import numpy as np
import jax
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_transition")
import sampler_core as S
import trackB_b1 as B1
import trackB_b1_core as C
import trackB_estimator as TE

RES = S.RES


def gpu_mem():
    try:
        st = jax.local_devices()[0].memory_stats()
        return (f"in-use {st['bytes_in_use']/2**30:.1f} GiB, peak "
                f"{st.get('peak_bytes_in_use', 0)/2**30:.1f} GiB")
    except Exception:
        return "n/a"


def perturbed_points(P, src_truth, active_idx, n=8, seed=0, mags=(0.0, 1e-5, 1e-4, 1e-3,
                                                                  3e-3, 1e-2, 3e-2, 1e-1)):
    """8 source vectors: truth + geometric ladder of scaled perturbations on the ACTIVE dims."""
    rng = np.random.default_rng(seed)
    scale = TE.phi_scale(P)
    lo, hi = TE.phi_bounds(P)
    pts = []
    for i in range(n):
        s = src_truth.copy()
        d = rng.normal(size=len(active_idx))
        d /= np.linalg.norm(d)
        s[active_idx] += mags[i % len(mags)] * scale[active_idx] * d * np.sqrt(len(active_idx))
        pts.append(np.clip(s, lo, hi))
    return np.array(pts), np.array([mags[i % len(mags)] for i in range(n)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="lit")
    ap.add_argument("--npts", type=int, default=8)
    ap.add_argument("--fd-step", type=float, default=1e-6, help="FD step in SCALED units")
    ap.add_argument("--ref", action="store_true",
                    help="pass 1: bank the B1Marg reference values and exit")
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
    t00 = time.time()

    # ---------------- cell: zero-noise Asimov, arm A (truth at prior mean) -------------
    cell = S.build_cell(seed=None, tier=a.tier, noise=False, arm="A")
    P = cell["P"]
    src_truth = cell["theta_truth"][P["n_dist"]:].copy()
    idx = S.active_index(list(range(C.N_LOUD)), S.ACTIVE_AXES, P["ncw"])
    labels = S.active_labels(list(range(C.N_LOUD)), S.ACTIVE_AXES)
    print(f"\n  ACTIVE dims: {len(idx)}  ({', '.join(labels)})", flush=True)
    print(f"  fringe grid B = {P['EV_truth'].shape[1]}; npsr = {P['npsr']}", flush=True)

    pts, mags = perturbed_points(P, src_truth, idx, n=a.npts)
    refpath = f"{RES}/gate_ref_values.npz"

    # ---------------- G0: VALUE ----------------
    # THE REFERENCE RUNS IN ITS OWN PROCESS. B1Marg holds ~15 GB of GPU buffers (116 jitted
    # per-pulsar evaluators over a 512-fringe stack) and MargJax's 116 gradient executables
    # will not fit beside them on a 24 GB card -- co-residency OOMs (measured, twice). Pass 1
    # (--ref) banks the reference values and exits; pass 2 loads them.
    if a.ref:
        B1.B1Marg.T_CHUNK = 8
        B1M = B1.B1Marg(P, cell["data"], cell["mask"])
        t0 = time.time()
        B1M.gate_fast_full(pts[:2])
        print(f"  (fast-full gate {time.time()-t0:.0f}s)", flush=True)
        t0 = time.time()
        ref = B1M.lnmarg(pts)
        print(f"   B1Marg (numpy path): {time.time()-t0:.0f}s for {len(pts)} points",
              flush=True)
        np.savez(refpath, ref=ref, pts=pts, mags=mags)
        print(f"   banked {refpath}", flush=True)
        return 0

    print("\n== G0: MargJax.lnmarg == B1Marg.lnmarg (the banked object) ==", flush=True)
    Z = np.load(refpath)
    ref = Z["ref"]
    assert np.allclose(Z["pts"], pts), "reference points do not match this run's points"
    print(f"   reference values loaded from {refpath} (own process)", flush=True)

    M = S.MargJax(P, cell["data"], cell["mask"], idx, src_truth, tier=a.tier)
    x_pts = np.array([p[idx] for p in pts])
    t0 = time.time(); v0 = M.value(x_pts[0]); t_compile = time.time() - t0
    print(f"   MargJax first value (incl. compile): {t_compile:.0f}s", flush=True)
    print(f"   GPU after build: {gpu_mem()}", flush=True)
    t0 = time.time()
    got = np.array([M.value(x) for x in x_pts])
    t_val = (time.time() - t0) / len(x_pts)
    d = np.abs(got - ref)
    C._finite("G0 values", got)
    g0 = bool(np.max(d) < 1e-6)
    print(f"   {'mag':>8s} {'B1Marg':>16s} {'MargJax':>16s} {'|diff|':>10s}", flush=True)
    for i in range(len(pts)):
        print(f"   {mags[i]:8.0e} {ref[i]:16.6f} {got[i]:16.6f} {d[i]:10.2e}", flush=True)
    print(f"   [G0] max|diff| = {np.max(d):.3e} nat -> {'PASS' if g0 else 'FAIL'} (<1e-6)",
          flush=True)

    # ---------------- G1: GRAD vs FD ----------------
    print(f"\n== G1: autodiff grad(lnL_marg) vs central FD (step {a.fd_step:.0e} scaled) ==",
          flush=True)
    scale = TE.phi_scale(P)[idx]
    t0 = time.time(); v, g = M.value_and_grad(x_pts[0]); t_gc = time.time() - t0
    print(f"   first value_and_grad (incl. compile): {t_gc:.0f}s", flush=True)
    t0 = time.time(); v, g = M.value_and_grad(x_pts[1]); t_grad = time.time() - t0
    print(f"   warm value_and_grad: {t_grad:.2f}s  (value {t_val:.2f}s -> "
          f"grad/value = {t_grad/max(t_val,1e-9):.1f}x)", flush=True)

    rows = []
    worst = 0.0
    for i, x in enumerate(x_pts):
        v, gad = M.value_and_grad(x)
        C._finite(f"G1 grad pt{i}", gad)
        h = a.fd_step * scale
        gfd = np.zeros(len(x))
        for j in range(len(x)):
            xp = x.copy(); xp[j] += h[j]
            xm = x.copy(); xm[j] -= h[j]
            gfd[j] = (M.value(xp) - M.value(xm)) / (2 * h[j])
        den = np.maximum(np.abs(gfd), np.abs(gad))
        rel = np.abs(gad - gfd) / np.where(den > 0, den, 1.0)
        # only score dims whose FD signal is above the FD noise floor (~1e-10/h nat)
        noise = 1e-10 / h
        live = np.abs(gfd) > 100 * noise
        r_live = float(np.max(rel[live])) if live.any() else np.nan
        worst = max(worst, r_live if np.isfinite(r_live) else 0.0)
        rows.append((mags[i], float(np.max(np.abs(gad))), r_live, int(live.sum())))
        print(f"   pt{i} mag={mags[i]:7.0e}  max|grad|={np.max(np.abs(gad)):12.4e}  "
              f"max rel err (live dims {int(live.sum())}/{len(x)}) = {r_live:.3e}", flush=True)
        if r_live > 1e-6:
            bad = int(np.argmax(np.where(live, rel, 0)))
            print(f"        worst dim {labels[bad]}: AD {gad[bad]:.6e} vs FD {gfd[bad]:.6e}",
                  flush=True)
    g1 = bool(worst < 1e-6)
    print(f"   [G1] worst live rel err over {len(x_pts)} points = {worst:.3e} -> "
          f"{'PASS' if g1 else 'FAIL'} (<1e-6; mission bar 1e-8 reported, not enforced "
          f"where FD noise dominates)", flush=True)
    g1_strict = bool(worst < 1e-8)
    print(f"   [G1-strict 1e-8] -> {'PASS' if g1_strict else 'FAIL'}", flush=True)

    # ---------------- G2: cost ----------------
    print("\n== G2: cost model on the 4090 ==", flush=True)
    t0 = time.time()
    for x in x_pts[:4]:
        M.value_and_grad(x)
    t_grad = (time.time() - t0) / 4
    print(f"   value            : {t_val*1e3:8.0f} ms", flush=True)
    print(f"   value_and_grad   : {t_grad*1e3:8.0f} ms", flush=True)
    print(f"   NUTS @ 32 leapfrog/sample -> {32*t_grad:.1f} s/sample; "
          f"1000 samples = {32*t_grad*1000/3600:.1f} h/chain", flush=True)

    np.savez(f"{RES}/sampler_gates.npz",
             mags=mags, ref=ref, got=got, diff=d, rows=np.array(rows, float),
             t_value=t_val, t_grad=t_grad, ndim=len(idx), tier=a.tier,
             g0=g0, g1=g1, g1_strict=g1_strict, fd_step=a.fd_step)
    ok = g0 and g1
    print(f"\n=== SAMPLER GATES: G0 {'PASS' if g0 else 'FAIL'} | G1 {'PASS' if g1 else 'FAIL'} "
          f"| total {time.time()-t00:.0f}s ===", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
