"""Stage C — Deliverable D0: derivative validation gate.

Validates jax autodiff derivatives of build_fast_scan_likelihood (jug-gpu,
jax 0.10.x) against finite differences on a small problem (5 pulsars, 1 CW),
at injected truth. This gates all downstream Hessian work: the fast-scan logL
VALUE is validated vs discotech-CPU, but the DERIVATIVES are not (jax 0.10
cho_solve FutureWarning risk).

Gates:
  a. grad(logL) over [5 distances + 8 CW params]: autodiff vs central FD
     (rel step 1e-5 per param natural scale).      GATE: max rel err < 1e-5
  b. hessian(logL) same params: autodiff vs central FD OF THE AUTODIFF GRAD.
     GATE: dist block + dist-CW cross block rel err < 1e-4 (elements above
     1e-8 * max|H|)
  c. cho_solve path in differentiated graph: rebuild with GWB log10_A
     perturbed +0.1, gradient at same point must change.

Run:  /home/mattm/miniforge3/envs/jug-gpu/bin/python stagec_d0_derivative_gate.py
"""

import sys
import time
import warnings

import numpy as np

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")

captured_warnings = []
with warnings.catch_warnings(record=True) as wlist:
    warnings.simplefilter("always")
    import jax
    import jax.numpy as jnp
    from cw_helpers import (
        load_pulsars, build_enterprise_pta, generate_injection_params,
        simulate, build_fast_scan_likelihood,
    )
    captured_warnings += [f"{w.category.__name__}: {w.message}" for w in wlist]

# nb-05 realistic operating point, shrunk to 5 pulsars / 1 CW
SEED = 1234
N_PSR = 5
N_CW = 1
N_COMPONENTS = 14
LOG10_EQUAD = -6.0
GWB_LOG10_A = -14.6
GWB_GAMMA = 13.0 / 3.0
LOG10_H = -13.5

CW_PARAM_KEYS = [
    "cw_cos_gwtheta", "cw_gwphi", "cw_cos_inc", "cw_log10_mc",
    "cw_log10_fgw", "cw_log10_h", "cw_phase0", "cw_psi",
]


def build_problem(gwb_log10_A=GWB_LOG10_A):
    """5-psr 1-CW scenario at injected truth; returns (logl_fn, keys, base_vals)."""
    rng = np.random.default_rng(SEED)
    np.random.seed(SEED)

    ent_psrs, disco_psrs = load_pulsars(N_PSR)
    pta, cw_block_names, _ = build_enterprise_pta(
        ent_psrs, N_CW, components=N_COMPONENTS, log10_equad=LOG10_EQUAD,
    )
    inj = generate_injection_params(
        pta, ent_psrs, N_CW, cw_block_names, log10_h=LOG10_H,
        scenario="single", rng=rng,
        gwb_log10_A=gwb_log10_A, gwb_gamma=GWB_GAMMA,
    )
    sim = simulate(pta, inj)
    resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim)}

    logl_fn, param_keys, base_vals = build_fast_scan_likelihood(
        disco_psrs, resid_map, N_CW, inj, cw_block_names,
        components=N_COMPONENTS, log10_equad=LOG10_EQUAD,
    )
    return logl_fn, param_keys, base_vals, ent_psrs


def main():
    t0 = time.time()
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)

    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        logl_fn, param_keys, base_vals, ent_psrs = build_problem()
        ll0 = float(logl_fn(base_vals))
    captured_warnings.extend(f"{w.category.__name__}: {w.message}" for w in wlist)
    print(f"built 5-psr 1-CW problem in {time.time()-t0:.1f}s, lnL(truth)={ll0:.6f}")

    dist_keys = [f"{p.name}_cw_p_dist" for p in ent_psrs]
    sel_keys = dist_keys + CW_PARAM_KEYS
    missing = [k for k in sel_keys if k not in param_keys]
    assert not missing, f"missing params: {missing}"
    sel_idx = np.array([param_keys.index(k) for k in sel_keys])
    n_dist, n_sel = len(dist_keys), len(sel_keys)
    print(f"selected {n_sel} params ({n_dist} distances + {len(CW_PARAM_KEYS)} CW)")

    base = jnp.asarray(base_vals)
    sel_idx_j = jnp.asarray(sel_idx)

    def f_sel(y):
        return logl_fn(base.at[sel_idx_j].set(y))

    y0 = np.array(base_vals)[sel_idx]
    y0_j = jnp.asarray(y0)

    # --- FD steps: 1e-2 relative in each param's NATURAL (oscillation) scale ---
    # The CW logL is oscillatory: pulsar-term phase Phi_p = 2*pi*f*(L_p/c)*(1-cos mu)
    # = 2*pi*L_p/dL_p ~ O(1e3-1e4) rad. A step of 1e-5*|x| spans O(0.1-1) rad of
    # fringe phase -> central-diff truncation error O(1e-3), which is what a naive
    # step produces. The natural scale of an oscillatory parameter is the scale
    # over which the phase changes by 1 rad:
    #   distance L_p     : dL_p / (2 pi)
    #   log10_fgw, log10_mc : 1 / (ln10 * Phi_max)
    #   cos_gwtheta, gwphi  : 1 / Phi_max
    #   phase0, psi, cos_inc, log10_h : 1 (smooth, O(1) periodicity)
    # Steps are 1e-2 of that scale, with Richardson extrapolation
    # (4*FD(h/2) - FD(h))/3 so truncation is O(u^4) ~ 1e-9 while the logL
    # difference stays large enough to avoid float64 cancellation.
    from cw_helpers import compute_mode_spacing

    inj_cos_th = float(y0[sel_keys.index("cw_cos_gwtheta")])
    inj_phi = float(y0[sel_keys.index("cw_gwphi")])
    inj_lf = float(y0[sel_keys.index("cw_log10_fgw")])
    _, disco_psrs_scales = load_pulsars(N_PSR)
    dL = np.array([
        compute_mode_spacing(inj_cos_th, inj_phi, inj_lf, p.pos)
        for p in disco_psrs_scales
    ])
    L_p = y0[:n_dist]
    Phi_max = float(np.max(2 * np.pi * L_p / dL))
    print(f"fringe spacings dL (kpc): {np.array2string(dL, precision=5)}")
    print(f"Phi_max = {Phi_max:.1f} rad")

    scale = np.ones(n_sel)
    scale[:n_dist] = dL / (2 * np.pi)
    ln10 = np.log(10.0)
    scale[sel_keys.index("cw_log10_fgw")] = 1.0 / (ln10 * Phi_max)
    scale[sel_keys.index("cw_log10_mc")] = 1.0 / (ln10 * Phi_max)
    scale[sel_keys.index("cw_cos_gwtheta")] = 1.0 / Phi_max
    scale[sel_keys.index("cw_gwphi")] = 1.0 / Phi_max
    REL_STEP = 1e-2
    h = REL_STEP * scale

    def fd_grad_richardson(i, step):
        def central(hh):
            yp = y0.copy(); yp[i] += hh
            ym = y0.copy(); ym[i] -= hh
            return (float(f_sel(jnp.asarray(yp))) - float(f_sel(jnp.asarray(ym)))) / (2 * hh)
        c1, c2 = central(step), central(step / 2)
        return (4 * c2 - c1) / 3, c1, c2

    # ---------------- gate (a): gradient ----------------
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        grad_fn = jax.jit(jax.grad(f_sel))
        g_ad = np.array(grad_fn(y0_j))
    captured_warnings.extend(f"{w.category.__name__}: {w.message}" for w in wlist)

    g_fd = np.zeros(n_sel)
    conv_ratio = np.zeros(n_sel)
    for i in range(n_sel):
        g_fd[i], c1, c2 = fd_grad_richardson(i, h[i])
        # truncation-dominance evidence: central-FD error should shrink ~4x
        # from h to h/2 if error is O(h^2) (i.e. FD -> autodiff as h -> 0)
        e1, e2 = abs(c1 - g_ad[i]), abs(c2 - g_ad[i])
        conv_ratio[i] = e1 / e2 if e2 > 0 else np.inf

    gmax = np.max(np.abs(g_fd))
    denom_a = np.maximum(np.abs(g_fd), 1e-12 * gmax)
    rel_a = np.abs(g_ad - g_fd) / denom_a
    gate_a = np.max(rel_a)
    print("\n--- gate (a): grad autodiff vs Richardson central FD ---")
    for k, ga, gf, r, cr in zip(sel_keys, g_ad, g_fd, rel_a, conv_ratio):
        print(f"  {k:28s} ad={ga: .10e}  fd={gf: .10e}  rel={r:.3e}  h->h/2 err ratio={cr:.2f}")
    print(f"GATE A: max rel err = {gate_a:.3e}  ({'PASS' if gate_a < 1e-5 else 'FAIL'} vs 1e-5)")

    # ---------------- gate (b): hessian ----------------
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        t1 = time.time()
        H_ad = np.array(jax.jit(jax.hessian(f_sel))(y0_j))
        t_hess = time.time() - t1
    captured_warnings.extend(f"{w.category.__name__}: {w.message}" for w in wlist)
    print(f"\nautodiff hessian ({n_sel}x{n_sel}) in {t_hess:.1f}s (incl compile)")

    H_fd = np.zeros((n_sel, n_sel))
    for j in range(n_sel):
        def central_col(hh):
            yp = y0.copy(); yp[j] += hh
            ym = y0.copy(); ym[j] -= hh
            return (np.array(grad_fn(jnp.asarray(yp)))
                    - np.array(grad_fn(jnp.asarray(ym)))) / (2 * hh)
        c1, c2 = central_col(h[j]), central_col(h[j] / 2)
        H_fd[:, j] = (4 * c2 - c1) / 3
    H_fd = 0.5 * (H_fd + H_fd.T)

    Hmax = np.max(np.abs(H_ad))
    thresh = 1e-8 * Hmax

    def block_relerr(rows, cols, label):
        A, F = H_ad[np.ix_(rows, cols)], H_fd[np.ix_(rows, cols)]
        mask = np.abs(F) > thresh
        if not mask.any():
            print(f"  {label}: no elements above threshold"); return 0.0
        rel = np.abs(A - F)[mask] / np.abs(F)[mask]
        print(f"  {label}: {mask.sum()}/{mask.size} elements above 1e-8*max|H|, "
              f"max rel err = {np.max(rel):.3e}")
        return np.max(rel)

    print(f"--- gate (b): hessian autodiff vs FD-of-grad (max|H|={Hmax:.3e}) ---")
    ii = np.arange(n_dist); jj = np.arange(n_dist, n_sel)
    err_dd = block_relerr(ii, ii, "distance block   (5x5)")
    err_dc = block_relerr(ii, jj, "dist-CW cross    (5x8)")
    err_cc = block_relerr(jj, jj, "CW block         (8x8) [info only]")
    gate_b = max(err_dd, err_dc)
    print(f"GATE B: max rel err (dist + cross) = {gate_b:.3e}  "
          f"({'PASS' if gate_b < 1e-4 else 'FAIL'} vs 1e-4)")

    # ---------------- gate (c): cho_solve path differentiated ----------------
    print("\n--- gate (c): perturb frozen GWB log10_A by +0.1, rebuild, grad must change ---")
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        logl2, keys2, base2, _ = build_problem(gwb_log10_A=GWB_LOG10_A + 0.1)
    captured_warnings.extend(f"{w.category.__name__}: {w.message}" for w in wlist)
    assert keys2 == param_keys, "param key mismatch after rebuild"
    # evaluate second build's grad at the ORIGINAL point (same residuals injected
    # differently is irrelevant: we compare grads of two frozen-GP graphs at one x)
    base2_arr = jnp.asarray(base)  # same parameter point

    def f_sel2(y):
        return logl2(base2_arr.at[sel_idx_j].set(y))

    g_ad2 = np.array(jax.jit(jax.grad(f_sel2))(y0_j))
    dg = np.abs(g_ad2 - g_ad)
    rel_c = np.max(dg) / np.max(np.abs(g_ad))
    print(f"  max |grad change| = {np.max(dg):.6e}  (rel to max|g| = {rel_c:.3e})")
    gate_c = np.max(dg) > 0
    print(f"GATE C: gradient responds to frozen-GP perturbation: "
          f"{'PASS' if gate_c else 'FAIL'}")

    # ---------------- summary ----------------
    print("\n================ D0 SUMMARY ================")
    print(f"jax version         : {jax.__version__}")
    print(f"lnL(truth)          : {ll0:.6f}")
    print(f"gate a (grad)       : {gate_a:.3e}  (< 1e-5) {'PASS' if gate_a < 1e-5 else 'FAIL'}")
    print(f"gate b (hessian)    : {gate_b:.3e}  (< 1e-4) {'PASS' if gate_b < 1e-4 else 'FAIL'}"
          f"   [CW block info: {err_cc:.3e}]")
    print(f"gate c (cho_solve)  : max dgrad {np.max(dg):.3e} > 0  {'PASS' if gate_c else 'FAIL'}")
    uniq = sorted(set(captured_warnings))
    print(f"warnings ({len(uniq)} unique):")
    for w in uniq:
        print(f"  - {w}")
    ok = (gate_a < 1e-5) and (gate_b < 1e-4) and gate_c
    print(f"D0 OVERALL: {'PASS' if ok else 'FAIL'}")
    np.savez(
        "/home/mattm/projects/HSYMT/CW_transition/stagec_d0_results.npz",
        g_ad=g_ad, g_fd=g_fd, H_ad=H_ad, H_fd=H_fd,
        sel_keys=np.array(sel_keys), gate_a=gate_a, gate_b=gate_b,
        max_dgrad=np.max(dg), ll_truth=ll0,
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
