"""Stage C — Deliverable D1: Hessian infrastructure at scale.

stagec_hessian(...) returns the Hessian of the fast-scan discovery logL at
truth over [N_psr pulsar distances + 8*N_CW CW params].

Two paths:
  dense : jax.hessian (forward-over-reverse). Fine while the transient
          Jacobian stack fits (expected up to N_CW ~ 20 at N_psr=116).
  hvp   : Hessian-vector products (jax.jvp over jax.grad), one chunk of unit
          tangents at a time; assembles ONLY the distance block H_dd
          (N_psr x N_psr) and the distance-CW cross block H_dc
          (N_psr x 8*N_CW). Memory-bounded for large N_CW.

GATE (run `gate`): dense vs HVP blocks agree < 1e-8 on the D0 small case
(5 psr, 1 CW), and dense matches the D0 FD-validated Hessian.

Bench (run `bench --ncw N --method dense|hvp`): 116 psr, nb-05 realistic
config; reports wall times and peak device memory (jax live-buffer stats).

Run in jug-gpu:
  python stagec_hessian.py gate
  python stagec_hessian.py bench --ncw 4 --method hvp
"""

import os
# BFC allocator (not "platform") so device.memory_stats()['peak_bytes_in_use']
# tracks the true transient peak; preallocate off so the number is honest.
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
# cronus is a SHARED GPU (co-tenant jobs spike memory). XLA autotuning profiles
# many GEMM/Cholesky configs, each needing transient scratch — under a co-tenant
# spike that profiling OOMs and kills the compile. Disable it: default configs
# are picked instead (slightly slower kernels, but the run is robust to
# co-tenancy). Does NOT affect the per-process peak_bytes_in_use we report.
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import argparse
import sys
import time

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import (  # noqa: E402  (imports jax first; env above wins)
    load_pulsars, build_enterprise_pta, generate_injection_params,
    simulate, build_fast_scan_likelihood,
)

# nb-05 realistic operating point
SEED = 1234
N_COMPONENTS = 14
LOG10_EQUAD = -6.0
GWB_LOG10_A = -14.6
GWB_GAMMA = 13.0 / 3.0
LOG10_H_RANGE = (-14.0, -13.5)
LOG10_FGW_RANGE = (-8.0, -7.5)
LOG10_MC_RANGE = (8.5, 9.5)

CW_PARAM_BASE = [
    "cos_gwtheta", "gwphi", "cos_inc", "log10_mc",
    "log10_fgw", "log10_h", "phase0", "psi",
]


def device_mem():
    """(bytes_in_use, peak_bytes_in_use) on the first device."""
    stats = jax.local_devices()[0].memory_stats() or {}
    return stats.get("bytes_in_use", -1), stats.get("peak_bytes_in_use", -1)


def cw_param_keys(n_cw):
    suffixes = ["" if i == 0 else f"_{i+1}" for i in range(n_cw)]
    return [f"cw_{base}{suffix}" for suffix in suffixes for base in CW_PARAM_BASE]


def build_stagec_problem(n_psr, n_cw, seed=SEED, scenario="realistic",
                         log10_h=None, gwb_log10_A=GWB_LOG10_A):
    """nb-05-style problem at injected truth.

    Returns dict with logl_fn, param_keys, base_vals, ent_psrs, disco_psrs,
    inj_params, cw_block_names.
    """
    rng = np.random.default_rng(seed)
    np.random.seed(seed)

    ent_psrs, disco_psrs = load_pulsars(n_psr)
    pta, cw_block_names, _ = build_enterprise_pta(
        ent_psrs, n_cw, components=N_COMPONENTS, log10_equad=LOG10_EQUAD,
    )
    inj = generate_injection_params(
        pta, ent_psrs, n_cw, cw_block_names, log10_h=log10_h,
        scenario=scenario, rng=rng,
        log10_h_range=LOG10_H_RANGE, log10_fgw_range=LOG10_FGW_RANGE,
        log10_mc_range=LOG10_MC_RANGE,
        gwb_log10_A=gwb_log10_A, gwb_gamma=GWB_GAMMA,
    )
    sim = simulate(pta, inj)
    resid_map = {getattr(p, "name", p): y for p, y in zip(pta.pulsars, sim)}

    logl_fn, param_keys, base_vals = build_fast_scan_likelihood(
        disco_psrs, resid_map, n_cw, inj, cw_block_names,
        components=N_COMPONENTS, log10_equad=LOG10_EQUAD,
    )
    return dict(
        logl_fn=logl_fn, param_keys=param_keys, base_vals=base_vals,
        ent_psrs=ent_psrs, disco_psrs=disco_psrs, inj_params=inj,
        cw_block_names=cw_block_names,
    )


def select_params(problem, n_cw):
    """Indices of [N_psr distances + 8*N_CW CW params] in the param vector."""
    dist_keys = [f"{p.name}_cw_p_dist" for p in problem["ent_psrs"]]
    sel_keys = dist_keys + cw_param_keys(n_cw)
    pk = problem["param_keys"]
    missing = [k for k in sel_keys if k not in pk]
    if missing:
        raise KeyError(f"missing params: {missing}")
    sel_idx = np.array([pk.index(k) for k in sel_keys])
    return sel_keys, sel_idx, len(dist_keys)


def stagec_hessian(problem, n_cw, method="dense", chunk=16, at=None):
    """Hessian of the fast-scan logL at truth over [distances + CW params].

    method="dense": returns full (n_sel, n_sel) Hessian via jax.hessian.
    method="hvp"  : returns (n_dist, n_sel) array [H_dd | H_dc] — distance
                    rows only — via chunked Hessian-vector products.
    Also returns metadata dict (sel_keys, n_dist, timings).
    """
    sel_keys, sel_idx, n_dist = select_params(problem, n_cw)
    n_sel = len(sel_keys)
    base = jnp.asarray(problem["base_vals"])
    sel_idx_j = jnp.asarray(sel_idx)
    logl_fn = problem["logl_fn"]

    def f_sel(y):
        return logl_fn(base.at[sel_idx_j].set(y))

    y0 = jnp.asarray(np.array(problem["base_vals"])[sel_idx] if at is None else at)

    meta = dict(sel_keys=sel_keys, n_dist=n_dist, n_sel=n_sel, method=method)

    if method == "dense":
        hess_fn = jax.jit(jax.hessian(f_sel))
        t0 = time.time()
        H = np.array(hess_fn(y0))
        meta["t_cold"] = time.time() - t0
        t0 = time.time()
        H = np.array(hess_fn(y0))
        meta["t_warm"] = time.time() - t0
        return H, meta

    if method == "hvp":
        grad_fn = jax.grad(f_sel)

        def hvp(v):
            return jax.jvp(grad_fn, (y0,), (v,))[1]

        hvp_chunk = jax.jit(jax.vmap(hvp))
        eye = np.eye(n_sel)
        rows_dist = slice(0, n_dist)

        def assemble():
            cols = np.empty((n_dist, n_sel))
            for s in range(0, n_sel, chunk):
                e = min(s + chunk, n_sel)
                V = eye[s:e]
                if e - s < chunk:  # pad to fixed shape, avoid recompile
                    V = np.vstack([V, np.zeros((chunk - (e - s), n_sel))])
                out = np.array(hvp_chunk(jnp.asarray(V)))
                cols[:, s:e] = out[: e - s, rows_dist].T
            return cols

        t0 = time.time()
        Hd = assemble()
        meta["t_cold"] = time.time() - t0
        t0 = time.time()
        Hd = assemble()
        meta["t_warm"] = time.time() - t0
        return Hd, meta

    raise ValueError(method)


# ============================================================
# CLI
# ============================================================

def run_gate():
    """Dense vs HVP on the D0 small case; dense vs D0 FD-validated Hessian."""
    print("building D0 small case (5 psr, 1 CW, scenario=single)...", flush=True)
    t0 = time.time()
    problem = build_stagec_problem(5, 1, scenario="single", log10_h=-13.5)
    print(f"  built in {time.time()-t0:.1f}s, "
          f"lnL={float(problem['logl_fn'](problem['base_vals'])):.6f}")

    H_dense, md = stagec_hessian(problem, 1, method="dense")
    Hd_hvp, mh = stagec_hessian(problem, 1, method="hvp", chunk=8)
    n_dist, n_sel = md["n_dist"], md["n_sel"]

    ref = H_dense[:n_dist, :]  # [H_dd | H_dc] rows from dense
    Hmax = np.max(np.abs(H_dense))
    mask = np.abs(ref) > 1e-12 * Hmax
    rel = np.abs(Hd_hvp - ref)[mask] / np.abs(ref)[mask]
    gate = float(np.max(rel))
    print(f"dense vs HVP  [H_dd | H_dc] ({n_dist}x{n_sel}): "
          f"max rel err = {gate:.3e}  ({'PASS' if gate < 1e-8 else 'FAIL'} vs 1e-8)")

    # tie to the D0 FD-validated Hessian (identical build config)
    d0 = np.load("/home/mattm/projects/HSYMT/CW_transition/stagec_d0_results.npz")
    dev = np.max(np.abs(H_dense - d0["H_ad"]) / np.maximum(np.abs(d0["H_ad"]), 1e-12 * Hmax))
    print(f"dense vs D0 FD-validated H_ad: max rel dev = {dev:.3e}")
    return 0 if gate < 1e-8 else 1


def run_bench(n_cw, method, chunk):
    print(f"=== bench: N_psr=116, N_CW={n_cw}, method={method} ===", flush=True)
    t0 = time.time()
    problem = build_stagec_problem(116, n_cw)
    t_build = time.time() - t0
    ll = float(problem["logl_fn"](problem["base_vals"]))
    used0, peak0 = device_mem()
    print(f"build {t_build:.1f}s, lnL(truth)={ll:.3f}, "
          f"device mem after build: in_use={used0/2**30:.2f} GiB, peak={peak0/2**30:.2f} GiB",
          flush=True)

    H, meta = stagec_hessian(problem, n_cw, method=method, chunk=chunk)
    used1, peak1 = device_mem()
    n_dist, n_sel = meta["n_dist"], meta["n_sel"]
    print(f"H shape {H.shape} (n_dist={n_dist}, n_sel={n_sel})")
    print(f"t_cold={meta['t_cold']:.2f}s (incl compile)  t_warm={meta['t_warm']:.2f}s")
    print(f"device mem after hessian: in_use={used1/2**30:.2f} GiB, "
          f"PEAK={peak1/2**30:.2f} GiB")
    print(f"RESULT {n_cw} {method} t_warm={meta['t_warm']:.3f} "
          f"peak_GiB={peak1/2**30:.3f} build_s={t_build:.1f} cold_s={meta['t_cold']:.1f}")
    # sanity: diagonal of distance block should be finite, mostly negative
    diag = np.diag(H[:n_dist, :n_dist]) if method == "dense" else np.diag(H[:, :n_dist])
    print(f"dist-block diag: n_finite={np.isfinite(diag).sum()}/{n_dist}, "
          f"n_negative={(diag < 0).sum()}")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["gate", "bench"])
    ap.add_argument("--ncw", type=int, default=1)
    ap.add_argument("--method", default="dense", choices=["dense", "hvp"])
    ap.add_argument("--chunk", type=int, default=16)
    args = ap.parse_args()
    print(f"jax {jax.__version__}, devices {jax.devices()}", flush=True)
    if args.cmd == "gate":
        return run_gate()
    return run_bench(args.ncw, args.method, args.chunk)


if __name__ == "__main__":
    sys.exit(main())
