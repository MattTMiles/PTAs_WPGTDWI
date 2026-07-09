"""Track B — B0.0 cache verification.

Fresh process. Build the N_CW=16 amortised Fisher (the estimator's core graphs:
fisher_logL, fisher_logL_batched, fisher_hessian) and TIME the first (cache-miss-or-hit)
compile of each. PASS = seconds (persistent cache hit), not ~465 s (cold XLA compile).

Run in jug-gpu:  python trackB_b00_cache.py
"""
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")

import sys, time
import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")
jax.config.update("jax_persistent_cache_min_compile_time_secs", 10)
import jax.numpy as jnp

sys.path.insert(0, "/home/mattm/projects/HSYMT/CW_lnL_check")
from cw_helpers import load_pulsars, build_enterprise_pta
from stagec_fisher import (build_fisher_amortised, make_geometry_injection,
                           N_COMPONENTS, LOG10_EQUAD)
from stagec_d4 import assemble_hessian

EQUAL_H = -13.75
NCW = int(os.environ.get("B00_NCW", "16"))

print(f"jax {jax.__version__}, {jax.devices()}", flush=True)
print(f"cache dir = {jax.config.jax_compilation_cache_dir}", flush=True)
print(f"min_compile_time_secs = {jax.config.jax_persistent_cache_min_compile_time_secs}", flush=True)

t0 = time.time()
ent_psrs, disco_psrs = load_pulsars(116)
pta, cwb, _ = build_enterprise_pta(ent_psrs, NCW, components=N_COMPONENTS, log10_equad=LOG10_EQUAD)
inj0 = make_geometry_injection(pta, ent_psrs, NCW, cwb, seed=1000, h_range=(EQUAL_H, EQUAL_H))
amo = build_fisher_amortised(disco_psrs, NCW, inj0, cwb)
print(f"[build] load+pta+build_fisher_amortised (no compile yet): {time.time()-t0:.1f}s", flush=True)

th0 = amo["theta_truth"]
dt0 = amo["inject_data"](th0)

# time each graph's FIRST call (compile-or-cache-load)
t = time.time(); _ = float(amo["fisher_logL"](th0, dt0)); t_logL = time.time()-t
print(f"[compile] fisher_logL          (N_CW={NCW}): {t_logL:8.2f}s", flush=True)

t = time.time()
_ = np.asarray(amo["fisher_logL_batched"](jnp.zeros((4, amo["n_theta"])), dt0))
t_batch = time.time()-t
print(f"[compile] fisher_logL_batched  (N_CW={NCW}): {t_batch:8.2f}s", flush=True)

t = time.time(); H = assemble_hessian(amo, th0, dt0); t_hess = time.time()-t
print(f"[compile] assemble_hessian/HVP (N_CW={NCW}): {t_hess:8.2f}s", flush=True)

# warm re-eval (should be ~ms-100ms)
t = time.time(); _ = float(amo["fisher_logL"](th0, dt0)); t_warm = time.time()-t
print(f"[warm]    fisher_logL re-eval:  {t_warm*1000:.1f} ms", flush=True)

worst = max(t_logL, t_batch, t_hess)
verdict = "PASS (cache hit, seconds)" if worst < 60 else "FAIL (cold compile, ~minutes)"
print(f"\n[B0.0 cache] worst first-compile = {worst:.1f}s -> {verdict}", flush=True)
