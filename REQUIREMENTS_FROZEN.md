# REQUIREMENTS_FROZEN — reproducibility snapshot

*Frozen environment record for the banked Stage-C / Track-A / Track-B GPU results.*
*Generated 2026-07-07 by a hygiene pass (BROOM). Snapshot only — nothing rebuilt.*

The banked GPU compute (`CW_transition/stagec_*`, `trackA_*`, `trackB_*`) was produced in
the **`jug-gpu`** conda env on an **RTX 4090**. None of the versions below were pinned in
the repo before this file; reproducibility rested on prose in `project_progress.md` §9 plus
two external editable git checkouts. This file records the actual runtime state.

---

## 1. Editable forks (the load-bearing, unpinned dependencies)

| dep | checkout path | commit (canonical) | notes |
|---|---|---|---|
| **discovery** | `/home/mattm/soft/discovery` (editable, `src/discovery`) | **`136b270f`** (`136b270f1891c28ae6d9840930a9dfbcf41fd52d`) — merge commit "Merge branch 'main' of github.com/MattTMiles/discovery", 2026-03-17 01:01:48, branch `main`; parent = `81966f9` | ACTUAL installed/`main` checkout; imported as `discovery` / `discovery.matrix`. See note below. |
| **enterprise_extensions** | `/home/mattm/soft/enterprise_extensions` (editable) | `d43fef9` (`d43fef99d0aa786ec6b08dbc2fd56dcce10f26f5`) — "updates for det/CW signals", 2025-11-11 | reports `enterprise_extensions.__version__ == 3.0.3` |

Both checkouts were **clean** (no uncommitted changes) at snapshot time.

**discovery hash: `136b270f` installed, `81966f9` doc-cited (resolved with evidence 2026-07-07).**
The docs (`project_progress.md`, `MANIFEST.md`) cite `81966f9`; the actual installed checkout
(`pip freeze` editable pin + `git rev-parse main` + `git worktree list`) is at `136b270f`. These
are **two distinct commits**, NOT a replace-ref/ref-identity artifact (`git replace -l` is empty):

- `81966f9a056e13b92350b2ebb33f48de89c00496` — "updates to try and make things faster and more
  memory friendly", 2026-03-17 01:01:19.
- `136b270f1891c28ae6d9840930a9dfbcf41fd52d` — merge commit "Merge branch 'main' of
  github.com/MattTMiles/discovery", 2026-03-17 01:01:48 (29 s later).
- `git merge-base --is-ancestor 81966f9 136b270f` → **true**: `81966f9` is the parent; `main`
  drifted forward one merge commit past the doc-cited hash.

**`136b270f` is what actually ran** the banked Stage-C compute and is the hash to pin. `81966f9`
is its parent and effectively the same code tree (the merge added no divergent work — a
same-minute remote-`main` merge), so results are unaffected; recorded as the historical
doc-cited hash so the docs and this file reconcile.

Also editable in `jug-gpu` (present in freeze, not part of the Stage-C import path but recorded
for completeness): `jug @ 068098d`, `PINT @ ea7652f`.

---

## 2. Runtime env: `jug-gpu`

- conda root: `/home/mattm/miniforge3/envs/jug-gpu`
- **Python 3.12.13**
- **jax / jaxlib 0.10.1**, `jax-cuda12-plugin` / `jax-cuda12-pjrt` 0.10.1 (CUDA-12 wheels)
- numpy 2.4.6, scipy 1.17.1, enterprise-pulsar 3.4.4, jaxopt 0.8.5
- bundled NVIDIA CUDA-12 wheels: cublas 12.9.2.10, cudnn 9.23.1.3, cusolver 11.7.5.82,
  cufft 11.4.1.4, nccl 2.30.7, nvcc/nvrtc/runtime 12.9.x

Full `pip freeze` (227 pkgs) appended in **§6**. Regenerate with:

```
/home/mattm/miniforge3/envs/jug-gpu/bin/python -m pip freeze
```

(The `discotech` env — jax 0.9.0.1 — has a broken GPU stack per project memory; do NOT use it
for Stage C. Stage-A toy Fisher used a since-deleted ephemeral jax-0.4.28 scratch venv, see §5.)

---

## 3. Driver / hardware context

- **GPU:** NVIDIA GeForce RTX 4090, 24564 MiB (24 GB), shared multi-tenant card (`cronus`).
- **Driver:** 550.120.
- **CUDA runtime:** supplied by the jax CUDA-12 pip wheels (§2), not a system CUDA.

The 24 GB ceiling + shared co-tenancy drive the memory discipline in the Stage-C code (HVP
Hessian assembly, padded fixed-shape XLA graphs, autotune disabled — see §4).

---

## 4. XLA / jax conventions (set at import in every Stage-C script)

Verbatim from `CW_transition/stagec_fisher.py` (mirrored across the `stagec_*` / `trackA_*` /
`trackB_*` modules):

```python
import os
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"   # per-process BFC pool, no 90% grab
os.environ.pop("XLA_PYTHON_CLIENT_ALLOCATOR", None)
os.environ.setdefault("XLA_FLAGS", "--xla_gpu_autotune_level=0")  # autotune OFF: co-tenant
                                                                  # spikes OOM the compile
import jax
jax.config.update("jax_enable_x64", True)                # float64 THROUGHOUT — required;
                                                          # the pulsar-term phase 2πL/dL ~ 1.6e4
                                                          # rad needs f64 or the fringe eval is junk
```

- **`--xla_gpu_autotune_level=0`** (autotuning disabled): default kernels, robust to the
  shared-GPU co-tenant spikes that otherwise OOM the compile profiler. Peak-memory numbers
  are unaffected (per-process pool).
- **`XLA_PYTHON_CLIENT_PREALLOCATE=false`**: no upfront 90% VRAM grab; coexist with co-tenants.
- **`jax_enable_x64=True`**: **mandatory.** Distances / pulsar-term phases are evaluated at
  ~1.6e4 rad; f32 would destroy the ~1e-3-nat fringe-breaking signal.

---

## 5. Caches & regenerability

| cache | location | size | regenerable? |
|---|---|---|---|
| jax persistent compilation cache | `~/.cache/jax_stagec_cache` | **9.2 GB** | **YES** — set via `jax.config.update("jax_compilation_cache_dir", "/home/mattm/.cache/jax_stagec_cache")`. Pure compile cache; delete freely, it rebuilds on next cold compile (the ~350–500 s Stage-C compiles). Not required for correctness, only speed. |
| Stage-A toy jax-0.4.28 venv | ephemeral prior-session scratch | — | **LOST** — never in-repo. Rebuild recipe (prose only): jax 0.4.28 / CUDA 12.4 / cuDNN 8.9. Only affects re-running the toy Fisher (Stage A/A.1/A.2/A.3); banked `.npz`/`.png` outputs are on disk. |
| Cornell parallax fetch (`get_distance.py`) | live catalogue, cached outside repo | — | **NOT reproducible** — pulls the live catalogue; can drift. **Guard: `CW_transition/best_distances.txt` is the FROZEN canonical EM prior. Never re-run A0 for numbers.** |
| 116 pulsar feathers | `data_products/*.feather` | — | frozen by mtime (2025-10-14), not hashed; the array behind every banked result. |

---

## 6. Full `jug-gpu` pip freeze (2026-07-07)

```
absl-py==2.4.0
annotated-doc==0.0.4
annotated-types==0.7.0
astropy @ conda (astropy-base_1764120505923)
astropy-iers-data @ conda
colorama==0.4.6
-e git+https://github.com/MattTMiles/discovery.git@136b270f1891c28ae6d9840930a9dfbcf41fd52d#egg=discovery
emcee==(conda 1734122663166)
enterprise-pulsar==3.4.4
-e git+https://github.com/MattTMiles/enterprise_extensions.git@d43fef99d0aa786ec6b08dbc2fd56dcce10f26f5#egg=enterprise_extensions
ephem==4.2.1
fastapi==0.136.3
healpy==1.19.0
httptools==0.8.0
jax==0.10.1
jax-cuda12-pjrt==0.10.1
jax-cuda12-plugin==0.10.1
jaxlib==0.10.1
jaxopt==0.8.5
joblib==1.5.3
-e git+https://github.com/MattTMiles/jug.git@068098d74b8d59282154c2c410b3a133297fcb0f#egg=jug_timing
loguru==0.7.3
matplotlib==3.10.9
ml_dtypes==0.5.4
MultiHMCGibbs @ git+https://github.com/CKrawczyk/MultiHMCGibbs@3ffa31658753c82d1f5e9835aa1e532caaa98e26
multipledispatch==1.0.0
munkres==1.1.4
nestle==0.2.0
numdifftools==0.9.41
numpy==2.4.6 (conda wheel)
numpyro==0.21.0
nvidia-cublas-cu12==12.9.2.10
nvidia-cuda-cccl-cu12==12.9.27
nvidia-cuda-cupti-cu12==12.9.79
nvidia-cuda-nvcc-cu12==12.9.86
nvidia-cuda-nvrtc-cu12==12.9.86
nvidia-cuda-runtime-cu12==12.9.79
nvidia-cudnn-cu12==9.23.1.3
nvidia-cufft-cu12==11.4.1.4
nvidia-cusolver-cu12==11.7.5.82
nvidia-cusparse-cu12==12.5.10.65
nvidia-nccl-cu12==2.30.7
nvidia-nvjitlink-cu12==12.9.86
nvidia-nvshmem-cu12==3.7.0
opt_einsum==3.4.0
optax==0.2.8
-e git+https://github.com/MattTMiles/PINT.git@ea7652ff9e240e84b7b79fa3ed604581050b77a2#egg=pint_pulsar
pydantic==2.13.4
pydantic_core==2.46.4
pyqtgraph==0.14.0
PySide6==6.11.1
scikit-sparse==0.5.0
scipy==1.17.1 (conda wheel)
setuptools==80.10.2
shiboken6==6.11.1
starlette==1.3.1
uncertainties==3.2.3
uvicorn==0.49.0
uvloop==0.22.1
watchfiles==1.2.0
websockets==16.0
wheel==0.47.0
```

*Abridged to version-bearing lines; the ~150 conda-forge `@ file://…` transitive deps
(jupyter, arviz, dask, aiohttp, …) are omitted. Full verbatim output: rerun the `pip freeze`
command in §2 — the canonical source, not this file.*
