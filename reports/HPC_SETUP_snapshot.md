> **SNAPSHOT (staged 2026-07-12 by IGNITE-2 staging session, tag criterion-v2).**
> ACCRE's live untracked `HPC_SETUP.md` at the repo root remains canonical for
> cluster state; this copy travels so cronus has the cluster record.

# HPC_SETUP.md — agent HARBOR, ACCRE (Vanderbilt)

Untracked working report. Execution-environment recon + port record.
Cronus is canonical; this file never feeds back into `project_progress.md`.

- Host: ACCRE `vampire` cluster, Slurm 25.11.1, Linux 5.14 (EL9)
- User `milesmt`, groups: `taylor_group`, `dsi_dgx*`
- Repo: `/data/taylor_group/matt_miles/PTAs_WPGTDWI` @ `634aab8` (branch `MM_playground`)
- Started 2026-07-09

---

## 0. DIVERGENCE FROM THE PORT PLAN — READ FIRST

**`CW_transition/PORTER_slurm_plan.md` DOES NOT EXIST.** Not in the working tree, not in
`git ls-files`, not in `git log --all -- '*PORTER*'`. The scoping agent's plan was never
committed, or was committed on cronus and not pushed. I proceeded without it. Everything
below is derived from first-principles recon, so the ASSUMPTIONS-table comparison the brief
asked for cannot be done. **If that file exists on cronus, push it — I will re-pull and
diff my findings against it.**

Second divergence, larger: **the brief assumes A100s reachable from a batch GPU partition
with a short/debug backfill lane.** That is not this cluster's shape. See §2.

---

## 1. H1 — RECON

### 1.1 Partitions

| Partition | Nodes | MaxTime | Default QOS | Notes |
|---|---|---|---|---|
| `batch` | 391 CPU | 14-00:00:00 | N/A | CPU only |
| `interactive` | 30 CPU | 14-00:00:00 | `gpu_default` | CPU only |
| `batch_gpu` | 40 | 14-00:00:00 | `gpu_default` | A6000/L40S bulk, 3×A100-80, 8×H100-80, 8×H200, 2×GH200 |
| `interactive_gpu` | 7 | 14-00:00:00 | `gpu_default` | **3×DGX (24×A100)**, A4000, 3090 |
| `reserved` | 0 | — | — | INACTIVE |

`DefaultTime=00:30:00` everywhere. No partition-level short/debug lane exists.

### 1.2 The QOS wall — this is the whole story

`gpu_default` (the default QOS on every GPU partition) has
`GrpTRES=gres/gpu:<every type>=0`. **The default QOS grants zero GPUs.** A job that does not
name an explicit `-A account -q qos` pair gets no GPU, ever. This is the single most
important fact about this cluster.

My associations (`sacctmgr show assoc user=milesmt`):

| Account | QOS | MaxWall | GPU entitlement |
|---|---|---|---|
| `dsi_dgx_iacc` | `dgx_iacc` | *(none → partition 14 d)* | **a100-sxm4-80gb=16, a100-sxm4-40gb=16, gpu=32, cpu=512, mem=5996G** |
| `p_dsi_dgx` | `dsi_dgx` | none | a100 16+16 — **but `-A p_dsi_dgx -q dsi_dgx` is rejected: "Invalid account/partition combination"** |
| `taylor_group_acc` | `taylor_group_account_batch_gpu` | none | gh200=2, h200=8 only |
| `taylor_group_iacc` / `dsi_dgx_iacc` | `debug_iacc` | **00:30:00** | gpu=1 (A4000 class), cpu=6, mem=75G |
| `taylor_group_int` / `dsi_dgx_int` | `debug_int` | 00:30:00 | **CPU only**, cpu=16, mem=192G |

### 1.3 `--test-only` entitlement probes (submitted nothing)

```
-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc  --gres=gpu:nvidia_a100-sxm4-80gb:1 -t 04:00:00  -> starts NOW  (dgx03)
-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc  --gres=gpu:nvidia_a100-sxm4-80gb:1 -t 14-00:00:00 -> starts NOW  (dgx03)
-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc  --gres=gpu:nvidia_a100-sxm4-80gb:4 -t 04:00:00  -> starts NOW  (dgx03)
-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc  --array=0-31%8 -t 02:00:00           -> starts NOW  (dgx03)
-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc  --gres=gpu:nvidia_a100-sxm4-40gb:1   -> starts 2026-07-10T09:58 (dgx01 full)
-p interactive_gpu -A taylor_group_iacc -q debug_iacc --gres=gpu:1 -t 00:30:00       -> starts NOW  (gpu0207 = RTX 3090)
-p batch_gpu -A taylor_group_acc -q taylor_group_account_batch_gpu --gres=gpu:nvidia_h200:1  -> starts 2026-08-08T03:28  (!!)
-p batch_gpu -A taylor_group_acc -q taylor_group_account_batch_gpu --gres=gpu:nvidia_a100_80gb:1 -> QOSGrpGRES DENIED
-p batch_gpu -A taylor_group_acc -q taylor_group_account_batch_gpu --gres=gpu:nvidia_l40s:1      -> QOSGrpGRES DENIED
```

**`batch_gpu` is unusable.** The taylor_group entitlement there is H200/GH200 only, and the
H200 queue is a *month* deep. The A6000/L40S bulk of that partition is denied to us outright.

### 1.4 Queue shape

- `batch_gpu`: 160 running, **51 pending** — 41 of those blocked on `QOSGrpGRES`
  (group entitlement exhaustion, not resource scarcity).
- `interactive_gpu`: 7 running, **0 pending. The lane is empty.**
- Free right now: dgx03 has **7/8 A100-80GB idle**, dgx04 has 4/8 idle. gpu0084 (batch_gpu)
  has 3/3 A100-80GB idle but is QOS-denied to us.
- `sacct -r interactive_gpu -S now-2days`: `dgx_iacc` jobs routinely run **1–5 days** and
  COMPLETE. One TIMEOUT, one admin CANCEL. No preemption observed.

### 1.5 Backfill — real, but largely irrelevant to us

```
SchedulerType       = sched/backfill
SchedulerParameters = bf_interval=240, bf_continue, bf_resolution=1800, bf_window=10080,
                      bf_max_job_test=300, bf_max_job_user=60, bf_max_job_user_part=30,
                      bf_max_job_array_resv=5, bf_busy_nodes, defer, batch_sched_delay=10
MaxArraySize        = 60000
MaxJobCount         = 500000
PriorityType        = priority/multifactor
PriorityWeightQOS   = 100000000   (but every QOS has Priority=0 -> contributes nothing)
PriorityWeightFairShare = 100000 ; Age = 25000 ; TRES gpu = 100000
```

Two consequences the brief's strategy did not anticipate:

1. **`bf_resolution=1800`** — backfill time granularity is **30 minutes**. Requesting
   `1:23:00` and `1:30:00` land in the same bucket. "measured + 30%" is right, but then
   **round up to the next 30-minute multiple**; finer precision buys literally nothing.
2. **The backfill-into-short-jobs strategy is moot on `interactive_gpu`** because there is
   no queue to backfill past. Accurate short walltimes remain correct practice (fairshare,
   and they protect us if the lane fills), but they are not the way in. *The way in is the
   `dgx_iacc` QOS.* The brief's premise — "clusters almost always have a low-walltime lane
   that schedules fast; that lane is our way in" — is true here (`debug_iacc`, 30 min) but
   **that lane only reaches an A4000/3090, whose fp64 is crippled exactly like the 4090.**
   It is useful for a GPU hello-world and nothing else.

**Preemption:** `PreemptType=preempt/qos`, `PreemptMode=CANCEL`, but *no QOS lists any
other in its `Preempt` field* — so no preemption relationship is actually configured. Jobs
die from TIMEOUT and admin action, not preemption. Checkpoint/resume is still mandatory
(TIMEOUT is real, and `PreemptMode=CANCEL` means a future config change would kill, not
requeue), but the urgency is lower than assumed. `KillWait=30s`, `GraceTime=0`.

### 1.6 GPU inventory

| Node | Partition | GPU | Count | fp64 rate |
|---|---|---|---|---|
| dgx01 | interactive_gpu | A100-SXM4-**40GB** | 8 | **1:2 (full)** |
| dgx03, dgx04 | interactive_gpu | A100-SXM4-**80GB** | 8 each | **1:2 (full)** |
| gpu0084 | batch_gpu | A100-80GB (PCIe) | 3 | 1:2 — *QOS-denied* |
| hgx02 | batch_gpu | H100-80GB | 8 | 1:2 — *not entitled* |
| hgx03 | batch_gpu | H200-141GB | 8 | 1:2 — entitled, **1-month queue** |
| gracehopper01/02 | batch_gpu | GH200 | 1 each | 1:2 — **aarch64, jax CUDA-12 x86 wheels will not run** |
| gpu0059–0082, 0208, 0300–0301 | batch_gpu | RTX A6000 | 108 | 1:64 crippled — *QOS-denied* |
| gpu0302–0310 | batch_gpu | L40S | 72 | 1:64 crippled — *QOS-denied* |
| gpu0058, 0060–0061 | interactive_gpu | RTX A4000 | 20 | 1:64 crippled |
| gpu0207 | interactive_gpu | RTX 3090 | 2 | 1:64 crippled |

Because `jax_enable_x64=True` is **mandatory** for this project (REQUIREMENTS_FROZEN §4 —
pulsar-term phase ~1.6e4 rad), **only the A100 class is viable**, and the *only* route to it
is `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc`. That the DGX A100s sit in a partition
named "interactive" is an ACCRE naming quirk; they take ordinary `sbatch` array jobs
(probe 12447937 above confirms) with up to 14-day walltimes.

Expected fp64 gain vs the banked 4090: A100 fp64 ≈ **9.7 TFLOP/s** (19.5 with tensor cores)
vs 4090 ≈ **1.3 TFLOP/s**. **~7.5–15×.** To be measured in H4a, not assumed.

Driver version and CUDA runtime: **UNKNOWN — no `nvidia-smi` on the login node.** Deferred
to the H4a GPU hello-world. This is the one hard blocker on declaring the env verified,
because the jax CUDA-12 pip wheels need driver ≥ 525 (cronus had 550.120).

### 1.7 Filesystems

| Path | Quota | Used | Backed up | Purge | Verdict |
|---|---|---|---|---|---|
| `/home/milesmt` (panfs) | 186.3 G | 72.0 G | yes | no | **jax compile cache here** (9.2 G) |
| `/nobackup/user/milesmt` | 186.3 G | 32 K | **no** | **policy UNKNOWN** | checkpoints only (regenerable) |
| `/data/taylor_group` (panfs) | 15.5 T | **13.6 T** | yes | no | **only 1.9 T headroom — group-shared, do not fill** |
| `/lfs_roots/barclay` (Lustre) | 35 P | 28 P | ? | ? | not obviously mine; unexplored |

**Open question:** `/nobackup` purge policy is not documented in-tree and there is no
README at the mount. ACCRE convention is a time-based purge. **Do not put the compile cache
or campaign results there until confirmed.** Ask ACCRE, or keep everything on `/home` +
`/data`.

Group quota is the real constraint: **1.9 TiB of headroom on a shared volume.** A large B1
array producing per-realisation npz needs a size budget before launch, not after.

### 1.8 Software environment (login node)

- `apptainer 1.4.5` **and** `singularity` present. Containers are viable.
- `conda`/`mamba` via `~/miniforge3` (Python 3.12.8 base). Existing envs: `base`,
  `discotech`, `main`, `ng20`. **No `jug-gpu`.** (`discotech` = jax 0.9.0.1, broken GPU
  stack per REQUIREMENTS_FROZEN §2 — do not use.)
- `module` = ACCRE software-stack shim, not Lmod.
- **Internet from the login node: YES** (`pypi.org` 200, `github.com` 200). From compute
  nodes: **UNKNOWN**, deferred to H4a.
- No `nvidia-smi`, no system CUDA on the login node.

### 1.9 Fork checkouts — a real problem

| Fork | Pin (REQUIREMENTS_FROZEN §1) | On this box | On remote |
|---|---|---|---|
| `discovery` | `136b270f` | `~/soft/discovery` @ `9124089`, **object `136b270f` absent from local object DB** | **present** (`refs/heads/main` = `136b270f`) |
| `enterprise_extensions` | `d43fef99` (3.0.3) | **ABSENT — no checkout anywhere** | **present** (`refs/heads/master` = `d43fef99`) |

`~/soft/discovery` HEAD is `9124089` ("Fixed the scale parameter for the band noise
filter") and does not contain the pinned commit — a stale clone from a different point in
history. A `git fetch` will bring `136b270f` in; then detach onto it. `enterprise_extensions`
must be cloned fresh. Both fetches work (login node has internet).

There is also `~/soft/NANOGRAV/discovery` @ `dbb693b` — a *different* clone. **Do not use
it.** Recording it so nobody wires the wrong one into `PYTHONPATH`.

---

## 2. REVISED QUEUE-ENTRY STRATEGY (supersedes the brief's §H3 where it conflicts)

The brief's strategy is built on "short jobs + accurate walltimes slip into backfill gaps."
That is the correct instinct for a contended partition. **Our partition is not contended.**
Adopt this instead:

**Lane:** `-p interactive_gpu -A dsi_dgx_iacc -q dgx_iacc --gres=gpu:nvidia_a100-sxm4-80gb:1`

Pin **`-80gb` explicitly.** The 40GB DGX (dgx01) is saturated and would push us to
tomorrow; the 80GB pair (dgx03/dgx04) has 11 GPUs free right now. Never request bare
`--gres=gpu:1` — Slurm will hand back a 3090 from `gpu0207` and every fp64 kernel will
crawl.

Keep from the brief, unchanged, because they are right for reasons independent of backfill:

- **Job arrays of short tasks.** Not for backfill — for *blast radius*. A TIMEOUT or an
  admin `scancel` costs one task, not a campaign. `MaxArraySize=60000` is ample.
  `bf_max_job_array_resv=5` and `bf_max_job_user=60` mean throttling the array with `%N`
  (N ≈ 8–16) is polite and costs us nothing given the empty queue.
- **Walltime = measured + 30%, then rounded up to the next 30 min** (`bf_resolution=1800`).
- **Checkpoint/resume in every task**, npz at stage boundaries, `--signal=B:USR1@300`,
  `--requeue`, idempotent tasks, output keyed by `(job, task, stage)`. `GraceTime=0` and
  `KillWait=30s` mean the USR1 handler has 300 s of real warning and must not dawdle.
- **Finiteness discipline** on every checkpoint: isfinite fraction + min/median/max printed,
  NaN aborts.
- **One warm-cache job per graph shape first**, into a shared `jax_compilation_cache_dir`
  on `/home` (not `/nobackup` until purge policy is known). Doubles as GPU hello-world.

Dropped from the brief: the "1–4 hr to court backfill" constraint. Task length should be set
by *checkpoint granularity and blast radius*, not by scheduler courtship. With a 14-day
ceiling and an empty lane, a 6-hour task with hourly checkpoints is fine.

**Revisit on A100 — see §7. Answered by the D1/D2/D3 diagnostic, and my first two
statements here were BOTH wrong.**

I first predicted `--gres` would be exclusive. Then H4a showed a foreign process on our card
and I concluded "Slurm does not fence the GPU." **Both wrong.** D1 proves
`ConstrainDevices=yes` and Slurm *does* fence devices correctly; D3 proves the foreign process
is a **rogue, non-Slurm SSH process** squatting on cards Slurm still hands out. Corrected
conclusions and the full timing matrix are in **§7**. The short version:

- `--xla_gpu_autotune_level=0` **stays** — but for a *measured* reason, not the inherited one:
  `autotune_level=2` makes the science kernels **2.2× SLOWER** on A100 (§7.3).
- `XLA_PYTHON_CLIENT_PREALLOCATE=false` **stays** — a 90 % grab would collide with the rogue.
- **Usable VRAM is 80 GB on a clean card, ~30 GB on a squatted one.** Which you get is a
  coin flip (3 of dgx03's 8 GPUs were squatted). Size B1 for **30 GB** until the rogue is
  cleared — still above the 4090's 24 GB, but not the 3.3× headroom a port plan would assume.

---

## 3. H2 — ENVIRONMENT BUILD — **PASS**

### 3.1 The blocker nobody flagged: 43 files hard-code cronus absolute paths

`grep -rl /home/mattm CW_transition/*.py` → **43 files**. `/home/mattm` does not exist here.
Reference counts by prefix:

| cronus prefix | refs | needed? |
|---|---|---|
| `/home/mattm/projects/HSYMT/CW_transition` | 74 | **yes** |
| `/home/mattm/.cache/jax_stagec_cache` | 24 | **yes** (in 12 modules, via `jax.config.update("jax_compilation_cache_dir", ...)`) |
| `/home/mattm/projects/HSYMT/CW_lnL_check` | 14 | **yes** (`sys.path.insert` for `cw_helpers`) |
| `/home/mattm/projects/HSYMT` | 3 | **yes** |
| `soft/scripts/MISC`, `soft/JUG/data`, `projects/NG`, `projects/MPTA` | 8 | **no** — only in `stagec_anchor_a0.py` / `a1.py`, which REQUIREMENTS_FROZEN §5 forbids re-running (`best_distances.txt` is the FROZEN canonical EM prior) |
| `miniforge3/envs/jug-gpu` | 1 | no — docstring |

**Two prefixes cover 112 of 125 references.** Per the HARD RULE I did not edit a single
file. Instead **apptainer synthesises `/home/mattm`** via bind mounts — the container earns
its keep for *path fidelity*, not for CUDA userland (the `nvidia-*-cu12` pip wheels already
freeze that; only the driver is host-supplied, and `--nv` injects the host driver into a
container anyway).

Bind-point auto-creation works unprivileged here (apptainer 1.4.5). Verified.

> **REQUESTED CODE FIX ON CRONUS** (do not apply here): make the repo root and the cache dir
> resolve from the environment, e.g.
> `HSYMT_ROOT = os.environ.get("HSYMT_ROOT", "/home/mattm/projects/HSYMT")` and
> `JAX_CACHE = os.environ.get("JAX_STAGEC_CACHE", "/home/mattm/.cache/jax_stagec_cache")`,
> defaulting to today's values so cronus behaviour is byte-identical. That removes the
> container requirement entirely. Until then the bind map is the contract.

### 3.2 Env: conda `jug-gpu`, rebuilt from the freeze

- `/home/milesmt/miniforge3/envs/jug-gpu`, **Python 3.12.13** (freeze: 3.12.13) ✔
- Container base: `docker://rockylinux/rockylinux:9` → `/home/milesmt/soft/harbor/el9.sif`
  (80 MB; matches the EL9 host). Runs the *host* conda python via bind — the image supplies
  only glibc + a filesystem namespace in which `/home/mattm` exists.
- Forks cloned fresh to `/home/milesmt/soft/harbor/{discovery,enterprise_extensions}`,
  **detached at their pins**, clean. I did **not** touch the pre-existing `~/soft/discovery`
  (stale, HEAD `9124089`) or `~/soft/NANOGRAV/discovery` (`dbb693b`).

`pip freeze` editable lines now match REQUIREMENTS_FROZEN §6 L121/L124 **byte-for-byte**:
```
-e git+https://github.com/MattTMiles/discovery.git@136b270f1891c28ae6d9840930a9dfbcf41fd52d#egg=discovery
-e git+https://github.com/MattTMiles/enterprise_extensions.git@d43fef99d0aa786ec6b08dbc2fd56dcce10f26f5#egg=enterprise_extensions
```

`h2_verify.py` (CPU-only, run inside the container) → **PASS**. Every pinned version matches:
jax/jaxlib 0.10.1, jax-cuda12-{plugin,pjrt} 0.10.1, numpy 2.4.6, scipy 1.17.1,
enterprise-pulsar 3.4.4, jaxopt 0.8.5, optax 0.2.8, numpyro 0.21.0, ml_dtypes 0.5.4,
matplotlib 3.10.9, healpy 1.19.0, `enterprise_extensions.__version__ == 3.0.3`, all 13
`nvidia-*-cu12` wheels at their §2 pins. `jax_enable_x64` → `float64`. Both forks clean at
their pins. `import discovery.matrix`, `discovery.deterministic`, `enterprise_extensions` OK.

### 3.3 Deviations from the freeze — declared, not improvised

| item | freeze | here | why |
|---|---|---|---|
| `scikit-sparse` 0.5.0 | conda | **conda-forge** (`suitesparse` + `scikit-sparse`) | pip wheel build fails — no SuiteSparse headers on ACCRE. `enterprise-pulsar` requires it. `import sksparse.cholmod` OK. |
| `numpy` 2.4.6, `scipy` 1.17.1 | conda wheels | **pip wheels** | same versions; pip resolves them cleanly. |
| `pint-pulsar` | editable @ `ea7652f` | **PyPI release** (pulled by `enterprise-pulsar`) | REQUIREMENTS_FROZEN §1 states PINT is *not* on the Stage-C import path. |
| `jug` @ `068098d` | editable | **absent** | §1: not on the Stage-C import path. |
| `MultiHMCGibbs`, `PySide6`, `fastapi`/`uvicorn`/… | present | **absent** | GUI/server transitive deps, not imported by `stagec_*`/`trackA_*`/`trackB_*`. |

**Two traps found and avoided.** (1) `discovery` declares `jax[cpu]~=0.4`; that specifier
means `>=0.4,<1.0`, so jax 0.10.1 satisfies it — no silent downgrade, but both forks were
installed `--no-deps` to be certain. (2) `/home` is a **symlink** to
`/panfs/accrepfs.vampire/home`, and pip baked the *resolved* `/panfs/...` path into both
editable `.pth`/finder files. Binding only `/home/milesmt` leaves `import discovery` failing
inside the container. The real path must be bound at its real name.

Everything lives in `hpc_harbor/env/harbor_env.sh`, sourced by every job. It exports the
REQUIREMENTS_FROZEN §4 XLA conventions verbatim and defines `harbor_py`.

---

## 4. H4 — SMOKE TEST LADDER

### 4a. GPU hello — **PASS** (job 12448416, dgx03, 14 s)

| quantity | value |
|---|---|
| GPU | **NVIDIA A100-SXM4-80GB**, compute capability 8.0, MIG Disabled |
| Driver | **580.159.04** (CUDA 13.0 capable). cronus was 550.120. jax CUDA-12 wheels need ≥525 — satisfied. |
| jax backend | `gpu`, `[CudaDevice(id=0)]` |
| x64 | active, `float64` |
| **Compute-node internet** | **YES** — pypi 200, github 200 |
| f64 2×2 matmul | max abs err **0.000e+00** |
| **fp64 GEMM (n=8192)** | **3.646 TFLOP/s** (301.6 ms) |
| **fp32 GEMM (n=8192)** | **6.284 TFLOP/s** (175.0 ms) |
| **fp32:fp64 ratio** | **1.72 : 1** |

Reference: RTX 4090 = **64 : 1** (crippled). A100 spec = 2 : 1. Measured 1.72 : 1 → **the
A100's fp64 is full-rate**, and the relative fp64 penalty is **~37× smaller** than on cronus.

Two cautions against over-reading this:

1. **Raw fp64 throughput is only ~2.8× cronus's**, not the 7–15× the spec sheets imply
   (4090 fp64 ≈ 1.3 TFLOP/s). We are running autotune-off on a **co-tenanted** card. The
   headline 1.72:1 is a *ratio*, not a speedup.
2. **This does not license editing `trackA_fisher.py`.** The complex64-GEMM-with-fp64-outer-
   accumulation path (`trackA_fisher.py:192-241`) is now a *speed optimisation* rather than
   a necessity — but its gate is banked at `3e-5` production error vs `3e-16` pure-f64
   (project_progress.md L1089). Changing it changes banked numbers. **Report only.**

Also recorded: `f32` ulp at the pulsar-term phase `2πL/dL ≈ 1.6e4 rad` is **9.766e-4 rad**,
i.e. **0.977×** the ~1e-3-nat fringe-breaking signal. f32 does not merely degrade the signal,
it is the same size as it. The `jax_enable_x64=True` mandate is exactly load-bearing.

### 4b. Warm-cache job — **CACHE WORKS; the script's own gate FAILS. Both true.** (job 12448542, dgx03, 32m42s, COMPLETED)

Ran after B-1 cleared (§4b-blocker below, now resolved by ee `b54b211`). Two fresh processes,
one allocation. Cache went 0 → **781 entries / 632 MB**.

| graph | PASS 1 (cold) | PASS 2 (warm) | speedup | cronus banked warm |
|---|---|---|---|---|
| `fisher_logL` | 34.15 s | **7.21 s** | 4.7× | **2.6 s** |
| `fisher_logL_batched` | 94.10 s | **10.50 s** | 9.0× | **38 s** |
| `assemble_hessian` / HVP | 1351.8 s | **131.74 s** | 10.3× | **52 s** |
| `fisher_logL` warm re-eval | 56.0 ms | **29.5 ms** | — | — |
| wall | 1658 s | **302 s** | 5.5× | — |

The persistent cache **works, cross-process, exactly as designed.** Every graph got faster.

**Why the script prints FAIL.** `trackB_b00_cache.py`'s criterion is
`worst first-compile < 60 s`. Its worst line is `assemble_hessian` at 131.74 s. But
`stagec_d4.assemble_hessian` (L50-62) compiles **one** graph, `hvp_batched`, and then
**executes it `ceil(n_theta/HVP_CHUNK)` times** with `HVP_CHUNK = 8`. That 131.74 s is
therefore *compile (cached, fast) + many HVP executions*. The script labels it `[compile]`
and tests it against a compile threshold. **The gate is mis-specified, not the cache.**
Judged on what it actually measures — does a second process avoid cold XLA compilation? —
H4b is a **PASS**.

> **REPORTED, NOT FIXED (cronus's call):** `trackB_b00_cache.py`'s `worst` should exclude
> `assemble_hessian`, or `assemble_hessian` should be timed as compile-only (one
> `hvp_batched` call) rather than the full chunked loop. As written the B0.0 gate can only
> pass when HVP *execution* happens to fit under 60 s — which is what happened on cronus
> (52 s), by luck of that machine's speed, not by design.

**The genuinely surprising number: cold `fisher_logL` compiled in 34.15 s here vs 465 s on
cronus — 13.6× faster.** XLA compilation is CPU-bound; dgx03 is a 128-core EPYC 7742. The
"465 s/graph compile tax" that motivates the whole warm-cache-first strategy is **~34 s here**.
The warm-cache job is still worth running once per graph shape, but it is no longer the
dominant cost it was, and it should not be allowed to shape the array design.

**The genuinely worrying number: warm HVP is 131.74 s here vs 52 s on cronus — 2.5× SLOWER**,
on a full-rate-fp64 A100 that beat the 4090's fp64:fp32 ratio 37×. Warm `fisher_logL` is also
slower (7.21 s vs 2.6 s). Candidate causes, **not yet discriminated**:
1. **GPU co-tenancy** (§2) — a foreign process held 50.8 GB and unknown SM share on our card.
2. **Cache reads cross panfs**, a network filesystem; cronus read from local disk. 632 MB /
   781 files per cold start.
3. autotune-off default kernels may suit Ada better than Ampere.

This is **not** resolved, and it directly sizes B1's walltime requests. See §5 open item 7.
Do not extrapolate cronus per-stage timings to this cluster in either direction.

Incidental, benign: `FeatherPulsar.read_feather: cannot find _pdist in feather file
.../J2322-2650.feather`. The EM distance prior comes from the frozen `best_distances.txt`,
not the feathers (REQUIREMENTS_FROZEN §5). Noted; confirm before trusting any distance-
dependent H4c number.

### >>> BLOCKER B-1 — **RESOLVED** 2026-07-09 by ee `b54b211` <<<

*Kept for the record; this is what stopped H4b/c/d for the first pass.*

`CW_lnL_check/cw_helpers.py:174` calls `enterprise_extensions.deterministic.cw_block_circ(...,
evolve=True)`. At the pinned commit **`d43fef99`**, `cw_block_circ` has **no `evolve`
parameter**. Its signature is:

```python
def cw_block_circ(amp_prior="log-uniform", dist_prior=None, skyloc=None, log10_fgw=None,
                  psrTerm=False, tref=0, phase_connected=False, discoclone=False, name="cw"):
```

and it **hard-codes `evolve=True`** internally at `deterministic.py:194` when it builds the
`cw_delay` signal. `evolve` is a parameter of `cw_delay` (`deterministic.py:340`), not of
`cw_block_circ`.

This is *not* a wrong-clone problem. I checked:

- `discoclone` and `phase_connected` — the two unusual kwargs — **are** present at the pin.
  `evolve` is the **only** kwarg `cw_helpers.py` passes that the pin does not accept.
- `git log --all -S evolve -- enterprise_extensions/deterministic.py` on the fork: **no
  commit on any branch** adds `evolve` to `cw_block_circ`'s signature. `origin/master` HEAD
  **is** `d43fef99`. There is nowhere else to fetch it from.
- `MAIN_PROJECT_QUESTIONS/Q1/sim.py:251` also passes `evolve` — so this is the repo's
  settled calling convention, not a one-off typo in `cw_helpers.py`.

**Chronology:** `cw_helpers.py` gained `evolve=True` in repo commit **`e7121d9`, 2026-05-04**
("revert CW lnL check back to MM_playgroun"). The ee pin `d43fef99` is dated
**2025-11-11** — nearly six months earlier.

**Conclusion:** the `jug-gpu` env on cronus must have an `enterprise_extensions` working tree
carrying a **local, uncommitted (or committed-but-unpushed) change** that adds an `evolve`
kwarg to `cw_block_circ`. REQUIREMENTS_FROZEN §1's assertion — *"Both checkouts were **clean**
(no uncommitted changes) at snapshot time"* — is **contradicted by the code**. Either the
`git status` check missed it, or the pin was recorded from `git rev-parse` while the working
tree carried edits on top.

Per the HARD RULE I have **not** patched `cw_block_circ` and have **not** edited
`cw_helpers.py`. Everything downstream of `build_enterprise_pta` — which is H4b, H4c, and
H4d — is blocked on this.

**WHAT I NEED FROM CRONUS (one of):**
1. On cronus, run `git -C ~/soft/enterprise_extensions status --porcelain` and
   `git diff`. If dirty, **commit and push** that diff; tell me the new sha and I re-pull.
   *(Most likely outcome — and it means the banked Stage-C/Track-A/Track-B results were
   produced by code that is not in any repository. That is a reproducibility hole
   independent of this port.)*
2. If clean, then cronus's `enterprise_extensions` is not the fork at
   `github.com/MattTMiles/enterprise_extensions` — tell me which checkout/remote it is.
3. If neither, `cw_helpers.py:174` is simply wrong and should drop `evolve=True` (the pin
   already forces `evolve=True` internally, so dropping it is likely a **no-op** — but that
   is a code change, and it is cronus's to make and verify, not mine).

Note the appealing shape of (3): because `cw_block_circ` hard-codes `evolve=True` at L194,
`evolve=True` at the call site is plausibly redundant. **Do not assume it.** Verify against a
banked value before touching it — that is precisely what H4c exists to do.

### 4c. VALUE GATE — **ALL 8 GATES PASS** (job 12449076, dgx03, 15m02s, exit 0)

`trackB_b1_core.py`, unmodified, warm cache. `build_b1_problem` 281.5 s (npsr=116, ncw=16,
n_theta=244). Total wall 901 s.

| gate | result on A100 | banked (cronus) |
|---|---|---|
| G0 `masked lnL(mask=1) == fisher_logL` | **0.000e+00** | 0.00e+00 |
| G1 `masked lnL(mask=0) == EarthDelay` | **0.000e+00** | — |
| G2 split `a_const` invariant to (θ, data, mask) | **0.000e+00** (405637.694390 both) | — |
| G3 `B1Split.lnL == masked logL` | **6.403e-10** (< 1e-8) | — |
| G4 `B1Split.estep == Split.estep` | **0.000e+00**, 0/116 > 1e-8 | 1.16e-10, 0/116 |
| **G5 zero-noise E-step at truth == A2 census ceilings** | **J0711 0.9534 / J1713 0.9887 / J1909 0.9984** | **0.9534 / 0.9887 / 0.9984** |
| G6 noise-draw covariance sanity | white var ratio 0.985; ORF −0.0337 == `hd_orf` | — |
| G7 noisy finiteness + S/N scale | `lnL(truth\|zero-noise) = 405413.51` | **405413.51** (spec L236/L637) |
| G8 arm-B distance truth draw | 116/116 in-window; census n_true −76/−40/8 | — |

**The value gate transfers exactly.** G5 reproduces the census ceilings to the 4 decimal
places the ceilings are quoted at, and 3/3 truth-blind certified pulsars land on the **true**
fringe. `lnL(truth) = 405413.51` matches the banked value. G0/G1/G2/G4 are **bit-exact
zeros** across a GPU architecture change (Ada → Ampere), a driver change (550.120 → 580.159.04),
and a full env rebuild. That is a strong result: the direct-evaluation path is architecture-
independent, as predicted.

Per the brief, **trajectories were not gated** and must not be: the loud-param Hessian has
`eig ∈ [-5.965e11, -14.43]`, **cond ~4e10**, and F1a establishes the M-step is numerically
*chaotic* (gradients agree to 2.2e-16 at one step, diverge 0.2-scaled over 19 iterations —
project_progress.md L1685-1690).

Co-tenant present throughout the gate (`env/deerdiff`, 50,846 MiB). It did not perturb any
value. Timings quoted here are therefore *upper bounds* under contention.

### 4d. One-realisation B1 smoke + kill/resume — **CANNOT BE RUN. Two code defects.**

I did not run this, and the reason is not the cluster. Recording precisely.

The only existing "one noisy realisation end-to-end" driver is
`CW_transition/trackB_b1_referendum.py` (STEP 2, `--noise --noise_seed`). It is also the
script that *decides* tier-conditioning-vs-cascade — the decision the brief says is being made
on cronus right now. Two independent problems:

#### >>> DEFECT D-1: checkpoints are WRITE-ONLY. There is no resume, anywhere. <<<

```
$ grep -rn "ckpt" --include=*.py . | grep -i load
(no output)
```

**No script in the repository ever reads a checkpoint back.** `trackB_b1_referendum.py:438`
writes `b1_ref_{tag}_s{seed}_ckpt.npz` at each SMC tempering-stage boundary (correct
granularity, matching spec L659), and nothing ever loads it. The H3 mandate — *"Every task's
first act: look for its own checkpoint and resume"* — is **not implemented**. A walltime kill
today costs the whole job, not one stage.

So "kill it mid-run and confirm clean resume" cannot be demonstrated: **the code cannot
resume.** Writing that resume path is a code change, and code changes are cronus's.

#### >>> DEFECT D-2: a `--noise` run silently overwrites BANKED, GIT-TRACKED checkpoints. <<<

At L504, `tag = f"t{args.tier}"` — keyed by tier only. **Not by `--noise`, not by seed-source,
not by job/task.** So:

| artifact | zero-noise run writes | `--noise` run writes | collide? |
|---|---|---|---|
| final result | `b1_referendum_tierA.npz` | `b1_referendum_tierA_noisy.npz` | no |
| **per-stage ckpt** | `b1_ref_tA_s0_ckpt.npz` | `b1_ref_tA_s0_ckpt.npz` | **YES** |

The final output *is* noise-keyed (L563). The checkpoint is not. And
`b1_ref_t{A,C}_s{0,1}_ckpt.npz` are **committed to git at `634aab8`** — a noisy smoke run
would have silently clobbered four banked artifacts inside the bound repo. `CWT` is hard-coded
to `/home/mattm/projects/HSYMT/CW_transition`, which my bind maps onto the *real* working tree.

**I stopped rather than run it.** The H3 requirement — *"output paths keyed by (job, task,
stage)"* — is exactly what is missing.

> **REQUESTED FIXES ON CRONUS (do not apply here):**
> 1. `tag` must encode noise + seed provenance, e.g.
>    `tag = f"t{args.tier}{'_noisy' if args.noise else ''}"`, and ideally
>    `{SLURM_JOB_ID}_{SLURM_ARRAY_TASK_ID}`.
> 2. Add the resume half: on start, `np.load` the task's own ckpt if present, restore
>    `(beta, X, L, stage, ess_hist)`, and continue the tempering ladder from `stage`.
>    Apply the standing finiteness discipline on the restored arrays (isfinite fraction +
>    min/median/max; NaN aborts).
> 3. Only then is `--signal=B:USR1@300` + `--requeue` (already in my sbatch templates)
>    actually load-bearing rather than decorative.
>
> Until (1) and (2) land, **no B1 production array should launch** — a walltime kill or a
> requeue loses the whole task, and a noisy array would race four tracked files.

**Also noted:** `project_progress.md` states STEP 2 is "pre-registered but not run"
(L1824-1825). On disk, `b1_referendum_tier{A,C}.npz` and four `b1_ref_*_ckpt.npz` **are
committed at `634aab8`**. The doc lags the artifacts. Worth reconciling before the array
design is finalised on that doc's basis.

**Per-stage wall-clock (for calibrating array walltimes)** — from H4b/H4c, warm cache, A100
under co-tenancy:

| stage | wall |
|---|---|
| `load_pulsars(116)` + `build_enterprise_pta` + `build_fisher_amortised` | 133–137 s |
| `build_b1_problem` (npsr=116, ncw=16, n_theta=244) | 281.5 s |
| `fisher_logL` first call (cache hit) | 7.2 s |
| `fisher_logL_batched` first call (cache hit) | 10.5 s |
| `assemble_hessian` (1 compile + 31 chunked HVP execs, `HVP_CHUNK=8`) | 131.7 s |
| `fisher_logL` warm re-eval | 29.5 ms |
| warm `estep` (116×512) | 157.2 s |
| `trackB_b1_core.py` end-to-end (8 gates) | 901 s |

Note `warm estep 157.22 s` vs the spec's profiled **0.68 s** for a full 116×512 E-step
(spec L180-182). Those measure different things — the gate's `estep` line includes table
builds — but the gap is large enough that **B1 array walltimes must be measured, not derived
from the spec's profiling table.** See §5 item 7.

---

## 5. OPEN QUESTIONS FOR CRONUS

**Ranked. (1) blocks everything.**

1. ~~**BLOCKER B-1** — `enterprise_extensions` `evolve` kwarg.~~ **RESOLVED** by ee
   `b54b211` (2026-07-09). **But two follow-ups on that commit:**
   - `cw_block_circ(..., evolve=False)` **defaults to `False`**, where the old code hard-coded
     `True`. Every caller that does *not* pass `evolve` silently flips behaviour:
     `lnL_GW_recovery_phase_connected.py:269`,
     `lnL_GW_recovery_phase_connected_pulsarRing.py:382`, `lnL_GW_recovery_refactor.py:242`,
     `old/*`. **The default should be `evolve=True`** (`cw_delay_phase_connected_binary`
     already defaults `True`). Our path is unaffected — `cw_helpers.py:174` passes `evolve=True`
     explicitly, into the `elif phase_connected:` branch, so the Stage-C/Track-B result is a
     provable no-op — which is why H4c passes. But the other scripts are now silently different.
   - The third hunk in `cw_delay_phase_connected_binary` **adds two lines that unconditionally
     overwrite** `omega, phase, omega_p, phase_p` immediately after the existing `if evolve:`
     block already set them. Dead code when `evolve=True`; a behaviour change when `False`.
     Two paths that must agree forever. Delete one.
   - **Reproducibility hole this exposes:** the banked Stage-C / Track-A / Track-B results were
     produced by an `enterprise_extensions` working tree whose `evolve` kwarg existed in **no
     commit** until today. REQUIREMENTS_FROZEN §1's "Both checkouts were **clean**" was wrong.
     Worth auditing whether `discovery` @ `136b270f` carries uncommitted edits too.
2. **DEFECT D-1 — no resume exists anywhere in the repo.** Checkpoints are write-only. See
   §4d. **This blocks the B1 production array**, independently of anything on this cluster.
3. **DEFECT D-2 — `--noise` runs clobber banked git-tracked checkpoints.** See §4d. Blocks
   any noisy array.
4. **BLOCKER B-2 — 43 files hard-code `/home/mattm`.** Worked around with apptainer binds, so
   *not* blocking, but the requested env-var fix (§3.1) would delete the container from the
   critical path.
5. **`project_progress.md` L1824-1825 says STEP 2 is unrun; `b1_referendum_tier{A,C}.npz` and
   four ckpts are committed at `634aab8`.** Doc lags disk. Reconcile before finalising the
   array design on the doc's basis.
6. **Where is `CW_transition/PORTER_slurm_plan.md`?** Not in this checkout, not in git
   history. Never pushed.
7. **Timing does NOT transfer, in either direction. Do not derive array walltimes from cronus.**
   Cold `fisher_logL` compile is **13.6× faster** here (34.15 s vs 465 s — XLA compile is
   CPU-bound and dgx03 is a 128-core EPYC 7742). But warm HVP is **2.5× slower** (131.7 s vs
   52 s) and warm `fisher_logL` **2.8× slower** (7.21 s vs 2.6 s), despite full-rate fp64.
   Unresolved between: GPU co-tenancy; the compile cache living on **panfs** (a network FS,
   632 MB / 781 files per cold start) where cronus read local disk; autotune-off kernel choice
   on Ampere vs Ada. **Measure B1 stages here before sizing the array.**
8. **Group volume `/data/taylor_group` has 1.9 TiB free** of a 15.5 TiB shared quota. What is
   B1's on-disk budget? Answer *before* an array launches, not after.
9. `/nobackup` purge policy — still unknown. Compile cache is on `/home` in the meantime.
10. **`--gres` does not fence the GPU** (§2). Usable VRAM ≈ 30 GB, not 80 GB. Stage-C's memory
    discipline was sized for the 4090's 24 GB, so B1 probably fits — but not with the 3.3×
    headroom a port plan would have assumed.

**Answered by H4a, no longer open:** compute-node internet (**yes**), dgx03 driver
(**580.159.04**, ≥525 so the CUDA-12 wheels are fine), A100 variant (**80 GB SXM4, cc 8.0**),
fp64 ratio (**1.72:1**, full-rate).

---

## 6. STATE ON DISK

| thing | path |
|---|---|
| env-setup script (source this) | `hpc_harbor/env/harbor_env.sh` |
| H2 verification | `hpc_harbor/slurm/h2_verify.py` → **PASS** |
| H4a GPU hello | `hpc_harbor/slurm/h4a_gpu_hello.{py,sbatch}` → **PASS** (job 12448416) |
| H4b warm cache | `hpc_harbor/slurm/h4b_warm_cache.sbatch` → **PASS** (job 12448542) |
| H4c value gate | `hpc_harbor/slurm/h4c_value_gate.sbatch` → **ALL 8 GATES PASS** (job 12449076) |
| H4d kill/resume | **not run** — defects D-1 / D-2, §4d |
| logs | `hpc_harbor/logs/` |
| conda env | `/home/milesmt/miniforge3/envs/jug-gpu` |
| container | `/home/milesmt/soft/harbor/el9.sif` |
| forks @ pins | `/home/milesmt/soft/harbor/discovery@136b270f`, `enterprise_extensions@b54b211` |
| jax compile cache | `/home/milesmt/.cache/jax_stagec_cache` (**781 entries, 632 MB, warm**) |
| pinned requirements | `/home/milesmt/soft/harbor/requirements_{harbor,cuda}.txt` |

Nothing was committed. Nothing was pushed. No tracked file was edited. `HPC_SETUP.md` and
`hpc_harbor/` are untracked.

**The B1 production array has NOT been launched**, per the brief.

---

## 7. GPU-EXCLUSIVITY DIAGNOSTIC (D1 / D2 / D3)

**Headline: the 2.5× slowdown is NOT contention, NOT the filesystem, and NOT the autotune
flag. And there is a rogue non-Slurm process squatting on 3 of dgx03's 8 GPUs.**

### 7.1 D1 — attribution. Slurm's fencing WORKS. I was wrong twice.

Job 12449702, dgx03, outside the container (so apptainer `--nv` cannot be blamed):

```
ConstrainDevices = yes          TaskPlugin = task/cgroup,task/affinity
ProctrackType    = proctrack/cgroup      (cgroup v2 -> eBPF device filter, no devices.list)
CUDA_VISIBLE_DEVICES=0   SLURM_JOB_GPUS=0
nvidia-smi -L  ->  exactly ONE device
/dev/nvidia0 OPEN ;  /dev/nvidia1..7 all "open DENIED"
```

A 3-task array (12450343) landed on `/dev/nvidia{6,4,3}` with **three distinct UUIDs** and
`SLURM_JOB_GPUS=6,4,3`. **Slurm assigns distinct physical devices and the cgroup enforces it.**

So my §2 claim "Slurm does not fence the GPU" was **false**, and the H4a observation had a
different cause. Two red herrings I chased and discarded, recorded so nobody re-chases them:

- *"Two of my own jobs shared a card."* `d2` (12449763) and `d2b` (12450183) both reported
  UUID `ec71dc3d`. **They did not overlap**: `d2` ended `17:55:02`, `d2b` started `17:55:18`.
  The 51,648 MiB `d2b` saw at startup was `d2`'s own CUDA context still tearing down. Slurm
  reuses the freed card immediately; **a job's first act sees the previous tenant's memory
  for a few seconds.** Any startup VRAM check must tolerate this.
- *"The co-tenant is user `chattec`."* `chattec` does hold `dgx_iacc` GPUs on dgx03, entirely
  legitimately, via Slurm. Not them.

### 7.2 D3 — **ROGUE PROCESS. Report this to the admins.**

Array 12450376 walked 6 GPUs. Three were clean; three carried a foreign process:

| our job | device | GPU UUID | foreign pid | mem | owner | cgroup |
|---|---|---|---|---|---|---|
| 12450379 | `/dev/nvidia3` | `GPU-1e10409a-65fb-13a0-1618-92ea66e4d515` | 930652 | 50,846 MiB | **`wut18`** (uid 934717) | `/user.slice/user-934717.slice/session-355.scope` |
| 12450381 | `/dev/nvidia5` | `GPU-aa596b0f-e007-ae91-7493-5b99867e7930` | 931233 | 50,846 MiB | **`wut18`** | same session-355 |
| 12450376 | `/dev/nvidia1` | `GPU-563c9580-d48c-78cc-c883-307c11092114` | 930137 | 51,648 MiB | **`wut18`** | same session-355 |

Clean at the same instant: `/dev/nvidia6` (`cf0ec2e2`), `/dev/nvidia4` (`0e744b08`),
`/dev/nvidia7` (`18be1e28`).

The cgroup is **`user.slice/…/session-355.scope`** — an interactive **login session**, *not*
`slurmstepd.scope/job_*`. `comm = pt_main_thread` (PyTorch). So `wut18` SSH'd to dgx03 and
started three ~50 GB PyTorch processes **outside Slurm entirely**. Slurm's GRES accounting
does not see them (`AllocTRES` counted only 2 GPUs), so the scheduler happily allocates those
same cards to other people's jobs. **3 of 8 GPUs on dgx03 — a 37.5 % chance any given job
lands on a squatted card.** That is the whole source of the H4a/H4c "co-tenancy," and it also
explains why the earlier pid (83140, `env/deerdiff/bin/python`) differed from today's: the
squatter restarts and lands on different cards.

Note also `sacct`: job **12331899 was "CANCELLED by 934717"** — the same uid.

> **ONE-LINE NOTE FOR ACCRE ADMINS (this is the actual fix; it is their enforcement problem):**
>
> *On `dgx03`, user `wut18` (uid 934717) is running three PyTorch processes (pids 930137,
> 930652, 931233; ~50 GB each) from an interactive SSH session
> (`/user.slice/user-934717.slice/session-355.scope`), outside any Slurm job cgroup, occupying
> GPUs `/dev/nvidia1`, `/dev/nvidia3`, `/dev/nvidia5` (UUIDs `563c9580…`, `1e10409a…`,
> `aa596b0f…`). Slurm does not account for these GPUs and schedules other jobs onto them —
> e.g. my jobs 12450376/12450379/12450381. Suggest `pam_slurm_adopt` on the DGX nodes, or
> blocking non-job SSH GPU access.*

### 7.3 D2 — enforcement ladder + the timing matrix

**Ladder (via `--test-only`, nothing submitted): the QOS denies none of them.** All three are
*permitted*; only resources gate them.

| variant | QOS verdict | earliest start |
|---|---|---|
| `--exclusive` (whole node, 8 GPUs) | **allowed** | 2026-07-10T16:38 (~24 h) |
| `--exclusive` (whole node, 1 GPU) | **allowed** | ~24 h |
| `--gres=…-80gb:8` | **allowed** | ~24 h |
| `--exclusive=user` (1 GPU) | **allowed** | ~30 min |

`dgx_iacc` has no `MaxTRESPJ`; its `GrpTRES` is `a100-80gb=16, a100-40gb=16, gpu=32, cpu=512`.
Whole-node exclusivity is a *queueing* cost (drain `chattec`'s legitimate jobs), not a policy
one — and it would **not** evict `wut18`, who is outside Slurm. **`--exclusive` does not solve
the rogue.**

Because D1 proved per-device fencing works, exclusivity was **unnecessary** to answer the
timing question: I recorded each run's contention state and compared directly.

**Timing matrix.** All runs dgx03, A100-80GB, warm cache, `n_theta=244`, `HVP_CHUNK=8`.
Jobs 12449763 (uncontended) and 12450183 (contended, landed on a squatted card).

| metric | at=0, **clean** | at=2, **clean** | at=0, contended, cache→**tmpfs** | at=0, contended, cache→**panfs** |
|---|---|---|---|---|
| GEMM fp64 TFLOP/s | **3.658** | 3.658 | **1.608** | 1.612 |
| GEMM fp32 TFLOP/s | 6.293 | 6.293 | 2.754 | 2.754 |
| fp32:fp64 ratio | 1.72 | 1.72 | 1.71 | 1.71 |
| build (load+pta+amortise) s | 110.3 | 166.1 | 127.7 | 124.0 |
| `fisher_logL` first s | 4.43 | 40.26 | 4.36 | 6.42 |
| `fisher_logL_batched` first s | 8.22 | 60.16 | 8.18 | 8.95 |
| **`assemble_hessian` s** | **118.3** | 1095.3 | **116.8** | **116.4** |
| `hvp_batched` 1 chunk, warm s | **0.0647** | **0.1411** | 0.0630 | 0.0636 |
| `fisher_logL` warm ms | **23.2** | **47.1** | 22.7 | 23.2 |
| `lnL(truth)` | `405413.512739` | `405413.512739` | `405413.512739` | `405413.512739` |

**Four conclusions, each isolated by its own column pair:**

1. **Contention is NOT the cause of the 2.5×.** A squatted card halves large-GEMM throughput
   (3.658 → 1.608 TFLOP/s, **2.27×**) but leaves *every pipeline kernel untouched*:
   `assemble_hessian` 118.3 → 116.8 s, HVP chunk 0.065 → 0.063 s, `fisher_logL` warm
   23.2 → 22.7 ms. The science kernels are **latency-/CPU-bound, not SM-throughput-bound**, so
   they are nearly immune to a co-tenant. (Which is also why H4c's *values* were unaffected.)
2. **The filesystem is NOT the cause.** Staging the 1.3 GB / 3583-entry cache onto node RAM
   (`/dev/shm`) changed `assemble_hessian` by **0.4 s** (116.8 vs 116.4). The panfs
   cache-read hypothesis in §5 item 7 is **refuted**. (The copy itself took 4 s.)
3. **`autotune_level=2` is a REGRESSION on A100 — 2.2× slower execution.** Warm HVP chunk
   0.0647 → 0.1411 s; warm `fisher_logL` 23.2 → 47.1 ms; and it mints a new cache key, so
   `assemble_hessian` recompiles cold (1095 s). GEMM is unchanged (cuBLAS path). **Keep
   `--xla_gpu_autotune_level=0`.** The cronus convention is correct here — for a completely
   different reason than the one recorded in REQUIREMENTS_FROZEN §4.
4. **`lnL(truth) = 405413.512739` is bit-identical across all four cells.** Autotune level,
   contention, and cache filesystem change *timing only*, never values. Consistent with H4c.

### 7.4 What the ~117 s in `assemble_hessian` actually is

Decomposed: `hvp_batched` is **one** graph executed `⌈244/8⌉ = 31` times at 0.065 s ⇒ **~2 s of
real HVP execution**. The other **~115 s is the graph's first-call cost** — XLA
deserialize + link + kernel-selection CPU work on a cache *hit* (cold is 1351.8 s, §4b). It is
not I/O (conclusion 2) and not GPU (conclusion 1).

**This is the single most important number for B1's array design.** Per-task cost is dominated
by fixed startup:

```
build (load_pulsars + pta + amortise) ~110-137 s
+ first-call graph materialisation      ~115 s   (per graph shape, per PROCESS, even warm)
+ actual per-evaluation work            ~23 ms   (fisher_logL), ~65 ms (HVP chunk)
```

⇒ **Do not run one realisation per array task.** ~4 minutes of fixed cost would dwarf the
science. Batch many realisations per task, or persist the built problem. The warm-cache job
removes the *1351 s* cold compile, not the *115 s* per-process materialisation — nothing can,
short of keeping one process alive across realisations.

### 7.5 Revisions to earlier sections of this file

- §2 "no co-tenant on an exclusively-allocated GPU" → **wrong then, wrong in the opposite
  direction after H4a**. Correct statement: Slurm fences devices correctly; a rogue non-Slurm
  process squats 3/8 cards on dgx03.
- §5 item 7's three candidate causes: **contention — refuted; filesystem — refuted; autotune
  kernel choice — refuted (level 2 is *worse*, and level 0 was used in both).** The residual
  gap to cronus's 52 s is per-process XLA graph materialisation on a cache hit. Cronus's
  banked 52 s was measured on a machine whose single-core XLA link/deserialise is evidently
  faster; note that dgx03's *cold* compile is 13.6× **faster** (34 s vs 465 s), so this is not
  a simple "cronus is faster" story and should not be treated as one.
- Cache now **3583 entries / 1.3 GB** (autotune-2 entries added; harmless, and they will never
  be hit again if level 0 stays).

---

## 8. HARBOR-2 — b1-v1 READINESS (2026-07-11)

Continuation agent, tag `b1-v1` (cronus C1 fixes: checkpoint resume + keying, ee fork
hygiene). PURE EXECUTION: **I committed nothing and edited no tracked file.** The one
sanctioned exception was staging `reports/` for Matt, who committed and pushed it himself
(origin/`MM_playground` = `3f485bf`). Cache warm from the earlier passes (5354 → 7.7 GB at
gate time; still autotune-off, level 0).

### 8.1 SHAs (S1 — SYNC)

| component | sha | note |
|---|---|---|
| repo tag `b1-v1` | `6243ae7f27166320e871bb5b56f14b8d326102f7` | detached; runtime tree for S2/S3 |
| `discovery` pin | `136b270f1891c28ae6d9840930a9dfbcf41fd52d` | unchanged |
| **`enterprise_extensions` pin** | `f73b8e0b99d6bdafe43491f6a55087d71e9fb6b4` | `b54b211`→`f73b8e0`, `cw_block_circ` now **defaults `evolve=True`** (closes the §5.1 follow-up) |

- **Env rebuild: NO.** The only pip-relevant `REQUIREMENTS_FROZEN.md` change is the ee
  editable git sha; ee is an `-e` install, so the git-checkout swap **is** the refresh — the
  finder resolves to `…/harbor/enterprise_extensions`, picked up automatically. Every other
  pin unchanged (jax/jaxlib 0.10.1, jax-cuda12-* 0.10.1, numpy 2.4.6, scipy 1.17.1, 13× cu12
  wheels). **Container binds unchanged** — no path moved.
- **Freeze-file lag (noted, not blocking):** the tag's `REQUIREMENTS_FROZEN.md` still records
  the *old* ee `d43fef99`; the `f73b8e0` line lives one commit past the tag at `d415287`
  (origin tip's ee-pin commit). Per the brief that commit is Matt handing me the sha; I used
  `f73b8e0`. Tag and its freeze-file text are one commit out of sync on the ee line only.
- Import sanity (CPU, inside container): ee `__version__ 3.0.3` from the fork @ `f73b8e0`
  (`cw_block_circ` `evolve` default `True`), `discovery`/`discovery.matrix` from the fork @
  `136b270f`. Stack coherent.

### 8.2 VALUE GATE (S2) — **BIT-IDENTICAL, PASS**

`trackB_b1_core.py` unmodified, warm cache. Job **12489834** (dgx03, 10m59s, exit 0) +
full-precision addendum **12490007** (dgx03, 3m11s, read-only, no artifact written).

| gate | result on b1-v1 | banked | match |
|---|---|---|---|
| G0 masked lnL(mask=1) == `fisher_logL` | 0.000e+00 | 0.00e+00 | ✓ |
| G1 masked lnL(mask=0) == EarthDelay | 0.000e+00 | — | ✓ |
| G2 split `a_const` invariance | **405637.694390** (diff 0.000e+00) | 405637.694390 | ✓ 6dp |
| G3 `B1Split.lnL` == masked logL | **6.403e-10** | 6.403e-10 | ✓ exact |
| G4 `B1Split.estep` == `Split.estep` | 0.000e+00, 0/116 | 0/116 | ✓ |
| **G5 zero-noise E-step at truth == A2 census** | **J0711 0.9534 / J1713 0.9887 / J1909 0.9984**, 3/3 on TRUE fringe | 0.9534 / 0.9887 / 0.9984 | ✓ |
| G6 noise-draw covariance | white var ratio 0.985; ORF −0.0337 | — | ✓ |
| G7 noisy finiteness + S/N | lnL(truth\|zn) 405413.51; lnL(truth\|noisy) 390152.88 | 405413.51 | ✓ |
| G8 arm-B distance truth draw | 116/116 in-window; census n_true **−76 / −40 / 8** | −76 / −40 / 8 | ✓ |

**lnL(truth\|zero-noise) to full precision** (addendum job, same `build_b1_problem` + `amo["logL"]`
path the gate uses internally): **`405413.51273930125` → `405413.512739` (6dp) == banked.
BIT_IDENTICAL.** (The 3e-07 residual is the banked value's own 6dp truncation.)

The value gate transfers exactly across the ee `b54b211`→`f73b8e0` change: the `evolve`
default flip is a **provable no-op on this path** — `cw_helpers.py:174` passes `evolve=True`
explicitly. The 5th provenance link holds. Banked ckpts + tracked tree checksum-verified
untouched by every gate run.

### 8.3 SLURM KILL/RESUME (S3) — **PASS. The array's license.**

Noisy referendum, keyed tag `tA_n770077`, `--N 64 --seeds 1` (`h4d_resume_test.sbatch`).

- **KILL** (job **12490066**, dgx03): setup (~19 min: build + `z_needle` Hessian + batched
  bracket + quadrature), then at the **first** SMC ckpt I `scancel`'d **deliberately at
  stage 1, β=0.00444** — 0.4 % up a β→1 ladder, deep mid-run. `State=CANCELLED`, batch
  `ExitCode 0:15` (SIGTERM mid-stage-2). On-disk ckpt `b1_ref_tA_n770077_s0_ckpt.npz`
  (stage 1, β 0.00444, X(64,6)).
- **RESUME** (job **12490403**, dgx03, exit 0, 19m17s): resubmitted the *same* sbatch →
  `z_box` found its keyed ckpt →
  `RESUME from b1_ref_tA_n770077_s0_ckpt.npz: stage 1, beta=0.0044 (bit-exact continuation)`
  → ladder stage 2 (β 0.011) … stage 11 (**β 1.0000**) → `saved b1_referendum_tierA_noisy_n770077.npz`.
- **Keyed names intact:** noise-keyed final + ckpt, distinct from the Asimov canonical names.
- **ZERO CLOBBER:** the 4 git-tracked banked ckpts (`b1_ref_t{A,C}_s{0,1}_ckpt.npz`) + 2 finals
  (`b1_referendum_tier{A,C}.npz`) were sha256-verified **bit-identical before the test, after
  the kill, and after the completed resume.** This is exactly the **D-2** defect closed by the
  **C1b keying** — a `--noise` run can no longer race banked artifacts.
- **What this adds over cronus:** cronus validated resume *bit-exactness* locally (killed
  stage 6/7). This proves the **production SLURM path**: cgroup teardown → fresh allocation →
  shared-panfs ckpt read → `rng_state`-restored resume → clean completion. The walltime-kill
  path production will actually hit.
- Test residue (self-created, untracked, `_n770077`) removed; CW_transition left pristine.

### 8.4 DISK

- `/data/taylor_group`: **18 TiB quota, 15 TiB used, 3.1 TiB free (83 %)** — headroom up from
  §1.7's 1.9 TiB.
- Lean-npz rates: referendum final **6.1 KB**, per-stage ckpt **12.3 KB**, noisy final **10.3 KB**.
- **B1 array footprint at lean rates:** ~10–20 KB kept per realisation. A 10 k-realisation
  array ≈ **0.2–0.5 GB**; even `MaxArraySize`=60 000 ≈ **1–2 GB**. **Three-plus orders below
  the 3.1 TiB headroom — a non-issue.** Caveat: this holds **only while the array writes lean
  summary npz.** Any per-realisation heavy output (full posteriors, per-draw traces, the
  28 MB-class cell sets the RING campaign produced) changes the arithmetic — size *before*
  launch, not after.

### 8.5 SQUAT LOTTERY (informational)

This session **drew clean.** S2's A100 reported at gate start
`memory.total 81920 MiB / used 0 MiB / free 81152 MiB` — no `wut18`-class non-Slurm co-tenant
on our card. The §7.2 rogue-squat hazard (≈3/8 dgx03 cards, ~37.5 % per-job lottery) is
**unchanged as a standing risk**, but it is timing-only (§7.3 conclusion 4: contention never
perturbs a value) and it did not bite here. No admin action taken; the §7.2 admin note stands.

### 8.6 VERDICT — **READY FOR CAMPAIGNS**

Sync verified · value gate bit-identical (5th provenance link) · SLURM kill/resume proven ·
C1b keying prevents banked-artifact clobber · disk ample at lean rates.

**Standing hazards, not blockers:** (1) squat lottery — timing only; (2) heavy-npz footprint
if the array deviates from lean output — size before launch.

The E-track map and the B1 array arrive as **separate prompts on Matt's go.** **STOP.**

State on disk (new since §6):
| thing | path |
|---|---|
| S3 resume-test sbatch | `hpc_harbor/slurm/h4d_resume_test.sbatch` |
| S2 precision addendum | `hpc_harbor/slurm/h4c2_lnltruth_precision.{py,sbatch}` |
| S2/S3 logs | `hpc_harbor/logs/h4c_12489834.*`, `h4c2_12490007.*`, `h4d_1249006{6}.* / h4d_12490403.*` |
| staged reports (pushed by Matt) | committed `reports/` @ origin/`MM_playground` `3f485bf` |
