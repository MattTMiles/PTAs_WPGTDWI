#!/usr/bin/env bash
#SBATCH --job-name=blfan
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-40gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/blfan_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/blfan_%A_%a.err
#
# BASELINE -- THE FAN.  usage:  sbatch --array=0-<N-1>%4 bl_fan.sbatch <ARM> <NSHARD>
#
# LANE: the BASELINE claim (GLACIER_capstone, 2026-07-29) -- interactive_gpu A100-40GB on
# dgx01 under dsi_dgx_iacc, capped at %4. PHOENIX holds the A100-80GB share on dgx03;
# this fan never requests that type. batch_gpu H200 is not entered.
#
# Every task is IDEMPOTENT (skip-on-exist), so a requeue or a re-submit resumes.
set -uo pipefail
ARM="${1:?arm C1|C2}"
NSHARD="${2:?nshard}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
# apptainer runs with --cleanenv, which STRIPS the thread vars. Forward them explicitly so
# the convention is DECLARED inside the container, not merely implied by taskset. (With a
# BANKED L_gwb the draw is thread-independent -- G-B0 -- but the convention is provenance.)
export APPTAINERENV_OMP_NUM_THREADS=8 APPTAINERENV_OPENBLAS_NUM_THREADS=8
export APPTAINERENV_MKL_NUM_THREADS=8 APPTAINERENV_NUMEXPR_NUM_THREADS=8
export APPTAINERENV_VECLIB_MAXIMUM_THREADS=8
# PER-TASK jax compile cache -- the harbor default, kept.
# I first set this to a SHARED cache on the reasoning that the documented race (concurrent
# array tasks colliding on xla_gpu_per_fusion_autotune_cache_dir/tmp/*.textproto) cannot
# fire with XLA_FLAGS=--xla_gpu_autotune_level=0. THAT REASONING IS WRONG: the smoke
# (12833821) created and kept writing that very directory under the pinned flag. So the
# race is live, and each task gets its own cache. Cost: the expensive first compile of the
# vmapped Adam graph is paid once per task instead of once per arm. Compile-cache isolation
# is NUMERICALLY inert (it caches XLA executables, not results).
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_baseline/${ARM}_t${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$HARBOR_JAXCACHE"
# ...but SEED it from the arm's warm smoke cache, so a private directory does not mean a
# private recompile. The vmapped Adam F_e graph is a very expensive first compile (the C1
# smoke spent >10 min on it); the smoke pays it once per arm and every fan task starts warm.
# The seed is HARDLINKED, not copied: /home/milesmt/.cache already holds 133 GB against a
# ~186 GiB quota, and 4 tasks x 4.9 GB of real copies per arm would put other campaigns'
# caches at risk. jax writes cache entries as new files (temp + rename) and never modifies
# one in place, so hardlinking existing entries is safe and costs zero space. Only top-level
# FILES are linked; xla_gpu_per_fusion_autotune_cache_dir/ is left to be created fresh per
# task, which is the directory the documented race is about.
WARM="/home/milesmt/.cache/jax_baseline/warm_${ARM}"
if [ -d "$WARM" ] && [ -z "$(ls -A "$HARBOR_JAXCACHE" 2>/dev/null)" ]; then
  N=$(find "$WARM" -maxdepth 1 -type f | wc -l)
  echo "[blfan] hardlink-seeding compile cache from $WARM ($N entries, 0 bytes copied)"
  find "$WARM" -maxdepth 1 -type f -exec cp -al {} "$HARBOR_JAXCACHE"/ \; 2>/dev/null || true
fi
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
set +e
SET=$(python3 -c "
import os
cpus = sorted(os.sched_getaffinity(0))
picked, seen = [], set()
for c in cpus:
    core = open(f'/sys/devices/system/cpu/cpu{c}/topology/core_id').read().strip()
    pkg = open(f'/sys/devices/system/cpu/cpu{c}/topology/physical_package_id').read().strip()
    if (pkg, core) not in seen:
        seen.add((pkg, core)); picked.append(c)
print(','.join(map(str, picked[:8])))")
echo "[pin] logical cpus -> $SET"
echo "[blfan] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} arm=${ARM} nshard=${NSHARD} start=$(date -Is)"
echo "[blfan] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/baseline/baseline.py run --arm $ARM --shard ${SLURM_ARRAY_TASK_ID} --nshard $NSHARD"
RC=$?
echo "[blfan] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
