#!/usr/bin/env bash
#SBATCH --job-name=gllloop
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllloop_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllloop_%A_%a.err
#
# GLACIER-LITE PART 1 -- THE LOOP DEMO.
# 8 tasks = {w50, w20} x realisations {0,1,2,3}; 8 laps each; 100-draw white-only floor
# re-cut EVERY lap (the banked floors were cut at 1414 ns white + real RN and do NOT
# apply at these noise levels -- the re-cut is mandatory and precedes all scoring).
#
# Task -> (white level, realisation):
#   0..3 = w50 r0..r3        4..7 = w20 r0..r3
#
# BUDGET-CUT ORDER if the 25 GPU-hr cap binds (pre-registered, from the brief):
#   1. realisations 4 -> 2 (cancel tasks 2,3 and 6,7)  -- lap depth is the point
#   2. drop the 50 ns level (cancel tasks 0..3) BEFORE touching 20 ns
# Every cut is reported. No silent truncation.
#
# Each lap checkpoints; a requeue resumes from the last banked lap.
# Kill step is SCORE-AND-LOG ONLY (capstone S4.23.1) -- it never executes.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
# per-task jax cache: defuses the shared-cache NOT_FOUND race across concurrent tasks
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_lite_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
export GLACIER_OUT=/home/mattm/projects/HSYMT/GLACIER_LITE_results
# harbor_py runs apptainer with --cleanenv, which STRIPS GLACIER_OUT. Apptainer does
# forward APPTAINERENV_*-prefixed vars through --cleanenv, so this is the one that
# actually reaches the container (job 12809208 banked into GLACIER_results without it).
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/GLACIER_LITE_results
set +e   # harbor_env.sh sets -e; the RC capture below needs it OFF
# dgx SMT: --cpus-per-task=8 hands back 16 LOGICAL cpus (8 physical x 2 SMT threads) on
# these nodes, which trips check_affinity's declared==8 convention and aborts the run
# (job 12803522). Bind one logical CPU per physical core -- gl_d2_a100.sbatch's fix,
# verbatim. CPU affinity is part of the provenance, so this is pinned, never relaxed.
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
echo "[pin] logical cpus -> $SET  (was $(python3 -c 'import os;print(len(os.sched_getaffinity(0)))'))"

WK=(w50 w50 w50 w50 w20 w20 w20 w20)
RL=(0 1 2 3 0 1 2 3)
i=${SLURM_ARRAY_TASK_ID}
echo "[gllloop] task=$i white=${WK[$i]} real=${RL[$i]}"
echo "[gllloop] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[gllloop] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lite.py loop --white ${WK[$i]} --real ${RL[$i]} --laps 8 --nnull 100"
RC=$?
echo "[gllloop] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
