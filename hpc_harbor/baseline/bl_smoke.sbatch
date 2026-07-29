#!/usr/bin/env bash
#SBATCH --job-name=blsmoke
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-40gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/blsmoke_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/blsmoke_%j.err
#
# BASELINE -- SMOKE (standing process rule: no array submit before a smoke).
# ONE venue build + the G-B0..G-B5 gate chain + one signal and one null realisation
# timed end to end, so the fan is sized on a MEASURED per-realisation wall, not a guess.
#
# LANE: the BASELINE lane claim (GLACIER_capstone, 2026-07-29) -- interactive_gpu A100
# under dsi_dgx_iacc, capped at 4 concurrent. batch_gpu H200 is not entered (its
# --test-only start date is 2026-11-04, still behind the nodeupgrade drain).
# cpus-per-task=8 is the pinned convention (noise-draw thread hazard); with a BANKED
# L_gwb the draw is thread-independent, but the convention is kept and printed anyway.
set -uo pipefail
ARM="${1:-C1}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
# apptainer runs --cleanenv and STRIPS these; forward them so the pinned 8-thread
# convention is DECLARED inside the container (job 12833821 logged OPENBLAS=None).
export APPTAINERENV_OMP_NUM_THREADS=8 APPTAINERENV_OPENBLAS_NUM_THREADS=8
export APPTAINERENV_MKL_NUM_THREADS=8 APPTAINERENV_NUMEXPR_NUM_THREADS=8
export APPTAINERENV_VECLIB_MAXIMUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_baseline/warm_${ARM}"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
set +e
# dgx SMT: --cpus-per-task=8 hands back 16 LOGICAL cpus; bind one logical cpu per
# physical core (gl_d2_a100.sbatch's fix, verbatim -- affinity is part of the provenance).
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
echo "[blsmoke] host=$(hostname) job=${SLURM_JOB_ID} arm=${ARM} start=$(date -Is)"
echo "[blsmoke] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/baseline/baseline.py smoke --arm $ARM ${2:-}"
RC=$?
echo "[blsmoke] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
