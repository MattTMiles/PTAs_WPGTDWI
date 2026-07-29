#!/usr/bin/env bash
#SBATCH --job-name=svt6b
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svt6b_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svt6b_%j.err
#
# SIEVE T6b -- common-random-numbers demo. GENERAL CPU LANE, no GPU entitlement used:
# the allocated 0.5 GPU-hr is NOT spent. Every claimed GPU lane (SUMMIT / PHOENIX
# dgx03-80GB / BASELINE dgx01-40GB / the frozen arm) is untouched.
# Banks to SIEVE_results/ ONLY (GLACIER_OUT is repointed) -- never GLACIER_results/,
# never FROZEN_results/.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export APPTAINERENV_OMP_NUM_THREADS=8 APPTAINERENV_OPENBLAS_NUM_THREADS=8
export APPTAINERENV_MKL_NUM_THREADS=8 APPTAINERENV_NUMEXPR_NUM_THREADS=8
export APPTAINERENV_VECLIB_MAXIMUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_sieve/t6b"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export JAX_PLATFORMS=cpu
mkdir -p /data/taylor_group/matt_miles/PTAs_WPGTDWI/SIEVE_results

# `--hint=nomultithread` is NOT sufficient on every node here: cn1486 handed back an
# affinity mask of 16 logical CPUs for cpus-per-task=8 and glacier_pop.check_affinity
# refused the run (correctly -- thread count is provenance). Pin to one logical CPU per
# PHYSICAL core, exactly as BASELINE's bl_fan.sbatch does.
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

export HARBOR_SIF HARBOR_ENV HARBOR_BINDS
harbor_cpu() {
  apptainer exec --cleanenv \
    --bind "$HARBOR_BINDS" \
    --env XLA_PYTHON_CLIENT_PREALLOCATE=false \
    --env JAX_PLATFORMS=cpu \
    --env PYTHONUNBUFFERED=1 \
    --env GLACIER_OUT=/home/mattm/projects/HSYMT/SIEVE_results \
    --pwd /home/mattm/projects/HSYMT/CW_transition \
    "$HARBOR_SIF" "$HARBOR_ENV/bin/python" "$@"
}

echo "[svt6b] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
t0=$SECONDS
export -f harbor_cpu
taskset -c "$SET" bash -c "harbor_cpu /home/mattm/projects/HSYMT/hpc_harbor/sieve/sv_t6b_crn.py --nsrc ${NSRC:-16} --iters ${ITERS:-3} --seeds ${SEEDS:-6}"
RC=$?
echo "[svt6b] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
