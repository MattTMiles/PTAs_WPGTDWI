#!/usr/bin/env bash
#SBATCH --job-name=svt2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svt2_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svt2_%j.err
#
# SIEVE T2 -- even/odd cross-validation. GENERAL CPU LANE. No GPU entitlement used, so
# every claimed GPU lane (SUMMIT / PHOENIX dgx03-80GB / BASELINE dgx01-40GB / frozen)
# is untouched. Three venue builds at ~640 s each (measured, probe 12834315) plus the
# noise drawer; 12 h is generous headroom.
# READS GLACIER_results/ (read-only). WRITES SIEVE_results/ only.
#
# usage: sbatch sv_t2.sbatch <rung> <arm> <sky> <iter>
set -uo pipefail
RUNG="${1:?rung}"; ARM="${2:?arm}"; SKY="${3:?sky}"; ITER="${4:?iter}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export APPTAINERENV_OMP_NUM_THREADS=8 APPTAINERENV_OPENBLAS_NUM_THREADS=8
export APPTAINERENV_MKL_NUM_THREADS=8 APPTAINERENV_NUMEXPR_NUM_THREADS=8
export APPTAINERENV_VECLIB_MAXIMUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_sieve/t2_${RUNG}_${ARM}_s${SKY}_i${ITER}"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export JAX_PLATFORMS=cpu
mkdir -p /data/taylor_group/matt_miles/PTAs_WPGTDWI/SIEVE_results

# one logical CPU per PHYSICAL core -- `--hint=nomultithread` alone returned a 16-CPU
# affinity mask on cn1486 and glacier_pop.check_affinity refused the run (correctly).
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
export -f harbor_cpu

echo "[svt2] host=$(hostname) job=${SLURM_JOB_ID} cell=${RUNG}/${ARM}/s${SKY}/i${ITER} start=$(date -Is)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_cpu /home/mattm/projects/HSYMT/hpc_harbor/sieve/sv_t2_xval.py --rung $RUNG --arm $ARM --sky $SKY --iter $ITER --max-fed ${MAXFED:-8}"
RC=$?
echo "[svt2] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
