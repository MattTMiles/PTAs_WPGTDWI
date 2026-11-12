#!/usr/bin/env bash
#SBATCH --job-name=svv4
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=10:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svv4_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svv4_%j.err
#
# SIEVE ADDENDUM V4 -- sigma-point E-step at working-regime belief widths
# (SIEVE-C T5 working-regime arm), GENERALISE arm A-SKY survivor cell, all 8 skies.
# GENERAL CPU LANE, no GPU entitlement: every claimed GPU lane (SUMMIT / PHOENIX /
# BASELINE / frozen) is untouched. ONE venue build (n_ecc=1, T=30, ncw=16 -- far
# smaller than GLACIER's ncw=287) serves all 8 skies.
# READS GENERALISE_results/ read-only. WRITES SIEVE_results/ only.
#
# usage: sbatch sv_v2.sbatch [cell] [nreal]
set -uo pipefail
CELL="${1:-1}"; NBEL="${2:-3}"; BUDGET="${3:-18000}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export APPTAINERENV_OMP_NUM_THREADS=8 APPTAINERENV_OPENBLAS_NUM_THREADS=8
export APPTAINERENV_MKL_NUM_THREADS=8 APPTAINERENV_NUMEXPR_NUM_THREADS=8
export APPTAINERENV_VECLIB_MAXIMUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_sieve/v4_c${CELL}"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export JAX_PLATFORMS=cpu
mkdir -p /data/taylor_group/matt_miles/PTAs_WPGTDWI/SIEVE_results

# one logical CPU per PHYSICAL core -- `--hint=nomultithread` alone returned a 16-CPU
# affinity mask on cn1486 and the campaign affinity check refused the run (correctly).
# 8 physical cores is ALSO the banked-noise-draw convention (NoiseDrawer eigh rotates
# the GWB draw with BLAS threads), so this pin is load-bearing twice over.
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

echo "[svv4] host=$(hostname) job=${SLURM_JOB_ID} cell=${CELL} nbelief=${NBEL} budget=${BUDGET} start=$(date -Is)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_cpu /home/mattm/projects/HSYMT/hpc_harbor/sieve/sv_v4_sigma.py --cell $CELL --nbelief $NBEL --budget-s $BUDGET"
RC=$?
echo "[svv4] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
