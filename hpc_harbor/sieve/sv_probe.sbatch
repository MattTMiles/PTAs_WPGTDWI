#!/usr/bin/env bash
#SBATCH --job-name=svprobe
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svprobe_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/svprobe_%j.err
#
# SIEVE-A cost probe. CPU lane only (general partition, no GPU entitlement used).
# No banks written. Times the venue build so T2 can be scoped before it is queued.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_sieve/probe"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export JAX_PLATFORMS=cpu

harbor_cpu() {
  apptainer exec --cleanenv \
    --bind "$HARBOR_BINDS" \
    --env XLA_PYTHON_CLIENT_PREALLOCATE=false \
    --env JAX_PLATFORMS=cpu \
    --env PYTHONUNBUFFERED=1 \
    --pwd /home/mattm/projects/HSYMT/CW_transition \
    "$HARBOR_SIF" "$HARBOR_ENV/bin/python" "$@"
}

echo "[svprobe] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
t0=$SECONDS
harbor_cpu /home/mattm/projects/HSYMT/hpc_harbor/sieve/sv_probe.py
RC=$?
echo "[svprobe] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
