#!/usr/bin/env bash
#SBATCH --job-name=gllgates
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllgates_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllgates_%j.err
#
# GLACIER-LITE -- BUILD GATES ON THE CPU LANE (no GPU entitlement needed).
# The A100 group entitlement is saturated by another group member and the H200 pool is
# queue-deep to late October, so the GPU smoke is waiting. This run fast-fails the
# riskiest NEW code -- the venue rebuild -- without spending a GPU slot:
#   1. frozen-population gate (the strain ladder is the banked POP, not a fresh draw)
#   2. white-noise gate       (likelihood N_diag rms == target -> the idealization
#                              reached the COVARIANCE, not just the data)
#   3. banked T=15 L_gwb loads STRICT against the rebuilt venue
#   4. BackgroundFit constructs (>=1 in-band GP mode; identity checks pass)
#   5. drain profile finite at nothing-fed
# Deliberately does NOT call estep_per_target (the XLA-CPU vm.max_map_count hazard).
# Produces NO banks -- provenance is unaffected by the CPU lane.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/cpu_lite_gates"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
export JAX_PLATFORMS=cpu

# harbor_py without --nv (no GPU on this partition)
harbor_cpu() {
  apptainer exec --cleanenv \
    --bind "$HARBOR_BINDS" \
    --env XLA_PYTHON_CLIENT_PREALLOCATE=false \
    --env JAX_PLATFORMS=cpu \
    --env PYTHONUNBUFFERED=1 \
    --pwd /home/mattm/projects/HSYMT/CW_transition \
    "$HARBOR_SIF" "$HARBOR_ENV/bin/python" "$@"
}

echo "[gllgates] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
t0=$SECONDS
harbor_cpu /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lite.py gates --white w50
RC=$?
echo "[gllgates] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
