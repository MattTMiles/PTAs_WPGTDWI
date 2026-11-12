#!/usr/bin/env bash
#SBATCH --job-name=ldga2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/ldga2_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/ldga2_%j.err
#
# LEDGER A2 -- RN MONITOR, REPLAYED on banked fan cells. CPU LANE ON PURPOSE.
#   * `rn_occupancy` never calls `estep_per_target`, so the XLA-CPU vm.max_map_count
#     hazard (the per-pulsar E-step thunk exhaustion) does not apply to this path.
#   * no GRES is requested, so this does not touch either GPU entitlement and cannot
#     contend with PHOENIX's frozen arm or BASELINE's share.
#   * the banked checkpoints are READ ONLY. Output goes to new ledger_a2_rn_* banks.
#
# Leg 1 is the CHEAP venue and is the real gate: it proves the monitor runs end to end on
# genuine machinery. Leg 2 is a production fan cell (Stage-1 none/s1, which feeds 2 members
# and carries the campaign's largest Stage-1 drain bite) and is allowed to fail on
# resources without taking leg 1 down -- a T=30/n=256 build on CPU is not something this
# programme has measured, so it is attempted, not promised.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/ledger_a2_cpu"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
export JAX_PLATFORMS=cpu

harbor_cpu() {
  apptainer exec --cleanenv --bind "$HARBOR_BINDS" \
    --env XLA_PYTHON_CLIENT_PREALLOCATE=false \
    --env JAX_PLATFORMS=cpu --env PYTHONUNBUFFERED=1 \
    --pwd /home/mattm/projects/HSYMT/CW_transition \
    "$HARBOR_SIF" "$HARBOR_ENV/bin/python" "$@"
}

echo "[LEDGER-A2] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
t0=$SECONDS

echo "--- leg 0: synthetic gates (must PASS before any venue number is read) ---"
harbor_cpu /home/mattm/projects/HSYMT/hpc_harbor/ledger/ledger_a2_rn_monitor.py gate
G=$?
echo "[LEDGER-A2] gates exit=$G"
if [ "$G" -ne 0 ]; then echo "[LEDGER-A2] gates FAILED -- not proceeding to any venue"; exit "$G"; fi

echo "--- leg 2: production fan cell, Stage-1 none/s1 (T=30, n=256) -- ATTEMPTED ---"
harbor_cpu /home/mattm/projects/HSYMT/hpc_harbor/ledger/ledger_a2_rn_monitor.py replay \
  --arm none --sky 1 --t 30 --n 256
R2=$?
echo "[LEDGER-A2] fan-cell replay exit=$R2 (non-zero is ALLOWED: CPU build cost is unmeasured)"

echo "[LEDGER-A2] wall=$((SECONDS-t0))s end=$(date -Is)"
exit 0
