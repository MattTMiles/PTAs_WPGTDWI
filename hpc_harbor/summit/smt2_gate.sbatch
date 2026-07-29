#!/usr/bin/env bash
#SBATCH --job-name=smt2gate
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt2gate_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt2gate_%j.err
#
# SUMMIT D2 GATES (one job, before the fan):
#   SG-D2b/c  CPU  tier prior columns + tiered truth-draw stream alignment
#   SG-D2a    GPU  wrapper inertness: stock lit path via the wrapper, i0 only,
#                  compared column-wise against the banked GLACIER ladder cell
# 1-h backfill class under the nodeupgrade drain (COMPASS pattern).
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_summit/d2_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/SUMMIT_results
echo "[SMT2GATE] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[SMT2GATE] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
SL=/home/mattm/projects/HSYMT/hpc_harbor/summit/summit_loop.py
t0=$SECONDS
harbor_py $SL gate      || exit 2
harbor_py $SL gatea     || exit 2
harbor_py $SL gateacmp  || exit 2
echo "=== SUMMIT D2 GATES: ALL PASS === wall=$((SECONDS-t0))s end=$(date -Is)"
