#!/usr/bin/env bash
#SBATCH --job-name=smt2fan
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt2fan_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt2fan_%A_%a.err
#
# SUMMIT D2 -- the PRIOR-TIER single-dial ladder at the feasible rung (r13p9,
# circular arm, T=30, w=1; the boundary-mapping center). lit INHERITS from the
# banked GLACIER ladder; this fan runs the two NEW tiers + the scrambled arm at
# the most aggressive tier (subpc -- manufacture risk is highest where the prior
# is tightest).
#
# ARRAY INDEX = FLAT CELL INDEX (11 cells):
#   0-3   cell2 --tier vlbi  --sky 0..3
#   4-7   cell2 --tier subpc --sky 0..3
#   8-10  null2 --tier subpc: (s0,r0), (s0,r1), (s1,r0)
# Floors: full on skies {0,1} + all nulls (run_cell policy unchanged).
#   sbatch --array=0-10%4 smt2_fan.sbatch   [after smt2_gate is green]
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_summit/d2_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/SUMMIT_results
sleep 5
echo "[SMT2FAN] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[SMT2FAN] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
SL=/home/mattm/projects/HSYMT/hpc_harbor/summit/summit_loop.py
CELLS=(
  "cell2 --tier vlbi --sky 0"  "cell2 --tier vlbi --sky 1"
  "cell2 --tier vlbi --sky 2"  "cell2 --tier vlbi --sky 3"
  "cell2 --tier subpc --sky 0" "cell2 --tier subpc --sky 1"
  "cell2 --tier subpc --sky 2" "cell2 --tier subpc --sky 3"
  "null2 --tier subpc --sky 0 --real 0"
  "null2 --tier subpc --sky 0 --real 1"
  "null2 --tier subpc --sky 1 --real 0"
)
k=$SLURM_ARRAY_TASK_ID
[ $k -lt ${#CELLS[@]} ] || { echo "[SMT2FAN] index $k out of range"; exit 2; }
echo "[SMT2FAN] cell $k: ${CELLS[$k]}"
t0=$SECONDS
harbor_py $SL ${CELLS[$k]}
RC=$?
echo "[SMT2FAN] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
