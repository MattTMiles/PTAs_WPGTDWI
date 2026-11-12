#!/usr/bin/env bash
#SBATCH --job-name=smt3fan
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt3fan_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt3fan_%A_%a.err
#
# SUMMIT D3 -- the N_psr single-dial ladder at the feasible rung (r13p9, circular
# arm, T=30, w=1, lit priors). 116 INHERITS from the banked GLACIER ladder.
# ARRAY INDEX = FLAT CELL INDEX (11 cells):
#   0-3   cell2 --next 30  --sky 0..3     (SUMMIT_results/d3_ext30)
#   4-7   cell2 --next 100 --sky 0..3     (SUMMIT_results/d3_ext100)
#   8-10  null2 --next 100: (s0,r0), (s0,r1), (s1,r0)   [scrambled at the extreme rung]
#   sbatch --array=0-10%4 smt3_fan.sbatch   [after smt3_gate is green x2]
set -uo pipefail
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
sleep 5
CELLS=(
  "30 cell2 --sky 0"  "30 cell2 --sky 1"  "30 cell2 --sky 2"  "30 cell2 --sky 3"
  "100 cell2 --sky 0" "100 cell2 --sky 1" "100 cell2 --sky 2" "100 cell2 --sky 3"
  "100 null2 --sky 0 --real 0" "100 null2 --sky 0 --real 1" "100 null2 --sky 1 --real 0"
)
k=$SLURM_ARRAY_TASK_ID
[ $k -lt ${#CELLS[@]} ] || { echo "[SMT3FAN] index $k out of range"; exit 2; }
read -r NEXT MODE REST <<< "${CELLS[$k]}"
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_summit/d3_shared_$([ "$NEXT" = 30 ] && echo 0 || echo 1)"
mkdir -p "$HARBOR_JAXCACHE"
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/SUMMIT_results/d3_ext${NEXT}
echo "[SMT3FAN] host=$(hostname) job=${SLURM_JOB_ID} task=$k ext=+${NEXT} args=$MODE $REST start=$(date -Is)"
echo "[SMT3FAN] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
SD=/home/mattm/projects/HSYMT/hpc_harbor/summit/summit_d3.py
t0=$SECONDS
harbor_py $SD $MODE --next $NEXT $REST
RC=$?
echo "[SMT3FAN] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
