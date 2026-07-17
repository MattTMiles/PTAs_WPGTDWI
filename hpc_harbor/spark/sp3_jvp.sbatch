#!/usr/bin/env bash
#SBATCH --job-name=sp3j
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3j_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3j_%A_%a.err
#
# ARM (c): THE BUDGETED CHUNKED-JVP. 16 marginal Hessians = 8 crossing-eligible units x rungs{0,8}.
# Grouped 2-4 per process (build ONCE per venue-pair, loop the share -- SPARK-2 cache lesson).
# 45-min cap PER HESSIAN enforced inside the process (spark3.py mode_jvp, cap_s=2700). A capped
# Hessian -> unit straddle-UNRESOLVED, banked as such; "unaffordable even chunked" over enough
# units is a pre-registered acceptable (AVALANCHE-shaped) outcome. Do NOT column-split a Hessian
# across tasks -- the cap bounds it and splitting buys coordination cost for nothing.
#   sbatch --array=0-4%4 sp3_jvp.sbatch     (5 processes, <=4 Hessians each)
set -uo pipefail
VEN=(A     A     A   B     B);   REALS=("0,1" "2,3" "4" "0,1" "2")
I=${SLURM_ARRAY_TASK_ID}
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/sp3_h200_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 8
echo "[SPARK3] host=$(hostname) job=${SLURM_JOB_ID} task=$I JVP venue=${VEN[$I]} reals=${REALS[$I]} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark3.py jvp \
    --venue "${VEN[$I]}" --reals "${REALS[$I]}"
RC=$?
echo "[SPARK3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
