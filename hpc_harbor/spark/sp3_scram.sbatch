#!/usr/bin/env bash
#SBATCH --job-name=sp3s
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3s_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3s_%A_%a.err
#
# ARM (a): THE SCRAMBLED-COUNTERPART ARM -- now LOAD-BEARING, not hygiene.
# The pessimistic-bound anatomy showed a poorly-constrained reservoir model ACCOMMODATES NOISE
# (floor 118->744). The optimistic edge's 11 crossings are exactly the kind of result that
# mechanism could manufacture. This arm scrambles the recovery counterpart so EVERY cert is FALSE
# by construction, at the rungs where the edge lives (0 baseline, 8 the edge), on the 5 crossing
# units. Manufacture (soft@opt raises false-cert count above unmodelled) -> STOP + anatomy.
# GOVERNS INTERPRETATION: no JVP reading ships until this reports.
set -uo pipefail
V=(A A A A A B B B B B); R=(0 0 1 1 4 4 0 0 1 1); K=(0 8 0 8 0 8 0 8 0 8)
# (5 crossing units A0 A1 A4 B0 B1, each at rungs 0 and 8)
I=${SLURM_ARRAY_TASK_ID}
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/sp3_h200_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 8
echo "[SPARK3] host=$(hostname) job=${SLURM_JOB_ID} task=$I SCRAM venue=${V[$I]} real=${R[$I]} rung=${K[$I]} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark3.py scram \
    --venue "${V[$I]}" --geo 3 --real "${R[$I]}" --rung "${K[$I]}" --n-null 100 --scr-mode L --self
RC=$?
echo "[SPARK3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
