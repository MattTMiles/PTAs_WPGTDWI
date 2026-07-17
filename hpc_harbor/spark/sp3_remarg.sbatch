#!/usr/bin/env bash
#SBATCH --job-name=sp3m
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3m_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3m_%A_%a.err
#
# THE STRADDLE SPLIT: re-score the soft state at the MARGINAL width (mode remarg), reading the
# JVP Hessian. 10 jobs = 5 crossing units (A0,A1,A4,B0,B1) x rungs {0,8}. Skips any whose JVP
# Hessian capped (straddle-unresolved for that unit). Idempotent (skip-on-exist).
#   sbatch --array=0-9%3 sp3_remarg.sbatch
set -uo pipefail
V=(A A A A A A B B B B); R=(0 0 1 1 4 4 0 0 1 1); K=(0 8 0 8 0 8 0 8 0 8)
I=${SLURM_ARRAY_TASK_ID}
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/sp3_h200_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 8
echo "[SPARK3] host=$(hostname) task=$I REMARG venue=${V[$I]} real=${R[$I]} rung=${K[$I]} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark3.py remarg \
    --venue "${V[$I]}" --geo 3 --real "${R[$I]}" --rung "${K[$I]}" --n-null 100 --self
RC=$?
echo "[SPARK3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
