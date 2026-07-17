#!/usr/bin/env bash
#SBATCH --job-name=sp3g
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --exclude=dgx03
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3g_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3g_%A_%a.err
#
# THE ARROW-2 GRID. 12 units = 2 venues x sky g3 x 6 realisations (every one rung-8-reachable).
#   sbatch --array=0-11%2 sp3_grid.sbatch
# THROTTLE %2 IS LOAD-BEARING, NOT POLITENESS. The per-target E-step is 3.62s solo but was
# MEASURED at 20.0s with 4 of these sharing dgx04 (5.5x) -- co-tenancy of my own jobs has now
# corrupted a wallclock number three times in this campaign (the "5-6x build" that was really
# 174s solo; gate1's absolutes; and a 4h walltime that would have timed out every unit at 9.3h).
# Units are idempotent PER RUNG, so a throttled queue costs wall-clock, never work.
set -uo pipefail
V=(A A A A A A B B B B B B); R=(0 1 2 3 4 5 0 1 2 3 4 5)
I=${SLURM_ARRAY_TASK_ID}
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/sp3_grid_${I}"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 8
echo "[SPARK3] device: $(nvidia-smi --query-gpu=uuid --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[SPARK3] host=$(hostname) job=${SLURM_JOB_ID} task=$I venue=${V[$I]} real=${R[$I]} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark3.py unit --venue "${V[$I]}" --geo 3 --real "${R[$I]}" --n-null 100
RC=$?
echo "[SPARK3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
