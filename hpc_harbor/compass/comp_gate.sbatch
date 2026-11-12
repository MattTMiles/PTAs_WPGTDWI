#!/usr/bin/env bash
#SBATCH --job-name=comp_gate
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:30:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/comp_gate_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/comp_gate_%j.err
#
# COMPASS gates gA (iso OFF, wrapper vs banked A-SKY) + g0 (iso ON, the mission gate:
# floor machinery in class at zero signal) + array audit + header estimate. One small
# job, STOP on failure, BEFORE any wide array. UNPACKED: cpus-per-task=8 is the pinned
# drawing convention (BLAS thread count is part of the noise seed); gates run at the
# exact convention. General H200 lane only; short walltime to ride backfill under the
# nodeupgrade drain.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_generalise_shared"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
MODE="${1:-gate}"
echo "[COMP] cache: $(du -sh $HARBOR_JAXCACHE 2>/dev/null | cut -f1)"
sleep 8
echo "[COMP] device: $(nvidia-smi --query-gpu=name,uuid,memory.total,memory.used --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[COMP] host=$(hostname) job=${SLURM_JOB_ID:-none} mode=$MODE cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/compass/compass.py "$MODE"
RC=$?
echo "[COMP] $MODE exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
