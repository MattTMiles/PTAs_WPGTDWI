#!/usr/bin/env bash
#SBATCH --job-name=comp_run
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=150G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/comp_%x_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/comp_%x_%A_%a.err
#
# COMPASS production fan. Usage:
#   sbatch --array=0-(NTASKS-1)%N --job-name=comp_<arm> comp_job.sbatch <ARM> <NTASKS>
# TWO-WAY GPU PACKING at the 8-cpu drawing convention (GENERALISE/TURBINE precedent):
# each task runs TWO worker processes, shards 2*TID and 2*TID+1 of 2*NTASKS.
# cpus-per-task=16 = 8 BLAS threads per drawing process, pinned via
# OPENBLAS/OMP/MKL_NUM_THREADS=8. Floors + signal for a unit ride the SAME process.
# 1-h walltime = the A-SKY backfill convention under the nodeupgrade drain.
# RESUME-BY-SKIP on keyed npz; a walltime kill costs the in-flight realisation only.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_generalise_shared"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
ARM="${1:?arm (I1|I2|S1)}"
NTASKS="${2:?ntasks}"
TID="${SLURM_ARRAY_TASK_ID:-0}"
NSH=$((2 * NTASKS))
sleep 8
echo "[COMP] device: $(nvidia-smi --query-gpu=name,uuid,memory.total,memory.used --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[COMP] host=$(hostname) job=${SLURM_JOB_ID:-none} arm=$ARM task=$TID/$NTASKS shards=$((2*TID)),$((2*TID+1))/$NSH cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"

comp_worker() {  # comp_worker <shard>
  apptainer exec --nv --cleanenv \
    --bind "$HARBOR_BINDS" \
    --env XLA_PYTHON_CLIENT_PREALLOCATE=false \
    --env XLA_FLAGS="$XLA_FLAGS" \
    --env PYTHONUNBUFFERED=1 \
    --env OPENBLAS_NUM_THREADS=8 \
    --env OMP_NUM_THREADS=8 \
    --env MKL_NUM_THREADS=8 \
    --pwd /home/mattm/projects/HSYMT/CW_transition \
    "$HARBOR_SIF" "$HARBOR_ENV/bin/python" \
    /home/mattm/projects/HSYMT/hpc_harbor/compass/compass.py run \
    --arm "$ARM" --shard "$1" --nshards "$NSH"
}

t0=$SECONDS
comp_worker $((2 * TID))     > >(sed "s/^/[w0] /") 2>&1 &
P0=$!
comp_worker $((2 * TID + 1)) > >(sed "s/^/[w1] /") 2>&1 &
P1=$!
wait $P0; R0=$?
wait $P1; R1=$?
RC=$(( R0 > R1 ? R0 : R1 ))
echo "[COMP] arm=$ARM task=$TID exit=$R0/$R1 wall=$((SECONDS-t0))s banked=$(ls /data/taylor_group/matt_miles/PTAs_WPGTDWI/COMPASS_results/comp_*.npz 2>/dev/null | wc -l) end=$(date -Is)"
exit $RC
