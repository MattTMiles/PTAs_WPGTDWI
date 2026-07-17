#!/usr/bin/env bash
#SBATCH --job-name=sp3h
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3h_%x_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3h_%x_%j.err
#
# THE H200 LANE -- GATED, NOT ASSUMED.
# All three A100 nodes (dgx01/03/04) went DOWN ("Node unexpectedly rebooted") 2026-07-17, so the
# entire Ampere lane this campaign is gated on is unavailable. My group's batch_gpu QOS grants
# EXACTLY gres/gpu:nvidia_gh200=2 + gres/gpu:nvidia_h200=8 and NO A100 -- which is the measured
# reason gpu0084's A100s are denied (QOSGrpGRES), confirming HPC_SETUP §1.6's "QOS-denied" note
# by test rather than by reading. GH200 is aarch64 (jax's CUDA-12 x86 wheels will not run), so
# the H200 is the only entitled card with full 1:2 fp64 that can run this stack at all.
#
# WHY THIS MUST BE GATED BEFORE ANY GRID RUNS ON IT. H200 is Hopper (sm_90); every banked number
# in this campaign was produced on A100 Ampere (sm_80). XLA emits DIFFERENT KERNELS per
# architecture, and reductions may reorder -- so bit-identity is NOT free across the change.
# g3a is exactly the right instrument: it re-draws a banked SURFACE realisation and compares the
# raw columns. If g3a is 0.000e+00 on H200 as it was on A100, the two lanes are interchangeable
# and the grid may run here. If it is not, the H200 lane is REFUSED for this campaign -- because
# the certified sets are drawn from SURFACE's banks and a lane that cannot reproduce them would
# cohere a certified set onto a realisation that is not there.
# (autotune is already OFF -- XLA_FLAGS=--xla_gpu_autotune_level=0 -- which is what makes kernel
# selection deterministic enough for this test to mean anything.)
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/sp3_h200_${SLURM_JOB_NAME}"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 8
echo "[SPARK3] device: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[SPARK3] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} args=$SP3_ARGS start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark3.py $SP3_ARGS
RC=$?
echo "[SPARK3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
