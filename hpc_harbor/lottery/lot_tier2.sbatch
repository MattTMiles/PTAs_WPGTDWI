#!/usr/bin/env bash
#SBATCH --job-name=lot_t2
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --signal=B:USR1@300
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/lottery_%x_%A.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/lottery_%x_%A.err
#
# LOTTERY tier-2: GPU validation of the channel-budget switch-on proxy on the GENERAL queue.
# One node, one A100-80GB, cpus-per-task=8 (the reproducibility lever -> OpenBLAS 8 threads,
# so nulls and signals for every cell share ONE GWB basis on ONE host). RESUME-BY-SKIP on the
# keyed lot_*.npz => a preemption costs only the in-flight realisation. Reserved dgx_iacc share
# and SPARK-3 jobs are NOT touched.
set -uo pipefail
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
sleep 8
echo "[LOTTERY] device: $(nvidia-smi --query-gpu=uuid,memory.total,memory.used --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory,process_name --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[LOTTERY] host=$(hostname) job=${SLURM_JOB_ID:-none} cpus=${SLURM_CPUS_PER_TASK:-?} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/lottery/lottery_tier2.py run
RC=$?
echo "[LOTTERY] run exit=$RC wall=$((SECONDS-t0))s banked=$(ls /home/mattm/projects/HSYMT/LOTTERY_results/lot_*.npz 2>/dev/null | wc -l) end=$(date -Is)"
exit $RC
