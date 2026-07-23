#!/usr/bin/env bash
#SBATCH --job-name=lot_band
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/lottery_%x_%A.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/lottery_%x_%A.err
#
# LOTTERY re-cut band: ONE basis state of the m1e03/m2e03 ensemble.
# --cpus-per-task is deliberately NOT set here -- it IS the basis knob and must be passed on the
# sbatch command line together with LOT_BASIS, e.g.
#     sbatch --cpus-per-task=4  --export=ALL,LOT_BASIS=b04 lot_band.sbatch
#     sbatch --cpus-per-task=16 --export=ALL,LOT_BASIS=b16 lot_band.sbatch
#     sbatch --cpus-per-task=2  --export=ALL,LOT_BASIS=b02 lot_band.sbatch
# The driver REFUSES to run if cpus-per-task disagrees with the ensemble spec for LOT_BASIS, so a
# mislabelled basis cannot enter the bank. General queue only; reserved dgx_iacc share and every
# SPARK-3 job untouched.
set -uo pipefail
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
sleep 8
echo "[BAND] device: $(nvidia-smi --query-gpu=uuid,memory.total,memory.used --format=csv,noheader)"
echo "[BAND] host=$(hostname) job=${SLURM_JOB_ID:-none} cpus=${SLURM_CPUS_PER_TASK:-?} basis=${LOT_BASIS:-unset} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/lottery/lottery_recut_band.py run "${LOT_BASIS:?LOT_BASIS must be exported: sbatch --export=ALL,LOT_BASIS=b04 ...}"
RC=$?
echo "[BAND] exit=$RC wall=$((SECONDS-t0))s basis=${LOT_BASIS:-unset} end=$(date -Is)"
exit $RC
