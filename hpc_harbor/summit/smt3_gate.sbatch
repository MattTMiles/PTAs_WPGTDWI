#!/usr/bin/env bash
#SBATCH --job-name=smt3gate
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt3gate_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smt3gate_%A_%a.err
#
# SUMMIT D3 PREP + GATES, one array task per extension rung:
#   task 0: lgwb ext+30  (array-keyed bank cut, cpus=8)  then gate battery +30
#   task 1: lgwb ext+100 then gate battery +100
# The gate battery: SG-D3b TWIN (clone at donor position reproduces the donor's
# scoreboard row -- every name-keyed inheritance verified at once), SG-D3c position
# audit, SG-D3e bank assert.   sbatch --array=0-1 smt3_gate.sbatch
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_summit/d3_shared_$SLURM_ARRAY_TASK_ID"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export APPTAINERENV_OMP_NUM_THREADS=8
export APPTAINERENV_OPENBLAS_NUM_THREADS=8
export APPTAINERENV_MKL_NUM_THREADS=8
NEXT=$([ "$SLURM_ARRAY_TASK_ID" -eq 0 ] && echo 30 || echo 100)
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/SUMMIT_results/d3_ext${NEXT}
mkdir -p /data/taylor_group/matt_miles/PTAs_WPGTDWI/SUMMIT_results/d3_ext${NEXT}
echo "[SMT3GATE] host=$(hostname) job=${SLURM_JOB_ID} ext=+${NEXT} start=$(date -Is)"
echo "[SMT3GATE] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
SD=/home/mattm/projects/HSYMT/hpc_harbor/summit/summit_d3.py
t0=$SECONDS
harbor_py $SD lgwb --next $NEXT || exit 2
harbor_py $SD gate --next $NEXT || exit 2
echo "=== SUMMIT D3 PREP+GATES ext+${NEXT}: ALL PASS === wall=$((SECONDS-t0))s end=$(date -Is)"
