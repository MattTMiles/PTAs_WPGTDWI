#!/usr/bin/env bash
#SBATCH --job-name=glfrnA
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfrnA_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfrnA_%j.err
#
# GLACIER g2a-ii REMEDY-A STUDY -- the band-extended population (flo -> -8.7, i.e. down
# to ~2 nHz, where the array's GWB sensitivity actually lives) at the MAIN-GATE venue
# (T=30, n=256). Measures: does the fitted background amplitude recover A_eff, and at
# what sigma -- the number Matt's remedy decision needs. Zero-noise Asimov: gate-class.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
du -sh /home/milesmt/.cache/jax_glacier 2>/dev/null || true
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
echo "[GLACIER-frnA] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-frnA] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_g2_forensic.py \
  --t 30 --n 256 --flo -8.7
RC=$?
echo "[GLACIER-frnA] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
