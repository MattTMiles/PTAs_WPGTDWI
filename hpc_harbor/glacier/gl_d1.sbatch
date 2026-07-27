#!/usr/bin/env bash
#SBATCH --job-name=gld1
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gld1_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gld1_%j.err
#
# GLACIER Stage-1 DRIVER GATES (mode gate): G-d1 plumbing + G-d2 per-target E-step
# scoreboard (NEWLY WIRED -- CertScoreboard on spark3's audited estep_per_target /
# rung_masks / score_from_LNL) + G-d3 carried-statement identity. Zero-noise Asimov at
# T=15/n=32: gate-class, no campaign spend, no L_gwb.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
du -sh /home/milesmt/.cache/jax_glacier 2>/dev/null || true
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
echo "[GLACIER-dgate] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-dgate] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_d1_mechanism.py
RC=$?
echo "[GLACIER-dgate] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
