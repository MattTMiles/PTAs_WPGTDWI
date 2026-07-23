#!/usr/bin/env bash
#SBATCH --job-name=gldrill
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=05:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gldrill_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gldrill_%j.err
#
# GLACIER RESUME DRILL (hold iii) -- *** DO NOT SUBMIT UNTIL FORGE-G2's runtime-SMASK IS
# IN THE PULLED TREE AND RE-GATED *** (the drill's bit-identity baseline must be cut on
# the FINAL tree; a pre-pull baseline is invalidated by the pull).
#
# LEG A1: drillcell --iters 2   (writes i0, i1 checkpoints)           } the "session loss"
# LEG A2: drillcell --iters 3   (i0, i1 must RESUME; i2 runs live)    } and the resume
# LEG B : drillcell --iters 3 --sub drill_straight   (uninterrupted control, own dir)
# CMP   : drillcmp -- the resumed leg's i2 bank must equal the control's, column-exact.
#
# Small census (n=32), carried floors (32 nulls at refit iters), T=30 venue, e07 arm.
# Separate processes per leg: clean memory, and A2's resume exercises COLD readback.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
du -sh /home/milesmt/.cache/jax_glacier 2>/dev/null || true
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
echo "[GLACIER-drill] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-drill] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
GL=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_loop.py

harbor_py "$GL" drillcell --iters 2 --n 32 --sub drill_resume \
  || { echo "[GLACIER-drill] LEG A1 FAILED"; exit 1; }
echo "[GLACIER-drill] leg A1 done t=$((SECONDS-t0))s"
harbor_py "$GL" drillcell --iters 3 --n 32 --sub drill_resume \
  || { echo "[GLACIER-drill] LEG A2 (resume) FAILED"; exit 1; }
echo "[GLACIER-drill] leg A2 done t=$((SECONDS-t0))s"
harbor_py "$GL" drillcell --iters 3 --n 32 --sub drill_straight \
  || { echo "[GLACIER-drill] LEG B (control) FAILED"; exit 1; }
echo "[GLACIER-drill] leg B done t=$((SECONDS-t0))s"
harbor_py "$GL" drillcmp --iters 3 --sub drill_resume --sub-b drill_straight
RC=$?
echo "[GLACIER-drill] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
