#!/usr/bin/env bash
#SBATCH --job-name=glg2ba
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg2ba_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg2ba_%j.err
#
# GLACIER g2 RE-GATE AT THE REMEDY-A BAND (flo -8.7, ~2 nHz). Authorisation item 1
# (2026-07-23): "the conservation gate must hold over a band the array hears; the Ahat
# measurement licenses the remedy, the re-gate certifies it -- pass tolerance as
# pre-registered (0.15 dex), STOP if it misses."
#
# Same fail-fast ladder as gl_g2.sbatch, band-extended (stems are band-tagged, so the
# incumbent FAIL banks are not overwritten -- they are evidence):
#   1. gate  --flo -8.7            -- CPU gates g2a-i/g2b/g2c/g2d at the extended band
#   2. fitgate --t 15 --n 32 --flo -8.7   -- SMOKE: plumbing only (edge REPORT-ONLY;
#                                            the smoke venue is information-poor)
#   3. fitgate --t 30 --n 256 --flo -8.7  -- THE GATE: iteration-0 fitted background
#                                            amplitude vs drawn summed power, tol 0.15 dex
#      (forensic job 12706074 measured this venue: Ahat -14.629+-0.021 vs A_eff -14.583)
#
# CACHE: the campaign dir. The g0b smask battery runs concurrent on ITS OWN dir; the
# dgate re-run + drill are chained AFTER this job (afterok) so campaign-dir users are
# serialized (the EMBER first-compile race).
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
du -sh /home/milesmt/.cache/jax_glacier 2>/dev/null || true
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLACIER-g2ba] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-g2ba] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
echo "[GLACIER-g2ba] foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
t0=$SECONDS
GL=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_pop.py
harbor_py "$GL" gate --flo -8.7 || { echo "[GLACIER-g2ba] CPU gates FAILED"; exit 1; }
echo "[GLACIER-g2ba] smoke (T=15 n=32) t=$((SECONDS-t0))s"
harbor_py "$GL" fitgate --t 15 --n 32 --flo -8.7 || { echo "[GLACIER-g2ba] SMOKE fitgate FAILED"; exit 1; }
echo "[GLACIER-g2ba] main (T=30 n=256) t=$((SECONDS-t0))s"
harbor_py "$GL" fitgate --t 30 --n 256 --flo -8.7
RC=$?
echo "[GLACIER-g2ba] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
