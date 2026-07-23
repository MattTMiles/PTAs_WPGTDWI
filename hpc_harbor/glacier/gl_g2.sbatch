#!/usr/bin/env bash
#SBATCH --job-name=glg2
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg2_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg2_%j.err
#
# GLACIER g2 -- the drawn population + the drainable background, GATED (g2a/b/c/d).
#
# WHAT RUNS, IN ORDER (fail-fast):
#   1. mode gate            -- pure-numpy gates g2a-i/g2b/g2c/g2d (seconds)
#   2. fitgate --t 15 --n 32  -- SMOKE: the full fit path on the cheap base build, so a
#                                trivial plumbing bug costs minutes, not the T=30 build
#   3. fitgate --t 30 --n 256 -- THE GATE: iteration-0 fitted background amplitude vs the
#                                drawn population's summed power (tol 0.15 dex, pre-stated)
#
# LANE: batch_gpu / H200 (the declared floor host). The profile itself is a factorisation
# loop (f64 Cholesky mandate); zero-noise Asimov throughout, so no L_gwb noise basis is
# touched -- but the affinity refusal (cpus=8) is enforced by the module anyway.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
du -sh /home/milesmt/.cache/jax_glacier 2>/dev/null || true
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLACIER-g2] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-g2] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
echo "[GLACIER-g2] foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
t0=$SECONDS
GL=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_pop.py
harbor_py "$GL" gate || { echo "[GLACIER-g2] CPU gates FAILED"; exit 1; }
echo "[GLACIER-g2] smoke (T=15 n=32) t=$((SECONDS-t0))s"
harbor_py "$GL" fitgate --t 15 --n 32 || { echo "[GLACIER-g2] SMOKE fitgate FAILED"; exit 1; }
echo "[GLACIER-g2] main (T=30 n=256) t=$((SECONDS-t0))s"
harbor_py "$GL" fitgate --t 30 --n 256
RC=$?
echo "[GLACIER-g2] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
