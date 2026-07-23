#!/usr/bin/env bash
#SBATCH --job-name=glfrn
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfrn_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfrn_%j.err
#
# GLACIER g2a-ii FORENSIC -- H-BUG vs H-FAINT on the pegged fit gate (see
# glacier_g2_forensic.py). Zero-noise Asimov: gate-class, no campaign GPU, no L_gwb.
#
# LEG 1: the smoke venue (T=15, n=32) -- the venue that failed. Fast.
# LEG 2: the MAIN-GATE venue (T=30, n=256) -- runs ONLY if leg 1 refutes H-BUG
#        (measuring a broken pipe's information content twice is waste).
#
# LANE: batch_gpu / H200 (the declared floor host), cpus=8 (affinity refusal in-module).
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
du -sh /home/milesmt/.cache/jax_glacier 2>/dev/null || true
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
echo "[GLACIER-frn] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-frn] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
echo "[GLACIER-frn] foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
t0=$SECONDS
FR=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_g2_forensic.py

harbor_py "$FR" --t 15 --n 32 | tee /tmp/glfrn_leg1.$$ ; RC1=${PIPESTATUS[0]}
echo "[GLACIER-frn] leg1 exit=$RC1 t=$((SECONDS-t0))s"
if [ "$RC1" -ne 0 ]; then echo "[GLACIER-frn] leg 1 errored; stopping"; exit "$RC1"; fi
if grep -q "H-BUG CONFIRMED" /tmp/glfrn_leg1.$$; then
  echo "[GLACIER-frn] H-BUG confirmed at the smoke venue -- leg 2 skipped (fix the pipe first)"
  rm -f /tmp/glfrn_leg1.$$; exit 0
fi
rm -f /tmp/glfrn_leg1.$$

echo "[GLACIER-frn] leg 2 (T=30 n=256) t=$((SECONDS-t0))s"
harbor_py "$FR" --t 30 --n 256
RC2=$?
echo "[GLACIER-frn] exit=$RC2 wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC2
