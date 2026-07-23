#!/usr/bin/env bash
#SBATCH --job-name=glg0
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg0_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg0_%j.err
#
# GLACIER g0 -- FORGE-G's Bg5a/b/c on THIS cluster's GPUs.
#
# WHY: FORGE_G_build.md 3.3 ran the full B5 suite on cronus CPU (all 0.000e+00) and
# auto-deferred the GPU path (that box's GPU DNN library fails to init). The brief makes the
# GPU re-run a BLOCKING gate before any campaign spend.
#
# THE LANE: batch_gpu / H200. Deliberate, not incidental -- SPARK-3's EDGE-POSITIVE verdict was
# cut on the H200 (SELF) lane (SPARK3 4.5, "Which host's floors the verdict is cut on"), and
# GLACIER pins its floor host to the same lane so its scoreboard is comparable to the verdict
# that licensed it. The reserved dgx share (interactive_gpu) is SPARK-3C's and is NOT touched.
#
# WHAT "BIT-AGREEMENT" MEANS HERE, STATED BEFORE THE RUN (this is a scope line, not a hedge):
#   Bg5a/b/c are WITHIN-MACHINE CONTRASTS -- always-fed vs incumbent, truth vs one cycle,
#   w_s==0 vs junk params. Every one is a DIFFERENCE of two paths evaluated on the SAME device.
#   * PASS/STOP criterion: every such difference is 0.000e+00 on the H200, exactly as on cronus.
#     These are the gates' own terms and a non-zero here is a real defect -> STOP.
#   * REPORT-ONLY: the ABSOLUTE lnL values will not be bit-identical to cronus's CPU values.
#     CPU and GPU reduce in different orders; a ~1e-12 relative offset is arithmetic, not a
#     defect. Quoting that as a STOP would be mistaking a gate for a goal (SPARK3 4.3).
#   Zero-noise throughout, so NONE of this touches the L_gwb basis hazard -- which is exactly
#   why these gates are lane-portable and the campaign's noisy floors are not.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/g0"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLACIER-g0] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-g0] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
echo "[GLACIER-g0] foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/CW_transition/trackB_b1_softsource.py gate
RC=$?
echo "[GLACIER-g0] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
