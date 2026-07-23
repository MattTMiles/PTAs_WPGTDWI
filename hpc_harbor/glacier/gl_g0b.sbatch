#!/usr/bin/env bash
#SBATCH --job-name=glg0b
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg0b_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glg0b_%j.err
#
# GLACIER g0b -- the FORGE-G2 runtime-SMASK battery + the Bg soft-source gates, ON THE
# H200 LANE, on the pulled tree (a2fd4df). Authorisation item 2 (2026-07-23).
#
# Two legs, fail-fast:
#   1. gl_smask_gate_h200.py -- SG battery in reference-absent mode (the cronus capture is
#      not bit-portable; see the wrapper docstring). Gated: SG2 bit-exact 0.0 on THIS host,
#      SG3 warm flip < 10 ms + cache-size constancy, SG5 grads, SG6 flip stability. The
#      one-time None->array structure retrace is paid in THIS job (t_first printed).
#   2. trackB_b1_softsource.py gate -- Bg2/Bg4 + Bg5a/b/c on the runtime-smask path
#      (same battery gl_g0 ran on the baked path, job 12635733, all 0.000e+00).
#
# CACHE: its own dir (g0) -- runs CONCURRENT with the g2 band re-gate, which owns the
# campaign dir; concurrent first-compiles on a shared dir are the EMBER race.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/g0"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLACIER-g0b] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-g0b] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
echo "[GLACIER-g0b] foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_smask_gate_h200.py \
  || { echo "[GLACIER-g0b] SMASK battery FAILED"; exit 1; }
echo "[GLACIER-g0b] smask battery PASS t=$((SECONDS-t0))s"
harbor_py /home/mattm/projects/HSYMT/CW_transition/trackB_b1_softsource.py gate
RC=$?
echo "[GLACIER-g0b] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
