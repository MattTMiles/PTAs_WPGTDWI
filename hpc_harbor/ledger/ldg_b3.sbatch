#!/usr/bin/env bash
#SBATCH --job-name=ldgb3
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-40gb:1
#SBATCH --nodelist=dgx01
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --array=0-3%2
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/ldgb3_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/ldgb3_%A_%a.err
#
# LEDGER B3 -- FEED-ORDER (M-step coordinate-order) DEPENDENCE.
# One rung (r13p5), one arm (e07), 4 skies, 3 permutations + identity, REPLAYED from the
# banked i0 checkpoints. Reads GLACIER_results/, writes ONLY ledger_b3_* banks there.
#
# *** PARKED 2026-07-29: NOT SUBMITTED. ***
# `sbatch --test-only` returns QOSGrpGRES on BOTH live entitlements -- dgx_iacc (dgx03
# A100-80 is 8/8, dgx01 A100-40 carries BASELINE's %4) and taylor_group_account_batch_gpu.
# There is no GPU headroom for a third session without displacing PHOENIX's frozen arm or
# BASELINE's share, and neither is LEDGER's to take (LANE_CLAIM.md). Charter: POSTED AND
# PARKED -- a stopped arm never idles another lane. Submit this file unchanged the moment
# either entitlement frees; nothing in it needs re-deciding.
#
# LANE CLAIM (LEDGER), for the moment it can be honoured:
#   * one card at a time in practice (%2 cap, 4 tasks, <=2 h each) -- under 1 GPU-hr total
#   * job prefix ldg*; banks ledger_b3_* only; no write to any existing bank
#   * yields immediately to PHOENIX (frozen arm) and BASELINE; never enters `reserved`
#   * host-agnostic BY DESIGN: B3 is a WITHIN-RUN permutation contrast, so it carries no
#     cross-host comparison and does NOT need the dgx03 pin the three-way arm needs.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/ledger_b3_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache

SKY=${SLURM_ARRAY_TASK_ID}
echo "[LEDGER-B3] host=$(hostname) job=${SLURM_JOB_ID} sky=${SKY} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[LEDGER-B3] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/ledger/ledger_b3_feedorder.py run \
  --rung r13p5 --arm e07 --sky "${SKY}" --iter 0 --t 30
RC=$?
echo "[LEDGER-B3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
