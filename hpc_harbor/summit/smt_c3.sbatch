#!/usr/bin/env bash
#SBATCH --job-name=smtc3
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smtc3_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/smtc3_%A_%a.err
#
# SUMMIT PHASE 2, CELL BLOCK C3 -- the first COMPOSED cell: prior tier x array capability.
#   vlbi (1 pc on union-18)  x  w0p5 (rms/2)  at r13p9 / circular / T=30, skies 0-3.
#
# WHY C3 FIRST: it composes the two dials whose SINGLE-dial results are clean and whose
# controls are clean -- D2 vlbi (S1.3: moves TRUST, saturates at 1 pc; 3 subpc nulls clean)
# and the array ladder's w0p5 (S0.1/S4.24; w0p25 nulls clean). It is the only C-block that
# is NOT blocked by STOP #S3 (which parks the +100 rung used by C1/C2/N).
#
# THE QUESTION: every single dial leaves CERT at exactly zero at the feasible rung. D2 moves
# trust; w0p5 moves reach (sky 3 crosses FEED). Does trust x reach compose -- does holding
# the drain while the feed set widens put anything on the certification board, or is the
# 0/24 cert wall untouched by composition too?
#
# ARRAY INDEX = SKY (0-3). Floors: full on skies {0,1} per run_cell's stock policy.
#
# LANE: dgx03 pin is NOT required here (this block is not compared bit-for-bit against a
# banked arm), but the SMT affinity pin IS -- dgx nodes are 2 threads/core and
# check_affinity STOPs at 16.
#
# BUDGET: 4 cells x ~0.5-1.4 h measured = ~2-5 GPU-hr. SUMMIT cumulative after C3 ~ 20
# GPU-hr against the 150 STOP. MUST DRAIN BY 2026-08-08 (AugustDowntime2026 08-10..08-13);
# at <=4 h/cell there is no exposure.
#
#   sbatch --array=0-3%4 smt_c3.sbatch      [after smtc3gate is green]
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_shared"
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/SUMMIT_results
sleep 5
echo "[SMTC3] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[SMTC3] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
SL=/home/mattm/projects/HSYMT/hpc_harbor/summit/summit_loop.py
k=$SLURM_ARRAY_TASK_ID
[ $k -ge 0 ] && [ $k -le 3 ] || { echo "[SMTC3] index $k out of range"; exit 2; }
SET=$(python3 -c "
import os
cpus = sorted(os.sched_getaffinity(0))
picked, seen = [], set()
for c in cpus:
    core = open(f'/sys/devices/system/cpu/cpu{c}/topology/core_id').read().strip()
    pkg = open(f'/sys/devices/system/cpu/cpu{c}/topology/physical_package_id').read().strip()
    if (pkg, core) not in seen:
        seen.add((pkg, core)); picked.append(c)
print(','.join(map(str, picked[:8])))")
echo "[SMTC3] C3 composed cell: tier=vlbi wscale=0.5 T=30 sky=$k (cpus $SET)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py $SL cell2 --tier vlbi --wscale 0.5 --t 30 --sky $k"
RC=$?
echo "[SMTC3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
