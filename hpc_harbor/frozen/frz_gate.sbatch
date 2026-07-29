#!/usr/bin/env bash
#SBATCH --job-name=frzgate
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodelist=dgx03
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:30:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/frzgate_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/frzgate_%j.err
#
# PHOENIX -- SG-F1: frozen-wrapper INERTNESS gate (GLACIER S4.15.1 item 2).
# Freeze DISABLED; one iteration of the stock path THROUGH the wrapper, column-compared
# to the banked live cell gl2_r13p9_cell_none_s0_T30_lit_i0__dgx03_NVIDIAA100-SXM4-80GB.
# Bar (SG-D2a's amended bar, same form): deterministic columns EXACT; Fisher-derived
# columns decision-invariant across the 0.5 feed bar + a 0.1 sanity cap.
#
# LANE IS NOT A CONVENIENCE: the live arm is banked on dgx03 A100-80GB. The host
# systematic is 1.72 sigma and not common-mode; the wander under measurement is
# 0.01-0.44 dex. The driver REFUSES any other lane.
#
# BUDGET: ~0.3 GPU-hr. Cumulative PHOENIX frozen-arm spend after gate+fan <= 3 GPU-hr,
# against GLACIER's standing cap. MUST DRAIN BY 2026-08-08 (AugustDowntime2026 takes the
# whole cluster 08-10..08-13). Walltime raised 01:00 -> 02:30 after attempt 12833771 TIMED
# OUT at 1:00 on a cold cache; still no exposure to the downtime window.
#
#   sbatch frz_gate.sbatch
set -uo pipefail
# SAME warm cache the banked live arm compiled against (gl_fan2b.sbatch verbatim):
# a cold cache cost the first gate attempt a TIMEOUT (build_b1_amortised alone 653.7s).
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_shared"
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/FROZEN_results
sleep 5
echo "[FRZGATE] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[FRZGATE] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
FZ=/home/mattm/projects/HSYMT/hpc_harbor/frozen/frozen_mstep.py
# dgx SMT: one logical CPU per physical core (affinity convention == 8).
# VERBATIM from gl_fan2b.sbatch -- the array ladder this arm is compared against ran
# under exactly this pin, and cpus-per-task is part of the provenance (thread-count
# hazard). Without it dgx03's 2 threads/core give affinity 16 and check_affinity STOPs.
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
echo "[FRZGATE] cpus $SET"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py $FZ gate --sky 0"
RC=$?
echo "[FRZGATE] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
