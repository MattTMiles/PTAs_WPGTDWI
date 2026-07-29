#!/usr/bin/env bash
#SBATCH --job-name=frzfan
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodelist=dgx03
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/frzfan_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/frzfan_%A_%a.err
#
# PHOENIX -- FROZEN-vs-LIVE M-step at the feasible rung (GLACIER S4.15.1 item 2).
# The FROZEN arm only: 4 cells = r13p9, circular arm, T=30, w=1, skies 0-3.
# The LIVE arm is INHERITED, never re-run (24/24 banked, dgx03 A100-80GB).
#
# ARRAY INDEX = SKY:  0 -> --sky 0 ... 3 -> --sky 3
# Floors: full on skies {0,1} per run_cell's stock policy (unchanged).
#
# LANE: PINNED to dgx03 (--nodelist). lane_tag() is host_GPU, NOT gpu-model-only, so
# the banked live arm's tag dgx03_NVIDIAA100-SXM4-80GB is reproducible on dgx03 ALONE.
# dgx04 frees 8 A100-80GB at 2026-07-29T13:02 but is a DIFFERENT HOST, and the measured
# host systematic (1.72 sigma) is not common-mode -- it would be a cross-lane comparison
# wearing a same-gres disguise. The driver refuses it; this pin makes the refusal moot.
# If dgx03 capacity ever binds, dgx04 can be LICENSED by re-running the SG-F1 inert gate
# there (0.3 GPU-hr) to measure the dgx03-vs-dgx04 offset directly -- not assumed.
#
# BUDGET: 4 cells x ~0.5 h measured feasible-rung cell wall = ~2 GPU-hr. dgx03 currently
# has 1 of 8 GPUs free (co-tenant holds 7), so expect ~2-2.5 h wallclock serialized
# rather than ~0.6 h at 4-wide; the %4 throttle takes any GPUs that free up.
# With the SG-F1 gate (~0.3) the frozen arm costs <= 3 GPU-hr.
# MUST DRAIN BY 2026-08-08 -- AugustDowntime2026 takes the whole cluster 08-10..08-13.
# At 0.6 h wallclock there is no exposure to that window.
#
# FRESH-DIR DISCIPLINE: the driver REFUSES a cell whose banks already exist; it does
# NOT skip-on-exist (SPARK-3C S5.4). Relaunching after a partial fan needs the partial
# cell's banks moved aside, deliberately.
#
#   sbatch --array=0-3%4 frz_fan.sbatch     [after frz_gate is SG-F1 PASS]
set -uo pipefail
# SAME warm cache the banked live arm compiled against (gl_fan2b.sbatch verbatim):
# a cold cache cost the first gate attempt a TIMEOUT (build_b1_amortised alone 653.7s).
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_shared"
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/FROZEN_results
sleep 5
echo "[FRZFAN] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[FRZFAN] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
FZ=/home/mattm/projects/HSYMT/hpc_harbor/frozen/frozen_mstep.py
k=$SLURM_ARRAY_TASK_ID
[ $k -ge 0 ] && [ $k -le 3 ] || { echo "[FRZFAN] index $k out of range"; exit 2; }
# dgx SMT: one logical CPU per physical core (affinity convention == 8).
# VERBATIM from gl_fan2b.sbatch -- the banked live arm ran under exactly this pin, and
# cpus-per-task is part of the provenance (thread-count hazard).
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
echo "[FRZFAN] frozen cell sky=$k (cpus $SET)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py $FZ cell --sky $k"
RC=$?
echo "[FRZFAN] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
