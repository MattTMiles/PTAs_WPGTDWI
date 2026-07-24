#!/usr/bin/env bash
#SBATCH --job-name=glfan2b
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfan2b_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfan2b_%A_%a.err
#
# GLACIER STAGE-2a ADDENDUM -- THE ARRAY-SIDE RUNGS (the dual ladder, Matt 2026-07-24).
# At the FEASIBLE sky rung only (r13p9, loudest NG15-consistent class): white-noise rms
# x{1, 1/2, 1/4} (catalogue RN + GWB held at measured levels -- they don't radiometer
# away), T in {30, 40}. Read out BOTH thresholds per cell: first-bite (N_resolved > 0)
# and certification (N_cert > 0) -- both banked already as per-iteration columns.
# DELIVERABLE: the ignition boundary in the (sky loudness x array capability) plane.
# Pre-registered question: does DSA-class rms (~1/4) bring the loudest NG15-consistent
# sky across the first-bite frontier at T=30-40?
#
# The (w=1, T=30) corner is NOT here -- it IS the sky ladder's r13p9 block (fan2 cells
# 0-7), reused. ARRAY INDEX = FLAT CELL INDEX:
#   0-39  signal: (w, T) in ((1,40), (0.5,30), (0.5,40), (0.25,30), (0.25,40))
#                 for arm in (e07, none) for sky in 0..3
#   40-51 null:   (w=0.25, T=40) for arm in (e07, none) for sky in (0,1) for real 0..2
#                 [scrambled arm at ONE array rung -- the likeliest-ignition corner]
# REQUIRES the T=40 L_gwb bank + its cross-host gate green (submit afterok on that).
#   sbatch --array=0-51%24 gl_fan2b.sbatch
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLFAN2B] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLFAN2B] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
GL=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_loop.py
CELLS=( $(python3 - <<'EOF'
combos = [(1.0, 40), (0.5, 30), (0.5, 40), (0.25, 30), (0.25, 40)]
cells = [f"cell2 --rung r13p9 --wscale {w} --t {t} --arm {a} --sky {s}"
         for (w, t) in combos for a in ("e07", "none") for s in range(4)]
cells += [f"null2 --rung r13p9 --wscale 0.25 --t 40 --arm {a} --sky {s} --real {k}"
          for a in ("e07", "none") for s in (0, 1) for k in range(3)]
print("\n".join(c.replace(" ", "@") for c in cells))
EOF
) )
k=$SLURM_ARRAY_TASK_ID
[ $k -lt ${#CELLS[@]} ] || { echo "[GLFAN2B] index $k out of range"; exit 2; }
ARGS=$(echo "${CELLS[$k]}" | tr '@' ' ')
# dgx SMT: one logical CPU per physical core (affinity convention == 8)
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
echo "[GLFAN2B] cell $k: $ARGS (cpus $SET)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py $GL $ARGS"
RC=$?
echo "[GLFAN2B] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
