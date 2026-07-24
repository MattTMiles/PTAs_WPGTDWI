#!/usr/bin/env bash
#SBATCH --job-name=glfanA
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfanA_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfanA_%A_%a.err
#
# GLACIER Stage-1 fan, A100 LANE (dgx use + host-freedom authorized by Matt 2026-07-23).
# VALID ONLY with the banked T=30 L_gwb (glacier_loop loads it strict; the cross-host
# gate gl_lgwb_xhost.py must be PASS on this lane before this file is submitted).
#
# ONE cell per task (no 2-way packing: 80 GB vs the H200's 141 GB -- memory headroom over
# walltime). ARRAY INDEX = FLAT CELL INDEX (same order as gl_fan.sbatch's CELLS list):
#   0-15  signal: (arm, sky) for arm in (e07, none) for sky in 0..7
#   16-35 null:   (arm, sky, real) for arm in (e07, none) for sky in (0,1) for real 0..4
# Submit a SUBSET complementary to the H200 fan -- the two arrays must not share cells
# (checkpoint collision), e.g.:  sbatch --array=10-35%24 gl_fan_a100.sbatch
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_shared"
mkdir -p "$HARBOR_JAXCACHE"
echo "[GLFANA] cache: $(du -sh $HARBOR_JAXCACHE 2>/dev/null | cut -f1) quota: $(df -h /home/milesmt | tail -1)"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLFANA] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLFANA] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
GL=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_loop.py
CELLS=( $(python3 - <<'EOF'
cells = [f"cell --arm {a} --sky {s}" for a in ("e07","none") for s in range(8)]
cells += [f"null --arm {a} --sky {s} --real {r}" for a in ("e07","none")
          for s in (0,1) for r in range(5)]
print("\n".join(c.replace(" ", "@") for c in cells))
EOF
) )
k=$SLURM_ARRAY_TASK_ID
[ $k -lt ${#CELLS[@]} ] || { echo "[GLFANA] index $k out of range"; exit 2; }
ARGS=$(echo "${CELLS[$k]}" | tr '@' ' ')
# dgx nodes run SMT: 8 cores arrive as 16 logical CPUs and the driver's affinity
# refusal (== 8, declared convention) fires. --hint=nomultithread is not honoured on
# this partition; bind to the first logical CPU of each core instead (one per core,
# affinity count 8 -- same as hgx03 where SMT is off).
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
echo "[GLFANA] cell $k: $ARGS (cpus $SET)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py $GL $ARGS"
RC=$?
echo "[GLFANA] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
