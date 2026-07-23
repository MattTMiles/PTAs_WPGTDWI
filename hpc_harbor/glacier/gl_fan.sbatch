#!/usr/bin/env bash
#SBATCH --job-name=glfan
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfan_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfan_%A_%a.err
#
# GLACIER Stage-1 fan -- from TURBINE's T3 template (gl_fan_TEMPLATE.sbatch), with the
# 2-way GPU packing adoption: 16 cpus/task, TWO drawing processes of 8 cpus each on ONE
# H200 (measured aggregate 1.33x, convention-preserving; each process ~1.48x slower solo,
# walltime sized accordingly). cpus-per-task=8 PER DRAWING PROCESS is load-bearing.
#
# *** HELD *** -- submit ONLY after (i) FORGE-G2 runtime-SMASK is pulled + re-gated on this
# lane (all Bg5 PASS, bit-exact pattern flip 0.0, warm flip <10ms), (ii) driver gate G-d2
# (per-target E-step scoreboard wiring) lands + passes, (iii) the resume drill passes.
# The driver itself REFUSES mode cell/null until the holds clear (glacier_loop.py main).
#
# CELL MAP: task i drives cells (2i, 2i+1) of the flattened list
#   [ (arm, sky) for arm in (e07, none) for sky in 0..7 ] + null cells
#   [ (arm, sky, real) for arm in (e07, none) for sky in (0, 1) for real in 0..4 ]
# = 16 signal + 20 null cells = 36 cells -> array 0-17%6.
#
#   sbatch --array=0-17%6 gl_fan.sbatch
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign_shared"
mkdir -p "$HARBOR_JAXCACHE"
# cache-prune + quota prologue (TURBINE SS7: 52GB/41 entries measured on a 186GB quota)
echo "[GLFAN] cache: $(du -sh $HARBOR_JAXCACHE 2>/dev/null | cut -f1) quota: $(df -h /home/milesmt | tail -1)"
find "$HARBOR_JAXCACHE" -maxdepth 1 -atime +14 -delete 2>/dev/null || true
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLFAN] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLFAN] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
echo "[GLFAN] foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
GL=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_loop.py
CELLS=( $(python3 - <<'EOF'
cells = [f"cell --arm {a} --sky {s}" for a in ("e07","none") for s in range(8)]
cells += [f"null --arm {a} --sky {s} --real {r}" for a in ("e07","none")
          for s in (0,1) for r in range(5)]
print("\n".join(c.replace(" ", "@") for c in cells))
EOF
) )
t0=$SECONDS
# split the allocation into two 8-cpu sets -- each drawing process MUST see affinity == 8
# (the module's check_affinity refusal; BLAS thread count is part of the provenance)
ALLCPUS=($(python3 -c "import os;print(' '.join(map(str,sorted(os.sched_getaffinity(0)))))"))
SET0=$(IFS=,; echo "${ALLCPUS[*]:0:8}")
SET1=$(IFS=,; echo "${ALLCPUS[*]:8:8}")
i=$((2 * SLURM_ARRAY_TASK_ID))
PIDS=(); RC=0
for j in 0 1; do
  k=$((i + j))
  [ $k -lt ${#CELLS[@]} ] || continue
  ARGS=$(echo "${CELLS[$k]}" | tr '@' ' ')
  SET=$([ $j -eq 0 ] && echo "$SET0" || echo "$SET1")
  echo "[GLFAN] launching cell $k: $ARGS (cpus $SET)"
  taskset -c "$SET" bash -c "harbor_py $GL $ARGS" &
  PIDS+=($!)
done
for pid in "${PIDS[@]}"; do
  wait "$pid" || RC=1
done
echo "[GLFAN] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
