#!/usr/bin/env bash
#SBATCH --job-name=glfan2h
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfan2h_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glfan2h_%A_%a.err
#
# GLACIER STAGE-2a -- THE CONDITIONING LADDER (Option C, Matt 2026-07-24).
# The deliverable is THE CURVE: ignition (and where ignited, drain slope +
# COMPOUNDING/CONVERGENT) vs sky loudness above NG15.
#
# ARRAY INDEX = FLAT CELL INDEX:
#   0-23  signal: (rung, arm, sky) for rung in (r13p9, r13p5, r13p25)
#                 for arm in (e07, none) for sky in 0..3
#   24-35 null:   (r13p25, arm, sky, real) for arm in (e07, none)
#                 for sky in (0, 1) for real in 0..2   [scrambled arm at ONE rung]
# Floors: full on skies {0,1} per (rung, arm) + all nulls (run_cell policy unchanged).
#   sbatch --array=0-35%24 gl_fan2.sbatch
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 5
echo "[GLFAN2H] host=$(hostname) job=${SLURM_JOB_ID} task=${SLURM_ARRAY_TASK_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLFAN2H] gpu: $(nvidia-smi --query-gpu=name,uuid,memory.total --format=csv,noheader)"
GL=/home/mattm/projects/HSYMT/hpc_harbor/glacier/glacier_loop.py
CELLS=( $(python3 - <<'EOF'
cells = [f"cell2 --rung {r} --arm {a} --sky {s}"
         for r in ("r13p9", "r13p5", "r13p25")
         for a in ("e07", "none") for s in range(4)]
cells += [f"null2 --rung r13p25 --arm {a} --sky {s} --real {k}"
          for a in ("e07", "none") for s in (0, 1) for k in range(3)]
print("\n".join(c.replace(" ", "@") for c in cells))
EOF
) )
k=$SLURM_ARRAY_TASK_ID
[ $k -lt ${#CELLS[@]} ] || { echo "[GLFAN2H] index $k out of range"; exit 2; }
ARGS=$(echo "${CELLS[$k]}" | tr '@' ' ')
echo "[GLFAN2H] cell $k: $ARGS"
t0=$SECONDS
harbor_py $GL $ARGS
RC=$?
echo "[GLFAN2H] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
