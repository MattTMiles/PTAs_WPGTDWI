#!/usr/bin/env bash
#SBATCH --job-name=gllng15
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllng15_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllng15_%A_%a.err
#
# GLACIER-LITE SECONDARY ARM -- the NG15-renormalised contrast (Matt, 2026-07-27).
# 20 ns only, 2 realisations, LOOP ONLY (no anchor cell, no completeness axis).
#
# WHY. The primary arm runs the frozen 16-source POP, which carries 572.8x the NG15 band
# power (A_equivalent -13.22). If the loop ignites there, the reading is ambiguous: is it
# the loop, or is it the loudness? This arm rescales every member's strain by ONE common
# -1.38 dex so the summed band power equals the NG15 target exactly (A_equivalent -14.60),
# leaving geometry, frequencies, phases, sky and census rank order untouched.
#
# PRE-LOGGED PREDICTION (banked verbatim into every cell as PREDICTION_RENORM):
#   "ignition here would bound the harvest's structure-limited conclusion to sub-SKA
#    noise; non-ignition strengthens it."
#
# Banked under *_ng15 stems. NEVER merged with the primary arm's cells.
# BUDGET CUT (pre-registered): this arm's realisations go 2 -> 1 BEFORE anything in the
# primary arm is touched.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_ng15_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
export GLACIER_OUT=/home/mattm/projects/HSYMT/GLACIER_LITE_results
# harbor_py runs apptainer with --cleanenv, which STRIPS GLACIER_OUT. Apptainer does
# forward APPTAINERENV_*-prefixed vars through --cleanenv, so this is the one that
# actually reaches the container (job 12809208 banked into GLACIER_results without it).
export APPTAINERENV_GLACIER_OUT=/home/mattm/projects/HSYMT/GLACIER_LITE_results
set +e   # harbor_env.sh sets -e; the RC capture below needs it OFF
# dgx SMT: --cpus-per-task=8 hands back 16 LOGICAL cpus (8 physical x 2 SMT threads) on
# these nodes, which trips check_affinity's declared==8 convention and aborts the run
# (job 12803522). Bind one logical CPU per physical core -- gl_d2_a100.sbatch's fix,
# verbatim. CPU affinity is part of the provenance, so this is pinned, never relaxed.
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
echo "[pin] logical cpus -> $SET  (was $(python3 -c 'import os;print(len(os.sched_getaffinity(0)))'))"

i=${SLURM_ARRAY_TASK_ID}
echo "[gllng15] task=$i white=w20 real=$i RENORMALISED to NG15"
echo "[gllng15] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[gllng15] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lite.py loop --white w20 --real $i --laps 8 --nnull 100 --renorm"
RC=$?
echo "[gllng15] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
