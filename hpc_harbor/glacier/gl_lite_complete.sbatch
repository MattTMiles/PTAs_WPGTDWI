#!/usr/bin/env bash
#SBATCH --job-name=gllcmpl
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
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllcmpl_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllcmpl_%A_%a.err
#
# GLACIER-LITE PART 3 -- THE COMPLETENESS AXIS (amendment 2026-07-27, cap now 35 GPU-hr).
# 20 ns only, Part-1 conditioning (3 loudest at T2), prior widths scaled by
# {1, 1/3, 1/10, 1/30} with the HARD FLOOR sigma_a >= 3*dL[a] (never sub-fringe).
# Loop to saturation: up to 12 laps, stop after 2 consecutive laps adding 0 new certs.
#
# Task -> (prior-width factor, realisation):
#   0,1 = pw 1     r0,r1        4,5 = pw 1/10  r0,r1
#   2,3 = pw 1/3   r0,r1        6,7 = pw 1/30  r0,r1
#
# BUDGET CUT ORDER (pre-registered): drop the pw=1/3 factor (tasks 2,3) BEFORE cutting
# realisations or laps -- leaves {1, 1/10, 1/30}. Reported, never silent.
#
# Per cell the bank carries: the K distribution over the 116 pulsars (min/median/max),
# C per lap and at saturation, the saturation lap, and the per-pulsar autopsy of the
# never-certified set (K, pulsar-term SNR, and the pre-stated audibility-vs-ambiguity
# label). Kill step SCORE-AND-LOG ONLY.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_cmpl_${SLURM_ARRAY_TASK_ID}"
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

PW=(1.0 1.0 0.3333333333333333 0.3333333333333333 0.1 0.1 0.03333333333333333 0.03333333333333333)
RL=(0 1 0 1 0 1 0 1)
i=${SLURM_ARRAY_TASK_ID}
echo "[gllcmpl] task=$i pw=${PW[$i]} real=${RL[$i]}"
echo "[gllcmpl] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[gllcmpl] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lite.py complete --pw ${PW[$i]} --real ${RL[$i]} --laps 12 --nnull 100"
RC=$?
echo "[gllcmpl] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
