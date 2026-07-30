#!/usr/bin/env bash
#SBATCH --job-name=gllanch
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllanch_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllanch_%j.err
#
# GLACIER-LITE PART 2 -- THE ANCHOR CELL (literature bridge, mandatory).
# SINGLE-SHOT, no loop: 20 ns white, 5 sources conditioned at FULL truth (T3, every
# parameter), per-pulsar distance posterior widths reported in both house conventions:
#   sigma_L_cond -- within-fringe conditional width (CW params KNOWN). This is the
#                   quantity the published single-shot forecasts report and the one
#                   comparable to their sub-pc claims.
#   sigma_L_marg -- fringe-comb sd, retaining the REAL fringe ambiguity, which this
#                   campaign is forbidden to idealize away.
# Compared against arXiv:2503.23017 (Yu & Pan) and arXiv:2512.10729 (Xiao et al.).
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_lite_anchor"
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
echo "[gllanch] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[gllanch] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lite.py anchor --white w20 --nsrc 5 --real 0"
RC=$?
echo "[gllanch] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
