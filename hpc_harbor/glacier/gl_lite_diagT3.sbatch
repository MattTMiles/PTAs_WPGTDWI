#!/usr/bin/env bash
#SBATCH --job-name=glldiag
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glldiag_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glldiag_%j.err
#
# GLACIER-LITE -- 1-CELL SMOKE (standing process rule: no array submit before a smoke).
# ONE venue build at 50 ns + ONE lap with a SHORT null set (32 draws).
# What it must prove before anything wide runs:
#   1. the frozen-population gate  (log10_h ladder == POP (3,-13.25,-14.25))
#   2. the WHITE-NOISE gate        (likelihood N_diag rms == 50 ns, i.e. the
#                                   idealization actually reached the covariance, not
#                                   just the data -- the deviation-(b) claim)
#   3. the banked T=15 L_gwb loads STRICT against the rebuilt venue (proves the GWB
#                                   block is untouched by the white-noise change)
#   4. BackgroundFit constructs     (>=1 in-band GP mode at band (-8.0,-7.5): the drain
#                                   has a basis; its two identity checks pass)
#   5. one full lap end-to-end     (Fisher -> frontier-v2 -> promote -> drain -> M-step
#                                   -> floor re-cut -> criterion v2.2 -> kill score)
# Stems carry '_smoke' so these banks can never satisfy the production resume check.
# A100 lane: the free lane on this cluster (batch_gpu H200 is queue-deep, HPC_SETUP 1.3).
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8 NUMEXPR_NUM_THREADS=8
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_shared"
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
echo "[glldiag] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[glldiag] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lite.py smoke --white w50 --cond T3"
RC=$?
echo "[glldiag] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
