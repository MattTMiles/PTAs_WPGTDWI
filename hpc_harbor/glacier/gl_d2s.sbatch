#!/usr/bin/env bash
#SBATCH --job-name=gld2s
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gld2s_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gld2s_%j.err
#
# A100 relocation of gl_d2.sbatch (hgx03 drained 2026-07-25 "update pending", no ETA;
# host freedom per Matt 2026-07-23 + the cross-host bit-identity gates).
# S4.23 STOP-MATERIAL KILL PASS: frozen S4.21 gate over the 7 STOPANAT banks.
# Report-only / non-pre-registered; decision rule recorded in S4.23.
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
set +e   # harbor_env.sh sets -e; RC1/RC2 capture below needs it OFF
sleep 5
echo "[GLACIER-d2s] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-d2s] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
# dgx SMT: bind one logical CPU per physical core (affinity convention == 8)
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
echo "[GLACIER-d2s] cpus $SET"
t0=$SECONDS
taskset -c "$SET" bash -c "harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_d2_stops.py"
RC=$?
echo "[GLACIER-d2s] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
