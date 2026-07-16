#!/usr/bin/env bash
#SBATCH --job-name=sp_s2
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:30:00
#SBATCH --exclude=dgx03
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp_s2_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp_s2_%j.err
set -uo pipefail
# Per-task JAX compile-cache isolation. MUST be set BEFORE sourcing harbor_env.sh
# (EMBER 7(a): the shared autotune cache races -> JaxRuntimeError NOT_FOUND, 11 kills).
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/s2"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 12
echo "[SPARK] device: $(nvidia-smi --query-gpu=uuid,memory.total,memory.used --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory,process_name --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[SPARK] host=$(hostname) job=${SLURM_JOB_ID:-none} cpus=${SLURM_CPUS_PER_TASK:-?} mode=s2 cache=$HARBOR_JAXCACHE start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark.py s2
RC=$?
echo "[SPARK] s2 exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
