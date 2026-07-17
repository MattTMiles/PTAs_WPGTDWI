#!/usr/bin/env bash
#SBATCH --job-name=sp3
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --exclude=dgx03
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3_%x_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3_%x_%j.err
set -uo pipefail
# cpus-per-task=8 is PART OF THE SEED (NoiseDrawer's eigh basis rotates with BLAS threads);
# the venue is T=40 -> no canonical L_gwb bank -> the RECOMPUTE branch -> gate g3a proves the
# reproduction of SURFACE's banked realisations rather than assuming it.
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/sp3_${SLURM_JOB_NAME}"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 8
echo "[SPARK3] device: $(nvidia-smi --query-gpu=uuid,memory.total,memory.used --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory,process_name --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[SPARK3] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} args=$SP3_ARGS start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark3.py $SP3_ARGS
RC=$?
echo "[SPARK3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
