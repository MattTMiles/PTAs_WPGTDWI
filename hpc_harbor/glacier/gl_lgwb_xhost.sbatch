#!/usr/bin/env bash
#SBATCH --job-name=glxhost
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glxhost_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/glxhost_%j.err
#
# Cross-host L_gwb identity gate ON THE dgx A100 LANE (dgx use authorized by Matt
# 2026-07-23, "don't worry about the not-ours thing"). PASS -> draws are host-free with
# the banked L; the fan may spread to the A100 pool.
set -uo pipefail
# APPTAINERENV_* survives harbor_py's --cleanenv (plain exports do not)
for v in OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS \
         VECLIB_MAXIMUM_THREADS NUMEXPR_NUM_THREADS; do
  export APPTAINERENV_${v}=8
done
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/a100"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
echo "[GL-xhost] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GL-xhost] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lgwb_xhost.py --t "${1:-30}"
RC=$?
echo "[GL-xhost] exit=$RC end=$(date -Is)"
exit $RC
