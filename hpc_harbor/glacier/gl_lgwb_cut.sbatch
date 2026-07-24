#!/usr/bin/env bash
#SBATCH --job-name=gllgwb
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllgwb_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gllgwb_%j.err
#
# Cut the canonical T=30 GLACIER L_gwb bank on the hgx03 lane at the cpus=8 convention.
# BLAS threads pinned BEFORE python starts -- the script refuses otherwise.
set -uo pipefail
# harbor_py runs apptainer --cleanenv: plain exports DO NOT reach the container.
# APPTAINERENV_* survives --cleanenv and lands as the bare name inside.
for v in OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS \
         VECLIB_MAXIMUM_THREADS NUMEXPR_NUM_THREADS; do
  export APPTAINERENV_${v}=8
done
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
echo "[GL-lgwb] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_lgwb_t30.py --t "${1:-30}"
RC=$?
echo "[GL-lgwb] exit=$RC end=$(date -Is)"
exit $RC
