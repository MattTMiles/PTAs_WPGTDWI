#!/usr/bin/env bash
#SBATCH --job-name=gld2
#SBATCH --partition=batch_gpu
#SBATCH --account=taylor_group_acc
#SBATCH --qos=taylor_group_account_batch_gpu
#SBATCH --gres=gpu:nvidia_h200:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gld2_%j.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/gld2_%j.err
#
# S4.20b + S4.21: (1) FRONTIER-v2 first validation case (replay of the reference
# artifact gl1_null0_none_s1; pre-registered expectation: wrong-sky copy refused at
# the data-support term). (2) D2 RIGIDITY-GATE KILL CHECK on the cascade cells (gate
# pre-stated in S4.21 from SPARK's banked floors BEFORE this job existed).
set -uo pipefail
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_glacier/campaign"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
set +e   # harbor_env.sh sets -e; RC1/RC2 capture below needs it OFF
echo "[GLACIER-d2] host=$(hostname) job=${SLURM_JOB_ID} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
echo "[GLACIER-d2] gpu: $(nvidia-smi --query-gpu=name,uuid --format=csv,noheader)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_v2_frontier.py
RC1=$?
echo "[GLACIER-d2] v2-frontier exit=$RC1 wall=$((SECONDS-t0))s"
t1=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/glacier/gl_d2_rigidity.py
RC2=$?
echo "[GLACIER-d2] d2-kill exit=$RC2 wall=$((SECONDS-t1))s end=$(date -Is)"
[ $RC1 -eq 0 ] && [ $RC2 -eq 0 ]
