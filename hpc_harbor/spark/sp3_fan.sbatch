#!/usr/bin/env bash
#SBATCH --job-name=sp3f
#SBATCH --partition=interactive_gpu
#SBATCH --account=dsi_dgx_iacc
#SBATCH --qos=dgx_iacc
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --exclude=gpu0058,gpu0060,gpu0061,gpu0207
#SBATCH --requeue
#SBATCH --output=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3f_%A_%a.out
#SBATCH --error=/data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/logs/sp3f_%A_%a.err
#
# THE ARROW-2 GRID, FANNED OUT BY RUNG.  48 jobs = 2 venues x sky g3 x 6 realisations x 4 rungs.
#   sbatch --array=0-47%12 sp3_fan.sbatch
#
# WHY RUNG GRANULARITY. The 4 rungs of a unit are INDEPENDENT (fresh pmask/Fisher/soft draw/null
# set at the same realisation); running them sequentially in one process buys only a shared
# build. One job per (unit,rung) is 4x the parallelism for the price of a rebuild. Below rung
# granularity the ~5-8 min build+compile stops amortising, so this is as fine as it pays to cut.
#
# THE LANE: gpu:1 + -N1 + exclude(non-A100, dgx03) = dgx01+dgx04 = 16 A100s. Every constraint
# is measured. (--nodelist=dgx01,dgx04 is WRONG here: it demands BOTH nodes per task, NODES=2.)
#  * NOT a bare gpu:1 without the exclude: interactive_gpu also holds RTX A4000 (gpu0058/60/61)
#    and a 3090 (gpu0207), which are 1:64-CRIPPLED in fp64 (HPC_SETUP §1.6). jax_enable_x64 is
#    MANDATORY here (pulsar-term phase ~1.6e4 rad; f32 destroys the ~1e-3-nat fringe-breaking
#    signal), so landing on one would be a correctness/perf disaster, not a slow job.
#  * NOT gpu:nvidia_a100-sxm4-80gb:1, which is what this campaign asked for until it was
#    MEASURED: a running rung job holds 988 MiB -- under 1 GB. The 80gb string bought nothing and
#    silently HALVED the lane by excluding dgx01's A100-40GB, whose fp64 is equally 1:2 full.
#  * dgx03 IS NOW ALLOWED, and the change is measured rather than inherited. SPARK barred it for
#    a rogue non-Slurm process squatting cards (SIREN: 51 GB resident; HPC_SETUP §7.2). Two facts
#    retire that bar FOR THIS WORKLOAD: (1) this job holds 988 MiB, so even a 51 GB squatter on
#    an 80 GB card leaves ~29 GB -- the OOM hazard the bar existed for is ~30x covered; and
#    (2) HPC_SETUP §7.3 conclusion 4 establishes that contention "never perturbs a value" -- it
#    costs timing only, and every number here is gated against the banks regardless. The sbatch
#    prints any foreign resident at startup, so the trail records what we actually landed beside.
#    NOTE 2026-07-17: dgx01 AND dgx04 both went DOWN ("Node unexpectedly rebooted"), so dgx03 is
#    the ONLY A100 capacity left. Barring it on inheritance would have idled the grid completely.
#
# THROTTLE %8: node contention is REAL and MEASURED (per-target E-step 3.62s at 1-way vs 20.0s
# with 4 of these on dgx04 = 5.5x) and it is a NODE effect, not a GPU-count effect -- those 4
# jobs held 4 distinct A100s. %8 across two nodes is the starting point; RE-MEASURE and adjust
# rather than trusting this number (co-tenancy has already corrupted three wallclock figures in
# this campaign, and there are foreign users on these nodes too).
set -uo pipefail
V=(); R=(); K=()
for v in A B; do for r in 0 1 2 3 4 5; do for k in 0 2 5 8; do V+=($v); R+=($r); K+=($k); done; done; done
I=${SLURM_ARRAY_TASK_ID}
# SHARED compile cache, pre-warmed by the gates. harbor_env: "Compile-cache isolation is
# NUMERICALLY inert: it caches XLA executables, not results." EMBER §7(a)'s race was on the
# autotune textprotos, and autotune is OFF here (XLA_FLAGS=--xla_gpu_autotune_level=0), so the
# race that motivated per-job caches is already defused. Sharing turns 48 cold compiles into 1.
export HARBOR_JAXCACHE="/home/milesmt/.cache/jax_spark/sp3_shared"
mkdir -p "$HARBOR_JAXCACHE"
source /data/taylor_group/matt_miles/PTAs_WPGTDWI/hpc_harbor/env/harbor_env.sh
export HARBOR_JAXCACHE_IN=/home/mattm/.cache/jax_stagec_cache
sleep 8
echo "[SPARK3] device: $(nvidia-smi --query-gpu=uuid --format=csv,noheader) | foreign: $(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader | tr '\n' ';' || echo none)"
echo "[SPARK3] host=$(hostname) job=${SLURM_JOB_ID} task=$I venue=${V[$I]} real=${R[$I]} rung=${K[$I]} cpus=${SLURM_CPUS_PER_TASK} start=$(date -Is)"
t0=$SECONDS
harbor_py /home/mattm/projects/HSYMT/hpc_harbor/spark/spark3.py unit \
    --venue "${V[$I]}" --geo 3 --real "${R[$I]}" --rung "${K[$I]}" --n-null 100
RC=$?
echo "[SPARK3] exit=$RC wall=$((SECONDS-t0))s end=$(date -Is)"
exit $RC
