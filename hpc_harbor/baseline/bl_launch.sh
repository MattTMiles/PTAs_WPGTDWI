#!/usr/bin/env bash
# BASELINE -- launch the fan. Run ONLY after both arms' smokes report GATES: 6/6.
#
# C1 runs first at %4; C2 is chained behind it with --dependency=afterany so the two never
# exceed the lane claim's 4-concurrent cap (GLACIER_capstone, 2026-07-29) and never contend
# with PHOENIX, which is on the A100-80GB share on dgx03.
#
# Every task is idempotent (skip-on-exist), so re-running this after a requeue, a preemption
# or a partial failure resumes rather than redoes.
#
#   C1: 130 realisations (30 sig + 100 null)  -> 33/33/32/32
#   C2: 128 realisations (120 sig + 8 xnull)  -> 32/32/32/32
#   measured wall 105 s/realisation + ~7 min venue+detector per task
set -euo pipefail
cd /data/taylor_group/matt_miles/PTAs_WPGTDWI

for ARM in C1 C2; do
  W="/home/milesmt/.cache/jax_baseline/warm_${ARM}"
  N=$(find "$W" -maxdepth 1 -type f 2>/dev/null | wc -l)
  echo "[launch] warm cache $ARM: $N entries"
  [ "$N" -ge 4 ] || { echo "*** REFUSING: $ARM has no warm compile cache ($W). Run the smoke first."; exit 1; }
done

JC1=$(sbatch --parsable --array=0-3%4 hpc_harbor/baseline/bl_fan.sbatch C1 4)
echo "[launch] C1 fan  -> $JC1  (array 0-3%4)"
JC2=$(sbatch --parsable --dependency=afterany:"$JC1" --array=0-3%4 \
      hpc_harbor/baseline/bl_fan.sbatch C2 4)
echo "[launch] C2 fan  -> $JC2  (array 0-3%4, after $JC1)"
echo
echo "score when both are COMPLETED:"
echo "  BL_HSYMT=/data/taylor_group/matt_miles/PTAs_WPGTDWI \\"
echo "    /home/milesmt/miniforge3/envs/jug-gpu/bin/python hpc_harbor/baseline/bl_score.py"
