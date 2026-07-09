#!/bin/bash
# B1 STEP 2 driver: strictly SERIAL (one process on the 24 GB card at a time).
# Concurrency was tried and failed: T_CHUNK=32 makes the per-pulsar activation
# (T, B, n_src, n_toa) ~2.6 GB for the 1222-TOA pulsar, and _rank_update then asks for
# 4.95 GiB -- tier B died RESOURCE_EXHAUSTED and tier C could not init CUDA.
set -u
cd /home/mattm/projects/HSYMT/CW_transition
PY=/home/mattm/miniforge3/envs/jug-gpu/bin/python
export XLA_FLAGS=--xla_gpu_autotune_level=0

# 1. tier A already banked ln Z_box from the first (broken-quadrature) run.
#    Recompute ONLY its needle with the bracketed quadrature.
echo "### tier A needle (bracketed quadrature)"
$PY -u trackB_b1_referendum.py --tier A --needle_only > logs/b1_ref_tierA_needle.log 2>&1
echo "### tier A needle rc=$?"

# 2. tiers B and C in full, one after the other.
for T in B C; do
  echo "### tier $T full"
  $PY -u trackB_b1_referendum.py --tier $T --N 192 --seeds 2 --mcmc 3 \
      > logs/b1_ref_tier${T}.log 2>&1
  echo "### tier $T rc=$?"
done
echo "### STEP 2 driver done"
