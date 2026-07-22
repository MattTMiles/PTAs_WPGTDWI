#!/bin/bash
# FORGE-G2 gate battery: runtime-SMASK gates, then the FORGE-G soft-source gates
# (Bg2/Bg4 logic + Bg5a/b/c GPU) on the new path. Sequential -- shared GPU.
set -u
cd /home/mattm/projects/HSYMT/CW_transition
PY=/home/mattm/miniforge3/envs/jug-gpu/bin/python
export XLA_FLAGS=--xla_gpu_autotune_level=0

echo "### FORGE-G2 runtime-smask gates"
$PY -u trackB_b1_smask_gate.py > logs/forge_g2_smask_gate.log 2>&1
echo "### smask gates rc=$?"

echo "### FORGE-G soft-source gates (Bg2/Bg4 + Bg5a/b/c) on the runtime-smask path"
$PY -u trackB_b1_softsource.py gate > logs/forge_g2_softsource_gate.log 2>&1
echo "### softsource gates rc=$?"
