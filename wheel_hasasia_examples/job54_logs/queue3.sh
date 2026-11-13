#!/bin/bash
cd /home/mattm/projects/HSYMT/wheel_hasasia_examples
export MPLBACKEND=Agg OMP_NUM_THREADS=4; unset DISPLAY
PY=/home/mattm/miniforge3/envs/sculpt/bin/python
while kill -0 2090564 2>/dev/null; do sleep 5; done
$PY job54_loc_sky.py main 10 10 4 3e7 _maps 1.0 maps
S10=10,11,12,13,14,15,16,17,18,19
$PY job54_loc_sky.py main 10 20 $S10 1e7
$PY job54_loc_sky.py main 20 20 $S10 1e7
echo QUEUE3_DONE
