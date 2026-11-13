#!/bin/bash
cd /home/mattm/projects/HSYMT/wheel_hasasia_examples
export MPLBACKEND=Agg OMP_NUM_THREADS=4; unset DISPLAY
PY=/home/mattm/miniforge3/envs/sculpt/bin/python
S20=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19
$PY job54_loc_sky.py main 10 10 $S20 3e7
$PY job54_loc_sky.py main 10 10 2 3e7 _half 0.5
$PY job54_loc_sky.py main 20 10 $S20 4e7
$PY job54_loc_sky.py partial 10 10 0,1,2,3,4,5,6,7 4e6
$PY job54_loc_sky.py main 10 20 $S20 1.5e7
$PY job54_loc_sky.py main 20 20 $S20 1.5e7
echo QUEUE_DONE
