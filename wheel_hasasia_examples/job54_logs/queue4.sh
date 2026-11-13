#!/bin/bash
cd /home/mattm/projects/HSYMT/wheel_hasasia_examples
export MPLBACKEND=Agg OMP_NUM_THREADS=4; unset DISPLAY
while kill -0 2104613 2>/dev/null; do sleep 5; done
/home/mattm/miniforge3/envs/sculpt/bin/python job54_loc_sky.py partial 10 10 8,9,10,11,12,13,14,15,16,17,18,19 4e6
echo QUEUE4_DONE
