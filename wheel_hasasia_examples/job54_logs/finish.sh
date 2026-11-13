#!/bin/sh
# Job 54 finisher: waits for the corrected rerun (pid 2191391), then regenerates tables and Parts 3-4 numbers. CPU only.
cd /home/mattm/projects/HSYMT/wheel_hasasia_examples
while kill -0 2191391 2>/dev/null; do sleep 30; done
export MPLBACKEND=Agg OMP_NUM_THREADS=1; unset DISPLAY
python job54_tables.py > job54_logs/tables.log 2>&1
JOB54_NW=8 python job54_parts34.py > job54_logs/parts34.log 2>&1
python talk_data/check.py > job54_logs/check.log 2>&1
date > job54_logs/FINISHED
