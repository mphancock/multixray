#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out
# SYS_DIR=/salilab/park1/matthew/xray

for JOB_DIR in $(ls -d  $SYS_DIR/271_native_2_wxray/*)
do
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done