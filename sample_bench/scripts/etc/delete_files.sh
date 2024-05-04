#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf
# SYS_DIR=/salilab/park1/matthew/xray

for JOB_DIR in $(ls -d  $SYS_DIR/182_bench/*)
do
    # echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done