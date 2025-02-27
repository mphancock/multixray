#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out

for JOB_DIR in $(ls -d  $SYS_DIR/285_2_state_3_cond/*)
do
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done