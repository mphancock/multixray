#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out

for JOB_DIR in $(ls -d  $SYS_DIR/283_2_cond_ref/*)
do
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done