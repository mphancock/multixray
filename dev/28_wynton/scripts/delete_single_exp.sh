#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out

for JOB_DIR in $(ls -d  $SYS_DIR/277_native_5_ref/*)
do
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done