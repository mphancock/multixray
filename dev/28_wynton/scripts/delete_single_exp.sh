#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out
# SYS_DIR=/salilab/park1/matthew/xray

for JOB_DIR in $(ls -d  $SYS_DIR/280_exp_all_2_phenix_ref/*)
do
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done