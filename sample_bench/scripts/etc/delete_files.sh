#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf
# SYS_DIR=/salilab/park1/matthew/xray

for JOB_DIR in $(ls -d  $SYS_DIR/178_wxray/*)
# for EXP_NAME in 126_2_state_wxray 127_2_cond_wxray 128_1_state_wxray 135_4_state_1_cond_wxray 136_4_state_2_cond_wxray 144_1_state_2_cond_wxray 149_N8_J1_wxray 150_N8_J2_wxray
do
    # JOB_DIR="$SYS_DIR/178_wxray/$JOB_DIR"
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done