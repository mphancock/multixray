#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf
# SYS_DIR=/salilab/park1/matthew/xray

for JOB_DIR in $(ls -d  $SYS_DIR/167_N2/*)
# for EXP_NAME in 121_native_decoys 122_native_decoys_1_state 123_natives_2_state 124_natives_2_cond 125_natives_1_state 126_2_state_wxray 127_2_cond_wxray 128_1_state_wxray 129_2_state 130_2_cond 131_1_state 132_native_2_state_wxray 133_native_2_cond_wxray 134_native_1_state_wxray 135_4_state_1_cond_wxray 136_4_state_2_cond_wxray 137_4_state_1_cond 138_4_state_2_cond 139_native_4_state_1_cond_wxray 140_native_4_state_2_cond_wxray 141_native_4_state_1_cond 142_native_4_state_2_cond
do
    # JOB_DIR="$SYS_DIR/$EXP_NAME"
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done