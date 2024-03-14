#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf
for EXP_NAME in 126_2_state_wxray 127_2_cond_wxray 128_1_state_wxray 129_2_state 130_2_cond 131_1_state
#  132_native_2_state_wxray 133_native_2_cond_wxray 134_native_1_state_wxray 135_4_state_1_cond_wxray 136_4_state_2_cond_wxray 137_4_state_1_cond 138_4_state_2_cond 139_native_4_state_1_cond_wxray 140_native_4_state_2_cond_wxray 141_native_4_state_1_cond 142_native_4_state_2_cond
do
    EXP_DIR="$SYS_DIR/$EXP_NAME"
    mkdir -p /salilab/park1/matthew/xray/$EXP_NAME

    for JOB_DIR in $(ls -d $EXP_DIR/*)
    do
        echo $JOB_DIR
        JOB_NAME=$(basename $JOB_DIR)

        DEST_DIR="/salilab/park1/matthew/xray/$EXP_NAME/$JOB_NAME"
        echo $DEST_DIR

        nohup cp -r $JOB_DIR $DEST_DIR &
    done
    # echo $JOB_DIR
    # nohup cp -r $JOB_DIR $DEST_DIR &
done