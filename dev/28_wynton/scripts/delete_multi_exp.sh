#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out
# SYS_DIR=/salilab/park1/matthew/xray

EXP_NAMES=(262_sb_base 263_sb_sa 264_sb_sa_w 265_no_w 266_decoys 267_full 268_natives_1_state 269_minor 270_fix_w)

for EXP_NAME in "${EXP_NAMES[@]}"
do
    EXP_DIR=$SYS_DIR/$EXP_NAME
    echo $EXP_DIR
    nohup rm -r $EXP_DIR &
done