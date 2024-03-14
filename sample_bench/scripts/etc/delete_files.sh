#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7

# SYS_DIR=/salilab/park1/matthew/xray

# for JOB_DIR in $(ls -d  $EXP_DIR/*)
for EXP_NAME in 127_short 128_short_2_cif 129_wxray 130_low_res 131_low_res_2_cif 132_low_res_weights 133_low_res_2_cif_weights 134_res_15 135_res_25 136_noise 137_native_1_cif 138_native_2_cif 139_1_state 140_control 141_control_2 142_control_3 143_native_1_cif_2 144_native_2_cif_2 145_native_1_state_2 146_native_1_state_2 147_native_2_cif_2 148_native_1_state_2 152_native_1_cif 153_native_2_cif 154_native_1_state
do
    JOB_DIR="$SYS_DIR/$EXP_NAME"
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &
done