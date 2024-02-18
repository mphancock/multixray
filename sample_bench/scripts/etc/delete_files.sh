#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7
# for JOB_DIR in $(ls -d  $EXP_DIR/*)
for EXP_NAME in 115_2_state 116_2_state_2_cif 117_2_state_no_weight 118_2_state_2_cif_no_weight 119_2_state_2_cif_w_1 120_2_state_2_cif_w_2 121_2_state_2_cif_w_3 122_52 123_52_2_cif 124_S_2_cif 125_2_cif_no_weights 126_test
do
    JOB_DIR="$SYS_DIR/$EXP_NAME"
    echo $JOB_DIR
    nohup rm -r $JOB_DIR &

    # for OUT_DIR in "$JOB_DIR"/output_*; do
    #     # Extract the Z value from the subfolder name using string manipulation
    #     folder_name=$(basename "$OUT_DIR")  # Get the base folder name without the path
    #     z_value=${folder_name#output_}        # Removes "output_" from the beginning of the folder name

    #     # Check if Z is greater than or equal to 50 using an if statement
    #     if [ "$z_value" -ge 50 ]; then
    #         echo "$OUT_DIR"
    #         rm -r "$OUT_DIR"
    #     fi
    # done

done