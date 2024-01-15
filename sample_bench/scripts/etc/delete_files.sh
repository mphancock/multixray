#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf
# for JOB_DIR in $(ls -d  $EXP_DIR/*)
for EXP_NAME in 55_7mhf_1 56_7mhf_2 57_7mhf_4 58_7mhf_8 59_7mhj_1 60_7mhj_2 61_7mhj_4 62_7mhj_8 63_overfit
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