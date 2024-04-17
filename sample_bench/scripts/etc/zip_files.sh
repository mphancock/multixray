#!/bin/bash


DATA_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf
DEST_DIR=/salilab/park1/matthew/xray

# zip -r $DEST_DIR/121_125.zip $DATA_DIR/121_native_decoys $DATA_DIR/122_native_decoys_1_state $DATA_DIR/123_natives_2_state $DATA_DIR/124_natives_2_cond $DATA_DIR/125_natives_1_state
# nohup zip -r /salilab/park1/matthew/xray/126_2_state_wxray/0.zip 0 &

for EXP_NAME in 126_2_state_wxray 127_2_cond_wxray 128_1_state_wxray
do
    EXP_DIR="$DATA_DIR/$EXP_NAME"
    mkdir -p /salilab/park1/matthew/xray/$EXP_NAME

    cd $EXP_DIR

    for JOB_DIR in $(ls -d $EXP_DIR/*)
    do
        JOB_NAME=$(basename $JOB_DIR)

        DEST_FILE="$DEST_DIR/$EXP_NAME/$JOB_NAME.zip"
        echo $JOB_NAME $DEST_FILE

        nohup zip -r $DEST_FILE $JOB_NAME &
    done
done