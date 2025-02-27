#!/bin/bash

# Define the source and destination directories
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/287_3_state_2_cond/20"

EXP_DIR="${JOB_DIR%/*}"
JOB_ID="${JOB_DIR##*/}"

echo $EXP_DIR
echo $JOB_ID

run=1
if [ -f "$EXP_DIR"_phenix_ref/$JOB_ID/output_"1"/log.csv ]; then
    echo "not running"
    exit 0
fi