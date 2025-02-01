#!/bin/bash

# Define the source and destination directories
OUT_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/280_exp_all_2_phenix_ref/0/output_0"

run=1
if [ -f "$OUT_DIR/log.csv" ]; then
    run=0
fi

if [ "$run" -eq 0 ]; then
    echo "not running"
    exit 0
fi
echo "running"