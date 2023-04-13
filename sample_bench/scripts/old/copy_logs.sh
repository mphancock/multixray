#!/bin/bash


JOB_NAME="single_md_3000_exp"
JOB_ID="346686"
for i in {0..99}
do
   LOG_FILE="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID/output_$i/log.csv"
   NEW_LOG_FILE="$HOME/xray/sample_bench/out/3ca7/$JOB_NAME/logs_$JOB_ID/log_$i.csv"

  cp "$LOG_FILE" "$NEW_LOG_FILE"
done