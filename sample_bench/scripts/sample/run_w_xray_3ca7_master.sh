#!/bin/bash


T=3000
RES=4
JOB_NAME="32_wxray_3000"
for JOB_ID in {0..81}
do
  if [[ $JOB_ID -lt 51 ]]
  then
    W_XRAY=$((( JOB_ID-0)*1000 ))
    DYN_W_XRAY=0
  else
    JOB_ID_DIFF=$(( JOB_ID-51 ))
    W_XRAY=$( echo "$JOB_ID_DIFF*.1" | bc )
    DYN_W_XRAY=1
  fi
  echo "$JOB_ID, $W_XRAY, $DYN_W_XRAY"

  qsub "$HOME/xray/sample_bench/scripts/sample/run_w_xray_3ca7.sh" "$JOB_NAME" "$JOB_ID" "$W_XRAY" "$DYN_W_XRAY" "$RES" "$T"
done