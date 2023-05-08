#!/bin/bash


T=300
RES=0
JOB_NAME="46_w_xray"
for JOB_ID in {0..80}
do
  if [[ $JOB_ID -lt 50 ]]
  then
    W_XRAY=$(( (JOB_ID+1)*1000 ))
    DYN_W_XRAY=0
  else
    JOB_ID_DIFF=$(( JOB_ID-50 ))
    W_XRAY=$( echo "$JOB_ID_DIFF*.1" | bc )
    DYN_W_XRAY=1
  fi
  echo "$JOB_ID, $W_XRAY, $DYN_W_XRAY"

  qsub "$HOME/xray/sample_bench/scripts/sample/run_w_xray_3ca7.sh" "$JOB_NAME" "$JOB_ID" "$W_XRAY" "$DYN_W_XRAY" "$RES" "$T"
done