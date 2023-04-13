#!/bin/bash


for JOB_ID in {0..106}
do
  if [[ $JOB_ID -lt 61 ]]
  then
    W_XRAY=$((( JOB_ID-0)*500 ))
    DYN_W_XRAY=0
  elif [[ $JOB_ID -lt 66 ]]
  then
    W_XRAY=$((( JOB_ID-61)*5000+30000 ))
    DYN_W_XRAY=0
  else
    JOB_ID_DIFF=$(( JOB_ID-66 ))
    W_XRAY=$( echo "$JOB_ID_DIFF*.1" | bc )
    DYN_W_XRAY=1
  fi
  echo "$JOB_ID", "$W_XRAY", "$DYN_W_XRAY"

  qsub run.sh "$JOB_ID" "$W_XRAY" "$DYN_W_XRAY"
done