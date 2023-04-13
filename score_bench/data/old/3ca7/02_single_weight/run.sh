#!/bin/bash


JOB_NAME="single_weight"
DECOY_DIR="3ca7_N_1000_x1"
NATIVE_CIF_NAME="3ca7_clean.cif"
MIN_RES=-1
W_XRAYS=(0 1 2 3 4 5)

for i in {0..5}
do
  JOB_ID="$i"
  python ~/xray/score_bench/scripts/score_bench_3ca7.py --job_name "$JOB_NAME" --job_id "$JOB_ID" --decoy_dir "$DECOY_DIR" --native_cif_name "$NATIVE_CIF_NAME" --min_res "$MIN_RES" --w_xray "${W_XRAYS[$i]}"
done