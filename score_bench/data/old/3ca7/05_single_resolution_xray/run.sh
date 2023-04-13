#!/bin/bash


JOB_NAME="05_single_resolution_xray"
DECOY_DIR="$HOME/xray/decoys/data/3ca7/3ca7_N_1000_x1_xray"
CIF_FILE="$HOME/xray/data/reflections/3ca7/3ca7_clean.cif"
MIN_RES=(0 2 3 4)
W_XRAY=30000

for i in {0..3}
do
  JOB_ID="$i"
  python ~/xray/score_bench/scripts/score_bench_3ca7.py --job_name "$JOB_NAME" --job_id "$JOB_ID" --decoy_dir "$DECOY_DIR" --cif_file "$CIF_FILE" --min_res "${MIN_RES[$i]}" --w_xray "$W_XRAY"
done