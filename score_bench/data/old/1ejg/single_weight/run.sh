#!/bin/bash


JOB_NAME="single_weight"
DECOY_DIR="1ejg_N_1000_x1"
NATIVE_CIF_NAME="1ejg_heavy.cif"
MIN_RES=1.5

W_XRAYS=(0 1 2 3 4 5)

for i in {0..5}
do
  JOB_ID="$i"
  python ~/xray/score_bench/scripts/score_bench_1ejg.py "$JOB_NAME" "$JOB_ID" "$DECOY_DIR" "$NATIVE_CIF_NAME" "$MIN_RES" "$((10**${W_XRAYS[$i]}))"
done