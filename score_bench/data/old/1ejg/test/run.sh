#!/bin/bash


conda activate imp_218_cctbx


JOB_NAME="test"
DECOY_DIR="1ejg_N_1000_x1"
NATIVE_CIF_NAME="1ejg_heavy.cif"
MIN_RES=
W_XRAY=1

JOB_ID=0
python ~/xray/score_bench/scripts/score_bench_1ejg.py "$JOB_NAME" "$JOB_ID" "$DECOY_DIR" "$NATIVE_CIF_NAME" "$MIN_RES" "$W_XRAY"