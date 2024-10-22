#!/bin/bash

TMPDIR="/wynton/group/sali/mhancock/xray/sample_bench/out/test_tmp"
rm -r "$TMPDIR"
mkdir -p "$TMPDIR"

JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/test"
rm -r "$JOB_DIR"
mkdir -p "$JOB_DIR"

RUN_ID=0
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"


cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .

python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --job_csv_file /wynton/home/sali/mhancock/xray/sample_bench/data/params/260.csv --job_id 0 --write

cd "$HOME/xray/sample_bench/scripts/sample"