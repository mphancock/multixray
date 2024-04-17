#!/bin/bash

TMPDIR="/wynton/home/sali/mhancock/xray/sample_bench/out/test"
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

python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id 0 --w_xray 0.5 --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 16 --init_weights rand --sa "{step3000,T300,dofA,pdb1,w1,res0}" --steps 2

cd "$HOME/xray/sample_bench/scripts/sample"