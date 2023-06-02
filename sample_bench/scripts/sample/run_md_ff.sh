#!/bin/bash


JOB_NAME="20_decoys"
JOB_ID="0"
JOB_DIR="/home/matthew/xray/sample_bench/out/7mhk/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=0
OUT_DIR="$JOB_DIR/output_$RUN_ID"
rm -r "$OUT_DIR"
mkdir "$OUT_DIR"

COM=os
START_PDB_FILE="/home/matthew/xray/data/pdbs/7mhk/7mhk_clean.pdb"
N_STATE=1
REF_PDB_FILE="/home/matthew/xray/data/pdbs/7mhk/7mhk_clean.pdb"
T=1000
SA=0
STEPS=10000
LOG_FILE="$OUT_DIR/log.csv"

python run_md_multi.py --out_dir "$OUT_DIR" --com "$COM" --start_pdb_file "$START_PDB_FILE" --n_state "$N_STATE" --ref_pdb_file "$REF_PDB_FILE" --T "$T" --sa "$SA" --steps "$STEPS" --log_file "$LOG_FILE"
