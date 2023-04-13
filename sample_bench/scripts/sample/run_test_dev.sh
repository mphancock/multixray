#!/bin/bash

TMPDIR="$HOME/xray/sample_bench/out/test"
rm -r "$TMPDIR"
mkdir "$TMPDIR"

JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/test"
rm -r "$JOB_DIR"
mkdir -p "$JOB_DIR"

RUN_ID=0
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

CIF_FILE="$HOME/xray/data/reflections/3ca7/3ca7.cif"
RES=0
W_XRAY=1
DYN_W_XRAY=1
COM=os
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_clean.pdb"
REF_PDB_FILE="$HOME/xray/data/pdbs/3ca7/3ca7_clean.pdb"
T=300
SA=0
LOG_FILE="$OUT_DIR/log.csv"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/md_3ca7_multi.py .
python md_3ca7_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_file "$CIF_FILE" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --com "$COM" --start_pdb_file "$START_PDB_FILE" --ref_pdb_file "$REF_PDB_FILE" --T "$T" --sa "$SA" --log_file "$LOG_FILE"

cd "$HOME/xray/sample_bench/scripts/sample"