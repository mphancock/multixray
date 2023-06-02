#!/bin/bash

# TMPDIR="/wynton/home/sali/mhancock/xray/sample_bench/out/test"
TMPDIR="/home/matthew/xray/sample_bench/tmp"
rm -r "$TMPDIR"
mkdir "$TMPDIR"

# JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/test"
JOB_DIR="/home/matthew/xray/sample_bench/out/test"
rm -r "$JOB_DIR"
mkdir -p "$JOB_DIR"

RUN_ID=0
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

# CIF_FILE="/wynton/home/sali/mhancock/xray/data/reflections/7mhk/7mhk.cif"
CIF_FILE="/home/matthew/xray/data/reflections/7mhk/7mhk.cif"
RES=0
W_XRAY=.5
DYN_W_XRAY=1
COM=os
# START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhk/7mhk_clean_h20.pdb"
START_PDB_FILE="/home/matthew/xray/data/pdbs/7mhk/7mhk_clean.pdb"
# CIF_FILE="/wynton/home/sali/mhancock/xray/data/reflections/7mhk/7mhk.cif"
CIF_FILE="/home/matthew/xray/data/reflections/7mhk/7mhk.cif"
N_STATE=2
WEIGHTS=1
# REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhk/7mhk_clean_h20.pdb"
REF_PDB_FILE="/home/matthew/xray/data/pdbs/7mhk/7mhk_clean.pdb"
UC_DIM="114.300 54.290 44.970 90.00 102.12 90.00"
SG_SYMBOL="C 1 2 1"
T=300
SA=0
LOG_FILE="$OUT_DIR/log.csv"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_file "$CIF_FILE" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --com "$COM" --start_pdb_file "$START_PDB_FILE" --n_state "$N_STATE" --weights "$WEIGHTS" --ref_pdb_file "$REF_PDB_FILE" --uc_dim "$UC_DIM" --sg_symbol="$SG_SYMBOL" --T "$T" --sa "$SA" --log_file "$LOG_FILE"

cd "$HOME/xray/sample_bench/scripts/sample"