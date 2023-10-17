#!/bin/bash

TMPDIR="/wynton/home/sali/mhancock/xray/sample_bench/data/out/test"
rm -r "$TMPDIR"
mkdir -p "$TMPDIR"

JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/data/out/test"
rm -r "$JOB_DIR"
mkdir -p "$JOB_DIR"

RUN_ID=0
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

CIF_FILES="/wynton/home/sali/mhancock/xray/dev/17_synthetic_native/data/cifs/2_state_ref/11.cif"
W_XRAY=1.0
DYN_W_XRAY=1
N_STATE=2
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/17_synthetic_native/data/pdbs/2_state_ref/11.pdb"
T=300
RES=0
LOG_FILE="$OUT_DIR/log.csv"

cd "$TMPDIR"
cp ~/xray/dev/18_sample_minor_state/scripts/run_18.py .
python run_18.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files "$CIF_FILES" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray --n_state "$N_STATE" --ref_pdb_file "$REF_PDB_FILE" --T "$T"

cd "$HOME/xray/dev/18_sample_minor_state/scripts"