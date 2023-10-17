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


CIF_FILES="/wynton/home/sali/mhancock/xray/data/cifs/3ca7/3ca7.cif"
W_XRAY=1.0
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
N_STATE=2
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
SA="{step1000,T300,dofA,pdb1,w1,res0};{step1000,T1000,dofS,pdb0,w0,res-1};{step1000,T1000,dofS,pdb0,w0,res3}"
LOG_FILE="$OUT_DIR/log.csv"
INIT_WEIGHTS="rand"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files $CIF_FILES --w_xray $W_XRAY --dyn_w_xray --start_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights "$INIT_WEIGHTS" --ref_pdb_file $REF_PDB_FILE --sa "$SA" --u_aniso_file "/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/3ca7/3ca7_refine_001.pdb"


cd "$HOME/xray/sample_bench/scripts/sample"