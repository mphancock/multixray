#!/bin/bash


JOB_DIR="/wynton/home/sali/mhancock/xray/dev/27_sample_side_chains/data/out"
rm -r "$JOB_DIR"
mkdir -p "$JOB_DIR"

OUT_DIR="$JOB_DIR/output_0"
mkdir "$OUT_DIR"

START_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/27_sample_side_chains/data/0_state_0.pdb"
REF_PDB_FILE="$START_PDB_FILE"
SA="1000,1000,S,-1,1"
INIT_WEIGHTS="rand"
B_FACTOR=15

echo $SA
python $HOME/xray/sample_bench/scripts/sample/run_md_multi.py --out_dir $OUT_DIR --start_pdb_file "$START_PDB_FILE" --n_state 1 --init_weights "$INIT_WEIGHTS" --ref_pdb_file "$REF_PDB_FILE" --sa "$SA" --steps 1000