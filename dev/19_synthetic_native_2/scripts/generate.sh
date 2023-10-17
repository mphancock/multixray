#!/bin/bash


N_STATE=4

for RUN_ID in {0..9}
do
    JOB_DIR="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/out/"$N_STATE"_state/$RUN_ID"
    rm -r "$JOB_DIR"
    mkdir -p "$JOB_DIR"

    OUT_DIR="$JOB_DIR/output_0"
    mkdir "$OUT_DIR"

    START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
    REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
    STEPS1=$((($RUN_ID+1)*100))
    SA="100,-1,$STEPS1,A;1000,-1,1000,S"
    INIT_WEIGHTS="rand"
    B_FACTOR=15

    echo $SA
    python $HOME/xray/sample_bench/scripts/sample/run_md_multi.py --out_dir $OUT_DIR --start_pdb_file "$START_PDB_FILE" --n_state "$N_STATE" --init_weights "$INIT_WEIGHTS" --ref_pdb_file "$REF_PDB_FILE" --T 300 --sa "$SA" --steps 100
done