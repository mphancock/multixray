#!/bin/bash


JOB_NAME="152_native_1_cif"
SAMPLE_BENCH_DIR="/wynton/home/sali/mhancock/xray/sample_bench/data/3ca7/$JOB_NAME"

mkdir -p "$SAMPLE_BENCH_DIR"

for JOB_ID in {0..9}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID"
    REF_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_0/$JOB_ID.pdb"

    # FIELD="xray_0+xray_1"
    FIELD="xray_0"
    FILE="$SAMPLE_BENCH_DIR/xray_volume_bench_$JOB_ID.csv"
    BONUS_FIELDS="rmsd_avg_0,pdb"
    python sample_volume_bench.py --job_dir "$JOB_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --file "$FILE"

    FIELD="rmsd_avg_0"
    FILE="$SAMPLE_BENCH_DIR/rmsd_volume_bench_$JOB_ID.csv"
    BONUS_FIELDS=""
    python sample_volume_bench.py --job_dir "$JOB_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --file "$FILE"
done