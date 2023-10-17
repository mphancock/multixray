#!/bin/bash


JOB_NAME="100_natives_4x"
SAMPLE_BENCH_DIR="/wynton/home/sali/mhancock/xray/sample_bench/data/3ca7/$JOB_NAME"

mkdir -p "$SAMPLE_BENCH_DIR"

for JOB_ID in {0..9}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID"
    REF_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/pdbs/4_state/$JOB_ID.pdb"

    FIELD="xray_0"
    FILE="$SAMPLE_BENCH_DIR/xray_volume_bench_$JOB_ID.csv"
    BONUS_FIELDS="rmsd_avg,pdb"
    python sample_volume_bench.py --job_dir "$JOB_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --file "$FILE"

    FIELD="rmsd_avg"
    FILE="$SAMPLE_BENCH_DIR/rmsd_volume_bench_$JOB_ID.csv"
    BONUS_FIELDS=""
    python sample_volume_bench.py --job_dir "$JOB_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --file "$FILE"
done