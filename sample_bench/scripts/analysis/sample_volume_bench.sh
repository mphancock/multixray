#!/bin/bash


JOB_NAME="60_sb_2"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine_4.pdb"
JOB_ID="134569"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID"
SAMPLE_BENCH_DIR="/wynton/home/sali/mhancock/xray/sample_bench/data/3ca7/$JOB_NAME"


FIELD="xray_0"
BONUS_FIELDS="pdb,rmsd,rmsd_avg,r_free_0"
python sample_volume_bench.py --job_dir "$JOB_DIR" --sample_bench_dir "$SAMPLE_BENCH_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS"

FIELD="rmsd_avg"
BONUS_FIELDS="pdb"
python sample_volume_bench.py --job_dir "$JOB_DIR" --sample_bench_dir "$SAMPLE_BENCH_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT