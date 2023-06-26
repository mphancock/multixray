#!/bin/bash


JOB_NAME="50_synth_sa_5"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine.pdb"
JOB_ID="9479402"
MAX_N="1000"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"
SAMPLE_BENCH_DIR="/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/$JOB_NAME"


FIELD="xray_0"
BONUS_FIELDS="pdb,rmsd"
# python sample_volume_bench.py --job_dir "$JOB_DIR" --sample_bench_dir "$SAMPLE_BENCH_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --max_n "$MAX_N"

FIELD="rmsd"
BONUS_FIELDS="pdb"
python sample_volume_bench.py --job_dir "$JOB_DIR" --sample_bench_dir "$SAMPLE_BENCH_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --max_n "$MAX_N"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT