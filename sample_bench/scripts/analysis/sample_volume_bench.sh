#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.o
#$ -j y
#$ -l h_rt=1:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N vol_bench
#$ -pe smp 25
#$ -l hostname='qb3-id*'


# module load CBI conda-stage
# conda activate imp_218_cctbx

JOB_NAME="46_synth_sa_3"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine.pdb"
JOB_ID="9417275_xray"
MAX_N="1"
FIELD="xray_0"
BONUS_FIELDS="pdb,rmsd"
# FIELD="rmsd"
# BONUS_FIELDS="pdb"

JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"
SAMPLE_BENCH_DIR="/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/$JOB_NAME"
# OUT_FILE="/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/$JOB_NAME/sample_bench_xray.csv"

python sample_volume_bench.py --job_dir "$JOB_DIR" --sample_bench_dir "$SAMPLE_BENCH_DIR" --ref_pdb_file "$REF_PDB_FILE" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --max_n "$MAX_N"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT