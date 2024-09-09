#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=8:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N c194
#$ -t 1-25
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_220_cctbx


JOB_NAME="194_bench_free_control"
TARGET="7mhf"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$TARGET/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --job_csv_file /wynton/home/sali/mhancock/xray/sample_bench/data/params/194_bench_free_control.csv --job_id 0

# If run_md_multi.py throws an error, then the contents will be copied (for debugging purposes)
cp -r "$TMP_OUT_DIR/*" "$OUT_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT