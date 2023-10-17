#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_219_cctbx


JOB_NAME="$1"
JOB_DIR="$2"
OFFSET="$4"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1+$OFFSET))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"

# Delete the output dir.
rm -r "$OUT_DIR"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .

# Arguement $3 is referenced without quotes to prevent it from being treated as a single string.
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" $3

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT