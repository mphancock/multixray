#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -t 1-50
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx


JOB_NAME="$1"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$2"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"

# Delete the job dir
rm -r "$OUT_DIR"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .

# Arguement $3 is referenced without quotes to prevent it from being treated as a single string.
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" $3

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT