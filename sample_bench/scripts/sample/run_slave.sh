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
conda activate imp_221_cctbx


JOB_DIR="$1"
OFFSET="$2"
# mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1+$OFFSET))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"

## PATCH
# ## check if the pdb directory exists
# run=1
# if [ -d "$OUT_DIR/pdbs" ]; then
#     file_count=$(find "$OUT_DIR/pdbs" -maxdepth 1 -type f | wc -l)
#     if [ "$file_count" -gt 250 ]; then
#         run=0
#     fi
# fi

# if [ "$run" -eq 0 ]; then
#     echo "not running"
#     exit 0
# fi
# echo "running"

# Delete the output dir.
# rm -r "$OUT_DIR"
# mkdir "$TMP_OUT_DIR"
# mkdir "$OUT_DIR"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .

# Arguement $3 is referenced without quotes to prevent it from being treated as a single string.
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" $3

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT