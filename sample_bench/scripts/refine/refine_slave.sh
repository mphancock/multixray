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
conda activate imp_220_cctbx


JOB_DIR="$1"
OFFSET="$2"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1+OFFSET))
OUT_DIR="$JOB_DIR/output_$RUN_ID"

python ~/xray/sample_bench/scripts/refine/refine_all_models.py --out_dir "$OUT_DIR" $3
# python refine_output.py --out_dir "$OUT_DIR" $3

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT