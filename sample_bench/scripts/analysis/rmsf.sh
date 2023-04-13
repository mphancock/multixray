#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.o
#$ -j y
#$ -l h_rt=1:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N rmsf
#$ -pe smp 50
#$ -l hostname='qb3-id*'


module load CBI conda-stage
conda activate imp_218

JOB_NAME="07_3000_4_exp"
JOB_NUM="890914"
N="10"

echo "JOB_NAME: $JOB_NAME"
echo "JOB_NUM: $JOB_NUM"
echo "N: $N"

python rmsf.py --job_name "$JOB_NAME" --job_num "$JOB_NUM" --n "$N"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT