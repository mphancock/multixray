#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.o
#$ -j y
#$ -l h_rt=1:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N vol_bench
#$ -pe smp 50
#$ -l hostname='qb3-id*'

#
#module load CBI conda-stage
#conda activate imp_218_cctbx

JOB_NAME="18_300_exp"
JOB_ID="1984501"
MAX_N="50"

python sample_volume_bench.py --job_name "$JOB_NAME" --job_id "$JOB_ID" --max_n "$MAX_N"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT