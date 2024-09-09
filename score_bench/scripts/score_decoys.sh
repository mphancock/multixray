#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/xray/score_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -N score
#$ -j y
#$ -l h_rt=02:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -pe smp 1
#$ -t 1-1
#$ -l hostname='qb3-id*'

eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_220_cctbx


python ~/xray/score_bench/scripts/score_decoys.py --job 122_native_decoys_1_state

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"

