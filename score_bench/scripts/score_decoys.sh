#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/xray/score_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -N score
#$ -j y
#$ -l h_rt=24:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -pe smp 1
#$ -t 1-1
#$ -l hostname='qb3-id*'

eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_221_cctbx


python /wynton/home/sali/mhancock/xray/score_bench/scripts/score_decoys_xray.py --decoy_file ~/xray/score_bench/data/266_decoys/rand1000_rand_w.csv --N 2 --J 2

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"

