#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=00:30:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N score
#$ -pe smp 8
#$ -t 1-1
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_219_cctbx

python ~/xray/sample_bench/scripts/analysis/score_analysis.py --n_cond 3 --job_name 162_J3_ijk