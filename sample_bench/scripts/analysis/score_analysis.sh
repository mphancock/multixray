#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=01:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N score_analysis
#$ -pe smp 8
#$ -t 1-1
#$ -l hostname='qb3-id*'

#!/bin/bash
module load CBI conda-stage
conda activate imp_219_cctbx

python ~/xray/sample_bench/scripts/analysis/score_analysis.py