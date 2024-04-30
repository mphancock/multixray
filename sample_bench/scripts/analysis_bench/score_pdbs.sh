#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=0:30:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N refine
#$ -pe smp 1
#$ -t 1-1
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_220_cctbx

# SGE_TASK_ID=48
EXP_NAME="182_bench"

python ~/xray/sample_bench/scripts/analysis_bench/score_pdbs.py --sample_file /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/$EXP_NAME/sample.csv