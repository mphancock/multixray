#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=06:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N refine
#$ -pe smp 1
#$ -t 1-1
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_220_cctbx

python refine_sample.py --exp_name 169_N8 --n_step 15