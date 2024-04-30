#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=00:30:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N score
#$ -pe smp 1
#$ -t 1-480
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_220_cctbx

START=$((($SGE_TASK_ID - 1)*100))

python ~/xray/sample_bench/scripts/analysis_exp/score_pdbs.py --sample_file /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/179_exp/ref_25000_$START.csv