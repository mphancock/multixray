#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=04:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N score
#$ -pe smp 1
#$ -t 1-1
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_220_cctbx

python ~/xray/sample_bench/scripts/analysis_exp/score_pdbs.py --pdb_df_file ~/xray/sample_bench/data/7mhf/169_N8/ref_15.csv --skip 0