#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=9:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N w_xray_rmsd_7mhk
#$ -pe smp 50
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx

python w_xray_rmsd_analysis.py