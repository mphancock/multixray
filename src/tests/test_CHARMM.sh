#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.o
#$ -j y
#$ -l h_rt=24:00:00
#$ -l mem_free=4G
#$ -l scratch=4G
#$ -N CHARMM
#$ -t 1-1
#$ -l hostname='qb3-id*'


eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_221_cctbx

python ~/xray/src/tests/test_CHARMMRestraint.py

