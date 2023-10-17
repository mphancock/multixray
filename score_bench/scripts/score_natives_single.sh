#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/xray/score_bench/tmp
#$ -N score
#$ -j y
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -pe smp 48
#$ -l hostname='qb3-id*'

# eval "$(conda shell.bash hook)"
# module load CBI conda-stage
# conda activate imp_218_cctbx

python ~/xray/score_bench/scripts/score_pdb_dir.py $1