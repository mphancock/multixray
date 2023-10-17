#$ -N dev_23
#$ -l h_rt=00:30:00
#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/xray/dev/tmp
#$ -j y
#$ -l mem_free=48G
#$ -l scratch=1G
#$ -pe smp 1
#$ -t 1-100
#$ -l hostname='qb3-id*'

eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx

python ~/xray/dev/23_crystal_sim/scripts/hetero_vs_degenerate.py --trial_id "$SGE_TASK_ID" --max_N 100
# python ~/xray/dev/23_crystal_sim/scripts/hetero_vs_degenerate.py --trial_id 0 --max_N 50
