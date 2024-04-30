#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=0:30:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N refine
#$ -pe smp 1
#$ -t 1-90
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_220_cctbx

# SGE_TASK_ID=48
EXP_NAME="181_wxray_control"
START=$((($SGE_TASK_ID - 1)*100))
N_PDB=100
MAX_FF=10000

echo $START, $N_PDB

python refine_sample.py --exp_name $EXP_NAME --max_ff $MAX_FF --start $START --n_pdb $N_PDB

python ~/xray/sample_bench/scripts/analysis_exp/score_pdbs.py --sample_file /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/$EXP_NAME/ref_"$MAX_FF"_"$START".csv