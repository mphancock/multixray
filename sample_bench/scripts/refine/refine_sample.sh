#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=01:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N refine
#$ -pe smp 1
#$ -t 1-10
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp_220_cctbx

# SGE_TASK_ID=48
EXP_NAME="182_bench"
START=$((($SGE_TASK_ID - 1)*1000))
N_PDB=1000
MAX_FF=10000

echo $START, $N_PDB

python /wynton/home/sali/mhancock/xray/sample_bench/scripts/refine/refine_sample.py --sample_file /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/182_bench/sample_per_out.csv --max_ff $MAX_FF --start $START --n_pdb $N_PDB

python /wynton/home/sali/mhancock/xray/sample_bench/scripts/analysis_bench/score_pdbs.py --sample_file /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/182_bench/ref_"$MAX_FF"_"$START".csv

# python refine_sample.py --exp_name $EXP_NAME --max_ff $MAX_FF --start $START --n_pdb $N_PDB

# python ~/xray/sample_bench/scripts/analysis_exp/score_pdbs.py --sample_file /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/$EXP_NAME/ref_"$MAX_FF"_"$START".csv