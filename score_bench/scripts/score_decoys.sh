#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/xray/score_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -N score
#$ -j y
#$ -l h_rt=00:20:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -pe smp 8
#$ -t 1-10
#$ -l hostname='qb3-id*'

eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_219_cctbx


I=$((SGE_TASK_ID-1))
python ~/xray/score_bench/scripts/score_decoys.py --i "$I" --n_cond 2 --job 122_native_decoys_1_state

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"

