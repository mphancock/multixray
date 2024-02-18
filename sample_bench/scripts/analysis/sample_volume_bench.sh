#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=00:20:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N vol_bench
#$ -pe smp 8
#$ -t 1-1
#$ -l hostname='qb3-id*'

#!/bin/bash
module load CBI conda-stage
conda activate imp_219_cctbx

JOB_NAME="124_natives_2_cond"
SYSTEM="7mhf"
SAMPLE_BENCH_DIR="/wynton/home/sali/mhancock/xray/sample_bench/data/$SYSTEM/$JOB_NAME"

mkdir -p "$SAMPLE_BENCH_DIR"

for JOB_ID in {0..9}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$SYSTEM/$JOB_NAME/$JOB_ID"

    # FIELD="xray_0"
    # BONUS_FIELDS="rmsd_0"
    FIELD="xray_0+xray_1"
    BONUS_FIELDS="rmsd_0+rmsd_1,pdb"

    FILE="$SAMPLE_BENCH_DIR/xray_volume_bench_$JOB_ID.csv"
    python sample_volume_bench.py --job_dir "$JOB_DIR" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --file "$FILE"

    # FIELD="rmsd_0"
    FIELD="rmsd_0+rmsd_1"
    FILE="$SAMPLE_BENCH_DIR/rmsd_volume_bench_$JOB_ID.csv"
    BONUS_FIELDS="pdb"
    python sample_volume_bench.py --job_dir "$JOB_DIR" --field "$FIELD" --bonus_fields "$BONUS_FIELDS" --file "$FILE"
done