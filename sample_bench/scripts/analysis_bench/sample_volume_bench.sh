#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=00:20:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N volume
#$ -pe smp 8
#$ -t 1-1
#$ -l hostname='qb3-id*'

#!/bin/bash
module load CBI conda-stage
conda activate imp_220_cctbx

JOB_NAME="152_native_N8_J2"
N_COND=2

SAMPLE_BENCH_DIR="/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/$JOB_NAME"
mkdir -p "$SAMPLE_BENCH_DIR"


for JOB_ID in {0..9}
# for JOB_ID in 0
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"

    FIELD="xray_0+xray_1"
    echo $FIELD
    FILE="$SAMPLE_BENCH_DIR"/volume_"$FIELD"_"$JOB_ID".csv
    python sample_volume_bench.py --job_dir "$JOB_DIR" --field "$FIELD" --bonus_fields  "$BONUS_FIELDS" --file "$FILE"

    FIELD="rmsd"
    fi
    BONUS_FIELDS="pdb"

    echo $FIELD
    FILE="$SAMPLE_BENCH_DIR"/volume_"$FIELD"_"$JOB_ID".csv
    python sample_volume_bench.py --job_dir "$JOB_DIR" --field "$FIELD" --bonus_fields "$BONUS_FIELDS" --file "$FILE"
done