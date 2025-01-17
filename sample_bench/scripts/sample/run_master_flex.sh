#! /bin/bash


EXP_ID=273
EXP_NAME="273_native_3_wxray"
N_JOBS="1-1000"
OFFSET="0"
H_RT="12:00:00"
JOB_FILE="/wynton/home/sali/mhancock/xray/sample_bench/data/params/$EXP_ID.csv"


for JOB_ID in {0..8}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$EXP_NAME/$JOB_ID"

    PARAMS="--job_csv_file $JOB_FILE --job_id $JOB_ID"

    qsub -N f"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_DIR" "$OFFSET" "$PARAMS"
done