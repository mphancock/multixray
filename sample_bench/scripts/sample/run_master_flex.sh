#! /bin/bash


EXP_ID=271
EXP_NAME="271_native_2_wxray"
N_JOBS="1-25"
OFFSET="0"
H_RT="24:00:00"
JOB_FILE="/wynton/home/sali/mhancock/xray/sample_bench/data/params/271.csv"


for JOB_ID in {0..17}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$EXP_NAME/$JOB_ID"

    PARAMS="--job_csv_file $JOB_FILE --job_id $JOB_ID"

    qsub -N f"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_DIR" "$OFFSET" "$PARAMS"
done