#! /bin/bash


EXP_ID=261
EXP_NAME="261_2_state_test"
N_JOBS="1-10"
OFFSET="0"
H_RT="24:00:00"
JOB_FILE="/wynton/home/sali/mhancock/xray/sample_bench/data/params/261.csv"
# JOB_FILE="/wynton/home/sali/mhancock/xray/sample_bench/data/params/$EXP_ID.csv"


for JOB_ID in {0..19}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$EXP_NAME/$JOB_ID"

    PARAMS="--job_csv_file $JOB_FILE --job_id $JOB_ID"

    qsub -N f"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_DIR" "$OFFSET" "$PARAMS"
done