#! /bin/bash


EXP_ID=202
EXP_NAME="202_no_wxray_auto"
N_JOBS="1-10"
OFFSET="0"
H_RT="24:00:00"

for JOB_ID in {0..13}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_ID"
    JOB_FILE="/wynton/home/sali/mhancock/xray/sample_bench/data/params/$EXP_ID.csv"

    PARAMS="--job_csv_file $JOB_FILE --job_id $JOB_ID"

    qsub -N f"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_DIR" "$OFFSET" "$PARAMS"
done