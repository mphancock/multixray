#! /bin/bash


EXP_ID=186
EXP_NAME="186_w_xray_bench"
N_JOBS="1-5"
OFFSET="0"
H_RT="6:00:00"

for JOB_ID in {0..439}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_ID"

    PARAMS="--job_csv_file /wynton/home/sali/mhancock/xray/sample_bench/data/params/w_xray_bench.csv --job_id $JOB_ID"

    qsub -N w"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_DIR" "$OFFSET" "$PARAMS"
done