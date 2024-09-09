#! /bin/bash


EXP_ID=188
EXP_NAME="188_w_xray_exp"
N_JOBS="1-5"
OFFSET="0"

H_RT="24:00:00"

for JOB_ID in {0..3464}
do
    JOB_NAME="$JOB_ID"
    echo "$JOB_NAME"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_NAME"

    PARAMS="--job_csv_file /wynton/home/sali/mhancock/xray/sample_bench/data/params/w_xray_exp.csv --max_ff 10000"

    qsub -N r"$EXP_ID"_"$JOB_NAME" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/refine/refine_slave.sh" "$JOB_DIR" 0 "$PARAMS"
done