#! /bin/bash


EXP_ID=187
EXP_NAME="187_bench"
H_RT="06:00:00"
N_JOBS="1-250"
OFFSET="0"

for JOB_ID in {0..39}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_ID"

    PARAMS="--job_csv_file /wynton/home/sali/mhancock/xray/sample_bench/data/params/bench.csv --job_id $JOB_ID"

    qsub -N b"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_DIR" "$OFFSET" "$PARAMS"
done