#! /bin/bash


EXP_ID=263
EXP_NAME="263_sb_sa"
N_JOBS="1-10"
OFFSET="0"

H_RT="1:00:00"

for JOB_ID in {0..35}
do
    JOB_NAME="$JOB_ID"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$EXP_NAME/$JOB_NAME"

    PARAMS="--job_csv_file /wynton/home/sali/mhancock/xray/sample_bench/data/params/263.csv --end_only"

    qsub -N r"$EXP_ID"_"$JOB_NAME" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/refine/refine_slave.sh" "$JOB_DIR" 0 "$PARAMS"
done