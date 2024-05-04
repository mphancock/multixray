#! /bin/bash


EXP_ID=180
EXP_NAME="180_wxray_bench"
N_JOBS="1-10"
OFFSET="0"

H_RT="6:00:00"

for n in {0..1}
do
    for j in {0..19}
    do
        for w in {0..10}
        do
            # echo $W_XRAY_ID
            JOB_NAME=n"$n"_j"$j"_w"$w"
            echo "$JOB_NAME"
            JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_NAME"

            qsub -N r"$JOB_NAME" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/refine/refine_slave.sh" "$JOB_NAME" "$JOB_DIR"
        done
    done
done