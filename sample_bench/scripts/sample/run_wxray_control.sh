#! /bin/bash


EXP_ID=181
EXP_NAME="181_wxray_control"
N_JOBS="1-25"
OFFSET="0"

H_RTS=("08:00:00" "12:00:00" "16:00:00" "24:00:00" "56:00:00")
N_STATES=(1 2 4 8 16)
W_XRAY_MULTS=(1 2 3 4 5 6)

for N_STATE_ID in {0..4}
do
    N_STATE=${N_STATES[$N_STATE_ID]}
    H_RT=${H_RTS[$N_STATE_ID]}
    for JOB_ID in {0..5}
    do
        for W_XRAY_MULT_ID in {0..5}
        do
            W_XRAY_MULT=${W_XRAY_MULTS[$W_XRAY_MULT_ID]}

            JOB_NAME=N"$N_STATE_ID"_J"$JOB_ID"_M"$W_XRAY_MULT_ID"
            echo "$JOB_NAME"
            JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_NAME"

            PARAMS="--input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id $JOB_ID --w_xray /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/178_wxray/best_wxray.csv --w_xray_multiplier $W_XRAY_MULT --n_state $N_STATE --init_weights rand --sa {step3000,T300,dofA,pdb1,w1,res0} --steps 2"

            echo $PARAMS
            qsub -N c"$JOB_NAME" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
        done
    done
done