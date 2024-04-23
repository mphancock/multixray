#! /bin/bash


EXP_ID=178
EXP_NAME="178_wxray"
N_JOBS="1-5"
OFFSET="0"

H_RTS=("08:00:00" "12:00:00" "16:00:00" "24:00:00" "56:00:00")
N_STATES=(1 2 4 8 16)
W_XRAYS=(0.03125 0.0625 .125 .25 .5 1 2 4 8 16 32)

for N_STATE_ID in {0..4}
do
    H_RT=${H_RTS[$N_STATE_ID]}
    for JOB_ID in {0..62}
    do
        for W_XRAY_ID in 10
        do
            # echo $W_XRAY_ID
            JOB_NAME=N"$N_STATE_ID"_J"$JOB_ID"_W"$W_XRAY_ID"
            echo "$JOB_NAME"
            JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_NAME"

            W_XRAY=${W_XRAYS[$W_XRAY_ID]}
            N_STATE=${N_STATES[$N_STATE_ID]}

            PARAMS="--input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id $JOB_ID --w_xray /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/178_wxray/best_wxray.csv --n_state $N_STATE --init_weights rand --sa {step3000,T300,dofA,pdb1,w1,res0} --steps 2"

            qsub -N w"$JOB_NAME" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
        done
    done
done