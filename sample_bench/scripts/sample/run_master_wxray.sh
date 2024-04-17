#! /bin/bash


JOB_NAME="177_wxray_7mhf"
EXP_ID=177
H_RT="08:00:00"
N_JOBS="1-50"
OFFSET="0"

N_STATES=(1 2 4 8 16)
W_XRAYS=(0.5 1.0 1.5 2.0 2.5 3.0)

for N_STATE_ID in {0..4}
do
    for W_XRAY_ID in {0..5}
    do
        JOB_ID=N"$N_STATE_ID"_W"$W_XRAY_ID"
        echo "$JOB_ID"
        JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"

        W_XRAY=${W_XRAYS[$W_XRAY_ID]}
        N_STATE=${N_STATES[$N_STATE_ID]}

        PARAMS="--input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id 0 --w_xray $W_XRAY --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state $N_STATE --init_weights rand --sa {step3000,T300,dofA,pdb1,w1,res0} --steps 2"

        qsub -N w"$EXP_ID""$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
    done
done