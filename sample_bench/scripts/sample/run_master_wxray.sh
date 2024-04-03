#! /bin/bash


JOB_NAME="175_w_xray_N16"
EXP_ID=175
H_RT="08:00:00"
N_JOBS="1-100"
OFFSET="0"
N_STATE=16

W_XRAYS=(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, .6, .7, .8, .9, 1.0)

for JOB_ID in {0..10}
do
    echo "$JOB_ID"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"

    W_XRAY=${W_XRAYS[$JOB_ID]}

    PARAMS="--input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id 0 --w_xray $W_XRAY --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state $N_STATE --init_weights rand --sa {step3000,T300,dofA,pdb1,w1,res0} --steps 2"

    qsub -N slav"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
done
