#! /bin/bash


EXP_ID=170
H_RT="28:00:00"
N_JOBS="1-500"
OFFSET="0"
N_STATE=16

for JOB_ID in {0..62}
do
    JOB_NAME="$EXP_ID"_N"$N_STATE"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"
    echo "$JOB_DIR"
    PARAMS="--input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id $JOB_ID --w_xray 0.5 --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state $N_STATE --init_weights rand --sa {step3000,T300,dofA,pdb1,w1,res0} --steps 2"

    qsub -N slav"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
done