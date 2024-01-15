#! /bin/bash


JOB_NAME="153_native_2_cif"
EXP_ID=153
N_STATE=2
H_RT="3:00:00"
N_JOBS="1-1000"
OFFSET="0"
INIT_WEIGHTS="rand"
SA="{step99999,T300,dofA,pdb1,w1,res2}"

START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
W_XRAY=1.0

for JOB_ID in {0..9}
do
    echo "$JOB_ID"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID"

    PARAMS="--cif_files "/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/0/$JOB_ID.cif,/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/1/$JOB_ID.cif" --w_xray $W_XRAY --dyn_w_xray --start_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights "$INIT_WEIGHTS" --ref_pdb_file "/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/scores/natives.csv" --sa $SA --bfactor 15"

    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
done
