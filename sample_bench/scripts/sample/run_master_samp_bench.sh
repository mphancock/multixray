#! /bin/bash


JOB_NAME="118_2_state_2_cif_no_weight"
EXP_ID=118
N_STATE=2
H_RT="12:00:00"
N_JOBS="1-500"
OFFSET="0"
INIT_WEIGHTS="ref"
SA="{step1000,T300,dofA,pdb1,w0,res0};{step1000,T1000,dofS,pdb0,w0,res-1};{step1000,T1000,dofS,pdb0,w0,res3}"

START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
W_XRAY=1.0
DYN_W_XRAY=1

for JOB_ID in {0..9}
do
    CIF_FILES="/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/2_state_0/$JOB_ID.cif,/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/2_state_1/$JOB_ID.cif"
    REF_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_0/$JOB_ID.pdb,/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_1/$JOB_ID.pdb"

    PARAMS="--cif_files $CIF_FILES --w_xray $W_XRAY --dyn_w_xray --start_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights "$INIT_WEIGHTS" --ref_pdb_file $REF_PDB_FILE --sa $SA --bfactor 15"

    echo "$JOB_ID"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID"
    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
done
