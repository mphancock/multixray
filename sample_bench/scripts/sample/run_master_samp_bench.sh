#! /bin/bash


JOB_NAME="100_natives_4x"
EXP_ID=100
N_STATE=4
PDB_DIR="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/pdbs/4_state"
CIF_DIR="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/cifs/4_state"
H_RT="24:00:00"
N_JOBS="1-1000"
OFFSET="0"
INIT_WEIGHTS="rand"
SA="{step1000,T300,dofA,pdb1,w1,res0};{step1000,T1000,dofS,pdb0,w0,res-1};{step1000,T1000,dofS,pdb0,w0,res3}"

START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
W_XRAY=1.0
DYN_W_XRAY=1

for JOB_ID in {0..9}
do
    CIF_FILES="$CIF_DIR/$JOB_ID.cif"
    REF_PDB_FILE="$PDB_DIR/$JOB_ID.pdb"

    PARAMS="--cif_files $CIF_FILES --w_xray $W_XRAY --dyn_w_xray --start_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights "$INIT_WEIGHTS" --ref_pdb_file $REF_PDB_FILE --sa $SA"

    echo "$JOB_ID"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID"
    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"

    # echo "$PARAMS"
done

# W_XRAY=$( echo "$JOB_ID*.1" | bc )