#! /bin/bash


JOB_NAME="82_native_4x_decoys"
EXP_ID=82
N_STATE=4

PDB_DIR=/wynton/home/sali/mhancock/xray/dev/17_synthetic_native/data/pdbs/"$N_STATE"_state_ref

INIT_WEIGHTS="ref"
T=100
STEPS=1000
H_RT="00:30:00"
N_JOBS="1-100"

for JOB_ID in {0..39}
do
    START_PDB_FILE="$PDB_DIR/$JOB_ID.pdb"
    REF_PDB_FILE="$PDB_DIR/$JOB_ID.pdb"


    PARAMS="--start_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights $INIT_WEIGHTS --ref_pdb_file $REF_PDB_FILE --T $T --steps $STEPS"

    echo "$JOB_ID"
    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_ID" "$PARAMS"
done

# W_XRAY=$( echo "$JOB_ID*.1" | bc )