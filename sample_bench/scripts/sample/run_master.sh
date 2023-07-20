#! /bin/bash


JOB_NAME="64_native_4x"
EXP_ID=64
N_STATE=4
PDB_DIR="$HOME/xray/dev/17_synthetic_native/data/pdbs/4_state_ref"
CIF_DIR="$HOME/xray/dev/17_synthetic_native/data/cifs/4_state_ref"

for JOB_ID in {0..39}
do
    CIF_FILES="$CIF_DIR/$JOB_ID.cif"
    RES=0
    W_XRAY=1.0
    DYN_W_XRAY=1
    START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
    REF_PDB_FILE="$PDB_DIR/$JOB_ID.pdb"
    T=300

    PARAMS="--cif_files $CIF_FILES --res $RES --w_xray $W_XRAY --dyn_w_xray --start_pdb_file $START_PDB_FILE --n_state $N_STATE --ref_pdb_file $REF_PDB_FILE --T $T --weights"

    echo "$JOB_ID"
    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt="03:00:00" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_ID" "$PARAMS"
done

# W_XRAY=$( echo "$JOB_ID*.1" | bc )
# DYN_W_XRAY=1

# echo "$JOB_ID, $W_XRAY, $DYN_W_XRAY"

# qsub "$HOME/xray/sample_bench/scripts/sample/run_w_xray.sh" "$JOB_NAME" "$JOB_ID" "$W_XRAY" "$DYN_W_XRAY" "$RES" "$T" "$SA" "$N_STATE"