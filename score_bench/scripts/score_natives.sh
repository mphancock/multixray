#!/bin/bash


DECOY_DIR_NAME="100_natives_4x"
JOB_NAME="100_natives_4x"

NATIVE_DIR="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/pdbs/4_state"
CIF_DIR="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/cifs/4_state"
SCORE_DIR="/wynton/home/sali/mhancock/xray/score_bench/data/3ca7/$JOB_NAME"

MIN_RES=0
SCORE_FS="ml,rmsd_avg,ff"

TIME="1:15:00"
for NATIVE_ID in {0..9}
# Natives 2x, 4x that scored poorly
# for NATIVE_ID in 3 4 10 11 18 19 30 34
# for NATIVE_ID in 1 3 8 9 14 18 25 29 30 34 36
# Natives 4x that scored very poorly
# for NATIVE_ID in 29
do
    DECOY_DIR="/wynton/group/sali/mhancock/xray/decoys/data/3ca7/$DECOY_DIR_NAME/$NATIVE_ID"

    SCORE_FILE_NAME="$NATIVE_ID"
    REF_PDB_FILE="$NATIVE_DIR"/"$NATIVE_ID".pdb
    CIF_FILE="$CIF_DIR"/"$NATIVE_ID".cif

    echo "$SCORE_FILE_NAME"
    echo "$REF_PDB_FILE"
    echo "$CIF_FILE"

    SCORE_FILE="$SCORE_DIR/$SCORE_FILE_NAME.csv"
    PARAMS_FILE="$SCORE_DIR/$SCORE_FILE_NAME.txt"
    rm "$SCORE_FILE"
    rm "$PARAMS_FILE"

    PARAMS="--pdb_dir $DECOY_DIR --cif_file $CIF_FILE --min_res $MIN_RES --ref_pdb_file $REF_PDB_FILE --param_file $PARAMS_FILE --score_file $SCORE_FILE --score_fs $SCORE_FS --add_native"

    # Not sure if I should run these scripts on Wynton compute nodes. The variance in wall clock is so high so it makes resource allocation challenging.
    # qsub -N score_"$NATIVE_ID" -l h_rt="$TIME" ~/xray/score_bench/scripts/score_natives_single.sh "$PARAMS"
    source ~/xray/score_bench/scripts/score_natives_single.sh "$PARAMS"
done