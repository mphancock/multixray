#!/bin/bash


PDB_DIR="$HOME/xray/data/pdbs"
PDB_FILE="$PDB_DIR/3ca7/remove_S_side_chains_3ca7_clean.pdb"

ERR_PCT=0

MOD_CIF_FILE_NAME="3ca7/remove_S_side_chains_3ca7_clean.cif"
EXP_MTZ_FILE_NAME="3ca7.mtz"

MAPS_DIR="$HOME/xray/data/maps"
REF_DIR="$HOME/xray/data/reflections"
EXP_MTZ_FILE="$MAPS_DIR/$EXP_MTZ_FILE_NAME"
MOD_MTZ_FILE=./"$MOD_CIF_FILE_NAME".tmp.mtz
MOD_CIF_FILE="$REF_DIR/$MOD_CIF_FILE_NAME"
rm "$MOD_CIF_FILE"

phenix.fmodel "$PDB_FILE" "$EXP_MTZ_FILE" add_sigmas=True add_random_error_to_amplitudes_percent="$ERR_PCT" data_column_label="FP,SIGFP" type=real file_name="$MOD_MTZ_FILE"

phenix.mtz_as_cif "$MOD_MTZ_FILE" output_file="$MOD_CIF_FILE" mtz_labels="FMODEL SIGFMODEL" cif_labels="_refln.F_meas_au _refln.F_meas_sigma_au"

rm "$MOD_MTZ_FILE"