#!/bin/bash


PDB_FILE="/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_0/0.pdb"
TMP_DIR="/wynton/home/sali/mhancock/xray/tmp"
MOD_CIF_FILE="/wynton/home/sali/mhancock/xray/tmp/0_2.cif"
EXP_MTZ_FILE="/wynton/home/sali/mhancock/xray/data/maps/3ca7.mtz"
MOD_MTZ_FILE=/wynton/home/sali/mhancock/xray/tmp/0_2.mtz

ERR_PCT=0

# MAPS_DIR="$HOME/xray/data/maps"
# REF_DIR="$HOME/xray/data/reflections"
# rm "$MOD_CIF_FILE"

# phenix.fmodel "$PDB_FILE" "$EXP_MTZ_FILE" add_sigmas=True add_random_error_to_amplitudes_percent="$ERR_PCT" data_column_label="FP,SIGFP" type=real file_name="$MOD_MTZ_FILE"
phenix.fmodel "$PDB_FILE" add_sigmas=True add_random_error_to_amplitudes_percent="$ERR_PCT" data_column_label="FP,SIGFP" type=real file_name="$MOD_MTZ_FILE" high_resolution=2.0

phenix.mtz_as_cif "$MOD_MTZ_FILE" output_file="$MOD_CIF_FILE" mtz_labels="FMODEL SIGFMODEL" cif_labels="_refln.F_meas_au _refln.F_meas_sigma_au"

rm "$MOD_MTZ_FILE"