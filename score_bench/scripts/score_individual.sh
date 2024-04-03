#!/bin/bash

# PDB_FILE=/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/143_native_1_cif_2/0/output_88/pdbs/2155.pdb
PDB_FILE=/wynton/home/sali/mhancock/xray/tmp/tmp.pdb
CIF_FILE=/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhi.cif

# phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" high_resolution=2.0
phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.F_meas_au" r_free_flags_label="_refln.status"