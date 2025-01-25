#!/bin/bash

PDB_FILE="/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/tmp/411_mod_refine_001.pdb"
CIF_FILE="$HOME/xray/dev/38_standard_flags/data/7mhh.cif"

# phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.F_meas_au" r_free_flags_label="_refln.status"
phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.F_meas_au" r_free_flags_label="_refln.status"