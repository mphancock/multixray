#!/bin/bash

PDB_FILE=~/xray/dev/26_phenix_refine/data/0/0_refine_001.pdb
CIF_FILE=$HOME/xray/data/cifs/3k0m/3k0m.cif

# phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.F_meas_au" r_free_flags_label="_refln.status"
phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.intensity" r_free_flags_label="_refln.status"