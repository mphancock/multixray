#!/bin/bash

PDB_FILE=~/xray/dev/42_traj_analysis/data/242/pdbs/11/1.pdb
CIF_FILE=$HOME/xray/data/cifs/3k0m/3k0m.cif

# phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.F_meas_au" r_free_flags_label="_refln.status"
phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.intensity" r_free_flags_label="_refln.status"