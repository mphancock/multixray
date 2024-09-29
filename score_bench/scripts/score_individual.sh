#!/bin/bash

PDB_FILE=~/xray/data/pdbs/3k0m/3k0m.pdb
CIF_FILE=$HOME/xray/data/cifs/3k0m/3k0n.cif

# phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.F_meas_au" r_free_flags_label="_refln.status"
phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.intensity" r_free_flags_label="_refln.status"