#!/bin/bash

PDB_FILE=~/xray/tmp/6_500_align.pdb
CIF_FILE=$HOME/xray/dev/39_bench_ensemble/data/cifs/7mhl.cif

# phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" high_resolution=2.0
phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" f_obs_label="_refln.F_meas_au" r_free_flags_label="_refln.status"