#!/bin/bash


PDB_FILE="$HOME/xray/dev/21_degeneracy/data/degens/34.pdb"
MTZ_OUT_FILE="$HOME/tmp/34.mtz"
CIF_OUT_FILE="$HOME/tmp/34.cif"


phenix.fmodel "$PDB_FILE" add_sigmas=True generate_fake_p1_symmetry=True high_resolution=1 low_resolution=5 data_column_label="FP,SIGFP" type=real file_name="$MTZ_OUT_FILE"

phenix.mtz_as_cif "$MTZ_OUT_FILE" output_file="$CIF_OUT_FILE" mtz_labels="FMODEL SIGFMODEL" cif_labels="_refln.F_meas_au _refln.F_meas_sigma_au"