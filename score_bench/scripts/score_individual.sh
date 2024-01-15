#!/bin/bash

# PDB_FILE=/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/143_native_1_cif_2/0/output_88/pdbs/2155.pdb
PDB_FILE=/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_0/0.pdb
CIF_FILE=/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/2_state_0_noise/0.cif

phenix.model_vs_data "$PDB_FILE" "$CIF_FILE" high_resolution=2.0