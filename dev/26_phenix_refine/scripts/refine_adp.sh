#!/bin/bash


PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7.pdb"
CIF_FILE="/wynton/home/sali/mhancock/xray/data/cifs/3ca7/3ca7.cif"

OUT_DIR="/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/3ca7"

cd "$OUT_DIR"

phenix.refine "$CIF_FILE" "$PDB_FILE" strategy=tls base_output_dir="$OUT_DIR" write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false

cd "/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/scripts"