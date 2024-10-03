PDB_FILE="/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/2_state_occs_tmp.pdb"
CIF_FILE="/wynton/home/sali/mhancock/xray/data/cifs/3k0m/3k0m.cif"
OUT_DIR="/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/tmp"

cd "$OUT_DIR" && rm *

phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=none ordered_solvent=true ordered_solvent.mode=every_macro_cycle refinement.input.xray_data.labels="_refln.intensity_meas,_refln.intensity_sigma" refinement.input.xray_data.r_free_flags.label="_refln.status" write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false write_mtz_file=false

cd ~/xray/dev/26_phenix_refine/scripts