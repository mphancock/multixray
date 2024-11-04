PDB_FILE="$HOME/xray/data/pdbs/3k0m/3k0m_clean.pdb"
# CIF_FILE="/home/matthew/xray/data/cifs/7mhf/7mhf.cif"
CIF_FILE="$HOME/xray/dev/45_synthetic_native_4/data/cifs/native_1.cif"
OUT_DIR="$HOME/xray/dev/26_phenix_refine/data/tmp"

cd "$OUT_DIR" && rm *

# phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=individual_sites ordered_solvent=true ordered_solvent.mode=every_macro_cycle refinement.input.xray_data.labels="_refln.F_meas_au,_refln.F_meas_sigma_au" refinement.input.xray_data.r_free_flags.label="_refln.status" write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
# phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=individual_sites ordered_solvent=true ordered_solvent.mode=every_macro_cycle strategy=individual_adp adp.individual.anisotropic="not element H" adp.individual.isotropic="element H" refinement.input.xray_data.labels="_refln.intensity_meas,_refln.intensity_sigma" refinement.input.xray_data.r_free_flags.label="_refln.status" write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=individual_sites refinement.input.xray_data.labels="_refln.F_meas_au,_refln.F_meas_sigma_au" simulated_annealing=true

cd ~/xray/dev/26_phenix_refine/scripts