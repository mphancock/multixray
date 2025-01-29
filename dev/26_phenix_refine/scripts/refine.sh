PDB_FILE="$HOME/xray/tmp/tmp.updated.pdb"
CIF_FILE="$HOME/xray/dev/38_standard_flags/data/7mhh.cif"

# PDB_FILE="$HOME/xray/dev/26_phenix_refine/data/3k0m_4_state_mod.pdb"
# CIF_FILE="$HOME/xray/data/cifs/3k0m/3k0m.cif"
OUT_DIR="$HOME/xray/dev/26_phenix_refine/data"

cd "$OUT_DIR" && rm *

# phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=individual_sites ordered_solvent=true ordered_solvent.mode=every_macro_cycle refinement.input.xray_data.labels="_refln.F_meas_au,_refln.F_meas_sigma_au" refinement.input.xray_data.r_free_flags.label="_refln.status" write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
# phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=individual_sites ordered_solvent=true ordered_solvent.mode=every_macro_cycle strategy=individual_adp adp.individual.anisotropic="not element H" adp.individual.isotropic="element H" refinement.input.xray_data.labels="_refln.intensity_meas,_refln.intensity_sigma" refinement.input.xray_data.r_free_flags.label="_refln.status" write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
# phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=individual_sites refinement.input.xray_data.labels="_refln.F_meas_au,_refln.F_meas_sigma_au" simulated_annealing=true refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" strategy=individual_sites+individual_adp ordered_solvent=true ordered_solvent.mode=every_macro_cycle  refinement.input.xray_data.labels="_refln.F_meas_au,_refln.F_meas_sigma_au" refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None crystal_symmetry.unit_cell=115.023,54.358,44.970,90.00,101.50,90.00 crystal_symmetry.space_group="C 1 2 1"

cd ~/xray/dev/26_phenix_refine/scripts