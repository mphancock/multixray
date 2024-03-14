PDB_FILE="/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/n2.pdb"
CIF_FILE="/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhi.cif"
OUT_DIR="/wynton/home/sali/mhancock/xray/dev/26_phenix_refine/data/n2"

cd "$OUT_DIR"
rm *

phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$OUT_DIR" ordered_solvent=true refinement.input.xray_data.labels="_refln.intensity_meas,_refln.intensity_sigma" refinement.input.xray_data.r_free_flags.label="_refln.status"

cd ~/xray/dev/26_phenix_refine/scripts