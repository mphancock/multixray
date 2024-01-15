# Setup the conda environment.
OUT_DIR="/wynton/home/sali/mhancock/xray/dev/30_ff_decoys/data/S/output_0"
rm -r "$OUT_DIR"
mkdir -p "$OUT_DIR"

START_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_0/0.pdb"
N_STATE=2
REF_PDB_FILE="$START_PDB_FILE"
SA="{step1000,T300,dofS,pdb0,w0,res-1}"
LOG_FILE="$OUT_DIR/log.csv"
INIT_WEIGHTS="uni"
STEPS=1000

python /wynton/home/sali/mhancock/xray/sample_bench/scripts/sample/run_md_multi.py --out_dir "$OUT_DIR" --start_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights "$INIT_WEIGHTS" --ref_pdb_file $REF_PDB_FILE --sa "$SA" --bfactor 15 --steps "$STEPS"