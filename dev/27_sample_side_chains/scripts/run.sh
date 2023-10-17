#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=12:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N dev_27
#$ -t 1-100
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_219_cctbx


JOB_DIR="/wynton/group/sali/mhancock/xray/dev/27_sample_side_chains/data/out"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
rm -r "$OUT_DIR"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"


CIF_FILES="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/cifs/4_state/0.cif"
W_XRAY=1.0
START_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/27_sample_side_chains/data/out/output_0/pdbs/101.pdb"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/19_synthetic_native_2/data/pdbs/4_state/0.pdb"
SA="1000,300,S,1,1;1000,1000,S,2,1;1000,2000,S,3,1;1000,3000,S,4,1;1000,2000,S,3,1;1000,1000,S,2,1"
LOG_FILE="$OUT_DIR/log.csv"
INIT_WEIGHTS="rand"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files $CIF_FILES --w_xray $W_XRAY --dyn_w_xray --start_pdb_file $START_PDB_FILE --n_state 1 --init_weights "$INIT_WEIGHTS" --ref_pdb_file $REF_PDB_FILE --weights --sa "$SA"


[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT