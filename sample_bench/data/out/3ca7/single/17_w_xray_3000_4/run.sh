#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=12:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N md_17
#$ -t 1-25
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/md_3ca7.py .

JOB_ID="$1"
W_XRAY="$2"
DYN_W_XRAY="$3"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/single_md/17_w_xray_3000_4/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
CIF_FILE="$HOME/xray/data/reflections/3ca7/3ca7.cif"
RES=4

START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_clean.pdb"
REF_PDB_FILE="$HOME/xray/data/pdbs/3ca7/3ca7_clean.pdb"
T=3000
STEPS=-1
LOG_FILE="$OUT_DIR/log.csv"

python md_3ca7.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_file "$CIF_FILE" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --start_pdb_file "$START_PDB_FILE" --ref_pdb_file "$REF_PDB_FILE" --T "$T" --steps "$STEPS" --log_file "$LOG_FILE"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT