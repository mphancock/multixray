#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=06:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N w_xray
#$ -t 1-50
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx


JOB_NAME="$1"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$2"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
rm -r "$OUT_DIR"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

CIF_FILE="$HOME/xray/data/reflections/3ca7/3ca7.cif"
RES="$5"
W_XRAY=$3
DYN_W_XRAY="$4"
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_clean_x16.pdb"
REF_PDB_FILE="$HOME/xray/data/pdbs/3ca7/3ca7_clean.pdb"
T="$6"
LOG_FILE="$OUT_DIR/log.csv"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/md_3ca7_multi.py .
python md_3ca7_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_file "$CIF_FILE" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --start_pdb_file "$START_PDB_FILE" --ref_pdb_file "$REF_PDB_FILE" --T "$T" --log_file "$LOG_FILE"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT