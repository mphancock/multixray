#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=06:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N w_xray_7mhk
#$ -t 1-50
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx


JOB_NAME="$1"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhk/$JOB_NAME/$2"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

CIF_FILE="$HOME/xray/data/reflections/7mhk/7mhk.cif"
RES=$5
W_XRAY=$3
DYN_W_XRAY=$4
COM=os
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhk/7mhk_clean.pdb"
N_STATE=1
REF_PDB_FILE="$HOME/xray/data/pdbs/7mhk/7mhk_clean.pdb"
UC_DIM="114.300 54.290 44.970 90.00 102.12 90.00"
SG_SYMBOL="C 1 2 1"
T=$6
SA=$7
LOG_FILE="$OUT_DIR/log.csv"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_file "$CIF_FILE" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --com "$COM" --start_pdb_file "$START_PDB_FILE" --n_state "$N_STATE" --ref_pdb_file "$REF_PDB_FILE" --uc_dim "$UC_DIM" --sg_symbol="$SG_SYMBOL" --T "$T" --sa "$SA" --log_file "$LOG_FILE"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT