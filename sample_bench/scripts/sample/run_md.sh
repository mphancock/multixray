#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=24:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N md_36_7mhk
#$ -t 1-1000
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx


JOB_NAME="36_iso_b_16"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhk/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

CIF_FILES="/wynton/home/sali/mhancock/xray/data/reflections/7mhf/7mhf.cif"
RES=0
W_XRAY=.5
DYN_W_XRAY=1
COM=os
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_no_alt_H_ion.pdb"
N_STATE=16
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_no_alt_H_ion.pdb"
T=300
SA=0
LOG_FILE="$OUT_DIR/log.csv"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files "$CIF_FILES" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --com "$COM" --start_pdb_file "$START_PDB_FILE" --n_state "$N_STATE" --ref_pdb_file "$REF_PDB_FILE" --T "$T" --sa "$SA" --log_file "$LOG_FILE"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT