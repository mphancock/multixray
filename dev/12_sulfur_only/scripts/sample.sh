#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=12:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N dev_14
#$ -t 1-500
#$ -l hostname='qb3-id*'

# Setup the conda environment.
# eval "$(conda shell.bash hook)"
# module load CBI conda-stage
# conda activate imp_218_cctbx

JOB_ID=0
JOB_NAME="14_s_sa"
JOB_DIR="/wynton/group/sali/mhancock/xray/dev/12_sulfur_only/data/out/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=0
TMP_OUT_DIR="$HOME/xray/tmp"
# RUN_ID=$((SGE_TASK_ID-1))
# TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

CIF_FILE="/wynton/home/sali/mhancock/xray/data/reflections/7mhf/7mhf_refine.cif"
RES=0
W_XRAY=1.0
DYN_W_XRAY=1
START_PDB_FILE="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/43_synth/9337679/output_0/pdbs/189.pdb"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine.pdb"
DOF="S"
T=300
SA="300,0,10,A;300,-1,90,S"
LOG_FILE="$OUT_DIR/log.csv"

python /wynton/home/sali/mhancock/xray/dev/12_sulfur_only/scripts/sample.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_file "$CIF_FILE" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --start_pdb_file "$START_PDB_FILE" --ref_pdb_file "$REF_PDB_FILE" --dof "$DOF" --T "$T" --sa "$SA" --log_file "$LOG_FILE"
# python /wynton/home/sali/mhancock/xray/dev/12_sulfur_only/scripts/sample.py --out_dir "$OUT_DIR" --start_pdb_file "$START_PDB_FILE" --ref_pdb_file "$REF_PDB_FILE" --dof "$DOF" --T "$T" --sa "$SA" --log_file "$LOG_FILE"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT