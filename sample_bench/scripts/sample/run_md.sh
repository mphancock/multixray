#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=12:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N md_61
#$ -t 1-2000
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx


JOB_NAME="61_sb_4"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

# CIF_FILES="/wynton/home/sali/mhancock/xray/data/reflections/3ca7/3ca7_refine_2.cif"
# RES=0
# W_XRAY=1.0
# DYN_W_XRAY=1
# COM=os
# START_PDB_FILE="/wynton/group/sali/mhancock/xray/decoys/data/3ca7/53_100/rand_1000_2x_53_54/832.pdb"
# N_STATE=2
# REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine_2.pdb"
# T=300
# SA="300,0,2000,A;1000,-1,100,A"

# cd "$TMPDIR"
# cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
# python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files "$CIF_FILES" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --com "$COM" --start_pdb_file "$START_PDB_FILE" --n_state "$N_STATE" --ref_pdb_file "$REF_PDB_FILE" --T "$T" --sa "$SA" --weights

CIF_FILES="/wynton/home/sali/mhancock/xray/data/reflections/3ca7/3ca7_refine_4.cif"
RES=0
W_XRAY=1.0
DYN_W_XRAY=1
COM=os
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
N_STATE=4
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine_4.pdb"
T=300

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files "$CIF_FILES" --res "$RES" --w_xray "$W_XRAY" --dyn_w_xray "$DYN_W_XRAY" --com "$COM" --start_pdb_file "$START_PDB_FILE" --n_state "$N_STATE" --ref_pdb_file "$REF_PDB_FILE" --T "$T"


[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT