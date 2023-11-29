#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=6:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N md_116
#$ -t 1-1000
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_219_cctbx


JOB_NAME="116"
TARGET="3ca7"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$TARGET/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"


CIF_FILES="/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/2_state_0/0.cif,/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/2_state_1/0.cif"
W_XRAY=1.0
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
N_STATE=2
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_0/0.pdb,/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/2_state_1/0.pdb"
SA="{step1000,T300,dofA,pdb1,w1,res0};{step1000,T1000,dofS,pdb0,w0,res-1};{step1000,T1000,dofS,pdb0,w0,res3}"
LOG_FILE="$OUT_DIR/log.csv"
INIT_WEIGHTS="ref"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files $CIF_FILES --w_xray $W_XRAY --dyn_w_xray --start_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights "$INIT_WEIGHTS" --ref_pdb_file $REF_PDB_FILE --sa "$SA" --bfactor 15


[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT