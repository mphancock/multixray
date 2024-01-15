#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=00:10:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N md_54
#$ -t 1-5000
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_219_cctbx


JOB_NAME="54_1000"
TARGET="3ca7"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$TARGET/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --start_pdb_file "/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb" --n_state 2 --init_weights "rand" --ref_pdb_file "/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb" --sa "{step1000,T1000,dofA,pdb1,w0,res0}" --bfactor 15 --steps 1000


[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT