#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=06:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N md_183
#$ -t 1-250
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_220_cctbx


JOB_NAME="183_test"
TARGET="7mhf"
JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/$TARGET/$JOB_NAME/$JOB_ID"
mkdir -p "$JOB_DIR"

RUN_ID=$((SGE_TASK_ID-1))
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"

cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --input_csv /wynton/home/sali/mhancock/xray/tmp/test.csv --job_id 0 --w_xray 0.5 --n_state 2 --init_weights rand --sa "{step3000,T300,dofA,pdb1,w1,res0}" --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine.pdb --steps 2

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT