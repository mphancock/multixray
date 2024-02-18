#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=03:00:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -N md_131
#$ -t 1-2000
#$ -l hostname='qb3-id*'

# Setup the conda environment.
eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_219_cctbx


JOB_NAME="131_1_state"
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
python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files /wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhi.cif --dyn_w_xray --w_xray 0.5 --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 1 --init_weights rand --n_cond 1 --ref_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --ref_occs "1" --sa "{step99999,T300,dofA,pdb1,w1,res0}" --bfactor 15


[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT