#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.o
#$ -j y
#$ -l h_rt=1:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N vol_bench
#$ -pe smp 50
#$ -l hostname='qb3-id*'


module load CBI conda-stage
conda activate imp_218_cctbx

JOB_NAME="14_32_h20"
JOB_ID="3608088"
MAX_N="1000"

JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhk/$JOB_NAME/$JOB_ID"
OUT_FILE="/wynton/home/sali/mhancock/xray/sample_bench/data/7mhk/$JOB_NAME/sample_bench.csv"

python run_sample_volume_bench.py --job_dir "$JOB_DIR" --out_file "$OUT_FILE" --max_n "$MAX_N"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
trap 'conda deactivate' EXIT