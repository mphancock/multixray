#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.o
#$ -j y
#$ -l h_rt=2:00:00
#$ -l mem_free=5G
#$ -l scratch=5G
#$ -N fix_46
#$ -pe smp 25
#$ -l hostname='qb3-id*'


module load CBI conda-stage
conda activate imp_218_cctbx

JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/46_synth_sa_3/9417275"
COPY_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/46_synth_sa_3/9417275_xray"

python ~/xray/dev/13_fix_sa_output/scripts/run.py --job_dir "$JOB_DIR" --copy_dir "$COPY_DIR"

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"