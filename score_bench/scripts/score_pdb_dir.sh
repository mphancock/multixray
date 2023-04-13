#$ -S /bin/bash
#$ -o /wynton/home/sali/mhancock/xray/score_bench/scripts
#$ -N score
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -pe smp 25
#$ -l hostname='qb3-id*'

eval "$(conda shell.bash hook)"
module load CBI conda-stage
conda activate imp_218_cctbx

JOB_NAME=00_ff
DECOY_NAME=3ca7_N_1000_x5

PDB_DIR="/wynton/group/sali/mhancock/xray/decoys/data/$JOB_NAME/$DECOY_NAME"
CIF_FILE="$HOME/xray/data/reflections/3ca7/3ca7.cif"
MIN_RES=4
W_XRAY=30000
INCLUDE_NATIVE=0

SCORE_DIR="$HOME/xray/score_bench/data/$JOB_NAME"
mkdir -p "$SCORE_DIR"
SCORE_FILE="$SCORE_DIR/$DECOY_NAME"_3ca7_$MIN_RES.csv

python ~/xray/score_bench/scripts/score_bench_3ca7.py --pdb_dir "$PDB_DIR" --cif_file "$CIF_FILE" --min_res "$MIN_RES" --w_xray "$W_XRAY" --score_file "$SCORE_FILE" --include_native "$INCLUDE_NATIVE"