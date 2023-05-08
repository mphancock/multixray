#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/xray/score_bench/tmp
#$ -N score
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -pe smp 25
#$ -l hostname='qb3-id*'

# eval "$(conda shell.bash hook)"
# module load CBI conda-stage
# conda activate imp_218_cctbx

JOB_NAME=51_16
DECOY_NAME=best_1000
MIN_RES=0
W_XRAY=30000
CIF_FILE="$HOME/xray/data/reflections/3ca7/3ca7.cif"

SCORE_FILE_NAME="$DECOY_NAME"_3ca7_"$MIN_RES"
PDB_DIR="/wynton/group/sali/mhancock/xray/decoys/data/3ca7/$JOB_NAME/$DECOY_NAME"

SCORE_DIR="$HOME/xray/score_bench/data/3ca7/$JOB_NAME"
mkdir -p "$SCORE_DIR"
SCORE_FILE="$SCORE_DIR/$SCORE_FILE_NAME.csv"
PARAMS_FILE="$SCORE_DIR/$SCORE_FILE_NAME.txt"
rm "$SCORE_FILE"
rm "$PARAMS_FILE"

python ~/xray/score_bench/scripts/score_bench_3ca7.py --pdb_dir "$PDB_DIR" --cif_file "$CIF_FILE" --min_res "$MIN_RES" --w_xray "$W_XRAY" --score_file "$SCORE_FILE" --param_file "$PARAMS_FILE" --add_native
