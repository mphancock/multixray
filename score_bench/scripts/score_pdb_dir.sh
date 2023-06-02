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

JOB_NAME=14_32_h20
DECOY_NAME=best_1000
MIN_RES=0
CIF_FILE="$HOME/xray/data/reflections/7mhk/7mhk.cif"
UC_DIM="114.300 54.290 44.970 90.00 102.12 90.00"
SG_SYMBOL="C 1 2 1"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhk/7mhk_clean_h20.pdb"

SCORE_FILE_NAME="$DECOY_NAME"_7mhk_"$MIN_RES"
PDB_DIR="/wynton/group/sali/mhancock/xray/decoys/data/7mhk/$JOB_NAME/$DECOY_NAME"

SCORE_DIR="$HOME/xray/score_bench/data/7mhk/$JOB_NAME"
mkdir -p "$SCORE_DIR"
SCORE_FILE="$SCORE_DIR/$SCORE_FILE_NAME.csv"
PARAMS_FILE="$SCORE_DIR/$SCORE_FILE_NAME.txt"
rm "$SCORE_FILE"
rm "$PARAMS_FILE"

# Can set the --add_native
python ~/xray/score_bench/scripts/score_pdb_dir.py --pdb_dir "$PDB_DIR" --cif_file "$CIF_FILE" --min_res "$MIN_RES" --uc_dim "$UC_DIM" --sg_symbol "$SG_SYMBOL" --ref_pdb_file "$REF_PDB_FILE" --param_file "$PARAMS_FILE" --score_file "$SCORE_FILE"
