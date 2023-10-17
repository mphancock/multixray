#$ -S /bin/bash
#$ -o /wynton/group/sali/mhancock/xray/score_bench/tmp
#$ -N score
#$ -j y
#$ -l h_rt=2:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -pe smp 25
#$ -l hostname='qb3-id*'

# eval "$(conda shell.bash hook)"
# module load CBI conda-stage
# conda activate imp_218_cctbx

DECOY_NAME=0
JOB_NAME=100_natives_4x
MIN_RES=0
CIF_FILE="$HOME/xray/data/cifs/3ca7/3ca7_refine.cif"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
SCORE_FS="ml,rmsd_avg"
PDB_DIR="/wynton/group/sali/mhancock/xray/decoys/data/3ca7/$JOB_NAME/$DECOY_NAME"

REF_PDB_FILE_NAME="$(basename $REF_PDB_FILE .pdb)"
SCORE_FILE_NAME="$REF_PDB_FILE_NAME"_"$DECOY_NAME"
echo $SCORE_FILE_NAME

SCORE_DIR="$HOME/xray/score_bench/data/3ca7/$JOB_NAME"
mkdir -p "$SCORE_DIR"
SCORE_FILE="$SCORE_DIR/$SCORE_FILE_NAME.csv"
PARAMS_FILE="$SCORE_DIR/$SCORE_FILE_NAME.txt"
rm "$SCORE_FILE"
rm "$PARAMS_FILE"

# Can set the --add_native
python ~/xray/score_bench/scripts/score_pdb_dir.py --pdb_dir "$PDB_DIR" --cif_file "$CIF_FILE" --min_res "$MIN_RES" --ref_pdb_file "$REF_PDB_FILE" --param_file "$PARAMS_FILE" --score_file "$SCORE_FILE" --score_fs "$SCORE_FS" --add_native
