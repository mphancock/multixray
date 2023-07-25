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

DECOY_NAME=rand_1000_2x
SCORE_FILE_NAME="$DECOY_NAME"_to_1x

JOB_NAME=53_100
MIN_RES=0
CIF_FILE="$HOME/xray/data/reflections/3ca7/3ca7_refine_2.cif"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb"
SCORE_FS="ml,rmsd_order,rmsd_avg,weight_delta,rmsd_dom,dom_weight"
PDB_DIR="/wynton/group/sali/mhancock/xray/decoys/data/3ca7/$JOB_NAME/$DECOY_NAME"

SCORE_DIR="$HOME/xray/score_bench/data/3ca7/$JOB_NAME"
mkdir -p "$SCORE_DIR"
SCORE_FILE="$SCORE_DIR/$SCORE_FILE_NAME.csv"
PARAMS_FILE="$SCORE_DIR/$SCORE_FILE_NAME.txt"
rm "$SCORE_FILE"
rm "$PARAMS_FILE"

# Can set the --add_native
python ~/xray/score_bench/scripts/score_pdb_dir.py --pdb_dir "$PDB_DIR" --cif_file "$CIF_FILE" --min_res "$MIN_RES" --ref_pdb_file "$REF_PDB_FILE" --param_file "$PARAMS_FILE" --score_file "$SCORE_FILE" --score_fs "$SCORE_FS" --add_native
