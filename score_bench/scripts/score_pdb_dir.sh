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

JOB_NAME=38_7mhf_decoys
DECOY_NAME=merge_2000
MIN_RES=0
CIF_FILE="$HOME/xray/data/reflections/7mhf/7mhf_no_H20_alt_H_ions.cif"
REF_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_no_H20_alt_H_ion.pdb"

SCORE_FILE_NAME="$DECOY_NAME"_7mhf_no_H20_alt_H_ions_"$MIN_RES"
PDB_DIR="/wynton/group/sali/mhancock/xray/decoys/data/7mhf/$JOB_NAME/$DECOY_NAME"

SCORE_DIR="$HOME/xray/score_bench/data/7mhk/$JOB_NAME"
mkdir -p "$SCORE_DIR"
SCORE_FILE="$SCORE_DIR/$SCORE_FILE_NAME.csv"
PARAMS_FILE="$SCORE_DIR/$SCORE_FILE_NAME.txt"
rm "$SCORE_FILE"
rm "$PARAMS_FILE"

# Can set the --add_native
python ~/xray/score_bench/scripts/score_pdb_dir.py --pdb_dir "$PDB_DIR" --cif_file "$CIF_FILE" --min_res "$MIN_RES" --ref_pdb_file "$REF_PDB_FILE" --param_file "$PARAMS_FILE" --score_file "$SCORE_FILE" --add_native
