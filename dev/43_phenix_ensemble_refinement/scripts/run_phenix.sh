#$ -S /bin/bash
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.o
#$ -j y
#$ -l h_rt=4:00:00
#$ -l mem_free=4G
#$ -l scratch=4G
#$ -N ensmble
#$ -t 1-1
#$ -l hostname='qb3-id*'


module load Sali phenix
cd /wynton/home/sali/mhancock/xray/dev/43_phenix_ensemble_refinement/scripts

MODEL_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhk.pdb"
# MODEL_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_clean.pdb"
DATA_FILE="/wynton/home/sali/mhancock/xray/data/maps/7mhk.mtz"

phenix.ensemble_refinement $MODEL_FILE $DATA_FILE ensemble_refinement.output_file_prefix=7mhk zip_final_model=False

