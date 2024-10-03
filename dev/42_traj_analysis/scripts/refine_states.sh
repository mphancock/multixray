#$ -S /bin/bash
#$ -N refine
#$ -cwd
#$ -o /wynton/group/sali/mhancock/xray/sample_bench/tmp/$JOB_ID.$TASK_ID.o
#$ -j y
#$ -l h_rt=00:30:00
#$ -l mem_free=1G
#$ -l scratch=1G
#$ -l hostname='qb3-id*'
#$ -t 1-10

source /wynton/home/sali/mhancock/phenix-1.21.2-5419/phenix_env.sh

RUN_ID=$((SGE_TASK_ID-1))

for STATE in {0..1}
do
    cd $TMPDIR

    FILE_NAME="$RUN_ID"_"$STATE"

    PDB_FILE="/wynton/home/sali/mhancock/xray/dev/42_traj_analysis/data/246/pdbs_state/32/$FILE_NAME.pdb"
    CIF_FILE="/wynton/home/sali/mhancock/xray/data/cifs/3k0m/3k0n.cif"

    phenix.refine "$CIF_FILE" "$PDB_FILE" base_output_dir="$TMPDIR" strategy=none ordered_solvent=true ordered_solvent.mode=every_macro_cycle refinement.input.xray_data.labels="_refln.intensity_meas,_refln.intensity_sigma" refinement.input.xray_data.r_free_flags.label="_refln.status" write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false write_mtz_file=false

    ls "$TMPDIR"

    cp "$TMPDIR"/"$FILE_NAME"_refine_001.pdb /wynton/home/sali/mhancock/xray/dev/42_traj_analysis/data/246/pdbs_state/32/"$FILE_NAME"_ref.pdb

    rm $TMPDIR/*
done

[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
