#! /bin/bash


JOB_NAME="155_native_N4_decoys"
EXP_ID=155
H_RT="01:00:00"
N_JOBS="1-10000"
OFFSET="0"


for JOB_ID in {0..9}
do
    echo "$JOB_ID"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"

    # PARAMS="--cif_files /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/0/$JOB_ID.cif --w_xray .5 --dyn_w_xray --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 4 --init_weights rand --n_cond 1 --ref_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv --ref_id $JOB_ID --sa {step99999,T300,dofA,pdb1,w1,res2} --bfactor 15"

    # PARAMS="--cif_files /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/0/$JOB_ID.cif,/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/1/$JOB_ID.cif --w_xray .5 --dyn_w_xray --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 8 --init_weights rand --n_cond 2 --ref_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv --ref_id $JOB_ID --sa {step99999,T300,dofA,pdb1,w1,res2} --bfactor 15"

    # PARAMS="--start_pdb_file /wynton/home/sali/mhancock/xray/dev/33_best_scoring_1_state/data/pdbs/7mhf_30/$JOB_ID.pdb --n_state 1 --init_weights rand --n_cond 1 --ref_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv --ref_id $JOB_ID --sa {step250,T50,dofA,pdb1,w0,res0};{step250,T100,dofA,pdb1,w0,res0};{step1000,T1000,dofA,pdb1,w0,res0} --bfactor 15 --steps 1000"

    # PARAMS="--start_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30/$JOB_ID.pdb --n_state 2 --n_cond 2 --init_weights rand --ref_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv --ref_id $JOB_ID --sa {step250,T50,dofA,pdb1,w0,res0};{step250,T100,dofA,pdb1,w0,res0};{step1000,T1000,dofA,pdb1,w0,res0} --bfactor 15 --steps 1000"

    # NATIVES
    PARAMS="--start_pdb_file /wynton/home/sali/mhancock/xray/dev/33_best_scoring/data/4_state/$JOB_ID.pdb --n_state 4 --init_weights rand --n_cond 2 --ref_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/scores/7mhf_30.csv --ref_id $JOB_ID --sa {step250,T50,dofA,pdb1,w0,res0};{step250,T100,dofA,pdb1,w0,res0};{step250,T1000,dofA,pdb1,w0,res0} --bfactor 15 --steps 10"

    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
done
