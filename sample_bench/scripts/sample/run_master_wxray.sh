#! /bin/bash


JOB_NAME="127_2_cond_wxray"
EXP_ID=127
H_RT="03:00:00"
N_JOBS="1-500"
OFFSET="0"

W_XRAYS=(0.1 0.25 0.5 1 1.5 2 5 10)


for JOB_ID in {0..7}
do
    echo "$JOB_ID"
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"

    W_XRAY=${W_XRAYS[$JOB_ID]}

    # PARAMS="--cif_files /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/0/0.cif --w_xray $W_XRAY --dyn_w_xray --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 2 --init_weights rand --ref_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30/0.pdb --ref_occs 0.5105879489376552,0.4894120510623448 --sa {step99999,T300,dofA,pdb1,w0,res2} --bfactor 15"

    # PARAMS="--cif_files /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/0/0.cif,/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/1/0.cif --w_xray $W_XRAY --dyn_w_xray --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 2 --init_weights 0.5105879489376552,0.4894120510623448;0.5555655588277791,0.4444344411722208 --ref_pdb_file /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/7mhf_20/0.pdb --ref_occs 0.5105879489376552,0.4894120510623448;0.5555655588277791,0.4444344411722208 --sa {step99999,T300,dofA,pdb1,w0,res2} --bfactor 15"

    # PARAMS="--cif_files /wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhi.cif --w_xray $W_XRAY --dyn_w_xray --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 1 --init_weights rand --n_cond 1 --ref_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --ref_occs 1 --sa {step99999,T300,dofA,pdb1,w1,res0} --bfactor 15 --steps 1000"

    PARAMS="--cif_files /wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhi.cif,/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhk.cif --w_xray $W_XRAY --dyn_w_xray --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 2 --init_weights rand --n_cond 2 --ref_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --ref_occs 1;1 --sa {step99999,T300,dofA,pdb1,w1,res0} --bfactor 15 --steps 1000"

    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
done
