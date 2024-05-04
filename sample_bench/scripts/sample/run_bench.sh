#! /bin/bash


EXP_ID=182
EXP_NAME="182_bench"
H_RT="06:00:00"
N_JOBS="1-250"
OFFSET="0"

N_STATES=(1 2)

for N_STATE_ID in {0..1}
do
    N_STATE=${N_STATES[$N_STATE_ID]}
    for JOB_ID in {0..19}
    do
        JOB_NAME=n"$N_STATE_ID"_j"$JOB_ID"
        echo "$JOB_NAME"
        JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$EXP_NAME/$JOB_NAME"

        PARAMS="--input_csv /wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/csvs/7mhf_30.csv --job_id $JOB_ID --w_xray /wynton/home/sali/mhancock/xray/sample_bench/data/7mhf/180_wxray_bench_ref_10000/best_wxray.csv --n_state $N_STATE --init_weights rand --sa {step3000,T300,dofA,pdb1,w1,res2} --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine.pdb --steps 2"

        qsub -N b"$JOB_NAME" -l h_rt=$H_RT -l mem_free=1G -l scratch=1G -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS" "$OFFSET"
    done
done
