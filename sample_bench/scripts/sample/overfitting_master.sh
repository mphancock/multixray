#! /bin/bash


JOB_NAME="63_overfit"
EXP_ID=63
H_RT="24:00:00"
N_JOBS="1-1"

CIF_FILES="$CIF_DIR/$JOB_ID.cif"
START_PDB_FILE="/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine.pdb"
N_STATE=2
INIT_WEIGHTS="rand"
T=300
SA="300,0,2000,A;1000,-1,100,A"

for JOB_ID in {0..7}
do
    JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/$JOB_NAME/$JOB_ID"

    if [ $JOB_ID == 2 ]; then
        W_XRAY=0.25
    else
        W_XRAY=0.5
    fi

    if [ $JOB_ID == 5 ]; then
        CIF_FILES="/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhj.cif,/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhf.cif"
    else
        CIF_FILES="/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhj.cif"
    fi

    PARAMS="--cif_files $CIF_FILES --dyn_w_xray --w_xray $W_XRAY --start_pdb_file $START_PDB_FILE --ref_pdb_file $START_PDB_FILE --n_state $N_STATE --init_weights $INIT_WEIGHTS --T $T --weights --sa $SA --bfactor $B_FACTOR"

    if [ $JOB_ID == 1 ]; then
    # If the condition is met, append the strings together
        PARAMS="$PARAMS --dropout "
    elif [ $JOB_ID == 2 ]; then
        PARAMS="$PARAMS --d_min 2.5 "
    elif [ $JOB_ID == 4 ]; then
        PARAMS="$PARAMS --rand_noise "
    elif [ $JOB_ID == 6 ]; then
        PARAMS="$PARAMS --main_chain "
    fi

    echo "$JOB_ID"
    echo "$PARAMS"
    qsub -N slave"$EXP_ID"_"$JOB_ID" -l h_rt=$H_RT -t "$N_JOBS" "$HOME/xray/sample_bench/scripts/sample/run_slave.sh" "$JOB_NAME" "$JOB_DIR" "$PARAMS"
done

# W_XRAY=$( echo "$JOB_ID*.1" | bc )