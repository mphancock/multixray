#!/bin/bash

TMPDIR="/wynton/home/sali/mhancock/xray/sample_bench/out/test"
rm -r "$TMPDIR"
mkdir -p "$TMPDIR"

JOB_DIR="/wynton/group/sali/mhancock/xray/sample_bench/out/test"
rm -r "$JOB_DIR"
mkdir -p "$JOB_DIR"

RUN_ID=0
TMP_OUT_DIR="$TMPDIR/output_$RUN_ID"
OUT_DIR="$JOB_DIR/output_$RUN_ID"
mkdir "$TMP_OUT_DIR"
mkdir "$OUT_DIR"


cd "$TMPDIR"
cp ~/xray/sample_bench/scripts/sample/run_md_multi.py .

python run_md_multi.py --out_dir "$OUT_DIR" --tmp_out_dir "$TMP_OUT_DIR" --cif_files /wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhi.cif,/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhk.cif,/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhf.cif --dyn_w_xray --w_xray 0.5 --start_pdb_file /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi_ref.pdb --n_state 2 --init_weights rand --n_cond 3 --ref_pdb_files /wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhi.pdb,/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhk.pdb,/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf.pdb --ref_occs "1;1;1" --sa "{step99999,T300,dofA,pdb1,w1,res2}" --bfactor 15

cd "$HOME/xray/sample_bench/scripts/sample"