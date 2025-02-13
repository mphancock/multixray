#!/bin/bash

cd ~/xray/tmp
rm *
cp ~/xray/sample_bench/scripts/refine/refine_all_models_phenix.py .

python refine_all_models_phenix.py --out_dir=/wynton/group/sali/mhancock/xray/sample_bench/out/283_2_cond_ref/0/output_0 --job_csv_file=~/xray/sample_bench/data/params/283.csv --tmp_dir=/wynton/home/sali/mhancock/xray/tmp --log_phenix