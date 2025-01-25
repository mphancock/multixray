#!/bin/bash

cd ~/xray/tmp
rm *
cp ~/xray/sample_bench/scripts/refine/refine_all_models_phenix.py .

python refine_all_models_phenix.py --out_dir=/wynton/group/sali/mhancock/xray/sample_bench/out/280_exp_all_2/9/output_0     --job_csv_file=/wynton/home/sali/mhancock/xray/sample_bench/data/params/280.csv --tmp_dir=/wynton/home/sali/mhancock/xray/tmp --log_phenix