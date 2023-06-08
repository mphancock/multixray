#!/bin/bash


FILES="$HOME/xray/dev/07_r_free/data/7mhf/*"
for f in $FILES
do
    echo "$f"
    phenix.model_vs_data "$f" ~/xray/data/reflections/7mhf/7mhf.cif f_obs_label="F_meas_au" r_free_flags_label="status"
    echo "\n\n\n"
done