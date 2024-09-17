#!/bin/bash

# Define the source and destination directories
src_base="/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/217_no_ref_com/5/output_"
dest_base="/wynton/home/sali/mhancock/xray/dev/42_traj_analysis/data/217/pdbs/5/"

# Loop from 1 to 25 for the NUM value
for num in {0..25}; do
  # Define the source and destination file paths
  src_file="${src_base}${num}/pdbs/500.pdb"
  dest_file="${dest_base}${num}.pdb"

  # Check if the source file exists before copying
  if [ -f "$src_file" ]; then
    # Copy the pdb file to the destination
    cp "$src_file" "$dest_file"
    echo "Copied $src_file to $dest_file"
  else
    echo "Source file $src_file not found, skipping."
  fi
done