#!/bin/bash

MTZ_FILE="/wynton/home/sali/mhancock/xray/data/maps/7mhf/7mhf_no_H2O_alt_H_ions.mtz"
CIF_FILE="/wynton/home/sali/mhancock/xray/data/reflections/7mhf/7mhf_no_H2O_alt_H_ions.cif"

phenix.mtz_as_cif "$MTZ_FILE" output_file="$CIF_FILE" mtz_labels="F FreeR_flag" cif_labels="_refln.F_meas_au _refln.status"