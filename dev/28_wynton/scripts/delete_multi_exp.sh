#!/bin/bash


SYS_DIR=/wynton/group/sali/mhancock/xray/sample_bench/out
# SYS_DIR=/salilab/park1/matthew/xray

EXP_NAMES=(186_w_xray_bench 187_bench 188_w_xray_exp 189_exp 190_free_control 191_orig_free 192_w_xray_free_control 193_bench_free_control 194_free_control_all_N 195_w_bench 196_w_bench 197_loop 198_all 199_refine 200_gradual 201_lower 202_no_wxray_auto 203_high_T 204_inflection 205_9E3_high_T 206_no_scale 207_inflection_2 208_8900 209_8900_low_res 210_update_freq 211_xray_freq 212_8900_high_T 213_singles 214_thermostat 215_no_thermostat 216_no_thermostat_no_com 217_no_ref_com 218_no_scale 219_no_k1_scale 220_mask 221_2_cif 222_wxray 223_wxray_thermo 224_all_states 225_all_states_freq 226_exp 227_ortho 228_thermo 229_3k0m 230_3k0m_freq 231_thermo 232_weights 233_weights 234_weights_10k 235_temp 236_2_state 237_all_states 238_res 239_waters 240_aniso 241_aniso_wxray 242_multi_wxray 243_3k0n 244_no_thermo 245_vel_cap 246_temp 247_res 248_temp_res 249_temp_4 250_refine 251_temp_er_state 252_auto_wxray 253_7mhf 254_r_scale 255_7mhf 256_ff 257_xray 258_wxray 259_sb_temp 260_1_state 261_2_state_test)

for EXP_NAME in "${EXP_NAMES[@]}"
do
    EXP_DIR=$SYS_DIR/$EXP_NAME
    echo $EXP_DIR
    nohup rm -r $EXP_DIR &
done