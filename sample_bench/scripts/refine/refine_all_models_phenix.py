from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import argparse
import time
import shutil
import os

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import read_pdb_and_refine, read_pdb_and_refine_to_max_ff, refine_posterior
from params import read_job_csv, build_weights_matrix
from align_imp import align_one_to_two
sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
from score import pool_score
from files import pdb_to_df, write_pdb_from_df, update_occs, update_model_based_on_altconf, update_alt_loc_by_model, insert_single_atom, renumber_hetero_residues, duplicate_heteroatoms_for_all_altlocs
from miller_ops import get_crystal_symmetry
from etc import get_zn_coords_and_occ_after_align


if __name__ == "__main__":
    ## 1. SETUP
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_dir")
    parser.add_argument("--job_csv_file")
    parser.add_argument("--log_phenix", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    tmp_dir = Path(args.tmp_dir)
    pdb_dir = Path(out_dir, "pdbs")

    if not pdb_dir.exists():
        raise RuntimeError("{} does not exist".format(pdb_dir))

    job_dir = out_dir.parents[0]
    job_name = job_dir.name

    input_csv = Path(args.job_csv_file)
    input_id = int(job_dir.name)

    params_dict = read_job_csv(input_csv, input_id)
    N = params_dict["N"]
    J = params_dict["J"]
    cif_files = params_dict["cifs"]
    ref_files = params_dict["ref_pdbs"]
    ref_w_mat = params_dict["ref_w_mat"]
    cif_names = [cif_file.stem for cif_file in cif_files]

    exp_dir = job_dir.parents[0]

    if not pdb_dir.exists():
        raise RuntimeError("PDB directory does not exist")

    log_file = Path(args.out_dir, "log.csv")
    log_df = pd.read_csv(log_file, index_col=0)
    log_df = log_df.loc[log_df["pdb"].notnull()]

    refined_exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/{}_phenix_ref".format(exp_dir.name))
    refined_out_dir = Path(refined_exp_dir, job_dir.name, out_dir.name)
    shutil.rmtree(refined_out_dir, ignore_errors=True)
    refined_pdb_dir = Path(refined_out_dir, "pdbs")
    refined_pdb_dir.mkdir(parents=True, exist_ok=True)

    refined_log_file = Path(refined_out_dir, "log.csv")
    refined_log_df = log_df.copy()

    refined_log_df["ff"] = np.nan
    # refined_log_df["pdb"] = np.nan
    refined_log_df["time"] = np.nan
    for cond in range(J):
        cif_name = cif_names[cond]
        for col_type in ["r_free", "r_work", "xray"]:
            refined_log_df["{}_{}".format(col_type, cif_name)] = np.nan

    # ## drop any rows that have pdb as nan
    refined_log_df.reset_index(drop=True, inplace=True)

    # ## pick 10 random rows to refine
    n_ref = 10
    refined_log_df = refined_log_df.sample(n=n_ref)
    # refined_log_df = refined_log_df.loc[[len(refined_log_df)-1]]
    refined_log_df.reset_index(drop=True, inplace=True)
    pdb_files = list(refined_log_df["pdb"])

    ## now add duplicate rows to represent the different conditions
    ref_indices = list()
    for i in range(len(refined_log_df)):
        ref_indices.extend([i]*J)
    refined_log_df = refined_log_df.loc[ref_indices]
    refined_log_df.reset_index(drop=True, inplace=True)

    print(len(refined_log_df))
    print(refined_log_df.head())

    ## 2. REFINE ALL MODELS AGAINST ALL CIFS
    ## for each pdb file we will refine against all cif files
    for i in range(len(pdb_files)):
        ## convert the pdb_file to altconf but store in local scratch
        pdb_file = Path(pdb_files[i])
    # for pdb_file in [Path("/wynton/group/sali/mhancock/xray/sample_bench/out/282_test_w/0/output_0/pdbs/3.pdb")]:
        tmp_pdb_file = Path(tmp_dir, pdb_file.name)
        decoy_w_mat = build_weights_matrix(refined_log_df, i, "w", N, cif_names)

        ## refine against all cif files and save as PDB_CIF.pdb
        for j in range(len(cif_files)):
            cif_file = cif_files[j]
            cif_name = cif_names[j]

            ## get the pdb csv
            pdb_df = pdb_to_df(pdb_file)

            ## update confs and occs
            pdb_df = update_occs(pdb_df, decoy_w_mat[:,j])
            pdb_df = update_alt_loc_by_model(pdb_df)

            ## get zn coordinates and occs and add to pdb csv
            zn_coords, zn_occ = get_zn_coords_and_occ_after_align(
                pdb_file_1=Path(Path.home(), "xray/data/pdbs/7mhf/{}.pdb".format(cif_name)),
                pdb_file_2=pdb_file
            )
            pdb_df = insert_single_atom(
                df=pdb_df,
                atom_data={"model": 1, "record": "HETATM","atom_serial": 3000, "atom_name": "ZN", "alt_loc": "", "residue_name": "ZN", "chain_id": "A","residue_seq": 401, "insertion": "", "x": zn_coords[0], "y": zn_coords[1], "z": zn_coords[2], "occupancy": zn_occ, "temp_factor": 53.43, "element": "ZN", "charge": ""},
                index=None
            )

            ## save as altconf
            write_pdb_from_df(
                df=pdb_df,
                out_pdb_file=tmp_pdb_file,
                single_model=True
            )

            ## convert the orig multistate pdb file to altconf tmp file (jth column of the decoy_w_mat)
            cif_file = cif_files[j]
            crystal_symmetry = get_crystal_symmetry(cif_file)
            refined_pdb_file = Path(refined_pdb_dir, "{}_{}.pdb".format(pdb_file.stem, cif_names[j]))
            print(pdb_file, tmp_pdb_file, refined_pdb_file)
            sg = crystal_symmetry.space_group_info().group().info()
            uc = crystal_symmetry.unit_cell()

            refine_command = "phenix.refine {} {} strategy=individual_sites+individual_adp ordered_solvent=true ordered_solvent.mode=every_macro_cycle  refinement.input.xray_data.labels=_refln.F_meas_au,_refln.F_meas_sigma_au refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None crystal_symmetry.unit_cell={},{},{},{},{},{} crystal_symmetry.space_group='{}' write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false".format(cif_file, tmp_pdb_file, uc.parameters()[0], uc.parameters()[1], uc.parameters()[2], uc.parameters()[3], uc.parameters()[4], uc.parameters()[5], sg)

            if not args.log_phenix:
                refine_command += "> /dev/null 2>&1"

            print(refine_command)
            os.system(refine_command)
            phenix_out_pdb_file = Path(tmp_dir, "{}_refine_001.pdb".format(pdb_file.stem))

            # ## convert all models back to multistate
            df = pdb_to_df(phenix_out_pdb_file)
            df = duplicate_heteroatoms_for_all_altlocs(df)
            df = update_model_based_on_altconf(df)
            df = renumber_hetero_residues(df)
            write_pdb_from_df(
                df=df,
                out_pdb_file=refined_pdb_file,
                single_model=False
            )
            refined_log_df.loc[i*len(cif_files)+j, "pdb"] = refined_pdb_file

            ## clean up temporary files from phenix
            os.system("rm *.pdb")
            os.system("rm *.mtz")
            os.system("rm *.log")

    ## 3. SCORE ALL REFINED MODELS AGAINST CORRESPONDING CIF
    print(refined_log_df.head())
    for i in range(len(refined_log_df)):
        ## a bit confusing because we are running each refined pdb file against a single cif file so we need to find correct cif, correct weights, etc
        refined_pdb_file = Path(refined_log_df.loc[i, "pdb"])

        cif_id = i % J
        cif_file = cif_files[cif_id]
        cif_name = cif_names[cif_id]

        decoy_w_mat = build_weights_matrix(refined_log_df, i, "w", N, cif_names)

        param_dict = dict()
        param_dict["decoy_files"] = [refined_pdb_file]
        param_dict["decoy_w_mat"] = decoy_w_mat[:, cif_id].reshape([-1,1])

        ## only 1 ref file for synthetic benchmark
        param_dict["ref_file"] = ref_files[0]
        param_dict["ref_w_mat"] = ref_w_mat
        param_dict["score_fs"] = ["xray_0", "ff"]
        param_dict["cif_files"] = [cif_file]
        param_dict["scale_k1"] = True
        param_dict["scale"] = True
        param_dict["remove_outliers"] = True
        param_dict["res"] = 0

        score_dict = pool_score(param_dict)
        refined_log_df.loc[i, "xray_{}".format(cif_name)] = score_dict["xray_0"]
        refined_log_df.loc[i, "r_free_{}".format(cif_name)] = score_dict["r_free_0"]
        refined_log_df.loc[i, "r_work_{}".format(cif_name)] = score_dict["r_work_0"]
        # refined_log_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_{}".format(cif_name)]
        refined_log_df.loc[i, "ff"] = score_dict["ff"]

        print(score_dict)

    print(refined_log_df.head())
    refined_log_df.to_csv(refined_log_file)

