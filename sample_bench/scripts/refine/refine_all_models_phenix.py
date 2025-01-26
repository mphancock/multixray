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
sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
from score import pool_score
from files import multi_to_altconf, altconf_to_multi


if __name__ == "__main__":
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

    print(cif_files)

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
    refined_log_df = refined_log_df.sample(n=10)
    # refined_log_df = refined_log_df.loc[[len(refined_log_df)-1]]
    refined_log_df.reset_index(drop=True, inplace=True)

    pdb_files = list(refined_log_df["pdb"])
    for i in range(len(pdb_files)):
        ## convert the pdb_file to altconf but store in local scratch
        pdb_file = Path(pdb_files[i])
    # for pdb_file in [Path(out_dir, "pdbs", "500.pdb")]:
        tmp_pdb_file = Path(tmp_dir, pdb_file.name)
        refined_pdb_file = Path(refined_pdb_dir, "{}.pdb".format(pdb_file.stem))
        print(pdb_file, tmp_pdb_file, refined_pdb_file)

        decoy_w_mat = build_weights_matrix(refined_log_df, i, "w", N, cif_names)
        multi_to_altconf(
            in_pdb_file=pdb_file,
            occs=decoy_w_mat[:,0],
            out_pdb_file=tmp_pdb_file
        )

        refine_command = "phenix.refine {} {} strategy=individual_sites+individual_adp ordered_solvent=true ordered_solvent.mode=every_macro_cycle  refinement.input.xray_data.labels=_refln.F_meas_au,_refln.F_meas_sigma_au refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None crystal_symmetry.unit_cell=115.023,54.358,44.970,90.00,101.50,90.00 crystal_symmetry.space_group='C 1 2 1' write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false".format(cif_files[0], tmp_pdb_file)

        if not args.log_phenix:
            refine_command += "> /dev/null 2>&1"

        print(refine_command)
        os.system(refine_command)

        phenix_out_pdb_file = Path(tmp_dir, "{}_refine_001.pdb".format(pdb_file.stem))

        ## convert all models back to multistate
        altconf_to_multi(
            in_pdb_file=phenix_out_pdb_file,
            out_pdb_file=refined_pdb_file,
            n_state=N
        )
        refined_log_df.loc[i, "pdb"] = refined_pdb_file

        ## clean up temporary files from phenix
        os.system("rm *.pdb")
        os.system("rm *.mtz")
        os.system("rm *.log")

    print(refined_log_df.head())
    for i in range(len(refined_log_df)):
        refined_pdb_file = Path(refined_log_df.loc[i, "pdb"])
        refined_log_df.loc[i, "pdb"] = refined_pdb_file

        cif_name = cif_names[0]
        decoy_w_mat = build_weights_matrix(refined_log_df, i, "w", N, cif_names)

        param_dict = dict()
        param_dict["decoy_files"] = [refined_pdb_file]
        param_dict["decoy_w_mat"] = decoy_w_mat

        ## only 1 ref file for synthetic benchmark
        param_dict["ref_file"] = ref_files[0]
        param_dict["ref_w_mat"] = ref_w_mat
        param_dict["score_fs"] = ["xray_0", "ff"]
        param_dict["cif_files"] = cif_files
        param_dict["scale_k1"] = True
        param_dict["scale"] = True
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

