from pathlib import Path
import sys
import argparse

import pandas as pd
import numpy as np

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from align_imp import align_ones_to_twos
from files import pdb_to_df, write_pdb_from_df
from params import read_job_csv, build_weights_matrix
from utility import get_n_state_from_pdb_file
from score import pool_score


### "refine" structure by copying atoms from reference structure
if __name__ == "__main__":
    ## turn off logging for IMP bc will print a lot of warnings for CHARMM
    IMP.set_log_level(IMP.SILENT)

    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir")
    parser.add_argument("--tmp_dir")
    parser.add_argument("--job_csv_file")
    args = parser.parse_args()

    log_file = Path(args.out_dir, "log.csv")
    score_df = pd.read_csv(log_file)
    for col in score_df.loc[:, score_df.columns.str.contains('xray|r_free|r_work', case=False)]:
        score_df[col] = np.nan

    pdb_file = Path(score_df.loc[0, "pdb"])
    out_dir = pdb_file.parents[1]
    job_dir = pdb_file.parents[2]
    input_id = int(job_dir.name)
    exp_dir = pdb_file.parents[3]
    exp_name = exp_dir.name
    new_exp_dir = Path(exp_dir.parents[0], exp_name + "_copy")
    new_pdb_dir = Path(new_exp_dir, job_dir.name, out_dir.name, "pdbs")
    new_log_file = Path(new_exp_dir, job_dir.name, out_dir.name, "log.csv")
    new_pdb_dir.mkdir(parents=True, exist_ok=True)

    input_csv = Path(args.job_csv_file)
    params_dict = read_job_csv(input_csv, input_id)
    n_state = params_dict["N"]

    for i in range(len(score_df)):
        pdb_file = Path(score_df.loc[i, "pdb"])
        out_pdb_file = Path(new_pdb_dir, pdb_file.name)
        print(out_pdb_file)

        cif_name = pdb_file.stem.split("_")[1]
        cif_file = Path("/wynton/home/sali/mhancock/xray/dev/38_standard_flags/data/{}.cif".format(cif_name))
        ref_pdb_file = Path("/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/{}.pdb".format(cif_name))

        ## 1 -- first align the pdb file to the ref pdb file
        ref_m, m = IMP.Model(), IMP.Model()
        ref_h = IMP.atom.read_pdb(str(ref_pdb_file), ref_m, IMP.atom.ATOMPDBSelector())
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.ATOMPDBSelector())

        align_ones_to_twos(hs, [ref_h, ref_h])
        IMP.atom.write_multimodel_pdb(hs, str(out_pdb_file))

        ## 2 -- then add the atoms from the ref pdb file
        ref_pdb_df = pdb_to_df(ref_pdb_file)
        pdb_df = pdb_to_df(out_pdb_file)

        for state in range(n_state):
            ref_het_atoms = ref_pdb_df[ref_pdb_df['record'] == "HETATM"].copy()
            ref_het_atoms['model'] = state+1
            pdb_df = pd.concat([pdb_df, ref_het_atoms], ignore_index=True)

        write_pdb_from_df(pdb_df, out_pdb_file)

        ## 3 -- score the pdb files
        score_df.loc[i, "pdb"] = str(out_pdb_file)
        decoy_w_mat = build_weights_matrix(score_df, i, "w", n_state, [cif_name])

        param_dict = dict()
        param_dict["decoy_files"] = [out_pdb_file]
        param_dict["decoy_w_mat"] = decoy_w_mat[:, 0].reshape([-1,1])

        ## only 1 ref file for synthetic benchmark
        param_dict["ref_file"] = ref_pdb_file

        ## not a good assumption at the moment
        param_dict["ref_w_mat"] = np.array([1]).reshape([-1,1])
        param_dict["score_fs"] = ["xray_0", "ff"]
        param_dict["cif_files"] = [cif_file]
        param_dict["scale_k1"] = True
        param_dict["scale"] = True
        param_dict["remove_outliers"] = True
        param_dict["res"] = 0

        score_dict = pool_score(param_dict)
        score_df.loc[i, "xray_{}".format(cif_name)] = score_dict["xray_0"]
        score_df.loc[i, "r_free_{}".format(cif_name)] = score_dict["r_free_0"]
        score_df.loc[i, "r_work_{}".format(cif_name)] = score_dict["r_work_0"]
        # refined_log_df.loc[i, "rmsd_{}".format(cif_name)] = score_dict["rmsd_{}".format(cif_name)]
        score_df.loc[i, "ff"] = score_dict["ff"]

        print(score_dict)

    score_df.to_csv(new_log_file)
    print(score_df.head())