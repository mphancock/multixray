from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
# from refine import refine, pool_refine
from utility import pool_read_pdb
from params import read_job_csv
from stat_df import pool_get_stat_info_df


def get_best_from_out_dir(
    params
):
    out_dir = params["out_dir"]
    field = params["field"]
    bonus_fields = params["bonus_fields"]
    log_files = [Path(out_dir, "log.csv")]

    params = dict()
    params["log_files"] = log_files
    params["equil"] = 0
    params["field"] = field
    params["bonus_fields"] = bonus_fields
    params["N"] = 1
    params["pdb_only"] = True
    params["max_rmsd"] = None
    params["max_ff"] = None
    stat_df = pool_get_stat_info_df(params)

    return stat_df


if __name__ == "__main__":
    exp_name = "267_full_ref"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out", exp_name)
    job_id = 1
    job_dir = Path(exp_dir, str(job_id))

    job_csv_file = Path(Path.home(), "xray/sample_bench/data/params/267.csv")
    sample_file = Path(Path.home(), "xray/sample_bench/data/analysis/{}/all_outs_{}.csv".format(exp_name, job_id))

    job_df = read_job_csv(job_csv_file=job_csv_file, job_id=job_id)

    cif_files = job_df["cifs"]
    cif_names = [Path(cif).stem for cif in cif_files]
    J = job_df["J"]
    N = job_df["N"]

    job_dir = Path(exp_dir, str(job_id))
    out_dirs = [out_dir for out_dir in job_dir.glob("output*")]

    # field = "xray_native_0+xray_native_1"
    field = "xray_native_0"
    bonus_fields = ["ff", "pdb", "rmsd"]

    for cif_name in cif_names:
        bonus_fields.append("rmsd_{}".format(cif_name))
        bonus_fields.append("xray_{}".format(cif_name))
        bonus_fields.append("r_work_{}".format(cif_name))
        bonus_fields.append("r_free_{}".format(cif_name))

    w_cols = list()
    for state in range(N):
        for cond in range(J):
            w_cols.append("w_{}_{}".format(state, cif_names[cond]))
    bonus_fields.extend(w_cols)

    pool_params = list()
    for out_dir in out_dirs:
        params_dict = dict()
        params_dict["out_dir"] = out_dir
        params_dict["field"] = field
        params_dict["bonus_fields"] = bonus_fields
        pool_params.append(params_dict)

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(get_best_from_out_dir, pool_params)

    sample_df = pd.DataFrame()
    for stat_df in pool_results:
        if len(stat_df) == 0:
            continue

        out_dir = Path(stat_df["pdb"].iloc[0]).parents[1]
        out_dir_name = out_dir.stem
        out_id = int(out_dir_name.split("_")[1])

        stat_df["out_id"] = out_id

        sample_df = pd.concat([sample_df, stat_df])

    sample_df.reset_index(drop=True, inplace=True)
    sample_df.to_csv(sample_file)


