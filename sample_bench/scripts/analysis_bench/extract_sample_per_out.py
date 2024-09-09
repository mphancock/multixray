from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df, pool_get_stat_info_df
sys.path.append(str(Path(Path.home(), "xray/src")))
# from refine import refine, pool_refine
from utility import pool_read_pdb
from params import read_job_csv


def get_best_from_out_dir(
    params
):
    out_dir = params["out_dir"]
    field = params["field"]
    bonus_fields = params["bonus_fields"]
    log_files = [Path(out_dir, "log.csv")]

    params = dict()
    params["log_files"] = log_files
    params["equil"] = 1
    params["field"] = field
    params["bonus_fields"] = bonus_fields
    params["N"] = 1
    params["pdb_only"] = True
    params["max_rmsd"] = None
    params["max_ff"] = None
    stat_df = pool_get_stat_info_df(params)

    return stat_df


if __name__ == "__main__":
    exp_name = "187_bench_ref_10000"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    # field_type = "r_free"

    field_1 = "xray_0"
    field_2 = "xray_0+xray_1"

    job_csv_file = Path(Path.home(), "xray/sample_bench/data/params/bench.csv")
    sample_file = Path(Path.home(), "xray/sample_bench/data/7mhf/{}/sample_per_out.csv".format(exp_name))

    pool_params = list()

    param_dicts = list()
    for job_id in range(40):
        param_dict = read_job_csv(job_csv_file=job_csv_file, job_id=job_id)
        param_dicts.append(param_dict)

    for job_id in range(40):
        param_dict = param_dicts[job_id]

        cif_files = param_dict["cifs"]
        cif_names = [Path(cif).stem for cif in cif_files]
        J = param_dict["J"]
        N = param_dict["N"]

        job_dir = Path(exp_dir, str(job_id))
        out_dirs = [out_dir for out_dir in job_dir.glob("output*")]

        bonus_fields = ["ff", "pdb", "r_free_0"]
        if J == 1:
            field = field_1
        else:
            field = field_2
            bonus_fields.append("r_free_1")

        w_cols = list()
        for state in range(N):
            w_cols.append("w_{}_{}".format(state, 0))
        bonus_fields.extend(w_cols)

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

        # print(len(stat_df))
        out_dir = Path(stat_df["pdb"].iloc[0]).parents[1]
        out_dir_name = out_dir.stem
        out_id = int(out_dir_name.split("_")[1])

        job_dir = out_dir.parents[0]
        job_id = int(job_dir.stem)

        # print(job_name, out_dir_name)
        # print("job_dir", job_dir)

        param_dict = param_dicts[job_id]
        N = param_dict["N"]
        J = param_dict["J"]
        stat_df["N"] = N
        stat_df["job_id"] = job_id
        # stat_df["cif_name"] = cif_name

        # job_cif_files = cif_map_df.loc[j].values[0].split(",")
        # job_cif_names = [Path(job_cif_file).stem for job_cif_file in job_cif_files]
        # job_cif_str = ",".join(job_cif_names)
        # stat_df["J"] = len(job_cif_names)
        # stat_df["cifs"] = job_cif_str

        stat_df["out_id"] = out_id

        # stat_df.drop(columns=["index"], inplace=True)

        if J == 1:
            stat_df.rename(columns={field_1: "field"}, inplace=True)
        else:
            stat_df.rename(columns={field_2: "field"}, inplace=True)
            # print(stat_df.columns)
            # stat_df["field"] = stat_df["field"] / 2

        # stat_df.rename(columns={"rmsd_0": "rmsd"}, inplace=True)

        for state in range(N):
            stat_df.rename(columns={"w_{}_{}".format(state, 0): "w_{}".format(state)}, inplace=True)

        # print(stat_df.columns)
        sample_df = pd.concat([sample_df, stat_df])

    sample_df.reset_index(drop=True, inplace=True)
    sample_df.to_csv(sample_file)


