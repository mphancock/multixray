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


def get_best_from_out_dir(
    params
):
    out_dir = params["out_dir"]
    field = params["field"]
    bonus_fields = params["bonus_fields"]
    log_files = [Path(out_dir, "log.csv")]

    params_2 = dict()
    params_2["log_files"] = log_files
    params_2["equil"] = 50
    params_2["field"] = field
    params_2["bonus_fields"] = bonus_fields
    params_2["N"] = 1
    params_2["pdb_only"] = True
    params_2["max_rmsd"] = None
    stat_df = pool_get_stat_info_df(params_2)

    return stat_df


if __name__ == "__main__":
    exp_name = "182_bench"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)

    n_states = [1, 2]

    cif_map_df = pd.read_csv(Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/csvs/7mhf_30.csv"), index_col=0)
    sample_file = Path(Path.home(), "xray/sample_bench/data/7mhf/{}/sample_per_out.csv".format(exp_name))

    pool_params = list()
    for n in range(len(n_states)):
        N = n_states[n]
        for j in range(20):
            cif_name = "0"
            job_cif_files = cif_map_df.loc[j, "cifs"].split(",")
            job_cif_names = [Path(cif).stem for cif in job_cif_files]
            J = len(job_cif_names)

            job_dir = Path(exp_dir, "N{}_J{}".format(n, j))
            out_dirs = [out_dir for out_dir in job_dir.glob("output*")]

            field = "r_free_{}".format(cif_name)
            bonus_fields = ["ff", "rmsd_0", "pdb"]

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
        out_dir = Path(stat_df["pdb"].iloc[0]).parents[1]
        out_dir_name = out_dir.stem
        out_dir_id = int(out_dir_name.split("_")[1])

        job_dir = out_dir.parents[0]
        # print("job_dir", job_dir)
        n = int(job_dir.stem.split("_")[0][1:])
        N = n_states[n]
        j = int(job_dir.stem.split("_")[1][1:])

        print(n, j, out_dir_id)

        stat_df["cif_name"] = cif_name
        stat_df["n"] = n
        stat_df["N"] = N
        stat_df["j"] = j

        job_cif_files = cif_map_df.loc[j].values[0].split(",")
        job_cif_names = [Path(job_cif_file).stem for job_cif_file in job_cif_files]
        job_cif_str = ",".join(job_cif_names)
        stat_df["J"] = len(job_cif_names)
        stat_df["cifs"] = job_cif_str

        stat_df["out_dir_id"] = out_dir_id

        # stat_df.drop(columns=["index"], inplace=True)

        stat_df.rename(columns={"r_free_0": "r_free"}, inplace=True)
        stat_df.rename(columns={"rmsd_0": "rmsd"}, inplace=True)

        for state in range(N):
            stat_df.rename(columns={"w_{}_{}".format(state, cif_name): "w_{}".format(state)}, inplace=True)

        print(stat_df.columns)
        sample_df = pd.concat([sample_df, stat_df])

    sample_df.reset_index(drop=True, inplace=True)
    sample_df.to_csv(sample_file)


