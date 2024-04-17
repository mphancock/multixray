from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
# from refine import refine, pool_refine
from utility import pool_read_pdb
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis_exp")))
from score_analysis_exp import field_id_mapper, map_std_field_to_field


if __name__ == "__main__":
    exp_name = "176_w_xray_7mhh"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    n_state = 8

    ref_exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/{}_ref".format(exp_name))

    cif_map_df = pd.read_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"), index_col=0)

    columns = ["cif_name"]
    cif_names = ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]
    for state in range(n_state):
        columns.append("w_{}".format(state))

    analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", exp_name)
    analysis_dir.mkdir(exist_ok=True)
    exp_stat_file = Path(analysis_dir, "sample.csv".format(exp_name))
    exp_stat_df = pd.DataFrame(columns=columns)

    Ns = [1, 2, 4, 8, 16]
    w_xrays = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

    job_id = 2

    # for i in range(len(Ns)):
    for N_id in [3]:
        for w_xray_id in range(len(w_xrays)):
            job_dir = Path(exp_dir, "N{}_W{}".format(N_id, w_xray_id))

            stat_dfs = list()
            for cif_name in ["7mhh"]:
                bonus_fields = ["pdb", "ff"]

                for state in range(n_state):
                    bonus_fields.append("w_{}_{}".format(state, 0))

                # for state in range(n_state):
                #     bonus_fields.append(map_std_field_to_field("w_{}_{}".format(state, cif_name), job_id))

                # field = map_std_field_to_field("r_free_{}".format(cif_name), job_id)

                field = "r_free_0"
                print(job_dir)
                print(N_id, w_xray_id, bonus_fields)

                # if not field:
                #     continue


                stat_df = get_stat_df(
                    log_files=[Path(out_dir, "log.csv") for out_dir in job_dir.glob("*/")],
                    field=field,
                    N=100,
                    bonus_fields=bonus_fields,
                    equil=1,
                    pdb_only=True,
                    max_rmsd=None
                )

                # print(i, w_xrays[j], stat_df.loc[0, "r_free_0"])

                print(len(stat_df))
                stat_dfs.append(stat_df)

                for i in range(len(stat_df)):
                    row = len(exp_stat_df)
                    exp_stat_df.loc[row, "w_xray_id"] = w_xray_id
                    exp_stat_df.loc[row, "r_free_0"] = stat_df.loc[i, "r_free_0"]
                    exp_stat_df.loc[row, "pdb"] = stat_df.loc[i, "pdb"]
                    exp_stat_df.loc[row, "cif_name"] = cif_name

                    for state in range(n_state):
                        w_field = "w_{}_{}".format(state, 0)
                        exp_stat_df.loc[row, "w_{}".format(state)] = stat_df.loc[i, w_field]

    exp_stat_df.to_csv(exp_stat_file)

