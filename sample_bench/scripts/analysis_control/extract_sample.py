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


if __name__ == "__main__":
    exp_name = "181_wxray_control"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)

    n_states = [1, 2, 4, 8, 16]

    cif_map_df = pd.read_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"), index_col=0)

    cif_names = ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]
    sample_file = Path(Path.home(), "xray/sample_bench/data/7mhf/{}/sample.csv".format(exp_name))
    sample_df = pd.DataFrame()

    for state_id in range(len(n_states)):
        N = n_states[state_id]
        for job_id in range(6):
            job_cif_files = cif_map_df.loc[job_id, "cifs"].split(",")
            job_cif_names = [Path(cif).stem for cif in job_cif_files]
            cif_name = job_cif_names[0]

            for M_id in range(6):
                M = M_id + 1
                print(state_id, job_id, M_id)
                print(cif_name)

                job_dir = Path(exp_dir, "N{}_J{}_M{}".format(state_id, job_id, M_id))
                log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*")]

                field = "r_free_{}".format(cif_name)
                bonus_fields = ["ff", "pdb"]

                for state in range(N):
                    bonus_fields.append("w_{}_{}".format(state, cif_name))

                stat_df = get_stat_df(
                    log_files=log_files,
                    field=field,
                    N=50,
                    bonus_fields=bonus_fields,
                    equil=50,
                    pdb_only=True,
                    max_rmsd=None
                )
                stat_df["cif_name"] = cif_name
                stat_df["N"] = N
                stat_df["job_id"] = job_id
                stat_df["M"] = M

                stat_df["cifs"] = cif_name
                stat_df.drop(columns=["index"], inplace=True)

                stat_df.rename(columns={"r_free_{}".format(cif_name): "r_free"}, inplace=True)

                for state in range(N):
                    stat_df.rename(columns={"w_{}_{}".format(state, cif_name): "w_{}".format(state)}, inplace=True)

                print(stat_df.columns)
                sample_df = pd.concat([sample_df, stat_df])

    sample_df.reset_index(drop=True, inplace=True)
    sample_df.to_csv(sample_file)


