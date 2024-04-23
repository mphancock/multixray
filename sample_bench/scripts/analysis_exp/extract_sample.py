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
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis")))
from score_analysis_exp import field_id_mapper, map_std_field_to_field


if __name__ == "__main__":
    exp_name = "179_exp"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)

    n_states = [1, 2, 4, 8, 16]

    cif_map_df = pd.read_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"), index_col=0)

    cif_names = ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]
    sample_file = Path(Path.home(), "xray/sample_bench/data/7mhf/{}/sample.csv".format(exp_name))
    sample_df = pd.DataFrame()

    for state_id in range(len(n_states)):
        N = n_states[state_id]
        for job_id in range(63):
            print(state_id, job_id)
            job_cif_files = cif_map_df.loc[job_id, "cifs"].split(",")
            job_cif_names = [Path(cif).stem for cif in job_cif_files]
            J = len(job_cif_names)

            job_dir = Path(exp_dir, "N{}_J{}".format(state_id, job_id))
            log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*")]

            for cif_name in cif_names:
                if cif_name not in job_cif_names:
                    continue

                field = "r_free_{}".format(cif_name)
                bonus_fields = ["ff", "pdb"]

                for state in range(N):
                    bonus_fields.append("w_{}_{}".format(state, cif_name))

                bonus_cif_fields = ["r_free_{}".format(cif_name) for cif_name in job_cif_names]
                bonus_cif_fields.remove(field)
                # bonus_fields.extend(bonus_cif_fields)

                print(cif_name)

                stat_df = get_stat_df(
                    log_files=log_files,
                    field=field,
                    N=100,
                    bonus_fields=bonus_fields,
                    equil=50,
                    pdb_only=True,
                    max_rmsd=None
                )
                print(cif_name)
                stat_df["cif_name"] = cif_name
                stat_df["N"] = N
                stat_df["J"] = J
                stat_df["job_id"] = job_id

                job_cif_files = cif_map_df.loc[job_id].values[0].split(",")
                job_cif_names = [Path(job_cif_file).stem for job_cif_file in job_cif_files]
                job_cif_str = ",".join(job_cif_names)
                stat_df["cifs"] = job_cif_str

                stat_df.rename(columns={"r_free_{}".format(cif_name): "r_free"}, inplace=True)

                for state in range(N):
                    stat_df.rename(columns={"w_{}_{}".format(state, cif_name): "w_{}".format(state)}, inplace=True)

                print(stat_df.columns)
                sample_df = pd.concat([sample_df, stat_df])

    sample_df.reset_index(drop=True, inplace=True)
    sample_df.to_csv(sample_file)


