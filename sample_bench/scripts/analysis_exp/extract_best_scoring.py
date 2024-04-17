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
from refine import refine, pool_refine
from utility import pool_read_pdb
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis")))
from score_analysis_exp import field_id_mapper, map_std_field_to_field


if __name__ == "__main__":
    exp_name = "169_N8"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    n_state = 8

    ref_exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/166_N1_ref")

    cif_map_df = pd.read_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"), index_col=0)

    columns = ["job_id", "cif_name"]
    cif_names = ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]
    for state in range(n_state):
        columns.append("w_{}".format(state))

    exp_stat_file = Path(Path.home(), "xray/sample_bench/data/7mhf/{}/ref_dataset.csv".format(exp_name))
    exp_stat_df = pd.DataFrame(columns=columns)
    for job_id in range(63):
        job_dir = Path(exp_dir, str(job_id))

        job_id = int(job_dir.name)
        n_cif = len(cif_map_df.loc[job_id, "cifs"].split(","))
        print(job_id, n_state)

        # bonus_fields = ["pdb"]

        stat_dfs = list()
        for cif_name in ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]:
            bonus_fields = ["pdb", "ff"]

            for state in range(n_state):
                bonus_fields.append(map_std_field_to_field("w_{}_{}".format(state, cif_name), job_id))

            field = map_std_field_to_field("r_free_{}".format(cif_name), job_id)

            if not field:
                continue

            print(field, bonus_fields)

            stat_df = get_stat_df(
                log_files=[Path(out_dir, "log.csv") for out_dir in job_dir.glob("*/")],
                field=field,
                N=100,
                bonus_fields=bonus_fields,
                equil=1,
                pdb_only=True,
                max_rmsd=None
            )

            # print(len(stat_df))
            # stat_dfs.append(stat_df)

            for i in range(len(stat_df)):
                row = len(exp_stat_df)
                exp_stat_df.loc[row, "pdb"] = stat_df.loc[i, "pdb"]
                exp_stat_df.loc[row, "job_id"] = job_id
                exp_stat_df.loc[row, "cif_name"] = cif_name

                for state in range(n_state):
                    w_field = map_std_field_to_field("w_{}_{}".format(state, cif_name), job_id)
                    exp_stat_df.loc[row, "w_{}".format(state)] = stat_df.loc[i, w_field]

    exp_stat_df.to_csv(exp_stat_file)

