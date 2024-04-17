from pathlib import Path
import sys
import pandas as pd
import numpy as np
import argparse

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df


if __name__ == "__main__":
    exp_name = "175_w_xray"
    n_cond = 1

    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)

    # min_field = "r_work_0"
    # bonus_fields = ["r_free_0", "ff", "pdb"]

    min_field = "r_free_0"
    bonus_fields = ["r_work_0", "ff", "pdb"]

    fields = [min_field] + bonus_fields
    columns = ["N", "W"] + fields
    w_xray_df = pd.DataFrame(columns=columns)

    for i in range(5):
        for j in range(11):
            job_dir = Path(exp_dir, "N{}_W{}".format(i, j))

            log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*/")]

            # print(log_files)

            stat_df = get_stat_df(
                log_files=log_files,
                field=min_field,
                N=1,
                bonus_fields=bonus_fields,
                equil=1,
                pdb_only=True,
                max_rmsd=None
            )

            job_num = i*11+j
            print(job_dir, job_num)

            w_xray_df.loc[job_num, "N"] = i
            w_xray_df.loc[job_num, "W"] = j
            for field in fields:
                w_xray_df.loc[job_num,field] = stat_df.iloc[0][field]

            # print(w_xray_df.head())

    analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", exp_name)
    analysis_dir.mkdir(exist_ok=True)

    log_file = Path(analysis_dir, "wxray_{}.csv".format(min_field))
    w_xray_df.to_csv(log_file)

