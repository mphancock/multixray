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
from params import read_job_csv


if __name__ == "__main__":
    exp_name = "253_7mhf"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out", exp_name)

    job_csv_file = Path(Path.home(), "xray/sample_bench/data/params/253.csv")
    sample_file = Path(Path.home(), "xray/sample_bench/data/analysis/{}/sample.csv".format(exp_name))
    sample_df = pd.DataFrame()

    for job_id in range(5):
        print(job_id)
        param_dict = read_job_csv(job_csv_file=job_csv_file, job_id=job_id)
        cif_files = param_dict["cifs"]
        J = param_dict["J"]
        N = param_dict["N"]

        job_cif_names = [Path(cif).stem for cif in cif_files]

        job_dir = Path(exp_dir, str(job_id))

        print(job_dir)
        log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*")]

        for cif_name in job_cif_names:
            field = "r_free_{}".format(cif_name)
            bonus_fields = ["ff", "pdb", "xray_{}".format(cif_name)]

            for state in range(N):
                bonus_fields.append("w_{}_{}".format(state, cif_name))

            stat_df = get_stat_df(
                log_files=log_files,
                field=field,
                N=1,
                bonus_fields=bonus_fields,
                equil=5,
                pdb_only=True,
                max_ff=None
            )

            stat_df["N"] = N
            stat_df["J"] = J
            stat_df["job_id"] = job_id

            stat_df.drop(columns=["index"], inplace=True)

            print(stat_df.columns)
            sample_df = pd.concat([sample_df, stat_df])

    sample_df.reset_index(drop=True, inplace=True)
    sample_df.to_csv(sample_file)


