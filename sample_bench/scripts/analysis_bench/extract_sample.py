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
    exp_name = "187_bench_ref_10000"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)

    job_csv_file = Path(Path.home(), "xray/sample_bench/data/params/bench.csv")
    sample_file = Path(Path.home(), "xray/sample_bench/data/7mhf/{}/sample.csv".format(exp_name))
    sample_df = pd.DataFrame()

    for job_id in range(40):
        param_dict = read_job_csv(job_csv_file=job_csv_file, job_id=job_id)
        cif_files = param_dict["cifs"]
        J = param_dict["J"]
        N = param_dict["N"]

        job_cif_names = [Path(cif).stem for cif in cif_files]

        job_dir = Path(exp_dir, str(job_id))
        log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*")]

        bonus_fields = ["ff", "pdb"]
        if J == 1:
            field = "xray_0"
            cif_names = [0]
            bonus_fields.extend(["r_free_0"])
        else:
            field = "xray_0+xray_1"
            cif_names = [0, 1]
            bonus_fields.extend(["xray_0", "xray_1", "r_free_0", "r_free_1"])

        for cif_name in cif_names:
            for state in range(N):
                bonus_fields.append("w_{}_{}".format(state, cif_name))

        # bonus_cif_fields = ["r_free_{}".format(cif_name) for cif_name in job_cif_names]
        # bonus_cif_fields.remove(field)

        stat_df = get_stat_df(
            log_files=log_files,
            field=field,
            N=1,
            bonus_fields=bonus_fields,
            equil=5,
            pdb_only=True,
            max_ff=None
        )
        # print(cif_name, stat_df["r_free_{}".format(cif_name)][0])
        # stat_df["cif_name"] = cif_name
        stat_df["N"] = N
        stat_df["J"] = J
        stat_df["job_id"] = job_id

        # job_cif_str = ",".join(job_cif_names)
        # stat_df["cifs"] = job_cif_str
        stat_df.drop(columns=["index"], inplace=True)

        # stat_df.rename(columns={"r_free_0": "r_free"}, inplace=True)

        # if J == 1:
        #     stat_df.rename(columns={"r_free_0".format(cif_name): "r_free"}, inplace=True)
        # else:
        #     stat_df.rename(columns={"r_free_0+r_free_1": "r_free"}, inplace=True)
        #     stat_df["r_free"] = stat_df["r_free"] / 2

        # for state in range(N):
        #     stat_df.rename(columns={"w_{}_{}".format(state, cif_name): "w_{}".format(state)}, inplace=True)

        print(stat_df.columns)
        sample_df = pd.concat([sample_df, stat_df])

    sample_df.reset_index(drop=True, inplace=True)
    sample_df.to_csv(sample_file)


