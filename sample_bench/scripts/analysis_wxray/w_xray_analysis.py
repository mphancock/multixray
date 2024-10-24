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
from utility import pool_read_pdb
from params import read_job_csv
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis_exp")))
# from score_analysis_exp import field_id_mapper, map_std_field_to_field


if __name__ == "__main__":
    exp_name = "243_3k0n"
    job_csv_file = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/params/243.csv")
    params_df = pd.read_csv(job_csv_file, index_col=0)

    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out", exp_name)
    analysis_dir = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/analysis", exp_name)
    analysis_dir.mkdir(exist_ok=True)

    w_xray_df = pd.DataFrame()
    for i in range(len(params_df)):
        job_dir = Path(exp_dir, str(i))
        print(job_dir)
        param_dict = read_job_csv(job_csv_file=job_csv_file, job_id=i)

        cif_files = param_dict["cifs"]
        cif_names = [cif_file.stem for cif_file in cif_files if cif_file]

        if len(cif_names) == 0:
            continue

        ## Perform the analysis based on the first r free
        field = "r_free_{}".format(cif_names[0])
        bonus_fields = ["r_work_{}".format(cif_names[0])]
        for j in range(1, len(cif_names)):
            bonus_fields.append("r_free_{}".format(cif_names[j]))
            bonus_fields.append("r_work_{}".format(cif_names[j]))
        bonus_fields.extend(["ff", "pdb"])

        N, J = param_dict["N"], param_dict["J"]
        w_cols = list()
        for state in range(N):
            w_cols.append("w_{}_{}".format(state, cif_names[0]))
        bonus_fields.extend(w_cols)

        # log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*/")]
        # print(log_files)
        log_files = list()
        for out_dir in job_dir.glob("*/"):
            log_file = Path(out_dir, "log.csv")
            if log_file.exists():
                log_files.append(log_file)

        if len(log_files) == 0:
            continue

        print(i, len(log_files))

        stat_df = get_stat_df(
            log_files=log_files,
            field=field,
            N=1,
            bonus_fields=bonus_fields,
            equil=1,
            pdb_only=True,
            max_rmsd=None
        )

        if len(stat_df) == 0:
            continue

        # df_id = len(w_xray_df)
        w_xray_df.loc[i, "N"] = N
        w_xray_df.loc[i, "J"] = J
        w_xray_df.loc[i, "w_xray"] = param_dict["w_xray"]
        w_xray_df.loc[i, "r_free"] = stat_df.loc[0, field]

        for w_col in w_cols:
            w_xray_df.loc[i, w_col] = stat_df.loc[0, w_col]

        for bonus_field in bonus_fields:
            w_xray_df.loc[i, bonus_field] = stat_df.loc[0, bonus_field]

        for j in range(len(cif_files)):
            w_xray_df.loc[i, "cif_{}".format(j)] = cif_files[j]

    w_xray_df.to_csv(Path(analysis_dir, "w_xray.csv"))
