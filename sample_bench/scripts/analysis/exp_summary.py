from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from stat_df import get_stat_df
from utility import pool_read_pdb
from params import read_job_csv
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis_exp")))
# from score_analysis_exp import field_id_mapper, map_std_field_to_field


if __name__ == "__main__":
    exp_name = "281_exp_4_test"
    job_csv_file = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/params/281.csv")
    params_df = pd.read_csv(job_csv_file, index_col=0)

    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out", exp_name)
    analysis_dir = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/analysis", exp_name)
    analysis_dir.mkdir(exist_ok=True)

    summary_df = pd.DataFrame()

    ## 1. ITERATE THROUGH EACH JOB AND FIND BEST MODEL
    ## if joint then the analysis is performed for collective satisfaction else the analysis is performed for each individual cif file for each job
    joint = False
    row = 0
    for i in range(len(params_df)):
        job_dir = Path(exp_dir, str(i))
        print(job_dir)
        param_dict = read_job_csv(job_csv_file=job_csv_file, job_id=i)

        cif_files = param_dict["cifs"]
        cif_names = [cif_file.stem for cif_file in cif_files if cif_file]

        if len(cif_names) == 0:
            continue

        ## fields is either every single individual r free or the collective r free
        fields = list()
        if not joint:
            for cif_name in cif_names:
                fields.append("r_free_{}".format(cif_name))
        else:
            field = "r_free_{}".format(cif_names[0])
            for j in range(1, len(cif_names)):
                field += "+r_free_{}".format(cif_names[j])

        ## 2. ITERATE THROUGH EACH CIF FILE (OR COMBO OF CIF FILES) TO FIND BEST MODEL
        ## for all fields perform the analysis
        for field in fields:
            ## Perform the analysis based on the first r free
            bonus_fields = list()
            bonus_fields.append("rmsd")
            for j in range(len(cif_names)):
                bonus_fields.append("r_free_{}".format(cif_names[j]))
                bonus_fields.append("r_work_{}".format(cif_names[j]))
                bonus_fields.append("xray_{}".format(cif_names[j]))
                # bonus_fields.append("rmsd_{}".format(cif_names[j]))

            bonus_fields.extend(["ff", "pdb"])
            if field in bonus_fields:
                bonus_fields.remove(field)

            N, J = param_dict["N"], param_dict["J"]
            w_cols = list()
            for state in range(N):
                for cond in range(J):
                    w_cols.append("w_{}_{}".format(state, cif_names[cond]))

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

            print(row, len(log_files))
            print(field, bonus_fields)
            stat_df = get_stat_df(
                log_files=log_files,
                field=field,
                N=1,
                bonus_fields=bonus_fields,
                equil=0,
                pdb_only=True,
                max_rmsd=None
            )

            if len(stat_df) == 0:
                continue

            # df_id = len(w_xray_df)
            summary_df.loc[row, "N"] = N
            summary_df.loc[row, "J"] = J
            summary_df.loc[row, "w_xray"] = param_dict["w_xray"]
            summary_df.loc[row, field] = stat_df.loc[0, field]

            for bonus_field in bonus_fields:
                summary_df.loc[row, bonus_field] = stat_df.loc[0, bonus_field]

            for j in range(len(cif_files)):
                summary_df.loc[row, "cif_{}".format(j)] = cif_files[j]

            row += 1
            print("len", len(summary_df))

    summary_df.to_csv(Path(analysis_dir, "summary.csv"))
