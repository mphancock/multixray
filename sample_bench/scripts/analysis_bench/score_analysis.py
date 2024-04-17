from pathlib import Path
import sys
import pandas as pd
import numpy as np
import argparse

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df



if __name__ == "__main__":
    columns = ["N", "J", "job", "xray_0", "xray_1", "xray_0+xray_1", "r_free_0", "r_free_1", "rmsd_0", "rmsd_1", "rmsd_0+rmsd_1", "w_0_0", "w_0_1", "w_1_0", "w_1_1", "w_2_0", "w_2_1", "w_3_0", "w_3_1", "pdb", "ff"]
    analysis_df = pd.DataFrame(columns=columns)
    analysis_df.astype({'job': 'int', 'N': 'int', 'J': 'int'})

    score_type = "rmsd"
    analysis_file = Path(Path.home(), "xray/sample_bench/data/7mhf/123_natives_2_state/score_analysis_{}.csv".format(score_type))
    analysis_file.parent.mkdir(exist_ok=True)

    params = [("125_natives_1_state", 1, 1), ("145_native_1_state_2_cond", 1, 2), ("123_natives_2_state", 2, 1), ("124_natives_2_cond", 2, 2), ("141_native_4_state_1_cond", 4, 1), ("142_native_4_state_2_cond", 4, 2)]
    # params = [("125_natives_1_state", 1, 1)]
    for exp_name, n_state, n_cond in params:
        exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
        job_dirs = list(Path(exp_dir).glob("*"))

        if n_cond == 1:
            min_field = "{}_0".format(score_type)
        elif n_cond == 2:
            min_field = "{}_0+{}_1".format(score_type, score_type)

        bonus_fields = ["xray_0", "xray_1", "xray_0+xray_1", "r_free_0", "r_free_1", "rmsd_0", "rmsd_1", "rmsd_0+rmsd_1", "pdb", "ff"]
        bonus_fields.remove(min_field)
        if n_cond == 1:
            bonus_fields.remove("xray_1")
            bonus_fields.remove("xray_0+xray_1")
            bonus_fields.remove("r_free_1")
            bonus_fields.remove("rmsd_1")
            bonus_fields.remove("rmsd_0+rmsd_1")

        for state in range(n_state):
            for cond in range(n_cond):
                bonus_fields.append("w_{}_{}".format(state, cond))

        # break

        for job_dir in job_dirs:
            print(job_dir, min_field, bonus_fields)
            job_log_files = list(job_dir.glob("output_*/log.csv"))

            stat_df = get_stat_df(
                log_files=job_log_files,
                field=min_field,
                N=1,
                bonus_fields=bonus_fields,
                equil=1,
                pdb_only=True
            )
            stat_df.loc[0, "job"] = int(job_dir.name)
            stat_df.loc[0, "N"] = n_state
            stat_df.loc[0, "J"] = n_cond
            stat_df = stat_df.astype({'job': 'int', 'N': 'int', 'J': 'int'})

            analysis_df = pd.concat([analysis_df, stat_df], join="outer")

    analysis_df.drop(columns=["index"], inplace=True)
    analysis_df = analysis_df.reset_index()
    analysis_df["rmsd_0+rmsd_1"] = analysis_df["rmsd_0+rmsd_1"] / 2
    analysis_df = analysis_df.sort_values(by=["N","J","job"])

    analysis_df.to_csv(analysis_file)
