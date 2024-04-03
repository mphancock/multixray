from pathlib import Path
import sys
import pandas as pd
import numpy as np
import argparse

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df


def get_score_analysis_df(
        exp_dir,
        job_ids,
        field,
        bonus_fields
):
    log_df = pd.DataFrame(index=job_ids)
    for job_id in job_ids:
        job_dir = Path(exp_dir, str(job_id))
        print(job_dir)
        out_dirs = list(job_dir.glob("output_*"))
        log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]

        # Get the average min R free and min R free from all output logs for a given w_xray.
        # R_free_0 is the first R_free reported.

        stat_df = get_stat_df.get_stat_df(
            log_file_groups=[[log_file] for log_file in log_files],
            main_field=field,
            main_stat="min",
            bonus_fields=bonus_fields,
            N=10,
            equil=1,
            pdb_only=True,
            rmsd_filter=10,
            test=False
        )

        # stat_df.to_csv(Path(Path.home(), "xray/tmp/stat_df.csv"))
        stat_df.dropna(inplace=True)
        print(len(stat_df))

        log_df.loc[job_id, field] = stat_df["{}_min_0".format(field)].min()

        # stat_df["{}_min_0".format(field)] = pd.to_numeric(stat_df["{}_min_0".format(field)])

        field_min_entry = stat_df.loc[stat_df["{}_min_0".format(field)].idxmin()]

        for bonus_field in bonus_fields:
            log_df.loc[job_id, bonus_field] = field_min_entry["{}_min_0_{}".format(field, bonus_field)]

    return log_df



if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--n_cond", type=int)
    # parser.add_argument("--job_name")
    # args = parser.parse_args()

    # params = [("165_J1_i", 1), ("164_J2_ik", 2), ("162_J3_ijk", 3), ("163_J4_fijk", 4)]
    # n_cond = args.n_cond
    # job_name = args.job_name


    for job_name, n_cond in [("143_native_1_state_2_cond_wxray", 2)]:
        fields = list()
        # fields = ["r_free_{}".format(i) for i in range(n_cond)]
        field = ""
        for i in range(n_cond):
            field = field + "xray_{}".format(i)
            if i < n_cond-1:
                field = field + "+"
        fields.append(field)

        bonus_fields = ["pdb", "ff"]

        for i in  range(n_cond):
            bonus_fields.append("xray_{}".format(i))
            bonus_fields.append("r_free_{}".format(i))
            bonus_fields.append("rmsd_{}".format(i))

        for field in fields:
            bonus_fields_tmp = bonus_fields.copy()
            if field in bonus_fields_tmp:
                bonus_fields_tmp.remove(field)

            exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", job_name)
            analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", job_name)
            analysis_dir.mkdir(exist_ok=True)

            print(field, bonus_fields_tmp)
            log_df = get_score_analysis_df(
                exp_dir=exp_dir,
                job_ids=list(range(8)),
                field=field,
                bonus_fields=bonus_fields_tmp
            )

            log_file = Path(analysis_dir, "score_{}.csv".format(field))
            log_df.to_csv(log_file)
