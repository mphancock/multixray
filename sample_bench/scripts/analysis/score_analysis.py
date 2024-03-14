from pathlib import Path
import sys
import pandas as pd
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df


def get_score_analysis_df(
        exp_dir,
        n_jobs,
        field,
        bonus_fields
):
    log_df = pd.DataFrame(index=list(range(n_jobs)))
    for job_id in range(n_jobs):
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
    n_jobs = 10
    # n_cond = 1
    # job_name = "123_natives_2_state"

    # params = [("123_natives_2_state", 1), ("125_natives_1_state", 1), ("141_native_4_state_1_cond", 1), ("151_native_N8_J1", 1)]
    # params = [("124_natives_2_cond", 2), ("142_native_4_state_2_cond", 2), ("145_native_1_state_2_cond", 2), ("152_native_N8_J2", 2)]
    for job_name, n_cond in [("141_native_4_state_1_cond", 1)]:
        for stat_type in ["rmsd", "xray"]:
            exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", job_name)
            analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", job_name)
            analysis_dir.mkdir(exist_ok=True)

            bonus_fields = ["pdb", "ff"]
            field = ""
            for i in range(n_cond):
                field = field + stat_type + "_{}".format(i)
                if i < n_cond-1:
                    field = field + "+"

                bonus_fields.append("xray_{}".format(i))
                bonus_fields.append("r_free_{}".format(i))
                bonus_fields.append("rmsd_{}".format(i))

            ## BUG FIX: Requesting the field as a bonus field causes issues with get_stat_df.
            if field in bonus_fields:
                bonus_fields.remove(field)

            log_df = get_score_analysis_df(
                exp_dir=exp_dir,
                n_jobs=n_jobs,
                field=field,
                bonus_fields=bonus_fields
            )

            log_file = Path(analysis_dir, "score_analysis_{}.csv".format(field))
            log_df.to_csv(log_file)
