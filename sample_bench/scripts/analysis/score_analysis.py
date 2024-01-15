from pathlib import Path
import sys
import pandas as pd
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df


if __name__ == "__main__":
    job_name = "152_native_1_cif"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7", job_name)
    analysis_dir = Path(Path.home(), "xray/sample_bench/data/3ca7", job_name)
    analysis_dir.mkdir(exist_ok=True)
    log_file = Path(analysis_dir, "score_analysis.csv".format(job_name))
    n_jobs = 10

    n_cif = 1
    field = "xray_0"
    bonus_fields = []
    for i in range(1,n_cif):
        field += "+xray_{}".format(i)
        bonus_fields.append("xray_{}".format(i))

    for i in range(n_cif):
        bonus_fields.append("r_free_{}".format(i))

    for i in range(n_cif):
        bonus_fields.append("rmsd_avg_{}".format(i))

    bonus_fields.append("pdb")

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

        log_df.loc[job_id, "min_{}".format(field)] = stat_df["{}_min_0".format(field)].min()
        # log_df.loc[job_id, "avg_min_{}".format(field)] = stat_df["{}_min_0".format(field)].mean()

        stat_df["{}_min_0".format(field)] = pd.to_numeric(stat_df["{}_min_0".format(field)])

        field_min_entry = stat_df.loc[stat_df["{}_min_0".format(field)].idxmin()]

        for bonus_field in bonus_fields:
            log_df.loc[job_id, "{}_min_0_{}".format(field, bonus_field)] = field_min_entry["{}_min_0_{}".format(field, bonus_field)]

    log_df.to_csv(log_file)