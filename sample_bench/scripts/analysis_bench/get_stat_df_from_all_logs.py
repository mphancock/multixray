from pathlib import Path
import argparse
import sys
import pandas as pd
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df


if __name__ == "__main__":
    job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/62_7mhj_8/1018130")
    n_states = 8

    job_name = job_dir.parents[0].name
    target_name = job_dir.parents[1].name

    analysis_dir = Path(Path.home(), "xray/sample_bench/data", target_name, job_name)

    print(analysis_dir)
    analysis_dir.mkdir(exist_ok=True)

    output_dirs = list(job_dir.glob("output*"))
    log_files = [Path(output_dir, "log.csv") for output_dir in output_dirs]

    main_field = "r_free_0"

    if main_field == "xray_0":
        bonus_fields=["r_free_0", "r_work_0", "pdb"]
    if main_field == "r_free_0":
        bonus_fields=["xray_0", "r_work_0", "pdb"]

    bonus_fields.extend(["occ_{}".format(i) for i in range(n_states)])

    stat_df = get_stat_df.get_stat_df(
        log_file_groups=[[log_file] for log_file in log_files],
        main_field=main_field,
        main_stat="min",
        bonus_fields=bonus_fields,
        N=1,
        equil=1,
        pdb_only=True
    )
    stat_df.to_csv(Path(analysis_dir, "stat_df_{}.csv".format(main_field)))

