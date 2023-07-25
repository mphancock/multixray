from pathlib import Path
import argparse
import sys
import pandas as pd
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df


if __name__ == "__main__":
    job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/57_sb_2x/9692409")
    output_dirs = list(job_dir.glob("output*"))
    log_files = [Path(output_dir, "log.csv") for output_dir in output_dirs]

    stat_df = get_stat_df.get_stat_df(
        log_file_groups=[[log_file] for log_file in log_files],
        main_field="xray_0",
        main_stat="min",
        bonus_fields=["rmsd", "r_free_0", "occ_0", "occ_1", "pdb"],
        N=1,
        equil=1,
        pdb_only=True
    )
    stat_df.to_csv(Path(Path.home(), "xray/tmp/stat_df.csv"))

