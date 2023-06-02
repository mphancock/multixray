from pathlib import Path
import pandas as pd
import numpy as np
import shutil
import sys

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df


if __name__ == "__main__":
    job_dir = Path(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhk/09_1_h20/3608081"))
    out_dirs = list(job_dir.glob("output*"))

    log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]
    stat_df = get_stat_df.get_stat_df(
        log_file_groups=[log_files],
        fields=["r_free"],
        stats=["min"],
        N=1000,
        offset=1,
        equil=1000,
        test=False
    )

    r_frees = list()
    for i in range(1000):
        r_frees.append(stat_df["r_free_min_{}".format(i)].iloc[0])

    print(np.min(r_frees))
