from pathlib import Path
import random
import pandas as pd
random.seed(0)
import IMP
import IMP.atom
import sys
import multiprocessing
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df


if __name__ == "__main__":
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/183_test/812183")
    out_dirs = list(exp_dir.glob("output*"))
    log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]
    print(log_files)

    stat_df = get_stat_df(
        log_files=log_files,
        field="r_free_0",
        N=100,
        bonus_fields=["xray_0", "pdb", "ff", "rmsd_0", "w_0_0", "w_1_0"],
        equil=50,
        pdb_only=True
    )

    stat_df.to_csv(Path(Path.home(), "xray/tmp/stat_df_2.csv"))
