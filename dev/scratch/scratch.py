from pathlib import Path
import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import seaborn as sns
sns.set_theme()
import math
import pandas as pd
import random
import numpy as np
import sys

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import refine_hs_to_max_ff


if __name__ == "__main__":
    log_dfs = list()
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/194_free_control_all_N/4")

    for out_dir in exp_dir.glob("output*"):
        log_df = pd.read_csv(Path(out_dir, "log.csv"))
        log_df = log_df[log_df['pdb'].notna()]
        log_dfs.append(log_df)

    print(len(log_dfs))

    all_log_df = pd.concat(log_dfs)
    # print(all_log_df.columns)

    print(all_log_df["r_free_7mhi"].min())
    print(all_log_df["r_work_7mhi"].min())
