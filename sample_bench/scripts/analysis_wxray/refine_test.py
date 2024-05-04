from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import argparse

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import read_pdb_and_refine, read_pdb_and_refine_to_max_ff
from utility import pool_read_pdb, get_n_state_from_pdb_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    w_xray_df = pd.read_csv(Path(Path.home(), "xray/sample_bench/data/7mhf/180_wxray_bench/best_wxray.csv"), index_col=0)

    pdb_files = [Path(pdb_file) for pdb_file in w_xray_df["pdb"]]
    print(len(pdb_files))
    pool_params = list()
    for i in range(len(pdb_files)):
        params_dict = dict()
        params_dict["pdb_file"] = pdb_files[i]
        params_dict["out_pdb_file"] = Path(Path.home(), "xray/tmp/refine/{}.pdb".format(i))
        params_dict["max_ff"] = 10000
        params_dict["log_file"] = None
        pool_params.append(params_dict)

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(read_pdb_and_refine_to_max_ff, pool_params)

    pdb_id = 0
    for pool_result in pool_results:
        # Make sure that the pdb file existed
        print(pdb_id, pool_result)
        pdb_id += 1

