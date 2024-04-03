import pandas as pd
from pathlib import Path
import multiprocessing
import numpy as np


def pool_get_max_time(log_file):
    try:
        log_df = pd.read_csv(log_file, index_col=0)
    except FileNotFoundError:
        print("File not found: {}".format(log_file))
        return None

    speed = log_df.iloc[-1]["time"] / log_df.iloc[-1]["step"]
    return speed


if __name__ == "__main__":
    pool_params = list()
    for i in range(500):
        log_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/167_N2/62/output_{}/log.csv".format(i))
        pool_params.append(log_file)

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(pool_get_max_time, pool_params)

    speeds = list()
    for speed in pool_results:
        if speed:
            speeds.append(speed)

    print(np.mean(speed))
    pool_obj.close()
