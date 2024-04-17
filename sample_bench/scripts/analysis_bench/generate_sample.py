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
import get_stat_df


def get_sample(
        log_files,
        N,
        field,
        bonus_fields
):

    for i in range(n_state):
        for j in range(n_cond):
            bonus_fields.append("w_{}_{}".format(i, j))

    stat_df = get_stat_df.get_stat_df(
        log_file_groups=[log_files],
        main_field=field,
        main_stat="min",
        N=N,
        bonus_fields=bonus_fields,
        equil=1,
        pdb_only=True,
        rmsd_filter=10,
        test=False
    )

    best_sample_df = pd.DataFrame(index=range(N), columns=[field]+bonus_fields)
    n_fields = len(bonus_fields)+1
    for i in range(N):
        scores = stat_df.iloc[0,i*n_fields:(i+1)*n_fields]
        best_sample_df.iloc[i][field] = scores.iloc[0]

        for j in range(len(bonus_fields)):
            best_sample_df.iloc[i][bonus_fields[j]] = scores.iloc[j+1]

    return best_sample_df

if __name__ == "__main__":
    target = "7mhf"

    n_cond = 3
    N = 100
    stat_type="xray"

    params = [("162_J3_ijk", 0, 1), ("162_J3_ijk", 1, 2), ("162_J3_ijk", 2, 4), ("162_J3_ijk", 3, 8)]
    # params = [("162_J3_ijk", 1, 2)]
    for job_name, job_id, n_state in params:
        pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/{}/{}/{}".format(target, job_name, job_id))
        print(pdb_dir)
        out_dirs = list(pdb_dir.glob("output_*"))
        log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]
        # log_files = [Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/162_J3_ijk/1/output_757/log.csv")]
        print(log_files)

        data_dir = Path(Path.home(), "xray/sample_bench/data", target, job_name)
        data_dir.mkdir(exist_ok=True)

        bonus_fields = ["pdb", "ff"]
        # field = "r_free_0"
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

        print(field, bonus_fields)

        best_sample_file = Path(data_dir, "sample_{}_{}.csv".format(job_id, field))
        best_sample_df = get_sample(
            log_files=log_files,
            N=N,
            field=field,
            bonus_fields=bonus_fields
        )

        best_sample_df.to_csv(best_sample_file)
