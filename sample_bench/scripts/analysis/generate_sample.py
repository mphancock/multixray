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


if __name__ == "__main__":
    target = "7mhf"

    for exp_num, job_num in [("112_test_30", 4839702)]:
        pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out", target, "{}/{}".format(exp_num, job_num))
        print(pdb_dir)
        out_dirs = list(pdb_dir.glob("output_*"))
        log_files = [Path(out_dir, "log.csv") for out_dir in out_dirs]
        print(log_files)

        job_dir = pdb_dir.parents[0]
        exp_name = job_dir.stem
        N = 100
        n_cif = 1

        field = "xray_0"
        bonus_fields = ["pdb", "r_free_0"]
        # for i in range(1,n_cif):
        #     field += "+xray_{}".format(i)

        # if n_cif > 1:
        #     for i in range(n_cif):
        #         bonus_fields.append("xray_{}".format(i))

        # for i in range(n_cif):
        #     bonus_fields.append("r_free_{}".format(i))

        # for i in range(n_cif):
        #     bonus_fields.append("rmsd_avg_{}".format(i))

        # field = "r_free_0+r_free_1"
        # bonus_fields = ["r_free_0", "r_free_1"]

        data_dir = Path(Path.home(), "xray/sample_bench/data", target, exp_name)
        data_dir.mkdir(exist_ok=True)
        best_sample_file = Path(data_dir, "sample_{}.csv".format(field))

        stat_df = get_stat_df.get_stat_df(
            log_file_groups=[log_files],
            main_field=field,
            main_stat="min",
            N=N,
            bonus_fields=bonus_fields,
            equil=10,
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

        best_sample_df.to_csv(best_sample_file)
