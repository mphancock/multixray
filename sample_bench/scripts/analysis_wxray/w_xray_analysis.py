from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
from utility import pool_read_pdb
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis_exp")))
# from score_analysis_exp import field_id_mapper, map_std_field_to_field


if __name__ == "__main__":
    exp_name = "180_wxray_bench"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    analysis_dir = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf", exp_name)
    analysis_dir.mkdir(exist_ok=True)
    n_states = [1, 2]
    # job_ids = list(range(63))
    job_ids = list(range(20))
    w_xrays = [0.03125, 0.0625, .125, .25, .5, 1, 2, 4, 8, 16, 32]
    # cif_df = pd.read_csv(Path("/wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv"), index_col=0)
    cif_df = pd.read_csv(Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/csvs/7mhf_20.csv"), index_col=0)

    w_xray_df = pd.DataFrame()

    for state_id in range(len(n_states)):
        for job_id in job_ids:
            cif_files = [Path(cif_file) for cif_file in cif_df.loc[job_id, "cifs"].split(",")]
            cif_names = [Path(cif_file).stem for cif_file in cif_files]
            field = "r_free_{}".format(cif_names[0])
            bonus_fields = ["r_free_{}".format(cif_names[0])]
            for i in range(1, len(cif_names)):
                # field = field+"+r_free_{}".format(cif_names[i])
                bonus_fields.append("r_free_{}".format(cif_names[i]))

            bonus_fields.extend(["ff", "pdb"])
            if field in bonus_fields:
                bonus_fields.remove(field)

            for w_xray_id in range(len(w_xrays)):
                job_dir = Path(exp_dir, "N{}_J{}_W{}".format(state_id, job_id, w_xray_id))
                log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("*/")]

                if len(log_files) == 0:
                    continue

                print(state_id, job_id, w_xray_id, len(log_files))
                print(field, bonus_fields)

                stat_df = get_stat_df(
                    log_files=log_files,
                    field=field,
                    N=1,
                    bonus_fields=bonus_fields,
                    equil=100,
                    pdb_only=True,
                    max_rmsd=None
                )

                if len(stat_df) == 0:
                    continue

                df_id = len(w_xray_df)
                w_xray_df.loc[df_id, "state_id"] = state_id
                w_xray_df.loc[df_id, "job_id"] = job_id
                w_xray_df.loc[df_id, "w_xray_id"] = w_xray_id
                w_xray_df.loc[df_id, "field"] = stat_df.loc[0, field]

                if field in ["r_free_{}".format(cif_name) for cif_name in cif_names]:
                    w_xray_df.loc[df_id, field] = stat_df.loc[0, field]

                for bonus_field in bonus_fields:
                    w_xray_df.loc[df_id, bonus_field] = stat_df.loc[0, bonus_field]

    w_xray_df.to_csv(Path(analysis_dir, "w_xray.csv"))
