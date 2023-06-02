from pathlib import Path
import sys
import multiprocessing
import pandas as pd
import os
import numpy as np
import time

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_rmsd_df


if __name__ == "__main__":
    job_name = "01_wxray_sa"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhk", job_name)
    log_file = Path(Path.home(), "xray/sample_bench/data/7mhk", job_name, "rmsd.csv")
    n_jobs = 80

    log_df = pd.DataFrame(index=list(range(n_jobs)), columns=["avg_rmsd_mean", "all_traj_rmsd_mean"])

    for job_id in range(n_jobs+1):
        t0 = time.time()
        job_dir = Path(exp_dir, str(job_id))
        print(job_dir)

        out_dirs = list(job_dir.glob("output_*"))
        pdb_dir_groups = [[Path(out_dir, "pdbs")] for out_dir in out_dirs]
        all_traj_pdb_dir_groups = [[]]
        for pdb_dir_group in pdb_dir_groups:
            all_traj_pdb_dir_groups[0].extend(pdb_dir_group)

        equil, offset = 50, 10
        # Get average rmsd for each trajectory for a given w_xray.
        rmsd_df = get_rmsd_df.get_avg_rmsd_df(
            pdb_dir_groups=pdb_dir_groups,
            equil=equil,
            offset=offset
        )

        # Get the average rmsd across all trajectories for a given w_xray.
        all_traj_rmsd_df = get_rmsd_df.get_avg_rmsd_df(
            pdb_dir_groups=all_traj_pdb_dir_groups,
            equil=equil,
            offset=offset
        )

        avg_rmsd_mean = rmsd_df["rmsd_mean"].mean()
        all_traj_rmsd_mean = all_traj_rmsd_df["rmsd_mean"].iloc[0]

        log_df.loc[job_id, "avg_rmsd_mean"] = avg_rmsd_mean
        log_df.loc[job_id, "all_traj_rmsd_mean"] = all_traj_rmsd_mean

        log_df.to_csv(log_file)
        print(time.time() - t0)

        # break





