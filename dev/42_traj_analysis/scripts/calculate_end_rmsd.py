from pathlib import Path
import pandas as pd

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
from align_imp import compute_rmsd_between_average_pdb


if __name__ == "__main__":
    job_ids = list()
    out_ids = list()
    end_rmsds = list()

    data_dir = Path("../data/209")
    data_dir.mkdir(exist_ok=True)

    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/209_9000_low_res")

    for job_id in range(26):
        for out_id in range(10):
            out_path = Path(exp_dir, "{}/output_{}".format(job_id, out_id))

            log_file = Path(out_path, "log.csv")
            log_df = pd.read_csv(log_file, index_col=0)

            if False:
                rmsd = log_df.loc[len(log_df)-1, "rmsd_0"]
            else:
                rmsd = log_df.loc[len(log_df)-1, "rmsd_7mhl"]

            job_ids.append(job_id)
            out_ids.append(out_id)
            end_rmsds.append(rmsd)

            print(job_id, out_id, rmsd)

    end_rmsd_df = pd.DataFrame()
    end_rmsd_df["job_id"] = job_ids
    end_rmsd_df["out_id"] = out_ids
    end_rmsd_df["end_rmsd"] = end_rmsds
    end_rmsd_df.to_csv(Path(data_dir, "end_rmsd.csv"))



