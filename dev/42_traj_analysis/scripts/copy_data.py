from pathlib import Path
import pandas as pd
import shutil

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
from align_imp import compute_rmsd_between_average_pdb


if __name__ == "__main__":
    job_ids = list()
    out_ids = list()

    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/241_aniso_wxray")
    dest_dir = Path("../data/241")

    # for job_id in range(14):
    for job_dir in exp_dir.glob("*"):
        print(job_dir)
        job_id = job_dir.stem

        dest_pdb_dir = Path(dest_dir, "pdbs/{}".format(job_id))
        dest_pdb_dir.mkdir(exist_ok=True, parents=True)

        dest_log_dir = Path(dest_dir, "logs/{}".format(job_id))
        dest_log_dir.mkdir(exist_ok=True, parents=True)

        for out_id in range(10):
            out_dir = Path(exp_dir, "{}/output_{}".format(job_id, out_id))
            pdb_dir = Path(out_dir, "pdbs")

            pdb_files = list(pdb_dir.glob("*.pdb"))
            pdb_file_nums = set([int(p.stem) for p in pdb_files])

            if len(pdb_file_nums) == 0:
                continue

            max_pdb_file = Path(pdb_dir, "{}.pdb".format(max(pdb_file_nums)))
            # print(max_pdb_file)
            shutil.copy(max_pdb_file, Path(dest_pdb_dir, "{}.pdb".format(out_id)))

        # for out_id in range(10):
            out_dir = Path(exp_dir, "{}/output_{}".format(job_id, out_id))

            log_file = Path(out_dir, "log.csv")
            if log_file.exists():
                shutil.copy(log_file, Path(dest_log_dir, "{}.csv".format(out_id)))

            print(job_id, out_id)




