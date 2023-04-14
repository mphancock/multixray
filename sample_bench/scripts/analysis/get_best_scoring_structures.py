from pathlib import Path
import pandas as pd
import numpy as np
import shutil


if __name__ == "__main__":
    job_name = "24_300_exp_equil_com"
    job_id = "2288531"

    job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/single_md/{}/{}".format(job_name, job_id))
    out_dirs = list(job_dir.glob("output_*"))
    # pdb_stat_df = pd.DataFrame(columns=["out_id", "pdb_file", "tot_score", "rmsd"])
    for out_dir in out_dirs:
        print(out_dir)
        out_id = out_dir.stem.split("_")[1]

        log_file = Path(out_dir, "log.csv")
        pdb_dir = Path(out_dir, "pdbs")

        pdb_files = list(pdb_dir.glob("*.pdb"))
        pdb_ids = [pdb_file.stem for pdb_file in pdb_files]
        pdb_ids.remove("-1")

        log_df = pd.read_csv(log_file, index_col=0)
        log_df["tot"] = log_df["xray"]*30000+log_df["ff"]

        min_tot_score = np.infty
        min_tot_score_pdb_id = np.infty
        for pdb_id in pdb_ids:
            row_id = int(pdb_id)*10+10
            tot_score = log_df.loc[row_id, "tot"]
            if tot_score < min_tot_score:
                min_tot_score = tot_score
                min_tot_score_pdb_id = pdb_id

        print(min_tot_score_pdb_id, min_tot_score)
        orig_pdb_file = Path(pdb_dir, "{}.pdb".format(min_tot_score_pdb_id))

        new_pdb_file = Path(Path.home(), "xray/sample_bench/data/analysis/{}/best_score/pdbs/{}.pdb".format(job_name, out_id))
        shutil.copy(orig_pdb_file, new_pdb_file)
