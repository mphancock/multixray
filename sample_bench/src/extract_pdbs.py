from pathlib import Path
import random
import shutil
import pandas as pd


def pool_get_n_random_files_from_pdb_dirs(
        params_dict
):
    job_id = params_dict["job_id"]
    sample_files = get_n_random_files_from_pdb_dirs(
        pdb_dirs=params_dict["pdb_dir"],
        n=params_dict["n"],
        equil=params_dict["equil"]
    )

    return job_id, sample_files


def get_n_random_files_from_pdb_dirs(
        pdb_dirs,
        n,
        equil
):
    all_pdb_files = list()
    for pdb_dir in pdb_dirs:
        # print(pdb_dir)
        pdb_files = list(pdb_dir.glob("*.pdb"))

        pdb_file_nums = [int(pdb_file.stem) for pdb_file in pdb_files]
        max_pdb_file_num = max(pdb_file_nums)

        # equil_pdb_files = list()
        for i in range(equil, max_pdb_file_num+1):
            all_pdb_files.append(Path(pdb_dir, "{}.pdb".format(i)))

    # print(len(all_pdb_files))
    sample_pdb_files = random.sample(all_pdb_files, n)
    return sample_pdb_files


def get_n_best_scoring_structures_from_out_dirs(
          out_dirs,
          N
):
    pdb_log_dfs = list()
    for out_dir in out_dirs:
        print(out_dir)
        out_id = out_dir.stem.split("_")[1]

        log_file = Path(out_dir, "log.csv")
        pdb_dir = Path(out_dir, "pdbs")

        pdb_files = list(pdb_dir.glob("*.pdb"))
        pdb_ids = [pdb_file.stem for pdb_file in pdb_files]

        log_df = pd.read_csv(log_file, index_col=0)
        pdb_log_df = log_df.iloc[::10]
        pdb_log_df = pdb_log_df.iloc[:len(pdb_files)]
        # There may be more pdbs than log entries.
        pdb_log_df["tot"] = pdb_log_df["xray"]*30000+pdb_log_df["ff"]
        pdb_log_df["pdb_file"] = [Path(pdb_dir, "{}.pdb".format(pdb_id)) for pdb_id in range(len(pdb_files))]
        pdb_log_dfs.append(pdb_log_df)

    all_pdb_log_df = pd.concat(pdb_log_dfs)
    best_df = all_pdb_log_df.nsmallest(N, "tot")
    print(best_df)
    best_pdf_files = list(best_df["pdb_file"])

    return best_pdf_files


if __name__ == "__main__":
    job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/34_1/2406616")
    pdb_files = get_n_best_scoring_structures_from_out_dirs(
        out_dirs=list(job_dir.glob("output_*"))[:10],
        N=10
    )
    print(pdb_files)


