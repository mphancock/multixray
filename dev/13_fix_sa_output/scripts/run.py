from pathlib import Path
import pandas as pd
import shutil
import multiprocessing
import argparse


def pool_fix_dir(
    params_dict
):
    out_dir = params_dict["out_dir"]
    copy_out_dir = params_dict["copy_out_dir"]

    print(out_dir, copy_out_dir)

    log_df = pd.read_csv(Path(out_dir, "log.csv"), index_col=0)

    drop_ids = list()
    for i in range(len(log_df)):
        copy = True
        if i > 0:
            if i % 1500 == 0 or i % 1500 > 500:
                drop_ids.append(i)
                copy = False

        if copy and i % 10 == 0:
            pdb_id = i // 10
            pdb_file = Path(out_dir, "pdbs", "{}.pdb".format(pdb_id))
            copy_pdb_file = Path(copy_out_dir, "pdbs", "{}.pdb".format(pdb_id))
            if pdb_file.exists():
                shutil.copy(pdb_file, copy_pdb_file)

    log_df.drop(drop_ids, inplace=True)
    log_df.to_csv(Path(copy_out_dir, "log.csv"))

    return copy_out_dir


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--job_dir")
    parser.add_argument("--copy_dir")
    args = parser.parse_args()
    print(args.job_dir)
    print(args.copy_dir)

    job_dir = Path(args.job_dir)
    copy_dir = Path(args.copy_dir)
    shutil.rmtree(copy_dir, ignore_errors=True)
    copy_dir.mkdir()

    pool_params = list()
    for out_dir in job_dir.glob("output*"):
        copy_out_dir = Path(copy_dir, out_dir.name)
        copy_out_dir.mkdir()

        pdb_copy_dir = Path(copy_out_dir, "pdbs")
        pdb_copy_dir.mkdir()

        params_dict = dict()
        params_dict["out_dir"] = out_dir
        params_dict["copy_out_dir"] = copy_out_dir
        pool_params.append(params_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_fix_dir,
        pool_params
    )

    for result in pool_results:
        continue







