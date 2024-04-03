from pathlib import Path
import pandas as pd
import shutil
import multiprocessing
import argparse


def copy_output_dir(out_dir):
    print(out_dir)
    log_file = Path(out_dir, "log.csv")
    file_parts = list(log_file.parts)
    output_name = file_parts[-2]

    new_out_dir = Path(new_dir, output_name)
    new_out_dir.mkdir(exist_ok=True, parents=True)

    new_pdb_dir = Path(new_out_dir, "pdbs")
    new_pdb_dir.mkdir(exist_ok=True)

    shutil.copy(Path(out_dir, "params.txt"), Path(new_out_dir, "params.txt"))

    try:
        log_df = pd.read_csv(log_file, index_col=0)
    except pd.errors.EmptyDataError:
        print("FAILED: {}".format(out_dir))
        return

    pdb_indices = log_df[~log_df['pdb'].isna()].index

    for id in pdb_indices:
        pdb_file = Path(log_df.loc[id, "pdb"])
        # print(pdb_file)
        new_pdb_file = Path(new_pdb_dir, pdb_file.name)
        log_df.loc[id, "pdb"] = str(new_pdb_file)

        shutil.copy(pdb_file, new_pdb_file)

    log_df.to_csv(Path(new_out_dir, "log.csv"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--orig")
    parser.add_argument("--new")
    args = parser.parse_args()

    data_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf")
    orig_dir = Path(data_dir, args.orig)
    new_dir = Path(data_dir, args.new)

    print(orig_dir)

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())

    out_dirs = list(orig_dir.glob("output_*"))
    # out_dirs = [Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/138_4_state_2_cond/131217/output_713")]

    # for out_dir in out_dirs:
    #     copy_output_dir(out_dir)

    pool_results = pool_obj.imap(copy_output_dir, out_dirs)
    for result in pool_results:
        continue


