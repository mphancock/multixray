from pathlib import Path
import sys
import pandas as pd
import random
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/decoys/scripts")))
import merge_pdbs


def pool_make_multi_state_pdb_files(
        param_dict
):
    return make_multi_state_pdb_file(
        id=param_dict["id"],
        pdb_files=param_dict["pdb_files"],
        out_file=param_dict["out_file"],
        n=param_dict["n"]
    )


def make_multi_state_pdb_file(
        id,
        pdb_files,
        out_file,
        n
):
    print(id)

    sample_pdb_files = random.sample(pdb_files, n)
    merge_pdbs.merge_pdbs(
        pdb_files=sample_pdb_files,
        weights=[1]*n,
        out_file=out_file,
        n=-1,
        order=False
    )

    return id, sample_pdb_files


if __name__ == "__main__":
    single_pdb_dir = Path(Path.home(), "xray/sample_bench/analysis/24_300_exp_equil_com/best_score/pdbs")
    pdb_files = list(single_pdb_dir.glob("*.pdb"))
    n = 32

    multi_pdb_dir = Path(Path.home(), "xray/sample_bench/analysis/24_300_exp_equil_com/multi_state_best_score/pdbs_{}".format(n))
    multi_pdb_dir.mkdir()

    meta_df = pd.DataFrame(index=list(range(1000)), columns=list(range(n)))

    params = list()
    for i in range(1000):
        param_dict = dict()
        param_dict["id"] = i
        param_dict["pdb_files"] = pdb_files
        param_dict["out_file"] = Path(multi_pdb_dir, "{}.pdb".format(i))
        param_dict["n"] = n

        params.append(param_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_make_multi_state_pdb_files,
        params
    )

    for id, sample_pdb_files in pool_results:
        meta_df.loc[id, list(range(n))] = sample_pdb_files

    meta_df.to_csv(Path(Path.home(), "xray/sample_bench/analysis/24_300_exp_equil_com/multi_state_best_score/meta/pdbs_{}.csv".format(n)))