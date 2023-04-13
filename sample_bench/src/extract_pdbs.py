from pathlib import Path
import random
import shutil


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



