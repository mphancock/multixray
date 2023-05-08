from pathlib import Path
import numpy as np
import multiprocessing
import pandas as pd
import sys

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp
import average_structure


def get_avg_rmsd_df(
        pdb_dir_groups,
        equil,
        offset
):
    indices = [str(group) for group in pdb_dir_groups]
    avg_rmsd_df = pd.DataFrame(index=indices, columns=["rmsd_avg"])

    # Need to use a different multiprocessing strategy because of how expensive the calculations are. get_stat_df computes the value for each group with a single process, while here we want to compute the value for each group with many processes.

    # It may be necesary to only compute these for a subset of pdb files.
    for pdb_dir_group in pdb_dir_groups:
        pdb_files = list()
        for pdb_dir in pdb_dir_group:
            pdb_dir_files = list(pdb_dir.glob("*.pdb"))
            pdb_dir_files_sort = list()
            for i in range(len(pdb_dir_files)):
                pdb_file = Path(pdb_dir, "{}.pdb".format(i))
                pdb_dir_files_sort.append(pdb_file)

            # Only include pdb files post equilibration with some offset.
            pdb_files.extend(pdb_dir_files_sort[equil::offset])

        if len(pdb_files) == 0:
            avg_rmsd = np.nan
        else:
            # The slower step is computing the average pdb file.
            avg_pdb_file = average_structure.get_average_pdb_file_from_pdb_files(pdb_files=pdb_files)
            print(avg_pdb_file)

            pool_params = list()
            for pdb_file in pdb_files:
                params_dict = dict()
                params_dict["pdb_file_0"] = pdb_file
                params_dict["pdb_file_1"] = avg_pdb_file
                pool_params.append(params_dict)
            print(len(pool_params))

            pool_obj = multiprocessing.Pool(
                multiprocessing.cpu_count()
            )

            pool_results = pool_obj.imap(
                align_imp.pool_compute_rmsd_aligning_first_to_second,
                pool_params
            )

            avg_rmsd = 0
            n_rmsd = 0
            for rmsd in pool_results:
                # print(rmsd)
                avg_rmsd = avg_rmsd + rmsd
                n_rmsd = n_rmsd + 1

            avg_rmsd = avg_rmsd / n_rmsd
            avg_rmsd_df.loc[str(pdb_dir_group), "rmsd_mean"] = avg_rmsd

    return avg_rmsd_df


if __name__ == "__main__":
    # pdb_dir_groups=[[Path(Path.home(), "xray/sample_bench/unit_tests/data")]]
    pdb_dir_groups=[list(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/46_w_xray/0").glob("*/pdbs"))]
    avg_rmsd_df = get_avg_rmsd_df(
        pdb_dir_groups=pdb_dir_groups,
        equil=100,
        offset=50
    )

    # print(avg_rmsd_df.head())