from pathlib import Path
import random
import numpy as np
import math
import multiprocessing
import pandas as pd
import argparse
import sys

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import average_structure



def pool_compute_res_rmsf(
        params_dict
):
    pdb_files = params_dict["pdb_files"]
    ref_pdb_file = params_dict["ref_pdb_file"]
    res_id = params_dict["res_id"]

    ref_s = IMP.atom.AllPDBSelector()
    ref_m = IMP.Model()
    ref_h = IMP.atom.read_pdb(str(ref_pdb_file), ref_m, ref_s)

    hs = list()
    ms = list()
    for pdb_file in pdb_files:
        s = IMP.atom.AllPDBSelector()
        m = IMP.Model()
        h = IMP.atom.read_pdb(str(pdb_file), m, s)

        hs.append(h)
        ms.append(m)

    res_rmsf = compute_res_rmsf(
        hs=hs,
        ref_h=ref_h,
        res_id=res_id
    )

    print(res_id, res_rmsf)

    return res_id, res_rmsf


def compute_res_rmsf(
        hs,
        ref_h,
        res_id
):
    ref_pids = IMP.atom.Selection(ref_h, residue_index=res_id, atom_type=IMP.atom.AT_CA).get_selected_particle_indexes()
    ref_xyzs = [IMP.core.XYZ(ref_h.get_model(), pid) for pid in ref_pids]

    squares = list()
    for h in hs:
        pids = IMP.atom.Selection(h, residue_index=res_id, atom_type=IMP.atom.AT_CA).get_selected_particle_indexes()
        xyzs = [IMP.core.XYZ(h.get_model(), pid) for pid in pids]

        for i in range(len(xyzs)):
            ref_xyz = ref_xyzs[i]
            xyz = xyzs[i]

            ref_coord = ref_xyz.get_coordinates()
            coord = xyz.get_coordinates()

            square = 0
            for i in range(3):
                square = square + (coord[i] - ref_coord[i])**2

            squares.append(square)

    res_rmsf = math.sqrt(np.mean(squares))
    return res_rmsf


def get_res_rmsf_df(
        pdb_files,
        ref_pdb_file,
        start_res_id,
        stop_res_id
):
    pool_params = list()
    for res_id in range(start_res_id,stop_res_id+1):
        params_dict = dict()
        params_dict["pdb_files"] = pdb_files
        params_dict["ref_pdb_file"] = ref_pdb_file
        params_dict["res_id"] = res_id
        pool_params.append(params_dict)

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_compute_res_rmsf,
        pool_params
    )

    rmsf_df = pd.DataFrame(index=list(range(48,98)), columns=["rmsf"])

    for pool_result in pool_results:
        res_id, res_rmsf = pool_result
        rmsf_df.loc[res_id, "rmsf"] = res_rmsf

    return rmsf_df


def get_avg_rmsf_df(
        pdb_dir_groups,
        equil
):
    indices = [str(group) for group in pdb_dir_groups]
    avg_rmsf_df = pd.DataFrame(index=indices, columns=["rmsf_avg"])

    # Need to use a different multiprocessing strategy because of how expensive the calculations are. get_stat_df computes the value for each group with a single process, while here we want to compute the value for each group with many processes.

    # It may be necesary to only compute these for a subset of pdb files.
    for pdb_dir_group in pdb_dir_groups:
        pdb_files = list()
        for pdb_dir in pdb_dir_group:
            pdb_files.extend(list(pdb_dir.glob("*.pdb")))

        avg_pdb_file = average_structure.get_average_pdb_file_from_pdb_files(pdb_files=pdb_files)
        print(avg_pdb_file)

        res_rmsf_df = get_res_rmsf_df(
            pdb_files=pdb_files,
            ref_pdb_file=avg_pdb_file,
            start_res_id=48,
            stop_res_id=97
        )

        rmsf_mean = res_rmsf_df["rmsf"].mean()
        avg_rmsf_df.loc[str(pdb_dir_group), "rmsf_avg"] = rmsf_mean
        print(avg_rmsf_df.head())

    return avg_rmsf_df


if __name__ == "__main__":
    pdb_dir_groups=[[Path(Path.home(), "xray/sample_bench/unit_tests/data")]]
    get_avg_rmsf_df(
        pdb_dir_groups=pdb_dir_groups,
        equil=0
    )


    # # parser = argparse.ArgumentParser()
    # # parser.add_argument("--job_name")
    # # parser.add_argument("--job_num", type=int)
    # # parser.add_argument("--n", type=int)
    # # args = parser.parse_args()
    # # print(args.job_name)
    # # print(args.job_num)
    # # print(args.n)
    # #
    # # job_name = args.job_name
    # # job_num = args.job_num
    # # n = args.n

    # # job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/single_md", job_name, str(job_num))
    # rmsf_file = Path(Path.home(), "xray/sample_bench/data/analysis/24_300_exp_equil_com/best_score/rmsf.csv")

    # # From each output dir, randomly select 10 pdb files for rmsf calculations.
    # # pdb_files = get_n_files_from_all_output_dirs(
    # #     job_dir=job_dir,
    # #     n=n,
    # #     equil=100
    # # )
    # pdb_dir = Path(Path.home(), "xray/sample_bench/data/analysis/24_300_exp_equil_com/best_score/pdbs")
    # pdb_files = list(pdb_dir.glob("*.pdb"))

    # # We run 3 analysis. First, rmsf of all structures to the native structure. Second, rmsf of all structures to the average structure. Third, average rmsf of all structures from each run to the average structure of the run.
    # native_pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    # avg_pdb_file = sample_average.get_average_pdb_file_from_pdb_files(
    #     pdb_files=pdb_files
    # )

    # all_res_rmsf_df = pd.DataFrame(index=list(range(48,97+1)))
    # for ref_pdb_file, name in [(avg_pdb_file, "avg")]:
    #     # Compute the rmsf for all the randomly selected pdb files relative to the reference pdb file (native).
    #     rmsf_df = get_res_rmsf_df(
    #         pdb_files=pdb_files,
    #         ref_pdb_file=ref_pdb_file,
    #         start_res_id=48,
    #         stop_res_id=97
    #     )
    #     # res_rmsf_df = res_rmsf_df.rename(columns={"rmsf": name})
    #     print(rmsf_df.head())
    #     all_res_rmsf_df = all_res_rmsf_df.join(rmsf_df)

    # # all_res_rmsf_df = all_res_rmsf_df.reset_index()
    # all_res_rmsf_df.to_csv(rmsf_file)