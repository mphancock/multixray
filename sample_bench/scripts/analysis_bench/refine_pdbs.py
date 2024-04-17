from pathlib import Path
import sys
import pandas as pd
import multiprocessing

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import refine, pool_refine


if __name__ == "__main__":
    pdb_df = pd.read_csv(Path(Path.home(), "xray/sample_bench/data/7mhf/166_N1/best_r_free.csv"), index_col=0)

    job_dir = Path(Path.home(), "xray/sample_bench/data/7mhf/166_N1")
    new_pdb_name = "best_r_free_ref"
    pdb_dir = Path(job_dir, new_pdb_name)

    ref_pdb_df = pd.DataFrame(index=list(range(24)),columns=pdb_df.columns)

    pool_params = list()
    n_step = 10

    all_hs = list()
    ms = list()
    for i in range(len(pdb_df)):
        pdb_file = pdb_df.loc[i, "pdb"]
        print(pdb_file)
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        all_hs.append(hs)
        ms.append(m)

        params_dict = dict()
        params_dict["hs"] = hs
        params_dict["n_step"] = n_step
        params_dict["log_file"] = None
        pool_params.append(params_dict)

        # print(IMP.core.XYZ(m, IMP.atom.Selection(hs[0]).get_selected_particle_indexes()[0]).get_x())
        # break

    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    pool_results = pool_obj.imap(pool_refine, pool_params)

    all_ref_hs = list()
    ref_ms = list()
    for pool_result in pool_results:
        m, hs = pool_result
        all_ref_hs.append(hs)
        ref_ms.append(m)

        # print(IMP.core.XYZ(m, IMP.atom.Selection(hs[0]).get_selected_particle_indexes()[0]).get_x())


    #     # refine(hs=hs, n_step=n_step)

    for i in range(len(all_ref_hs)):
        hs = all_ref_hs[i]
        # print(IMP.core.XYZ(hs[0].get_model(), IMP.atom.Selection(hs[0]).get_selected_particle_indexes()[0]).get_x())

        N = pdb_df.loc[i, "N"]
        cif = pdb_df.loc[i, "cif"]
        ref_pdb_df.loc[i, "N"] = N
        ref_pdb_df.loc[i, "cif"] = cif

        for j in range(N):
            ref_pdb_df.loc[i, "w_{}".format(j)] = pdb_df.loc[i, "w_{}".format(j)]

        ref_pdb_file = Path(pdb_dir, "{}.pdb".format(i))
        ref_pdb_df.loc[i, "pdb"] = ref_pdb_file

        print(ref_pdb_file)
        IMP.atom.write_multimodel_pdb(hs, str(ref_pdb_file))

    ref_pdb_df.to_csv(Path(job_dir, "{}.csv".format(new_pdb_name)))



