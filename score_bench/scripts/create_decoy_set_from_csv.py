from pathlib import Path
import pandas as pd
import shutil
import multiprocessing

import IMP
import IMP.atom


def write_multistate_pool(
        params_dict
):
    decoy_file = params_dict["decoy_file"]
    pdb_files = params_dict["pdb_files"]
    ids = params_dict["ids"]
    weights = params_dict["weights"]

    n_state = len(pdb_files)

    m_decoys = list()
    h_decoys = list()
    for i in range(n_state):
        pdb_file = pdb_files[i]
        id = ids[i]
        w = weights[i]

        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        h = hs[id]
        set_w(h, w)

        m_decoys.append(m)
        h_decoys.append(h)

    IMP.atom.write_multimodel_pdb(h_decoys,str(decoy_file))

    print(decoy_file)
    return decoy_file


def set_w(
        h,
        w
    ):
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(h.get_model(), pid).set_occupancy(w)

if __name__ == "__main__":
    for i in range(10):
        decoy_meta_file = Path(Path.home(), "xray/decoys/data/3ca7/100_natives_4x/{}.csv".format(i))

        decoy_df = pd.read_csv(decoy_meta_file, index_col=0)

        pdb_dir = Path(decoy_df.loc[0, "decoy_file"]).parents[0]
        shutil.rmtree(pdb_dir)
        pdb_dir.mkdir(exist_ok=True)

        n_state = (len(decoy_df.columns)-1)//3
        print(n_state)

        pool_params = list()
        for i in range(len(decoy_df)):
            params_dict = dict()
            decoy_file = Path(decoy_df.loc[i, "decoy_file"])

            params_dict["decoy_file"] = decoy_file
            params_dict["pdb_files"] = list()
            params_dict["ids"] = list()
            params_dict["weights"] = list()

            m_decoys = list()
            h_decoys = list()
            for j in range(n_state):
                pdb_file = Path(decoy_df.loc[i, "pdb_{}".format(j)])
                id = int(decoy_df.loc[i, "id_{}".format(j)])
                w = decoy_df.loc[i, "w_{}".format(j)]

                params_dict["pdb_files"].append(pdb_file)
                params_dict["ids"].append(id)
                params_dict["weights"].append(w)

            pool_params.append(params_dict)

        pool_obj = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        pool_results = pool_obj.imap(
            write_multistate_pool,
            pool_params
        )

        for pool_result in pool_results:
            continue

        pool_obj.close()







