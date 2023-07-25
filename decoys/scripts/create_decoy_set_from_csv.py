from pathlib import Path
import pandas as pd
import shutil

import IMP
import IMP.atom


def set_w(
        h,
        w
):
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    for pid in pids:
        IMP.atom.Atom(h.get_model(), pid).set_occupancy(w)

if __name__ == "__main__":
    # decoy_meta_file = Path(Path.home(), "xray/decoys/data/7mhf/38_7mhf_decoys/merge_2000.csv")
    decoy_meta_file = Path(Path.home(), "xray/dev/17_synthetic_native/data/1x_decoy_df.csv")

    decoy_df = pd.read_csv(decoy_meta_file, index_col=0)

    pdb_dir = Path(decoy_df.loc[0, "decoy_file"]).parents[0]
    pdb_dir.mkdir(exist_ok=True)

    n_state = (len(decoy_df.columns)-1)//3
    print(n_state)
    for i in range(len(decoy_df)):
        decoy_file = Path(decoy_df.loc[i, "decoy_file"])

        m_decoys = list()
        h_decoys = list()
        for j in range(n_state):
            pdb_file = Path(decoy_df.loc[i, "pdb_{}".format(j)])
            state = decoy_df.loc[i, "id_{}".format(j)]
            w = decoy_df.loc[i, "w_{}".format(j)]

            m = IMP.Model()
            hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
            h = hs[state]
            set_w(h, w)

            m_decoys.append(m)
            h_decoys.append(h)

        IMP.atom.write_multimodel_pdb(h_decoys,str(decoy_file))





