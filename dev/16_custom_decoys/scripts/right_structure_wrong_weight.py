from pathlib import Path
import random
random.seed(0)
import pandas as pd

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_refine_2.pdb")
    decoy_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data/3ca7/weight/3ca7_refine_2")
    decoy_meta_file = Path(Path.home(), "xray/decoys/data/3ca7/weight/3ca7_refine_2.csv")
    decoy_dir.mkdir(exist_ok=True)

    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    n_state = len(hs)

    cols = ["decoy"]
    for i in range(n_state):
        cols.extend(["structure_{}".format(i), "w_{}".format(i)])

    decoy_df = pd.DataFrame(index=list(range(1000)), columns=cols)

    for i in range(1000):
        decoy_file = Path(decoy_dir, "{}.pdb".format(i))
        decoy_df.loc[i, "decoy"] = decoy_file
        print(decoy_file)

        ws = [random.random() for j in range(n_state)]

        for j in range(n_state):
            w = ws[j]/sum(ws)

            decoy_df.loc[i, "structure_{}".format(j)] = pdb_file
            decoy_df.loc[i, "w_{}".format(j)] = w

            pids = IMP.atom.Selection(hs[j]).get_selected_particle_indexes()
            for pid in pids:
                IMP.atom.Atom(m, pid).set_occupancy(w)

        IMP.atom.write_multimodel_pdb(hs, str(decoy_file))

    decoy_df.to_csv(decoy_meta_file)
