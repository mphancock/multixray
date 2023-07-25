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
    name = "wrong_structure_right_weight"

    pdb_file_0 = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_refine_2.pdb")
    m_0 = IMP.Model()
    hs_0 = IMP.atom.read_multimodel_pdb(str(pdb_file_0), m_0, IMP.atom.AllPDBSelector())
    n_state = len(hs_0)
    occs_0 = list()
    for i in range(n_state):
        pids = IMP.atom.Selection(hs_0[i]).get_selected_particle_indexes()
        occs_0.append(IMP.atom.Atom(m_0, pids[0]).get_occupancy())

    orig_decoy_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data/3ca7/53_100/rand_1000_2x_53_54")

    decoy_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data/3ca7/weight", name)
    decoy_meta_file = Path(Path.home(), "xray/decoys/data/3ca7/weight/{}.csv".format(name))
    decoy_dir.mkdir(exist_ok=True)

    cols = ["decoy"]
    for i in range(n_state):
        cols.extend(["structure_{}".format(i), "state_{}".format(i), "w_{}".format(i)])

    decoy_df = pd.DataFrame(index=list(range(1000)), columns=cols)

    orig_decoy_files = list(orig_decoy_dir.glob("*.pdb"))
    for i in range(1000):
        pdb_file = orig_decoy_files[i]

        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        flip_rmsd = False
        rmsd = align_imp.compute_rmsd(
            h_0s=hs,
            h_1s=hs_0,
            ca_only=True
        )

        rmsd_tmp = align_imp.compute_rmsd(
            h_0s=[hs[1],hs[0]],
            h_1s=hs_0,
            ca_only=True
        )
        print(rmsd, rmsd_tmp)

        if rmsd_tmp < rmsd:
            flip_rmsd = True
            rmsd = rmsd_tmp

        decoy_df.loc[i, "structure_0"] = pdb_file
        decoy_df.loc[i, "structure_1"] = pdb_file

        decoy_df.loc[i, "w_0"] = occs_0[0]
        decoy_df.loc[i, "w_1"] = occs_0[1]

        if flip_rmsd:
            decoy_df.loc[i, "state_0"] = 1
            decoy_df.loc[i, "state_1"] = 0
            hs = [hs[1],hs[0]]
        else:
            decoy_df.loc[i, "state_0"] = 0
            decoy_df.loc[i, "state_1"] = 1

        for j in range(2):
            pids = IMP.atom.Selection(hs[j]).get_selected_particle_indexes()

            for pid in pids:
                IMP.atom.Atom(m, pid).set_occupancy(decoy_df.loc[i, "w_{}".format(j)])

        decoy_file = Path(decoy_dir, "{}.pdb".format(i))
        decoy_df.loc[i, "decoy"] = decoy_file
        print(decoy_file)

        IMP.atom.write_multimodel_pdb(hs, str(decoy_file))

    decoy_df.to_csv(decoy_meta_file)












