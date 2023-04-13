from pathlib import Path
import sys

import IMP
import IMP.atom
import IMP.core
import IMP.algebra

sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


if __name__ == "__main__":
    decoy_dir = Path(Path.home(), "xray/decoys/data/decoy_sets/3ca7/3ca7_N_1000_x1")

    ref_pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)

    pids_0 = IMP.atom.Selection(h_0).get_selected_particle_indexes()

    for pdb_file in decoy_dir.glob("*"):
        m = IMP.Model()
        s = IMP.atom.AllPDBSelector()
        h = IMP.atom.read_pdb(str(pdb_file), m, s)

        rmsd = align_imp.compute_rmsd(
            h=h,
            h_0=h_0
        )

        pids = IMP.atom.Selection(h).get_selected_particle_indexes()

        dists = list()
        for i in range(len(pids_0)):
            pid_0 = pids_0[i]
            pid = pids[i]

            xyz_0 = IMP.core.XYZ(m_0, pid_0)
            xyz = IMP.core.XYZ(m, pid)

            dist = IMP.algebra.get_distance(xyz_0.get_coordinates(), xyz.get_coordinates())

            dists.append(dist)

        print(pdb_file.stem, rmsd, max(dists))
