from pathlib import Path

import IMP
import IMP.core
import IMP.atom
import IMP.algebra


if __name__ == "__main__":
    pdb_dir = Path(Path.home(), "xray/dev/21_degeneracy/data/degens")

    ref_ca_cb_dist = 1.51

    for pdb_file in pdb_dir.glob("*.pdb"):
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        ca_cb_dist_delta = 0

        for h in hs:
            for res in [71, 72]:
                ca = IMP.atom.Selection(h, residue_index=res, atom_type=IMP.atom.AT_CA).get_selected_particles()[0]
                cb = IMP.atom.Selection(h, residue_index=res, atom_type=IMP.atom.AT_CB).get_selected_particles()[0]

                ca_xyz = IMP.core.XYZ(m, ca).get_coordinates()
                cb_xyz = IMP.core.XYZ(m, cb).get_coordinates()

                ca_cb_dist = (ca_xyz - cb_xyz).get_magnitude()

                ca_cb_dist_delta  =ca_cb_dist_delta + abs(ca_cb_dist - ref_ca_cb_dist)

        print(pdb_file.stem, ca_cb_dist_delta)






