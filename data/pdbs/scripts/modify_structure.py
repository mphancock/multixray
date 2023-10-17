from pathlib import Path

import IMP
import IMP.atom


def set_b_factors(
    hs,
    b
):
    for h in hs:
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        for pid in pids:
            at = IMP.atom.Atom(m, pid)
            at.set_temperature_factor(b)


def set_occupancies(
        hs,
        occ
):
    m = hs[0].get_model()
    for h in hs:
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(m, pid).set_occupancy(occ)


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_refine.pdb")
    new_pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_refine_b_factor.pdb")

    m = IMP.Model()
    sel = IMP.atom.AllPDBSelector()

    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, sel)

    set_b_factors(hs, 15)

    IMP.atom.write_multimodel_pdb(hs, str(new_pdb_file))
