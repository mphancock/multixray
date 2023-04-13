from pathlib import Path

import IMP
import IMP.atom


def set_b_factors():
    for pid in m.get_particle_indexes():
        at = IMP.atom.Atom(m, pid)
        at.set_temperature_factor(20)


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
    pdb_file = Path(Path.home(), "xray/tmp/9.pdb")
    new_pdb_file = Path(Path.home(), "xray/tmp/9_25.pdb")

    m = IMP.Model()
    sel = IMP.atom.AllPDBSelector()

    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, sel)

    set_occupancies(
        hs=hs,
        occ=.25
    )

    IMP.atom.write_multimodel_pdb(hs, str(new_pdb_file))
