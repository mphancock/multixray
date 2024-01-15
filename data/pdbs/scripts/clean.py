from pathlib import Path
import sys

import IMP
import IMP.atom


def normalize_coordinates(
        m, pids
):
    x_coords = list()
    y_coords = list()
    z_coords = list()

    for pid in pids:
        xyz = IMP.core.XYZ(m, pid)
        x_coords.append(xyz.get_x())
        y_coords.append(xyz.get_y())
        z_coords.append(xyz.get_z())

    x_max, x_min = max(x_coords), min(x_coords)
    y_max, y_min = max(y_coords), min(y_coords)
    z_max, z_min = max(z_coords), min(z_coords)

    for pid in pids:
        xyz = IMP.core.XYZ(m, pid)
        x_coord = xyz.get_x()
        y_coord = xyz.get_y()
        z_coord = xyz.get_z()

        xyz.set_x(x_coord - x_min)
        xyz.set_y(y_coord - y_min)
        xyz.set_z(z_coord - z_min)

    print(x_max - x_min)
    print(y_max - y_min)
    print(z_max - z_min)


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhi.pdb")
    save_pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhi_clean.pdb")

    m = IMP.Model()
    # sel = IMP.atom.NonWaterPDBSelector()
    # sel_1 = IMP.atom.NonHydrogenPDBSelector()
    # sel = IMP.atom.AllPDBSelector()
    sel_2 = IMP.atom.NonAlternativePDBSelector()
    sel_1 = IMP.atom.NonWaterNonHydrogenPDBSelector()
    # sel_2 = IMP.atom.WaterPDBSelector()
    sel = IMP.atom.AndPDBSelector(sel_1, sel_2)

    h = IMP.atom.read_pdb(
        str(pdb_file),
        m,
        sel
    )

    # Add any missing atoms via CHARMM force field.
    # charmm_restraints(m, h)
    atoms = IMP.atom.Selection(hierarchy=h)
    pids = atoms.get_selected_particle_indexes()
    print(len(pids))

    # Set all atomic b factors to 0 and occupancies to 1.
    # for pid in pids:
    #     atom = IMP.atom.Atom(m, pid)
    #     atom.set_occupancy(1)
        # atom.set_temperature_factor(0)

    IMP.atom.write_pdb(
        mhd=h,
        out=str(save_pdb_file)
    )



