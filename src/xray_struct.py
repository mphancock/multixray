from pathlib import Path
import math

import IMP.core
import IMP.atom

import cctbx.xray
import cctbx.crystal
import cctbx.array_family.flex
from mmtbx.tls import tools
import mmtbx.model
import iotbx.pdb
from cctbx.array_family import flex


"""
Given an IMP model, the unit cell dimensions, and the space group; the function converts the IMP structure representation into the cctbx structure representation and returns the corresponding cctbx xray structure.
"""
def get_xray_structure(
        hs,
        occs,
        pids,
        crystal_symmetry,
        u_anisos=None,
        delta=None
):
    m = hs[0].get_model()

    scatterers = cctbx.array_family.flex.xray_scatterer()
    unit_cell = crystal_symmetry.unit_cell()

    ## do a check that the names are equal across all the states
    names = list()
    cnt = 0
    for i in range(len(hs)):
        h = hs[i]

        if i == 0:
            names = [m.get_particle_name(pid) for pid in IMP.atom.Selection(h).get_selected_particle_indexes()]

        # Get only subset of pids that are in th state corresponding to h.
        pids_state = IMP.atom.Selection(h).get_selected_particle_indexes()
        # pids_state_subset = list(set(pids).intersection(set(pids_state)))

        for j in range(len(pids_state)):
            ## only check the name of the first particle in the state because there may be a different number of solvent atoms
            if j == 0 and names[j] != m.get_particle_name(pids_state[j]):
                raise RuntimeError("The names of the particles across states don't match")

            pid = pids_state[j]
            d = IMP.core.XYZR(m, pid)
            coords_cart = (d.get_x(), d.get_y(), d.get_z())
            if delta:
                coords_cart = (coords_cart[0]+delta[0], coords_cart[1]+delta[1], coords_cart[2]+delta[2])

            coords_frac = unit_cell.fractionalize(coords_cart)

            element_table = IMP.atom.ElementTable()
            a = IMP.atom.Atom(m, pid)
            e = a.get_element()
            e_name = element_table.get_name(e)

            occ = occs[i]
            b_factor = a.get_temperature_factor()

            # print(i, cnt, m.get_particle_name(pid))
            # cnt = cnt + 1

            scatterer = cctbx.xray.scatterer(
                label=m.get_particle_name(pid),
                site=coords_frac,
                b=b_factor,
                occupancy=occ,
                scattering_type=e_name,
                fp=0,
                fdp=0
            )

            # print(i)
            if u_anisos:
                scatterer.set_use_u(False, True)
                u_star = u_anisos[pid]
                scatterer.u_star = u_star
            else:
                scatterer.set_use_u(True, False)

            scatterers.append(scatterer)
            cnt += 1

    xray_structure = cctbx.xray.structure(
        scatterers=scatterers,
        crystal_symmetry=crystal_symmetry
    )

    return xray_structure


"""
Given a set of single structure pdb files, the corresponding weights, and the cctbx crystal symmetry; the function returns a cctbx xray structure containing multiple structures.
"""
def get_cctbx_multi_structure_from_pdbs(
        pdb_files,
        weights,
        crystal_symmetry
):
    if len(weights) != len(pdb_files):
        raise RuntimeError("Length of pdb_files must match the length of weights")

    w_sum = 0
    for w in weights:
        w_sum = w_sum + w
    if not (w_sum > .999) and not (w_sum < 1.001):
        raise RuntimeError("Sum of weights must be 1")

    multi_structure = iotbx.pdb.input(str(pdb_files[0])).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )
    multi_structure.adjust_occupancy(
        occ_max=weights[0],
        occ_min=weights[0],
        selection=multi_structure.all_selection()
    )

    for i in range(1, len(pdb_files)):
        pdb_file = pdb_files[i]
        w = weights[i]
        tmp_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
            crystal_symmetry=crystal_symmetry
        )

        tmp_structure.adjust_occupancy(
            occ_max=w,
            occ_min=w,
            selection=tmp_structure.all_selection()
        )

        tmp_scatterers = tmp_structure.scatterers()
        multi_structure.add_scatterers(
            scatterers=tmp_scatterers
        )

    return multi_structure


if __name__ == "__main__":
    u_aniso_file = Path(Path.home(), "xray/dev/26_phenix_refine/data/7mhj_heavy/7mhj_heavy_refine_001.pdb")
    # get_u_aniso_from_file(u_aniso_file)

    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhj_heavy.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=(114.880, 54.736, 45.240, 90.00, 101.42, 90.00),
        space_group_symbol="P1"
    )

    xray_struct = iotbx.pdb.input(str(u_aniso_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    # xray_struct = get_xray_structure(
    #     m,
    #     pids,
    #     crystal_symmetry=crystal_symmetry,
    #     u_aniso_file=u_aniso_file
    # )

    # print(xray_struct.scatterers()[0].u_star)
    print(xray_struct.show_scatterers())
    print(xray_struct.show_scatterer_flags_summary())
