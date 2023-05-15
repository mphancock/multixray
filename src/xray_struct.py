import math
import IMP.core
import IMP.atom
import cctbx.xray
import cctbx.crystal
import cctbx.array_family.flex
import iotbx.pdb


"""
Given an IMP model, the unit cell dimensions, and the space group; the function converts the IMP structure representation into the cctbx structure representation and returns the corresponding cctbx xray structure.
"""
def get_xray_structure(
        m,
        uc_dim,
        sg_symbol
):
    scatterers = cctbx.array_family.flex.xray_scatterer()
    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=uc_dim,
        space_group_symbol=sg_symbol
    )

    unit_cell = crystal_symmetry.unit_cell()

    b_factors = cctbx.array_family.flex.double()
    for pid in m.get_particle_indexes():
        if IMP.atom.Atom.get_is_setup(m, pid):

            d = IMP.core.XYZR(m, pid)
            d.set_coordinates_are_optimized(True)
            coords_cart = (d.get_x(), d.get_y(), d.get_z())
            coords_frac = unit_cell.fractionalize(coords_cart)

            element_table = IMP.atom.ElementTable()
            a = IMP.atom.Atom(m, pid)
            e = a.get_element()
            e_name = element_table.get_name(e)

            b_factor = a.get_temperature_factor()
            # b_factors.append(b_factor)
            u = b_factor / (8 * math.pi**2)

            occ = a.get_occupancy()

            scatterer = cctbx.xray.scatterer(
                label=m.get_particle_name(pid),
                site=coords_frac,
                u=u,
                occupancy=occ,
                scattering_type=e_name,
                fp=0,
                fdp=0
            )
            scatterer.set_use_u(True, False)
            scatterers.append(scatterer)

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