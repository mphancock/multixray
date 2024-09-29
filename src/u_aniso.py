import mmtbx.f_model
import mmtbx.model
import cctbx.crystal
import cctbx.xray
from iotbx import pdb

import miller_ops


## we need to read in only non altloc atoms to get u aniso becasue that is what IMP will be doing
def read_non_altloc_structure(
    pdb_file
):
    pdb_inp = pdb.input(file_name=str(pdb_file))
    hierarchy = pdb_inp.construct_hierarchy()

    asc = hierarchy.atom_selection_cache()

    # Create a selection for atoms where altloc is not 'B'
    sel = asc.selection("not altloc 'B'")

    # Alternatively, if you want to select atoms where altloc is blank or any character except 'B'
    # sel = asc.selection("altloc '?' and not altloc 'B'")

    # Apply the selection to get a new hierarchy with only the selected atoms
    hierarchy_not_altloc_B = hierarchy.select(sel)

    atoms = hierarchy_not_altloc_B.atoms()

    crystal_symmetry = pdb_inp.crystal_symmetry()

    xray_structure = hierarchy_not_altloc_B.extract_xray_structure(crystal_symmetry=crystal_symmetry)

    return xray_structure


def get_u_anisos_from_file(
    u_aniso_file
):
    # The particles will line up based on the order they are read in from the pdb file.
    xray_struct = read_non_altloc_structure(u_aniso_file)
    scatterers = xray_struct.scatterers()

    u_anisos = list()
    for i in range(len(scatterers)):
        scatterer = scatterers[i]
        print(scatterer.flags.use_u_iso(), scatterer.u_star)
        u_aniso = scatterer.u_star
        u_anisos.append(u_aniso)

    return u_anisos


