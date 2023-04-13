import iotbx.pdb


def get_crystal_info(pdb_file):
    pdb_input = iotbx.pdb.hierarchy.input(file_name=str(pdb_file))
    structure = pdb_input.xray_structure_simple()
    symmetry = structure.crystal_symmetry()
    uc = symmetry.unit_cell()
    uc_dimensions = uc.parameters()
    sg_info = symmetry.space_group_info()
    sg_symbol = sg_info.symbol_and_number().split("(")[0]

    return uc_dimensions, sg_symbol