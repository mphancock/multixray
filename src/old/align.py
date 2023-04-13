from pathlib import Path
import Bio.PDB


def test_input(align_input):
    # pdb_file, ref_pdb_file, align_file = align_input

    aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    ref_atoms = list()
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_model = pdb_parser.get_structure("reference", str(align_input["pdb"]))[0]
    for chain in ref_model.get_chains():
        for res in chain:
            if res.resname in aas:
                ref_atoms.append(res['CA'])

    sample_atoms = list()
    try:
        sample_structure = pdb_parser.get_structure("sample", str(align_input["pdb"]))
    except Bio.PDB.PDBExceptions.PDBConstructionException:
        return Bio.PDB.PDBExceptions.PDBConstructionException
    except ValueError:
        return ValueError

    try:
        sample_model = [model for model in sample_structure.get_models()][0]
    except IndexError:
        return IndexError

    for res in sample_structure.get_residues():
        if res.resname in aas:
            try:
                sample_atoms.append(res['CA'])
            except KeyError:
                return KeyError

    if len(ref_atoms) != len(sample_atoms):
        return AssertionError

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    if super_imposer.rms > 15:
        return AssertionError

    return None


def align(align_input):
    # pdb_file, ref_pdb_file, align_file, backbone = align_input
    error = test_input(align_input)
    if error:
        return align_input["pdb"], error

    # pdb_file_0, code = test_input(align_input)
    # if code:
    #     raise RuntimeError("Invalid input: {}".format(pdb_file))

    aas = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    ref_atoms = list()
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_model = pdb_parser.get_structure("reference", str(align_input["pdb"]))[0]
    for chain in ref_model.get_chains():
        for res in chain:
            if res.resname in aas:
                if align_input["ca"]:
                    ref_atoms.append(res['CA'])
                else:
                    for atom in res:
                        ref_atoms.append(atom)

    sample_atoms = list()
    sample_structure = pdb_parser.get_structure("sample", str(align_input["pdb"]))
    # sample_model = sample_structure.get_models()[0]
    # for chain in sample_model.get_chains():
    #     for res in chain:
    #         if res.resname in aas:
    #             sample_atoms.append(res['CA'])
    for res in sample_structure.get_residues():
        if res.resname in aas:
            if align_input["ca"]:
                sample_atoms.append(res['CA'])
            else:
                for atom in res:
                    sample_atoms.append(atom)

    sample_model = [model for model in sample_structure.get_models()][0]
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    if align_input["save"]:
        io = Bio.PDB.PDBIO()
        io.set_structure(sample_structure)
        io.save(str(align_input["save"]))

    # print(pdb_file, super_imposer.rms)
    return align_input["pdb"], super_imposer.rms

