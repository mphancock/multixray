import IMP
import IMP.atom


def get_n_state_from_pdb_file(pdb_file):
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
    return len(hs)
