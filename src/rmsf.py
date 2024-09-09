import numpy as np

import IMP
import IMP.atom


def get_rmsfs(
    pdb_files,
    n_res
):
    ## read in coords
    mat = np.zeros((len(pdb_files), n_res, 3))
    for i in range(len(pdb_files)):
        pdb_file = pdb_files[i]
        m = IMP.Model()
        h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.CAlphaPDBSelector())

        for j in range(1,n_res+1):
            pid = IMP.atom.Selection(h, residue_index=j, atom_type=IMP.atom.AtomType("CA")).get_selected_particles()[0]
            xyz = IMP.core.XYZR(m, pid)
            mat[i,j-1,:] = xyz.get_coordinates()

    ## normalize
    for i in range(n_res):
        for j in range(3):
            mat[:,i,j] = mat[:,i,j] - np.mean(mat[:,i,j])

    ## calculate fluctuations
    rmsfs = list()
    for i in range(n_res):
        rmsf = 0
        for j in range(len(pdb_files)):
            rmsf += np.linalg.norm(mat[j,i,:])

        rmsf /= len(pdb_files)

        # print(i, rmsf)
        rmsfs.append(rmsf)

    return rmsfs



