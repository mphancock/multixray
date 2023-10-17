from pathlib import Path
import sys
import pandas as pd
import multiprocessing
import numpy as np
import pandas as pd

import IMP
import IMP.atom

from pyRMSD.utils.proteinReading import Reader
from pyRMSD.matrixHandler import MatrixHandler


if __name__ == "__main__":
    multi_pdb_file = Path(Path.home(), "xray/tmp/sample.pdb")

    pdb_meta_file = Path(Path.home(), "xray/sample_bench/data/7mhf/59_7mhj_1/stat_df_xray_0.csv")
    pdb_file_df = pd.read_csv(pdb_meta_file)
    pdb_files = list(pdb_file_df["xray_0_min_0_pdb"])

    # reader = Reader().readThisFile(str(multi_pdb_file)).gettingOnlyCAs()
    # coordinates = reader.read()
    # num_of_atoms = reader.numberOfAtoms
    # num_of_frames = reader.numberOfFrames

    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(multi_pdb_file), m, IMP.atom.CAlphaPDBSelector())

    n_struct = len(hs)
    n_atoms = len(IMP.atom.Selection(hs[0]).get_selected_particle_indexes())
    all_coords = np.ndarray(shape=(n_struct, n_atoms, 3))
    for i in range(n_struct):
        h = hs[i]
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()

        for j in range(n_atoms):
            pid = pids[j]
            xyz = IMP.core.XYZ(m, pid)
            for k in range(3):
                all_coords[i, j, k] = xyz.get_coordinates()[k]

    rmsd_matrix = MatrixHandler()\
    .createMatrix(all_coords, 'QCP_OMP_CALCULATOR')

    avg_struct_id = 0
    min_total_rmsd = np.infty

    for i in range(len(rmsd_matrix)):
        sum = 0
        for j in range(len(rmsd_matrix)):
            if i != j:
                sum = sum + rmsd_matrix[i,j]
        print(i, sum)

        if sum < min_total_rmsd:
            min_total_rmsd = sum
            avg_struct_id = i

    print(pdb_files[avg_struct_id], min_total_rmsd)
