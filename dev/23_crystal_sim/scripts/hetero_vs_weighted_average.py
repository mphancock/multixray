from pathlib import Path
import sys
import numpy as np
import pandas as pd
import random
# random.seed(0)
import time
import argparse

sys.path.append(str(Path(Path.home(), "xray/dev/23_crystal_sim/src")))
import crystal
import benchmark


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--trial_id", type=float)
    parser.add_argument("--max_N", type=int)
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()

    if args.test:
        trial_id = -1
    else:
        trial_id = args.trial_id
    uc_dim = np.array([10,10,10])
    n_scatt = 2
    n_state = 2

    # Get all hkls.
    hkls = list()
    if args.test:
        hkls.append(np.array([0,0,0]))
    else:
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    hkls.append(np.array([i,j,k]))

    print("Trial: {}".format(trial_id))
    unit_cells = benchmark.get_n_random_unit_cells(
        n_state=n_state,
        n_scatt=n_scatt,
        uc_dim=uc_dim,
        rand_type=True
    )
    uc_occs = [random.random() for unit_cell in unit_cells]
    uc_occs = [occ/sum(uc_occs) for occ in uc_occs]

    all_trial_df = pd.DataFrame(index=["{}{}{}".format(hkl[0],hkl[1],hkl[2]) for hkl in hkls])
    for N in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
        print(N)
        crystal_het = crystal.Crystal(
            crystal_dim=np.array([N,N,N]),
            unit_cells=unit_cells,
            unit_cell_occs=uc_occs
        )

        # Compute the hkls.
        t0 = time.time()
        F_hets = crystal_het.get_reflection_fft(hkls=hkls)
        print(time.time() - t0)

        F_arr_weight = np.zeros((uc_dim[0],uc_dim[1],uc_dim[2]), dtype=np.complex128)
        for i in range(len(unit_cells)):
            uc_arr = unit_cells[i].build_array()
            F_arr = np.fft.fftn(uc_arr)
            F_arr_weight = F_arr_weight + uc_occs[i]*F_arr

        trial_df = benchmark.get_trial_df(
            F_1s=F_hets,
            F_2s=[F_arr_weight[hkl[0],hkl[1],hkl[2]] for hkl in hkls],
            hkls=hkls
        )

        trial_df.rename(columns={"delta_real": "delta_real_{}".format(N), "delta_imag": "delta_imag_{}".format(N)}, inplace=True)
        all_trial_df = pd.merge(all_trial_df, trial_df, left_index=True, right_index=True)

    # Write the trial_df to output_file.
    if not args.test:
        output_file = Path(Path.home(), "xray/dev/23_crystal_sim/data/hetero_vs_weighted_average/trial_{}.csv".format(trial_id))
        all_trial_df.to_csv(output_file)