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
import unit_cell
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
    # NEED EVEN UC_OCCS.
    unit_cells = benchmark.get_n_random_unit_cells(
        n_state=n_state,
        n_scatt=n_scatt,
        uc_dim=uc_dim,
        rand_type=True
    )
    uc_occs = [.5,.5]

    # print(unit_cells[0].build_array()[0,0,0])
    # print(unit_cells[1].build_array()[0,0,0])

    # print(np.sum(unit_cells[0].build_array()))
    # print(np.sum(unit_cells[1].build_array()))


    unit_cell_inters = list()
    i = 0
    for scatt_1, scatt_2 in [(unit_cells[0].get_scatterers()[0], unit_cells[1].get_scatterers()[1]), (unit_cells[1].get_scatterers()[0], unit_cells[0].get_scatterers()[1])]:
        uc = unit_cell.UnitCell(
            name="uc_{}".format(i),
            uc_dim=uc_dim,
            scatterers=[scatt_1, scatt_2]
        )
        unit_cell_inters.append(uc)
        i = i+1

    for uc in unit_cells:
        uc.show_scatterers()

    for uc in unit_cell_inters:
        uc.show_scatterers()

    # print(unit_cells[0].get_reflection_fft(hkls=[np.array([0,0,0])]))
    # print(unit_cells[1].get_reflection_fft(hkls=[np.array([0,0,0])]))

    # print(unit_cell_inters[0].get_reflection_fft(hkls=[np.array([0,0,0])]))
    # print(unit_cell_inters[1].get_reflection_fft(hkls=[np.array([0,0,0])]))

    all_trial_df = pd.DataFrame(index=["{}{}{}".format(hkl[0],hkl[1],hkl[2]) for hkl in hkls])
    for N in range(10, args.max_N+1, 10):
        print("N: {}".format(N))
        crystal_het = crystal.Crystal(
            crystal_dim=np.array([N,N,N]),
            unit_cells=unit_cells,
            unit_cell_occs=uc_occs
        )

        # Compute the hkls.
        t0 = time.time()
        F_hets = crystal_het.get_reflection_fft(hkls=hkls)
        print("time: {}".format(time.time() - t0))

        crystal_inter = crystal.Crystal(
            crystal_dim=np.array([N,N,N]),
            unit_cells=unit_cell_inters,
            unit_cell_occs=uc_occs
        )
        F_inters = crystal_inter.get_reflection_fft(hkls=hkls)

        trial_df = benchmark.get_trial_df(
            F_1s=F_hets,
            F_2s=F_inters,
            hkls=hkls
        )

        trial_df.rename(columns={"delta_real": "delta_real_{}".format(N), "delta_imag": "delta_imag_{}".format(N)}, inplace=True)
        all_trial_df = pd.merge(all_trial_df, trial_df, left_index=True, right_index=True)

    # Write the trial_df to output_file.
    if not args.test:
        output_file = Path(Path.home(), "xray/dev/23_crystal_sim/data/hetero_vs_intermixed/trial_{}.csv".format(trial_id))
        all_trial_df.to_csv(output_file)