from pathlib import Path
import sys
import numpy as np
import pandas as pd
import random
import time
import argparse

sys.path.append(str(Path(Path.home(), "xray/dev/23_crystal_sim/src")))
import crystal
import unit_cell
import benchmark
import scatterer


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
    n_scatt = 1
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
        scatt_types=["C"],
        force_even=True
    )

    # We have to make weights .5,.5 else the divide will not be clean.
    uc_occs = [.5,.5]

    avg_xyz = np.array([0,0,0])
    scatterers = list()
    for i in range(len(unit_cells)):
        xyz = unit_cells[i].get_scatterers()[0].get_xyz()
        scatt = scatterer.Scatterer(
            xyz=xyz,
            scatt_type='C',
            occ=uc_occs[i]
        )
        scatterers.append(scatt)

    uc_avg = unit_cell.UnitCell(
        name="uc_avg",
        uc_dim=uc_dim,
        scatterers=scatterers
    )

    # for uc in unit_cells:
    #     avg_xyz = avg_xyz + uc.get_scatterers()[0].get_xyz()

    # avg_xyz = avg_xyz / len(unit_cells)

    # avg_xyz = avg_xyz.astype(int)
    # print("avg_xyz: {}".format(avg_xyz))

    # scatt_avg = scatterer.Scatterer(
    #     xyz=avg_xyz,
    #     scatt_type='C',
    #     occ=1.0
    # )

    # uc_avg = unit_cell.UnitCell(
    #     name="uc_{}".format(i),
    #     uc_dim=uc_dim,
    #     scatterers=[scatt_avg]
    # )

    # for uc in unit_cells:
    #     print(uc.get_reflection_fft(hkls=[np.array([1,0,0])]))

    # print(uc_avg.get_reflection_fft(hkls=[np.array([1,0,0])]))


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

        F_avgs = uc_avg.get_reflection_fft(hkls=hkls)

        trial_df = benchmark.get_trial_df(
            F_1s=F_hets,
            F_2s=F_avgs,
            hkls=hkls
        )

        trial_df.rename(columns={"delta_real": "delta_real_{}".format(N), "delta_imag": "delta_imag_{}".format(N)}, inplace=True)
        all_trial_df = pd.merge(all_trial_df, trial_df, left_index=True, right_index=True)

    # Write the trial_df to output_file.
    if not args.test:
        output_file = Path(Path.home(), "xray/dev/23_crystal_sim/data/hetero_vs_degenerate/trial_{}.csv".format(trial_id))
        all_trial_df.to_csv(output_file)