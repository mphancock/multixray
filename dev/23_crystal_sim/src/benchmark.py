import pandas as pd
import random
from pathlib import Path
import sys
import numpy as np

sys.path.append(str(Path(Path.home(), "xray/dev/23_crystal_sim/src")))
import scatterer
import unit_cell


def get_n_random_unit_cells(
        n_state,
        n_scatt,
        uc_dim,
        scatt_types,
        force_even=False
):
    unit_cells = list()
    for i in range(n_state):
        scatts = list()
        for j in range(n_scatt):
            scatt_type = np.random.choice(scatt_types)
            # if rand_type:
            #     scatt_type = random.choice(['H', 'C', 'O', 'N'])
            # else:
            #     scatt_type = 'H'

            xyz = np.zeros(3)
            for k in range(3):
                xyz[k] = random.randint(0, uc_dim[k]-1)
                if force_even:
                    while xyz[k] % 2 != 0:
                        xyz[k] = random.randint(0, uc_dim[k]-1)
            xyz = xyz.astype(int)

            scatt = scatterer.Scatterer(
                xyz=xyz,
                scatt_type=scatt_type,
                occ=1.0
            )
            print("scatt xyz: {}".format(scatt.get_xyz()))
            scatts.append(scatt)

        uc = unit_cell.UnitCell(
            name="uc_{}".format(i),
            uc_dim=uc_dim,
            scatterers=scatts
        )

        unit_cells.append(uc)

    return unit_cells

def get_trial_df(
        F_1s,
        F_2s,
        hkls
):
    trial_df = pd.DataFrame()
    for i in range(len(hkls)):
        hkl = hkls[i]
        F_1 = F_1s[i]
        F_2 = F_2s[i]

        delta_real = (F_1-F_2).real
        delta_imag = (F_1-F_2).imag

        if (hkl == np.array([0,0,0])).all() or (hkl == np.array([1,0,0])).all():
            print(hkl, F_1, F_2)

        trial_df.loc["{}{}{}".format(hkl[0],hkl[1],hkl[2]), "delta_real"] = delta_real
        trial_df.loc["{}{}{}".format(hkl[0],hkl[1],hkl[2]), "delta_imag"] = delta_imag

    return trial_df