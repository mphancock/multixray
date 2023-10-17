from pathlib import Path
import sys
import numpy as np
import random
random.seed(0)

sys.path.append(str(Path(Path.home(), "xray/crystal_sim/src")))
import ft


if __name__ == "__main__":
    ucs = list()
    ijks = list()
    cryst = np.zeros([4,4,4])

    for i in range(2):
        for j in range(2):
            for k in range(2):
                uc = np.zeros([2,2,2])
                uc[i,j,k] = random.choice([0,1])
                ucs.append(uc)
                ijks.append((i,j,k))

                cryst[i*2:(i+1)*2, j*2:(j+1)*2, k*2:(k+1)*2] = uc

    fft_cryst = np.fft.fftn(cryst)

    uc_sum_dft = np.zeros([4,4,4],dtype=complex)
    for i in range(len(ucs)):
        uc = ucs[i]
        ijk = ijks[i]
        x_0, y_0, z_0 = ijk[0]*2, ijk[1]*2, ijk[2]*2
        uc_dft = ft.dft3D_with_offset(uc, x_0, y_0, z_0, x_N=4, y_N=4, z_N=4)

        uc_sum_dft = uc_sum_dft + uc_dft

    print(cryst[0,0,0], ucs[0][0,0,0])

    for i in range(4):
        for j in range(4):
            for k in range(4):
                print(fft_cryst[i,j,k], uc_sum_dft[i,j,k])

    # rec_latt_cryst = np.fft.fftn(test_cryst)
    # print(rec_latt_uc_sum[0,0], rec_latt_cryst[0,0])
    # print(rec_latt_uc_sum[0,2], rec_latt_cryst[0,2])
    # print(rec_latt_uc_sum[2,0], rec_latt_cryst[2,0])
    # print(rec_latt_uc_sum[2,2], rec_latt_cryst[2,2])

