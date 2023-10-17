import numpy as np
import random
from distributions import gaussian, find_min_std
from math_func import is_within_ellipsoid, correct_bounds
from aux_functions import get_coords_lists, get_color
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/crystal_sim/src")))
import ft


class Crystal:
    def __init__(
            self,
            crystal_dim,
            unit_cells,
            unit_cell_occs
    ):
        # The number of unit cells in each direction.
        self.crystal_dim = crystal_dim
        self.unit_cells = unit_cells
        self.uc_dim = unit_cells[0].get_uc_dim()
        for uc in self.unit_cells:
            if (uc.get_uc_dim() != self.uc_dim).any():
                raise Exception("Unit cell dimensions {} do not match crystal dimensions {}".format(uc.get_uc_dim(), self.uc_dim))

        self.occs = unit_cell_occs

    # def get_largest_scatt_radius(self):
    #     max_radius = np.zeros(3)
    #     for i in range(len(self.unit_cells)):
    #         radius = self.unit_cells[i].get_max_scatt_radius()
    #         if radius[0] > max_radius[0]:
    #             max_radius[0] = radius[0]
    #         if radius[1] > max_radius[1]:
    #             max_radius[1] = radius[1]
    #         if radius[2] > max_radius[2]:
    #             max_radius[2] = radius[2]

    #     return max_radius.astype(int)

    ## scattering coords must be 2D numpy array with dim = (# of scatt, 3)
    def add_unit_cell(
        self,
        uc,
        occ
    ):
        if (uc.get_uc_dim() != self.uc_dim).any():
            raise Exception("Unit cell dimensions {} do not match crystal dimensions {}".format(uc.get_uc_dim(), self.uc_dim))

        self.unit_cells.append(uc)
        self.occs.append(occ)

    # def get_total_occ(self):
    #     total_occ = 0
    #     for i in range(len(self.unit_cells)):
    #         total_occ = total_occ + self.unit_cells[i].get_occupancy()

    #     return total_occ

    # def sample_unit_cells(self):
    #     random = np.random.random()
    #     total_occ = self.get_total_occ()
    #     cum_occ = 0
    #     i = 0
    #     while cum_occ < random:
    #         cum_occ = cum_occ + self.unit_cells[i].get_occupancy()
    #         i = i + 1

    #     return i - 1


    def get_reflection(
            self,
            hkl
    ):
        for i in range(3):
            if hkl[i] < 0 or hkl[i] > self.crystal_dim[i]:
                raise RuntimeError("hkl must be within crystal dimensions")

        F = 0 + 0j

        pool_params = list()
        for i in range(self.crystal_dim[0]):
            for j in range(self.crystal_dim[1]):
                for k in range(self.crystal_dim[2]):
                    uc = random.choices(self.unit_cells, self.occs, k=1)[0]
                    x = uc.build_array()

                    params_dict = {}
                    params_dict["x"] = x
                    params_dict["x_0"] = i*self.uc_dim[0]
                    params_dict["y_0"] = j*self.uc_dim[1]
                    params_dict["z_0"] = k*self.uc_dim[2]
                    params_dict["x_N"] = self.crystal_dim[0]*self.uc_dim[0]
                    params_dict["y_N"] = self.crystal_dim[1]*self.uc_dim[1]
                    params_dict["z_N"] = self.crystal_dim[2]*self.uc_dim[2]
                    params_dict["h"] = hkl[0]*self.crystal_dim[0]
                    params_dict["k"] = hkl[1]*self.crystal_dim[1]
                    params_dict["l"] = hkl[2]*self.crystal_dim[2]

                    pool_params.append(params_dict)

        print("CPUs: {}".format(multiprocessing.cpu_count()))
        pool_obj = multiprocessing.Pool(
            multiprocessing.cpu_count()
        )

        pool_results = pool_obj.imap(
            ft.dft3D_with_offset_for_hkl_pool,
            pool_params
        )

        for F_uc in pool_results:
            F = F + F_uc

        return F

    def get_reflection_fft(
            self,
            hkls
    ):
        cryst = np.ndarray([self.crystal_dim[0]*self.uc_dim[0], self.crystal_dim[1]*self.uc_dim[1], self.crystal_dim[2]*self.uc_dim[2]])
        for i in range(self.crystal_dim[0]):
            for j in range(self.crystal_dim[1]):
                for k in range(self.crystal_dim[2]):
                    uc = random.choices(self.unit_cells, self.occs, k=1)[0]
                    if i == 0 and j == 0 and k == 0:
                        print(uc.get_name())
                    x = uc.build_array()
                    cryst[i*self.uc_dim[0]:(i+1)*self.uc_dim[0], j*self.uc_dim[1]:(j+1)*self.uc_dim[1], k*self.uc_dim[2]:(k+1)*self.uc_dim[2]] = x

        Fs = list()

        cryst_ft = np.fft.fftn(cryst)
        for hkl in hkls:
            F = cryst_ft[hkl[0]*self.crystal_dim[0], hkl[1]*self.crystal_dim[1], hkl[2]*self.crystal_dim[2]] / (self.crystal_dim[0]*self.crystal_dim[1]*self.crystal_dim[2])
            Fs.append(F)

        return Fs



        # crystal_dim = (self.unit_dim * self.n_unit + 2 * self.get_largest_scatt_radius()).astype(int)
        # crystal = np.zeros(crystal_dim)
        # # print(crystal.shape)

        # for i in range(self.n_unit):
        #     for j in range(self.n_unit):
        #         for k in range(self.n_unit):
        #             latt_pt = np.array([i * self.unit_dim[0], j * self.unit_dim[1], k * self.unit_dim[2]])
        #             latt_pt = latt_pt + self.get_largest_scatt_radius()
        #             # rand = random.random()

        #             index = self.sample_unit_cells()
        #             curr_unit_cell = self.unit_cells[index]

        #             padding = curr_unit_cell.get_max_scatt_radius()

        #             # print(latt_pt[0]-padding[0], latt_pt[0]+self.unit_dim[0]+padding[0])

        #             present_val = crystal[latt_pt[0]-padding[0]:latt_pt[0]+self.unit_dim[0]+padding[0], \
        #                 latt_pt[1]-padding[1]:latt_pt[1]+self.unit_dim[1]+padding[1], \
        #                 latt_pt[2]-padding[2]:latt_pt[2]+self.unit_dim[2]+padding[2]].copy()


        #             crystal[latt_pt[0]-padding[0]:latt_pt[0]+self.unit_dim[0]+padding[0], \
        #                 latt_pt[1]-padding[1]:latt_pt[1]+self.unit_dim[1]+padding[1], \
        #                 latt_pt[2]-padding[2]:latt_pt[2]+self.unit_dim[2]+padding[2]] \
        #                 = present_val + curr_unit_cell.build_array()

        # print(crystal.shape)
        # print(np.sum(crystal))
        # freq = np.fft.fftn(crystal)

        # return freq

#     def compute_sf(self):
# #         hkl_corr = hkl * self.n_unit
# #        sf = 0
#         # for i in range(self.n_unit):
#         #     for j in range(self.n_unit):
#         #         for k in range(self.n_unit):
#         #             latt_vec = np.array([i,j,k]) * self.unit_dim
#         #             uc = self.unit_cells[self.sample_unit_cells()]
#         #             sf = sf + uc.compute_ft(latt_vec, hkl_corr, self.unit_dim * self.n_unit)

#         freq = self.FT()
#         # print(freq.shape)
#         sf = freq[::self.get_n_unit(),::self.get_n_unit(),::self.get_n_unit()]
#         # print(sf.shape)

#         return sf


        # return None

    # def FT_danielson_lanczos(self, n=1):
    #     for i in range()
