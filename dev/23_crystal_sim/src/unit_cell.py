import numpy as np


class UnitCell:
    def __init__(
            self,
            name,
            uc_dim,
            scatterers
    ):
        self.name = name
        self.uc_dim = uc_dim
        self.scatterers = scatterers

    def get_name(self):
        return self.name

    def get_uc_dim(self):
        return self.uc_dim

    def get_scatterers(self):
        return self.scatterers

    def set_uc_dim(self, uc_dim):
        self.uc_dim = uc_dim

    ## add a single scatterer to the unit cell
    def add_scatterer(self, scatterer):
        self.scatterers.append(scatterer)

    ## add a list of scatterers to the unit cell
    def add_scatterers(self, scatterers):
        for scatt in scatterers:
            self.add_scatterer(scatt)

    ## get_scatterer_in_array returns array that is larger than unit dim
    ## returned array is unit dim + radius
    ## determine what the largest array returned will be for unit cell
    # def get_max_scatt_radius(self):
    #     max_radius = np.zeros(3)
    #     for scatterer in self.scatterers:
    #         radius = scatterer.get_radius()
    #         if radius[0] > max_radius[0]:
    #             max_radius[0] = radius[0]
    #         if radius[1] > max_radius[1]:
    #             max_radius[1] = radius[1]
    #         if radius[2] > max_radius[2]:
    #             max_radius[2] = radius[2]

    #     return max_radius.astype(int)

    def build_array(self):
        arr = np.zeros(self.uc_dim)
        for scatterer in self.scatterers:
            xyz = scatterer.get_xyz()
            arr[xyz[0], xyz[1], xyz[2]] = scatterer.get_form_factor() * scatterer.get_occ()

        return arr

    def get_reflection_fft(
            self,
            hkls
    ):
        Fs = list()
        uc_rs = self.build_array()
        # print(sum(uc_rs))
        uc_fft = np.fft.fftn(uc_rs)
        for hkl in hkls:
            F = uc_fft[hkl[0], hkl[1], hkl[2]]
            Fs.append(F)

        return Fs

    def show_scatterers(
            self
    ):
        for i in range(len(self.get_scatterers())):
            print("scatt {}: {} {}".format(i, self.get_scatterers()[i].get_xyz(), self.get_scatterers()[i].get_type()))


    # def compute_ft(self, freq):
    #     arr = self.build_array()
    #     ft = 0
    #     for i in range(arr.shape[0]):
    #         for j in range(arr.shape[1]):
    #             for k in range(arr.shape[2]):
    #                 ft = ft + arr[i,j,k] * np.e**(-2 * np.pi * np.array([1j]) * np.dot(freq, np.array([i,j,k])/self.get_dim()))

    #     return ft

    ## compute first num_hkl structure factors in all 3 dimensions
    # def compute_sf_ana(self, num_hkl):
    #     sf = np.ndarray(num_hkl, dtype='complex')
    #     for h in range(num_hkl[0]):
    #         for k in range(num_hkl[1]):
    #             for l in range(num_hkl[2]):
    #                 sf[h,k,l] = self.compute_ft(np.array([h,k,l]))

    #     return sf

    # def compute_sf_fft(self):
    #     arr = self.build_array()
    #     sf = np.fft.fftn(arr)
    #     return sf

    # def get_scatterers_by_type(self):
    #     type_dict = dict()
    #     for scatt in self.scatterers:
    #         scatt_type = scatt.get_type()
    #         if scatt_type in list(type_dict.keys()):
    #             type_dict[scatt_type].append(scatt)
    #         else:
    #             type_dict[scatt_type] = [scatt]

    #     return type_dict

    # def show_graph(self, proj=False, color=False, bonds=False):
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')

    #     if not color:
    #         xs, ys, zs = get_coords_lists(self.scatterers)
    #         ax.scatter(xs, ys, zs, c='r', marker='o')
    #         # ax.scatter(xs_proj, ys_proj, zs_proj, c='b', marker='o')
    #     else:
    #         type_dict = self.get_scatterers_by_type()
    #         for scatt_type in list(type_dict.keys()):
    #             xs, ys, zs = get_coords_lists(type_dict[scatt_type])
    #             ax.scatter(xs, ys, zs, c=get_color(scatt_type), edgecolors='k', marker='o')

    #     ax.set_xlim([0, self.get_dim()[0]])
    #     ax.set_ylim([0, self.get_dim()[1]])
    #     ax.set_zlim([0, self.get_dim()[2]])

    #     if bonds:
    #         for scatt_1 in self.scatterers:
    #             for scatt_2 in self.scatterers:
    #                 pos_1 = scatt_1.get_pos_mean()
    #                 pos_2 = scatt_2.get_pos_mean()
    #                 if not (pos_1 == pos_2).all() and scatt_1.get_type() != 'H' and scatt_2.get_type() != 'H':
    #                     if np.sum((pos_1 - pos_2)**2) < 3:
    #                         ax.plot([pos_1[0], pos_2[0]], [pos_1[1], pos_2[1]], [pos_1[2], pos_2[2]], c='k')

    #     # ax.plot([0,0,0], [5,5,5])
    #     ax.set_xlabel('x')
    #     ax.set_ylabel('y')
    #     ax.set_zlabel('z')

    #     return fig
