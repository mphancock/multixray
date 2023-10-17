import numpy as np


form_fac_dict = {None: 1, 'H': 1, 'C': 6, 'O': 8, 'N': 14}

class Scatterer():
    def __init__(
            self,
            xyz,
            occ,
            scatt_type=None
        ):
        self.xyz = xyz
        self.occ = occ
        self.scatt_type = scatt_type
        self.form_factor = form_fac_dict[scatt_type]

    def set_xyz(
            self,
            xyz
    ):
        self.xyz = xyz

    def set_occ(
            self,
            occ
    ):
        self.occ = occ

    def set_scatt_type(
            self,
            scatt_type
    ):
        self.type = type
        self.form_factor = form_fac_dict[scatt_type]


    ## set the radius (3D array) of the scatterer ellipsoid
    # def set_radius(self, radius):
    #     scatt_sigma = find_min_std(radius, .999)
    #     self.radius = radius
    #     self.scatt_sigma = scatt_sigma

    def get_xyz(self):
        return self.xyz

    def get_occ(self):
        return self.occ

    def get_type(self):
        return self.scatt_type

    def get_form_factor(self):
        return self.form_factor

    ## sample the scatterer center
    ## check that sampled center does not exceed dim
    # def sample_center(self, dim):
    #     rand_sample = np.random.normal(0, self.pos_sigma, 3)
    #     for i in range(rand_sample.shape[0]):
    #         rand_sample[i] = np.round(rand_sample[i])

    #     center = self.pos_mean + rand_sample
    #     center = correct_bounds(center, dim)

    #     return center.astype(int)

    ## return an array with unit cell dimension with the scatterer as a gaussian randomly sampled around mean position
    ## max radius is the maximum scatterer radius from entire unit cell
    ## informs the dimensions of the returned array
    # def get_scatterer_in_array(self, unit_dim, max_radius):
    #     arr = np.zeros(unit_dim.astype(int) + max_radius.astype(int) * 2)
    #     center = self.sample_center(unit_dim)

    #     ## check if any direction of scattering radius is 0
    #     if (self.radius == np.array([0,0,0])).any():
    #         arr[center[0], center[1], center[2]] = 1 * self.form_fac
    #         return arr

    #     for i in range(arr.shape[0]):
    #         for j in range(arr.shape[1]):
    #             for k in range(arr.shape[2]):
    #                 ## check if coordinate i,j,k is within the ellipsoid around center with radius
    #                 ## then compute scattering density for coordinate i,j,k from guassian centered at center
    #                 if is_within_ellipsoid(np.array([i,j,k]), center + max_radius, self.radius):
    #                     ## gaussian takes 3 params (arrays): X, mu, sigma (covariance matrix)
    #                     ## multiply by atomic form factor
    #                     arr[i,j,k] = arr[i,j,k] + self.form_fac * gaussian(np.array([i,j,k]), center + max_radius, np.identity(3)*self.scatt_sigma**2)

    #     return arr
