from xtal_simulator import Scatterer, Unit_Cell, Crystal 
import numpy as np 
import matplotlib.pyplot as plt 
import os 
import sys 


def compute_sf_rmsd(sf_A, sf_B, n_unit):
    rmsd = np.absolute(sf_A - sf_B)**2
    rmsd = np.sum(rmsd) 
    rmsd = rmsd / (sf_A.shape[0] * sf_A.shape[1] * sf_A.shape[2]) 
    rmsd = rmsd**(1/2)
    return rmsd 

def compute_sf_avg(sf_A, sf_B): 
    mean = np.sum(sf_A - sf_B) 
    mean = np.absolute(mean) 
    mean = mean / (sf_A.shape[0] * sf_A.shape[1] * sf_A.shape[2]) 

    return mean 

def single_scatt(plot_path, pos_std, max_n_unit): 
    unit_dim = np.array([10,10,10])
    rmsd_mu = {0: list()} 
    rmsd_sig = {0: list()} 

    mean_mu = list() 
    mean_sig = list() 
    N = 20 
    for i in range(1, max_n_unit): 
        n_unit = i   

        mean_n_unit = list() 
        rmsd_n_unit = {0: list()}
        for j in range(N): 
            pos_A = np.array(np.floor(np.random.rand(3) * unit_dim))

            scatt_1 = Scatterer() 
            scatt_1.set_pos_mean(pos_A)
            scatt_1.set_pos_sigma(0)

            xtal = Crystal()
            xtal.set_n_unit(n_unit) 
            xtal.set_unit_dim(unit_dim)

            uc_1 = Unit_Cell() 
            uc_1.set_dim(unit_dim)
            uc_1.add_scatterer(scatt_1)

            sf_1 = uc_1.compute_sf_fft()

            uc_1.scatterers[0].set_pos_sigma(pos_std)
            xtal.create_unit_cell(uc_1, 1) 
            sf_mix = xtal.compute_sf() / xtal.get_n_unit()**3

            rmsd_n_unit[0].append(compute_sf_rmsd(sf_1, sf_mix, n_unit))
            mean_n_unit.append(compute_sf_avg(sf_1, sf_mix))

        rmsd_mu[0].append(np.mean(rmsd_n_unit[0]))
        rmsd_sig[0].append(np.std(rmsd_n_unit[0]))

        mean_mu.append(np.mean(mean_n_unit)) 
        mean_sig.append(np.std(mean_n_unit))

    p1 = plt.errorbar(range(1, max_n_unit), rmsd_mu[0], yerr=rmsd_sig[0]/np.sqrt(N), fmt='bo', capsize=5) 
    # p1 = plt.errorbar(range(1, max_n_unit), mean_mu, yerr=mean_sig/np.sqrt(N), fmt='bo', capsize=5) 
    plt.title('Fhkl RMSD disordered/mean order (N={}, sig={})'.format(N, pos_std))
    plt.legend([p1], ['mix vs mean avg'], loc='upper right')
    plt.xlabel('num UC') 
    plt.ylabel('Fhkl RMSD')
    plt.savefig(os.path.join(plot_path, 'fhkl_order_vs_dis_single_scatt_sigma_{}.png'.format(pos_std)))
    # plt.show()

def multi_scatt(plot_path, pos_std, max_n_unit): 
    unit_dim = np.array([10,10,10])
    rmsd_mu = {0: list(), 1: list(), 2:list()} 
    rmsd_sig = {0: list(), 1: list(), 2:list()} 
    N = 20 
    for i in range(1, max_n_unit): 
        n_unit = i   

        rmsd_n_unit = {0: list(), 1: list(), 2:list()}
        for j in range(N): 
            occ = np.random.rand() 

            pos_A = np.array(np.floor(np.random.rand(3) * unit_dim))
            pos_B = np.array(np.floor(np.random.rand(3) * unit_dim))

            while ((pos_A == pos_B).all() == True): 
                pos_B = np.array(np.floor(np.random.rand(3) * unit_dim))

            scatt_1 = Scatterer() 
            scatt_1.set_pos_mean(pos_A)
            scatt_1.set_pos_sigma(0)

            scatt_2 = Scatterer() 
            scatt_2.set_pos_mean(pos_B)
            scatt_2.set_pos_sigma(0)

            xtal = Crystal()
            xtal.set_n_unit(n_unit) 
            xtal.set_unit_dim(unit_dim)

            uc_1 = Unit_Cell() 
            uc_1.set_dim(unit_dim)
            uc_1.add_scatterer(scatt_1)

            uc_2 = Unit_Cell()
            uc_2.set_dim(unit_dim)
            uc_2.add_scatterer(scatt_2)

            sf_1 = uc_1.compute_sf_fft()
            sf_2 = uc_2.compute_sf_fft()
            sf_avg = sf_1 * occ + sf_2 * (1-occ) 

            uc_1.scatterers[0].set_pos_sigma(pos_std)
            uc_2.scatterers[0].set_pos_sigma(pos_std)

            xtal.create_unit_cell(uc_1, occ) 
            xtal.create_unit_cell(uc_2, 1-occ)
            sf_mix = xtal.compute_sf() / xtal.get_n_unit()**3

            rmsd_n_unit[0].append(compute_sf_rmsd(sf_mix, sf_avg, n_unit))
            rmsd_n_unit[1].append(compute_sf_rmsd(sf_mix, sf_1, n_unit))
            rmsd_n_unit[2].append(compute_sf_rmsd(sf_mix, sf_2, n_unit))

        rmsd_mu[0].append(np.mean(rmsd_n_unit[0]))
        rmsd_mu[1].append(np.mean(rmsd_n_unit[1]))
        rmsd_mu[2].append(np.mean(rmsd_n_unit[2]))

        rmsd_sig[0].append(np.std(rmsd_n_unit[0]))
        rmsd_sig[1].append(np.std(rmsd_n_unit[1]))
        rmsd_sig[2].append(np.std(rmsd_n_unit[2]))

    p1 = plt.errorbar(range(1, max_n_unit), rmsd_mu[0], yerr=rmsd_sig[0]/np.sqrt(N), fmt='bo', capsize=5) 
    p2 = plt.errorbar(range(1, max_n_unit), rmsd_mu[1], yerr=rmsd_sig[1]/np.sqrt(N), fmt='ro', capsize=5) 
    p3 = plt.errorbar(range(1, max_n_unit), rmsd_mu[2], yerr=rmsd_sig[2]/np.sqrt(N), fmt='go', capsize=5) 
    plt.title('Fhkl RMSD disordered/occ weighted order vs. UC dim (N={}, sig={})'.format(N, pos_std))
    plt.xlabel('num UC') 
    plt.ylabel('Fhkl RMSD')
    plt.legend([p1, p2, p3], ['dis vs weight avg', 'dis vs order 1', 'dis vs order 2'], loc='upper right')
    plt.savefig(os.path.join(plot_path, 'fhkl_order_vs_dis_multi_scat_sigma_{}.png'.format(pos_std)))
    # plt.show()

if __name__ == '__main__': 
    plot_path = '/wynton/home/sali/mhancock/xtal/results/xtal_simulation/plots'
    # plot_path = '/Volumes/hard_drive/scratch/'
    max_n_unit = 50 

    if int(sys.argv[1]) == 0: 
        single_scatt(plot_path, int(sys.argv[2]), max_n_unit)
    if int(sys.argv[1]) == 1: 
        multi_scatt(plot_path, int(sys.argv[2]), max_n_unit)