import sys 
sys.path.append('../')
from xtal_simulator import Crystal, Unit_Cell, Scatterer 
import numpy as np 
import matplotlib.pyplot as plt 
import sys  
import os 
 
def simulation(dim, scatt_coord, radius): 
    xtal = Crystal()
    xtal.set_n_unit(1)
    xtal.set_unit_dim(dim)

    scatt = Scatterer() 
    scatt.set_pos_mean(scatt_coord)
    scatt.set_pos_sigma(np.array([0,0,0]))
    scatt.set_radius(radius)

    uc = Unit_Cell()
    uc.add_scatterer(scatt)
    xtal.create_unit_cell(uc, 1)

    freq = xtal.FT()
    return freq 

def compute_rmsd(dim_num, scatt_coord, radius): 
    dim = np.array([dim_num]*3) 
    freq = simulation(dim, scatt_coord, radius) 

    rmsd = 0 
    for h in range(dim_num): 
        for k in range(dim_num): 
            for l in range(dim_num): 
                sf_num = freq[h,k,l] 
                sf_ana = np.e**(-2*np.pi*np.array([1j])*np.dot(np.array([h,k,l]), scatt_coord/dim))
                rmsd = rmsd + (sf_num - sf_ana)**2

    rmsd = rmsd / dim_num**3
    rmsd = rmsd ** (1/2)

    return rmsd 


if __name__ == '__main__': 
    dim_num = 100  
    N = 25 
    rmsd_mu = list() 
    rmsd_sig = list()

    rmsd_mu_0 = list() 
    rmsd_sig_0 = list() 
    for i in range(5,dim_num): 
        print(i)
        rmsd = list() 
        rmsd_0 = list() 
        for j in range(N): 
            rand = np.random.rand(3)

            dim = np.array([i]*3) 
            scatt_coord = np.floor(rand * dim).astype(int)
            rmsd.append(np.absolute(compute_rmsd(i,scatt_coord, np.array([2,2,2]))))
            rmsd_0.append(np.absolute(compute_rmsd(i,scatt_coord, np.array([0,0,0]))))

        rmsd_mu.append(np.mean(rmsd))
        rmsd_sig.append(np.std(rmsd)) 

        rmsd_mu_0.append(np.mean(rmsd_0))
        rmsd_sig_0.append(np.std(rmsd_0)) 

    plot_path = '/wynton/home/sali/mhancock/xtal/results/xtal_simulation/plots'
    # plot_path = '/Volumes/hard_drive/scratch/'

    gauss = plt.errorbar(range(5,dim_num), rmsd_mu, yerr=rmsd_sig/np.sqrt(N), fmt='bo', capsize=5) 
    point = plt.errorbar(range(5,dim_num), rmsd_mu_0, yerr=rmsd_sig_0/np.sqrt(N), fmt='ro', capsize=5) 
    plt.title('RMSD of numerical/analytical Fhkl (1 UC) vs. dimensions of UC (N = {})'.format(N))
    plt.xlabel('dimension of UC') 
    plt.ylabel('Fhkl RMSD')
    plt.legend([gauss, point], ['Gaussian (r=2)', 'point'])

    plt.savefig(os.path.join(plot_path, 'fhkl_num_vs_analytical_1_uc.png'))