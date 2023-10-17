import sys 
sys.path.append('../lib')

import os 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R 

from xtal_simulator import Crystal, Unit_Cell 
from amino_acid import add_threonine, add_leucine
from aux_functions import compute_sf_rmsd 

max_n_unit = 50 
N = 5  
uc_1 = Unit_Cell(np.array([20,10,10])) 
uc_2 = Unit_Cell(np.array([20,10,10]))

r_rot = R.from_quat([0, np.sin(np.pi/4), 0, np.cos(np.pi/4)])

uc_1 = add_leucine(np.array([5,5,5]), uc_1, 0, R.from_quat([0,0,0,1])) 
uc_1 = add_leucine(np.array([8,6,5]), uc_1, 1, r_rot)
uc_1 = add_threonine(np.array([11,5,5]), uc_1, 0)
uc_1 = add_leucine(np.array([14,6,5]), uc_1, 1, R.from_quat([0,0,0,1]))
uc_1 = add_threonine(np.array([17,5,5]), uc_1, 0)

uc_2 = add_leucine(np.array([5,5,5]), uc_2, 0, R.from_quat([0,0,0,1])) 
uc_2 = add_leucine(np.array([8,6,5]), uc_2, 1, R.from_quat([0,0,0,1]))
uc_2 = add_threonine(np.array([11,5,5]), uc_2, 0)
uc_2 = add_leucine(np.array([14,6,5]), uc_2, 1, R.from_quat([0,0,0,1]))
uc_2 = add_threonine(np.array([17,5,5]), uc_2, 0)

rmsd_mu = list() 
rmsd_sig = list() 
for i in range(1,max_n_unit): 
    rmsd = list() 
    for j in range(N): 
        xtal = Crystal(n_unit=i, unit_dim=np.array([9,7,7])) 
        xtal.create_unit_cell(uc_1, occ=.5)
        xtal.create_unit_cell(uc_2, occ=.5)

        sf_1 = uc_1.compute_sf_fft() 
        sf_2 = uc_2.compute_sf_fft() 
        sf_xtal = xtal.compute_sf() 

        rmsd.append(compute_sf_rmsd(sf_xtal/i**3, (sf_1 + sf_2)/2))

    rmsd_mu.append(np.mean(rmsd))
    rmsd_sig.append(np.std(rmsd))

    print(i, np.mean(rmsd), np.std(rmsd))


p1 = plt.errorbar(range(1, max_n_unit), rmsd_mu, yerr=rmsd_sig/np.sqrt(N), fmt='bo', capsize=5) 
# p2 = plt.errorbar(range(1, max_n_unit), rmsd_mu[1], yerr=rmsd_sig[1]/np.sqrt(N), fmt='ro', capsize=5) 
# p3 = plt.errorbar(range(1, max_n_unit), rmsd_mu[2], yerr=rmsd_sig[2]/np.sqrt(N), fmt='go', capsize=5) 
plt.title('Fhkl RMSD of disordered xtal & 2 occ weighted ordered xtal (N={})'.format(N))
plt.xlabel('number of unit cells') 
plt.ylabel('Fhkl RMSD')
plt.legend([p1], ['dis vs weight avg'], loc='upper right')

result_dir = '/wynton/home/sali/mhancock/xtal/results/xtal_simulation/plots'
# result_dir = '/Volumes/hard_drive/scratch'

plt.savefig(os.path.join(result_dir, 'fhkl_order_vs_dis_peptide.png'))
