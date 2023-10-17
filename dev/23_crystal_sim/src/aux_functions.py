import numpy as np 

## return list of coords 
def get_coords_lists(scatt_list): 
    xs = list() 
    ys = list() 
    zs = list() 

    xs_proj = list() 
    ys_proj = list() 
    zs_proj = list() 

    for scatt in scatt_list: 
        pos = scatt.get_pos_mean() 
        xs.append(pos[0]) 
        ys.append(pos[1]) 
        zs.append(pos[2])

        # if proj: 
        #     xs_proj.append(pos[0]) 
        #     ys_proj.append(pos[1]) 
        #     zs_proj.append(0)

        #     xs_proj.append(0) 
        #     ys_proj.append(pos[1]) 
        #     zs_proj.append(pos[2])

        #     xs_proj.append(pos[0]) 
        #     ys_proj.append(0)
        #     zs_proj.append(pos[2])

    return xs, ys, zs

def get_color(scatt_type): 
    if scatt_type == 'H': 
         return 'w' 
    elif scatt_type == 'O': 
        return 'r' 
    elif scatt_type == 'C': 
        return 'k' 
    elif scatt_type == 'N': 
        return 'b' 
    else: 
        return 'g' 


## compute 'rmsd' between 2 structure factor arrays 
def compute_sf_rmsd(sf_A, sf_B):
    rmsd = np.absolute(sf_A - sf_B)**2
    rmsd = np.sum(rmsd) 
    rmsd = rmsd / (sf_A.shape[0] * sf_A.shape[1] * sf_A.shape[2]) 
    rmsd = rmsd**(1/2)
    return rmsd