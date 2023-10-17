from xtal_simulator import Scatterer, Unit_Cell, Crystal 

from scipy.spatial.transform import Rotation as R 
import matplotlib.pyplot as plt 
import numpy as np 

def add_amino_acid_base(loc, uc, r, pep_bond_l=True, pep_bond_r=True):
    c_alpha = Scatterer(r.apply(np.array([0,0,0]))+loc, 'C') 
    c_alpha_h = Scatterer(r.apply(np.array([0,0,-1]))+loc, 'H')
    cooh_c = Scatterer(r.apply(np.array([-1,-1,0]))+loc, 'C') 
    cooh_o_1 = Scatterer(r.apply(np.array([-2,0,0]))+loc, 'O')
    cooh_o_2 = Scatterer(r.apply(np.array([-1,-2,0]))+loc, 'O')
    cooh_h = Scatterer(r.apply(np.array([-1,-3,0]))+loc, 'H')
    amino_n = Scatterer(r.apply(np.array([1,-1,0]))+loc, 'N') 
    amino_h_1 = Scatterer(r.apply(np.array([2,-1,0]))+loc, 'H')
    amino_h_2 = Scatterer(r.apply(np.array([1,-2,0]))+loc, 'H')
    
    amino_acid_base = [c_alpha, c_alpha_h, amino_n, cooh_c,\
                       cooh_o_2, cooh_h]
    if not pep_bond_l: 
        amino_acid_base.extend([cooh_o_1])
    if not pep_bond_r: 
        amino_acid_base.extend([amino_h_1, amino_h_2])


    uc.add_scatterers(amino_acid_base) 
    return uc 

def add_threonine(loc, uc, order): 
    if order == 0: 
        r_base = R.from_quat([-np.sin(np.pi/2),0,0,np.cos(np.pi/2)])
        r = R.from_quat([-np.sin(np.pi/4),0,0,np.cos(np.pi/4)])
    else: 
        r_base = R.from_quat([0,0,0,1])
        r = R.from_quat([np.sin(np.pi/4),0,0,np.cos(np.pi/4)])
 
    uc = add_amino_acid_base(loc, uc, r_base)

    c_1 = Scatterer(r.apply(np.array([0,1,0]))+loc, 'C') 
    c_2 = Scatterer(r.apply(np.array([1,2,0]))+loc, 'C') 
    o_1 = Scatterer(r.apply(np.array([-1,2,0]))+loc, 'O') 
    c_2_h_1 = Scatterer(r.apply(np.array([1,2,-1]))+loc, 'H') 
    c_2_h_2 = Scatterer(r.apply(np.array([1,3,0]))+loc, 'H')
    c_2_h_3 = Scatterer(r.apply(np.array([1,2,1]))+loc, 'H')
    o_1_h_1 = Scatterer(r.apply(np.array([-1,3,0]))+loc, 'H') 

    uc.add_scatterers([c_1, c_2, o_1, c_2_h_1, c_2_h_2, c_2_h_3, o_1_h_1])
    # uc = add_amino_acid_base(loc, uc)

    return uc 


## rotate amino acid around c alpha with r quaternion 
## rotate amino acid around c2 with r_rot (represents rotameric state)
def add_leucine(loc, uc, order, r_rot): 
    if order == 0: 
        r_base = R.from_quat([-np.sin(np.pi/2),0,0,np.cos(np.pi/2)])
        r = R.from_quat([0,0,-np.sin(np.pi/4),np.cos(np.pi/4)])
        r = r * R.from_quat([-np.sin(np.pi/4),0,0,np.cos(np.pi/4)])
    else: 
        r_base = R.from_quat([0,0,0,1])
        r = R.from_quat([0,0,np.sin(np.pi/4),np.cos(np.pi/4)])
        r = r * R.from_quat([np.sin(np.pi/4),0,0,np.cos(np.pi/4)]) 

    uc = add_amino_acid_base(loc, uc, r_base)

    c_1 = r.apply(np.array([0,1,0]))+loc
    c_2 = r.apply(np.array([1,2,0]))+loc
    c_3 = r_rot.apply((r.apply(np.array([2,1,0]))+loc - c_2)) + c_2
    c_4 = r_rot.apply(r.apply(np.array([1,3,0]))+loc - c_2) + c_2

    scatt = list() 
    for c in [c_1, c_2, c_3, c_4]: 
        scatt.append(Scatterer(c, 'C'))

    c_1_h_1 = r.apply(np.array([0,1,1]))+loc 
    c_1_h_2 = r.apply(np.array([0,1,-1]))+loc 
    c_2_h_1 = r.apply(np.array([1,2,1]))+loc 

    c_3_h_1 = r_rot.apply(r.apply(np.array([3,1,0]))+loc - c_2) + c_2
    c_3_h_2 = r_rot.apply(r.apply(np.array([2,1,1]))+loc - c_2) + c_2 
    c_3_h_3 = r_rot.apply(r.apply(np.array([2,1,-1]))+loc - c_2) + c_2 
    c_4_h_1 = r_rot.apply(r.apply(np.array([1,4,0]))+loc - c_2) + c_2 
    c_4_h_2 = r_rot.apply(r.apply(np.array([0,3,0]))+loc - c_2) + c_2 
    c_4_h_3 = r_rot.apply(r.apply(np.array([1,3,1]))+loc - c_2) + c_2  

    for c_h in [c_1_h_1,c_1_h_2,c_2_h_1,c_3_h_1,c_3_h_2,c_3_h_3,c_4_h_1,c_4_h_2,c_4_h_3]: 
        scatt.append(Scatterer(c_h, 'H'))

    uc.add_scatterers(scatt)

    return uc


# uc = Unit_Cell() 
# uc.set_dim(np.array([20,10,10]))
# r_rot = R.from_quat([0, np.sin(np.pi/4), 0, np.cos(np.pi/4)])

# uc = add_leucine(np.array([5,5,5]), uc, 0, R.from_quat([0,0,0,1])) 
# uc = add_leucine(np.array([8,6,5]), uc, 1, r_rot)
# uc = add_threonine(np.array([11,5,5]), uc, 0)
# uc = add_leucine(np.array([14,6,5]), uc, 1, R.from_quat([0,0,0,1]))
# uc = add_threonine(np.array([17,5,5]), uc, 0)
# # uc = add_leucine(np.array([8,6,5]), uc, r, r_rot) 

# # uc = add_threonine(np.array([5,5,5]), uc, [-np.sin(np.pi/4),0,0,np.cos(np.pi/4)])
# # uc = add_leucine(np.array([5,5,5]), uc, [-np.sin(np.pi/4),0,0,np.cos(np.pi/4)])
# fig = uc.show_graph(True, color=True, bonds=True)

# plt.show()