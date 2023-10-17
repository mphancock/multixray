import numpy as np 
from math import radians

def is_within_ellipsoid(point, center, radius):
    sqr_sum = np.sum(np.square((point-center) / radius))
    if np.sqrt(sqr_sum) > 1: 
        return False  
    else: 
        return True 


def correct_bounds(center, dim):
    if center[0] < 0: 
        center[0] = 0 
    if center[1] < 0: 
        center[1] = 0
    if center[2] < 0: 
        center[2] = 0 
    if center[0] >= dim[0]: 
        center[0] = dim[0]-1 
    if center[1] >= dim[1]: 
        center[1] = dim[1]-1 
    if center[2] >= dim[2]: 
        center[2] = dim[2]-1

    return center 


## take a 3-D numpy array theta that describes rotation in degrees 
def rotation_matrix(theta): 
    theta = np.array([radians(theta[0]), radians(theta[1]), radians(theta[2])])
    Rx = np.zeros([3,3])
    Rx[0,0] = 1 
    Rx[1,1] = np.cos(theta[0]) 
    Rx[1,2] = -np.sin(theta[0]) 
    Rx[2,1] = np.sin(theta[0]) 
    Rx[2,2] = np.cos(theta[0])

    Ry = np.zeros([3,3])
    Ry[0,0] = np.cos(theta[1])
    Ry[0,2] = np.sin(theta[1])
    Ry[1,1] = 1 
    Ry[2,0] = -np.sin(theta[1])
    Ry[2,2] = np.cos(theta[1])

    Rz = np.zeros([3,3]) 
    Rz[0,0] = np.cos(theta[2])
    Rz[0,1] = -np.sin(theta[2])
    Rz[1,0] = np.sin(theta[2])
    Rz[1,1] = np.cos(theta[2])
    Rz[2,2] = 1 

    rot_mat = np.matmul(Rz, np.matmul(Ry, Rx))
    return rot_mat 


def rotate_vector(vec, theta): 
    rot_mat = rotation_matrix(theta) 
    vec_rot = np.matmul(rot_mat, vec) 

    return vec_rot
    
