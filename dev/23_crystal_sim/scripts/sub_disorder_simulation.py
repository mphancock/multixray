import numpy as np 
from xtal_simulator import Crystal 
import gc 
import sys 
import os 
# import psutil

def run_sim(coords, occ, n_unit, unit_dim): 
    xtal = Crystal()
    xtal.set_n_unit(n_unit)
    xtal.set_unit_dim(unit_dim)

    for i in range(len(coords)): 
        xtal.create_unit_cell(coords[i], occ[i])

    freq = xtal.FT()

    return freq


def analysis(tmp_dir): 
    ft_order_1 = np.load(os.path.join(tmp_dir, 'ordered_1.npy'))
    ft_order_2 = np.load(os.path.join(tmp_dir, 'ordered_2.npy'))

    ft_mean_order = np.add(ft_order_1, ft_order_2) / 2 

    del ft_order_1 
    del ft_order_2 
    gc.collect()

    ft_disord = np.load(os.path.join(tmp_dir, 'disorder.npy'))
    print(ft_disord[10,0,0])
    print(ft_mean_order[10,0,0])
    # print(ft_disord[3,3,3])

    # print('non zero elements: {}'.format(np.count_nonzero(ft_disord)))
    # print('non zero elements: {}'.format(np.count_nonzero(ft_mean_order)))

    diff = np.subtract(ft_disord, ft_mean_order)

    del ft_disord 
    del ft_mean_order
    gc.collect()

    rmsd = np.sqrt( np.sum(np.square(np.abs(diff))) / diff.size)
    # print('rmsd: {}'.format(rmsd))


if __name__ == '__main__': 
    tmp_dir = sys.argv[1] 
    n_unit = int(sys.argv[2])
    unit_dim = np.array([7,7,7])
    n_trial = 10 
    results = np.ndarray(10)
    # result_dir = '/wynton/home/sali/mhancock/xtal/results/xtal_simulation/sub_disorder'
    result_dir = '/Volumes/hard_drive/scratch'

    for i in range(n_trial): 
        if i == 0: 
            freq = run_sim([np.array([[3,3,3]])], [1], n_unit, unit_dim)
            np.save(os.path.join(tmp_dir, 'ordered_1.npy'), freq)
            del freq 
            gc.collect()
        
            freq = run_sim([np.array([[2,3,3]])], [1], n_unit, unit_dim)
            np.save(os.path.join(tmp_dir, 'ordered_2.npy'), freq)
            del freq
            gc.collect()

        freq = run_sim([np.array([[2,3,3]]), np.array([[3,3,3]])], [.5,.5], n_unit, unit_dim)
        np.save(os.path.join(tmp_dir, 'disorder.npy'), freq)
        del freq
        gc.collect()

        rmsd = analysis(tmp_dir)
        results[i] = rmsd

    np.save(os.path.join(result_dir, 'rmsd_{}.npy'.format(n_unit)), results)

    # xtal = Crystal()
    # xtal.set_n_unit(10)
    # xtal.set_unit_dim(np.array([7,7,7]))

    # scatterers = np.array([[3,3,3]])
    # xtal.create_unit_cell(scatterers, .5)

    # scatterers = np.array([[2,3,3]])
    # xtal.create_unit_cell(scatterers, .5)

    # xtal.FT()



# xtal = Crystal()
# xtal.set_n_unit(200)
# xtal.set_unit_dim(np.array([7,7,7]))

# scatterers = np.array([[3,3,3]])
# xtal.create_unit_cell(scatterers, .5)

# scatterers = np.array([[2,3,3]])
# xtal.create_unit_cell(scatterers, .5)

# freq = xtal.FT()
# np.save('../results/xtal_simulation/disordered_ft.npy', freq)

# del xtal 
# del freq 
# gc.collect()

# xtal = Crystal()
# xtal.set_n_unit(200)
# xtal.set_unit_dim(np.array([7,7,7]))

# scatterers = np.array([[3,3,3]])
# xtal.create_unit_cell(scatterers, 1)

# freq = xtal.FT()
# np.save('../results/xtal_simulation/ordered_ft_1.npy', freq)

# del xtal 
# del freq 
# gc.collect()

# xtal = Crystal()
# xtal.set_n_unit(200)
# xtal.set_unit_dim(np.array([7,7,7]))

# scatterers = np.array([[2,3,3]])
# xtal.create_unit_cell(scatterers, 1)

# freq = xtal.FT()
# np.save('../results/xtal_simulation/ordered_ft_2.npy', freq)

# del xtal 
# del freq 
# gc.collect()