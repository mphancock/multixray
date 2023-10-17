import numpy as np
import os 
import gc
import sys 

def analysis(ft_dir): 
    ft_order_1 = np.load(os.path.join(ft_dir, 'ordered_ft_1.npy'))
    ft_order_2 = np.load(os.path.join(ft_dir, 'ordered_ft_2.npy'))

    ft_mean_order = np.add(ft_order_1, ft_order_2) / 2 

    print('shape: {}'.format(ft_mean_order.shape))

    del ft_order_1 
    del ft_order_2 
    gc.collect()

    ft_disord = np.load(os.path.join(results_path, 'disordered_ft.npy'))

    print('non zero elements: {}'.format(np.count_nonzero(ft_disord)))
    print('non zero elements: {}'.format(np.count_nonzero(ft_mean_order)))

    diff = np.subtract(ft_disord, ft_mean_order)

    del ft_disord 
    del fft_mean_order
    gc.collect()

    rmsd = np.sqrt( np.sum(np.square(np.abs(diff))) / diff.size)
    print('\nrmsd: {}'.format(rmsd))


if __name__ == '__main__': 
    dim = sys.argv[0]
    disorder_analysis(dim)

