import numpy as np


def dft(arr, freq):
    fourier_term = 0
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            for k in range(arr.shape[2]):
                sum = sum + arr[i,j,k] * \
                np.exp(np.array([0-1j]) * 2 * np.pi * np.dot(np.array([i,j,k], freq)))


def dft3D_with_offset(x, x_0, y_0, z_0, x_N, y_N, z_N):
    M, N, P = x.shape
    X = np.zeros((x_N, y_N, z_N),dtype=complex)

    for u in range(x_N):
        for v in range(y_N):
            for w in range(z_N):
                X[u, v, w] = dft3D_with_offset_for_hkl(x, x_0, y_0, z_0, x_N, y_N, z_N, u, v, w)

    return X


def dft3D_with_offset_for_hkl(x, x_0, y_0, z_0, x_N, y_N, z_N, h, k, l):
    M, N, P = x.shape
    X = np.zeros((x_N, y_N, z_N),dtype=complex)

    X = 0
    for m in range(x_0, M+x_0):
        for n in range(y_0, N+y_0):
            for p in range(z_0, P+z_0):
                X += x[m-x_0, n-y_0, p-z_0] * np.exp(-2j * np.pi * (h * m / (x_N) + k * n / (y_N) + l * p / (z_N)))

    return X


def dft3D_with_offset_for_hkl_pool(
        params_dict
):
    x = params_dict["x"]
    x_0, y_0, z_0 = params_dict["x_0"], params_dict["y_0"], params_dict["z_0"]
    x_N, y_N, z_N = params_dict["x_N"], params_dict["y_N"], params_dict["z_N"]
    h, k, l = params_dict["h"], params_dict["k"], params_dict["l"]

    return dft3D_with_offset_for_hkl(x, x_0, y_0, z_0, x_N, y_N, z_N, h, k, l)