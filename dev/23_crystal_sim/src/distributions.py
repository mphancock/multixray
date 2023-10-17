import numpy as np 
from integration import midpoint_quadrature

## probability density function of N dimensional Gaussian at X 
## pdf is defined by N dimensional mean, mu, and N,N dimensional covariance matrix, sigma
def gaussian(X, mu, sigma):
    dim = X.shape[0] 
    if dim != mu.shape[0] or dim != sigma.shape[0]: 
        return -1 

    if dim == 1: 
        f = np.exp(-(1/2) * np.dot(np.dot((X-mu), 1/sigma), (X-mu)))
        f = f / ((2 * np.pi) ** (dim/2) * sigma[0]**(1/2))
    else: 
        f = np.exp(-(1/2) * np.dot(np.dot((X-mu), np.linalg.inv(sigma)), (X-mu)))
        f = f / ((2 * np.pi) ** (dim/2) * np.linalg.det(sigma)**(1/2))
 
    return f 


## cdf calculator for univariate normal distribution 
## take sigma (standard deviation)
## gaussian function requires sigma matrix (covariance)
def gaussian_cdf(mu, sigma, a, b): 
    cdf = midpoint_quadrature(lambda x: gaussian(np.array([x]), np.array([mu]), np.array([sigma**2])), 1000, a, b)    
    return cdf 


## return minimum sigma of Gaussian that satisfies following inequality 
## F(a) - F(-a) > min  
## a is a N-D vector 
def find_min_std(a, min_val): 
    sig_vec = np.ndarray(a.shape)
    for i in range(a.shape[0]):
        cdf = .5 
        sigma = 0
        while cdf > min_val/2: 
            sigma = sigma + .1 
            ## gaussian_cdf takes mu, sigma (std dev), 0, and bound 
            cdf = gaussian_cdf(0, sigma, 0, a[i])

        sig_vec[i] = sigma - .1 

    # print('sig_vec', sig_vec)
    return sig_vec  

