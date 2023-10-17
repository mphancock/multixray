## numeric integration using midpoint quadrature algorithm 
## integrates 1 dimensional function f between a and b with n_pts equally spaced quadrature points 
def midpoint_quadrature(f, n_pts, a, b): 
    d = (b-a)/n_pts 
    # print(d)
    F = 0  

    for i in range(n_pts): 
        x = a + i * d + d / 2 
        F = F + f(x) * d
    
    return F  

## numeric integration using midpoint quadrature algorithm 
## integrates 3 dimensional function f between a and b with n_pts equally spaced quadrature points 
def midpoint_quadrature_3(f, n_pts, a, b):
    d = ((b[0]-a[0])/n_pts[0], (b[1]-a[1])/n_pts[1], (b[2]-a[2])/n_pts[2])
    F = 0 

    for i in range(n_pts[0]): 
        x = a[0] + i * d[0] + d[0] / 2
        for j in range(n_pts[1]): 
            y = a[1] + j * d[1] + d[1] / 2
            for k in range(n_pts[2]): 
                z = a[2] + k * d[2] + d[2] / 2

                F = F + f(x, y, z) * d[0] * d[1] * d[2]
                 
    return F 