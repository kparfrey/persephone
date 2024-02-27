import numpy as np

x1 =  -1.0
x2 =   1.0

eps = 1e-15
n   = 8  # No. of interior points needed = N_soln - 1

m = (n + 1)//2

xm = 0.5 * (x2 + x1)
xl = 0.5 * (x2 - x1)

x = np.ndarray((n))

niterations = 0

for i in range(m):
    z = np.cos(np.pi * (i + 0.75)/(n + 0.5))
    
    while True:
        p1 = 1.0
        p2 = 0.0
        for j in range(n):
            p3 = p2
            p2 = p1
            p1 = ((2.0*j + 1.0) * z * p2 - j * p3) / (j + 1.0)
        pp = n * (z * p1 - p2)/(z*z - 1.0)
        z1 = z
        z = z1 - p1/pp

        niterations += 1
        if (np.abs(z-z1) < eps):
            break

    x[i] = xm - xl * z
    x[n - 1 - i] = xm + xl * z

print(x)
print("No. of iterations: ", niterations)
