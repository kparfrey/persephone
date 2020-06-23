import matplotlib.pyplot as plt
import numpy as np

pi = np.pi

def disc_to_physical(r, theta, Rm, Zm):
    R = 0.0
    Z = 0.0
    M = Rm.shape[0]
    for m in range(1, M):
        R += Rm[m] * np.cos(m*theta)
        Z += Zm[m] * np.sin(m*theta)
    Rfinal = r*R + Rm[0]
    Zfinal = r*Z
    return (Rfinal, Zfinal)


def square(theta, h):
    N = theta.shape[0]
    r = np.ndarray((N,))
    for n in range(N):
        t = theta[n]
        if (t < pi/4) or (t > 3*pi/4 and t < 5*pi/4) or (t > 7*pi/4):
            r[n] = h / np.abs(np.cos(t))
        else:
            r[n] = h / np.abs(np.sin(t))
    return r


def corner_lines(t_outer, h, Rm, Zm):
    Nr = 70
    R0s = np.array([h, -h, -h,  h])
    Z0s = np.array([h,  h, -h, -h])
    theta_add = np.array([0, pi, pi, 0])
    Rcl = np.ndarray((4,Nr)) # These are in physical space
    Zcl = np.ndarray((4,Nr))
    for i in range(4):
        R0 = R0s[i]  # Cyl coords on unit disc for inner (square) corners
        Z0 = Z0s[i]
        R1 = np.cos(t_outer[i]) # Cyl coords on unit disc for outer corners
        Z1 = np.sin(t_outer[i])
        slope = (Z1-Z0)/(R1-R0)
        R = np.linspace(R0, R1, Nr)
        Z = Z0 + slope*(R - R0)
        r = np.sqrt(R*R + Z*Z)
        theta = np.arctan(Z/R) + theta_add[i]
        Rcl[i,:], Zcl[i,:] = disc_to_physical(r, theta, Rm, Zm)
    return Rcl, Zcl



def circle(theta):
    N = theta.shape[0]
    r = np.ones((N,))
    return r



##########################################################################



M = 5
Rm = np.zeros(M)
Zm = np.zeros(M)

Rm[0] = 2.0
Rm[1] = 1.0
Zm[1] = 2.0
Rm[2] = 0.2
Zm[2] = 0.2
Rm[3] = 0.0
Zm[3] = 0.0

h = 0.4  # Inner square half-width on unit disc

Ntheta = 97  # (Multiple of 8) + 1 to place points exactly at pi/4, 3*pi/4 etc
theta = np.linspace(0, 2*np.pi, Ntheta)

Ro, Zo = disc_to_physical(circle(theta),    theta, Rm, Zm)
Rs, Zs = disc_to_physical(square(theta, h), theta, Rm, Zm)

### Quad corners on the outer boundary
outer_corners  = np.array([pi/4, 3*pi/4, 5*pi/4, 7*pi/4]) ### Default: same as square
#outer_corners  = np.array([pi/4, 0.85*pi, 1.15*pi, 7*pi/4])

### Outer quad corners in physical space
Rc, Zc = disc_to_physical(circle(outer_corners), outer_corners, Rm, Zm)

### Lines connecting inner square to outer quad corners, physical space
Rcl, Zcl = corner_lines(outer_corners, h, Rm, Zm)



plt.figure()
ax = plt.gca()
ax.set_aspect('equal')

### Plot unit disc without transformation ###
Rm0 = Zm0 = np.array([0,1])
Ro0, Zo0 = disc_to_physical(circle(theta),    theta, Rm0, Zm0)
Rs0, Zs0 = disc_to_physical(square(theta, h), theta, Rm0, Zm0)
Rc0, Zc0 = disc_to_physical(circle(outer_corners), outer_corners, Rm0, Zm0)
Rcl0, Zcl0 = corner_lines(outer_corners, h, Rm0, Zm0)
plt.plot(Ro0, Zo0, 'C0')
plt.plot(Rs0, Zs0, 'C1')
for i in [0,1,2,3]:
    plt.plot(Rcl0[i,:], Zcl0[i,:], 'C2')
plt.plot(Rc0, Zc0, '.', color='C3')

### Plot transformed cross-section ##########
plt.plot(Ro, Zo, 'C0')
plt.plot(Rs, Zs, 'C1')

for i in [0,1,2,3]:
    plt.plot(Rcl[i,:], Zcl[i,:], 'C2')
plt.plot(Rc, Zc, '.', color='C3')

plt.show()
#plt.savefig('/Users/kyle/Desktop/mapping_simple.pdf')