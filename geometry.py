import matplotlib.pyplot as plt
import numpy as np

pi = np.pi

def disc_to_physical(r, theta, Rm, Zm):
    R = 0.0
    Z = 0.0
    M = Rm.shape[0]
    for m in range(M):
        R += Rm[m] * np.cos(m*theta)
        Z += Zm[m] * np.sin(m*theta)
    return (r*R,r*Z)


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


def corner_line(theta, h):
    N = 70


def circle(theta):
    N = theta.shape[0]
    r = np.ones((N,))
    return r


M = 5
Rm = np.zeros(M)
Zm = np.zeros(M)

Rm[1] = 1.0
Zm[1] = 2.0
Rm[2] = 0.4
Zm[2] = 0.2
#Rm[3] = -0.1
Zm[3] = 0.0

h = 0.4  # Inner square half-width on unit disc

Ntheta = 97  # (Multiple of 8) + 1 to place points exactly at pi/4, 3*pi/4 etc
theta = np.linspace(0, 2*np.pi, Ntheta)

Ro, Zo = disc_to_physical(circle(theta),      theta, Rm, Zm)
Rs, Zs = disc_to_physical(square(theta, h), theta, Rm, Zm)

### Quad corners on the outer boundary
### Simplest case --- all have the same pi/4-like angles as the square's corners
theta_corners = np.array([pi/4, 3*pi/4, 5*pi/4, 7*pi/4])
Rc, Zc = disc_to_physical(circle(theta_corners), theta_corners, Rm, Zm)

### Corner lines
### Simplest case --- will want to generalize to that all four corners on the disc
### are arbitrary, and in general each is an independent function of toroidal angle
Nr = 70
r0 = np.sqrt(2*h*h)
rcl = np.linspace(r0, 1.0, Nr)
unity = np.ones((Nr,))
Rcl = np.ndarray((4,Nr))
Zcl = np.ndarray((4,Nr))
for i in [0,1,2,3]:
    Rcl[i,:], Zcl[i,:] = disc_to_physical(rcl, theta_corners[i]*unity, Rm, Zm) 


plt.figure()
ax = plt.gca()
ax.set_aspect('equal')

plt.plot(Ro, Zo)
plt.plot(Rs, Zs)

for i in [0,1,2,3]:
    plt.plot(Rcl[i,:], Zcl[i,:], 'k')
#plt.plot(Rc, Zc, 'o', color='r')

plt.show()
