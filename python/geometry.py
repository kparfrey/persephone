import matplotlib.pyplot as plt
import numpy as np

pi = np.pi

### Apply the Fourier VMEC mapping
def disc_to_physical_VMEC(r, theta, modes):
    Rm, Zm = modes
    R = 0.0
    Z = 0.0
    M = Rm.shape[0]
    for m in range(1, M):
        R += Rm[m] * np.cos(m*theta)
        Z += Zm[m] * np.sin(m*theta)
    Rfinal = r*R + Rm[0]
    Zfinal = r*Z
    return (Rfinal, Zfinal)


def disc_to_physical_Garabedian(r, theta, modes):
    Deltam, w0 = modes
    M = Deltam.shape[0]

    cost = np.cos(theta)
    sint = np.sin(theta)

    R = 0.0
    Z = 0.0
    for m in range(0, M):
        cosmt = np.cos(m*theta)
        sinmt = np.sin(m*theta)
        R += Deltam[m] * ( cost*cosmt + sint*sinmt ) # using e^{-im0}
        Z += Deltam[m] * (-cost*sinmt + sint*cosmt )
    Rfinal = r*R + w0[0]
    Zfinal = r*Z + w0[1]

    '''
    # Same, but for boundary, r = 1
    Rb = 0.0
    Zb = 0.0
    for m in range(0, M):
        cosmt = np.cos(m*theta)
        sinmt = np.sin(m*theta)
        Rb += Deltam[m] * ( cost*cosmt + sint*sinmt )
        Zb += Deltam[m] * (-cost*sinmt + sint*cosmt )
    Rbfinal = Rb*r + w0[0]
    Zbfinal = Zb*r + w0[1]

    rp = np.sqrt(Rfinal**2 + Zfinal**2)
    rb = np.sqrt(Rbfinal**2 + Zbfinal**2)
    rf = rb/rp

    Rfinal = Rbfinal * rf/rb
    Zfinal = Zbfinal * rf/rb
    '''
    return (Rfinal, Zfinal)


### Given theta, find the corresponding radius coordinate 
### of a square of half-width h
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


def corner_lines(t_outer, h, transform, modes):
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
        Rcl[i,:], Zcl[i,:] = transform(r, theta, modes)
    return Rcl, Zcl



### Unit disc: return a radius=1 array the same length as the theta array
def circle(theta):
    N = theta.shape[0]
    r = np.ones((N,))
    return r



##########################################################################



### For VMEC mapping
M = 5
Rm = np.zeros(M)
Zm = np.zeros(M)

Rm[0] = 2.0
Rm[1] = 1.0
Zm[1] = 2.0
Rm[2] = 0.2
Zm[2] = 0.2
Rm[3] = -0.3
Zm[3] = -0.1


### For Garabedian mapping
# Delta_0 = minor radius, ie radius of the cross section
# Delta_1 = shifts the shape without deformation in the +ve/-ve R direction,
#           without moving the coordinate origin: gives a deformed interior
#           coord system even with a circular outer boundary. Should probably
#           locate the origin in phys space with separate w0, and use this for 
#           control over the created coordinate grid?
# Delta_2 = makes elliptical: +ve = wide, -ve = tall
# Delta_3 = makes triangular: +ve = outerward pointing, -ve = inward pointing
# Delta_4 = makes square
# Delta_5 = makes pentagonal etc. etc.
Deltam = 1.1 * np.array([1.0, 0.0, 0.3, 0.1, 0.05, 0.035, -0.02])
w0 = [3.0, 0.0] # physical-space R,Z location of the polar-coordinate origin


### Choose mapping type

#modes = [Rm, Zm]
#disc_to_physical = disc_to_physical_VMEC

modes = [Deltam, w0]
disc_to_physical = disc_to_physical_Garabedian



h = 0.4  # Inner square half-width on unit disc

Ntheta = 97  # (Multiple of 8) + 1 to place points exactly at pi/4, 3*pi/4 etc
theta = np.linspace(0, 2*np.pi, Ntheta)

Ro, Zo = disc_to_physical(circle(theta),    theta, modes)
Rs, Zs = disc_to_physical(square(theta, h), theta, modes)


### Theta angles of the four corners on the outer boundary
## Square
outer_corners  = np.array([pi/4, 3*pi/4, 5*pi/4, 7*pi/4])

## or:
## Move the two "left" points to improve the angles in physical space
#outer_corners  = np.array([pi/4, 0.85*pi, 1.15*pi, 7*pi/4])


### Outer quad corners in physical space
Rc, Zc = disc_to_physical(circle(outer_corners), outer_corners, modes)

### Lines connecting inner square to outer quad corners, physical space
Rcl, Zcl = corner_lines(outer_corners, h, disc_to_physical, modes)

#origin = disc_to_physical(1e-9, 0.0, modes)
#print("\n")
#print("Mapped origin R: %.5lf" % origin[0])
#print("Mapped origin Z: %.5lf" % origin[1])

plt.figure()


### Plot unit disc without transformation ###
if False:
    Rm0 = Zm0 = np.array([0,1])
    Ro0, Zo0 = disc_to_physical_VMEC(circle(theta),    theta, [Rm0, Zm0])
    Rs0, Zs0 = disc_to_physical_VMEC(square(theta, h), theta, [Rm0, Zm0])
    Rc0, Zc0 = disc_to_physical_VMEC(circle(outer_corners), outer_corners, [Rm0, Zm0])
    Rcl0, Zcl0 = corner_lines(outer_corners, h, disc_to_physical_VMEC, [Rm0, Zm0])


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

ax = plt.gca()
ax.set_aspect('equal')

plt.show()
#plt.savefig('~/Desktop/mapping_simple.pdf')
