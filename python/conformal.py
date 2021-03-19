import numpy as np

def fix_sqrt(z):
    for i in range(len(z)):
        if (z[i].imag < 0):
            z[i] = - z[i]
    return


# The first map, z -> z / (1 - z/b)
def f0(z, b):
    for i in range(len(z)):
        if (z[i] == np.inf):
            z[i] = -b
        elif (b != np.inf):
            if (z[i] == b):
                z[i] = np.inf
            else:
                z[i] /= 1.0 - z[i]/b;
    return


def phi_1(z, z0, z1):
    z[0]  = np.inf
    z[1]  = 0.0
    z[2:] = 1j * np.lib.scimath.sqrt((z[2:] - z1)/(z[2:] - z0))
    fix_sqrt(z[2:])
    return


# The i^th forward map, given a
def phi_i(z, a):
    asq = abs(a)**2
    b = asq / a.real
    c = asq / a.imag

    f0(z, b)

    z[...] = np.lib.scimath.sqrt(z**2 + c**2)
    fix_sqrt(z)

    return


# Assume points are in anti-clockwise order
def phi_np1(z, zeta_np1):
    f0(z[1:], zeta_np1)
    z[1:] = -(z[1:]**2)
    return


def phi_disc(z, a):
    a_conj = a.real - 1j*a.imag
    z[...] = (z - a)/(z - a_conj)
    return




class ConformalMap(object):
    zeta = None
    z    = None
    n    = None

    def __init__(self, z, origin=0. + 1j*0.):
        self.n = len(z) - 1
        self.zeta = np.zeros((self.n+2), dtype=complex)

        # Pack origin onto the end of z so it's mapped 
        # along with the z points
        z = np.concatenate((z, np.array((origin,))))

        z0 = z[0]
        z1 = z[1]
        phi_1(z, z0, z1)

        for i in range(2, self.n+1):
            a = z[i]
            self.zeta[i] = a # So [0], [1] not set?
            phi_i(z, a)

        self.zeta[self.n + 1] = z[0]
        phi_np1(z, self.zeta[self.n+1])

        phi_disc(z, z[-1])

        print(z)
        self.z = z
        return
