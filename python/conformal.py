import numpy as np

def fix_sqrt(z):
    for i in range(len(z)):
        if (z[i].imag < 0):
            z[i] = - z[i]
    return


def sqrt(z):
    w = np.lib.scimath.sqrt(z)
    fix_sqrt(w)
    return w


def sqrt_scalar(z):
    w = np.lib.scimath.sqrt(z)
    if w.imag < 0:
        w = - w
    return w



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


def f0_scalar(z, b):
    if (z == np.inf):
        z = -b
    elif (b != np.inf):
        if (z == b):
            z = np.inf
        else:
            z /= 1.0 - z/b;
    return z


def phi_1(z, z0, z1):
    z[0]  = np.inf
    z[1]  = 0.0
    z[2:] = 1j * np.lib.scimath.sqrt((z[2:] - z1)/(z[2:] - z0))
    fix_sqrt(z[2:])
    return


def phi_inv_1(z, z0, z1):
    zsq = z**2
    z[...] = (z1 + zsq*z0)/(zsq + 1.0)
    return


# The i^th forward map, given a
def phi_i(z, a):
    asq = abs(a)**2
    b = asq / a.real
    c = asq / a.imag

    f0(z, b)

    for i in range(len(z)):
        if z[i] == 0:
            z[i] = c
        else:
            z[i] *= np.lib.scimath.sqrt(1.0 + c**2/z[i]**2)

    ### The above gives much lower rounding error than
    ### than doing this directly:
    #z[...] = np.lib.scimath.sqrt(z**2 + c**2)
    
    fix_sqrt(z)

    return


def phi_inv_i(z, a):
    asq = abs(a)**2
    b = asq / a.real
    c = asq / a.imag

    z *= np.lib.scimath.sqrt(1.0 - c**2 / z**2)
    fix_sqrt(z)
    f0(z, -b)

    return


# Assume points are in anti-clockwise order
def phi_np1(z, zeta_np1):
    f0(z[1:], zeta_np1)
    z[1:] = -(z[1:]**2)
    return


def phi_inv_np1(z, zeta_np1):
    for i in range(len(z)):
        if z[i] == np.inf:
            z[i] = zeta_np1
        else:
            z[i] = 1j * np.lib.scimath.sqrt(z[i])
            if z[i].imag < 0:
                z[i] *= -1
            z[i] = f0_scalar(z[i], -zeta_np1)
    return


def phi_disc(z, a):
    a_conj = a.real - 1j*a.imag

    for i in range(len(z)):
        if z[i] == np.inf:
            z[i] =  1.0 + 0*1j
        elif z[i] == a_conj:
            z[i] = np.inf
        else:
            z[i] = (z[i] - a)/(z[i] - a_conj)

    return


def phi_inv_disc(z, a):
    a_conj = a.real - 1j*a.imag

    for i in range(len(z)):
        if z[i] == 1.0 + 0*1j:
            z[i] = np.inf
        else:
            z[i] = (z[i]*a_conj - a)/(z[i] - 1)
    return



def test_shape(N):
    Am = [0.0, 2.0, 0.2, 0.35, 0.1]
    M = len(Am)
    z = np.zeros((N,), dtype=np.clongdouble)
    for n in range(N):
        for m in range(M):
            z[n] += Am[m] * np.exp(1j*n*m*2*np.pi/N)

    return z


#T: no. of theta sectors, R: no. of radial sectors
def uds_grid(T=14, R=6):
    z = np.zeros((T*R,), dtype=np.clongdouble)
    for t in range(T):
        for r in range(1,R+1):
            z[t*R + r-1] += (r/R) * np.exp(1j*t*2*np.pi/T)

    return z


class ConformalMap(object):
    zeta = None
    origin_mapped = None
    z    = None
    n    = None

    def __init__(self, z, origin=0):
        self.n = len(z) - 1
        self.zeta = np.zeros((self.n+2), dtype=np.clongdouble)

        # Pack origin onto the end of z so it's mapped 
        # along with the z points
        z = np.concatenate((z, np.array((origin,))))

        z0 = z[0]
        z1 = z[1]

        self.zeta[0] = z0
        self.zeta[1] = z1
        phi_1(z, z0, z1)
        for i in range(2, self.n+1):
            a = z[i]
            self.zeta[i] = a
            phi_i(z, a)

        self.zeta[self.n + 1] = z[0]
        phi_np1(z, self.zeta[self.n+1])

        self.origin_mapped = z[-1]
        phi_disc(z, z[-1])

        self.z = z
        return


    def inverse(self, z_in):
        z = np.copy(z_in)

        phi_inv_disc(z, self.origin_mapped)
        phi_inv_np1(z, self.zeta[self.n+1])

        for i in reversed(range(2, self.n+1)):
            phi_inv_i(z, self.zeta[i])

        phi_inv_1(z, self.zeta[0], self.zeta[1])

        return z

