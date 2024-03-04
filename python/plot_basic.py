import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
#import fejer_quadrature
#from scipy.interpolate import griddata


### Returns weights for 1D Fejer (type 1) quadrature on npoints points
### Therefore polynomial_order = npoints - 1
def fejer_weights(npoints):
    npoints = int(npoints)
    if npoints == 1:
        return np.array([1.0])

    steps = np.arange(1, npoints, 2)
    length = len(steps)
    remains = npoints - length

    kappa = np.arange(remains)
    beta = np.hstack(
        [
            2 * np.exp(1j * np.pi * kappa / npoints) / (1 - 4 * kappa**2),
            np.zeros(length + 1),
        ]
    )
    beta = beta[:-1] + np.conjugate(beta[:0:-1])

    weights = np.fft.ifft(beta)
    assert max(weights.imag) < 1e-15
    weights = weights.real / 2.0

    return weights 


class Group(object):
    group = None

    r0 = None
    r1 = None
    r2 = None
    r  = None

    Jrdetg = None
    quadweight = None

    e0 = None
    e1 = None
    e2 = None
    e3 = None
    edges = None

    Nelem_tot = None
    Nelem = None
    Ns = None
    Nf = None

    def __init__(self, group_id):
        self.group = group_id

        return
    

    def elem_id_1d(self, i, j, k):
        return (k*self.Nelem[1] + j)*self.Nelem[0] + i



class Mesh(object):
    Ngroup    = None
    Nproc     = None
    Nfield    = None
    Nelem_tot = None
    Nelem     = None
    Ns        = None
    Nf        = None

    weights_found = False

    fileref = None # Handle for h5py file object
    
    g = dict() # Dictionary of Groups

    def __init__(self):
        print("Loading mesh file")
        filename = 'mesh.h5'
        m = h5py.File(filename, 'r')

        self.Ngroup = m['Ngroup'][()]
        self.Nproc  = m['Nproc'][()]
        self.Nfield = m['Nfield'][()]
        
        for ig in range(self.Ngroup):
            sg = str(ig)
            group = Group(ig)

            group.r0 = m[sg]['r']['0']
            group.r1 = m[sg]['r']['1']
            group.r2 = m[sg]['r']['2']
            group.r  = [group.r0, group.r1, group.r2]
            
            group.Jrdetg     = m[sg]['Jrdetg']
            #group.quadweight = m[sg]['quadweight']

            ## This seems to load the data to disk - may want to use list 
            #group.r  = np.array((group.r0,group.r1,group.r2))

            # Edges 0-3, the "front" poloidal edges
            group.e0 = m[sg]['edges']['0']
            group.e1 = m[sg]['edges']['1']
            group.e2 = m[sg]['edges']['2']
            group.e3 = m[sg]['edges']['3']
            # List rather than array since edges don't all have same length
            group.edges = [group.e0, group.e1, group.e2, group.e3]

            group.Nelem_tot = m[sg]['Nelem_tot'][()]
            group.Nelem = m[sg]['Nelem'][:] # No. elems in this group, in each dir
            group.Ns    = m[sg]['Ns'][:]    # Soln points per elem, in each dir
            group.Nf    = m[sg]['Nf'][:]    # Flux    "      "          "

            self.g[ig] = group

        # Don't close the file --- not loading arrays into memory
        # Save handle to the h5py object in case want to close later
        self.fileref = m
        return


    # Draw the "front" poloidal face of a slice of elements 
    def draw_edges(self, nphi = 0, width = 0.2):
        for ig in range(self.Ngroup):
            edges = self.g[ig].edges
            for i in range(self.g[ig].Nelem[0]):
                for j in range(self.g[ig].Nelem[1]):
                    elem = self.g[ig].elem_id_1d(i,j,nphi)
                    for edge_id in range(4):
                        plt.plot(edges[edge_id][elem,0,:], edges[edge_id][elem,1,:], 'k-', lw=width, zorder=1) 

        ax = plt.gca()
        ax.set_aspect('equal')
        return


    # Draw the points immediately inside the "front" poloidal face of the
    # point or element slice indicated by nphi
    def draw_soln_points(self, nphi = 0, size = 5.0, nphi_elemwise=False):
        # By default nphi indexes by soln point. Convert to by element.
        # Assume Ns[2] the same for all groups
        if nphi_elemwise:
            nphi *= self.g[0].Ns[2]
            
        for i in range(self.Ngroup):
            r0 = self.g[i].r0
            r1 = self.g[i].r1

            # Add c='k' to draw the points in black
            plt.scatter(r0[:,:,nphi], r1[:,:,nphi], marker='.', s=size)

        ax = plt.gca()
        ax.set_aspect('equal')
        return
    

    ### Find weights for 3D Fejer quadrature
    ### For now, assume all elements in a group are identical
    def setup_quadrature_weights(self):
        print("Finding quadrature weights")

        for ig in range(self.Ngroup):
            g  = self.g[ig]
            Nelem = g.Nelem
            Ns = g.Ns

            g.quadweight = np.ndarray(g.Jrdetg.shape)

            w0 = fejer_weights(Ns[0])
            w1 = fejer_weights(Ns[1])
            w2 = fejer_weights(Ns[2])

            # Not the most efficient, but at least it's clear...
            for ke in range(Nelem[2]):
                for je in range(Nelem[1]):
                    for ie in range(Nelem[0]):
                        for k in range(Ns[2]):
                            for j in range(Ns[1]):
                                for i in range(Ns[0]):
                                    il = ie * Ns[0] + i
                                    jl = je * Ns[1] + j
                                    kl = ke * Ns[2] + k
                                    g.quadweight[il][jl][kl] = w0[i] * w1[j] * w2[k]

        self.weights_found = True
        return







class Snapshot(object):
    m       = None # Mesh object
    dfile   = None # H5 file for data
    filenum = None # Number of loaded data file
    time    = None
    step    = None
    dt      = None

    Ngroup  = None

    # Use these to convert B to e.g. Teslas
    # The code writes B / sqrt{mu_0}
    mu0     = 4 * np.pi * 1e-7
    sqrt_mu0 = np.sqrt(mu0) 

    adiabatic_idx = 5.0 / 3.0;

    def __init__(self, filenum = 0):
        self.m = Mesh()
        self.Ngroup = self.m.Ngroup

        if (filenum >= 0):
            self.load_data(filenum)

        return


    def load_data(self, filenum = 0):
        print("Loading data from file number %d" % filenum)
        self.filenum = filenum

        if self.dfile is not None:
            self.dfile.close()

        filename = 'data%04d.h5' % filenum
        self.dfile = h5py.File(filename,'r')

        self.rho = [0] * self.Ngroup 
        self.p   = [0] * self.Ngroup 
        self.v0  = [0] * self.Ngroup 
        self.v1  = [0] * self.Ngroup 
        self.v2  = [0] * self.Ngroup 
        #self.v   = [0] * self.Ngroup
        self.B0  = [0] * self.Ngroup
        self.B1  = [0] * self.Ngroup 
        self.B2  = [0] * self.Ngroup 
        self.psi = [0] * self.Ngroup
        self.divB = [0] * self.Ngroup
        #self.B   = [0] * self.Ngroup

        self.vsq = [0] * self.Ngroup
        self.Bsq = [0] * self.Ngroup
        self.divB_B = [0] * self.Ngroup # |divB/B|
        self.beta = [0] * self.Ngroup
        self.over_beta = [0] * self.Ngroup
        self.ptot = [0] * self.Ngroup
        self.cs   = [0] * self.Ngroup   # sound speed
        self.mach = [0] * self.Ngroup
        
        self.R   = [0] * self.Ngroup
        self.Z   = [0] * self.Ngroup
        self.Phi = [0] * self.Ngroup

        for ig in range(self.Ngroup):
            sg = str(ig)
            self.rho[ig] = self.dfile[sg]['rho']
            self.p[ig]   = self.dfile[sg]['p']
            self.v0[ig]  = self.dfile[sg]['v0']
            self.v1[ig]  = self.dfile[sg]['v1']
            self.v2[ig]  = self.dfile[sg]['v2']
            #self.v[ig]   = [self.v0, self.v1, self.v2]
            self.B0[ig]  = self.dfile[sg]['B0']
            self.B1[ig]  = self.dfile[sg]['B1']
            self.B2[ig]  = self.dfile[sg]['B2']
            #self.B[ig]   = [self.B0, self.B1, self.B2]
            self.psi[ig] = self.dfile[sg]['psi']
            self.divB[ig] = self.dfile[sg]['divB']

            self.R[ig] = self.m.g[ig].r0
            self.Z[ig] = self.m.g[ig].r1
            self.Phi[ig] = self.m.g[ig].r2

            # These actually load data into memory - might want to reorganise
            self.vsq[ig] = self.v0[ig][...]**2 + self.v1[ig][...]**2 + self.v2[ig][...]**2
            self.Bsq[ig] = self.B0[ig][...]**2 + self.B1[ig][...]**2 + self.B2[ig][...]**2
            self.divB_B[ig] = np.abs(self.divB[ig][...])/(np.sqrt(self.Bsq[ig]) + 1e-15)
            self.beta[ig] = np.abs(self.p[ig]) / (0.5 * self.Bsq[ig][...] + 1e-15)
            self.over_beta[ig] = 1.0 / self.beta[ig][...]
            self.ptot[ig] = self.p[ig] + 0.5*self.Bsq[ig]
            self.cs[ig] = np.sqrt(self.adiabatic_idx * self.p[ig][...] / self.rho[ig][...])
            self.mach[ig] = np.sqrt(self.vsq[ig])/self.cs[ig]


        # Can't make "global" v array, since the groups can have
        # different numbers of points, elements etc.
        #self.v   = np.array((self.v0, self.v1, self.v2))

        self.time = self.dfile['time'][()]
        self.step = self.dfile['step'][()]
        self.dt   = self.dfile['dt'][()]

        # Don't close the file --- not loading arrays into memory
        return


    # k is the grid index in the azimuthal direction 
    def contour_plot(self, variable, width=0.7, levels=[None,], k=0, mag=False):
        if levels[0] == None:
            if mag:
                fmax = np.amax(np.abs(variable))
                fmin = np.amin(np.abs(variable))
            else:
                fmax = np.amax(variable)
                fmin = np.amin(variable)
            levels = np.linspace(fmin, fmax, 30)

        R = [0] * self.Ngroup
        Z = [0] * self.Ngroup
        f = [0] * self.Ngroup

        for ig in range(self.Ngroup):
            R[ig] = self.R[ig][:,:,k]
            Z[ig] = self.Z[ig][:,:,k]
            f[ig] = variable[ig][:,:,k]
            if mag:
                f[ig] = np.abs(f[ig])

        if self.Ngroup > 1:
            R = self.connect_groups(R)
            Z = self.connect_groups(Z)
            f = self.connect_groups(f)
            
        for ig in range(self.Ngroup):
            plt.contour(R[ig], Z[ig], f[ig], levels=levels, linewidths=width, zorder=5)

        plt.colorbar()
        plt.title('t = %.4lf' % self.time)

        ax = plt.gca()
        ax.set_aspect('equal')
        return


    # Add extra lines of fake solution points along the group boundaries, to avoid
    # white lines in the composite plot.
    def connect_groups(self, f):
        (N0cen, N1cen) = f[0].shape
        (N0out, N1out) = f[1].shape

        F = [0] * self.Ngroup
        F[0] = np.ndarray((N0cen+2, N1cen+2))
        for i in [1,2,3,4]:
            F[i] = np.ndarray((N0out+1, N1out+2))

        F[0][1:-1,1:-1] = f[0]
        for i in [1,2,3,4]:
            F[i][1:,1:-1] = f[i] 

        # Radial lines between outer groups
        F[1][1:,0] = F[4][1:,-1] = 0.5 * (f[1][:,0] + f[4][:,-1])
        F[2][1:,0] = F[1][1:,-1] = 0.5 * (f[2][:,0] + f[1][:,-1])
        F[3][1:,0] = F[2][1:,-1] = 0.5 * (f[3][:,0] + f[2][:,-1])
        F[4][1:,0] = F[3][1:,-1] = 0.5 * (f[4][:,0] + f[3][:,-1])

        # Lines between the central group and an outer group
        F[0][1:-1,0]  = F[1][0,1:-1]    = 0.5 * (f[0][:,0] + f[1][0,:])
        F[0][1:-1,-1] = F[3][0,1:-1][::-1] = 0.5 * (f[0][:,-1] + f[3][0,::-1])
        F[0][0,1:-1]  = F[4][0,1:-1][::-1] = 0.5 * (f[0][0,:]  + f[4][0,::-1])
        F[0][-1,1:-1] = F[2][0,1:-1]    = 0.5 * (f[0][-1,:] + f[2][0,:])

        # The four corner at which the central group meets two outer groups
        third = 1.0/3.0
        F[0][0,0] = F[1][0,0] = F[4][0,-1] = third * (f[0][0,0] + f[1][0,0] + f[4][0,-1])
        F[0][-1,0] = F[1][0,-1] = F[2][0,0] = third * (f[0][-1,0] + f[1][0,-1] + f[2][0,0])
        F[0][-1,-1] = F[2][0,-1] = F[3][0,0] = third * (f[0][-1,-1] + f[2][0,-1] + f[3][0,0])
        F[0][0,-1] = F[3][0,-1] = F[4][0,0] = third * (f[0][0,-1] + f[3][0,-1] + f[4][0,0])

        return F


    def global_integral(self, f):
        sum_all_groups = 0.0

        if self.m.weights_found == False:
            self.m.setup_quadrature_weights()

        for ig in range(self.Ngroup):
            g = self.m.g[ig]
            sum_all_groups += np.sum( f[ig][()] * g.Jrdetg[()] * g.quadweight[()] , axis=None )

        return sum_all_groups


