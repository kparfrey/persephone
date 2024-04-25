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


#### Convenience functions for basic operations on lists of arrays #####
#### d stands for "domain"

### Global minimum and maximum
def dmin(var):
    return np.amin([np.amin(x) for x in var])

def dmax(var):
    return np.amax([np.amax(x) for x in var])

def dabs(var):
    return [np.abs(x) for x in var]

def dsqrt(var):
    return [np.sqrt(x) for x in var]

def dsq(var):
    return [np.square(x) for x in var]

def dmult(a, b):
    return [np.multiply(x,y) for (x,y) in zip(a,b)]
########################################################################

class Group(object):

    def __init__(self, group_id):
        self.group = group_id
        self.r0 = None
        self.r1 = None
        self.r2 = None
        self.r  = None

        self.Jrdetg = None
        self.quadweight = None

        self.e0 = None
        self.e1 = None
        self.e2 = None
        self.e3 = None
        self.edges = None

        self.Nelem_tot = None
        self.Nelem = None
        self.Ns = None
        self.Nf = None

        return
    

    def elem_id_1d(self, i, j, k):
        return (k*self.Nelem[1] + j)*self.Nelem[0] + i



class Mesh(object):

    def __init__(self):
        print("Loading mesh file")
        filename = 'data/mesh.h5'
        m = h5py.File(filename, 'r')

        self.Ngroup = m['Ngroup'][()]
        self.Nproc  = m['Nproc'][()]
        self.Nfield = m['Nfield'][()]

        self.weights_found = False
        
        ### Total numbers for the whole mesh
        self.Nelem = 0
        self.Ns    = 0
        self.Nf    = 0

        self.g = dict() # Dictionary of Groups

        for ig in range(self.Ngroup):
            sg = str(ig)
            group = Group(ig)

            group.r0 = m[sg]['r']['0']
            group.r1 = m[sg]['r']['1']
            group.r2 = m[sg]['r']['2']
            group.r  = [group.r0, group.r1, group.r2]
            
            group.Jrdetg     = m[sg]['Jrdetg']

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

            self.Nelem += group.Nelem_tot
            self.Ns    += group.Ns[0] * group.Ns[1] * group.Ns[2]
            self.Nf    += group.Nf[0] * group.Nf[1] * group.Nf[2]

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

        try:
            self.dfile.close()
        except:
            pass

        filename = 'data/data%04d.h5' % filenum
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
        self.Etot = [0] * self.Ngroup
        
        self.R   = [0] * self.Ngroup
        self.Z   = [0] * self.Ngroup
        self.Phi = [0] * self.Ngroup

        self.PoloidalFlux = [0] * self.Ngroup
        self.MaxPFOuterSurface = None

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
            self.Etot[ig] = ( 0.5 * self.rho[ig][...] * self.vsq[ig] +
                              0.5 * self.Bsq[ig] +
                              self.p[ig][...] / (self.adiabatic_idx - 1.0) + 
                              0.5 * self.psi[ig][...] * self.psi[ig][...] )


        # Can't make "global" v array, since the groups can have
        # different numbers of points, elements etc.
        #self.v   = np.array((self.v0, self.v1, self.v2))

        self.time = self.dfile['time'][()]
        self.step = self.dfile['step'][()]
        self.dt   = self.dfile['dt'][()]

        # Don't close the file --- not loading arrays into memory
        return


    # k is the grid index in the azimuthal direction 
    def contour_plot(self, variable, width=0.7, levels=[None,], Nlevels=20, k=0, mag=False, balance=False, **kwargs):
        if mag:
            variable = dabs(variable)

        f = variable
        if (len(f[0].shape) == 3):
            f = [group_array[:,:,k] for group_array in variable]
        R = [group_array[:,:,k] for group_array in self.R]
        Z = [group_array[:,:,k] for group_array in self.Z]

        levels = np.array(levels) # Convert to array if not already an array

        if levels[0] == None:
            fmax = dmax(f)
            fmin = dmin(f)
            if balance:
                fmax = np.amax([fmax, np.abs(fmin)])
                fmin = - fmax
            levels = np.linspace(1.01*fmin, 0.99*fmax, Nlevels)
        
        print("Contour minimum: ", levels[0])
        print("Contour maximum: ", levels[-1])

        # Make separate plots for each group
        if self.Ngroup > 1:
            R = self.connect_groups(R)
            Z = self.connect_groups(Z)
            f = self.connect_groups(f)
            
        plots = [plt.contour(R[ig], Z[ig], f[ig], levels=levels, linewidths=width, zorder=5, **kwargs) for ig in range(self.Ngroup)]

        #plt.colorbar(mappable=plots[-1]) # Colorbar from last group only; hard to get colorbar applicable for all groups

        # Make standalone custom colorbar --- doesn't show contour levels but at least shows colors and limits
        Ncc = 20 # Nlevels at which you cease to do a "careful colorbar", labeling each contour level
        cmap = plots[0].cmap
        boundaries = np.ndarray((len(levels)+1,))
        boundaries[0]    = levels[0]
        boundaries[1:-1] = 0.5 * np.array(levels[:-1] + levels[1:])
        boundaries[-1]   = levels[-1]
        if len(levels) < Ncc:
            boundaries[0]  -= 0.5*(levels[1]-levels[0]) # Extend a half-step to get lower color boundary
            boundaries[-1] += 0.5 * (levels[-1] - levels[-2])
        norm = mpl.colors.BoundaryNorm(boundaries, cmap.N)
        cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap))
        if len(levels) < Ncc:
            cbar.set_ticks(levels)

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


    ### For now, assume dealing with a strictly axisymmetric simulation in torus geometry
    def find_poloidal_flux_function(self):
        for ig in range(1,5):
            BR = self.B0[ig][:,:,0]
            BZ = self.B1[ig][:,:,0]
            R  = self.R[ig][:,:,0]
            Z  = self.Z[ig][:,:,0]
            (N0,N1) = BR.shape # dir-0 is normal to outer boundary, dir-1 is tangential
            self.PoloidalFlux[ig] = np.zeros((N0,N1))
            BRa = np.zeros((N0,N1))
            BZa = np.zeros((N0,N1))
            BRa[-1,:] = BR[-1,:] # Outer surface -- not taking account of normal flux going to 
            BZa[-1,:] = BZ[-1,:] # zero on the boundary yet...
            BRa[:-1,:] = 0.5 * (BR[:-1,:] + BR[1:,:])
            BZa[:-1,:] = 0.5 * (BZ[:-1,:] + BZ[1:,:])

            RdR = np.zeros((N0,N1))
            RdZ = np.zeros((N0,N1))
            RdR[-1,:] = R[-1,:] * 0.2 * (R[-2,:] - R[-1,:]) # dR_outer ~ 0.2 dR of neighbouring cell?
            RdZ[-1,:] = R[-1,:] * 0.2 * (Z[-2,:] - Z[-1,:])
            RdR[:-1,:] = 0.5 * (R[:-1,:] + R[1:,:]) * (R[:-1,:] - R[1:,:])
            RdZ[:-1,:] = 0.5 * (R[:-1,:] + R[1:,:]) * (Z[:-1,:] - Z[1:,:])

            ### Now integrate inwards from outer boundary
            for j in range(N1): # tangential direction
                self.PoloidalFlux[ig][N0-1,j] = BZa[N0-1,j]*RdR[N0-1,j] - BRa[N0-1,j]*RdZ[N0-1,j]
                for i in reversed(range(N0-1)): # radial direction, inward from boundary
                    self.PoloidalFlux[ig][i,j] = (self.PoloidalFlux[ig][i+1,j] +
                                                    BZa[i,j]*RdR[i,j] - BRa[i,j]*RdZ[i,j])

        ### Store the maximum value of the poloidal flux function reached on the outermost
        ### line of solution points.
        self.MaxPFOuterSurface = np.amax([np.amax(x[-1,:]) for x in self.PoloidalFlux[1:]]) 

        ### Central group --- integrate from group 2 since their coordinate axes are aligned
        BR = self.B0[0][:,:,0]
        BZ = self.B1[0][:,:,0]
        R  = self.R[0][:,:,0]
        Z  = self.Z[0][:,:,0]
        (N0,N1) = BR.shape # dir-0 is normal to outer boundary, dir-1 is tangential
        self.PoloidalFlux[0] = np.zeros((N0,N1))
        BRa = np.zeros((N0,N1))
        BZa = np.zeros((N0,N1))
        ### The top of group 0 attaches to the bottom of group 2
        BRa[-1,:] = 0.5 * (BR[-1,:] + self.B0[2][0,:,0])
        BZa[-1,:] = 0.5 * (BZ[-1,:] + self.B1[2][0,:,0])
        BRa[:-1,:] = 0.5 * (BR[:-1,:] + BR[1:,:])
        BZa[:-1,:] = 0.5 * (BZ[:-1,:] + BZ[1:,:])

        RdR = np.zeros((N0,N1))
        RdZ = np.zeros((N0,N1))
        RdR[-1,:] = 0.5 * (R[-1,:] + self.R[2][0,:,0]) * (R[-1,:] - self.R[2][0,:,0])
        RdZ[-1,:] = 0.5 * (R[-1,:] + self.R[2][0,:,0]) * (Z[-1,:] - self.Z[2][0,:,0])
        RdR[:-1,:] = 0.5 * (R[:-1,:] + R[1:,:]) * (R[:-1,:] - R[1:,:])
        RdZ[:-1,:] = 0.5 * (R[:-1,:] + R[1:,:]) * (Z[:-1,:] - Z[1:,:])
        
        for j in range(N1): 
            self.PoloidalFlux[0][N0-1,j] = self.PoloidalFlux[2][0,j] + BZa[N0-1,j]*RdR[N0-1,j] - BRa[N0-1,j]*RdZ[N0-1,j]
            for i in reversed(range(N0-1)): 
                self.PoloidalFlux[0][i,j] = (self.PoloidalFlux[0][i+1,j] +
                                                BZa[i,j]*RdR[i,j] - BRa[i,j]*RdZ[i,j])
        
        ### Blend the three other edges of group 0 --- isn't usually necessary, and is a sign something's going on...
        ### Outermost set of solution points --- set equal to values from the other groups
        #self.PoloidalFlux[0][:,0]  = self.PoloidalFlux[1][0,:]       # with group 1
        #self.PoloidalFlux[0][:,-1] = self.PoloidalFlux[3][0,:][::-1] # group 3 -- dir-1 is reversed
        #self.PoloidalFlux[0][0,:]  = self.PoloidalFlux[4][0,:][::-1] # group 4 -- "      "   " 

        return


