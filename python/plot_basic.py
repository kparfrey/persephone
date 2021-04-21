import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
#from scipy.interpolate import griddata

class Group(object):
    group = None

    r0 = None
    r1 = None
    r  = None

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





class Snapshot(object):
    m       = None # Mesh object
    dfile   = None # H5 file for data
    filenum = None # Number of loaded data file
    time    = None
    step    = None
    dt      = None

    Ngroup  = None


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
        self.v   = [0] * self.Ngroup
        
        for ig in range(self.Ngroup):
            sg = str(ig)
            self.rho[ig] = self.dfile[sg]['rho']
            self.p[ig]   = self.dfile[sg]['p']
            self.v0[ig]  = self.dfile[sg]['v0']
            self.v1[ig]  = self.dfile[sg]['v1']
            self.v2[ig]  = self.dfile[sg]['v2']
            self.v[ig]   = [self.v0, self.v1, self.v2]

        # Can't make "global" v array, since the groups can have
        # different numbers of points, elements etc.
        #self.v   = np.array((self.v0, self.v1, self.v2))

        self.time = self.dfile['time'][()]
        self.step = self.dfile['step'][()]
        self.dt   = self.dfile['dt'][()]

        # Don't close the file --- not loading arrays into memory
        return


    # Assume 0-1 plane
    def contour_plot(self, var='rho', width=0.7, levels=[None,]):
        for ig in range(self.Ngroup):
            sg = str(ig)
            r0 = self.m.g[ig].r0
            r1 = self.m.g[ig].r1
            
            if levels[0] == None:
                levels = np.linspace(0.0, 0.99, 20)
            
            plt.contour(r0[:,:,0], r1[:,:,0], self.dfile[sg][var][:,:,0], levels=levels,
                                                             linewidths=width, zorder=5)

        plt.title('t = %.4lf' % self.time)

        ax = plt.gca()
        ax.set_aspect('equal')
        return
