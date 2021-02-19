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

    def __init__(self, group_id):
        self.group = group_id

        return



class Mesh(object):
    Ngroup    = None
    Nproc     = None
    Nfield    = None
    Nelem_tot = None
    Nelem     = None
    Ns        = None
    Nf        = None
    
    g = dict() # Dictionary of Groups

    def __init__(self):
        print("Loading mesh file")
        filename = 'mesh.h5'
        m = h5py.File(filename)

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

            group.e0 = m[sg]['edges']['0']
            group.e1 = m[sg]['edges']['1']
            group.e2 = m[sg]['edges']['2']
            group.e3 = m[sg]['edges']['3']
            # List rather than array since edges don't all have same length
            group.edges = [group.e0, group.e1, group.e2, group.e3]

            group.Nelem_tot = m[sg]['Nelem_tot'][()]
            group.Nelem = m[sg]['Nelem'][:]
            group.Ns    = m[sg]['Ns'][:]
            group.Nf    = m[sg]['Nf'][:]

            self.g[ig] = group

        # Don't close the file --- not loading arrays into memory
        return


    def draw_edges(self, width = 0.2):
        for ig in range(self.Ngroup):
            edges = self.g[ig].edges
            for i in range(self.g[ig].Nelem_tot):
                for j in range(4):
                    plt.plot(edges[j][i,0,:], edges[j][i,1,:], 'k-', lw=width, zorder=1) 

        ax = plt.gca()
        ax.set_aspect('equal')
        return


    def draw_soln_points(self, size = 5.0):
        for i in range(self.Ngroup):
            r0 = self.g[i].r0
            r1 = self.g[i].r1
            
            # Add c='k' to draw the points in black
            plt.scatter(r0[:,:,0], r1[:,:,0], marker='.', s=size)

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
        self.dfile = h5py.File(filename)

        self.rho = self.dfile['rho']
        self.p   = self.dfile['p']
        self.v0  = self.dfile['v0']
        self.v1  = self.dfile['v1']
        self.v2  = self.dfile['v2']
        self.v   = np.array((self.v0, self.v1, self.v2))

        self.time = self.dfile['time'][()]
        self.step = self.dfile['step'][()]
        self.dt   = self.dfile['dt'][()]

        # Don't close the file --- not loading arrays into memory
        return


    # Assume 0-1 plane
    """
    def contour_plot(self, var='rho', width=0.4, levels=[None,]):
        r0 = self.m.r0
        r1 = self.m.r1
        
        if levels[0] == None:
            levels = np.linspace(0.5, 0.99, 12)
        
        plt.contour(r0[:,:,0], r1[:,:,0], self.dfile[var][:,:,0], levels=levels,
                                                         linewidths=width, zorder=5)
        plt.title('t = %.2lf' % self.time)

        ax = plt.gca()
        ax.set_aspect('equal')
        return
    """

