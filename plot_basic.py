import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
#from scipy.interpolate import griddata

class Mesh(object):
    Nelem_tot = None
    Nelem     = None
    Ns        = None
    Nf        = None

    def __init__(self):
        print("Loading mesh file")
        filename = 'mesh.h5'
        m = h5py.File(filename)

        self.r0 = m['coords']['r']['0']
        self.r1 = m['coords']['r']['1']
        self.r2 = m['coords']['r']['2']
        self.r = np.array((self.r0,self.r1,self.r2))
        
        self.e0 = m['edges']['0']
        self.e1 = m['edges']['1']
        self.e2 = m['edges']['2']
        self.e3 = m['edges']['3']
        # List rather than array since edges don't all have same length
        self.edges = [self.e0, self.e1, self.e2, self.e3]

        self.Nelem_tot = m['Nelem_tot'][()]
        self.Nelem  = m['Nelem'][:]
        self.Ns     = m['Ns'][:]
        self.Nf     = m['Nf'][:]

        # Don't close the file --- not loading arrays into memory
        return


class Snapshot(object):
    m       = None # Mesh object
    dfile   = None # H5 file for data
    filenum = None # Number of loaded data file
    time    = None
    step    = None
    dt      = None

    def __init__(self, filenum = 0):
        self.m = Mesh()
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


    def draw_edges(self, width = 0.2):
        m = self.m

        for i in range(m.Nelem_tot):
            for j in range(4):
                plt.plot(m.edges[j][i,0,:], m.edges[j][i,1,:], 'k-', lw=width, zorder=1) 

        ax = plt.gca()
        ax.set_aspect('equal')
        return


    def draw_soln_points(self, size = 5.0):
        r0 = self.m.r0
        r1 = self.m.r1
        
        # Add c='k' to draw the points in black
        plt.scatter(r0[:,:,0], r1[:,:,0], marker='.', s=size)

        ax = plt.gca()
        ax.set_aspect('equal')
        return


    # Assume 0-1 plane
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

