import plot_basic as pb
import matplotlib.pyplot as plt
import numpy as np
import os

movie_dir = '/Users/kyle/Desktop/mesh_movie/'

os.chdir('data/')

m = pb.Mesh()

Nphi_per_elem = m.g[0].Ns[2]
Nphi = Nphi_per_elem * m.g[0].Nelem[2]

plt.figure()

for i in range(0,Nphi,1):
    plt.clf()
    phi = m.g[0].r2[0,0,i]
    
    m.draw_edges(i/Nphi_per_elem)
    m.draw_soln_points(i, size=3)
    plt.xlim(0.9,3.5)
    plt.xlabel(r"$R$")
    plt.ylabel(r"$Z$")
    plt.title(r'$\phi = %.3lf \pi, n_{\phi} = %d$' % (phi/np.pi,i))
    ax = plt.gca()
    ax.set_aspect('equal')

    plt.savefig(movie_dir+'%04d.png' % i, dpi=150)
    print("Done %d" % i)

print("All done!")

