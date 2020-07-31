import plot_basic as pb
import matplotlib.pyplot as plt
import os

Nfiles = 299
movie_dir = '/Users/kyle/Desktop/movie/'

os.chdir('data/')

s = pb.Snapshot()

plt.figure()

for i in range(Nfiles):
    s.load_data(i)
    plt.clf()
    
    s.draw_edges(0.2)
    s.contour_plot(var='rho', width=1.5)

    plt.savefig(movie_dir+'%04d.png' % i, dpi=150)

print("All done!")

    
