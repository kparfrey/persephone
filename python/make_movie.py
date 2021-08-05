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
    
    s.m.draw_edges(0.05)
    s.contour_plot(var='p', width=1.0)

    plt.savefig(movie_dir+'%04d.png' % i, dpi=150)

print("All done!")

    
