# For DCT details see 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.fftpack.dct.html

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct

N = 8 # Number of Gauss points

def f(x):
    return 2*x**3 - 5 * x**2 + 4*x + 0.5

M    = np.ndarray((N,N)) # Nodes to modes
Minv = np.ndarray((N,N)) # Modes to nodes
F    = np.zeros((N,N))   # Final total filtering matrix
L    = np.ndarray((N,)) # diagonal filter matrix

x = np.ndarray((N,))

for i in range(N):
    x[i] = 0.5 * (1.0 - np.cos(np.pi*(i+0.5)/N))


fx = f(x)

modes = dct(fx, type=2)

# Type 2 DCT
for k in range(N):
    for j in range(N):
        M[k,j] = 2.0 * np.cos(0.5 * np.pi * k * (2*j + 1)/N)

# Type 3 DCT, multiplied by 0.5 / N 
norm = 0.5/N
for j in range(N):
    Minv[j,0] = norm * 1.0
    for k in range(1,N):
        Minv[j,k] = norm * 2.0 * np.cos(0.5 * np.pi * k * (2*j + 1)/N)

alpha = 36.0
s = 4.0

for i in range(N):
    eta = i/(N-1.0)
    L[i] = np.exp(-alpha * eta**s)


LM = np.zeros((N,N))
for i in range(N):
    filtval = L[i]
    for j in range(N):
        LM[i,j] = filtval * M[i,j] # Since L is diagonal

for i in range(N):
    for j in range(N):
        for m in range(N):
            F[i,j] += Minv[i,m] * LM[m,j]


#########################################


modes_matrix = np.zeros((N,))
for i in range(N):
    for j in range(N):
        modes_matrix[i] += M[i,j] * fx[j]

fx1 = np.zeros((N,))
for i in range(N):
    for j in range(N):
        fx1[i] += Minv[i,j] * modes_matrix[j]

fx_filt = np.zeros((N,))
for i in range(N):
    for j in range(N):
        fx_filt[i] += F[i,j] * fx[j]

modes_filt = np.zeros((N,))
for i in range(N):
    for j in range(N):
        modes_filt[i] += LM[i,j] * fx[j]


plt.figure()
plt.plot(x, fx, linewidth=1.5)
#plt.plot(x, fx,'.')
#plt.plot(x, fx1, linewidth=0.5)
plt.plot(x, fx_filt, linewidth=1.5)

plt.figure()
plt.plot(modes, linewidth=1.5)
#plt.plot(modes,'.')
#plt.plot(modes_matrix,linewidth=0.5)
plt.plot(modes_filt,linewidth=1.5)

plt.figure()
plt.plot(L)

plt.show()
