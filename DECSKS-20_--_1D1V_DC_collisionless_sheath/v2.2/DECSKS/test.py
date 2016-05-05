import numpy as np

Nx = 25

D = np.zeros([Nx,Nx])
for i in range(Nx):
    if i == 0 or i == Nx-1:
        D[i,i] = 1
    else:
        D[i,i-1] = 1
        D[i,i] = -2
        D[i,i+1] = 1

