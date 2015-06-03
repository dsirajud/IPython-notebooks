import numpy as np
import matplotlib.pyplot as plt

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error = [0.0675141216017, # Nx24
         0.0161712472223, # Nx48
         0.00914501177604, # Nx96
         0.000551128433045, # Nx192
         6.74899383061e-06, # Nx384
         5.67282743459e-08, # Nx768
         4.58162827303e-10] # Nx1536 

Nx = [24,
48,
96,
192,
384,
768,
1536]


Nx = np.array(Nx)
dx = 2./Nx

for i in range(1,len(Nx)):
    order = np.log2(error[i-1] / error[i])
    print "Nx%d/Nx%d : order = % g" % (Nx[i-1], Nx[i], order)

fig, ax = plt.subplots()
ax.loglog(Nx, error, 'o')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()
plt.show()
