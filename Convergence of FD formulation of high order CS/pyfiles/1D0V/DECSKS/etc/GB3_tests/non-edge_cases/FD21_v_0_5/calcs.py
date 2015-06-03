import numpy as np
import matplotlib.pyplot as plt

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.0875272425701,
              0.0437463991949,
              0.00374114468444,
              1.13961974703e-05,
              6.22569862759e-10,
              1.77542201736e-10,
              8.19865444473e-11]

Nx = [24,
48,
96,
192,
384,
768,
1536]


Nx = np.array(Nx)
orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step:"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

fig, ax = plt.subplots()
ax.loglog(Nx, error_norm, 'o')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()
plt.show()
