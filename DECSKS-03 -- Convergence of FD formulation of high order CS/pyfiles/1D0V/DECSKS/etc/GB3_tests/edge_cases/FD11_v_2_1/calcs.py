import numpy as np
import matplotlib.pyplot as plt

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.117123029204,
              0.0960259305052,
              0.0186822543343,
              0.00236022919004,
              0.0114955963398,
              0.0358812363677,
              0.0343273342871]

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
