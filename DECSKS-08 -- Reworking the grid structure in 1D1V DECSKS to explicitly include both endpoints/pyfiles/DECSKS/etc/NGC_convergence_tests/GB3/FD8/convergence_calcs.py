#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

Nx = [24,
      48,
      96,
      192,
      384,
      768,
      1536]

L2 = [0.0545892017126,
      0.00626735872757,
      8.62608345024e-05,
      2.48298856624e-07,
      5.53991079128e-10,
      1.12575793363e-12,
      4.16018116531e-15]


print "Nx                L2 error                order"
print "--------------------------------------------------------"
for i in range(len(Nx)):
    if i == 0:
        print "%d                %2.5e" % (Nx[i], L2[i])
    else:
        order = np.log2(L2[i-1] / L2[i])
        print "%d                %2.5e                %g" % (Nx[i], L2[i], order)


order = -np.polyfit(np.log2(Nx[1:]), np.log2(L2[1:]), 1)[0]

fig, ax = plt.subplots()
ax.loglog(Nx, L2, '--o')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.ylabel('$L^2$ norm of $f_{num} - f_{exact}$', fontsize = 16)
plt.xlabel('number of gridpoints $N_x$', fontsize = 16)
plt.text(np.max(Nx), np.max(L2), r'order ($|slope|$) = %g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
plt.grid()
plt.show()
