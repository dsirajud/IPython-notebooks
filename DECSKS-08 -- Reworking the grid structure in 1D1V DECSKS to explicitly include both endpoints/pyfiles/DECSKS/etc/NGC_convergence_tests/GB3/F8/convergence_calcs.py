#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

Nx = [24,
      48,
      96,
      192,
      384,
      768]


L2 = [0.0334109577276,
      0.000138773646963,
      5.12971008836e-08,
      9.88269132188e-11,
      1.90117467797e-13,
      1.53600267929e-15]


order = -np.polyfit(np.log2(Nx), np.log2(L2), 1)[0]

print "Nx                L2 error                order"
print "--------------------------------------------------------"
for i in range(len(Nx)):
    if i == 0:
        print "%d                %2.5e" % (Nx[i], L2[i])
    else:
        order = np.log2(L2[i-1] / L2[i])
        print "%d                %2.5e                %g" % (Nx[i], L2[i], order)


fig, ax = plt.subplots()
ax.loglog(Nx, L2, '--o')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.ylabel('$L^2$ norm of $f_{num} - f_{exact}$', fontsize = 16)
plt.xlabel('number of gridpoints $N_x$', fontsize = 16)
plt.text(np.max(Nx), np.max(L2), r'order ($|slope|$) = %g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
plt.grid()
plt.show()
