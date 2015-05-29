import numpy as np
import matplotlib.pyplot as plt

# the following errors were outputted from successive CS tests with the corresponding mesh sizes

error_norm = [0.112811293016,
              0.0505443821755,
              0.00552345804204,
              0.0001298328338,
              1.6330421753e-07,
              9.90644408585e-11,
              5.14652228979e-14]

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

print "order calculations at each refinement step: LTE order = 12"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

fig, ax = plt.subplots()
order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{12})$, order = %g' % order)
ax.hold('on')

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.106422187398,
              0.0488785772965,
              0.00505269856632,
              9.51438114925e-05,
              6.96411113297e-08,
              2.20092669395e-11,
              1.71158396293e-14]

orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 13"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{13})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.106249913515,
              0.0485198076544,
              0.00475066758542,
              6.27180711852e-05,
              2.64531223041e-08,
              4.36603744294e-12,
              1.64944533995e-14]

orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 14"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{14})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.101612572432,
              0.047071802681,
              0.00446399010269,
              4.75496675788e-05,
              1.19766404584e-08,
              1.07334055089e-12,
              3.04422903091e-13]


orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 15"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{15})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.105373121488,
              0.04653695228,
              0.0042717525573,
              3.31014669033e-05,
              4.90168891124e-09,
              4.33789442272e-13,
              2.99567928258e-13]

orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 16"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{16})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.0875272425701,
              0.0437463991949,
              0.00374114468444,
              1.13961974703e-05,
              6.22569862759e-10,
              1.77542201736e-10,
              8.19865444473e-11]

orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 22"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{22})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()



ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

# shrink frame of the plot to make room for a legend on the right side
frame = ax.get_position() # position of plot center
ax.set_position([frame.x0, frame.y0,
                 frame.width * 0.8, frame.height])

# Place legend to the right of the shrunken frame
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
          fancybox=True, ncol=1)
plt.grid()
plt.show()
