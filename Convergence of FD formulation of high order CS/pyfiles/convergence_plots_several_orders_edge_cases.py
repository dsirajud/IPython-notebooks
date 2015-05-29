import numpy as np
import matplotlib.pyplot as plt

# the following errors were outputted from successive CS tests with the corresponding mesh sizes

error_norm = [0.0675141216017, # Nx24
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
orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 8"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

fig, ax = plt.subplots()
order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{8})$, order = %g' % order)
ax.hold('on')

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.0933791494565,
         0.0177063387756,
         0.0173742156515,
         0.000437814728139,
         2.03932489416e-06,
         9.9955975972e-09,
         4.0101019395e-11]

orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 9"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{9})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.101172411677,
              0.0202810036715,
              0.0114927064009,
              0.000300056419446,
              5.8438700153e-07,
              2.15424422863e-09,
              0.00164369969121]

orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 10"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{10})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.110412229688,
              0.0495531129451,
              0.0249962581766,
              0.000741144051204,
              8.43100699903e-05,
              0.0349792732794,
              0.0395441636386]

orders = np.zeros(len(Nx))
dx = 2./Nx

print "order calculations at each refinement step: LTE order = 11"

for n in range(len(Nx)):
    if n == 0:
        print "Nx%d        error = %g       ----" % (Nx[n], error_norm[n])
    else:
        orders[n] = np.log2(error_norm[n-1] / error_norm[n])
        print "Nx%d        error = %g       order = %g" % (Nx[n], error_norm[n], orders[n])
print '\n'

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{11})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.117123029204,
              0.0960259305052,
              0.0186822543343,
              0.00236022919004,
              0.0114955963398,
              0.0358812363677,
              0.0343273342871]

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

order = np.max(orders)
ax.loglog(Nx, error_norm, '-o', label = '$O(\Delta x^{12})$, order = %g' % order)
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.grid()

####################

# the following errors were outputted from successive CS tests with the corresponding mesh sizes
error_norm = [0.180565107006,
              0.151931643039,
              0.0472540483564,
              0.043555005061,
              0.0281880779251,
              0.0258890357207,
              0.0232929041963]

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
error_norm = [0.202096063037,
              0.142659964027,
              0.0448100740151,
              0.031347715831,
              0.0208353930884,
              0.011856210771,
              0.00850951657219]

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




# shrink frame of the plot to make room for a legend on the right side
frame = ax.get_position() # position of plot center
ax.set_position([frame.x0, frame.y0,
                 frame.width * 0.8, frame.height])

# Place legend to the right of the shrunken frame
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
          fancybox=True, ncol=1)
plt.grid()
plt.show()
