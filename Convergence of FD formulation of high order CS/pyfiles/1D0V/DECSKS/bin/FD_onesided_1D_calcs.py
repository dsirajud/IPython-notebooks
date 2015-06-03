# FD 5
import numpy as np
import math

Nx = 6
Nxs = [32, 64, 128, 256, 512, 1024]
L2 = np.zeros(Nx)
L2[0] = 7.0675680166488327738072e-02
L2[1] = 2.2476804189265533334696e-02
L2[2] = 4.0535113755227995813168e-03
L2[3] = 4.8341524068649076306681e-05
L2[4] = 2.4821603420530102303287e-04
L2[5] = 2.7715212113163948184913e-03

for q in range(1,Nx):
    intro_str = "for Nx%d/Nx%d: " % (Nxs[q-1],Nxs[q])
    order =  math.log(L2[q-1] / L2[q] , 2)
    order_str = 'order = %f' % order
    print intro_str + order_str
