
import numpy as np

# .................................................................#
# (L) linear increasing/decreasing piecewise velocity distribution #
def v_L(x):
    v_low = -0.25
    v_high = 0.75
    return np.where(x < 0.0, (v_low + 1*x/3.), (v_high - 2*x/3.))

# .................................................................#
# +/- (H) step function velocity distribution #
def v_H(x):
    v_low = -1
    v_high = 2
    return np.where(x < 0, v_low,v_high)

# .................................................................#
# +/- linear, const, sinusoid (LCS) velocity distribution #
def v_LCS(x):

    v = np.zeros(len(x))
    v_low = -1
    T      = 0.4
    v_const = 1

    for each in range(len(x)):
        if x[each] < 0:
            v[each] = v_low
        elif 0 <= x[each] < 0.10:
            v[each] = 1 + 2*x[each]/1.
        else:
            v[each] = 1.5 - 1.*np.sin(2*np.pi*x[each]/T)
            #else:
            #v[each] = 15 + 20*x[each]/30

    return v


# .................................................................#
# sinusoid (s) velocity distribution #

def v_S(x):

    return np.sin(2*np.pi*x / 1.0)

# .................................................................#
# (L) linear increasing/decreasing piecewise velocity distribution #
def v_0(x):
    v = 1.
    vv = v*np.ones(len(x))

    return vv

# .................................................................#
# (L) linear increasing/decreasing piecewise velocity distribution #
def v_0n(x): #n = negative
    v = 5.
    vv = -v*np.ones(len(x))

    return vv

# .................................................................#
# rotation field when using q = {x,-y}

def v_2piz(z): #n = negative
    v = 2*np.pi*z
    
    return v
