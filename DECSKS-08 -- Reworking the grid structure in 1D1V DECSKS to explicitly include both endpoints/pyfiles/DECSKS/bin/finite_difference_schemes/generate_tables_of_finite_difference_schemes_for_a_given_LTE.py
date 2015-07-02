#!/usr/bin/env python
#
# This script generates two types of outputs (see __main__ routine at bottom):
#
#    (1) decorated full summary of all families of explicit finite differences schemes
#        for derivative specified queried by user at given LTE
#
#    (2) individual files for each derivative which contains each family of
#        schemes (e.g. 'forward', 'backward', 'central') with varying degrees
#        of asymmetry.
#
#    *** (1) is generated as a presentable table that can be viewed by the user
#
#           default filepath:
#           ../etc/Table_of_finite_difference_schemes_at_const_LTE_decorated.dat
#
#    *** (2) is generated in a particular format that a read function knows about
#            so that the tables of [weights, stencil] can be read into and stored
#            in a dictionary FD_schemes in store_schemes(dn_max)
#
#           default filepaths for (2): 
#           ../etc/f1_FD_coefficients.dat
#           ../etc/f2_FD_coefficients.dat
#           , etc., one for each dn-th derivative dnf
#
#
# Quickstart:
#
# the script takes in a system arguments corresponding to the order of the
# LTE truncation error desired on the derivatives and the derivative desired dn
#
# $ python generate_tables_of_finite_difference_schemes_for_a_given_LTE LTE dn
#
# where, LTE and dn are both integers
#
# see __main__ routine below for execution details

def find_weights(_dn, _p, _stencil):
    """Assemble and solve the matrix system characterizing
    the explicit scheme of LTE order p for the derivative
    order dn according to stencil provided

    _ indicates variable is restricted to scope of this function

    inputs:
    _dn -- (int) order of derivative
    _p -- (int) order of truncation error desired
    _stencil -- (list, int) list of integers i in x + i(dx)
         the choice of scheme (asymmetric (forward, backward,
         unequal sampling on each side) vs. symmetric (centered)
         will determine this range), A[i,:] = x ** i

    outputs:
    _w -- (ndarray, ndim=2) finite difference coefficients from solving AC = b
    _stencil -- (list, int) the integer i in (x + ih) the weight applies to
         i.e. returns the input
    """
    from math import factorial
    import numpy as np
    import numpy.linalg as LA

    eps = np.finfo(np.float).eps # smallest floating point number
    _stencil = np.array(_stencil)
    A = np.zeros([_dn+_p, _dn+_p])
    b = np.zeros(len(_stencil))

    for i in range(_dn+_p):
        A[i,:] = _stencil ** i

        if i == _dn:
            b[i] = 1.
        else:
            b[i] = 0.

    _w = factorial(_dn)*LA.solve(A,b) # = w[0], w[1], .., w[p+d-1]
    _w[np.abs(_w) < 2*eps] = 0 # zero out any weight close (factor of 2) to machine epsilon
    return _w, _stencil

def main(LTE = 2, dn = 1, rel_path = './../../etc/finite_difference_schemes/'):
    import misc
    import os
    import sys
    import numpy as np # for modulo operation

    decorated_out, dn_schemes_out = misc.write_header(dn, LTE)
    # generate finite difference schemes
    stencil_size = LTE + dn
    stencil_center = stencil_size // 2
    jmax = stencil_size - 1

    for j in range(stencil_size):
        if j < (int( (LTE + dn)/2)):
            label = 'F' + str(j)
        elif j == stencil_center and np.mod(stencil_size, 2) == 1:
            label = 'C0'
        elif j == jmax:
            label = 'B0'
        else:
            label = 'B' + str(jmax - j)

        stencil = range(0 - j, LTE + dn - j)
        w, stencil = find_weights(dn, LTE, stencil)

        decorated_out, dn_schemes_out = misc.write_rest(decorated_out,
                                                 dn_schemes_out,
                                                 label,
                                                 LTE,
                                                 dn,
                                                 stencil,
                                                 w)

    decorated_out = misc.write_footer(decorated_out)
    decorated_out.close()
    dn_schemes_out.close()

    # files created in bin/, move to etc/ via os.rename()
    # move files to ../etc via os.rename()
    os.rename('./' + dn_schemes_out.name, rel_path + dn_schemes_out.name)
    os.rename('./' + decorated_out.name, rel_path + decorated_out.name)

if __name__ == '__main__':

    import misc
    import os

    import sys
    import numpy as np # for modulo operation

    LTE = int(sys.argv[1])  # LTE
    dn = int(sys.argv[2]) # order of derivative whose schemes are desired


    decorated_out, dn_schemes_out = misc.write_header(dn, LTE)
    # generate finite difference schemes
    stencil_size = LTE + dn
    stencil_center = stencil_size // 2
    jmax = stencil_size - 1

    for j in range(stencil_size):
        if j < (int( (LTE + dn)/2)):
            label = 'F' + str(j)
        elif j == stencil_center and np.mod(stencil_size, 2) == 1:
            label = 'C0'
        elif j == jmax:
            label = 'B0'
        else:
            label = 'B' + str(jmax - j)

        stencil = range(0 - j, LTE + dn - j)
        w, stencil = find_weights(dn, LTE, stencil)

        decorated_out, dn_schemes_out = misc.write_rest(decorated_out,
                                                 dn_schemes_out,
                                                 label,
                                                 LTE,
                                                 dn,
                                                 stencil,
                                                 w)

    decorated_out = misc.write_footer(decorated_out)
    decorated_out.close()
    dn_schemes_out.close()


    # files created in bin/, move to etc/ via os.rename()
    # move files to ../etc via os.rename()
    rel_path = './../etc/' 
    os.rename('./' + dn_schemes_out.name, rel_path + dn_schemes_out.name)
    os.rename('./' + decorated_out.name, rel_path + decorated_out.name)
