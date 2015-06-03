#!/usr/bin/env python

# ./DECSKS/bin/finite_difference_schemes/
# generate_tables_of_finite_difference_schemes.py
#
# This script generates two types of outputs (see __main__ routine at bottom):
#
#    (1) full summary of all families of explicit finite differences schemes
#        for all derivatives specified
#
#    (2) individual files for each derivative which contains each family of
#        schemes (e.g. 'forward', 'backward', 'central') with varying degrees
#        of asymmetry and also varying error on each derivative according to the
#        prescription on high order CS, see notes below
#
#    *** (1) is generated as a presentable table that can be viewed by the user
#
#   default filepath:
#   ./DECSKS/etc/finite_difference_schemes/Tables_of_finite_difference_schemes.dat
#
#    *** (2) is generated in a particular format that a read function knows about
#            so that the tables of [weights, stencil] can be read into and stored
#            in a dictionary FD_schemes in store_schemes(dn_max)
#
#   default filepath:
#   ./DECSKS/etc/finite_difference_schemes/f1_FD_coefficients.dat
#   ./DECSKS/etc/finite_difference_schemes/f2_FD_coefficients.dat
#   , etc., one for each derivative fdn
#
#
# Quickstart:
#
# the script takes in a system argument corresponding to the order of the global
# error N of the method
#
# $ python generate_table_of_finite_difference_schemes_required_for_a_chosen_GE_on_CS N
#
# where N is an integer
#
# Notes:
#
#         To correct CS to be of global error O(N), what is needed are derivative
#         coefficients d in Uf[i] = Sum_q c[q] d[q,i] for the qth derivative to have an
#         *LTE* of O(N), since an expansion of the CS update takes the form:
#
#         f_CS(t,z) = f(t,z) +
#                     \sum_p (-1)^p \frac{(\Delta z)^p}{p!}\frac{\partial^p}U(t,x)f(t,x)
#                     + O(\Delta z^{N+1})
#
#
#         Since terms dz ** q (q >= 1) multiply the flux Uf = Sum_q c[q] d[q,i],
#         then LTE[d] >= O(N) so that LTE[f_CS] = O(N+1) which gives
#         a GE[f_CS] = O(N)
#
#         Now, d itself is also multiply by powers of dz
#         i.e. d = dz ** dn * dnf, where dnf is the dn-th derivative of f, dn is the
#                             derivative number
#
#         Thus, for an Nth order CS scheme, this product needs to have LTE O(N)
#         which implies. Labeling the LTE on the derivative as p,
#
#         dn = 1, LTE(dnf) = p = O(N-1) so that LTE[dz * dnf] = O(N)
#         dn = 2, LTE(dnf) = p = O(N-2) so that LTE[dz**2 * dnf] = O(N)
#         dn = 3, LTE(dnf) = p = O(N-3) so that LTE[dz**3 * dnf] = O(N)
#         ... and so on
#
#         i.e. for each dn, truncation error on dnf, p is related to the
#         global error N = GE by
#
#         p = GE - dn, or p = N - dn for each dn

def HighPrecisionE(number):
    """Converts a number into a string object
    while retaining a chosen degree of precision. This
    is designed to evade the truncation that is involved
    with str() so that outputs can store numbers with high
    precision

    inputs:
    number -- (number)

    outputs:
    string object with chosen precision in scientific notation
    """

    return "%.22e" % number

def strings_for_writeout():
    """creates a dictionary of the strings below used in
    the outfile writing decoration

    inputs:
    None

    outputs:
    table_strings -- (dict) decorating strings used in output writing
    """
    newline = '\n'
    divider = '-----------------------------------'
    spc = ' '
    colon = ':'
    equaldivider = '========================='
    outfilename_suffix = '_FD_coefficients.dat'
    number_of_lines_to_read = 'number of lines to read: '

    table_strings = dict(newline = newline,
                         divider = divider,
                         spc = spc,
                         colon = colon,
                         equaldivider = equaldivider,
                         outfilename_suffix = outfilename_suffix,
                         number_of_lines_to_read = number_of_lines_to_read
                        )

    return table_strings

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
    return _w, _stencil

def write_scheme(outfile, _w, _stencil):
    """For a provided outfile by the user, write
    the scheme (no plural!) according to the provided finite difference
    weights (_w) and stencil (_stencil)

    inputs:
    outfile -- (file) an opened ('w' mode) file to be written to
    _w -- (list) weights of the finite difference scheme
    _stencil -- (list) corresponding stencil for the scheme

    outputs:
    outfile -- (file) return the same file after writing is complete
    """
    length = len(stencil)

    for i in range(length):
        outfile.write( HighPrecisionE(_w[i]) + ', '
                       + str(_stencil[i]) + '\n')

    return outfile


if __name__ == '__main__':

    import os
    import sys
    import numpy as np
    import numpy.linalg as LA

    N = int(sys.argv[1]) # specified global error in phase space for CS scheme
    dn_max = N - 1 # maximum derivative in range

    strings = strings_for_writeout()

    # store as compact names for brevity
    newline = strings['newline']
    divider = strings['divider']
    spc = strings['spc']
    colon = strings['colon']
    equaldivider = strings['equaldivider']
    outfilename_suffix = strings['outfilename_suffix']
    number_of_lines_to_read = strings['number_of_lines_to_read']

    # script intended to be run in directory
    # ./DECSKS/bin/finite_difference_schemes/ and writes to
    # ./DECSKS/etc/finite_difference_schemes/
    rel_path = './../../etc/finite_difference_schemes/'

    # open file for writing full table with all schemes for all derivatives
    fullsummary_out = open(rel_path + 'Tables_of_finite_difference_schemes.dat','w')

    # generate schemes and write to tables (output files) for each dn
    for dn in range(1,N):

        p = N - dn       # LTE error p requirement decreases with each derivative
                         # see notes in top matter, LTE[d] = O(N) is what is
                         # needed, LTE[dnf] = O(p) = O(N - dn), dnf = dn-th
                         # derivative

        # header for output
        derivative = 'f^(' + str(dn) + ')'
        fullsummary_out.write(derivative + colon + spc + 'LTE = O('
                      + str(p) + ')' + '\n')
        fullsummary_out.write(equaldivider + newline + newline)

        # open files for each derivative
        # f1 + outfilename_suffix, f2 + outfilename, ...
        # these files are read in by finite_difference_schemes.py
        # read(*args), and store(*kwargs) functions
        single_dn = 'f' + str(dn)
        single_dn_outname = single_dn + outfilename_suffix
        single_dn_schemes_out = open(single_dn_outname, 'w')

        nmax = p + dn - 1
        for n in range(int((p+dn))):
            if n < (int( (p+dn)/2)):
                label = 'forward ' + str(n)
            elif n == int((p + dn) / 2.) and np.mod(p+dn, 2) == 1:
                label = 'central 0'
            elif n == nmax:
                label = 'backward 0'
            else:
                label = 'backward ' + str(nmax - n)

            fullsummary_out.write(derivative + colon + spc +
                                    label + '\n')

            fullsummary_out.write(divider + newline)

            stencil = range(0 - n, p + dn - n)
            w, stencil = find_weights(dn,p,stencil)
            fullsummary_out.write('(weights, stencil)' + newline)
            fullsummary_out = write_scheme(fullsummary_out,
                                             w, stencil)
            single_dn_schemes_out.write(label + newline)
            single_dn_schemes_out.write(number_of_lines_to_read + str(p+dn) + newline)
            single_dn_schemes_out = write_scheme(single_dn_schemes_out,
                                         w, stencil)


            fullsummary_out.write(divider + newline)
        fullsummary_out.write(equaldivider + newline)


        single_dn_schemes_out.close()
        # move file to ./DECSKS/etc/finite_difference_schemes/ vis os.rename()
        os.rename('./' + single_dn_outname, rel_path + single_dn_outname)
    fullsummary_out.close()
