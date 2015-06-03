#!/usr/bin/env python
#
# This script generates two types of outputs (see __main__ routine at bottom):
#
#    (1) decorated full summary of all families of explicit finite differences schemes
#        for derivative specified queried by user at given LTE on CS
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
#           ../etc/f1_LTE_N-1_FD_coefficients.dat
#           ../etc/f2_LTE_N-2_FD_coefficients.dat
#           , etc., one for each dn-th derivative dnf
#
#
# Quickstart:
#
# the script takes in a system arguments corresponding to the order of the
# LTE truncation error desired on the derivatives and the derivative desired dn
#
# $ python generate_tables_of_finite_difference_schemes_for_a_given_LTE_on_CS LTE
#
# where, LTE is an int
#
# see __main__ routine below for execution details

import generate_tables_of_finite_difference_schemes_for_a_given_LTE

def main(LTE_CS = 5):
    N = LTE_CS - 1
    for dn in range(1,N):
        generate_tables_of_finite_difference_schemes_for_a_given_LTE.main(LTE = N - dn, dn = dn)

if __name__ == '__main__':
    import sys

    N = int(sys.argv[1])
    for dn in range(1,N):
        generate_tables_of_finite_difference_schemes_for_a_given_LTE.main(LTE = N - dn, dn = dn)
