Generates tables of finite difference weights and corresponding stencils for 
derivatives of any order at user-specified local truncation error order

Usage:

    Run these scripts to generate tables of finite difference coefficients
for a chosen order of accuracy (p, LTE) and according to how many derivatives
are needed for a given high order CS (= LTE - 2, e.g. FD5 has global error 
GE = 5, LTE = 6 which requires 4 derivatives, dn = 1, 2, 3, 4). 

finite difference weights and corresonding stencils can then be read by 
./DECSKS/lib/read.py to store all families of schemes in a dictionary 
FD_schemes.py that are retrieved during a CS simulation as needed

Contents:

(1) generate_tables_of_finite_difference_schemes.py ---

    sys.argv[1] = (int) p, order of LTE for each derivative
    sys.argv[2] = (int) dn_max, maximum derivative order needed

    FD5 uses p = 6 (LTE), dn_max = 4 (LTE - 2)

    Generates families of finite difference schemes for all 
    derivatives dn = 1, 2, ... dn_max and writes two types 
    of files:

    (a) Tables_of_finite_difference_schemes.dat ---

        Full summary of all schemes for every derivative,
        formatted in an easy to read table with dividers
        and labels

    (b) fi_FD_coefficients.dat --- (i = 1, 2, 3... dn_max)

        one file for every derivative order dn

        each file contains all the families of schemes formatted
        in a way that is known to reader/storing functions:

            ./finite_difference_schemes/finite_difference_schemes.py
                store_finite_difference_schemes(*args)
                read_finite_difference_schemes(*args, **kwargs)

            + TODO: include these two routines in ./DECSKS/lib/read.py
              

Example usage:

# calculate and create tables for p = 6 (LTE) and for dn = 1, 2, 3, 4, i.e. FD5

$ python generate_tables_of_finite_difference_schemes.py 6 4

# store all schemes in a dictionary (of a dictionary of a dictionary of...)

$ ipython

In [1]: from finite_difference_schemes import *
In [2]: FD_schemes = store_schemes(dn_max = 4)

# note: dn_max should be the same as max dn used in table generation, sys.argv[2]

# FD_schemes is a dictionary containing all schemes:
#
# FD_schemes['dn#'][handedness][asymmetry] = {'w' : [list of weights,
#                                        'stencil' : [list of stencil indices]}
#
# e.g., selecting the weights and stencil for a B2 scheme for f^(3):
#
# w = FD_schemes['dn3']['backward']['2']['w']
# stencil = FD_schemes['dn3']['backward']['2']['stencil']

# Example of stencils: LTE p = 6 

# p + dn = 6 grid points needed
# grid
# x-----x-----x-----x-----x-----x-----x

for different degrees of asymmetry choose point to be sampled the 0 below, 
then we have a family of schemes for each degree of asymmetry from each class of
scheme: 'forward', 'backward', 'central', which are measured deviations from the
base scheme, per example:

# forward + 0 = F0 scheme: (asymmetry = 0)
# 0-----x-----x-----x-----x-----x-----x

# forward - 1 = F1 scheme: (asymmetry = 1)
# x-----0-----x-----x-----x-----x-----x

# forward - 2 = F2 scheme: (asymmetry = 2)
# x-----x-----0-----x-----x-----x-----x

# central = C0 scheme: (asymmetry = 0)
# x-----x-----x-----0-----x-----x-----x

# backward + 2 = B2 scheme: (asymmetry = 2)
# x-----x-----x-----x-----0-----x-----x

# backward + 1 = B1 scheme: (asymmetry = 1)
# x-----x-----x-----x-----x-----0-----x

# backward + 0 = B0 scheme: (asymmetry = 0)
# x-----x-----x-----x-----x-----x-----0


# f'' O(6) for all regions
# grid
# x-----x-----x-----x-----x-----x-----x-----x

This time, we will label indices of the stencil
centered at the point 0

# forward + 0 = F0 scheme: (asymmetry = 0)
# 0-----1-----2-----3-----4-----5-----6-----7

# forward - 1 = F1 scheme: (asymmetry = 1)
# [-1]----0-----1-----2-----3-----4-----5-----6

# forward - 2 = F2 scheme: (asymmetry = 2)
# [-2]---[-1]---0-----1-----2-----3-----4-----5

# forward - 3 = F3 scheme: (asymmetry = 3)
# [-3]---[-2]--[-1]---0-----1-----2-----3-----4

# backward + 3 = B3 scheme: (asymmetry = 3)
# [-4]---[-3]--[-2]---[-1]-----0-----1-----2----3

# backward + 2 = B2 scheme: (asymmetry = 2)
# [-5]---[-4]--[-3]----[-2]---[-1]---0-----1----2

# backward + 1 = B1 scheme: (asymmetry = 1)
# [-6]---[-5]--[-4]----[-3]---[-2]---[-1]----0----1

# backward + 0 = B0 scheme (asymmetry = 0):
# [-7]---[-6]--[-5]---[-4]---[-3]---[-2]----[-1]----0


