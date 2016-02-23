import numpy as np
import scipy.misc

fact = scipy.misc.factorial

def beta_0(a):

    if a >= 0:
        return a
    else:
        return a

def beta_1(a):

    if a >= 0:
        return 1/2. * (a ** 2 - a)
    else:
        return 1/2. * (a ** 2 + a)

def beta_2(a):

    if a >= 0:
        return a ** 3 / 6. - a ** 2 / 4. + a / 12.
    else:
        return a ** 3 / 6. + a ** 2 / 4. + a / 12.

def beta_3(a):

    if a >= 0:
        return a ** 4 / 24. - a ** 3 / 12. + a ** 2 / 24
    else:
        return a ** 4 / 24. + a ** 3 / 12. + a ** 2 / 24

def beta_4(a):

    if a >= 0:
        return a ** 5 / 120. - a ** 4 / 48. + a ** 3 / 72. - a / 720.
    else:
        return a ** 5 / 120. + a ** 4 / 48. + a ** 3 / 72. -  a / 720.

def beta_5(a):

    if a >= 0:
        return a ** 6 / 720. - a ** 5 /240. + a ** 4 / 288. - a ** 2 / 1440
    else:
        return a ** 6 / 720. + a ** 5 /240. + a ** 4 / 288. - a ** 2 / 1440


