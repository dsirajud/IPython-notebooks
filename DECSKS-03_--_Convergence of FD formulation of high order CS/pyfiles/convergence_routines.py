import _mypath
import pyfiles
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA

def function(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_c = -0.25):

    return a*np.exp(-((_x + x_a) / w_a) ** 2) + np.exp(-((_x/w_b)**2)) + c*np.exp(-((_x + x_c)/w_c) ** 2)

def function_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a*np.exp(-((_x + x_a) / w_a) ** 2) + b*np.exp(-(((_x + x_b)/w_b)**2)) + c*np.exp(-((_x + x_c)/w_c) ** 2) \
      + d*np.exp(-((_x + x_d) / w_d) ** 2)

def df1(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0, x_c = -0.25):

    return a*(-2./w_a**2)*(_x + x_a)*np.exp(-((_x + x_a) / w_a)**2) \
      + (-2./w_b**2)*_x*np.exp(-(_x/w_b)**2) \
      + c*(-2./w_c**2)*(_x + x_c)*np.exp(-((_x + x_c)/w_c)**2)

def df2(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0, x_c = -0.25):

    return a*(-2./w_a**2) * (1 - 2./w_a**2 * (_x + x_a)**2) * np.exp(-( (_x + x_a) / w_a)**2) \
      + (-2./w_b**2) * (1 - 2./w_b**2 * _x ** 2) * np.exp(-(_x / w_b)**2) \
      + c*(-2./w_c**2) * (1 - 2/w_c**2 * (_x + x_c)**2) * np.exp(-( (_x + x_c) / w_c) ** 2)

def df3(_x, a = 3/4., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0, x_c = -0.25):

    return a * (-2. / w_a**2) ** 2 * ( 3*(_x + x_a) - 2./w_a**2 * (_x + x_a) ** 3) * np.exp(-( (_x + x_a) / w_a) ** 2) \
      + (-2./w_b**2)**2 * (3*_x - 2./w_b**2 * _x**3) * np.exp(-(_x / w_b)**2) \
      + c*(-2./w_c**2) ** 2 * (3*(_x + x_c) - 2./w_c**2 * (_x + x_c)**3) * np.exp(-( (_x + x_c) / w_c) ** 2)

def df4(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a ** 2) ** 2 * ( 3*(1 + 2*(-2./w_a**2)*(_x + x_a) ** 2) + (-2. / w_a**2)**2 * (_x + x_a) ** 4) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
            + b * (-2./w_b ** 2) ** 2 * ( 3*(1 + 2*(-2./w_b**2)*(_x + x_b) ** 2) + (-2. / w_b**2)**2 * (_x + x_b) ** 4) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
            + c * (-2./w_c ** 2) ** 2 * ( 3*(1 + 2*(-2./w_c**2)*(_x + x_c) ** 2) + (-2. / w_c**2)**2 * (_x + x_c) ** 4) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df5(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a ** 2) ** 3 * ( 15*(_x + x_a) + 10*(-2./w_a**2) * (_x + x_a) ** 3 + (-2./w_a**2) ** 2 *  (_x + x_a)**5) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
      + b * (-2./w_b ** 2) ** 3 * ( 15*(_x + x_b) + 10*(-2./w_b**2) * (_x + x_b) ** 3 + (-2./w_b**2) ** 2 *  (_x + x_b)**5) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
      + c * (-2./w_c ** 2) ** 3 * ( 15*(_x + x_c) + 10*(-2./w_c**2) * (_x + x_c) ** 3 + (-2./w_c**2) ** 2 *  (_x + x_c)**5) * np.exp(-( (_x + x_c) / w_c) ** 2 )

def df6(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a**2) ** 3 * (15
                                    + 45 * (-2./w_a**2) * (_x + x_a)**2
                                    + 15 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + (-2./w_a**2) ** 3 * (_x + x_a)**6) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 3 * (15
                                    + 45 * (-2./w_b**2) * (_x + x_b)**2
                                    + 15 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + (-2./w_b**2) ** 3 * (_x + x_b)**6) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 3 * (15
                                    + 45 * (-2./w_c**2) * (_x + x_c)**2
                                    + 15 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + (-2./w_c**2) ** 3 * (_x + x_c)**6) * np.exp(-( (_x + x_c) / w_c) ** 2 )

def df7(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a**2) ** 4 * (105 * (_x + x_a)
                                    + 105 * (-2./w_a**2) * (_x + x_a) ** 3
                                    + 21 * (-2./w_a**2) ** 2 * (_x + x_a) ** 5
                                    +      (-2./w_a**2) ** 3 * (_x + x_a) ** 7
                                    ) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 4 * (105 * (_x + x_b)
                                    + 105 * (-2./w_b**2) * (_x + x_b) ** 3
                                    + 21 * (-2./w_b**2) ** 2 * (_x + x_b) ** 5
                                    +      (-2./w_b**2) ** 3 * (_x + x_b) ** 7
                                    ) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 4 * (105 * (_x + x_c)
                                    + 105 * (-2./w_c**2) * (_x + x_c) ** 3
                                    + 21 * (-2./w_c**2) ** 2 * (_x + x_c) ** 5
                                    +      (-2./w_c**2) ** 3 * (_x + x_c) ** 7
                                    ) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df8(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a**2) ** 4 * (105
                                    + 420 * (-2./w_a**2) * (_x + x_a)**2
                                    + 210 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + 28 * (-2./w_a**2) ** 3 * (_x + x_a)**6
                                    + (-2./w_a**2) ** 4 * (_x + x_a)**8) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 4 * (105
                                    + 420 * (-2./w_b**2) * (_x + x_b)**2
                                    + 210 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + 28 * (-2./w_b**2) ** 3 * (_x + x_b)**6
                                    + (-2./w_b**2) ** 4 * (_x + x_b)**8) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 4 * (105
                                    + 420 * (-2./w_c**2) * (_x + x_c)**2
                                    + 210 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + 28 * (-2./w_c**2) ** 3 * (_x + x_c)**6
                                    + (-2./w_c**2) ** 4 * (_x + x_c)**8) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df9(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a**2) ** 5 * (945 * (_x + x_a)
                                    + 1260 * (-2./w_a**2) * (_x + x_a) ** 3
                                    + 378 * (-2./w_a**2) ** 2 * (_x + x_a) ** 5
                                    + 36 *  (-2./w_a**2) ** 3 * (_x + x_a) ** 7
                                    + (-2./w_a**2) ** 4 * (_x + x_a) ** 9
                                    ) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 5 * (945 * (_x + x_b)
                                    + 1260 * (-2./w_b**2) * (_x + x_b) ** 3
                                    + 378 * (-2./w_b**2) ** 2 * (_x + x_b) ** 5
                                    + 36 *  (-2./w_b**2) ** 3 * (_x + x_b) ** 7
                                    + (-2./w_b**2) ** 4 * (_x + x_b) ** 9
                                    ) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 5 * (945 * (_x + x_c)
                                    + 1260 * (-2./w_c**2) * (_x + x_c) ** 3
                                    + 378 * (-2./w_c**2) ** 2 * (_x + x_c) ** 5
                                    + 36 *  (-2./w_c**2) ** 3 * (_x + x_c) ** 7
                                    + (-2./w_c**2) ** 4 * (_x + x_c) ** 9
                                    ) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df10(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):


    return a * (-2./w_a**2) ** 5 * (945
                                    + 4725 * (-2./w_a**2) * (_x + x_a)**2
                                    + 3150 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + 630 * (-2./w_a**2) ** 3 * (_x + x_a)**6
                                    + 45 * (-2./w_a**2) ** 4 * (_x + x_a)**8
                                    + (-2./w_a**2) ** 5 * (_x + x_a) ** 10)* np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 4 * (945
                                    + 4725 * (-2./w_b**2) * (_x + x_b)**2
                                    + 3150 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + 630 * (-2./w_b**2) ** 3 * (_x + x_b)**6
                                    + 45 * (-2./w_b**2) ** 4 * (_x + x_b)**8
                                    + (-2./w_b**2) ** 5 * (_x + x_b) ** 10)* np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 4 * (945
                                    + 4725 * (-2./w_c**2) * (_x + x_c)**2
                                    + 3150 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + 630 * (-2./w_c**2) ** 3 * (_x + x_c)**6
                                    + 45 * (-2./w_c**2) ** 4 * (_x + x_c)**8
                                    + (-2./w_c**2) ** 5 * (_x + x_c) ** 10)* np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df11(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a**2) ** 6 * (10395 * (_x + x_a)
                                    + 17325 * (-2./w_a**2) * (_x + x_a) ** 3
                                    + 6930 * (-2./w_a**2) ** 2 * (_x + x_a) ** 5
                                    + 990 *  (-2./w_a**2) ** 3 * (_x + x_a) ** 7
                                    + 55 * (-2./w_a**2) ** 4 * (_x + x_a) ** 9
                                    + (-2./w_a**2) ** 5 * (_x + x_a) ** 11
                                    ) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 6 * (10395 * (_x + x_b)
                                    + 17325 * (-2./w_b**2) * (_x + x_b) ** 3
                                    + 6930 * (-2./w_b**2) ** 2 * (_x + x_b) ** 5
                                    + 990 *  (-2./w_b**2) ** 3 * (_x + x_b) ** 7
                                    + 55 * (-2./w_b**2) ** 4 * (_x + x_b) ** 9
                                    + (-2./w_b**2) ** 5 * (_x + x_b) ** 11
                                    ) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 6 * (10395 * (_x + x_c)
                                    + 17325 * (-2./w_c**2) * (_x + x_c) ** 3
                                    + 6930 * (-2./w_c**2) ** 2 * (_x + x_c) ** 5
                                    + 990 *  (-2./w_c**2) ** 3 * (_x + x_c) ** 7
                                    + 55 * (-2./w_c**2) ** 4 * (_x + x_c) ** 9
                                    + (-2./w_c**2) ** 5 * (_x + x_c) ** 11
                                    ) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df12(_x, a = 3/4., b = 1., c = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, x_a = 0.25, x_b = 0., x_c = -0.25):

    return a * (-2./w_a**2) ** 6 * (10395
                                    + 62370 * (-2./w_a**2) * (_x + x_a)**2
                                    + 51975 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + 13860 * (-2./w_a**2) ** 3 * (_x + x_a)**6
                                    + 1485 * (-2./w_a**2) ** 4 * (_x + x_a)**8
                                    + 66 * (-2./w_a**2) ** 5 * (_x + x_a) ** 10
                                    + (-2./w_a**2) ** 6 * (_x + x_a) ** 12)* np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 6 * (10395
                                    + 62370 * (-2./w_b**2) * (_x + x_b)**2
                                    + 51975 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + 13860 * (-2./w_b**2) ** 3 * (_x + x_b)**6
                                    + 1485 * (-2./w_b**2) ** 4 * (_x + x_b)**8
                                    + 66 * (-2./w_b**2) ** 5 * (_x + x_b) ** 10
                                    + (-2./w_b**2) ** 6 * (_x + x_b) ** 12)* np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 6 * (10395
                                    + 62370 * (-2./w_c**2) * (_x + x_c)**2
                                    + 51975 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + 13860 * (-2./w_c**2) ** 3 * (_x + x_c)**6
                                    + 1485 * (-2./w_c**2) ** 4 * (_x + x_c)**8
                                    + 66 * (-2./w_c**2) ** 5 * (_x + x_c) ** 10
                                    + (-2./w_c**2) ** 6 * (_x + x_c) ** 12)* np.exp(-( (_x + x_c) / w_c) ** 2 ) \

def df1_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a*(-2./w_a**2)*(_x + x_a)*np.exp(-((_x + x_a) / w_a)**2) \
      + b*(-2./w_b**2)*(_x + x_b)*np.exp(-((_x + x_b)/w_b)**2) \
      + c*(-2./w_c**2)*(_x + x_c)*np.exp(-((_x + x_c)/w_c)**2) \
      + d*(-2./w_d**2)*(_x + x_d)*np.exp(-((_x + x_d)/w_d)**2)

def df2_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a*(-2./w_a**2) * (1 - 2./w_a**2 * (_x + x_a)**2) * np.exp(-( (_x + x_a) / w_a)**2) \
      + b*(-2./w_b**2) * (1 - 2./w_b**2 * (_x + x_b) ** 2) * np.exp(-((_x + x_b) / w_b)**2) \
      + c*(-2./w_c**2) * (1 - 2/w_c**2 * (_x + x_c)**2) * np.exp(-( (_x + x_c) / w_c) ** 2)\
      + d*(-2./w_d**2) * (1 - 2/w_d**2 * (_x + x_d)**2) * np.exp(-( (_x + x_d) / w_d) ** 2)

def df3_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2. / w_a**2) ** 2 * ( 3*(_x + x_a) - 2./w_a**2 * (_x + x_a) ** 3) * np.exp(-( (_x + x_a) / w_a) ** 2) \
      + b*(-2./w_b**2)**2 * (3*(_x + x_b) - 2./w_b**2 * (_x + x_b)**3) * np.exp(-((_x + x_b) / w_b)**2) \
      + c*(-2./w_c**2) ** 2 * (3*(_x + x_c) - 2./w_c**2 * (_x + x_c)**3) * np.exp(-( (_x + x_c) / w_c) ** 2)\
      + d*(-2./w_d**2) ** 2 * (3*(_x + x_d) - 2./w_d**2 * (_x + x_d)**3) * np.exp(-( (_x + x_d) / w_d) ** 2)

def df4_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a ** 2) ** 2 * ( 3*(1 + 2*(-2./w_a**2)*(_x + x_a) ** 2) + (-2. / w_a**2)**2 * (_x + x_a) ** 4) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
            + b * (-2./w_b ** 2) ** 2 * ( 3*(1 + 2*(-2./w_b**2)*(_x + x_b) ** 2) + (-2. / w_b**2)**2 * (_x + x_b) ** 4) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
            + c * (-2./w_c ** 2) ** 2 * ( 3*(1 + 2*(-2./w_c**2)*(_x + x_c) ** 2) + (-2. / w_c**2)**2 * (_x + x_c) ** 4) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
            + d * (-2./w_d ** 2) ** 2 * ( 3*(1 + 2*(-2./w_d**2)*(_x + x_d) ** 2) + (-2. / w_d**2)**2 * (_x + x_d) ** 4) * np.exp(-( (_x + x_d) / w_d) ** 2 )

def df5_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a ** 2) ** 3 * ( 15*(_x + x_a) + 10*(-2./w_a**2) * (_x + x_a) ** 3 + (-2./w_a**2) ** 2 *  (_x + x_a)**5) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
      + b * (-2./w_b ** 2) ** 3 * ( 15*(_x + x_b) + 10*(-2./w_b**2) * (_x + x_b) ** 3 + (-2./w_b**2) ** 2 *  (_x + x_b)**5) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
      + c * (-2./w_c ** 2) ** 3 * ( 15*(_x + x_c) + 10*(-2./w_c**2) * (_x + x_c) ** 3 + (-2./w_c**2) ** 2 *  (_x + x_c)**5) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
      + d * (-2./w_d ** 2) ** 3 * ( 15*(_x + x_d) + 10*(-2./w_d**2) * (_x + x_d) ** 3 + (-2./w_d**2) ** 2 *  (_x + x_d)**5) * np.exp(-( (_x + x_d) / w_d) ** 2 )

def df6_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a**2) ** 3 * (15
                                    + 45 * (-2./w_a**2) * (_x + x_a)**2
                                    + 15 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + (-2./w_a**2) ** 3 * (_x + x_a)**6) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 3 * (15
                                    + 45 * (-2./w_b**2) * (_x + x_b)**2
                                    + 15 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + (-2./w_b**2) ** 3 * (_x + x_b)**6) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 3 * (15
                                    + 45 * (-2./w_c**2) * (_x + x_c)**2
                                    + 15 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + (-2./w_c**2) ** 3 * (_x + x_c)**6) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
                                    \
                                    + d * (-2./w_d**2) ** 3 * (15
                                    + 45 * (-2./w_d**2) * (_x + x_d)**2
                                    + 15 * (-2./w_d**2)**2 * (_x + x_d)**4
                                    + (-2./w_d**2) ** 3 * (_x + x_d)**6) * np.exp(-( (_x + x_d) / w_d) ** 2 )


def df7_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a**2) ** 4 * (105 * (_x + x_a)
                                    + 105 * (-2./w_a**2) * (_x + x_a) ** 3
                                    + 21 * (-2./w_a**2) ** 2 * (_x + x_a) ** 5
                                    +      (-2./w_a**2) ** 3 * (_x + x_a) ** 7
                                    ) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 4 * (105 * (_x + x_b)
                                    + 105 * (-2./w_b**2) * (_x + x_b) ** 3
                                    + 21 * (-2./w_b**2) ** 2 * (_x + x_b) ** 5
                                    +      (-2./w_b**2) ** 3 * (_x + x_b) ** 7
                                    ) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 4 * (105 * (_x + x_c)
                                    + 105 * (-2./w_c**2) * (_x + x_c) ** 3
                                    + 21 * (-2./w_c**2) ** 2 * (_x + x_c) ** 5
                                    +      (-2./w_c**2) ** 3 * (_x + x_c) ** 7
                                    ) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
                                    \
                                    + d * (-2./w_d**2) ** 4 * (105 * (_x + x_d)
                                    + 105 * (-2./w_d**2) * (_x + x_d) ** 3
                                    + 21 * (-2./w_d**2) ** 2 * (_x + x_d) ** 5
                                    +      (-2./w_d**2) ** 3 * (_x + x_d) ** 7
                                    ) * np.exp(-( (_x + x_d) / w_d) ** 2 ) \

def df8_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a**2) ** 4 * (105
                                    + 420 * (-2./w_a**2) * (_x + x_a)**2
                                    + 210 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + 28 * (-2./w_a**2) ** 3 * (_x + x_a)**6
                                    + (-2./w_a**2) ** 4 * (_x + x_a)**8) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 4 * (105
                                    + 420 * (-2./w_b**2) * (_x + x_b)**2
                                    + 210 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + 28 * (-2./w_b**2) ** 3 * (_x + x_b)**6
                                    + (-2./w_b**2) ** 4 * (_x + x_b)**8) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 4 * (105
                                    + 420 * (-2./w_c**2) * (_x + x_c)**2
                                    + 210 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + 28 * (-2./w_c**2) ** 3 * (_x + x_c)**6
                                    + (-2./w_c**2) ** 4 * (_x + x_c)**8) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
                                    + d * (-2./w_d**2) ** 4 * (105
                                    + 420 * (-2./w_d**2) * (_x + x_d)**2
                                    + 210 * (-2./w_d**2)**2 * (_x + x_d)**4
                                    + 28 * (-2./w_d**2) ** 3 * (_x + x_d)**6
                                    + (-2./w_d**2) ** 4 * (_x + x_d)**8) * np.exp(-( (_x + x_d) / w_d) ** 2 ) \

def df9_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a**2) ** 5 * (945 * (_x + x_a)
                                    + 1260 * (-2./w_a**2) * (_x + x_a) ** 3
                                    + 378 * (-2./w_a**2) ** 2 * (_x + x_a) ** 5
                                    + 36 *  (-2./w_a**2) ** 3 * (_x + x_a) ** 7
                                    + (-2./w_a**2) ** 4 * (_x + x_a) ** 9
                                    ) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 5 * (945 * (_x + x_b)
                                    + 1260 * (-2./w_b**2) * (_x + x_b) ** 3
                                    + 378 * (-2./w_b**2) ** 2 * (_x + x_b) ** 5
                                    + 36 *  (-2./w_b**2) ** 3 * (_x + x_b) ** 7
                                    + (-2./w_b**2) ** 4 * (_x + x_b) ** 9
                                    ) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 5 * (945 * (_x + x_c)
                                    + 1260 * (-2./w_c**2) * (_x + x_c) ** 3
                                    + 378 * (-2./w_c**2) ** 2 * (_x + x_c) ** 5
                                    + 36 *  (-2./w_c**2) ** 3 * (_x + x_c) ** 7
                                    + (-2./w_c**2) ** 4 * (_x + x_c) ** 9
                                    ) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
                                    \
                                    + d * (-2./w_d**2) ** 5 * (945 * (_x + x_d)
                                    + 1260 * (-2./w_d**2) * (_x + x_d) ** 3
                                    + 378 * (-2./w_d**2) ** 2 * (_x + x_d) ** 5
                                    + 36 *  (-2./w_d**2) ** 3 * (_x + x_d) ** 7
                                    + (-2./w_d**2) ** 4 * (_x + x_d) ** 9
                                    ) * np.exp(-( (_x + x_d) / w_d) ** 2 ) \

def df10_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a**2) ** 5 * (945
                                    + 4725 * (-2./w_a**2) * (_x + x_a)**2
                                    + 3150 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + 630 * (-2./w_a**2) ** 3 * (_x + x_a)**6
                                    + 45 * (-2./w_a**2) ** 4 * (_x + x_a)**8
                                    + (-2./w_a**2) ** 5 * (_x + x_a) ** 10)* np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 4 * (945
                                    + 4725 * (-2./w_b**2) * (_x + x_b)**2
                                    + 3150 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + 630 * (-2./w_b**2) ** 3 * (_x + x_b)**6
                                    + 45 * (-2./w_b**2) ** 4 * (_x + x_b)**8
                                    + (-2./w_b**2) ** 5 * (_x + x_b) ** 10)* np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 4 * (945
                                    + 4725 * (-2./w_c**2) * (_x + x_c)**2
                                    + 3150 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + 630 * (-2./w_c**2) ** 3 * (_x + x_c)**6
                                    + 45 * (-2./w_c**2) ** 4 * (_x + x_c)**8
                                    + (-2./w_c**2) ** 5 * (_x + x_c) ** 10)* np.exp(-( (_x + x_c) / w_c) ** 2 ) \
                                    \
                                    + d * (-2./w_d**2) ** 4 * (945
                                    + 4725 * (-2./w_d**2) * (_x + x_d)**2
                                    + 3150 * (-2./w_d**2)**2 * (_x + x_d)**4
                                    + 630 * (-2./w_d**2) ** 3 * (_x + x_d)**6
                                    + 45 * (-2./w_d**2) ** 4 * (_x + x_d)**8
                                    + (-2./w_d**2) ** 5 * (_x + x_d) ** 10)* np.exp(-( (_x + x_d) / w_d) ** 2 )



def df11_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a**2) ** 6 * (10395 * (_x + x_a)
                                    + 17325 * (-2./w_a**2) * (_x + x_a) ** 3
                                    + 6930 * (-2./w_a**2) ** 2 * (_x + x_a) ** 5
                                    + 990 *  (-2./w_a**2) ** 3 * (_x + x_a) ** 7
                                    + 55 * (-2./w_a**2) ** 4 * (_x + x_a) ** 9
                                    + (-2./w_a**2) ** 5 * (_x + x_a) ** 11
                                    ) * np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 6 * (10395 * (_x + x_b)
                                    + 17325 * (-2./w_b**2) * (_x + x_b) ** 3
                                    + 6930 * (-2./w_b**2) ** 2 * (_x + x_b) ** 5
                                    + 990 *  (-2./w_b**2) ** 3 * (_x + x_b) ** 7
                                    + 55 * (-2./w_b**2) ** 4 * (_x + x_b) ** 9
                                    + (-2./w_b**2) ** 5 * (_x + x_b) ** 11
                                    ) * np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 6 * (10395 * (_x + x_c)
                                    + 17325 * (-2./w_c**2) * (_x + x_c) ** 3
                                    + 6930 * (-2./w_c**2) ** 2 * (_x + x_c) ** 5
                                    + 990 *  (-2./w_c**2) ** 3 * (_x + x_c) ** 7
                                    + 55 * (-2./w_c**2) ** 4 * (_x + x_c) ** 9
                                    + (-2./w_c**2) ** 5 * (_x + x_c) ** 11
                                    ) * np.exp(-( (_x + x_c) / w_c) ** 2 ) \
                                    \
                                    + d * (-2./w_d**2) ** 6 * (10395 * (_x + x_d)
                                    + 17325 * (-2./w_d**2) * (_x + x_d) ** 3
                                    + 6930 * (-2./w_d**2) ** 2 * (_x + x_d) ** 5
                                    + 990 *  (-2./w_d**2) ** 3 * (_x + x_d) ** 7
                                    + 55 * (-2./w_d**2) ** 4 * (_x + x_d) ** 9
                                    + (-2./w_d**2) ** 5 * (_x + x_d) ** 11
                                    ) * np.exp(-( (_x + x_d) / w_d) ** 2 )

def df12_edge(_x, a = 3/4., b = 1., c = 1/2., d = 1/2., w_a = 0.03, w_b = 0.06, w_c = 0.1, w_d = 0.1,
             x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    return a * (-2./w_a**2) ** 6 * (10395
                                    + 62370 * (-2./w_a**2) * (_x + x_a)**2
                                    + 51975 * (-2./w_a**2)**2 * (_x + x_a)**4
                                    + 13860 * (-2./w_a**2) ** 3 * (_x + x_a)**6
                                    + 1485 * (-2./w_a**2) ** 4 * (_x + x_a)**8
                                    + 66 * (-2./w_a**2) ** 5 * (_x + x_a) ** 10
                                    + (-2./w_a**2) ** 6 * (_x + x_a) ** 12)* np.exp(-( (_x + x_a) / w_a) ** 2 ) \
                                    \
                                    + b * (-2./w_b**2) ** 6 * (10395
                                    + 62370 * (-2./w_b**2) * (_x + x_b)**2
                                    + 51975 * (-2./w_b**2)**2 * (_x + x_b)**4
                                    + 13860 * (-2./w_b**2) ** 3 * (_x + x_b)**6
                                    + 1485 * (-2./w_b**2) ** 4 * (_x + x_b)**8
                                    + 66 * (-2./w_b**2) ** 5 * (_x + x_b) ** 10
                                    + (-2./w_b**2) ** 6 * (_x + x_b) ** 12)* np.exp(-( (_x + x_b) / w_b) ** 2 ) \
                                    \
                                    + c * (-2./w_c**2) ** 6 * (10395
                                    + 62370 * (-2./w_c**2) * (_x + x_c)**2
                                    + 51975 * (-2./w_c**2)**2 * (_x + x_c)**4
                                    + 13860 * (-2./w_c**2) ** 3 * (_x + x_c)**6
                                    + 1485 * (-2./w_c**2) ** 4 * (_x + x_c)**8
                                    + 66 * (-2./w_c**2) ** 5 * (_x + x_c) ** 10
                                    + (-2./w_c**2) ** 6 * (_x + x_c) ** 12)* np.exp(-( (_x + x_c) / w_c) ** 2 ) \
                                    \
                                    + d * (-2./w_d**2) ** 6 * (10395
                                    + 62370 * (-2./w_d**2) * (_x + x_d)**2
                                    + 51975 * (-2./w_d**2)**2 * (_x + x_d)**4
                                    + 13860 * (-2./w_d**2) ** 3 * (_x + x_d)**6
                                    + 1485 * (-2./w_d**2) ** 4 * (_x + x_d)**8
                                    + 66 * (-2./w_d**2) ** 5 * (_x + x_d) ** 10
                                    + (-2./w_d**2) ** 6 * (_x + x_d) ** 12)* np.exp(-( (_x + x_d) / w_d) ** 2 ) \

def domain(_a = -0.5, _b = 1.5, _Nx = 21):

    _L = float(_b - _a)
    _dx = _L / _Nx

    _x = np.zeros(_Nx)
    for i in range(_Nx):
        _x[i] = _a + i*_dx

    return _x, _dx, _L


def FD_derivative_matrix_formulation(_dn = 1, _p = 1, _Nx = 128, edge = 'no',
                                     x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):

    x, dx, L = domain(_a = -0.5, _b = 1.5, _Nx = _Nx)
    if edge == 'yes':
        f = function_edge(x,
                           x_a = x_a, x_b = x_b, x_c = x_c, x_d = x_d)

        df_str = 'df' + str(_dn) + '_edge'
        df = eval(df_str)(x,
                       x_a = x_a, x_b = x_b, x_c = x_c, x_d = x_d)

    else:
        f = function(x)

        df_str = 'df' + str(_dn)
        df = eval(df_str)(x)

    # compute derivative approximations
    # extract subdictionary pertaining to scheme of derivative order dn
    dn_key = 'dn' + str(_dn)
    FD_schemes = pyfiles.lib.make_FD_schemes_dict.store(_dn, _p)

    FD_schemes = FD_schemes[dn_key] # dictionary containing the family of
                                    # FD schemes or order dn

    imax = len(x) - 1
    stencil_size = _p + _dn
    stencil_center = stencil_size // 2
    W = np.zeros([_Nx, _Nx])
    for i in range(_Nx):
        if i < stencil_center:
            handedness = 'forward'
            asymmetry = str(i)
        elif imax - i < stencil_center:
            handedness = 'backward'
            asymmetry = str(imax - i)
        else:
            if np.mod(stencil_size,2) == 1:
                handedness = 'central'
                asymmetry = str(0)
            else:
                handedness = 'forward'
                asymmetry = str(stencil_center - 1)
        FD_scheme = FD_schemes[handedness][asymmetry]
        w = FD_scheme['w']
        stencil = FD_scheme['stencil']

        W[i, i + np.array(stencil)] = w

    df_approx = W.dot(f)
    df_approx /= dx ** _dn
    error = df_approx - df
    L2_error_norm = np.sqrt(dx/L)* LA.norm(error, 2) # L2 norm of error
    Linf_error = np.max(error) # L-infinity norm of error

    return L2_error_norm, df_approx

def convergence_routine_nonperiodic(NumGrids = 10, Nx = 18, LTE = 3, dn = 1, plots = 'yes', edge = 'no',
                                     x_a = -1.025, x_b = -1.275, x_c = -1.525, x_d = 0.475):
    """executes FD_derivative_matrix_formulation(*args) on progressively finer grids.
    The first grid contains Nx points (must be at least as large as the stencil size
    = O(LTE) + order of derivative), the subsequent (NumGrids-1)
    grids double the number of points from the previous grid

    inputs:
    NumGrids -- (int) number of grids the derivative is to be evaluated on
    Nx -- (int) number of gridpoints on coarsest grid
    LTE -- (int) local truncation error of the derivatives being used
    dn -- (int) order of derivative to be evaluated

    NOTE: derivative tables must be generated for the required dn at chosen LTE
    by executing

   $ python generate_tables_of_finite_difference_schemes_for_a_given_LTE.py LTE dn

    from ./bin/ which stores the table in ./etc/
    """
    error_norm = np.zeros(NumGrids)
    grids = np.zeros(NumGrids)
    orders = np.zeros(NumGrids)

    # calculate first order then loop over remaining
    q = 0
    error_norm[q], df_approx = FD_derivative_matrix_formulation(dn, LTE, Nx, edge = edge,
                                                                x_a = x_a, x_b = x_b, x_c = x_c, x_d = x_d)
    grids[q] = Nx

    for q in range(1, NumGrids):
        Nx *= 2
        error_norm[q], df_approx = FD_derivative_matrix_formulation(dn, LTE, Nx, edge = edge,
                                                                    x_a = x_a, x_b = x_b, x_c = x_c, x_d = x_d)
        orders[q] = np.log2(error_norm[q-1] / error_norm[q])
        grids[q] = Nx

        #    order = -np.polyfit(np.log2(grids), np.log2(error_norm), 1)[0]
        #    print "slope of a linear fit = %g \n" % order

    if np.mod(dn + LTE, 2) == 1:
        print "a central differencing scheme exists: LTE between %d and %d expected if density far from edges\n" % (LTE, LTE + 1)
    else:
        print "no central differencing scheme exists: LTE of %d expected\n" % LTE

    print "order calculations at each refinement step:"

    for n in range(NumGrids):
        if n == 0:
            print "Nx%d        error = %g       ----" % (grids[n],error_norm[n])
        else:
            print "Nx%d        error = %g       order = %g" % (grids[n],error_norm[n],orders[n])

    print '\n'

    if plots.lower() == 'yes':
        fig, ax = plt.subplots()
        ax.loglog(grids, error_norm, '--o', markersize = 8, linewidth = 2)
        ax.set_xscale('log', basex=2)
        ax.set_yscale('log', basey=2)
        plt.xlabel(r'Number of grid points $N_x$ over $x\in [-0.5,1.5]$', fontsize = 16)
        plt.ylabel(r'$L^2$ error of $f_{\Delta x} - f_{exact}$', fontsize = 16)
        #        plt.text(np.max(grids), np.max(error_norm),'slope = -%g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
        plt.title(r'Error in $d^{%i}f/dx^{%i}$ using FD schemes at LTE = $O(\Delta x^{%i})$' % (dn,dn,LTE))
        plt.grid()

    return error_norm, grids


def convergence_for_several_derivatives_at_const_LTE(NumGrids = 10, _LTE = 3, _dn_min = 1, _dn_max = 4, plots = 'yes'):
    """runs a convergence routine on a chosen number of grids (NumGrids) that whose uniform
    spacing between grid points is halved in each subsequent run, for the range of derivatives
    dn_min <= dn <= dn_max, dn_min must be greater than or equal to 1"""

    error_histories = np.zeros([NumGrids, _dn_max+1]) # no values stored in [:,0] entry
    grid_histories = np.zeros([NumGrids, _dn_max+1])  # no values stored in [:,0] entry

    fig, ax = plt.subplots()

    for d in range(_dn_min, _dn_max + 1):
        print "derivative dn = %d: \n" % d

        error_histories[:,d], grid_histories[:,d] = convergence_routine_nonperiodic(NumGrids = 10,
                                                    Nx = 21, LTE = _LTE, dn = d, plots = plots)

        plot_label = '$dn = $' + str(d)
        ax.loglog(grid_histories[:,d], error_histories[:,d], '--o', markersize = 8, linewidth = 2, label = plot_label)
        ax.hold('on')
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    plt.xlabel(r'Number of grid points $N_x$ over $x\in [-1,1]$', fontsize = 16)
    plt.ylabel('$L^2$ error', fontsize = 16)


    # shrink frame of the plot to make room for a legend on the right side
    frame = ax.get_position() # position of plot center
    ax.set_position([frame.x0, frame.y0,
                     frame.width * 0.8, frame.height])

    # Place legend to the right of the shrunken frame
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
              fancybox=True, ncol=1)
    plt.grid()

    return None

def convergence_for_a_single_derivative_at_several_LTE(NumGrids = 10, _dn = 8, _LTE_min = 1, _LTE_max = 4, plots = 'yes'):
    """runs a convergence routine on a chosen number of grids (NumGrids) that whose uniform
    spacing between grid points is halved in each subsequent run, for the range of
    LTEs: _LTE_min <= LTE <= _LTE_max, _LTE_min must be greater than or equal to 1 for
    a single derivative or order dn"""

    error_histories = np.zeros([NumGrids, _LTE_max+1]) # no values stored in [:,0] entry
    grid_histories = np.zeros([NumGrids, _LTE_max+1])  # no values stored in [:,0] entry

    fig, ax = plt.subplots()

    for _LTE in range(_LTE_min, _LTE_max + 1):
        if np.mod(_LTE + _dn, 2) == 1:
            print "LTE = %d: central scheme exists (order between %d and %d expected) \n" % (_LTE, _LTE, _LTE + 1)
        else:
            print "LTE = %d: \n" % _LTE

        error_histories[:,_LTE], grid_histories[:, _LTE] = convergence_routine_nonperiodic(NumGrids = 10,
                                                    Nx = 21, LTE = _LTE, dn = _dn, plots = plots)

        plot_label = '$LTE = $' + str(_LTE)
        ax.loglog(grid_histories[:, _LTE], error_histories[:,_LTE], '--o', markersize = 8, linewidth = 2, label = plot_label)
        ax.hold('on')
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    plt.xlabel(r'Number of grid points $N_x$ over $x\in [-0.5,1.5]$', fontsize = 16)
    plt.ylabel('$L^2$ error', fontsize = 16)


    # shrink frame of the plot to make room for a legend on the right side
    frame = ax.get_position() # position of plot center
    ax.set_position([frame.x0, frame.y0,
                     frame.width * 0.8, frame.height])

    # Place legend to the right of the shrunken frame
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
              fancybox=True, ncol=1)
    plt.grid()

    return None

if __name__ == '__main__':
    """
    run from terminal as

        $ python main.py NumGrids Nx LTE dn

    sys.argv inputs:
    NumGrids -- (int) number of grids the derivative is to be evaluated on
    Nx -- (int) number of gridpoints on coarsest grid
    LTE -- (int) local truncation error of the derivatives being used
    dn -- (int) order of derivative to be evaluated

    NOTE: derivative tables must be generated for the required dn at chosen LTE
    by executing

   $ python generate_tables_of_finite_difference_schemes_for_a_given_LTE.py LTE dn

    from ./bin/ which stores the table in ./etc/
    """
    import sys

    NumGrids = int(sys.argv[1])
    Nx = int(sys.argv[2])
    LTE = int(sys.argv[3])
    dn = int(sys.argv[4])

    error_norm = np.zeros(NumGrids)
    grids = np.zeros(NumGrids)
    orders = np.zeros(NumGrids)
    p = LTE # relabeling is due to author's personal collocation with this symbol

    # calculate first order then loop over remaining
    q = 0
    error_norm[q] = FD_derivative_matrix_formulation(dn, p, Nx)
    grids[q] = Nx

    for q in range(1, NumGrids):
        Nx *= 2
        error_norm[q] = FD_derivative_matrix_formulation(dn, p, Nx)
        orders[q] = np.log2(error_norm[q-1] / error_norm[q])
        grids[q] = Nx

    order = -np.polyfit(np.log2(grids), np.log2(error_norm), 1)[0]
    print "slope of a linear fit = %g" % order

    print "order calculations at each refinement step:"

    for n in range(NumGrids):
        if n == 0:
            print "Nx%d        error = %g       ----" % (grids[n],error_norm[n])
        else:
            print "Nx%d        error = %g       order = %g" % (grids[n],error_norm[n],orders[n])

    fig, ax = plt.subplots()
    ax.loglog(grids, error_norm, '--o')
    #ax.hold('on')
    #ax.loglog(abscissa, order_line, '-b')
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    #ax.set_xlim(1, 1e4)
    #ax.set_ylim(1e-15, 1)
    plt.ylabel('$L^2$ norm of $f_{num} - f_{exact}$', fontsize = 14)
    plt.xlabel('number of gridpoints $N_x$', fontsize = 14)
    plt.text(np.max(grids), np.max(error_norm),'slope = -%g' % order, horizontalalignment = 'right', verticalalignment = 'top', fontsize = 16)
    plt.grid()
    plt.show()
