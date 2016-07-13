import DECSKS
import numpy as np



def setup(fe, fi, x, vx):
    """
    see notebook s29 for details

    Here we take the case that ghost points are required at both boundaries
    in x, but not in vx
    
    """
