import numpy as np

class phasespace_var:
    def __init__(self, f_old):
        self.MCs = np.zeros_like(f_old)
        self.prepointvaluemesh = np.random.randn(f_old.shape[0],f_old.shape[1])
        self.width = .5
        self.CFL = CourantNumber(f_old)

class CourantNumber:
    """Returns a CFL number instance of the phasespace variable z

    inputs:
    z -- (instance) phase space variable from class lib.domain.Setup

    outputs:
    self -- (instance) CFL number ascribed to variable z convection

    Note: broadcasting ensures self.numbers is the same shape as z.MCs
    """
    def __init__(self, f_old):

        self.numbers = np.zeros_like(f_old)
        self.frac = np.zeros_like(f_old)
        self.int = np.zeros_like(f_old)


    def compute_numbers(self, z, vz, dt):
        """Calculates the CFL numbers and corresponding integer and fractional
        parts, call by z.CFL.calculate_numbers(*args)

        inputs:
        self -- (lib.domain.CourantNumber instance) CFL instance with attribute containers
                containers CFL.numbers, CFL.int, CFL.frac.

                NOTE: the method is an attribute belonging to the subinstance z.CFL
                hence, self refers to z.CFL, not z.

        z -- (lib.domain.Setup instance) phasespace instance being advected
        vz -- (lib.domain.Setup instance) velocity for self
        dt -- (float) width of time step, can be a fraction of t.width for split
              methods

        outputs:
        None -- updates attributes
        """
        self.numbers = vz.prepointvaluemesh * dt / z.width

        # if >= 0 , self.int = floor(self.numbers), else ceil(self.numbers)
        # i.e. the sign of CFL.frac agrees with the sign of vz
        self.int = np.where(self.numbers >=0, np.floor(self.numbers),
                            np.ceil(self.numbers))

        # remaining portion is the fractional CFL number
        self.frac = self.numbers - self.int
